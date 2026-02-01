#!/bin/bash
#SBATCH --job-name=gatk_merge_gvcf
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err

set -euo pipefail

usage() {
  echo "Usage: $0 -g <gvcf_dir> -r <reference_fasta> [-l interval] [-c cutoff]"
  echo "  -g  Directory containing per-sample .gvcf.gz files"
  echo "  -r  Reference FASTA (used for GATK and indexing)"
  echo "  -l  Optional interval (e.g., chr1) passed to GenomicsDBImport/GenotypeGVCFs"
  echo "  -c  Optional indel length cutoff passed to dropSV.py (default: dropSV.py default)"
  exit 1
}

if ! command -v module >/dev/null 2>&1; then
  echo "ERROR: environment modules not available (module command not found)."
  exit 1
fi

LOG_DIR="${SLURM_SUBMIT_DIR:-.}/logs"
mkdir -p "$LOG_DIR"

module load picard || { echo "ERROR: failed to load picard module"; exit 1; }
module load tabix || { echo "ERROR: failed to load tabix module"; exit 1; }
module load gatk || { echo "ERROR: failed to load gatk module"; exit 1; }

echo "Tool versions:"
picard --version || true
tabix --version || true
gatk --version || true

GVCF_DIR=""
REF_FASTA=""
INTERVAL=""
DROP_CUTOFF=""

while getopts "g:r:l:c:" opt; do
  case "$opt" in
    g) GVCF_DIR="$OPTARG" ;;
    r) REF_FASTA="$OPTARG" ;;
    l) INTERVAL="$OPTARG" ;;
    c) DROP_CUTOFF="$OPTARG" ;;
    *) usage ;;
  esac
done

if [ -z "$GVCF_DIR" ] || [ -z "$REF_FASTA" ]; then
  usage
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DROP_SV="${SCRIPT_DIR}/dropSV.py"

if [ ! -x "$DROP_SV" ]; then
  echo "ERROR: dropSV.py not found or not executable at $DROP_SV"
  exit 1
fi

if [ ! -f "$REF_FASTA" ]; then
  echo "ERROR: reference FASTA not found: $REF_FASTA"
  exit 1
fi

command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH"; exit 1; }
command -v gatk >/dev/null 2>&1 || { echo "ERROR: gatk not found in PATH"; exit 1; }
command -v bgzip >/dev/null 2>&1 || { echo "ERROR: bgzip not found in PATH"; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "ERROR: tabix not found in PATH"; exit 1; }
command -v picard >/dev/null 2>&1 || { echo "ERROR: picard not found in PATH"; exit 1; }

# Ensure reference index and dictionary exist.
if [ ! -f "${REF_FASTA}.fai" ]; then
  samtools faidx "$REF_FASTA"
fi

DICT="${REF_FASTA%.*}.dict"
if [ ! -f "$DICT" ]; then
  picard CreateSequenceDictionary R="$REF_FASTA" O="$DICT"
fi

# Run dropSV.py on the input directory.
if [ -n "$DROP_CUTOFF" ]; then
  python3 "$DROP_SV" -d "$GVCF_DIR" -c "$DROP_CUTOFF"
else
  python3 "$DROP_SV" -d "$GVCF_DIR"
fi

CLEAN_DIR="${GVCF_DIR%/}/cleangVCF"
if [ ! -d "$CLEAN_DIR" ]; then
  echo "ERROR: expected cleaned gVCFs in $CLEAN_DIR"
  exit 1
fi

# Ensure all cleaned gVCFs are bgzipped + indexed.
for f in "$CLEAN_DIR"/*.gvcf; do
  [ -e "$f" ] || continue
  bgzip -c "$f" > "${f}.gz"
  tabix -p vcf "${f}.gz"
done

for f in "$CLEAN_DIR"/*.gvcf.gz; do
  [ -e "$f" ] || continue
  if [ ! -f "${f}.tbi" ]; then
    tabix -p vcf "$f"
  fi
done

GVCF_FILES=()
while IFS= read -r line; do
  GVCF_FILES+=("$line")
done < <(find "$CLEAN_DIR" -maxdepth 1 -type f -name "*.gvcf.gz" | sort)
if [ "${#GVCF_FILES[@]}" -eq 0 ]; then
  echo "ERROR: no .gvcf.gz files found in $CLEAN_DIR"
  exit 1
fi

WORKSPACE="genomicsdb_workspace"
OUT_GVCF="${PWD}/combined.gvcf.gz"

GVCF_ARGS=()
for f in "${GVCF_FILES[@]}"; do
  GVCF_ARGS+=( -V "$f" )
done

IMPORT_CMD=(gatk --java-options "-Xmx100g -Xms100g" GenomicsDBImport "${GVCF_ARGS[@]}" \
  --genomicsdb-workspace-path "$WORKSPACE")
GENO_CMD=(gatk --java-options "-Xmx100g -Xms100g" GenotypeGVCFs \
  -R "$REF_FASTA" -V "gendb://${WORKSPACE}" -O "$OUT_GVCF")

if [ -n "$INTERVAL" ]; then
  IMPORT_CMD+=( -L "$INTERVAL" )
  GENO_CMD+=( -L "$INTERVAL" )
fi

"${IMPORT_CMD[@]}"
"${GENO_CMD[@]}"

echo "Combined gVCF written to $OUT_GVCF"
