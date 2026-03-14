# ARG Pipeline (Direct MAF)

This repository provides a Snakemake workflow for processing AnchorWave MAFs directly into per-contig site outputs. The workflow emits all-sites VCFs, variant-only VCFs, and BED masks from the alignments. Written with the assistance of [Codex](https://openai.com/codex/).

## Requirements

- Conda
- The environment defined in `argprep.yml`

## Setup

```bash
conda env create -f argprep.yml
conda activate argprep
```

## Configure

Edit `config.yaml` and set:

- `maf_dir`: directory containing `*.maf` or `*.maf.gz`
- `reference_fasta`: reference FASTA path
- `results_dir`: output directory

Optional controls:

- `contigs`: explicit contigs to run
- `samples`: explicit sample list
- `direct_max_missing_count`
- `direct_max_missing_fraction`
- `direct_mask_indels`
- `direct_mask_indel_adjacent_snps`
- `direct_treat_n_as_missing`
- `direct_allow_multiallelic_snps`

Indel masking behavior:

- `direct_mask_indels: true` masks reference positions directly overlapped by deletions.
- `direct_mask_indel_adjacent_snps: true` additionally masks SNPs immediately adjacent to an insertion or deletion.

## Run

Local:

```bash
snakemake -j 8
```

SLURM:

```bash
snakemake --profile profiles/slurm
```

## Outputs

Outputs are written under `results/sites/` by default:

- `combined.<contig>.all_sites.vcf`
- `combined.<contig>.variants.vcf`
- `combined.<contig>.masked.bed`
- `combined.<contig>.coverage.txt`
- `combined.<contig>.site_summary.tsv`

## Testing

```bash
pytest -q
```

## Simulation Helper

The repository includes [scripts/simulate_msprime_indels.py](/Users/jeffreyross-ibarra/src/argprep/scripts/simulate_msprime_indels.py) for generating haploid test datasets with msprime SNP variation plus branch-based indels on the tree sequence.

Example:

```bash
python scripts/simulate_msprime_indels.py \
  --sequence-length 10000 \
  --num-samples 8 \
  --theta 0.01 \
  --rho 0.01 \
  --ne 10000 \
  --indel-rate 1e-9 \
  --indel-lambda 0.01 \
  --seed 7 \
  --out-prefix /tmp/sim/example
```

Outputs:

- `<prefix>.reference.fa`
- `<prefix>.samples.fa`
- `<prefix>.indels.tsv`
- `<prefix>.summary.tsv`
- `<prefix>.maf/`

Summary fields include:

- `seed`
- `sequence_length`
- `ancestral_bp_with_indel_in_ge1_sample`
- `total_snps`
- `snps_without_overlapping_indel`
