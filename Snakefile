import sys
from pathlib import Path

from snakemake.io import glob_wildcards

from scripts.common import normalize_contig, read_maf_contigs

if not workflow.configfiles:
    raise ValueError(
        "A config file is required. Run Snakemake with --configfile path/to/options.yaml."
    )

wildcard_constraints:
    contig="[^/]+"


def _config_bool(value, default=False):
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        return value != 0
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"1", "true", "t", "yes", "y", "on"}:
            return True
        if normalized in {"0", "false", "f", "no", "n", "off", ""}:
            return False
    raise ValueError(f"Invalid boolean config value: {value!r}")


REQUIRED_CONFIG_KEYS = ("maf_dir", "reference_fasta")
missing_config_keys = [key for key in REQUIRED_CONFIG_KEYS if key not in config]
if missing_config_keys:
    raise ValueError(
        "Missing required config keys: "
        + ", ".join(missing_config_keys)
        + ". Check the file passed via --configfile."
    )


MAF_DIR = Path(config["maf_dir"]).resolve()
ORIG_REF_FASTA = Path(config["reference_fasta"]).resolve()
RESULTS_DIR = Path(config.get("results_dir", "results")).resolve()

ALLOW_MULTIALLELIC = _config_bool(config.get("allow_multiallelic_snps", True))
MASK_INDELS = _config_bool(config.get("mask_indels", True))
MASK_INDEL_ADJACENT_SNPS = _config_bool(config.get("mask_indel_adjacent_snps", True))
TREAT_N_AS_MISSING = _config_bool(config.get("treat_n_as_missing", True))
MAX_MISSING_COUNT = config.get("max_missing_count")
MAX_MISSING_FRACTION = config.get("max_missing_fraction")

DEFAULT_MEM_MB = 48000
DEFAULT_THREADS = 2
DEFAULT_TIME = "24:00:00"
SUMMARY_WINDOW_BP = int(config.get("summary_window_bp", 100000))

DIRECT_REF_FASTA = RESULTS_DIR / "refs" / "reference_sites.fa"
REF_FAI = str(DIRECT_REF_FASTA) + ".fai"

def _maf_path_for_sample(sample: str) -> Path:
    maf = MAF_DIR / f"{sample}.maf"
    maf_gz = MAF_DIR / f"{sample}.maf.gz"
    if maf.exists():
        return maf
    if maf_gz.exists():
        return maf_gz
    return maf


def _discover_samples():
    if "samples" in config:
        return list(config["samples"])
    maf_pattern = str(MAF_DIR / "{sample}.maf")
    maf_gz_pattern = str(MAF_DIR / "{sample}.maf.gz")
    samples = set(glob_wildcards(maf_pattern).sample)
    samples.update(glob_wildcards(maf_gz_pattern).sample)
    return sorted(samples)


def _read_maf_contig_sets(samples: list[str]) -> dict[str, set[str]]:
    contigs_by_sample: dict[str, set[str]] = {}
    for sample in samples:
        contigs_by_sample[sample] = read_maf_contigs(_maf_path_for_sample(sample))
    return contigs_by_sample


def _read_fai_contigs(fai: Path) -> list[str]:
    contigs = []
    with fai.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            contigs.append(line.split("\t", 1)[0])
    return contigs


def _resolve_requested_contigs(
    requested: list[str], available: list[str]
) -> tuple[list[str], list[str], list[tuple[str, str]]]:
    available_set = set(available)
    available_norm: dict[str, list[str]] = {}
    for name in available:
        available_norm.setdefault(normalize_contig(name), []).append(name)

    kept: list[str] = []
    dropped: list[str] = []
    remapped: list[tuple[str, str]] = []
    seen: set[str] = set()
    for raw in requested:
        req = str(raw)
        mapped = req
        if req in available_set:
            mapped = req
        else:
            candidates = available_norm.get(normalize_contig(req), [])
            if len(candidates) == 1:
                mapped = candidates[0]
                remapped.append((req, mapped))
            else:
                dropped.append(req)
                continue
        if mapped not in seen:
            kept.append(mapped)
            seen.add(mapped)
    return kept, dropped, remapped


SAMPLES = _discover_samples()
if not SAMPLES:
    raise ValueError(f"No MAF files found in {MAF_DIR}")

MAF_CONTIGS_BY_SAMPLE = _read_maf_contig_sets(SAMPLES)
if MAF_CONTIGS_BY_SAMPLE:
    MAF_CONTIG_INTERSECTION = sorted(set.intersection(*MAF_CONTIGS_BY_SAMPLE.values()))
else:
    MAF_CONTIG_INTERSECTION = []


def _active_contig_resolution() -> tuple[list[str], list[str], list[str], list[tuple[str, str]]]:
    ckpt = checkpoints.index_reference.get()
    fai = Path(str(ckpt.output.fai))
    available = _read_fai_contigs(fai)
    if "contigs" in config:
        requested = [str(c) for c in config["contigs"]]
        kept, dropped, remapped = _resolve_requested_contigs(requested, available)
        if not kept:
            raise ValueError(
                "None of the configured contigs are present in reference .fai: "
                + ", ".join(requested[:10])
            )
        return kept, dropped, requested, remapped

    requested = list(MAF_CONTIG_INTERSECTION)
    if not requested:
        raise ValueError(
            "No contigs are shared across all MAF files. Set explicit 'contigs' in options.yaml to override."
        )
    kept, dropped, remapped = _resolve_requested_contigs(requested, available)
    if not kept:
        raise ValueError(
            "No shared MAF contigs are present in reference .fai. Set explicit 'contigs' in options.yaml to override."
        )
    return kept, dropped, requested, remapped


def _active_contigs() -> list[str]:
    return _active_contig_resolution()[0]


def _maf_input(sample):
    return str(_maf_path_for_sample(sample))


def _direct_prefix(contig):
    return RESULTS_DIR / "sites" / f"combined.{contig}"


def _direct_all_sites_out(contig):
    return Path(str(_direct_prefix(contig)) + ".all_sites.vcf")


def _direct_variants_out(contig):
    return Path(str(_direct_prefix(contig)) + ".variants.vcf")


def _direct_mask_out(contig):
    return Path(str(_direct_prefix(contig)) + ".masked.bed")


def _all_targets(_wc):
    contigs = _active_contigs()
    return (
        [str(_direct_all_sites_out(c)) for c in contigs]
        + [str(_direct_variants_out(c)) for c in contigs]
        + [str(_direct_mask_out(c)) for c in contigs]
        + [str(RESULTS_DIR / "summary.html")]
    )


rule all:
    input: _all_targets


rule prepare_reference:
    input:
        ref=str(ORIG_REF_FASTA),
    output:
        ref=str(DIRECT_REF_FASTA),
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{output.ref}")"
        cp "{input.ref}" "{output.ref}"
        """


checkpoint index_reference:
    input:
        ref=str(DIRECT_REF_FASTA),
    output:
        fai=REF_FAI,
    shell:
        """
        set -euo pipefail
        samtools faidx "{input.ref}"
        """


rule direct_maf_sites:
    threads: int(config.get("maf_threads", DEFAULT_THREADS))
    resources:
        mem_mb=int(config.get("maf_mem_mb", DEFAULT_MEM_MB)),
        time=str(config.get("maf_time", DEFAULT_TIME))
    input:
        mafs=lambda wc: [_maf_input(sample) for sample in SAMPLES],
        ref=str(DIRECT_REF_FASTA),
        fai=REF_FAI,
    output:
        all_sites=str(RESULTS_DIR / "sites" / "combined.{contig}.all_sites.vcf"),
        variants=str(RESULTS_DIR / "sites" / "combined.{contig}.variants.vcf"),
        mask=str(RESULTS_DIR / "sites" / "combined.{contig}.masked.bed"),
        summary=str(RESULTS_DIR / "sites" / "combined.{contig}.site_summary.tsv"),
    params:
        maf_dir=str(MAF_DIR),
        samples=" ".join(SAMPLES),
        max_missing_count=(
            None if MAX_MISSING_COUNT in (None, "") else int(MAX_MISSING_COUNT)
        ),
        max_missing_fraction=(
            None
            if MAX_MISSING_FRACTION in (None, "")
            else float(MAX_MISSING_FRACTION)
        ),
        allow_multiallelic=ALLOW_MULTIALLELIC,
        mask_indels=MASK_INDELS,
        mask_indel_adjacent_snps=MASK_INDEL_ADJACENT_SNPS,
        treat_n_as_missing=TREAT_N_AS_MISSING,
        out_prefix=lambda wc: str(_direct_prefix(wc.contig)),
    shell:
        """
        set -euo pipefail
        mkdir -p "{RESULTS_DIR}/sites"
        cmd=(python "{workflow.basedir}/scripts/maf_to_sites.py"
          --maf-dir "{params.maf_dir}"
          --reference-fasta "{input.ref}"
          --contig "{wildcards.contig}"
          --out-prefix "{params.out_prefix}"
          --samples {params.samples})
        if [ "{params.max_missing_count}" != "None" ]; then
          cmd+=(--max-missing-count "{params.max_missing_count}")
        fi
        if [ "{params.max_missing_fraction}" != "None" ]; then
          cmd+=(--max-missing-fraction "{params.max_missing_fraction}")
        fi
        if [ "{params.allow_multiallelic}" = "True" ]; then
          cmd+=(--allow-multiallelic-snps)
        fi
        if [ "{params.mask_indels}" = "True" ]; then
          cmd+=(--mask-indels)
        fi
        if [ "{params.mask_indel_adjacent_snps}" = "False" ]; then
          cmd+=(--keep-indel-adjacent-snps)
        fi
        if [ "{params.treat_n_as_missing}" = "True" ]; then
          cmd+=(--treat-n-as-missing)
        fi
        "${{cmd[@]}}"
        """


rule summary_report:
    input:
        all_sites=lambda wc: [str(_direct_all_sites_out(c)) for c in _active_contigs()],
        masks=lambda wc: [str(_direct_mask_out(c)) for c in _active_contigs()],
        summaries=lambda wc: [str(_direct_prefix(c)) + ".site_summary.tsv" for c in _active_contigs()],
        fai=REF_FAI,
    output:
        report=str(RESULTS_DIR / "summary.html"),
    params:
        window_bp=SUMMARY_WINDOW_BP,
    shell:
        """
        set -euo pipefail
        python "{workflow.basedir}/scripts/summary_report.py" \
          --fai "{input.fai}" \
          --window-bp "{params.window_bp}" \
          --report-out "{output.report}" \
          --all-sites {input.all_sites} \
          --masked-beds {input.masks} \
          --site-summaries {input.summaries}
        """
