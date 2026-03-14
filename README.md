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

- `max_missing_count`
- `max_missing_fraction`
- `mask_indels`
- `mask_indel_adjacent_snps`
- `treat_n_as_missing`
- `allow_multiallelic_snps`

Advanced override:

- `contigs`: restrict the run to specific contigs instead of using the shared MAF/reference contigs automatically
- `samples`: restrict the run to specific sample basenames instead of using all `*.maf` / `*.maf.gz` files in `maf_dir`

Example CLI override:

```bash
snakemake -j 8 --config samples='["sampleA","sampleB"]' contigs='["chr1","chr2"]'
```

Missingness thresholds:

- `max_missing_count` is an absolute number of missing samples allowed at a retained site.
- `max_missing_fraction` is a fraction of samples allowed to be missing.
- If both are set, the workflow uses the stricter threshold.
- The fraction is converted to a count with downward truncation. For example, with 10 samples, `0.15` allows `1` missing sample.

Indel masking behavior:

- `mask_indels: true` masks reference positions directly overlapped by deletions.
- `mask_indel_adjacent_snps: true` additionally masks SNPs immediately adjacent to an insertion or deletion.
- `mask_indels: false` disables indel-based masking entirely, so indel-overlapped and indel-adjacent sites are judged only by the remaining filters such as missingness.
- `mask_indel_adjacent_snps` only has an effect when `mask_indels: true`.

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
- `combined.<contig>.site_summary.tsv`
- `summary.html`

The pipeline still validates that retained sites plus the mask span each contig exactly, but that check is now internal and is no longer written as a separate `coverage.txt` file.

## Testing

```bash
pytest -q
```

## Simulation Helper

The repository includes [scripts/simulate_msprime_indels.py](https://github.com/RILAB/argprep/blob/main/scripts/simulate_msprime_indels.py) for generating haploid test datasets with msprime SNP variation plus branch-based indels on the tree sequence.

Example:

```bash
python scripts/simulate_msprime_indels.py \
  --sequence-length 1000000 \  # ancestral sequence length in bp
  --num-samples 8 \  # number of haploid samples
  --theta 0.01 \  # scaled mutation parameter, 4Ne*mu
  --rho 0.01 \  # scaled recombination parameter, 4Ne*r
  --ne 10000 \  # effective population size used to convert theta and rho
  --indel-rate 1e-8 \  # indel events per bp per generation on branches
  --indel-lambda 0.001 \  # exponential indel size rate; mean size = 1/lambda = 1000 bp
  --seed 8675309 \  # RNG seed for reproducibility
  --out-prefix example_mafs/example  # output prefix for FASTA, MAF, TSV files
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
