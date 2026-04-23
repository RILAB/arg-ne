# ARGprep Pipeline 

This repository provides a Snakemake workflow for processing AnchorWave MAFs directly into per-contig site outputs. The workflow emits all-sites VCFs, variant-only VCFs, and BED masks from the alignments. Written with the aid of [Codex](https://openai.com/codex/) and [Claude](https://claude.ai/). Note that v1.0 was a major rewrite from v0.4, and no longer uses Tassel or GATK.

If you use this please cite: 

Ross-Ibara, J. 2026. ARGprep: A pipeline to prepare pairwise whole-genome alignments for ancestral recombination graph estimation. [![DOI](https://zenodo.org/badge/1131411242.svg)](https://doi.org/10.5281/zenodo.19655050)


## Requirements

- Conda
- The environment defined in `argprep.yml`

## Setup

```bash
conda env create -f argprep.yml
conda activate argprep
```

## Configure

Create or edit a config file such as `options.yaml`.
The workflow requires `--configfile` and will fail if it is omitted.

Required keys:

- `maf_dir`: directory containing `*.maf` or `*.maf.gz`
- `reference_fasta`: reference FASTA path

Optional path keys:

- `results_dir`: output directory (default: `results`)

Optional controls (defaults shown):

- `max_missing_count` - no default; see missingness thresholds below
- `max_missing_fraction` - no default; see missingness thresholds below
- `mask_indels: true` - mask reference positions overlapped by deletions
- `mask_indel_adjacent_snps: true` - mask SNPs immediately flanking an indel (only applies when `mask_indels: true`)
- `treat_n_as_missing: false` - treat `N` bases as missing rather than as a call
- `allow_multiallelic_snps: true` - retain sites with more than two alleles
- `add_ref: false` - append a synthetic `REF` sample (genotype `0`) to both VCFs
- `summary_window_bp: 100000` - window size in bp for binned per-contig plots in `summary.html` (this does not affect the per-MAF tables)

SLURM resource overrides (for the `direct_maf_sites` rule):

- `maf_threads: 2`
- `maf_mem_mb: 48000`
- `maf_time: "24:00:00"`

SLURM profile keys (required when using `--profile profiles/slurm`):

- `slurm_account`
- `slurm_partition`

Advanced override:

- `contigs`: restrict the run to specific contigs instead of using the automatic shared-contig behavior
- `samples`: restrict the run to specific sample basenames instead of using all `*.maf` / `*.maf.gz` files in `maf_dir`

Contig and sample selection behavior:

- If `samples` is omitted, samples are auto-discovered from both `*.maf` and `*.maf.gz` in `maf_dir`.
- If both `<sample>.maf` and `<sample>.maf.gz` exist, `<sample>.maf` is used.
- If `contigs` is omitted, the workflow uses the intersection of contigs present in all selected MAFs.
- Requested contigs are matched to reference `.fai` contigs with normalization (for example `chr01` can map to `1` when unambiguous).
- Requested contigs that are unmatched or ambiguous after normalization are skipped.
- The workflow errors only if no contigs remain after resolution.

Example CLI override:

```bash
snakemake -j 8 --configfile options.yaml --config samples='["sampleA","sampleB"]' contigs='["chr1","chr2"]'
```

Missingness thresholds:

- `max_missing_count` is an absolute number of missing samples allowed at a retained site.
- `max_missing_fraction` is a fraction of samples allowed to be missing.
- If both are set, the workflow uses the stricter threshold.
- The fraction is converted to a count with downward truncation. For example, with 10 samples, `0.15` allows `1` missing sample.
- **If neither is set, the default is 0 - any site where even one sample is unaligned or missing is masked.** Set one of these options explicitly if you want to retain sites with partial coverage.
- A sample counts as missing at a site if it has no alignment block covering that position, carries a gap/indel character (`-`), carries an `N` base, or has any other non-ACGT character. In particular, **indel gap characters always count as missing** regardless of the `mask_indels` setting — `mask_indels` controls whether indel-overlapped sites are masked outright, but even when `mask_indels: false`, a `-` at a site still contributes to the missing-sample count.

Indel masking behavior:

- `mask_indels: true` masks reference positions directly overlapped by deletions.
- `mask_indel_adjacent_snps: true` additionally masks SNPs immediately adjacent to an insertion or deletion.
- `mask_indels: false` disables indel-based masking entirely, so indel-overlapped and indel-adjacent sites are judged only by the remaining filters such as missingness.
- `mask_indel_adjacent_snps` only has an effect when `mask_indels: true`.

Reference-sample behavior:

- `add_ref: true` appends a synthetic `REF` sample to both final VCFs.
- The added sample is emitted as genotype `0` at every retained site in `all_sites` and `variants`.

## Run

Local:

```bash
snakemake -j 8 --configfile options.yaml
```

SLURM:

```bash
snakemake --profile profiles/slurm --configfile options.yaml
```

When using the SLURM profile, set `slurm_account` and `slurm_partition` in your config file. Slurm defaults for other resources are defined in `profiles/slurm/config.yaml`. Parsing the MAFs is the most computationally expensive step in the pipeline, and direct-maf rule resources can be overridden in `options.yaml` (`maf_threads`, `maf_mem_mb`, `maf_time`).

## Outputs

Outputs are written under `results/` by default (or under `results_dir` if provided):

- `sites/combined.<contig>.all_sites.vcf`
- `sites/combined.<contig>.vcf`
- `sites/combined.<contig>.mask.bed`
- `sites/combined.<contig>.site_summary.tsv`
- `sites/combined.<contig>.<sample>.missing.bed` (per-sample missing regions used by per-MAF summary stats)
- `summary.html` (genome-wide overview plus per-MAF tables and per-contig per-MAF breakdowns)

The `site_summary.tsv` contains one metric per row with columns `metric` and `value`:

| metric | description |
|---|---|
| `contig` | contig name |
| `contig_length` | contig length in bp |
| `samples` | number of samples |
| `allowed_missing` | effective missing-sample threshold used |
| `all_sites` | retained sites (invariant + variant) |
| `variants` | retained variant sites |
| `invariant` | retained invariant sites |
| `masked_total` | total masked positions |
| `masked_intervals` | number of merged BED intervals in the mask |
| `masked_missingness` | positions masked due to too many missing samples |
| `masked_indel` | positions masked due to indel overlap or adjacency |
| `masked_multiallelic` | positions masked due to more than two alleles |
| `masked_no_alignment` | positions masked because at least one sample had no alignment |
| `masked_ref_non_acgt` | reference positions with non-ACGT bases (always masked) |

The pipeline still validates that retained sites plus the mask span each contig exactly, but that check is now internal and is no longer written as a separate `coverage.txt` file.

## Testing

```bash
pytest -q
```

## Simulation Helper

The repository includes [scripts/simulate_msprime_indels.py](https://github.com/RILAB/argprep/blob/main/scripts/simulate_msprime_indels.py) for generating haploid test datasets with msprime SNP variation plus branch-based indels on the tree sequence. Note that these simulations are not intended to be evolutionarily accurate, but simply to give a reasonable example data.

Example:

```bash
python scripts/simulate_msprime_indels.py \
  --sequence-length 1000000 \
  --num-samples 8 \
  --theta 0.01 \
  --rho 0.01 \
  --ne 10000 \
  --indel-rate 1e-8 \
  --indel-lambda 0.001 \
  --seed 8675309 \
  --out-prefix example_data/example
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
- `reference_bp_with_indel_in_ge1_sample`
- `total_snps`
- `snps_without_overlapping_indel`
