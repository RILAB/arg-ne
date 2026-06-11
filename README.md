# ARGprep Pipeline 

This repository provides a Snakemake workflow for processing AnchorWave MAFs directly into per-contig site outputs. The workflow emits all-sites VCFs, variant-only VCFs, and BED masks from the alignments. Written with the aid of [Codex](https://openai.com/codex/) and [Claude](https://claude.ai/). Note that v1.0 was a major rewrite from v0.4, and no longer uses Tassel or GATK. 

If you use this please cite: 

Ross-Ibarra, J. 2026. ARGprep: A pipeline to prepare pairwise whole-genome alignments for ancestral recombination graph estimation. [doi: 10.5281/zenodo.19655050](https://doi.org/10.5281/zenodo.19655050)

If your use case is pairwise variant discovery (SNPs, large indels, inversions) rather than ARG-ready all-sites output, [wgatools](https://github.com/wjwei-handsome/wgatools) is a potential alternative. See [WGATOOLS_COMPARISON.md](WGATOOLS_COMPARISON.md) for a detailed comparison of the two approaches.

> **Which version to use:** check out the most recent tagged commit (e.g. `git checkout v1.5`) rather than an older release. See [changelog.md](changelog.md) for a per-version breakdown of changes.

## Requirements

- Conda
- The environment defined in `argprep.yml`

## Setup

```bash
conda env create -f argprep.yml
conda activate argprep
```

## Container Setup (alternative)

Instead of a local Conda environment, you can use Singularity or Docker. Both
containers bundle all dependencies, the pipeline code, example data, and the
full test suite. Tests are run automatically during the build so the image is
verified before it is produced.

### Singularity

```bash
singularity build --fakeroot argprep.sif singularity.def
```

Run the pipeline on your data (from your project directory):

```bash
singularity exec argprep.sif snakemake -j 8 \
  --snakefile /opt/argprep/Snakefile \
  --configfile options.yaml
```

Run the bundled example:

```bash
singularity exec argprep.sif snakemake -j 4 \
  --snakefile /opt/argprep/Snakefile \
  --configfile /opt/argprep/example_data/options.yaml \
  --config results_dir=example_results
```

Run tests inside the container:

```bash
singularity exec argprep.sif bash -c 'cd /opt/argprep && pytest -q -p no:cacheprovider'
```

Our singularity setup defaults to `/tmp` but for HPC users that want to use some variable made available by their HPC admins like `$SCRATCH` then you could re-write your ENV variable:

```bash
singularity exec --env XDG_CACHE_HOME=$SCRATCH/argprep_cache argprep.sif
```

### Docker

```bash
docker build -t argprep .
```

Run the pipeline (mount your data directory):

```bash
docker run --rm -v $(pwd):/data -w /data argprep \
  -j 8 --snakefile /opt/argprep/Snakefile --configfile options.yaml
```

Run tests:

```bash
docker run --rm --entrypoint pytest argprep -q -p no:cacheprovider
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
- `mask_indel_adjacent_snps: false` - when `true`, mask SNPs immediately flanking an indel in any sample (see [NOTES.md](NOTES.md) for exact semantics)
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
- A sample counts as missing at a site if it has no alignment block covering that position, carries a gap (`-`), an `N`, or any other non-ACGT character. To drop every site overlapped by a deletion in any sample, set `max_missing_count: 0` (the default) — gaps always contribute to the missing-sample count.

Adjacent-SNP masking is documented in detail in [NOTES.md](NOTES.md).

Reference-sample behavior:

- `add_ref: true` appends a synthetic `REF` sample (genotype `0` at every retained site) to both final VCFs.

## Run

Local:

```bash
snakemake -j 8 --configfile options.yaml
```

SLURM (recommended — submit the controller as its own job):

```bash
sbatch profiles/slurm/run-controller.sbatch options.yaml
# watch progress:
tail -f logs/slurm/controller-<jobid>.out
```

This runs the long-lived Snakemake controller as a SLURM job on a
**non-preemptable** partition and lets it submit the per-rule jobs. The
controller must outlive every rule job, so it should never run on a preemptable
queue. The wrapper always passes `--rerun-incomplete`.

The two cluster-specific values in `run-controller.sbatch` are
`--partition=high` and `--account=jrigrp` (the UCD farm defaults). Override them
for another cluster without editing the file — `sbatch` CLI flags win over the
`#SBATCH` lines:

```bash
sbatch --partition=<your-nonpreemptable> --account=<your-acct> \
       profiles/slurm/run-controller.sbatch options.yaml
```

#### Running rule jobs on a preemptable queue

Set `slurm_partition: low` (or your cluster's preemptable partition) in your
config file to send the per-rule jobs to the cheap queue. Preemption is handled
safely: `profiles/slurm/config.yaml` wires in `profiles/slurm/status-sacct.sh`,
which maps a `PREEMPTED` job to *running* rather than *failed*. A preempted job
is auto-requeued by SLURM (requires `PreemptMode=REQUEUE`, the farm `low`
default) and reruns the rule from scratch; the controller waits for it instead
of aborting the run. The same status command also makes genuinely failed jobs
(TIMEOUT, OOM, NODE_FAIL, scancel) fail cleanly instead of hanging.

You can still launch the controller directly (e.g. on the head node for a quick
run), but it is then vulnerable to being killed:

```bash
snakemake --profile profiles/slurm --configfile options.yaml --rerun-incomplete
```

When using the SLURM profile, set `slurm_account` and `slurm_partition` in your config file. Slurm defaults for other resources are defined in `profiles/slurm/config.yaml`. Parsing the MAFs is the most computationally expensive step in the pipeline, and direct-maf rule resources can be overridden in `options.yaml` (`maf_threads`, `maf_mem_mb`, `maf_time`).

### Try it on the bundled example

`example_data/` ships with a small simulated dataset (`example.maf/`, `example.reference.fa`) and a matching `options.yaml`. From the repo root:

```bash
snakemake -j 4 --configfile example_data/options.yaml
```

Outputs land in `example_results/`. To regenerate the example data from scratch, see the [Simulation Helper](#simulation-helper) section.

## Outputs

Outputs are written under `results/` by default (or under `results_dir` if provided):

- `sites/combined.<contig>.all_sites.vcf` — every retained site (invariant + variant) that passed all filters; `INFO=SC=invariant|variant` distinguishes the two
- `sites/combined.<contig>.vcf` — variant-only subset of `all_sites.vcf`
- `sites/combined.<contig>.mask.bed` — merged BED intervals for masked positions
- `sites/combined.<contig>.site_summary.tsv` — per-contig counts (see table below)
- `sites/combined.<contig>.<sample>.missing.bed` — per-sample missing regions used by per-MAF summary stats; 4-column BED (`chrom`, `start`, `end`, `sample`)
- `summary.html` — genome-wide overview plus per-MAF tables and per-contig per-MAF breakdowns

Both VCFs share the same header and use a single haploid `GT` per sample (`0` for the REF allele, `1`/`2`/... for ALTs in `ALT` order, `.` for missing). `INFO` carries `NS` (non-missing samples), `MS` (missing samples), and `SC` (`invariant` or `variant`). All retained sites are emitted with `FILTER=PASS`; filtered-out positions appear in the BED mask, not the VCFs.

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
| `masked_indel_adjacent` | SNPs masked because they immediately flank an indel |
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
