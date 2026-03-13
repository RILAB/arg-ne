# ARG Pipeline (Snakemake)

This repository provides a Snakemake workflow for processing AnchorWave MAFs. The original path converts MAFs into per-contig merged gVCFs, splits them into clean/filtered/invariant sets, and produces mask bedfiles for downstream ARG inference. The newer direct path reads MAFs directly and emits all-sites VCFs, variant-only VCFs, and mask BEDs without going through gVCF reconstruction (see [logic.md](https://github.com/RILAB/argprep/blob/main/logic.md) for details). Written with the assistance of [Codex](https://openai.com/codex/).

## Requirements

- Conda (you may need to `module load conda` on your cluster)
- TASSEL (provided via the `tassel-5-standalone` submodule)
- GATK, Picard, htslib (installed in the conda env defined below)

## Setup

Clone with submodules so TASSEL is available:

```bash
git clone --recurse-submodules <repo>
```

Create and activate the environment (do this before running Snakemake):

```bash
module load conda
conda env create -f argprep.yml
conda activate argprep
```

## Configure

Edit `config.yaml` to point at your MAF directory and reference FASTA. At minimum you must set:

- `maf_dir`: directory containing `*.maf` or `*.maf.gz`
- `reference_fasta`: reference FASTA path (plain `.fa/.fasta` or bgzipped `.fa.gz`)
- `direct_maf_pipeline`: set to `true` to use the direct MAF-to-sites path

If your reference FASTA does **not** have an index (`.fai`), either create one (`samtools faidx`) or set `contigs:` explicitly in `config.yaml`.

### Direct MAF pipeline options

When `direct_maf_pipeline: true`, the workflow runs the direct per-contig site builder in `scripts/maf_to_sites.py`.

Key options:

- `direct_max_missing_count`: maximum number of missing samples allowed at a retained site
- `direct_max_missing_fraction`: alternate missingness threshold as a fraction of samples
- `direct_mask_indels`: if `true`, positions touched by indels are masked into the BED output
- `direct_treat_n_as_missing`: if `true`, `N` and other non-ACGT calls are treated as missing
- `direct_allow_multiallelic_snps`: if `false`, multiallelic SNP sites are masked instead of emitted

## Run 

The pipeline can be run one of two ways, both from the repo root. **It is recommended you first run on the example data provided to ensure the pipeline works on your system.**

### On Slurm 

A default SLURM profile is provided under `profiles/slurm/`. Edit `profiles/slurm/config.yaml` to customize sbatch options if needed.
Defaults for account/partition and baseline resources are set in `config.yaml` (`slurm_account`, `slurm_partition`, `default_mem_mb`, `default_time`, and `default_threads`). You will need to change these to match your own SLURM setup.
SLURM stdout/stderr logs are written to `logs/slurm/` by default.
The provided profile sets `jobs: 200` (Snakemake's max concurrent workflow jobs, equivalent to `-j 200`), which you can override on the command line if needed.

```bash
snakemake --profile profiles/slurm
```

### Locally

```bash
snakemake -j 8 
```

To run the direct MAF pipeline explicitly:

```bash
snakemake -j 8 --config direct_maf_pipeline=true
```

Common options:

- `-j <N>`: number of parallel Snakemake jobs (not per-job threads/CPUs)
- `--rerun-incomplete`: clean up partial outputs
- `--printshellcmds`: show executed commands
- `--notemp`: keep temporary intermediates (see Notes)

## Workflow Outputs

By default the workflow uses these locations (override in `config.yaml`):

### Direct MAF pipeline outputs

When `direct_maf_pipeline: true`, the main outputs are written under `results/sites/`:

- `results/sites/combined.<contig>.all_sites.vcf` : retained invariant and variant sites
- `results/sites/combined.<contig>.variants.vcf` : retained variant sites only
- `results/sites/combined.<contig>.masked.bed` : masked positions merged into BED intervals
- `results/sites/combined.<contig>.coverage.txt` : validation that all-sites VCF plus mask BED spans the full contig
- `results/sites/combined.<contig>.site_summary.tsv` : per-contig site counts by class and exclusion reason

### Legacy gVCF pipeline outputs

When `direct_maf_pipeline: false`, the workflow uses the gVCF-based path and writes:

- `gvcf/` : TASSEL gVCFs (`*.gvcf.gz`) from MAFs
- `gvcf/cleangVCF/` : cleaned gVCFs from `scripts/dropSV.py`
- `gvcf/cleangVCF/dropped_indels.bed` : bedfile of large indels
- `gvcf/cleangVCF/split_gvcf/` : per-contig gVCFs for merging
- `results/combined/combined.<contig>.gvcf.gz` : merged gVCF per contig
- `results/split/combined.<contig>.inv` : invariant sites
- `results/split/combined.<contig>.filtered` : filtered sites
- `results/split/combined.<contig>.clean.vcf.gz` : clean sites (bgzip-compressed by default)
- `results/split/combined.<contig>.missing.bed` : missing positions
- `results/split/combined.<contig>.clean.mask.bed` : merged mask bed
- `results/split/combined.<contig>.coverage.txt` : split coverage validation summary
- `results/split/combined.<contig>.accessible.npz` : boolean accessibility array (union of clean + invariant sites), for scikit-allel statistics
- `results/summary.html` : HTML summary of jobs run, outputs created, and warnings
  - Contig remaps (e.g., configured `chr1` remapped to reference contig `1`) are reported under the remap warning and are not repeated in the generic MAF-vs-reference mismatch warnings.

## Testing

Run from the repo root:

```bash
pytest -q
```

- Runs unit/fast tests by default.
- Some integration tests are gated and skipped unless enabled.

To run integration tests too:

```bash
RUN_INTEGRATION=1 pytest -q
```

- Requires the `argprep` conda environment and external tools used by the workflow (e.g., `snakemake`, `bcftools`, `tabix`, `bgzip`, `gatk`, `java`; some tests also need `samtools`).
- These tests take longer and may run parts of the Snakemake workflow.

For a quick SLURM-profile syntax/config sanity check without submitting jobs, use a dry run:

```bash
snakemake --profile profiles/slurm -n
```

## Simulation Helper

The repository includes [scripts/simulate_msprime_indels.py](/Users/jeffreyross-ibarra/src/argprep/scripts/simulate_msprime_indels.py) for generating small haploid test datasets with msprime SNP variation plus branch-based indels on the tree sequence.

Example:

```bash
python scripts/simulate_msprime_indels.py \
  --sequence-length 10000 \
  --num-samples 8 \
  --theta 0.01 \
  --rho 0.01 \
  --indel-rate 0.001 \
  --indel-lambda 0.5 \
  --seed 7 \
  --out-prefix /tmp/sim/example
```

Outputs:

- `<prefix>.reference.fa` : simulated reference sequence
- `<prefix>.samples.fa` : simulated sample sequences after SNPs and shared/private indels
- `<prefix>.indels.tsv` : per-sample indel event table with `shared_event_id` linking the same branch event across descendant samples
- `<prefix>.summary.tsv` : summary counts for ancestral bp affected by deletions in at least one sample and SNPs that do not overlap any deleted ancestral bp
- `<prefix>.maf/` : one pairwise MAF per sample against the ancestral sequence

Parameters:

- `--sequence-length` : reference length to simulate
- `--num-samples` : number of haploid samples
- `--theta` : scaled mutation parameter; the script uses `theta / 4` as the per-site mutation rate with `population_size=1`
- `--rho` : scaled recombination parameter; the script uses `rho / 4` as the per-site recombination rate with `population_size=1`
- `--indel-rate` : indel events per base per generation on tree-sequence branches; descendant samples inherit branch events according to the local genealogy
- `--indel-lambda` : rate parameter for the exponential indel size distribution; mean indel size is `1 / lambda`
- `--seed` : RNG seed for reproducibility

The simulator includes both insertions and deletions. Insertions appear in the sample FASTA/MAF outputs, while the ancestral-bp summary counts only deletion-overlapped ancestral positions because insertions do not consume reference bases.

## Notes

- You will need to change the SLURM account and partition in `config.yaml`.
- `direct_maf_pipeline: true` bypasses TASSEL/GATK merging and does not produce `results/summary.html` or the legacy split/accessibility outputs.
- The default memory setting of 20GB may be too small -- if your run is failing with OOM (out of memory), you may need to set this to a larger value.
- `bgzip_output` defaults to `true`; by default the `.inv`, `.filtered`, `.clean.vcf`, and `.missing.bed` files are bgzip-compressed (`.gz` suffix).
- All gzipped outputs in this pipeline use bgzip (required for `tabix`).
- `scripts/dropSV.py` removes indels larger than `drop_cutoff` (if set in `config.yaml`).
- `scripts/split.py` supports `--filter-multiallelic` and `--bgzip-output` (toggle via `config.yaml`).
- `scripts/filt_to_bed.py` merges `<prefix>.filtered`, `<prefix>.missing.bed`, and `dropped_indels.bed` into a final mask bed.
- `make_accessibility` builds a per-contig accessibility array from the union of `combined.<contig>.clean.vcf.gz` and `combined.<contig>.inv.gz` (default names) using the reference `.fai` to size the array. The output is a compressed NumPy archive containing a boolean array named `mask`, intended for scikit-allel statistics.
- Ploidy is inferred from MAF block structure by default (max non-reference `s` lines per block, typically `1` for pairwise MAFs). You can override with `ploidy` in `config.yaml`.
- Optional: enable `vt_normalize: true` in `config.yaml` to normalize merged gVCFs with `vt normalize` after `SelectVariants`.
- If GenomicsDBImport fails with a buffer-size error, increase `genomicsdb_vcf_buffer_size` and `genomicsdb_segment_size` in `config.yaml` (set them above your longest gVCF line length).
- Large intermediate files are marked as temporary and removed after a successful run (per-sample gVCFs, cleaned gVCFs, per-contig split gVCFs, and the GenomicsDB workspace). Use `snakemake --notemp` if you want to preserve them for debugging or reruns.
- Resource knobs (memory/threads/time) and GenomicsDB buffer sizes are configurable in `config.yaml` (e.g., `merge_contig_mem_mb`, `maf_to_gvcf_*`, `genomicsdb_*`).
- To cap concurrent contig-merge jobs on SLURM, set `merge_gvcf_max_jobs` in `config.yaml` (used by the profile as a global `merge_gvcf_jobs` resource limit).
- To cap the SLURM array concurrency for `scripts/maf_to_gvcf.sh`, set `maf_to_gvcf_array_max_jobs` in `config.yaml` (default 4).

## Downstream Uses

### ARG estimation

Use Nate Pope's SINGER Snakemake [pipeline](https://github.com/nspope/singer-snakemake) with `combined.<contig>.clean.vcf.gz` and `combined.<contig>.clean.mask.bed` as inputs.

### Population genetic statistics 

If you use scikit-allel, you can use the `combined.<contig>.clean.vcf.gz` VCF and load the accessibility mask like this:

```python
import numpy as np
mask = np.load("results/split/combined.1.accessible.npz")["mask"]
```
