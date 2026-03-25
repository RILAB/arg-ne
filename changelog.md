## Changes since v1.0

- Added optional `add_ref` output behavior so final `all_sites` and `variants` VCFs can include a synthetic `REF` sample with genotype `0` at every retained site.
- Fixed .gitignore not updating example_data

## New Pipeline Argrep v1.0

- Replaced the legacy TASSEL/gVCF/GATK workflow with a direct MAF-to-sites pipeline.
- The workflow now reads per-sample MAFs directly and emits three per-contig outputs:
  - `combined.<contig>.all_sites.vcf`
  - `combined.<contig>.variants.vcf`
  - `combined.<contig>.masked.bed`
- Site classification is now done in reference coordinates from the MAF alignment itself rather than from merged gVCF records.
- Missingness filtering is configurable with either `max_missing_count` or `max_missing_fraction`; if both are set, the stricter threshold is used.
- Indel handling is configurable:
  - `mask_indels` masks deleted reference positions
  - `mask_indel_adjacent_snps` optionally masks SNPs adjacent to insertions or deletions
- Multiallelic SNPs are retained by default and can be disabled with `allow_multiallelic_snps`.
- Coverage validation was simplified: `maf_to_sites.py` now performs the retained-sites-plus-mask contig-span check internally without a separate `coverage.txt` artifact.
- Added a direct-pipeline `summary.html` report with 100 kb window plots showing percent invariant, percent variable, and percent missing along each contig.
- Added an msprime-based simulation helper that generates:
  - reference FASTA
  - sample FASTA sequences
  - pairwise MAFs against the reference sequence
  - realized indel tables
  - summary truth counts for SNPs and indel-affected reference positions

## Changes since v0.4

- Clean VCF handling was tightened: `<NON_REF>` is now always stripped from `ALT`, and reference-donor `AD/PL/DP` values are copied where needed for cleaner downstream outputs.
- Summary warning collection was fixed to use only current-run logs passed by Snakemake, avoiding leakage of historical warnings from old `.snakemake` and SLURM logs.
- Pipeline defaults were updated in `config.yaml`/`Snakefile`, including example-data-oriented path defaults and associated resource/default cleanup.
- Integration tests were adjusted to be more portable (not relying on hardcoded `PATH` behavior), and test coverage was updated across split/integration/summary-related tests.
- Random-position integration test paths were aligned with example-output conventions and related test defaults were refreshed.
- README/changelog documentation was refreshed with expanded SLURM usage notes, formatting fixes, and changelog updates.

## Changes since v0.3

- Added/expanded contig handling across the workflow: reference renaming from gVCF/MAF contigs, numeric-suffix normalization, configured-contig remapping (e.g. `chr1` -> `1`), filtering contigs to the reference `.fai`, and defaulting to the shared MAF-contig intersection when `contigs:` is not set.
- Improved missing-contig behavior and reporting: split/merge steps now skip missing per-sample contigs more safely, preserve valid placeholders, and surface clearer warnings in the HTML summary.
- Expanded HTML summary reporting with stronger warning parsing, MAF-vs-reference contig mismatch warnings, ploidy reporting, and a missing-genotype exclusion histogram.
- Fixed summary warning de-duplication so contigs reported in "Configured contigs were remapped..." are not also repeated in generic MAF/reference mismatch warnings.
- Improved split/filter/mask robustness: refined split classification and filtered BED span handling, gzip/bgzip split-output support, gzipped missing-BED input support, and better empty-coverage handling in `check_split_coverage.py`.
- `.clean.vcf` output now always strips `<NON_REF>` from `ALT`, including GT-only (`outputJustGT=true`) runs.
- Added optional `add_reference` split output behavior: when enabled, `.clean.vcf` appends a synthetic `REF` sample with reference-only genotypes (`0`, `0/0`, etc.) and fallback ploidy matching the pipeline-resolved ploidy.
- Fixed boolean config parsing in `Snakefile` so CLI overrides like `--config bgzip_output=false add_reference=true` are interpreted correctly instead of treating non-empty strings as truthy.
- Added accessibility mask outputs (`combined.<contig>.accessible.npz`) for scikit-allel workflows and documented downstream use.
- Added ploidy inference from MAF blocks (with summary reporting) while retaining config override support.
- Added and refined resource/concurrency controls: configurable `maf_to_gvcf_*` and `merge_contig_*` resources, `merge_gvcf_max_jobs`, `maf_to_gvcf_array_max_jobs`, GenomicsDB buffer tuning, and `default_mem_mb` fallback behavior.
- Improved SLURM profile/runtime behavior: numeric default resource values in the profile, SLURM stdout/stderr written to `logs/slurm/`, and clearer documentation of profile `jobs` (Snakemake parallelism) vs per-job CPUs.
- Added substantial test coverage: pytest scaffolding, split/filt/dropSV unit tests, coverage checks, integration tests (including contig mismatch/default-contig behavior), random-position validation, summary-report regressions, Snakefile contig-resolution unit tests, and a SLURM-profile dry-run smoke test.
- Documentation updates: `logic.md`, expanded README notes/options/outputs/testing instructions, and a standalone `changelog.md`.

## Changes since v0.2

- Moved HTML summary generation into `scripts/summary_report.py` and simplified `Snakefile`.
- Corrected example MAF inputs so `example_data/*.maf.gz` are valid gzip files.
- Updated split classification so `ALT=<NON_REF>`-only records are treated as invariant.
- Updated SINGER clean-output formatting to strip `<NON_REF>` while preserving genotype/sample fields.
- Added `logic.md` with detailed site-routing/filtering logic and concrete examples.
- Added `results/split/combined.<contig>.coverage.txt` to documented workflow outputs.
- Added split-test coverage for invariant/nonref and genotype-preserving clean formatting.
- Updated SLURM profile default resources to numeric values to avoid resource conversion/submission errors.
- Added merge_gvcf_max_jobs pipeline concurrency control

## Changes since v0.1

- Added HTML summary report with embedded SVG histograms and expanded output details.
- Split logic tightened: clean sites now require all samples called; missing GTs are routed to filtered.
- Invariant/filtered/clean outputs are enforced as mutually exclusive per position; filtered BED spans now respect END/REF lengths and subtract inv/clean.
- Merged gVCFs are produced via GATK SelectVariants with genotype calling; TASSEL `outputJustGT` default set to `false` to retain likelihoods for calling.
- Added accessibility mask generation (`combined.<contig>.accessible.npz`) for scikit-allel workflows.
- New/expanded validation and tests: split coverage checks, filtered-bed tests, integration tests gated by `RUN_INTEGRATION=1`.
- Example data regenerated via msprime with indels and missing data and AnchorWave-style MAF formatting.
- `check_split_coverage.py` now reports overlap intervals with file names to aid debugging.
- `filt_to_bed.py` filters masks to the target contig, preventing cross-contig lines in `combined.<contig>.filtered.bed`.
- SLURM default resources now read `default_*` from `config.yaml` instead of hardcoded profile values.
