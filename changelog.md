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
- Added accessibility mask generation (`combined.<contig>.accessible.npz`) for scikit‑allel workflows.
- New/expanded validation and tests: split coverage checks, filtered‑bed tests, integration tests gated by `RUN_INTEGRATION=1`.
- Example data regenerated via msprime with indels and missing data and AnchorWave‑style MAF formatting.
- `check_split_coverage.py` now reports overlap intervals with file names to aid debugging.
- `filt_to_bed.py` filters masks to the target contig, preventing cross‑contig lines in `combined.<contig>.filtered.bed`.
- SLURM default resources now read `default_*` from `config.yaml` instead of hardcoded profile values.
