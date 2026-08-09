# ARGprep implementation plan — 2026-08-09

## Decisions

- `combined.<contig>.all_sites.vcf` remains a required scientific output.
- Contig-name normalization remains supported because MAF and reference names may legitimately differ.
- Explicitly requested contigs should fail when any name is missing or ambiguous; automatic discovery may skip unmatched contigs with a visible warning.
- Protection for implausible inputs may be removed when doing so materially simplifies the implementation. In particular, zero-length reference contigs do not need successful output support.
- Changes should be delivered as small, independently testable commits. Correctness fixes should not be bundled with cleanup or output redesign.

## Desired end state

The workflow continues to emit:

- the complete retained-site `all_sites.vcf`;
- the variant-only VCF;
- the global mask BED;
- per-sample missing BEDs;
- per-contig site summaries;
- optional ARGweaver `.sites` files;
- `summary.html`.

The main performance change is that `summary_report.py` no longer rereads every genotype in every `all_sites.vcf`. `maf_to_sites.py` will accumulate the report statistics while each NumPy call block is already in memory and write compact per-contig statistics. The VCF remains a first-class scientific artifact, but report generation becomes proportional to summary size rather than VCF size.

## Phase 1 — fix multiallelic configuration propagation

### Implementation

1. In the `direct_maf_sites` rule, pass `--allow-multiallelic-snps` when `allow_multiallelic_snps` is true.
2. Pass `--mask-multiallelic-snps` when it is false.
3. Keep the Python CLI defaults unchanged for backward compatibility with direct script users.

### Tests

- Add a workflow-level fixture containing a site with more than two ALT alleles.
- Run it with `allow_multiallelic_snps: true` and assert that the site is retained.
- Run it with `allow_multiallelic_snps: false` and assert that the site is absent from both VCFs, present in the mask, and counted under `masked_multiallelic`.

### Acceptance criteria

- Both Boolean YAML values reach distinct, tested Python behavior.
- The existing direct-CLI multiallelic tests continue to pass.
- This commit contains no unrelated cleanup.

## Phase 2 — remove verified dead state and code

### Implementation

Remove only items whose lack of production callers has been verified:

- `per_sample_missing_retained_by_contig` and its associated accumulation in `summary_report.py`;
- `read_fasta_contigs`, `extract_end`, and `extract_dp` from `scripts/common.py`;
- `apply_indel_events` and its isolated test if it is confirmed not to be a supported public helper.

Remove `maf_threads` from the documented configuration and the rule if there is no near-term plan to add actual multithreading. Otherwise force/document a value of 1. Update `example_data/options.yaml`, which currently requests two threads for a single-threaded program.

### Tests

- Run the complete suite.
- Confirm that repository-wide searches find no remaining references to removed symbols or `maf_threads`, if that option is removed.

### Acceptance criteria

- Scientific outputs are byte-identical for the example workflow.
- Resource requests no longer imply nonexistent CPU parallelism.
- Cleanup is kept separate from functional changes.

## Phase 3 — fail on unmatched explicitly requested contigs

### Implementation

Preserve current normalization rules, including case folding, leading `chr`, leading zeros, and recognized assembly prefixes.

Change only explicit `contigs:` behavior:

1. Resolve every requested name against the reference index.
2. Collect missing and ambiguous requested names.
3. Raise one actionable error listing all unresolved names if either collection is nonempty.
4. Continue deduplicating aliases that resolve to the same exact reference contig, but report the remapping.

For automatic contig discovery:

- keep the normalized intersection behavior;
- continue skipping contigs that cannot be mapped to the reference;
- preserve a clear stderr warning listing skipped names;
- fail if no contigs remain.

### Tests

- Explicit list with all exact names succeeds.
- Explicit aliases that map unambiguously succeed.
- An explicit list containing one valid and one missing name fails instead of running partially.
- An explicit ambiguous alias fails and lists its candidates.
- Automatic discovery retains its current skip-and-warn behavior.

### Acceptance criteria

- A typo in an explicit contig list cannot silently reduce the analysis.
- Existing useful MAF/reference alias handling remains available.

## Phase 4 — retain `all_sites.vcf` while eliminating its report rescan

### Statistics contract

Define a compact per-contig report-statistics file, preferably TSV. It should contain enough information to reproduce the current report without reading either VCF:

- window index or `[start, end)` coordinates;
- invariant retained-site count per window;
- variant retained-site count per window;
- masked-base count per window;
- per-sample retained called-site total;
- per-sample carried-variant total;
- sample order.

Per-sample full-contig missing-base totals may continue to come from the existing missing BEDs, or may be included in the statistics file if doing so removes another scan without duplicating policy.

### Implementation in `maf_to_sites.py`

1. Accept the configured `summary_window_bp` as an argument.
2. Allocate small window-count arrays for the current contig.
3. During the existing chunked classification loop:
   - increment invariant or variant counts for retained positions;
   - increment masked-window counts whenever a position is masked;
   - accumulate per-sample called and carried-variant totals from the current call block.
4. Write the compact report-statistics file after the normal coverage invariant succeeds.
5. Continue writing every current scientific output, including the complete `all_sites.vcf`.

The implementation should reuse the existing NumPy block and masks. It should not introduce a second per-position Python pass merely to produce statistics.

### Implementation in `summary_report.py`

1. Replace `read_all_sites_stats` with a reader for the compact statistics files.
2. Remove genotype parsing and the unused retained-missing structure.
3. Preserve the current HTML totals, window plots, per-sample tables, and sample order.
4. Continue reading mask or missing BEDs only where their interval-level information is still required. If all report values are available in the compact statistics, remove those report-only rescans as a follow-up.

### Workflow changes

- Add one statistics output per contig to `direct_maf_sites`.
- Make those files inputs to `summary_report`.
- Keep `all_sites.vcf` in `rule all` and as a declared scientific output.
- Remove `all_sites.vcf` from the report rule's inputs once the report no longer reads it. This prevents a report-only rebuild from unnecessarily coupling parsing to the large files while still retaining them as workflow outputs.

### Verification

Before switching the report reader:

1. Run the existing workflow on representative fixtures under multiple options:
   - default missingness;
   - partial missingness allowed;
   - multiallelic sites allowed and masked;
   - indel-adjacent masking on and off;
   - `add_ref` on and off;
   - quality masking enabled;
   - ARGweaver output enabled.
2. Generate reports through both the old VCF-scanning path and the new statistics path.
3. Compare all numeric tables and plotted data arrays exactly.
4. Confirm that the scientific VCFs and BEDs are byte-identical before and after the change.
5. Measure report runtime and peak memory on at least one production-scale contig.

HTML text need not remain byte-identical if internal metadata or formatting changes, but every reported numeric value must match.

### Acceptance criteria

- `all_sites.vcf` remains complete and unchanged in meaning.
- `summary_report.py` does not open or parse `all_sites.vcf`.
- Report statistics match the previous implementation exactly across the option matrix.
- Report-generation time is no longer proportional to the number of VCF genotype cells.

## Phase 5 — remove zero-length-contig output support

### Implementation

1. Delete `write_empty_contig_outputs` and its duplicated output construction.
2. After loading the reference contig, raise a clear error if its length is zero.
3. Keep 1 bp and other short contigs on the normal code path.

### Tests

- Replace successful zero-length output tests with one concise rejection test.
- Retain or add a 1 bp normal-path test only if short-contig behavior is otherwise untested.

### Acceptance criteria

- The alternate zero-length output path and ARGweaver special case are gone.
- Normal scientific inputs and outputs are unchanged.

## Phase 6 — simplify optional ARGweaver output lifecycle

### Decision required

Choose one contract:

1. **Always generate `.sites`:** simplest DAG and no stale-output problem, at the cost of an output some users do not need.
2. **Generate only when enabled and do not clean stale files:** simplest optional behavior; document that changing output modes in an existing results directory may leave obsolete files.
3. **Keep current cleanup behavior:** retain only if toggling the option in-place is an important supported workflow.

Preferred simplification: option 2. Output cleanup should not be a side effect of building `summary.html`.

### Implementation if option 2 is accepted

- Remove stale `.sites` deletion from the summary rule.
- Remove the no-op ARGweaver dependency from the report rule.
- Keep conditional `.sites` generation in `direct_maf_sites` and `rule all`.
- Document the stale-file behavior.

### Tests

- Test enabled and disabled fresh runs.
- Remove the test that requires disabling the option to delete an old file.

## Phase 7 — evaluate standalone helpers

Decide based on actual use rather than code aesthetics:

- Keep `window_to_fasta.py` if users construct reference-anchored windows from `all_sites.vcf`; retaining the VCF means its architecture remains valid.
- Delete `chr1_variable_plot.py` if the main report can provide an appropriate variable-site scale. Parsing generated SVG coordinates back into data is brittle and duplicates chart rendering.

If retained, add focused tests only for supported behavior. Do not expand edge-case handling for these standalone helpers.

## Phase 8 — tighten parsing policy last

This phase has the greatest chance of rejecting files that currently limp through, so it should follow the higher-value changes.

### Recommended strictness

- Propagate required-file `OSError` exceptions instead of converting them into empty contig collections.
- For selected MAF `s` rows, require both alignment strings to have equal lengths before iterating them.
- Raise on malformed numeric MAF fields with file/block context.
- Decide explicitly whether malformed quality-BED rows should fail or be skipped. Favor failure for non-comment data rows because silent masking omissions have scientific consequences.

Do not add elaborate recovery. The goal is a short fail-fast path with actionable errors.

### Tests

- One malformed MAF row test.
- One unequal alignment-string test.
- One malformed quality-BED row test if strict BED parsing is selected.
- One missing/unreadable required-file test.

## Commit sequence

Recommended commits:

1. `Fix multiallelic masking config propagation`
2. `Remove unused summary and common helpers`
3. `Fail on unresolved explicit contigs`
4. `Emit compact report statistics during site calling`
5. `Build summary report from compact statistics`
6. `Reject zero-length reference contigs`
7. `Simplify optional ARGweaver output lifecycle`
8. `Remove unused standalone helper scripts` (only after usage decisions)
9. `Fail fast on malformed required inputs`

Phases 4 and 5 of this commit list deliberately separate statistics production from report consumption. During development, the old report reader can remain available for exact cross-validation before it is removed.

## Completion checklist

- [ ] `allow_multiallelic_snps: false` is effective through Snakemake.
- [ ] Every Boolean workflow option has at least one integration test for its non-default state.
- [ ] Explicit unresolved contigs fail; automatic discovery remains compatible.
- [ ] `all_sites.vcf` remains a required, documented scientific output.
- [ ] All preexisting scientific outputs are unchanged except where a confirmed bug is fixed.
- [ ] Report generation no longer parses `all_sites.vcf`.
- [ ] Report values match the old implementation exactly.
- [ ] Zero-length contigs are rejected instead of receiving duplicate output handling.
- [ ] Single-threaded resource documentation matches actual behavior.
- [ ] Optional ARGweaver output has one clear lifecycle contract.
- [ ] Full tests pass in the `argprep` conda environment.
- [ ] Production-scale report runtime and memory measurements are recorded.
