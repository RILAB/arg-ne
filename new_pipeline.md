# New Pipeline Design

## Goal

Replace the current TASSEL + GATK-heavy workflow with a direct MAF-based pipeline that produces three final outputs:

1. A VCF containing both variant and invariant sites that satisfy a maximum missingness threshold.
2. A VCF containing only variant sites that satisfy the same missingness threshold.
3. A BED file marking regions where alignment is absent or missingness exceeds the threshold.

This is a better match to the biological intent than reconstructing gVCF semantics from MAFs and then filtering back down to alignment-aware site masks.


## Why Start Directly From MAF

The MAF files already contain the information needed to answer the real questions:

- Which reference positions are aligned in each sample?
- Which positions are missing in one or more samples?
- Which positions are invariant?
- Which positions are polymorphic?
- Which positions should be masked because too many samples are missing?

The current workflow uses multiple format conversions:

`MAF -> per-sample gVCF -> cleaned gVCF -> merged gVCF -> split/filter/mask`

That adds operational complexity, makes the logic harder to validate, and obscures the actual rule set being applied to alignment-derived positions.

A direct pipeline would instead look like:

`MAFs -> per-position multisample alignment table -> VCFs + BED mask`


## Proposed Outputs

### 1. All-sites VCF

Suggested name:

`all_sites.maxmissX.vcf.gz`

Contents:

- Includes both invariant and variant positions.
- Includes only positions where missingness is at or below the threshold.
- Uses the reference base as `REF`.
- Uses `ALT=.` for invariant sites.
- Uses standard ALT alleles for variant sites.


### 2. Variant-only VCF

Suggested name:

`variants.maxmissX.vcf.gz`

Contents:

- Subset of the all-sites VCF.
- Includes only positions that are polymorphic among retained calls.
- Uses the same missingness threshold as the all-sites VCF.


### 3. BED Mask

Suggested name:

`masked.maxmissX.bed.gz`

Contents:

- Merged intervals where missingness exceeds the threshold.
- Includes regions with no alignment.
- May also include sites removed for ambiguity or indel-related rules, depending on final decisions.


## Core Processing Model

For each reference contig:

1. Parse MAF blocks for all samples.
2. Project aligned bases onto reference coordinates.
3. Build a per-position table across samples.
4. For each reference position:
   - determine which samples have usable calls
   - count missing samples
   - determine whether the site is invariant or variant
   - decide whether the position is kept or masked
5. Merge adjacent masked bases into BED intervals.
6. Emit:
   - all-sites VCF
   - variant-only VCF
   - BED mask


## Recommended Site Logic

At each reference position:

- If missingness exceeds the threshold, the site goes to the BED mask and is excluded from both VCFs.
- If missingness is within threshold:
  - emit to the all-sites VCF
  - emit to the variant-only VCF only if polymorphic

This gives one consistent rule across all outputs.


## Key Decisions That Must Be Made

These need to be fixed before implementation.

### Missingness definition

Decide what counts as missing:

- no aligned MAF block for the sample at that reference position
- gap character `-`
- ambiguous base such as `N`
- other non-ACGT symbols
- masked/lowercase sequence, if present

Recommended default:

- treat no block, `-`, and non-ACGT bases as missing


### Threshold semantics

Decide whether missingness is:

- an absolute maximum number of missing samples
- a fraction of samples

Recommended default:

- support both
- internal logic should reduce both to an absolute missing-count threshold per run


### Variant definition

Decide whether to keep:

- only biallelic SNPs
- multiallelic SNPs
- non-SNP substitutions

Recommended default:

- keep SNPs only
- allow multiallelic SNPs if all ALT alleles are single-base A/C/G/T


### Indel treatment

Decide how to treat positions involved in insertions/deletions.

Options:

- exclude all indel-affected columns from both VCFs and place them in BED
- try to represent indels explicitly in the variant VCF

Recommended default:

- exclude indel-affected columns from both VCFs and mask them in BED

This keeps the outputs aligned with the goal of producing clean SNP-only site sets.


### Invariant-site definition

Decide what counts as invariant for the all-sites VCF:

- only positions where all retained non-missing calls match the reference
- positions where all retained calls are identical, even if they differ from the reference

Recommended default:

- only emit reference-consistent invariant sites


### Reference coordinate authority

Decide whether the reference FASTA or the reference sequence embedded in the MAF blocks is the source of truth.

Recommended default:

- reference FASTA should define contig names, lengths, and REF alleles
- MAF-derived reference sequence should be treated as supporting alignment data, not the final authority


### Sample set handling

Decide how sample membership is determined:

- explicit sample list
- inferred from files
- inferred from MAF block labels

Recommended default:

- allow explicit sample selection
- fail clearly if required samples are absent


## Suggested Implementation Shape

The simplest implementation is likely one main Python script:

`scripts/maf_to_sites.py`

Inputs:

- `--maf-dir`
- `--reference-fasta`
- `--samples`
- `--max-missing-count` and/or `--max-missing-fraction`
- optional flags such as:
  - `--mask-indels`
  - `--treat-n-as-missing`
  - `--allow-multiallelic-snps`

Outputs:

- `all_sites.vcf.gz`
- `variants.vcf.gz`
- `masked.bed.gz`

Optional helper outputs:

- site-level QC summary
- contig-level counts
- positions excluded for each reason


## Internal Data Model

One practical approach:

- parse each sample’s MAF into reference-coordinate intervals
- stream contig by contig
- maintain per-position calls across samples

For each reference position, store:

- contig
- 1-based reference position
- reference base
- per-sample call
- missing sample count
- site class:
  - masked
  - invariant
  - variant


## Validation Strategy

The new pipeline should be validated against small synthetic MAFs first.

Tests should cover:

- fully invariant aligned sites
- simple SNPs
- multiallelic SNPs
- positions with one or more missing samples
- threshold edge cases
- no-alignment regions
- indel-containing regions
- ambiguous bases
- multiple contigs

The most important property is that the BED mask and the VCF outputs are logically consistent:

- any site in either VCF must not be masked
- any masked position must be absent from both VCFs


## Migration Recommendation

Do not rewrite everything at once.

Recommended path:

1. Implement a standalone prototype that reads small MAF sets directly.
2. Lock down missingness and SNP rules with tests.
3. Compare output against the current workflow on example data.
4. Replace the current multi-stage path once the direct method is trusted.


## Summary

The desired outputs are fundamentally alignment-derived masks and SNP site tables. Because of that, a direct MAF-based implementation is likely:

- simpler
- easier to explain
- easier to validate
- less dependent on format-conversion machinery

The main remaining work is not tool selection. It is making the biological rules explicit and implementing them consistently.
