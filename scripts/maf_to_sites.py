#!/usr/bin/env python3
"""
Directly project per-sample MAF alignments onto reference coordinates and emit:
  - an all-sites VCF
  - a variant-only VCF
  - a BED mask of excluded positions
  - a per-contig QC summary
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

try:
    from scripts.common import merge_intervals, normalize_contig, open_text
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from scripts.common import merge_intervals, normalize_contig, open_text


NUC_TO_CODE = {
    "A": 1,
    "C": 2,
    "G": 3,
    "T": 4,
    "-": 5,
    "N": 6,
    "?": 7,
}
CODE_TO_BASE = {
    1: "A",
    2: "C",
    3: "G",
    4: "T",
}
VALID_BASES = {"A", "C", "G", "T"}
MISSING_CODES = {0, 5, 6, 7}


@dataclass
class MafRecord:
    src: str
    start: int
    size: int
    strand: str
    src_size: int
    text: str


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--maf-dir", required=True, help="Directory containing per-sample .maf/.maf.gz files")
    ap.add_argument("--reference-fasta", required=True, help="Reference FASTA path")
    ap.add_argument("--contig", required=True, help="Reference contig to process")
    ap.add_argument("--out-prefix", required=True, help="Output prefix for this contig")
    ap.add_argument(
        "--samples",
        nargs="*",
        default=None,
        help="Explicit sample names; defaults to inferring *.maf and *.maf.gz stems from --maf-dir",
    )
    ap.add_argument("--max-missing-count", type=int, default=None)
    ap.add_argument("--max-missing-fraction", type=float, default=None)
    ap.add_argument("--allow-multiallelic-snps", dest="allow_multiallelic_snps", action="store_true", default=True)
    ap.add_argument("--mask-multiallelic-snps", dest="allow_multiallelic_snps", action="store_false")
    ap.add_argument("--mask-indels", action="store_true", default=False)
    ap.add_argument(
        "--mask-indel-adjacent-snps",
        dest="mask_indel_adjacent_snps",
        action="store_true",
        default=True,
    )
    ap.add_argument(
        "--keep-indel-adjacent-snps",
        dest="mask_indel_adjacent_snps",
        action="store_false",
    )
    ap.add_argument("--treat-n-as-missing", action="store_true", default=False)
    ap.add_argument("--add-ref", action="store_true", default=False)
    return ap.parse_args()


def discover_samples(maf_dir: Path) -> list[str]:
    samples = set()
    for path in maf_dir.iterdir():
        if path.name.endswith(".maf.gz"):
            samples.add(path.name[: -len(".maf.gz")])
        elif path.name.endswith(".maf"):
            samples.add(path.name[: -len(".maf")])
    return sorted(samples)


def maf_path_for_sample(maf_dir: Path, sample: str) -> Path:
    plain = maf_dir / f"{sample}.maf"
    gz = maf_dir / f"{sample}.maf.gz"
    if plain.exists():
        return plain
    if gz.exists():
        return gz
    raise FileNotFoundError(f"Missing MAF for sample '{sample}' under {maf_dir}")


def read_contig_sequence(reference_fasta: Path, contig: str) -> str:
    seq_parts: list[str] = []
    current: str | None = None
    with open_text(reference_fasta, "rt", errors="ignore") as handle:
        for raw in handle:
            if raw.startswith(">"):
                name = raw[1:].strip().split()[0]
                current = name
                continue
            if current == contig:
                seq_parts.append(raw.strip())
    if not seq_parts:
        raise ValueError(f"Contig '{contig}' not found in {reference_fasta}")
    return "".join(seq_parts).upper()


def iter_maf_blocks(path: Path) -> Iterator[list[MafRecord]]:
    block: list[MafRecord] = []
    with open_text(path, "rt", errors="ignore") as handle:
        for raw in handle:
            stripped = raw.strip()
            if not stripped:
                if block:
                    yield block
                    block = []
                continue
            if stripped.startswith("#"):
                continue
            parts = stripped.split()
            if not parts:
                continue
            if parts[0] == "a":
                if block:
                    yield block
                    block = []
                continue
            if parts[0] != "s" or len(parts) < 7:
                continue
            block.append(
                MafRecord(
                    src=parts[1],
                    start=int(parts[2]),
                    size=int(parts[3]),
                    strand=parts[4],
                    src_size=int(parts[5]),
                    text=parts[6],
                )
            )
    if block:
        yield block


def choose_sample_record(block: list[MafRecord], contig: str) -> tuple[MafRecord, MafRecord] | None:
    if not block:
        return None
    ref = block[0]
    if normalize_contig(ref.src) != normalize_contig(contig):
        return None
    if len(block) >= 2:
        # Some pairwise MAFs use the same contig name for both the reference
        # and query rows, so the second alignment row is still the sample.
        return ref, block[1]
    return None


def _assign_code(calls: bytearray, idx: int, code: int) -> None:
    existing = calls[idx]
    if existing == 0:
        calls[idx] = code
        return
    if existing == code:
        return
    calls[idx] = NUC_TO_CODE["?"]


def load_sample_calls(
    maf_path: Path,
    contig: str,
    contig_len: int,
    *,
    mask_indels: bool,
    treat_n_as_missing: bool,
) -> tuple[bytearray, bytearray, bytearray]:
    calls = bytearray(contig_len)
    indel_flags = bytearray(contig_len)
    adjacent_indel_flags = bytearray(contig_len)

    for block in iter_maf_blocks(maf_path):
        chosen = choose_sample_record(block, contig)
        if chosen is None:
            continue
        ref_record, sample_record = chosen
        ref_pos = ref_record.start
        prev_ref_idx: int | None = None
        mark_next_ref_adjacent = False
        for ref_char, sample_char in zip(ref_record.text.upper(), sample_record.text.upper()):
            if ref_char == "-":
                if mask_indels and sample_char != "-" and prev_ref_idx is not None:
                    adjacent_indel_flags[prev_ref_idx] = 1
                    mark_next_ref_adjacent = True
                continue

            if ref_pos >= contig_len:
                break
            idx = ref_pos
            ref_pos += 1
            prev_ref_idx = idx

            if mask_indels and mark_next_ref_adjacent:
                adjacent_indel_flags[idx] = 1
                mark_next_ref_adjacent = False

            if sample_char == "-":
                _assign_code(calls, idx, NUC_TO_CODE["-"])
                if mask_indels:
                    indel_flags[idx] = 1
                    adjacent_indel_flags[idx] = 1
                    if idx > 0:
                        adjacent_indel_flags[idx - 1] = 1
                    mark_next_ref_adjacent = True
                continue

            if sample_char in VALID_BASES:
                _assign_code(calls, idx, NUC_TO_CODE[sample_char])
                continue

            if sample_char == "N" and not treat_n_as_missing:
                _assign_code(calls, idx, NUC_TO_CODE["N"])
                continue

            _assign_code(calls, idx, NUC_TO_CODE["?"])

    return calls, indel_flags, adjacent_indel_flags


def missing_threshold(sample_count: int, max_missing_count: int | None, max_missing_fraction: float | None) -> int:
    if max_missing_count is None and max_missing_fraction is None:
        return 0
    thresholds: list[int] = []
    if max_missing_count is not None:
        thresholds.append(max(0, max_missing_count))
    if max_missing_fraction is not None:
        if not 0 <= max_missing_fraction <= 1:
            raise ValueError("--max-missing-fraction must be between 0 and 1")
        thresholds.append(int(sample_count * max_missing_fraction))
    return min(thresholds)


def format_gt(call_code: int, alt_order: list[str]) -> str:
    if call_code in MISSING_CODES:
        return "."
    base = CODE_TO_BASE[call_code]
    if base not in alt_order:
        return "0"
    return str(alt_order.index(base) + 1)


def vcf_header(contig: str, contig_len: int, samples: list[str]) -> list[str]:
    return [
        "##fileformat=VCFv4.2",
        "##source=argprep.maf_to_sites",
        f"##contig=<ID={contig},length={contig_len}>",
        '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of non-missing samples">',
        '##INFO=<ID=MS,Number=1,Type=Integer,Description="Number of missing samples">',
        '##INFO=<ID=SC,Number=1,Type=String,Description="Site class: invariant or variant">',
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]


def write_lines(path: Path, lines: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        for line in lines:
            handle.write(line)
            handle.write("\n")


def intervals_from_positions(contig: str, masked_positions: list[int]) -> list[tuple[str, int, int]]:
    intervals: list[tuple[int, int]] = []
    for pos in masked_positions:
        intervals.append((pos, pos + 1))
    merged = merge_intervals(intervals)
    return [(contig, start, end) for start, end in merged]


def summarize_site_and_mask_coverage(
    contig: str,
    contig_len: int,
    site_positions_0based: list[int],
    masked_positions_0based: list[int],
) -> dict[str, int | str]:
    site_intervals = merge_intervals([(pos, pos + 1) for pos in site_positions_0based])
    mask_intervals = merge_intervals([(pos, pos + 1) for pos in masked_positions_0based])

    i = j = 0
    overlap = 0
    while i < len(site_intervals) and j < len(mask_intervals):
        site_start, site_end = site_intervals[i]
        mask_start, mask_end = mask_intervals[j]
        if site_end <= mask_start:
            i += 1
            continue
        if mask_end <= site_start:
            j += 1
            continue
        overlap += min(site_end, mask_end) - max(site_start, mask_start)
        if site_end <= mask_end:
            i += 1
        else:
            j += 1
    if overlap:
        raise ValueError(f"overlap detected for {contig}: site_vcf vs mask_bed={overlap}")

    total = sum(end - start for start, end in merge_intervals(site_intervals + mask_intervals))
    if total != contig_len:
        raise ValueError(f"coverage mismatch for {contig}: union={total}, chrom_len={contig_len}")

    return {
        "chrom": contig,
        "site_bp": sum(end - start for start, end in site_intervals),
        "mask_bed_bp": sum(end - start for start, end in mask_intervals),
        "total_bp": total,
        "chrom_len": contig_len,
    }


def main() -> None:
    args = parse_args()
    maf_dir = Path(args.maf_dir)
    reference_fasta = Path(args.reference_fasta)
    contig = args.contig
    out_prefix = Path(args.out_prefix)
    samples = sorted(args.samples) if args.samples else discover_samples(maf_dir)
    if not samples:
        raise ValueError(f"No samples found under {maf_dir}")
    output_samples = [*samples, "REF"] if args.add_ref else samples

    contig_seq = read_contig_sequence(reference_fasta, contig)
    contig_len = len(contig_seq)
    sample_arrays: list[bytearray] = []
    indel_flags = bytearray(contig_len)
    adjacent_indel_flags = bytearray(contig_len)

    for sample in samples:
        calls, sample_indels, sample_adjacent_indels = load_sample_calls(
            maf_path_for_sample(maf_dir, sample),
            contig,
            contig_len,
            mask_indels=args.mask_indels,
            treat_n_as_missing=args.treat_n_as_missing,
        )
        sample_arrays.append(calls)
        for idx, flag in enumerate(sample_indels):
            if flag:
                indel_flags[idx] = 1
        for idx, flag in enumerate(sample_adjacent_indels):
            if flag:
                adjacent_indel_flags[idx] = 1

    allowed_missing = missing_threshold(
        len(samples),
        args.max_missing_count,
        args.max_missing_fraction,
    )

    all_sites_path = out_prefix.with_suffix(out_prefix.suffix + ".all_sites.vcf")
    variants_path = out_prefix.with_suffix(out_prefix.suffix + ".vcf")
    mask_path = out_prefix.with_suffix(out_prefix.suffix + ".mask.bed")
    summary_path = out_prefix.with_suffix(out_prefix.suffix + ".site_summary.tsv")

    all_sites_path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(all_sites_path, "wt") as all_sites, open_text(variants_path, "wt") as variants:
        for line in vcf_header(contig, contig_len, output_samples):
            all_sites.write(f"{line}\n")
            variants.write(f"{line}\n")

        masked_positions: list[int] = []
        retained_positions: list[int] = []
        counts: Counter[str] = Counter()

        for idx, ref_base in enumerate(contig_seq):
            pos = idx + 1
            if ref_base not in VALID_BASES:
                masked_positions.append(idx)
                counts["masked_ref_non_acgt"] += 1
                continue

            alleles: list[str] = []
            allele_set: set[str] = set()
            missing = 0
            call_codes: list[int] = []
            has_unaligned_sample = False
            for calls in sample_arrays:
                code = calls[idx]
                call_codes.append(code)
                if code == 0:
                    has_unaligned_sample = True
                if code in MISSING_CODES:
                    missing += 1
                    continue
                base = CODE_TO_BASE[code]
                alleles.append(base)
                allele_set.add(base)

            if not alleles:
                masked_positions.append(idx)
                if indel_flags[idx]:
                    counts["masked_indel"] += 1
                elif has_unaligned_sample:
                    counts["masked_no_alignment"] += 1
                else:
                    counts["masked_missingness"] += 1
                continue

            alt_order = sorted(a for a in allele_set if a != ref_base)

            if indel_flags[idx] or (
                args.mask_indel_adjacent_snps and alt_order and adjacent_indel_flags[idx]
            ):
                masked_positions.append(idx)
                counts["masked_indel"] += 1
                continue

            if missing > allowed_missing:
                masked_positions.append(idx)
                if has_unaligned_sample:
                    counts["masked_no_alignment"] += 1
                else:
                    counts["masked_missingness"] += 1
                continue

            if len(alt_order) > 1 and not args.allow_multiallelic_snps:
                masked_positions.append(idx)
                counts["masked_multiallelic"] += 1
                continue

            info = f"NS={len(alleles)};MS={missing}"
            if not alt_order and allele_set == {ref_base}:
                counts["all_sites"] += 1
                counts["invariant"] += 1
                retained_positions.append(idx)
                genotypes = [format_gt(code, alt_order) for code in call_codes]
                if args.add_ref:
                    genotypes.append("0")
                record = (
                    f"{contig}\t{pos}\t.\t{ref_base}\t.\t.\tPASS\t{info};SC=invariant\tGT\t"
                    + "\t".join(genotypes)
                )
                all_sites.write(f"{record}\n")
                continue

            counts["all_sites"] += 1
            counts["variants"] += 1
            retained_positions.append(idx)
            genotypes = [format_gt(code, alt_order) for code in call_codes]
            if args.add_ref:
                genotypes.append("0")
            alt_field = ",".join(alt_order)
            record = (
                f"{contig}\t{pos}\t.\t{ref_base}\t{alt_field}\t.\tPASS\t{info};SC=variant\tGT\t"
                + "\t".join(genotypes)
            )
            all_sites.write(f"{record}\n")
            variants.write(f"{record}\n")

    bed_lines = [f"{chrom}\t{start}\t{end}" for chrom, start, end in intervals_from_positions(contig, masked_positions)]
    write_lines(mask_path, bed_lines)

    summarize_site_and_mask_coverage(
        contig,
        contig_len,
        retained_positions,
        masked_positions,
    )

    summary_lines = [
        "metric\tvalue",
        f"contig\t{contig}",
        f"contig_length\t{contig_len}",
        f"samples\t{len(samples)}",
        f"allowed_missing\t{allowed_missing}",
    ]
    for key in (
        "all_sites",
        "variants",
        "invariant",
        "masked_missingness",
        "masked_indel",
        "masked_multiallelic",
        "masked_no_alignment",
        "masked_ref_non_acgt",
    ):
        summary_lines.append(f"{key}\t{counts.get(key, 0)}")
    summary_lines.append(f"masked_total\t{len(masked_positions)}")
    summary_lines.append(f"masked_intervals\t{len(bed_lines)}")
    write_lines(summary_path, summary_lines)


if __name__ == "__main__":
    main()
