#!/usr/bin/env python3
"""
Check that a direct all-sites VCF and mask BED together cover the full contig.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple

try:
    from scripts.common import merge_intervals, open_text
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from scripts.common import merge_intervals, open_text


def load_fai_length(path: Path, chrom: str) -> int:
    with path.open("r", encoding="utf-8") as fin:
        for raw in fin:
            if not raw.strip():
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) >= 2 and cols[0] == chrom:
                return int(cols[1])
    raise ValueError(f"Chromosome '{chrom}' not found in {path}")


def read_site_intervals(path: Path) -> tuple[str | None, List[Tuple[int, int]]]:
    chrom = None
    intervals: List[Tuple[int, int]] = []
    with open_text(path, "rt") as fin:
        for raw in fin:
            if raw.startswith("#") or not raw.strip():
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            if chrom is None:
                chrom = cols[0]
            pos = int(cols[1])
            intervals.append((pos - 1, pos))
    return chrom, intervals


def read_bed_intervals(path: Path) -> tuple[str | None, List[Tuple[int, int]]]:
    chrom = None
    intervals: List[Tuple[int, int]] = []
    with path.open("r", encoding="utf-8") as fin:
        for raw in fin:
            if not raw.strip() or raw.startswith("#"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            if chrom is None:
                chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            if end > start:
                intervals.append((start, end))
    return chrom, intervals


def overlap_bp(a: List[Tuple[int, int]], b: List[Tuple[int, int]]) -> int:
    i = j = 0
    total = 0
    a = merge_intervals(a)
    b = merge_intervals(b)
    while i < len(a) and j < len(b):
        a_s, a_e = a[i]
        b_s, b_e = b[j]
        if a_e <= b_s:
            i += 1
            continue
        if b_e <= a_s:
            j += 1
            continue
        total += min(a_e, b_e) - max(a_s, b_s)
        if a_e <= b_e:
            i += 1
        else:
            j += 1
    return total


def summarize_site_and_mask_coverage(
    site_vcf: Path,
    mask_bed: Path,
    fai: Path,
    *,
    chrom_hint: str | None = None,
) -> Dict[str, int | str]:
    site_chrom, site_intervals = read_site_intervals(site_vcf)
    bed_chrom, bed_intervals = read_bed_intervals(mask_bed)

    chrom = site_chrom or bed_chrom or chrom_hint
    if chrom is None:
        raise ValueError("unable to determine chromosome from site VCF or mask BED")
    if site_chrom is not None and site_chrom != chrom:
        raise ValueError(f"mismatched chromosome in site VCF: {site_chrom} != {chrom}")
    if bed_chrom is not None and bed_chrom != chrom:
        raise ValueError(f"mismatched chromosome in mask BED: {bed_chrom} != {chrom}")

    chrom_len = load_fai_length(fai, chrom)
    overlap = overlap_bp(site_intervals, bed_intervals)
    if overlap:
        raise ValueError(f"overlap detected for {chrom}: site_vcf vs mask_bed={overlap}")

    merged = merge_intervals(site_intervals + bed_intervals)
    total = sum(e - s for s, e in merged)
    if total != chrom_len:
        raise ValueError(f"coverage mismatch for {chrom}: union={total}, chrom_len={chrom_len}")

    return {
        "chrom": chrom,
        "site_bp": sum(e - s for s, e in merge_intervals(site_intervals)),
        "mask_bed_bp": sum(e - s for s, e in merge_intervals(bed_intervals)),
        "total_bp": total,
        "chrom_len": chrom_len,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Check direct all-sites coverage against reference length.")
    ap.add_argument("--site-vcf", required=True, help="All-sites VCF path")
    ap.add_argument("--mask-bed", required=True, help="Mask BED path")
    ap.add_argument("--fai", required=True, help="Reference .fai path")
    ap.add_argument("--chrom", help="Optional explicit chromosome/contig name")
    ap.add_argument("--report-out", help="Optional report output path")
    args = ap.parse_args()

    try:
        summary = summarize_site_and_mask_coverage(
            Path(args.site_vcf),
            Path(args.mask_bed),
            Path(args.fai),
            chrom_hint=args.chrom,
        )
    except ValueError as exc:
        sys.exit(f"ERROR: {exc}")

    report = Path(args.report_out) if args.report_out else Path(args.site_vcf).with_suffix(".coverage.txt")
    report.write_text(
        f"chrom={summary['chrom']}\n"
        f"site_bp={summary['site_bp']}\n"
        f"mask_bed_bp={summary['mask_bed_bp']}\n"
        f"total_bp={summary['total_bp']}\n"
        f"chrom_len={summary['chrom_len']}\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
