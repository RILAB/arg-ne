#!/usr/bin/env python3
"""Build a reference-anchored FASTA alignment for a single window from
combined.<contig>.all_sites.vcf + per-sample missing BEDs.

Usage:
    python scripts/window_to_fasta.py \\
        --vcf admix_results/sites/combined.10.all_sites.vcf \\
        --bed-dir admix_results/sites \\
        --reference admix/b73.fa \\
        --contig 10 --start 96300001 --end 96400000 \\
        --out admix_results/window_chr10_96.3-96.4Mb.fa
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def fetch_reference(reference: Path, contig: str, start: int, end: int) -> str:
    region = f"{contig}:{start}-{end}"
    out = subprocess.run(
        ["samtools", "faidx", str(reference), region],
        check=True, capture_output=True, text=True,
    ).stdout
    seq = "".join(line.strip() for line in out.splitlines()[1:]).upper()
    expected = end - start + 1
    if len(seq) != expected:
        sys.exit(f"reference length {len(seq)} != expected {expected}")
    return seq


def parse_vcf(vcf: Path, contig: str, start: int, end: int):
    samples: list[str] = []
    records: list[tuple[int, str, list[str], list[str]]] = []
    with open(vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] != contig:
                continue
            pos = int(parts[1])
            if pos < start or pos > end:
                continue
            ref, alts = parts[3], parts[4].split(",")
            gts = parts[9:]
            records.append((pos, ref, alts, gts))
    return samples, records


def apply_missing_bed(bed: Path, seq: bytearray, contig: str, start: int, end: int) -> None:
    with open(bed) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] != contig:
                continue
            bs, be = int(f[1]), int(f[2])
            # BED is 0-based half-open; window in 1-based inclusive [start, end]
            # overlap range (1-based): [max(bs+1, start), min(be, end)]
            lo = max(bs + 1, start)
            hi = min(be, end)
            if lo > hi:
                continue
            for p in range(lo, hi + 1):
                seq[p - start] = ord("N")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", type=Path, required=True)
    ap.add_argument("--bed-dir", type=Path, required=True)
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--contig", required=True)
    ap.add_argument("--start", type=int, required=True, help="1-based inclusive")
    ap.add_argument("--end", type=int, required=True, help="1-based inclusive")
    ap.add_argument("--bed-prefix", default="combined")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    ref_seq = fetch_reference(args.reference, args.contig, args.start, args.end)
    samples, records = parse_vcf(args.vcf, args.contig, args.start, args.end)
    print(f"window length: {len(ref_seq)} bp; samples: {len(samples)}; variant records: {len(records)}", file=sys.stderr)

    out_seqs: dict[str, bytearray] = {s: bytearray(ref_seq, "ascii") for s in samples}

    for pos, ref, alts, gts in records:
        idx = pos - args.start
        alleles = [ref] + alts
        for s, gt in zip(samples, gts):
            if gt == ".":
                out_seqs[s][idx] = ord("N")
                continue
            try:
                ai = int(gt)
            except ValueError:
                out_seqs[s][idx] = ord("N")
                continue
            if ai < 0 or ai >= len(alleles):
                out_seqs[s][idx] = ord("N")
                continue
            allele = alleles[ai]
            # all-sites VCF emits single-base alleles (indels are masked)
            out_seqs[s][idx] = ord(allele[0].upper())

    for s in samples:
        bed = args.bed_dir / f"{args.bed_prefix}.{args.contig}.{s}.missing.bed"
        if not bed.exists():
            sys.exit(f"missing BED not found: {bed}")
        apply_missing_bed(bed, out_seqs[s], args.contig, args.start, args.end)

    with open(args.out, "w") as fh:
        fh.write(f">B73 {args.contig}:{args.start}-{args.end}\n")
        s = ref_seq
        for i in range(0, len(s), 80):
            fh.write(s[i:i + 80] + "\n")
        for name in samples:
            fh.write(f">{name}\n")
            t = out_seqs[name].decode("ascii")
            for i in range(0, len(t), 80):
                fh.write(t[i:i + 80] + "\n")
    print(f"wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
