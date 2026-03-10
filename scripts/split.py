#!/usr/bin/env python3
"""
VCF splitter / filter with END-span expansion.

Purpose
-------
Split a gVCF into three mutually exclusive outputs:
  <prefix>.inv       invariant sites (INFO="." or END=... spans)
  <prefix>.filtered  sites removed for quality/format reasons
  <prefix>.clean.vcf usable variant sites for downstream inference

The script also emits <prefix>.missing.bed, a mask of positions absent from
the input gVCF (gaps between covered positions). This helps track accessiblity.

Key behavior
------------
Records with END= in INFO are expanded so each base in [POS, END] is written
as a separate record in .inv. This can be large.
"""

from __future__ import annotations

import argparse
import datetime
import getpass
import math
import os
import platform
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO

try:
    from scripts.common import extract_dp, extract_end, open_text
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from scripts.common import extract_dp, extract_end, open_text

SCRIPT_NAME = "vcf-splitter"
SCRIPT_VERSION = "1.0.0"
VALID_BASES = {"A", "C", "G", "T"}


@dataclass(frozen=True)
class OutputPaths:
    prefix: str
    inv: str
    filtered: str
    clean: str
    missing: str
    missing_gt_stats: str
    inv_final: str
    filtered_final: str
    clean_final: str
    missing_final: str


@dataclass
class Counters:
    inv_bp: int = 0
    filtered_bp: int = 0
    clean_bp: int = 0
    missing_bp: int = 0
    missing_gt_snp_total: int = 0
    missing_gt_snp_by_sample: dict[str, int] = field(default_factory=dict)


@dataclass
class RecordDecision:
    line: str
    cols: list[str]
    is_inv: bool
    is_filtered: bool
    has_missing_gt: bool
    missing_gt_sample_idxs: list[int]


@dataclass
class MissingCoverageTracker:
    contig_lengths: dict[str, int]
    last_chrom: str | None = None
    last_end: int | None = None

    def update(self, cols: list[str], handle: TextIO, counters: Counters) -> None:
        if len(cols) < 2:
            return
        chrom = cols[0]
        try:
            pos = int(cols[1])
        except ValueError:
            return

        if chrom != self.last_chrom:
            self.finish_contig(handle, counters)
            self.last_chrom = chrom
            self.last_end = None

        ref = cols[3] if len(cols) >= 4 else "N"
        end_val = extract_end(cols[7]) if len(cols) >= 8 else None
        if end_val is not None and end_val >= pos:
            curr_start = pos
            curr_end = end_val
        else:
            curr_start = pos
            curr_end = pos + max(len(ref), 1) - 1

        if self.last_end is None:
            if curr_start > 1:
                handle.write(f"{chrom}\t0\t{curr_start - 1}\n")
                counters.missing_bp += curr_start - 1
            self.last_end = curr_end
            return

        if curr_start <= self.last_end + 1:
            self.last_end = max(self.last_end, curr_end)
            return

        gap_start = self.last_end + 1
        gap_end = curr_start - 1
        if gap_end >= gap_start:
            handle.write(f"{chrom}\t{gap_start - 1}\t{gap_end}\n")
            counters.missing_bp += gap_end - gap_start + 1
        self.last_end = curr_end

    def finish_contig(self, handle: TextIO, counters: Counters) -> None:
        if self.last_chrom is None or self.last_end is None:
            return
        chrom_len = self.contig_lengths.get(self.last_chrom)
        if chrom_len is not None and self.last_end < chrom_len:
            handle.write(f"{self.last_chrom}\t{self.last_end}\t{chrom_len}\n")
            counters.missing_bp += chrom_len - self.last_end

    def finish(self, handle: TextIO, counters: Counters) -> None:
        self.finish_contig(handle, counters)


def build_provenance_headers(
    in_path: str,
    out_inv: str,
    out_filt: str,
    out_clean: str,
    bgzip_output: bool,
) -> list[str]:
    ts = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    user = getpass.getuser()
    host = platform.node()
    cmdline = " ".join([os.path.basename(sys.argv[0])] + sys.argv[1:])
    lines = [
        f"##source={SCRIPT_NAME} v{SCRIPT_VERSION}",
        f"##run.user={user}",
        f"##run.host={host}",
        f"##run.timestamp={ts}",
        f"##run.commandline={cmdline}",
        f"##input.file={in_path}",
        f"##output.inv={out_inv}",
        f"##output.filtered={out_filt}",
        f"##output.clean={out_clean}",
        f"##parameters.bgzip_output={'true' if bgzip_output else 'false'}",
        "##notes=Records with valid INFO/END are expanded across POS..END and routed to .inv; other filters apply as documented.",
    ]
    return [line + "\n" for line in lines]


def require_tool(tool: str) -> None:
    if shutil.which(tool) is None:
        sys.exit(f"ERROR: required tool not found in PATH: {tool}")


def load_fai_lengths(path: str) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with open(path, "rt", encoding="utf-8") as fin:
        for raw in fin:
            if not raw.strip():
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            try:
                lengths[cols[0]] = int(cols[1])
            except ValueError:
                continue
    return lengths


def format_clean_record(cols: list[str]) -> str:
    out = list(cols)
    alts = [alt for alt in out[4].split(",") if alt and alt != "<NON_REF>"]
    if alts:
        out[4] = ",".join(alts)
    return "\t".join(out) + "\n"


def infer_reference_gt(sample_fields: list[str], default_ploidy: int = 2) -> str:
    for sample in sample_fields:
        gt = sample.split(":", 1)[0]
        if not gt or gt in (".", "./.", ".|."):
            continue
        if "|" in gt:
            return "|".join(["0"] * len(gt.split("|")))
        if "/" in gt:
            return "/".join(["0"] * len(gt.split("/")))
        return "0"
    ploidy = max(int(default_ploidy), 1)
    return "0" if ploidy == 1 else "/".join(["0"] * ploidy)


def build_reference_sample_value(
    format_field: str,
    alt_field: str,
    sample_fields: list[str],
    default_ploidy: int = 2,
) -> str:
    if format_field in ("", "."):
        return "."
    fmt_keys = format_field.split(":")
    gt_value = infer_reference_gt(sample_fields, default_ploidy=default_ploidy)
    gt_ploidy = max(gt_value.count("/") + gt_value.count("|") + 1, 1)
    donor_values: dict[str, str] = {}
    for sample in sample_fields:
        sample_parts = sample.split(":")
        if not sample_parts:
            continue
        sample_gt = sample_parts[0]
        if not sample_gt or "." in sample_gt:
            continue
        if sample_gt.replace("|", "/") != gt_value.replace("|", "/"):
            continue
        for idx, key in enumerate(fmt_keys):
            if key in {"AD", "PL", "DP"} and idx < len(sample_parts) and sample_parts[idx] not in ("", "."):
                donor_values[key] = sample_parts[idx]
        break

    alt_count = 0 if alt_field in ("", ".") else len([alt for alt in alt_field.split(",") if alt != ""])
    allele_count = 1 + alt_count
    out_fields: list[str] = []
    for key in fmt_keys:
        if key == "GT":
            out_fields.append(gt_value)
        elif key == "DP":
            out_fields.append(donor_values.get("DP", "1"))
        elif key == "AD":
            out_fields.append(donor_values.get("AD", ",".join(["1"] + (["0"] * alt_count))))
        elif key == "PL":
            pl_len = math.comb(allele_count + gt_ploidy - 1, gt_ploidy)
            out_fields.append(donor_values.get("PL", ",".join(["0"] + (["99"] * max(pl_len - 1, 0)))))
        else:
            out_fields.append(".")
    return ":".join(out_fields)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Split/filter a VCF into inv/filtered/clean, expanding END spans."
    )
    ap.add_argument(
        "--filter-multiallelic",
        action="store_true",
        help="Filter sites with multiple A/C/G/T alleles (e.g. ALT=T,G)",
    )
    ap.add_argument("vcf", help="Input VCF (.vcf or .vcf.gz)")
    ap.add_argument(
        "--out-prefix",
        default=None,
        help="Output prefix (default: input filename without .vcf/.vcf.gz)",
    )
    ap.add_argument(
        "--bgzip-output",
        action="store_true",
        help="Bgzip all output files (.gz)",
    )
    ap.add_argument(
        "--add-reference",
        action="store_true",
        help='Append a synthetic "REF" sample to .clean.vcf records with GT set to reference alleles.',
    )
    ap.add_argument(
        "--reference-ploidy",
        type=int,
        default=2,
        help="Fallback ploidy for synthetic REF GT when sample GT ploidy cannot be inferred at a site (default: 2).",
    )
    ap.add_argument(
        "--fai",
        default=None,
        help="Reference .fai to fill missing BED gaps at contig ends.",
    )
    ap.add_argument(
        "--missing-gt-stats-out",
        default=None,
        help=(
            "Optional TSV output path for per-sample counts of SNP sites excluded "
            "from clean due to missing genotype calls."
        ),
    )
    return ap.parse_args()


def resolve_output_paths(args: argparse.Namespace, in_path: str) -> OutputPaths:
    if args.out_prefix:
        prefix = args.out_prefix
    else:
        path = Path(in_path)
        name = path.name
        if name.endswith(".vcf.gz"):
            name = name[:-7]
        elif name.endswith(".vcf"):
            name = name[:-4]
        prefix = str(path.with_name(name))

    inv = prefix + ".inv"
    filtered = prefix + ".filtered"
    clean = prefix + ".clean.vcf"
    missing = prefix + ".missing.bed"
    missing_gt_stats = args.missing_gt_stats_out or (prefix + ".missing_gt_snp_by_sample.tsv")
    suffix = ".gz" if args.bgzip_output else ""
    return OutputPaths(
        prefix=prefix,
        inv=inv,
        filtered=filtered,
        clean=clean,
        missing=missing,
        missing_gt_stats=missing_gt_stats,
        inv_final=inv + suffix,
        filtered_final=filtered + suffix,
        clean_final=clean + suffix,
        missing_final=missing + suffix,
    )


def classify_record(cols: list[str], filter_multiallelic: bool) -> RecordDecision:
    line = "\t".join(cols)
    if len(cols) < 8:
        return RecordDecision(
            line=line,
            cols=cols,
            is_inv=False,
            is_filtered=True,
            has_missing_gt=False,
            missing_gt_sample_idxs=[],
        )

    ref = cols[3]
    alt_field = cols[4]
    info = cols[7]
    alts = [alt.strip() for alt in alt_field.split(",")] if alt_field != "." else []
    alts_no_nonref = [alt for alt in alts if alt != "<NON_REF>"]
    extract_dp(info)

    valid_base_alts = [alt for alt in alts_no_nonref if alt in VALID_BASES]
    has_non_acgt_nonstar = any(alt not in VALID_BASES and alt != "*" for alt in alts_no_nonref)
    is_filtered = (
        ("*" in alts_no_nonref)
        or (len(ref) > 1)
        or (filter_multiallelic and len(set(valid_base_alts)) > 1)
        or has_non_acgt_nonstar
    )
    is_inv = (alt_field == ".") or (len(alts) > 0 and len(alts_no_nonref) == 0)

    gts = [col.split(":", 1)[0] for col in cols[9:]] if len(cols) > 9 else ["."]
    missing_gt_sample_idxs = [idx for idx, gt in enumerate(gts) if ("." in gt) or (gt == "")]
    return RecordDecision(
        line=line,
        cols=cols,
        is_inv=is_inv,
        is_filtered=is_filtered,
        has_missing_gt=bool(missing_gt_sample_idxs),
        missing_gt_sample_idxs=missing_gt_sample_idxs,
    )


def write_headers(
    raw: str,
    *,
    header_buffer: list[str],
    provenance_lines: list[str],
    f_inv: TextIO,
    f_filt: TextIO,
    f_clean: TextIO,
    add_reference: bool,
) -> list[str]:
    header_cols = raw.rstrip("\n").split("\t")
    sample_names = header_cols[9:] if len(header_cols) > 9 else []
    if add_reference and "REF" in sample_names:
        sys.exit('ERROR: --add-reference requested, but input already contains a sample named "REF".')

    for handle in (f_inv, f_filt, f_clean):
        for line in provenance_lines:
            handle.write(line)
        for line in header_buffer:
            handle.write(line)

    f_inv.write(raw)
    f_filt.write(raw)
    f_clean.write(raw.rstrip("\n") + "\tREF\n" if add_reference else raw)
    return sample_names


def flush_buffered_headers(
    header_buffer: list[str],
    provenance_lines: list[str],
    f_inv: TextIO,
    f_filt: TextIO,
    f_clean: TextIO,
) -> None:
    for handle in (f_inv, f_filt, f_clean):
        for line in header_buffer:
            handle.write(line)
        for line in provenance_lines:
            handle.write(line)


def write_missing_gt_stats(
    path: str,
    sample_names: list[str],
    counters: Counters,
) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("sample\texcluded_snp_sites\ttotal_excluded_snp_sites\n")
        names = sample_names if sample_names else sorted(counters.missing_gt_snp_by_sample)
        for sample in names:
            handle.write(
                f"{sample}\t{counters.missing_gt_snp_by_sample.get(sample, 0)}\t"
                f"{counters.missing_gt_snp_total}\n"
            )


def bgzip_outputs(paths: list[str]) -> None:
    for path in paths:
        proc = subprocess.run(["bgzip", "-f", path], check=False)
        if proc.returncode != 0:
            sys.exit(f"ERROR: bgzip failed for {path}")


def expand_end_record(cols: list[str], handle: TextIO) -> int:
    try:
        start = int(cols[1])
    except ValueError:
        handle.write("\t".join(cols) + "\n")
        return 1

    end_val = extract_end(cols[7])
    if end_val is None or end_val <= start:
        handle.write("\t".join(cols) + "\n")
        return 1

    for pos in range(start, end_val + 1):
        cols[1] = str(pos)
        handle.write("\t".join(cols) + "\n")
    return end_val - start + 1


def record_progress(record_count: int, raw: str, start_time: float) -> None:
    if record_count % 100_000 != 0:
        return
    fields = raw.rstrip("\n").split("\t")
    chrom, pos = (fields[0], fields[1]) if len(fields) >= 2 else ("?", "0")
    elapsed = max(time.time() - start_time, 1e-6)
    sys.stderr.write(
        f"\r[progress] {record_count:,} records  @ chrom {chrom} Mb {int(pos)/1E6}  "
        f"({record_count / elapsed:,.0f} rec/s)"
    )
    sys.stderr.flush()


def emit_group(
    records: list[RecordDecision],
    *,
    sample_names: list[str],
    f_inv: TextIO,
    f_filt: TextIO,
    f_clean: TextIO,
    counters: Counters,
    add_reference: bool,
    reference_ploidy: int,
) -> None:
    if not records:
        return

    if any(record.is_inv for record in records):
        for record in records:
            f_inv.write(record.line + "\n")
            counters.inv_bp += 1
        records.clear()
        return

    has_filtered = any(record.is_filtered for record in records)
    has_missing_gt = any(record.has_missing_gt for record in records)
    if has_filtered or has_missing_gt:
        if has_missing_gt and not has_filtered:
            counters.missing_gt_snp_total += 1
            missing_samples: set[str] = set()
            for record in records:
                for idx in record.missing_gt_sample_idxs:
                    if idx < len(sample_names):
                        missing_samples.add(sample_names[idx])
                    else:
                        missing_samples.add(f"sample_{idx + 1}")
            for sample in missing_samples:
                counters.missing_gt_snp_by_sample[sample] = counters.missing_gt_snp_by_sample.get(sample, 0) + 1

        for record in records:
            f_filt.write(record.line + "\n")
            counters.filtered_bp += 1
        records.clear()
        return

    for record in records:
        clean_cols = list(record.cols)
        if add_reference and len(clean_cols) >= 9:
            clean_cols.append(
                build_reference_sample_value(
                    clean_cols[8],
                    clean_cols[4],
                    clean_cols[9:],
                    default_ploidy=reference_ploidy,
                )
            )
        f_clean.write(format_clean_record(clean_cols))
        counters.clean_bp += 1
    records.clear()


def main() -> None:
    args = parse_args()
    if args.bgzip_output:
        require_tool("bgzip")

    in_path = args.vcf
    outputs = resolve_output_paths(args, in_path)
    contig_lengths = load_fai_lengths(args.fai) if args.fai else {}
    counters = Counters()
    tracker = MissingCoverageTracker(contig_lengths)
    provenance_lines = build_provenance_headers(
        in_path=in_path,
        out_inv=outputs.inv_final,
        out_filt=outputs.filtered_final,
        out_clean=outputs.clean_final,
        bgzip_output=args.bgzip_output,
    )

    sample_names: list[str] = []
    header_buffer: list[str] = []
    headers_written = False
    group: list[RecordDecision] = []
    group_key: tuple[str, int] | None = None
    start_time = time.time()

    with open_text(in_path, "rt") as fin, \
         open(outputs.inv, "wt", encoding="utf-8") as f_inv, \
         open(outputs.filtered, "wt", encoding="utf-8") as f_filt, \
         open(outputs.clean, "wt", encoding="utf-8") as f_clean, \
         open(outputs.missing, "wt", encoding="utf-8") as f_missing:
        for record_count, raw in enumerate(fin, start=1):
            record_progress(record_count, raw, start_time)

            if raw.startswith("#"):
                if raw.startswith("##fileformat"):
                    f_inv.write(raw)
                    f_filt.write(raw)
                    f_clean.write(raw)
                elif raw.startswith("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"):
                    sample_names = write_headers(
                        raw,
                        header_buffer=header_buffer,
                        provenance_lines=provenance_lines,
                        f_inv=f_inv,
                        f_filt=f_filt,
                        f_clean=f_clean,
                        add_reference=args.add_reference,
                    )
                    headers_written = True
                else:
                    header_buffer.append(raw)
                continue

            if not headers_written:
                flush_buffered_headers(header_buffer, provenance_lines, f_inv, f_filt, f_clean)
                headers_written = True

            cols = raw.rstrip("\n").split("\t")
            tracker.update(cols, f_missing, counters)

            if len(cols) < 8:
                f_filt.write(raw.rstrip("\n") + "\n")
                counters.filtered_bp += 1
                continue

            if extract_end(cols[7]) is not None:
                emit_group(
                    group,
                    sample_names=sample_names,
                    f_inv=f_inv,
                    f_filt=f_filt,
                    f_clean=f_clean,
                    counters=counters,
                    add_reference=args.add_reference,
                    reference_ploidy=max(args.reference_ploidy, 1),
                )
                group_key = None
                counters.inv_bp += expand_end_record(cols, f_inv)
                continue

            decision = classify_record(cols, args.filter_multiallelic)
            try:
                current_key = (cols[0], int(cols[1]))
            except ValueError:
                f_filt.write(decision.line + "\n")
                counters.filtered_bp += 1
                continue

            if group_key is None:
                group_key = current_key
            elif current_key != group_key:
                emit_group(
                    group,
                    sample_names=sample_names,
                    f_inv=f_inv,
                    f_filt=f_filt,
                    f_clean=f_clean,
                    counters=counters,
                    add_reference=args.add_reference,
                    reference_ploidy=max(args.reference_ploidy, 1),
                )
                group_key = current_key
            group.append(decision)

        emit_group(
            group,
            sample_names=sample_names,
            f_inv=f_inv,
            f_filt=f_filt,
            f_clean=f_clean,
            counters=counters,
            add_reference=args.add_reference,
            reference_ploidy=max(args.reference_ploidy, 1),
        )
        tracker.finish(f_missing, counters)

        print("Output summary (non-header bp):", file=sys.stderr)
        print(f"  inv:      {counters.inv_bp:,}", file=sys.stderr)
        print(f"  filtered: {counters.filtered_bp:,}", file=sys.stderr)
        print(f"  clean:    {counters.clean_bp:,}", file=sys.stderr)
        print(f"  missing:  {counters.missing_bp:,}", file=sys.stderr)

    write_missing_gt_stats(outputs.missing_gt_stats, sample_names, counters)
    if args.bgzip_output:
        bgzip_outputs([outputs.inv, outputs.filtered, outputs.clean, outputs.missing])


if __name__ == "__main__":
    main()
