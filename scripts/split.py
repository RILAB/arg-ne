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
import gzip
import math
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import TextIO

SCRIPT_NAME = "vcf-splitter"
SCRIPT_VERSION = "1.0.0"

# Valid single-base SNP alleles
VALID_BASES = {"A", "C", "G", "T"}

# ------------------------------------------------------------
# Provenance header helper
# ------------------------------------------------------------
def build_provenance_headers(
    in_path: str,
    out_inv: str,
    out_filt: str,
    out_clean: str,
    bgzip_output: bool,
) -> list[str]:
    """
    Build extra VCF header lines describing how the outputs were produced.
    These must appear before the #CHROM header line.
    """
    ts = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    user = getpass.getuser()
    host = platform.node()

    # Record the exact command line for reproducibility.
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
    ]
    lines.append(
        "##notes=Records with valid INFO/END are expanded across POS..END and routed to .inv; other filters apply as documented."
    )
    return [ln + "\n" for ln in lines]



# ------------------------------------------------------------
# Utility: open plain-text or gzipped files transparently
# ------------------------------------------------------------
def open_maybe_gzip(path: str, mode: str) -> TextIO:
    """
    Open a file that may or may not be gzipped.

    Parameters
    ----------
    path : str
        Path to file (.vcf or .vcf.gz)
    mode : str
        Text mode: "rt" for reading, "wt" for writing

    Returns
    -------
    TextIO
    """
    if path.endswith(".gz"):
        return gzip.open(path, mode)  # type: ignore
    return open(path, mode, encoding="utf-8")


def require_tool(tool: str) -> None:
    if shutil.which(tool) is None:
        sys.exit(f"ERROR: required tool not found in PATH: {tool}")


def load_fai_lengths(path: str) -> dict[str, int]:
    """
    Load contig lengths from a .fai file.
    Expects: contig<TAB>length as the first two columns.
    """
    lengths: dict[str, int] = {}
    with open(path, "rt", encoding="utf-8") as fin:
        for raw in fin:
            if not raw.strip():
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            name = cols[0]
            try:
                length = int(cols[1])
            except ValueError:
                continue
            lengths[name] = length
    return lengths


# ------------------------------------------------------------
# INFO field parsers
# ------------------------------------------------------------
def extract_dp(info: str) -> int | None:
    """
    Extract integer DP value from the INFO column.

    Returns
    -------
    int : DP value if present and parseable
    None : if INFO=".", DP missing, or malformed
    """
    if info == ".":
        return None

    for field in info.split(";"):
        if field.startswith("DP="):
            try:
                return int(field.split("=", 1)[1])
            except ValueError:
                return None

    return None


def extract_end(info: str) -> int | None:
    """
    Extract integer END value from the INFO column.

    Returns
    -------
    int : END value if present and parseable
    None : if INFO="." or END missing/malformed
    """
    if info == ".":
        return None

    for field in info.split(";"):
        if field.startswith("END="):
            try:
                return int(field.split("=", 1)[1])
            except ValueError:
                return None

    return None


# ------------------------------------------------------------
# Clean-output ALT normalization
# ------------------------------------------------------------
def format_clean_record(cols: list[str]) -> str:
    """
    Remove <NON_REF> from ALT while preserving all other fields.
    """
    out = list(cols)
    alts = [a for a in out[4].split(",") if a and a != "<NON_REF>"]
    if alts:
        out[4] = ",".join(alts)
    return "\t".join(out) + "\n"


def infer_reference_gt(sample_fields: list[str], default_ploidy: int = 2) -> str:
    """
    Infer a reference-only GT string (e.g. 0/0, 0|0, or 0) from existing samples.
    Defaults to diploid unphased if no usable GT is present.
    """
    for sample in sample_fields:
        gt = sample.split(":", 1)[0]
        if not gt or gt in (".", "./.", ".|."):
            continue
        if "|" in gt:
            parts = gt.split("|")
            if parts:
                return "|".join(["0"] * len(parts))
        if "/" in gt:
            parts = gt.split("/")
            if parts:
                return "/".join(["0"] * len(parts))
        return "0"
    ploidy = max(int(default_ploidy), 1)
    if ploidy == 1:
        return "0"
    return "/".join(["0"] * ploidy)


def build_reference_sample_value(
    format_field: str,
    alt_field: str,
    sample_fields: list[str],
    default_ploidy: int = 2,
) -> str:
    """
    Build a REF sample column that matches the FORMAT layout with GT fixed to reference.
    Non-GT subfields are set to missing ('.').
    """
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
        sample_gt_norm = sample_gt.replace("|", "/")
        gt_value_norm = gt_value.replace("|", "/")
        if sample_gt_norm != gt_value_norm:
            continue
        for idx, key in enumerate(fmt_keys):
            if key in {"AD", "PL", "DP"} and idx < len(sample_parts) and sample_parts[idx] not in ("", "."):
                donor_values[key] = sample_parts[idx]
        break
    alt_count = 0 if alt_field in ("", ".") else len([a for a in alt_field.split(",") if a != ""])
    ploidy = gt_ploidy
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
            pl_len = math.comb(allele_count + ploidy - 1, ploidy)
            out_fields.append(donor_values.get("PL", ",".join(["0"] + (["99"] * max(pl_len - 1, 0)))))
        else:
            out_fields.append(".")
    return ":".join(out_fields)


def _update_info_dp(info: str, dp: int) -> str:
    """
    Replace or append INFO DP with the provided depth.
    """
    if info == "." or info == "":
        return f"DP={dp}"
    parts = info.split(";")
    updated = False
    for i, part in enumerate(parts):
        if part.startswith("DP="):
            parts[i] = f"DP={dp}"
            updated = True
            break
    if not updated:
        parts.append(f"DP={dp}")
    return ";".join(parts)
# ------------------------------------------------------------
# Main driver
# ------------------------------------------------------------
def main() -> None:
    """
    Main entry point.
    Handles argument parsing, file routing, and record classification.
    """
    # Counters track base pairs written to each output.
    inv_bp = 0
    filtered_bp = 0
    clean_bp = 0
    missing_bp = 0
    sample_names: list[str] = []
    missing_gt_snp_total = 0
    missing_gt_snp_by_sample: dict[str, int] = {}


    # ---------------- Argument parsing ----------------
    ap = argparse.ArgumentParser(
        description="Split/filter a VCF into inv/filtered/clean, expanding END spans."
    )
    ap.add_argument(
        "--filter-multiallelic",
        action="store_true",
        help="Filter sites with multiple A/C/G/T alleles (e.g. ALT=T,G)"
    )
    ap.add_argument(
        "vcf",
        help="Input VCF (.vcf or .vcf.gz)"
    )
    ap.add_argument(
        "--out-prefix",
        default=None,
        help="Output prefix (default: input filename without .vcf/.vcf.gz)"
    )
    ap.add_argument(
        "--bgzip-output",
        action="store_true",
        help="Bgzip all output files (.gz)"
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
    args = ap.parse_args()

    in_path = args.vcf
    contig_lengths = load_fai_lengths(args.fai) if args.fai else {}

    # ---------------- Output prefix logic ----------------
    if args.out_prefix:
        prefix = args.out_prefix
    else:
        # Strip .vcf or .vcf.gz from input name
        p = Path(in_path)
        name = p.name
        if name.endswith(".vcf.gz"):
            name = name[:-7]
        elif name.endswith(".vcf"):
            name = name[:-4]
        prefix = str(p.with_name(name))

    out_inv = prefix + ".inv"
    out_filt = prefix + ".filtered"
    out_clean = prefix + ".clean.vcf"
    out_missing = prefix + ".missing.bed"
    out_missing_gt_stats = args.missing_gt_stats_out or (prefix + ".missing_gt_snp_by_sample.tsv")
    bgzip_output = args.bgzip_output
    add_reference = args.add_reference
    reference_ploidy = max(args.reference_ploidy, 1)

    out_inv_final = out_inv + (".gz" if bgzip_output else "")
    out_filt_final = out_filt + (".gz" if bgzip_output else "")
    out_clean_final = out_clean + (".gz" if bgzip_output else "")
    out_missing_final = out_missing + (".gz" if bgzip_output else "")

    # ---------------- Open files ----------------
    if bgzip_output:
        require_tool("bgzip")

    with open_maybe_gzip(in_path, "rt") as fin, \
         open(out_inv, "wt", encoding="utf-8") as f_inv, \
         open(out_filt, "wt", encoding="utf-8") as f_filt, \
         open(out_clean, "wt", encoding="utf-8") as f_clean, \
         open(out_missing, "wt", encoding="utf-8") as f_missing:
    
        # Buffer header lines so we can inject provenance right before #CHROM.
        header_buffer: list[str] = []
        headers_written = False
    
        # Prebuild provenance lines (inserted before #CHROM).
        prov_lines = build_provenance_headers(
            in_path=in_path,
            out_inv=out_inv_final,
            out_filt=out_filt_final,
            out_clean=out_clean_final,
            bgzip_output=bgzip_output,
        )
    
        record_count = 0
        last_report_time = time.time()
        last_chrom: str | None = None
        last_end: int | None = None
        group: list[dict] = []
        group_chrom: str | None = None
        group_pos: int | None = None

        def flush_group(records: list[dict]) -> None:
            nonlocal inv_bp, filtered_bp, clean_bp
            nonlocal missing_gt_snp_total, missing_gt_snp_by_sample
            if not records:
                return
            has_inv = any(r["is_inv"] for r in records)
            has_filtered = any(r["is_filtered"] for r in records)
            has_missing_gt = any(r["has_missing_gt"] for r in records)
            if has_inv:
                for r in records:
                    f_inv.write(r["line"] + "\n")
                    inv_bp += 1
                records.clear()
                return
            if has_filtered or has_missing_gt:
                # Track SNP-like sites that were excluded only because sample GTs were missing.
                if has_missing_gt and not has_filtered:
                    missing_gt_snp_total += 1
                    missing_samples: set[str] = set()
                    for r in records:
                        for idx in r["missing_gt_sample_idxs"]:
                            if idx < len(sample_names):
                                missing_samples.add(sample_names[idx])
                            else:
                                missing_samples.add(f"sample_{idx + 1}")
                    for sample in missing_samples:
                        missing_gt_snp_by_sample[sample] = missing_gt_snp_by_sample.get(sample, 0) + 1
                for r in records:
                    f_filt.write(r["line"] + "\n")
                    filtered_bp += 1
                records.clear()
                return
            for r in records:
                clean_cols = list(r["cols"])
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
                clean_bp += 1
            records.clear()
        for raw in fin:
            record_count += 1
            # Periodic progress line to stderr for large files.
            if record_count % 100_000 == 0:
                fields = raw.rstrip("\n").split("\t")
                if len(fields) >= 2:
                    chrom, pos = fields[0], fields[1]
                else:
                    chrom, pos = "?", "?"

                now = time.time()
                rate = record_count / max(now - last_report_time, 1e-6)

                sys.stderr.write(
                    f"\r[progress] {record_count:,} records  @ chrom {chrom} Mb {int(pos)/1E6}  "
                    f"({rate:,.0f} rec/s)"
                )
                sys.stderr.flush()

            # ---------------- Header handling ----------------
            if raw.startswith("#"):
                # Preserve fileformat line in all outputs.
                if raw.startswith("##fileformat"):
                    f_inv.write(raw)
                    f_filt.write(raw)
                    f_clean.write(raw)
                # Buffer meta/header lines until we see the #CHROM header line.
                elif raw.startswith("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"):
                    header_cols = raw.rstrip("\n").split("\t")
                    sample_names = header_cols[9:] if len(header_cols) > 9 else []
                    if add_reference and "REF" in sample_names:
                        sys.exit('ERROR: --add-reference requested, but input already contains a sample named "REF".')
                    # First: inject provenance headers (##...).
                    for ln in prov_lines:
                        f_inv.write(ln)
                        f_filt.write(ln)
                        f_clean.write(ln)
                    
                    # Second: write buffered input meta-headers (##...).
                    for h in header_buffer:
                        f_inv.write(h)
                        f_filt.write(h)
                        f_clean.write(h)

                    # Third: write the #CHROM header line itself.
                    f_inv.write(raw)
                    f_filt.write(raw)
                    if add_reference:
                        f_clean.write(raw.rstrip("\n") + "\tREF\n")
                    else:
                        f_clean.write(raw)
                    headers_written = True
                else:
                    header_buffer.append(raw)
                continue

            # If input is malformed and lacks #CHROM, flush buffered headers now.
            if not headers_written:
                for h in header_buffer:
                    f_inv.write(h)
                    f_filt.write(h)
                    f_clean.write(h)
                for ln in prov_lines:
                    f_inv.write(ln)
                    f_filt.write(ln)
                    f_clean.write(ln)
                headers_written = True

            # ---------------- Data line processing ----------------
            line = raw.rstrip("\n")
            cols = line.split("\t")

            # Track missing bp: any gap between covered positions.
            # We infer covered span from END= when present, else from REF length.
            if len(cols) >= 2:
                chrom = cols[0]
                try:
                    pos_for_missing = int(cols[1])
                except ValueError:
                    pos_for_missing = None
                if pos_for_missing is not None:
                    # Reset gap tracking when chromosome changes.
                    if chrom != last_chrom:
                        # Close out the previous chromosome tail if we know its length.
                        if last_chrom is not None and last_end is not None:
                            chrom_len = contig_lengths.get(last_chrom)
                            if chrom_len is not None and last_end < chrom_len:
                                f_missing.write(f"{last_chrom}\t{last_end}\t{chrom_len}\n")
                                missing_bp += (chrom_len - last_end)
                        last_chrom = chrom
                        last_end = None

                    ref_for_missing = cols[3] if len(cols) >= 4 else "N"
                    end_for_missing = None
                    if len(cols) >= 8:
                        end_for_missing = extract_end(cols[7])

                    # Determine span covered by this record for gap tracking.
                    if end_for_missing is not None and end_for_missing >= pos_for_missing:
                        curr_start = pos_for_missing
                        curr_end = end_for_missing
                    else:
                        ref_len = len(ref_for_missing) if ref_for_missing else 1
                        curr_start = pos_for_missing
                        curr_end = pos_for_missing + max(ref_len, 1) - 1

                    if last_end is None:
                        # If the contig doesn't start at position 1, add initial missing span.
                        if curr_start > 1:
                            f_missing.write(f"{chrom}\t0\t{curr_start - 1}\n")
                            missing_bp += (curr_start - 1)
                        last_end = curr_end
                    elif curr_start <= last_end + 1:
                        if curr_end > last_end:
                            last_end = curr_end
                    else:
                        # Gap detected: write missing BED interval (0-based, half-open).
                        gap_start = last_end + 1
                        gap_end = curr_start - 1
                        if gap_end >= gap_start:
                            f_missing.write(f"{chrom}\t{gap_start - 1}\t{gap_end}\n")
                            missing_bp += (gap_end - gap_start + 1)
                        last_end = curr_end

            # Malformed line (no INFO column) => filtered
            if len(cols) < 8:
                f_filt.write(line + "\n")
                filtered_bp += 1
                continue

            # VCF columns (0-based): 0 CHROM, 1 POS, 3 REF, 4 ALT, 7 INFO
            ref = cols[3]
            alt_field = cols[4]
            info = cols[7]

            # ============================================================
            # 1) .inv (highest priority): ALT="." (reference) OR contains END=
            # If END= exists, expand across POS..END inclusive.
            # ============================================================
            end_val = extract_end(info)
            if end_val is not None:
                # Flush any buffered records at the previous position.
                flush_group(group)
                if end_val is None:
                    f_inv.write(line + "\n")
                    inv_bp += 1
                    continue

                try:
                    pos0 = int(cols[1])
                except ValueError:
                    f_inv.write(line + "\n")
                    inv_bp += 1
                    continue

                if end_val <= pos0:
                    f_inv.write(line + "\n")
                    inv_bp += 1
                    continue

                span = end_val - pos0 + 1
                inv_bp += span
                for pos in range(pos0, end_val + 1):
                    cols[1] = str(pos)
                    f_inv.write("\t".join(cols) + "\n")
                continue

            # ============================================================
            # 2) .filtered
            # ============================================================
            alts = [a.strip() for a in alt_field.split(",")] if alt_field != "." else []
            alts_no_nonref = [a for a in alts if a != "<NON_REF>"]

            # Depth filtering disabled: DP is parsed only for downstream metadata.
            dp = extract_dp(info)
            has_star = ("*" in alts_no_nonref)
            ref_long = (len(ref) > 1)

            valid_base_alts = [a for a in alts_no_nonref if a in VALID_BASES]
            multiple_valid_bases = (len(set(valid_base_alts)) > 1)

            has_non_acgt_nonstar = any(
                (a not in VALID_BASES) and (a != "*") for a in alts_no_nonref
            )

            is_filtered = (
                has_star
                or ref_long
                or (args.filter_multiallelic and multiple_valid_bases)
                or has_non_acgt_nonstar
            )
            # Treat ALT="." or ALT="<NON_REF>"-only (gVCF-style invariant records) as invariant.
            is_inv = (alt_field == ".") or (len(alts) > 0 and len(alts_no_nonref) == 0)

            # Group by chrom/pos to enforce mutual exclusivity.
            if group_pos is None:
                group_chrom = chrom
                group_pos = int(cols[1])
            if chrom != group_chrom or int(cols[1]) != group_pos:
                flush_group(group)
                group_chrom = chrom
                group_pos = int(cols[1])
            gts = [c.split(":", 1)[0] for c in cols[9:]] if len(cols) > 9 else ["."]
            missing_gt_sample_idxs = [
                i for i, gt in enumerate(gts)
                if ("." in gt) or (gt == "")
            ]
            has_missing_gt = len(missing_gt_sample_idxs) > 0
            group.append(
                {
                    "line": line,
                    "cols": cols,
                    "is_inv": is_inv,
                    "is_filtered": is_filtered,
                    "has_missing_gt": has_missing_gt,
                    "missing_gt_sample_idxs": missing_gt_sample_idxs,
                }
            )

        flush_group(group)

        print("Output summary (non-header bp):",file=sys.stderr)
        print(f"  inv:      {inv_bp:,}",file=sys.stderr)
        print(f"  filtered: {filtered_bp:,}",file=sys.stderr)
        print(f"  clean:    {clean_bp:,}",file=sys.stderr)
        # Close out the final chromosome tail if we know its length.
        if last_chrom is not None and last_end is not None:
            chrom_len = contig_lengths.get(last_chrom)
            if chrom_len is not None and last_end < chrom_len:
                f_missing.write(f"{last_chrom}\t{last_end}\t{chrom_len}\n")
                missing_bp += (chrom_len - last_end)

        print(f"  missing:  {missing_bp:,}",file=sys.stderr)

    with open(out_missing_gt_stats, "w", encoding="utf-8") as f_stats:
        f_stats.write("sample\texcluded_snp_sites\ttotal_excluded_snp_sites\n")
        names = sample_names if sample_names else sorted(missing_gt_snp_by_sample.keys())
        for sample in names:
            count = missing_gt_snp_by_sample.get(sample, 0)
            f_stats.write(f"{sample}\t{count}\t{missing_gt_snp_total}\n")

    if bgzip_output:
        for path in (out_inv, out_filt, out_clean, out_missing):
            proc = subprocess.run(["bgzip", "-f", path], check=False)
            if proc.returncode != 0:
                sys.exit(f"ERROR: bgzip failed for {path}")
    
if __name__ == "__main__":
    main()
