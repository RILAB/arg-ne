from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import Iterable, TextIO


def open_text(
    path: str | Path,
    mode: str,
    *,
    encoding: str = "utf-8",
    errors: str | None = None,
) -> TextIO:
    path_str = str(path)
    kwargs = {}
    if "b" not in mode:
        kwargs["encoding"] = encoding
        if errors is not None:
            kwargs["errors"] = errors
    if path_str.endswith(".gz"):
        return gzip.open(path_str, mode, **kwargs)  # type: ignore[arg-type]
    return open(path_str, mode, **kwargs)


def open_fasta(path: Path) -> TextIO:
    return open_text(path, "rt")


def normalize_contig(name: str) -> str:
    name = name.strip().lower()
    if "." in name and name.rsplit(".", 1)[-1].startswith("chr"):
        # Drop an assembly/genome prefix that precedes a "chr" contig token,
        # e.g. "Zm-B73v5.chr5" -> "chr5" or "Zx-TIL25.chr2" -> "chr2".
        # Gated on the "chr" token so accession-style names such as
        # "NC_050096.1" (where ".1" is a version suffix) are left untouched.
        name = name.rsplit(".", 1)[-1]
    name = re.sub(r"^chr", "", name, flags=re.IGNORECASE)
    m = re.match(r"^(.*?)(\d+)$", name)
    if m:
        prefix, num = m.groups()
        num = num.lstrip("0") or "0"
        name = f"{prefix}{num}"
    else:
        name = name.lstrip("0")
    return name if name else "0"


def read_fasta_contigs(path: Path) -> list[str]:
    contigs: list[str] = []
    try:
        with open_fasta(path) as handle:
            for line in handle:
                if line.startswith(">"):
                    contigs.append(line[1:].strip().split()[0])
    except OSError:
        pass
    return contigs


def _maf_reference_contigs(lines: Iterable[str]) -> set[str]:
    contigs: set[str] = set()
    first_src: str | None = None
    for line in lines:
        if not line or line.startswith("#"):
            continue
        stripped = line.strip()
        if not stripped:
            if first_src is not None:
                contigs.add(first_src)
                first_src = None
            continue
        parts = stripped.split()
        if not parts:
            continue
        if parts[0] == "a":
            if first_src is not None:
                contigs.add(first_src)
            first_src = None
            continue
        if parts[0] == "s" and len(parts) >= 2 and first_src is None:
            first_src = parts[1]
    if first_src is not None:
        contigs.add(first_src)
    return contigs


def read_maf_contigs(path: Path) -> set[str]:
    try:
        with open_text(path, "rt", errors="ignore") as handle:
            return _maf_reference_contigs(handle)
    except OSError:
        return set()


def extract_info_int(info: str, key: str) -> int | None:
    if info == ".":
        return None
    prefix = f"{key}="
    for field in info.split(";"):
        if field.startswith(prefix):
            try:
                return int(field.split("=", 1)[1])
            except ValueError:
                return None
    return None


def extract_end(info: str) -> int | None:
    return extract_info_int(info, "END")


def extract_dp(info: str) -> int | None:
    return extract_info_int(info, "DP")


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    sorted_intervals = sorted(intervals)
    merged: list[tuple[int, int]] = []
    cur_start, cur_end = sorted_intervals[0]
    for start, end in sorted_intervals[1:]:
        if start <= cur_end:
            cur_end = max(cur_end, end)
            continue
        merged.append((cur_start, cur_end))
        cur_start, cur_end = start, end
    merged.append((cur_start, cur_end))
    return merged
