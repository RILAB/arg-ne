#!/usr/bin/env python3
from __future__ import annotations

import argparse
import html
import re
import subprocess
from datetime import datetime
from pathlib import Path

try:
    from scripts.common import merge_intervals, open_text
except ModuleNotFoundError:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from scripts.common import merge_intervals, open_text


def read_fai_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                lengths[parts[0]] = int(parts[1])
    return lengths


def read_site_summary(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline()
        if not header:
            return values
        for line in handle:
            if not line.strip():
                continue
            key, value = line.rstrip("\n").split("\t", 1)
            values[key] = value
    return values


def window_count(length: int, window_bp: int) -> int:
    return (length + window_bp - 1) // window_bp


def window_midpoints(length: int, window_bp: int) -> list[int]:
    mids: list[int] = []
    for idx in range(window_count(length, window_bp)):
        start = idx * window_bp
        end = min(start + window_bp, length)
        mids.append(start + (end - start) // 2)
    return mids


def read_all_sites_stats(
    path: Path,
    lengths: dict[str, int],
    window_bp: int,
) -> tuple[
    dict[str, list[int]],
    dict[str, list[int]],
    list[str],
    dict[str, dict[str, int]],
    dict[str, dict[str, int]],
    dict[str, dict[str, int]],
]:
    """Single pass over one all_sites VCF.

    Returns:
      invariant_by_contig: window-binned invariant site counts
      variant_by_contig:   window-binned variant site counts
      samples:             sample column order from the VCF header
      per_sample_variant_by_contig: contig -> sample -> sites where sample carries non-ref
      per_sample_missing_retained_by_contig: contig -> sample -> retained sites where sample is missing (GT=.)
      per_sample_called_by_contig: contig -> sample -> retained sites where sample has a call (ref or alt)
    """
    invariant: dict[str, list[int]] = {}
    variant: dict[str, list[int]] = {}
    samples: list[str] = []
    per_sample_variant: dict[str, dict[str, int]] = {}
    per_sample_missing: dict[str, dict[str, int]] = {}
    per_sample_called: dict[str, dict[str, int]] = {}

    with open_text(path, "rt", errors="ignore") as handle:
        for line in handle:
            if not line.strip():
                continue
            if line.startswith("#CHROM"):
                header_parts = line.rstrip("\n").split("\t")
                if len(header_parts) > 9:
                    samples = header_parts[9:]
                continue
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            contig = parts[0]
            length = lengths.get(contig)
            if length is None:
                continue
            if contig not in invariant:
                invariant[contig] = [0] * window_count(length, window_bp)
                variant[contig] = [0] * window_count(length, window_bp)
            pos = int(parts[1])
            idx = max((pos - 1) // window_bp, 0)
            if idx >= len(invariant[contig]):
                continue
            is_variant = parts[4] != "."
            if is_variant:
                variant[contig][idx] += 1
            else:
                invariant[contig][idx] += 1

            if samples and len(parts) >= 9 + len(samples):
                sv = per_sample_variant.setdefault(contig, {s: 0 for s in samples})
                sm = per_sample_missing.setdefault(contig, {s: 0 for s in samples})
                sc = per_sample_called.setdefault(contig, {s: 0 for s in samples})
                for s_idx, sample in enumerate(samples):
                    gt = parts[9 + s_idx].partition(":")[0]
                    if gt == ".":
                        sm[sample] += 1
                    else:
                        sc[sample] += 1
                        if is_variant and gt != "0":
                            sv[sample] += 1

    return invariant, variant, samples, per_sample_variant, per_sample_missing, per_sample_called


def read_sample_missing_bp(
    bed_paths: list[Path],
    samples: list[str],
    lengths: dict[str, int],
) -> dict[str, dict[str, int]]:
    """Sum BED interval lengths per sample per contig.

    Sample name is inferred from each BED path's filename suffix
    ".{sample}.missing.bed". Paths that don't match any known sample are ignored.
    """
    # Sort samples by length desc so longer names match first (avoids a sample
    # name being a suffix of another's).
    ordered = sorted(samples, key=len, reverse=True)
    data: dict[str, dict[str, int]] = {}
    for bed_path in bed_paths:
        name = bed_path.name
        matched = None
        for s in ordered:
            if name.endswith(f".{s}.missing.bed"):
                matched = s
                break
        if matched is None:
            continue
        data.setdefault(matched, {})
        with bed_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                contig = parts[0]
                if contig not in lengths:
                    continue
                start = int(parts[1])
                end = int(parts[2])
                span = max(end - start, 0)
                if span:
                    data[matched][contig] = data[matched].get(contig, 0) + span
    return data


def read_mask_percentages(
    path: Path,
    lengths: dict[str, int],
    window_bp: int,
) -> dict[str, list[int]]:
    masked: dict[str, list[int]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            contig = parts[0]
            length = lengths.get(contig)
            if length is None:
                continue
            if contig not in masked:
                masked[contig] = [0] * window_count(length, window_bp)
            start = int(parts[1])
            end = int(parts[2])
            pos = start
            while pos < end:
                idx = pos // window_bp
                if idx >= len(masked[contig]):
                    break
                window_end = min((idx + 1) * window_bp, length)
                span_end = min(end, window_end)
                masked[contig][idx] += max(span_end - pos, 0)
                pos = span_end
    return masked


def to_percentages(counts: list[int], length: int, window_bp: int) -> list[float]:
    values: list[float] = []
    for idx, count in enumerate(counts):
        start = idx * window_bp
        end = min(start + window_bp, length)
        span = max(end - start, 1)
        values.append((count / span) * 100.0)
    return values


def fmt_int(s: str) -> str:
    try:
        return f"{int(s):,}"
    except (ValueError, TypeError):
        return s


def _anchor(contig: str) -> str:
    return "c-" + re.sub(r"[^\w-]", "_", contig)


def svg_combined_plot(
    xs: list[int],
    series: list[tuple[str, str, list[float]]],
    *,
    width: int = 900,
    height: int = 300,
) -> str:
    """Multi-series line plot; y-axis fixed 0–100% for comparability across contigs.

    Uses polylines only (no per-point circles) to keep SVG element count low
    even for large chromosomes with thousands of windows.
    """
    margin = {"left": 60, "right": 20, "top": 28, "bottom": 55}
    plot_w = width - margin["left"] - margin["right"]
    plot_h = height - margin["top"] - margin["bottom"]
    x_min = min(xs) if xs else 0
    x_max = max(xs) if xs else 1
    if x_max == x_min:
        x_max = x_min + 1
    y_max = 100.0

    def x_scale(x_val: int) -> float:
        return margin["left"] + (x_val - x_min) / (x_max - x_min) * plot_w

    def y_scale(y_val: float) -> float:
        return margin["top"] + plot_h - (y_val / y_max) * plot_h

    parts = [
        f'<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
        f'xmlns="http://www.w3.org/2000/svg" role="img">',
        '<rect width="100%" height="100%" fill="white"/>',
    ]

    # horizontal grid lines
    for i in range(5):
        frac = i / 4
        y = margin["top"] + plot_h - frac * plot_h
        parts.append(
            f'<line x1="{margin["left"]}" y1="{y:.2f}" '
            f'x2="{margin["left"] + plot_w}" y2="{y:.2f}" '
            f'stroke="#e0e0e0" stroke-width="1"/>'
        )

    # axes
    parts.append(
        f'<line x1="{margin["left"]}" y1="{margin["top"] + plot_h}" '
        f'x2="{margin["left"] + plot_w}" y2="{margin["top"] + plot_h}" stroke="#333" stroke-width="1"/>'
    )
    parts.append(
        f'<line x1="{margin["left"]}" y1="{margin["top"]}" '
        f'x2="{margin["left"]}" y2="{margin["top"] + plot_h}" stroke="#333" stroke-width="1"/>'
    )

    # axis labels
    parts.append(
        f'<text x="{width / 2}" y="{height - 8}" text-anchor="middle" '
        f'font-size="12" font-family="sans-serif">Window midpoint (bp)</text>'
    )
    parts.append(
        f'<text x="16" y="{height / 2}" text-anchor="middle" font-size="12" '
        f'font-family="sans-serif" transform="rotate(-90 16 {height / 2})">% of window</text>'
    )

    # y ticks
    for i in range(5):
        frac = i / 4
        y = margin["top"] + plot_h - frac * plot_h
        val = frac * y_max
        parts.append(
            f'<line x1="{margin["left"] - 4}" y1="{y:.2f}" '
            f'x2="{margin["left"]}" y2="{y:.2f}" stroke="#333" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{margin["left"] - 8}" y="{y + 4:.2f}" text-anchor="end" '
            f'font-size="10" font-family="sans-serif">{val:.0f}</text>'
        )

    # x ticks
    for i in range(5):
        frac = i / 4
        x = margin["left"] + frac * plot_w
        val = int(round(x_min + frac * (x_max - x_min)))
        parts.append(
            f'<line x1="{x:.2f}" y1="{margin["top"] + plot_h}" '
            f'x2="{x:.2f}" y2="{margin["top"] + plot_h + 4}" stroke="#333" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{x:.2f}" y="{margin["top"] + plot_h + 18}" text-anchor="middle" '
            f'font-size="10" font-family="sans-serif">{val:,}</text>'
        )

    # series lines + legend entries
    for idx, (label, color, values) in enumerate(series):
        if not values:
            continue
        points = " ".join(f"{x_scale(x):.2f},{y_scale(y):.2f}" for x, y in zip(xs, values))
        parts.append(
            f'<polyline fill="none" stroke="{color}" stroke-width="2" points="{points}"/>'
        )
        legend_x = margin["left"] + 10 + idx * 170
        legend_y = margin["top"] - 16
        parts.append(
            f'<line x1="{legend_x}" y1="{legend_y}" x2="{legend_x + 18}" y2="{legend_y}" '
            f'stroke="{color}" stroke-width="2"/>'
        )
        parts.append(
            f'<text x="{legend_x + 24}" y="{legend_y + 4}" '
            f'font-size="11" font-family="sans-serif">{html.escape(label)}</text>'
        )

    parts.append("</svg>")
    return "\n".join(parts)


def svg_genome_overview(
    contigs: list[str],
    lengths: list[int],
    retained: list[int],
    masked: list[int],
    *,
    width: int = 900,
    height: int = 300,
) -> str:
    """Stacked vertical bar chart: retained / masked fraction per contig."""
    if not contigs:
        return ""
    margin = {"left": 50, "right": 150, "top": 20, "bottom": 65}
    plot_w = width - margin["left"] - margin["right"]
    plot_h = height - margin["top"] - margin["bottom"]
    n = len(contigs)
    bar_step = plot_w / n
    bar_w = bar_step * 0.75

    colors = {"retained": "#4C78A8", "masked": "#E45756"}

    parts = [
        f'<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
        f'xmlns="http://www.w3.org/2000/svg" role="img">',
        '<rect width="100%" height="100%" fill="white"/>',
    ]

    # horizontal grid lines + y ticks
    for i in range(5):
        frac = i / 4
        y = margin["top"] + plot_h - frac * plot_h
        pct_label = int(frac * 100)
        parts.append(
            f'<line x1="{margin["left"]}" y1="{y:.2f}" '
            f'x2="{margin["left"] + plot_w}" y2="{y:.2f}" stroke="#e0e0e0" stroke-width="1"/>'
        )
        parts.append(
            f'<line x1="{margin["left"] - 4}" y1="{y:.2f}" '
            f'x2="{margin["left"]}" y2="{y:.2f}" stroke="#333" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{margin["left"] - 8}" y="{y + 4:.2f}" text-anchor="end" '
            f'font-size="10" font-family="sans-serif">{pct_label}</text>'
        )

    # axes
    parts.append(
        f'<line x1="{margin["left"]}" y1="{margin["top"] + plot_h}" '
        f'x2="{margin["left"] + plot_w}" y2="{margin["top"] + plot_h}" stroke="#333" stroke-width="1"/>'
    )
    parts.append(
        f'<line x1="{margin["left"]}" y1="{margin["top"]}" '
        f'x2="{margin["left"]}" y2="{margin["top"] + plot_h}" stroke="#333" stroke-width="1"/>'
    )

    # y-axis label
    parts.append(
        f'<text x="14" y="{height / 2}" text-anchor="middle" font-size="12" '
        f'font-family="sans-serif" transform="rotate(-90 14 {height / 2})">% of contig</text>'
    )

    for i, (contig, length, ret, mask) in enumerate(zip(contigs, lengths, retained, masked)):
        bar_x = margin["left"] + i * bar_step + (bar_step - bar_w) / 2
        label_x = bar_x + bar_w / 2
        y_bottom = margin["top"] + plot_h

        def bar_h(val: int) -> float:
            return (val / length) * plot_h if length > 0 else 0.0

        ret_h = bar_h(ret)
        mask_h = bar_h(mask)

        parts.append(
            f'<rect x="{bar_x:.2f}" y="{y_bottom - ret_h:.2f}" '
            f'width="{bar_w:.2f}" height="{ret_h:.2f}" fill="{colors["retained"]}"/>'
        )
        parts.append(
            f'<rect x="{bar_x:.2f}" y="{y_bottom - ret_h - mask_h:.2f}" '
            f'width="{bar_w:.2f}" height="{mask_h:.2f}" fill="{colors["masked"]}"/>'
        )

        # rotated contig label
        parts.append(
            f'<text x="{label_x:.2f}" y="{y_bottom + 8}" text-anchor="end" '
            f'font-size="11" font-family="sans-serif" '
            f'transform="rotate(-40 {label_x:.2f} {y_bottom + 8})">'
            f'{html.escape(contig)}</text>'
        )

    # legend (right side)
    legend_x = margin["left"] + plot_w + 16
    for idx, (color, label) in enumerate([
        (colors["retained"], "Retained"),
        (colors["masked"], "Masked"),
    ]):
        ly = margin["top"] + 20 + idx * 24
        parts.append(f'<rect x="{legend_x}" y="{ly}" width="14" height="14" fill="{color}"/>')
        parts.append(
            f'<text x="{legend_x + 20}" y="{ly + 11}" '
            f'font-size="12" font-family="sans-serif">{label}</text>'
        )

    parts.append("</svg>")
    return "\n".join(parts)


def read_options_yaml(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except OSError as exc:
        return f"# Unable to read options file: {path}\n# {exc}"


def get_argprep_version() -> str:
    repo_root = Path(__file__).resolve().parents[1]
    try:
        result = subprocess.run(
            ["git", "-C", str(repo_root), "describe", "--tags", "--always", "--dirty"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return "unknown"
    version = result.stdout.strip()
    return version if version else "unknown"


def build_report(
    all_sites_paths: list[Path],
    mask_paths: list[Path],
    summary_paths: list[Path],
    fai: Path,
    window_bp: int,
    report_out: Path,
    options_yaml: Path,
    sample_missing_beds: list[Path] | None = None,
) -> None:
    lengths = read_fai_lengths(fai)
    invariant_counts: dict[str, list[int]] = {}
    variant_counts: dict[str, list[int]] = {}
    masked_counts: dict[str, list[int]] = {}
    summaries: dict[str, dict[str, str]] = {}
    sample_order: list[str] = []
    sample_variant_by_contig: dict[str, dict[str, int]] = {}
    sample_missing_retained_by_contig: dict[str, dict[str, int]] = {}
    sample_called_by_contig: dict[str, dict[str, int]] = {}

    for path in all_sites_paths:
        inv, var, samples, sv, sm, sc = read_all_sites_stats(path, lengths, window_bp)
        if samples and not sample_order:
            sample_order = list(samples)
        for contig, values in inv.items():
            if contig not in invariant_counts:
                invariant_counts[contig] = [0] * len(values)
                variant_counts[contig] = [0] * len(values)
            invariant_counts[contig] = [a + b for a, b in zip(invariant_counts[contig], values)]
            variant_counts[contig] = [a + b for a, b in zip(variant_counts[contig], var[contig])]
        for contig, d in sv.items():
            sample_variant_by_contig.setdefault(contig, {}).update(d)
        for contig, d in sm.items():
            sample_missing_retained_by_contig.setdefault(contig, {}).update(d)
        for contig, d in sc.items():
            sample_called_by_contig.setdefault(contig, {}).update(d)

    for path in mask_paths:
        masked = read_mask_percentages(path, lengths, window_bp)
        for contig, values in masked.items():
            if contig not in masked_counts:
                masked_counts[contig] = [0] * len(values)
            masked_counts[contig] = [a + b for a, b in zip(masked_counts[contig], values)]

    for path in summary_paths:
        values = read_site_summary(path)
        contig = values.get("contig")
        if contig:
            summaries[contig] = values

    active_contigs = [
        c for c in sorted(lengths)
        if any(invariant_counts.get(c, []))
        or any(variant_counts.get(c, []))
        or any(masked_counts.get(c, []))
        or c in summaries
    ]

    # genome-wide totals (from summary TSVs only)
    total_length = sum(lengths[c] for c in active_contigs)
    total_all_sites = sum(int(summaries[c].get("all_sites", 0)) for c in active_contigs if c in summaries)
    total_variants = sum(int(summaries[c].get("variants", 0)) for c in active_contigs if c in summaries)
    total_masked = sum(int(summaries[c].get("masked_total", 0)) for c in active_contigs if c in summaries)

    # per-sample (per-MAF) aggregation. BED-derived missing bp is optional; when
    # provided it restricts the reported sample set to those with a matching
    # ".{sample}.missing.bed" file (so synthesized columns like a REF track from
    # add_ref are excluded automatically).
    sample_bed_missing = read_sample_missing_bp(
        list(sample_missing_beds or []), sample_order, lengths
    ) if sample_missing_beds else {}
    if sample_missing_beds:
        per_maf_samples = [s for s in sample_order if s in sample_bed_missing]
    else:
        per_maf_samples = list(sample_order)

    CSS = """\
body{font-family:sans-serif;margin:24px;color:#111;max-width:1000px}
table{border-collapse:collapse;margin:12px 0}
th,td{border:1px solid #ccc;padding:4px 10px;text-align:right}
th{background:#f0f0f0;text-align:center}
td:first-child{text-align:left}
h1,h2,h3{margin-top:1.4em}
pre{background:#f6f8fa;border:1px solid #d0d7de;border-radius:6px;padding:12px;overflow:auto;white-space:pre-wrap}
details{margin:8px 0;border:1px solid #d0d7de;border-radius:6px;padding:4px 12px}
summary{cursor:pointer;font-size:1.05em;font-weight:bold;padding:6px 0;list-style:revert}
summary:hover{color:#0969da}
.toc{columns:4;column-gap:1em;margin:12px 0}
.toc a{display:block;color:#0969da;text-decoration:none;padding:1px 0}
.toc a:hover{text-decoration:underline}
tr.warn td{background:#fff8dc}
tr.bad  td{background:#ffd7d7}
"""

    report_out.parent.mkdir(parents=True, exist_ok=True)
    with report_out.open("w", encoding="utf-8") as fh:
        w = fh.write
        w('<!doctype html>\n<html lang="en">\n<head>\n')
        w('<meta charset="utf-8" />\n<meta name="viewport" content="width=device-width, initial-scale=1" />\n')
        w("<title>ARGprep summary</title>\n")
        w(f"<style>{CSS}</style>\n")
        w("</head>\n<body>\n")
        w("<h1>ARGprep Summary</h1>\n")
        w(
            f'<p>Generated: <code>{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</code>'
            f" &nbsp;|&nbsp; Window size: <code>{window_bp:,}</code> bp</p>\n"
        )

        # ── genome-wide overview ──────────────────────────────────────────────
        w("<h2>Genome-wide overview</h2>\n")
        w(svg_genome_overview(
            active_contigs,
            [lengths[c] for c in active_contigs],
            [int(summaries[c].get("all_sites", 0)) if c in summaries else 0 for c in active_contigs],
            [int(summaries[c].get("masked_total", 0)) if c in summaries else 0 for c in active_contigs],
        ))
        w("\n")

        w("<table>\n")
        w(
            "<tr><th>Contig</th><th>Length (bp)</th>"
            "<th>Retained sites</th><th>% retained</th>"
            "<th>Variants</th><th>% variant</th>"
            "<th>Masked</th><th>% masked</th></tr>\n"
        )
        for c in active_contigs:
            s = summaries.get(c, {})
            clen = lengths[c]
            ret = int(s.get("all_sites", 0))
            var = int(s.get("variants", 0))
            mask = int(s.get("masked_total", 0))
            ret_pct = ret / clen * 100 if clen else 0.0
            var_pct = var / ret * 100 if ret else 0.0
            mask_pct = mask / clen * 100 if clen else 0.0
            row_cls = " class=\"bad\"" if mask_pct > 50 else (" class=\"warn\"" if mask_pct > 25 else "")
            w(
                f'<tr{row_cls}>'
                f'<td><a href="#{_anchor(c)}">{html.escape(c)}</a></td>'
                f"<td>{clen:,}</td>"
                f"<td>{ret:,}</td>"
                f"<td>{ret_pct:.1f}%</td>"
                f"<td>{var:,}</td>"
                f"<td>{var_pct:.1f}%</td>"
                f"<td>{mask:,}</td>"
                f"<td>{mask_pct:.1f}%</td>"
                f"</tr>\n"
            )
        if active_contigs:
            ret_pct_all = total_all_sites / total_length * 100 if total_length else 0.0
            var_pct_all = total_variants / total_all_sites * 100 if total_all_sites else 0.0
            mask_pct_all = total_masked / total_length * 100 if total_length else 0.0
            row_cls = " class=\"bad\"" if mask_pct_all > 50 else (" class=\"warn\"" if mask_pct_all > 25 else "")
            w(
                f'<tr{row_cls} style="font-weight:bold;border-top:2px solid #999">'
                f"<td>Total</td>"
                f"<td>{total_length:,}</td>"
                f"<td>{total_all_sites:,}</td>"
                f"<td>{ret_pct_all:.1f}%</td>"
                f"<td>{total_variants:,}</td>"
                f"<td>{var_pct_all:.1f}%</td>"
                f"<td>{total_masked:,}</td>"
                f"<td>{mask_pct_all:.1f}%</td>"
                f"</tr>\n"
            )
        w("</table>\n")

        # ── per-MAF (per-sample) summary ─────────────────────────────────────
        if per_maf_samples:
            w("<h2>Per-MAF summary</h2>\n")
            w(
                "<p>One row per input MAF / sample. Missing bp is computed across "
                "the full contig (from the per-sample <code>.missing.bed</code> masks). "
                "Variant and called counts are across retained sites in the "
                "<code>all_sites</code> VCFs.</p>\n"
            )
            w("<table>\n")
            w(
                "<tr><th>Sample</th>"
                "<th>Missing bp</th><th>% genome missing</th>"
                "<th>Retained sites called</th>"
                "<th>Variants carried</th><th>% variant of called</th>"
                "</tr>\n"
            )
            totals_missing = 0
            totals_called = 0
            totals_variant = 0
            for sample in per_maf_samples:
                miss_bp = (
                    sum(sample_bed_missing.get(sample, {}).get(c, 0) for c in active_contigs)
                    if sample_bed_missing else 0
                )
                called = sum(
                    sample_called_by_contig.get(c, {}).get(sample, 0) for c in active_contigs
                )
                carried = sum(
                    sample_variant_by_contig.get(c, {}).get(sample, 0) for c in active_contigs
                )
                miss_pct = miss_bp / total_length * 100 if total_length else 0.0
                var_of_called = carried / called * 100 if called else 0.0
                row_cls = (
                    ' class="bad"' if miss_pct > 50
                    else (' class="warn"' if miss_pct > 25 else "")
                )
                miss_cell = f"{miss_bp:,}" if sample_bed_missing else "&mdash;"
                miss_pct_cell = f"{miss_pct:.1f}%" if sample_bed_missing else "&mdash;"
                w(
                    f"<tr{row_cls}>"
                    f"<td>{html.escape(sample)}</td>"
                    f"<td>{miss_cell}</td>"
                    f"<td>{miss_pct_cell}</td>"
                    f"<td>{called:,}</td>"
                    f"<td>{carried:,}</td>"
                    f"<td>{var_of_called:.1f}%</td>"
                    f"</tr>\n"
                )
                totals_missing += miss_bp
                totals_called += called
                totals_variant += carried
            if len(per_maf_samples) > 1:
                n = len(per_maf_samples)
                mean_missing = totals_missing / n
                mean_called = totals_called / n
                mean_variant = totals_variant / n
                mean_miss_pct = mean_missing / total_length * 100 if total_length else 0.0
                mean_var_pct = mean_variant / mean_called * 100 if mean_called else 0.0
                mean_miss_cell = f"{mean_missing:,.0f}" if sample_bed_missing else "&mdash;"
                mean_miss_pct_cell = f"{mean_miss_pct:.1f}%" if sample_bed_missing else "&mdash;"
                w(
                    f'<tr style="font-weight:bold;border-top:2px solid #999">'
                    f"<td>Mean</td>"
                    f"<td>{mean_miss_cell}</td>"
                    f"<td>{mean_miss_pct_cell}</td>"
                    f"<td>{mean_called:,.0f}</td>"
                    f"<td>{mean_variant:,.0f}</td>"
                    f"<td>{mean_var_pct:.1f}%</td>"
                    f"</tr>\n"
                )
            w("</table>\n")

        # ── table of contents ─────────────────────────────────────────────────
        w("<h2>Contigs</h2>\n<div class=\"toc\">\n")
        for c in active_contigs:
            clen = lengths[c]
            s = summaries.get(c, {})
            ret = int(s.get("all_sites", 0))
            ret_pct = ret / clen * 100 if clen else 0.0
            w(f'<a href="#{_anchor(c)}">{html.escape(c)} ({ret_pct:.0f}%)</a>\n')
        w("</div>\n")

        # ── per-contig sections ───────────────────────────────────────────────
        for contig in active_contigs:
            length = lengths[contig]
            xs = window_midpoints(length, window_bp)
            inv_raw = invariant_counts.get(contig, [0] * len(xs))
            var_raw = variant_counts.get(contig, [0] * len(xs))
            miss_raw = masked_counts.get(contig, [0] * len(xs))
            inv_pct = to_percentages(inv_raw, length, window_bp)
            var_pct = to_percentages(var_raw, length, window_bp)
            miss_pct = to_percentages(miss_raw, length, window_bp)

            summary = summaries.get(contig, {})
            ret = int(summary.get("all_sites", 0))
            mask = int(summary.get("masked_total", 0))
            ret_pct_val = ret / length * 100 if length else 0.0
            mask_pct_val = mask / length * 100 if length else 0.0

            w(f'<details id="{_anchor(contig)}">\n')
            w(
                f"<summary>{html.escape(contig)}"
                f" &mdash; {length:,} bp"
                f" &nbsp;|&nbsp; {ret_pct_val:.1f}% retained"
                f" &nbsp;|&nbsp; {mask_pct_val:.1f}% masked"
                f"</summary>\n"
            )

            if summary:
                w("<table>\n<tr><th>Metric</th><th>Value</th><th>% of contig</th></tr>\n")
                for key, label in (
                    ("contig_length", "Length (bp)"),
                    ("samples", "Samples"),
                    ("allowed_missing", "Allowed missing"),
                    ("all_sites", "Retained sites"),
                    ("variants", "Variants"),
                    ("invariant", "Invariant"),
                    ("masked_total", "Masked total"),
                    ("masked_missingness", "Masked — missingness"),
                    ("masked_indel_adjacent", "Masked — indel-adjacent SNP"),
                    ("masked_multiallelic", "Masked — multiallelic"),
                    ("masked_no_alignment", "Masked — no alignment"),
                    ("masked_ref_non_acgt", "Masked — non-ACGT ref"),
                ):
                    if key not in summary:
                        continue
                    raw_val = summary[key]
                    fmt_val = fmt_int(raw_val)
                    pct_cell = ""
                    if key not in ("samples", "allowed_missing") and length:
                        try:
                            pct_cell = f"{int(raw_val) / length * 100:.1f}%"
                        except (ValueError, TypeError):
                            pass
                    row_cls = ""
                    if key == "masked_total":
                        try:
                            mp = int(raw_val) / length * 100 if length else 0.0
                            row_cls = ' class="bad"' if mp > 50 else (' class="warn"' if mp > 25 else "")
                        except (ValueError, TypeError):
                            pass
                    w(
                        f"<tr{row_cls}>"
                        f"<td>{html.escape(label)}</td>"
                        f"<td>{html.escape(fmt_val)}</td>"
                        f"<td>{html.escape(pct_cell)}</td>"
                        f"</tr>\n"
                    )
                w("</table>\n")

            w(svg_combined_plot(
                xs,
                [
                    ("Invariant (%)", "#4C78A8", inv_pct),
                    ("Variable (%)", "#F58518", var_pct),
                    ("Missing (%)", "#E45756", miss_pct),
                ],
            ))

            # per-MAF breakdown for this contig
            if per_maf_samples:
                contig_miss_map = {
                    s: sample_bed_missing.get(s, {}).get(contig, 0)
                    for s in per_maf_samples
                } if sample_bed_missing else {}
                contig_called_map = sample_called_by_contig.get(contig, {})
                contig_var_map = sample_variant_by_contig.get(contig, {})
                w("<h4>Per-MAF on this contig</h4>\n")
                w("<table>\n")
                w(
                    "<tr><th>Sample</th>"
                    "<th>Missing bp</th><th>% contig missing</th>"
                    "<th>Retained called</th><th>Variants carried</th>"
                    "</tr>\n"
                )
                for sample in per_maf_samples:
                    miss_bp = contig_miss_map.get(sample, 0)
                    called = contig_called_map.get(sample, 0)
                    carried = contig_var_map.get(sample, 0)
                    miss_pct = miss_bp / length * 100 if length else 0.0
                    row_cls = (
                        ' class="bad"' if miss_pct > 50
                        else (' class="warn"' if miss_pct > 25 else "")
                    )
                    miss_cell = f"{miss_bp:,}" if sample_bed_missing else "&mdash;"
                    miss_pct_cell = f"{miss_pct:.1f}%" if sample_bed_missing else "&mdash;"
                    w(
                        f"<tr{row_cls}>"
                        f"<td>{html.escape(sample)}</td>"
                        f"<td>{miss_cell}</td>"
                        f"<td>{miss_pct_cell}</td>"
                        f"<td>{called:,}</td>"
                        f"<td>{carried:,}</td>"
                        f"</tr>\n"
                    )
                w("</table>\n")

            w("\n</details>\n")

        # ── run configuration ─────────────────────────────────────────────────
        w("<h2>Run configuration</h2>\n")
        w(f"<p>ARGprep version: <code>{html.escape(get_argprep_version())}</code></p>\n")
        w(f"<p>Source options file: <code>{html.escape(str(options_yaml))}</code></p>\n")
        w("<pre><code>" + html.escape(read_options_yaml(options_yaml)) + "</code></pre>\n")
        w("</body>\n</html>\n")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Build a direct-pipeline HTML summary report.")
    ap.add_argument("--fai", required=True)
    ap.add_argument("--window-bp", type=int, default=100000)
    ap.add_argument("--report-out", required=True)
    ap.add_argument("--all-sites", nargs="+", required=True)
    ap.add_argument("--masked-beds", nargs="+", required=True)
    ap.add_argument("--site-summaries", nargs="+", required=True)
    ap.add_argument("--options-yaml", required=True)
    ap.add_argument(
        "--sample-missing-beds",
        nargs="*",
        default=[],
        help="Per-sample missing BED files (combined.{contig}.{sample}.missing.bed). "
        "Used to populate the per-MAF missingness table.",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    build_report(
        [Path(p) for p in args.all_sites],
        [Path(p) for p in args.masked_beds],
        [Path(p) for p in args.site_summaries],
        Path(args.fai),
        args.window_bp,
        Path(args.report_out),
        Path(args.options_yaml),
        [Path(p) for p in args.sample_missing_beds],
    )


if __name__ == "__main__":
    main()
