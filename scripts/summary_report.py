#!/usr/bin/env python3
from __future__ import annotations

import argparse
import html
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


def read_all_sites_percentages(
    path: Path,
    lengths: dict[str, int],
    window_bp: int,
) -> tuple[dict[str, list[int]], dict[str, list[int]]]:
    invariant: dict[str, list[int]] = {
        contig: [0] * window_count(length, window_bp) for contig, length in lengths.items()
    }
    variant: dict[str, list[int]] = {
        contig: [0] * window_count(length, window_bp) for contig, length in lengths.items()
    }
    with open_text(path, "rt", errors="ignore") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            contig = parts[0]
            if contig not in invariant:
                continue
            pos = int(parts[1])
            idx = max((pos - 1) // window_bp, 0)
            if parts[4] == ".":
                invariant[contig][idx] += 1
            else:
                variant[contig][idx] += 1
    return invariant, variant


def read_mask_percentages(
    path: Path,
    lengths: dict[str, int],
    window_bp: int,
) -> dict[str, list[int]]:
    masked: dict[str, list[int]] = {
        contig: [0] * window_count(length, window_bp) for contig, length in lengths.items()
    }
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            contig = parts[0]
            if contig not in masked:
                continue
            start = int(parts[1])
            end = int(parts[2])
            pos = start
            while pos < end:
                idx = pos // window_bp
                window_end = min((idx + 1) * window_bp, lengths[contig])
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


def svg_series_plot(
    xs: list[int],
    label: str,
    color: str,
    values: list[float],
    *,
    width: int = 900,
    height: int = 280,
) -> str:
    margin = {"left": 60, "right": 20, "top": 20, "bottom": 55}
    plot_w = width - margin["left"] - margin["right"]
    plot_h = height - margin["top"] - margin["bottom"]
    x_min = min(xs) if xs else 0
    x_max = max(xs) if xs else 1
    if x_max == x_min:
        x_max = x_min + 1
    y_max = max(values, default=1.0)
    y_max = max(y_max, 1.0)

    def x_scale(x_val: int) -> float:
        return margin["left"] + (x_val - x_min) / (x_max - x_min) * plot_w

    def y_scale(y_val: float) -> float:
        return margin["top"] + plot_h - (y_val / y_max) * plot_h

    parts = [
        f'<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" role="img">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<line x1="{margin["left"]}" y1="{margin["top"] + plot_h}" x2="{margin["left"] + plot_w}" y2="{margin["top"] + plot_h}" stroke="#333" stroke-width="1"/>',
        f'<line x1="{margin["left"]}" y1="{margin["top"]}" x2="{margin["left"]}" y2="{margin["top"] + plot_h}" stroke="#333" stroke-width="1"/>',
        f'<text x="{width / 2}" y="{height - 10}" text-anchor="middle" font-size="12" font-family="sans-serif">Window midpoint (bp)</text>',
        f'<text x="18" y="{height / 2}" text-anchor="middle" font-size="12" font-family="sans-serif" transform="rotate(-90 18 {height / 2})">Percent of window</text>',
    ]

    for i in range(5):
        frac = i / 4
        y = margin["top"] + plot_h - frac * plot_h
        val = frac * y_max
        parts.append(f'<line x1="{margin["left"] - 4}" y1="{y:.2f}" x2="{margin["left"]}" y2="{y:.2f}" stroke="#333" stroke-width="1"/>')
        parts.append(f'<text x="{margin["left"] - 8}" y="{y + 4:.2f}" text-anchor="end" font-size="10" font-family="sans-serif">{val:.1f}</text>')

    for i in range(5):
        frac = i / 4
        x = margin["left"] + frac * plot_w
        val = int(round(x_min + frac * (x_max - x_min)))
        parts.append(f'<line x1="{x:.2f}" y1="{margin["top"] + plot_h}" x2="{x:.2f}" y2="{margin["top"] + plot_h + 4}" stroke="#333" stroke-width="1"/>')
        parts.append(f'<text x="{x:.2f}" y="{margin["top"] + plot_h + 18}" text-anchor="middle" font-size="10" font-family="sans-serif">{val:,}</text>')

    legend_x = margin["left"] + 10
    legend_y = margin["top"] + 8
    points = " ".join(
        f"{x_scale(x):.2f},{y_scale(y):.2f}" for x, y in zip(xs, values)
    )
    parts.append(f'<polyline fill="none" stroke="{color}" stroke-width="2" points="{points}"/>')
    for x, y in zip(xs, values):
        parts.append(f'<circle cx="{x_scale(x):.2f}" cy="{y_scale(y):.2f}" r="2.5" fill="{color}"/>')
    parts.append(f'<line x1="{legend_x}" y1="{legend_y}" x2="{legend_x + 18}" y2="{legend_y}" stroke="{color}" stroke-width="2"/>')
    parts.append(f'<text x="{legend_x + 24}" y="{legend_y + 4}" font-size="11" font-family="sans-serif">{html.escape(label)}</text>')

    parts.append("</svg>")
    return "\n".join(parts)


def build_report(
    all_sites_paths: list[Path],
    mask_paths: list[Path],
    summary_paths: list[Path],
    fai: Path,
    window_bp: int,
    report_out: Path,
) -> None:
    lengths = read_fai_lengths(fai)
    invariant_counts: dict[str, list[int]] = {}
    variant_counts: dict[str, list[int]] = {}
    masked_counts: dict[str, list[int]] = {}
    summaries: dict[str, dict[str, str]] = {}

    for path in all_sites_paths:
        inv, var = read_all_sites_percentages(path, lengths, window_bp)
        for contig, values in inv.items():
            if contig not in invariant_counts:
                invariant_counts[contig] = [0] * len(values)
                variant_counts[contig] = [0] * len(values)
            invariant_counts[contig] = [a + b for a, b in zip(invariant_counts[contig], values)]
            variant_counts[contig] = [a + b for a, b in zip(variant_counts[contig], var[contig])]

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

    report_out.parent.mkdir(parents=True, exist_ok=True)
    with report_out.open("w", encoding="utf-8") as handle:
        handle.write("<!doctype html>\n<html lang=\"en\">\n<head>\n")
        handle.write('<meta charset="utf-8" />\n<meta name="viewport" content="width=device-width, initial-scale=1" />\n')
        handle.write("<title>ARGprep summary</title>\n")
        handle.write("<style>body{font-family:sans-serif;margin:24px;color:#111}table{border-collapse:collapse;margin:12px 0}th,td{border:1px solid #ccc;padding:4px 8px}h1,h2,h3{margin-top:1.4em}</style>\n")
        handle.write("</head>\n<body>\n")
        handle.write("<h1>ARGprep Summary</h1>\n")
        handle.write(f"<p>Window size: <code>{window_bp:,}</code> bp</p>\n")
        for contig in sorted(lengths):
            length = lengths[contig]
            xs = window_midpoints(length, window_bp)
            inv_pct = to_percentages(invariant_counts.get(contig, [0] * len(xs)), length, window_bp)
            var_pct = to_percentages(variant_counts.get(contig, [0] * len(xs)), length, window_bp)
            miss_pct = to_percentages(masked_counts.get(contig, [0] * len(xs)), length, window_bp)
            handle.write(f"<h2>{html.escape(contig)}</h2>\n")
            summary = summaries.get(contig)
            if summary:
                handle.write("<table>\n<tr><th>Metric</th><th>Value</th></tr>\n")
                for key in (
                    "contig_length",
                    "samples",
                    "allowed_missing",
                    "all_sites",
                    "variants",
                    "invariant",
                    "masked_total",
                    "masked_missingness",
                    "masked_indel",
                    "masked_multiallelic",
                    "masked_no_alignment",
                    "masked_ref_non_acgt",
                ):
                    if key in summary:
                        handle.write(
                            f"<tr><td>{html.escape(key)}</td><td>{html.escape(summary[key])}</td></tr>\n"
                        )
                handle.write("</table>\n")
            for label, color, values in (
                ("Invariant (%)", "#4C78A8", inv_pct),
                ("Variable (%)", "#F58518", var_pct),
                ("Missing (%)", "#E45756", miss_pct),
            ):
                handle.write(f"<h3>{html.escape(label)}</h3>\n")
                handle.write(svg_series_plot(xs, label, color, values))
        handle.write("</body>\n</html>\n")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Build a direct-pipeline HTML summary report.")
    ap.add_argument("--fai", required=True)
    ap.add_argument("--window-bp", type=int, default=100000)
    ap.add_argument("--report-out", required=True)
    ap.add_argument("--all-sites", nargs="+", required=True)
    ap.add_argument("--masked-beds", nargs="+", required=True)
    ap.add_argument("--site-summaries", nargs="+", required=True)
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
    )


if __name__ == "__main__":
    main()
