#!/usr/bin/env python3
"""Remake the chromosome-1 line plot from admix_results/summary.html showing
only the Variable (%) series, with an auto-scaled y-axis.

The variable-site percentages are recovered directly from the orange polyline
already rendered in summary.html (no need to rescan the multi-GB all_sites VCF).
The original plot locks the y-axis to 0-100%, which flattens the variable line
near the bottom; here we auto-scale so the signal is visible.
"""
from __future__ import annotations

import html
import math
import re
from pathlib import Path

# --- geometry of the original svg_combined_plot (see scripts/summary_report.py)
ORIG = dict(width=900, height=300, left=60, right=20, top=28, bottom=55, y_max=100.0)
ORIG_PLOT_H = ORIG["height"] - ORIG["top"] - ORIG["bottom"]  # 217

WINDOW_BP = 100_000
CHR1_LEN = 308_452_471

HERE = Path(__file__).resolve().parents[1] / "admix_results"
SUMMARY = HERE / "summary.html"
OUT = HERE / "chr1_variable.html"


def extract_chr1_variable_pcts() -> list[float]:
    text = SUMMARY.read_text(encoding="utf-8")
    # isolate the chromosome-1 <details id="c-1"> ... block
    start = text.index('<details id="c-1">')
    end = text.index('<details id="c-2">')
    block = text[start:end]
    m = re.search(
        r'<polyline fill="none" stroke="#F58518" stroke-width="2" points="([^"]*)"',
        block,
    )
    if not m:
        raise SystemExit("Could not find the chr1 Variable (#F58518) polyline")
    pts = m.group(1).split()
    # invert y_scale: y_px = top + plot_h - (val/y_max)*plot_h
    vals = []
    for p in pts:
        _, y_str = p.split(",")
        y_px = float(y_str)
        val = (ORIG["top"] + ORIG_PLOT_H - y_px) / ORIG_PLOT_H * ORIG["y_max"]
        vals.append(max(val, 0.0))
    return vals


def window_midpoints(length: int, window_bp: int) -> list[int]:
    mids = []
    n = (length + window_bp - 1) // window_bp
    for idx in range(n):
        s = idx * window_bp
        e = min(s + window_bp, length)
        mids.append(s + (e - s) // 2)
    return mids


def nice_ceiling(v: float) -> float:
    if v <= 0:
        return 1.0
    exp = math.floor(math.log10(v))
    base = 10 ** exp
    for mult in (1, 1.5, 2, 2.5, 3, 4, 5, 7.5, 10):
        if mult * base >= v:
            return mult * base
    return 10 * base


def render(xs: list[int], values: list[float], *, width=900, height=320) -> str:
    margin = {"left": 60, "right": 20, "top": 28, "bottom": 55}
    plot_w = width - margin["left"] - margin["right"]
    plot_h = height - margin["top"] - margin["bottom"]
    x_min, x_max = min(xs), max(xs)
    if x_max == x_min:
        x_max = x_min + 1
    y_max = nice_ceiling(max(values) if values else 1.0)

    def xs_(x):
        return margin["left"] + (x - x_min) / (x_max - x_min) * plot_w

    def ys_(y):
        return margin["top"] + plot_h - (y / y_max) * plot_h

    parts = [
        f'<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
        f'xmlns="http://www.w3.org/2000/svg" role="img">',
        '<rect width="100%" height="100%" fill="white"/>',
    ]
    # grid + y ticks
    for i in range(5):
        frac = i / 4
        y = margin["top"] + plot_h - frac * plot_h
        val = frac * y_max
        parts.append(
            f'<line x1="{margin["left"]}" y1="{y:.2f}" x2="{margin["left"] + plot_w}" '
            f'y2="{y:.2f}" stroke="#e0e0e0" stroke-width="1"/>'
        )
        parts.append(
            f'<line x1="{margin["left"] - 4}" y1="{y:.2f}" x2="{margin["left"]}" '
            f'y2="{y:.2f}" stroke="#333" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{margin["left"] - 8}" y="{y + 4:.2f}" text-anchor="end" '
            f'font-size="10" font-family="sans-serif">{val:.2f}</text>'
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
        f'font-family="sans-serif" transform="rotate(-90 16 {height / 2})">Variable sites (% of window)</text>'
    )
    # x ticks
    for i in range(5):
        frac = i / 4
        x = margin["left"] + frac * plot_w
        val = int(round(x_min + frac * (x_max - x_min)))
        parts.append(
            f'<line x1="{x:.2f}" y1="{margin["top"] + plot_h}" x2="{x:.2f}" '
            f'y2="{margin["top"] + plot_h + 4}" stroke="#333" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{x:.2f}" y="{margin["top"] + plot_h + 18}" text-anchor="middle" '
            f'font-size="10" font-family="sans-serif">{val:,}</text>'
        )
    # legend
    parts.append(
        f'<line x1="{margin["left"] + 10}" y1="{margin["top"] - 16}" '
        f'x2="{margin["left"] + 28}" y2="{margin["top"] - 16}" stroke="#F58518" stroke-width="2"/>'
    )
    parts.append(
        f'<text x="{margin["left"] + 34}" y="{margin["top"] - 12}" '
        f'font-size="11" font-family="sans-serif">Variable (%) &mdash; chromosome 1</text>'
    )
    # data line
    points = " ".join(f"{xs_(x):.2f},{ys_(y):.2f}" for x, y in zip(xs, values))
    parts.append(f'<polyline fill="none" stroke="#F58518" stroke-width="2" points="{points}"/>')
    parts.append("</svg>")
    return "\n".join(parts)


def main() -> None:
    values = extract_chr1_variable_pcts()
    xs = window_midpoints(CHR1_LEN, WINDOW_BP)
    if len(xs) != len(values):
        # the polyline has one point per window; align defensively
        n = min(len(xs), len(values))
        xs, values = xs[:n], values[:n]
    svg = render(xs, values)
    doc = (
        '<!doctype html>\n<html lang="en">\n<head>\n'
        '<meta charset="utf-8" />\n'
        "<title>ARGprep — chr1 variable sites</title>\n"
        "<style>body{font-family:sans-serif;margin:24px;color:#111;max-width:1000px}</style>\n"
        "</head>\n<body>\n"
        "<h1>Chromosome 1 — variable sites</h1>\n"
        f"<p>Window size: <code>{WINDOW_BP:,}</code> bp &nbsp;|&nbsp; "
        f"peak window: <code>{max(values):.3f}%</code> &nbsp;|&nbsp; "
        f"mean: <code>{sum(values) / len(values):.3f}%</code></p>\n"
        "<p>Variable-site percentage per 100 kb window, auto-scaled y-axis "
        "(the combined plot in <code>summary.html</code> is locked to 0–100%, "
        "which flattens this series).</p>\n"
        + svg
        + "\n</body>\n</html>\n"
    )
    OUT.write_text(doc, encoding="utf-8")
    print(f"wrote {OUT}  ({len(values)} windows)")


if __name__ == "__main__":
    main()
