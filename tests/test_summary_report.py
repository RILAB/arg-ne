import re
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.summary_report import nice_axis_max, svg_combined_plot  # noqa: E402


@pytest.mark.parametrize(
    ("peak", "expected"),
    [
        (0.0, 1.0),      # no data at all still yields a usable axis
        (0.0031, 0.005),
        (0.08, 0.1),
        (0.23, 0.25),
        (1.4, 2.0),
        (4.7, 5.0),
        (47.0, 50.0),
        (99.9, 100.0),
        (140.0, 100.0),  # a percentage axis never exceeds 100
    ],
)
def test_nice_axis_max_rounds_up_to_readable_bound(peak: float, expected: float):
    assert nice_axis_max(peak) == expected
    assert nice_axis_max(peak) >= min(peak, 100.0)


@pytest.mark.parametrize(
    ("y_max", "expected"),
    [
        (100.0, ["0", "25", "50", "75", "100"]),
        (5.0, ["0.00", "1.25", "2.50", "3.75", "5.00"]),      # not 0,1,2,4,5
        (2.0, ["0.0", "0.5", "1.0", "1.5", "2.0"]),
        (0.05, ["0.0000", "0.0125", "0.0250", "0.0375", "0.0500"]),
    ],
)
def test_axis_ticks_use_enough_decimals_to_be_exact(y_max: float, expected: list[str]):
    """Quartered axes give steps like 1.25; rounding those to "1" misreads the plot."""
    svg = svg_combined_plot([0, 1], [("Variable (%)", "#F58518", [0.0, y_max])], y_max=y_max)
    ticks = re.findall(r'text-anchor="end"[^>]*>([^<]+)</text>', svg)
    assert ticks == expected


def test_series_above_axis_max_is_clamped_into_the_plot():
    """Out-of-range values must not draw outside the axes."""
    svg = svg_combined_plot([0, 1], [("Variable (%)", "#F58518", [0.0, 99.0])], y_max=1.0)
    points = svg.split('points="')[1].split('"')[0]
    ys = [float(pair.split(",")[1]) for pair in points.split()]
    # y=1.0 (the axis max) maps to the top margin; nothing may sit above it.
    assert min(ys) == pytest.approx(28.0)


def test_summary_report_html_contains_separate_plots_per_contig(tmp_path: Path):
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t10\t0\t0\t0\nchr2\t8\t0\t0\t0\nscaffold99\t500\t0\t0\t0\n", encoding="utf-8")

    report_stats = tmp_path / "combined.report_stats.tsv"
    report_stats.write_text(
        "record_type\tcontig\tstart\tend\tsample\tinvariant\tvariant\tmasked\tcalled\tcarried_variant\n"
        "window\tchr1\t0\t4\t\t1\t1\t2\t\t\n"
        "window\tchr1\t4\t8\t\t1\t0\t3\t\t\n"
        "window\tchr1\t8\t10\t\t0\t0\t2\t\t\n"
        "sample\tchr1\t\t\ts1\t\t\t\t3\t1\n"
        "window\tchr2\t0\t4\t\t1\t1\t2\t\t\n"
        "window\tchr2\t4\t8\t\t0\t0\t4\t\t\n"
        "sample\tchr2\t\t\ts1\t\t\t\t2\t1\n",
        encoding="utf-8",
    )

    summary1 = tmp_path / "combined.chr1.site_summary.tsv"
    summary1.write_text(
        "metric\tvalue\n"
        "contig\tchr1\n"
        "contig_length\t10\n"
        "samples\t1\n"
        "allowed_missing\t0\n"
        "all_sites\t3\n"
        "variants\t1\n"
        "invariant\t2\n"
        "masked_total\t7\n",
        encoding="utf-8",
    )

    summary2 = tmp_path / "combined.chr2.site_summary.tsv"
    summary2.write_text(
        "metric\tvalue\n"
        "contig\tchr2\n"
        "contig_length\t8\n"
        "samples\t1\n"
        "allowed_missing\t0\n"
        "all_sites\t2\n"
        "variants\t1\n"
        "invariant\t1\n"
        "masked_total\t6\n",
        encoding="utf-8",
    )

    options_yaml = tmp_path / "options.yaml"
    options_yaml.write_text(
        "maf_dir: /tmp/maf\nreference_fasta: /tmp/ref.fa\nadd_ref: true\n",
        encoding="utf-8",
    )

    report = tmp_path / "summary.html"
    subprocess.run(
        [
            sys.executable,
            str(Path("scripts") / "summary_report.py"),
            "--fai",
            str(fai),
            "--window-bp",
            "4",
            "--report-out",
            str(report),
            "--report-stats",
            str(report_stats),
            "--site-summaries",
            str(summary1),
            str(summary2),
            "--options-yaml",
            str(options_yaml),
        ],
        cwd=Path.cwd(),
        check=True,
    )

    html = report.read_text(encoding="utf-8")
    assert "ARGprep Summary" in html
    assert "chr1" in html
    assert "chr2" in html
    assert "Invariant (%)" in html
    assert "Variable (%)" in html
    assert "Missing (%)" in html
    assert "Unaligned" not in html
    # Two SVGs per active contig (chr1, chr2) — no per-series <h3> headings:
    # a fixed 0-100% invariant/missing plot, and a separately scaled variable
    # plot, so the variable signal is not flattened against the 0-100% axis.
    assert html.count("<details") == 2
    assert html.count("<polyline") >= 6  # 3 series × 2 contigs
    assert html.count("% of window") == 2  # one per contig
    assert html.count("% variable of window") == 2  # one per contig
    assert "ARGprep version:" in html
    assert "Source options file:" in html
    assert "maf_dir: /tmp/maf" in html
    # Scaffolds with no data should be omitted.
    assert "scaffold99" not in html


def test_summary_report_per_maf_section(tmp_path: Path):
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t10\t0\t0\t0\n", encoding="utf-8")

    report_stats = tmp_path / "combined.chr1.report_stats.tsv"
    report_stats.write_text(
        "record_type\tcontig\tstart\tend\tsample\tinvariant\tvariant\tmasked\tcalled\tcarried_variant\n"
        "window\tchr1\t0\t4\t\t1\t2\t1\t\t\n"
        "window\tchr1\t4\t8\t\t1\t0\t3\t\t\n"
        "window\tchr1\t8\t10\t\t0\t0\t2\t\t\n"
        "sample\tchr1\t\t\ts1\t\t\t\t4\t2\n"
        "sample\tchr1\t\t\ts2\t\t\t\t3\t1\n",
        encoding="utf-8",
    )

    summary1 = tmp_path / "combined.chr1.site_summary.tsv"
    summary1.write_text(
        "metric\tvalue\n"
        "contig\tchr1\n"
        "contig_length\t10\n"
        "samples\t2\n"
        "allowed_missing\t1\n"
        "all_sites\t4\n"
        "variants\t2\n"
        "invariant\t2\n"
        "masked_total\t4\n",
        encoding="utf-8",
    )

    # s1 missing 2 bp on chr1, s2 missing 5 bp on chr1
    bed_s1 = tmp_path / "combined.chr1.s1.missing.bed"
    bed_s1.write_text("chr1\t8\t10\n", encoding="utf-8")
    bed_s2 = tmp_path / "combined.chr1.s2.missing.bed"
    bed_s2.write_text("chr1\t0\t5\n", encoding="utf-8")

    options_yaml = tmp_path / "options.yaml"
    options_yaml.write_text("maf_dir: /tmp/maf\nreference_fasta: /tmp/ref.fa\n", encoding="utf-8")

    report = tmp_path / "summary.html"
    subprocess.run(
        [
            sys.executable,
            str(Path("scripts") / "summary_report.py"),
            "--fai", str(fai),
            "--window-bp", "4",
            "--report-out", str(report),
            "--report-stats", str(report_stats),
            "--site-summaries", str(summary1),
            "--sample-missing-beds", str(bed_s1), str(bed_s2),
            "--options-yaml", str(options_yaml),
        ],
        cwd=Path.cwd(),
        check=True,
    )

    html = report.read_text(encoding="utf-8")
    assert "Per-MAF summary" in html
    assert "Per-MAF on this contig" in html
    # s1: 2 missing bp, 2 variants carried (positions 2 and 3), 4 called sites
    # s2: 5 missing bp, 1 variant carried (position 3), 3 called sites (pos 5 is missing in retained)
    assert "<td>s1</td><td>2</td><td>20.0%</td><td>4</td><td>2</td>" in html
    assert "<td>s2</td><td>5</td><td>50.0%</td><td>3</td><td>1</td>" in html
    assert "<td>REF</td>" not in html
