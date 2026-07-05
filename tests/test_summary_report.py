import subprocess
import sys
from pathlib import Path


def test_summary_report_html_contains_separate_plots_per_contig(tmp_path: Path):
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t10\t0\t0\t0\nchr2\t8\t0\t0\t0\nscaffold99\t500\t0\t0\t0\n", encoding="utf-8")

    all_sites = tmp_path / "combined.all_sites.vcf"
    all_sites.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
        "chr1\t1\t.\tA\t.\t.\tPASS\tNS=1;MS=0;SC=invariant\tGT\t0\n"
        "chr1\t2\t.\tC\tT\t.\tPASS\tNS=1;MS=0;SC=variant\tGT\t1\n"
        "chr1\t5\t.\tG\t.\t.\tPASS\tNS=1;MS=0;SC=invariant\tGT\t0\n"
        "chr2\t2\t.\tT\t.\t.\tPASS\tNS=1;MS=0;SC=invariant\tGT\t0\n"
        "chr2\t3\t.\tA\tG\t.\tPASS\tNS=1;MS=0;SC=variant\tGT\t1\n",
        encoding="utf-8",
    )

    masked = tmp_path / "combined.mask.bed"
    masked.write_text("chr1\t2\t4\nchr1\t5\t10\nchr2\t4\t8\n", encoding="utf-8")

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
            "--all-sites",
            str(all_sites),
            "--masked-beds",
            str(masked),
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
    # One combined SVG per active contig (chr1, chr2) — no per-series <h3> headings.
    assert html.count("<details") == 2
    assert html.count("<polyline") >= 6  # 3 series × 2 contigs
    assert "ARGprep version:" in html
    assert "Source options file:" in html
    assert "maf_dir: /tmp/maf" in html
    # Scaffolds with no data should be omitted.
    assert "scaffold99" not in html


def test_summary_report_per_maf_section(tmp_path: Path):
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t10\t0\t0\t0\n", encoding="utf-8")

    all_sites = tmp_path / "combined.chr1.all_sites.vcf"
    all_sites.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\tREF\n"
        "chr1\t1\t.\tA\t.\t.\tPASS\tNS=2;MS=0;SC=invariant\tGT:DP\t0:4\t0:4\t0:4\n"
        "chr1\t2\t.\tC\tT\t.\tPASS\tNS=2;MS=0;SC=variant\tGT:DP\t1:3\t0:3\t0:3\n"
        "chr1\t3\t.\tG\tA\t.\tPASS\tNS=2;MS=0;SC=variant\tGT:DP\t1:3\t1:3\t0:3\n"
        "chr1\t5\t.\tG\t.\t.\tPASS\tNS=1;MS=1;SC=invariant\tGT:DP\t0:3\t.:0\t0:3\n",
        encoding="utf-8",
    )

    masked = tmp_path / "combined.chr1.mask.bed"
    masked.write_text("chr1\t6\t10\n", encoding="utf-8")

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
            "--all-sites", str(all_sites),
            "--masked-beds", str(masked),
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
