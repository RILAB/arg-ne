import subprocess
import sys
from pathlib import Path


def test_summary_report_html_contains_separate_plots_per_contig(tmp_path: Path):
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t10\t0\t0\t0\nchr2\t8\t0\t0\t0\n", encoding="utf-8")

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

    masked = tmp_path / "combined.masked.bed"
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
    assert html.count("<h3>Invariant (%)</h3>") == 2
    assert html.count("<h3>Variable (%)</h3>") == 2
    assert html.count("<h3>Missing (%)</h3>") == 2
