import gzip
import runpy
from pathlib import Path
from types import SimpleNamespace


def _run_summary_report(tmp_path: Path, monkeypatch) -> str:
    repo = Path(__file__).resolve().parents[1]
    monkeypatch.chdir(tmp_path)

    # The report scans only explicitly provided current-run warning logs.
    (tmp_path / "logs").mkdir()
    warning_log = tmp_path / "logs" / "test.log"
    warning_log.write_text(
        "Shell command:\n"
        '  echo "WARNING: from shell command block should be ignored"\n'
        'echo "WARNING: explicit echo line should be ignored"\n'
        "WARNING: real log warning should be kept\n",
        encoding="utf-8",
    )

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    (maf_dir / "plain.maf").write_text(
        "##maf version=1\n"
        "a score=0\n"
        "s chr1 0 4 + 4 ACGT\n"
        "s sampleA 0 4 + 4 ACGT\n",
        encoding="utf-8",
    )
    with gzip.open(maf_dir / "extra.maf.gz", "wt", encoding="utf-8") as handle:
        handle.write(
            "##maf version=1\n"
            "a score=0\n"
            "s maf_only_contig 0 4 + 4 ACGT\n"
            "s sampleB 0 4 + 4 ACGT\n"
        )

    ref_fa = tmp_path / "ref.fa"
    ref_fa.write_text(">1\nACGTACGTAC\n>fai_only_contig\nACGT\n", encoding="utf-8")
    ref_fai = tmp_path / "ref.fa.fai"
    ref_fai.write_text("1\t10\t0\t0\t0\n", encoding="utf-8")

    status_tsv = tmp_path / "split.status.tsv"
    status_tsv.write_text(
        "sampleToRef\tchr1\tmissing\n"
        "sampleToRef\tchr1\tmissing\n"
        "sampleToRef\tchr2\tpresent\n",
        encoding="utf-8",
    )

    report_path = tmp_path / "results" / "summary.html"
    fake_snakemake = SimpleNamespace(
        output=SimpleNamespace(report=str(report_path)),
        input=SimpleNamespace(beds=[], invs=[], cleans=[], missing_gt_stats=[]),
        params=SimpleNamespace(
            contigs=["1"],
            jobs=[],
            temp_paths=[],
            arg_outputs=[],
            split_prefixes={},
            split_status_files=[str(status_tsv)],
            dropped_contigs_not_in_ref=[],
            requested_contigs=["chr1"],
            remapped_contigs=[("chr1", "1")],
            contigs_not_in_all_mafs=[],
            ploidy=1,
            ploidy_source="test",
            ploidy_file_values={},
            ploidy_warnings=[],
            maf_dir=str(maf_dir),
            orig_ref_fasta=str(ref_fa),
            ref_fai=str(ref_fai),
            warning_logs=[str(warning_log)],
        ),
    )

    runpy.run_path(
        str(repo / "scripts" / "summary_report.py"),
        init_globals={"snakemake": fake_snakemake},
        run_name="__main__",
    )
    return report_path.read_text(encoding="utf-8", errors="ignore")


def test_summary_report_warning_filters_and_deduplicates(monkeypatch, tmp_path: Path):
    html = _run_summary_report(tmp_path, monkeypatch)

    assert "Configured contigs were remapped to renamed-reference contigs" in html
    assert "chr1-&gt;1" in html

    assert "MAF contigs not present in reference" in html
    assert "maf_only_contig" in html
    assert "MAF contigs not present in reference (showing up to 5): chr1" not in html

    assert "Reference contigs not present in MAFs" in html
    assert "fai_only_contig" in html
    assert "Reference contigs not present in MAFs (showing up to 5): 1" not in html

    assert "gVCF sampleToRef.gvcf.gz is missing configured contigs: chr1" in html
    assert "gVCF sampleToRef.gvcf.gz is missing configured contigs: chr1, chr1" not in html

    assert "real log warning should be kept" in html
    assert "from shell command block should be ignored" not in html
    assert "explicit echo line should be ignored" not in html
