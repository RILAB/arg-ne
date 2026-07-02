import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.maf_to_sites import (
    discover_samples,
    load_quality_mask,
    maf_path_for_sample,
    missing_threshold,
    read_contig_sequence,
    summarize_site_and_mask_coverage,
)


def _run(cmd, cwd):
    subprocess.run(cmd, cwd=cwd, check=True)


def _write_pairwise_maf(path: Path, contig: str, ref_seq: str, sample: str, sample_seq: str) -> None:
    path.write_text(
        "##maf version=1\n"
        "a score=0\n"
        f"s {contig} 0 {len(ref_seq.replace('-', ''))} + {len(ref_seq.replace('-', ''))} {ref_seq}\n"
        f"s {sample} 0 {len(sample_seq.replace('-', ''))} + {len(sample_seq.replace('-', ''))} {sample_seq}\n",
        encoding="utf-8",
    )


def _write_pairwise_maf_strand(
    path: Path,
    contig: str,
    ref_seq: str,
    sample: str,
    sample_seq: str,
    strand: str,
    start: int,
    src_size: int,
) -> None:
    path.write_text(
        "##maf version=1\n"
        "a score=0\n"
        f"s {contig} 0 {len(ref_seq.replace('-', ''))} + {len(ref_seq.replace('-', ''))} {ref_seq}\n"
        f"s {sample} {start} {len(sample_seq.replace('-', ''))} {strand} {src_size} {sample_seq}\n",
        encoding="utf-8",
    )


def _write_same_src_pairwise_maf(path: Path, contig: str, ref_seq: str, sample_seq: str) -> None:
    path.write_text(
        "##maf version=1\n"
        "a score=0\n"
        f"s {contig} 0 {len(ref_seq.replace('-', ''))} + {len(ref_seq.replace('-', ''))} {ref_seq}\n"
        f"s {contig} 0 {len(sample_seq.replace('-', ''))} + {len(sample_seq.replace('-', ''))} {sample_seq}\n",
        encoding="utf-8",
    )


def _read_vcf_records(path: Path) -> list[list[str]]:
    records = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            records.append(line.rstrip("\n").split("\t"))
    return records


def _read_vcf_header_line(path: Path) -> str:
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM\t"):
                return line.rstrip("\n")
    raise AssertionError(f"No VCF header line found in {path}")


def _read_bed(path: Path) -> list[tuple]:
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if len(fields) >= 4:
                rows.append((chrom, start, end, fields[3]))
            else:
                rows.append((chrom, start, end))
    return rows


def test_maf_to_sites_emits_expected_records_and_mask(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t8\t6\t8\t9\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGTACGT", "s1", "ACGTACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "ACGTACGT", "s2", "AT-TCCGN")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "chr1",
            "--out-prefix",
            str(out_prefix),
            "--samples",
            "s1",
            "s2",
            "--max-missing-count",
            "0",
            "--mask-indel-adjacent-snps",
        ],
        cwd=Path.cwd(),
    )

    all_sites = Path(str(out_prefix) + ".all_sites.vcf")
    variants = Path(str(out_prefix) + ".vcf")
    masked = Path(str(out_prefix) + ".mask.bed")

    all_records = _read_vcf_records(all_sites)
    variant_records = _read_vcf_records(variants)
    bed = _read_bed(masked)

    assert [record[1] for record in all_records] == ["1", "4", "5", "6", "7"]
    assert [record[1] for record in variant_records] == ["5"]
    assert bed == [("chr1", 1, 3), ("chr1", 7, 8)]


def test_maf_to_sites_emits_per_sample_missing_masks(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t8\t6\t8\t9\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    # s1: fully aligned, no missing data
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGTACGT", "s1", "ACGTACGT")
    # s2: deletion at ref pos 2 (idx 2), N at ref pos 7 (idx 7) treated as missing
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "ACGTACGT", "s2", "AT-TCCGN")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir", str(maf_dir),
            "--reference-fasta", str(ref),
            "--contig", "chr1",
            "--out-prefix", str(out_prefix),
            "--samples", "s1", "s2",
            "--max-missing-count", "0",
        ],
        cwd=Path.cwd(),
    )

    s1_mask = _read_bed(Path(str(out_prefix) + ".s1.missing.bed"))
    s2_mask = _read_bed(Path(str(out_prefix) + ".s2.missing.bed"))

    assert s1_mask == []
    assert s2_mask == [("chr1", 2, 3, "s2"), ("chr1", 7, 8, "s2")]


def test_maf_to_sites_preserves_explicit_sample_order(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nAC\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t2\t6\t2\t3\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "AC", "s1", "AT")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "AC", "s2", "AC")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir", str(maf_dir),
            "--reference-fasta", str(ref),
            "--contig", "chr1",
            "--out-prefix", str(out_prefix),
            "--samples", "s2", "s1",
            "--max-missing-count", "0",
            "--keep-indel-adjacent-snps",
        ],
        cwd=Path.cwd(),
    )

    all_sites = Path(str(out_prefix) + ".all_sites.vcf")
    assert _read_vcf_header_line(all_sites).endswith("\ts2\ts1")
    records = _read_vcf_records(all_sites)
    assert [record[1] for record in records] == ["1", "2"]
    assert records[0][9:] == ["0", "0"]
    assert records[1][9:] == ["0", "1"]


def test_maf_to_sites_per_sample_missing_mask_includes_unaligned_positions(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    # s1: only covers the first base; positions 1-3 have no alignment (code 0)
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "A", "s1", "A")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "ACGT", "s2", "ACGT")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir", str(maf_dir),
            "--reference-fasta", str(ref),
            "--contig", "chr1",
            "--out-prefix", str(out_prefix),
            "--samples", "s1", "s2",
            "--max-missing-count", "0",
        ],
        cwd=Path.cwd(),
    )

    s1_mask = _read_bed(Path(str(out_prefix) + ".s1.missing.bed"))
    s2_mask = _read_bed(Path(str(out_prefix) + ".s2.missing.bed"))

    assert s1_mask == [("chr1", 1, 4, "s1")]  # positions 1, 2, 3 unaligned → merged interval
    assert s2_mask == []


def test_maf_to_sites_keeps_multiallelic_sites_by_default(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nA\n", encoding="utf-8")
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t1\t6\t1\t2\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "A", "s1", "A")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "A", "s2", "C")
    _write_pairwise_maf(maf_dir / "s3.maf", "chr1", "A", "s3", "G")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "chr1",
            "--out-prefix",
            str(out_prefix),
            "--samples",
            "s1",
            "s2",
            "s3",
            "--max-missing-count",
            "0",
        ],
        cwd=Path.cwd(),
    )

    all_sites = Path(str(out_prefix) + ".all_sites.vcf")
    variants = Path(str(out_prefix) + ".vcf")
    masked = Path(str(out_prefix) + ".mask.bed")

    all_records = _read_vcf_records(all_sites)
    variant_records = _read_vcf_records(variants)

    assert len(all_records) == 1
    assert len(variant_records) == 1
    assert all_records[0][3] == "A"
    assert all_records[0][4] == "C,G"
    assert all_records[0][9:] == ["0", "1", "2"]
    assert variant_records == all_records
    assert _read_bed(masked) == []


def test_maf_to_sites_accepts_pairwise_blocks_with_same_src_names(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nA\n", encoding="utf-8")
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t1\t6\t1\t2\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_same_src_pairwise_maf(maf_dir / "s1.maf", "chr1", "A", "A")
    _write_same_src_pairwise_maf(maf_dir / "s2.maf", "chr1", "A", "C")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "chr1",
            "--out-prefix",
            str(out_prefix),
            "--samples",
            "s1",
            "s2",
            "--max-missing-count",
            "0",
        ],
        cwd=Path.cwd(),
    )

    all_records = _read_vcf_records(Path(str(out_prefix) + ".all_sites.vcf"))
    variant_records = _read_vcf_records(Path(str(out_prefix) + ".vcf"))
    summary = Path(str(out_prefix) + ".site_summary.tsv").read_text(encoding="utf-8")

    assert len(all_records) == 1
    assert len(variant_records) == 1
    assert all_records[0][3] == "A"
    assert all_records[0][4] == "C"
    assert "all_sites\t1" in summary
    assert "masked_no_alignment\t0" in summary


def test_maf_to_sites_masks_snp_adjacent_to_insertion(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "AC-GT", "s1", "ATCGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "AC-GT", "s2", "AC-GT")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "chr1",
            "--out-prefix",
            str(out_prefix),
            "--samples",
            "s1",
            "s2",
            "--max-missing-count",
            "0",
            "--mask-indel-adjacent-snps",
        ],
        cwd=Path.cwd(),
    )

    all_sites = Path(str(out_prefix) + ".all_sites.vcf")
    variants = Path(str(out_prefix) + ".vcf")
    masked = Path(str(out_prefix) + ".mask.bed")

    all_records = _read_vcf_records(all_sites)
    assert [record[1] for record in all_records] == ["1", "3", "4"]
    assert [record[1] for record in _read_vcf_records(variants)] == []
    assert _read_bed(masked) == [("chr1", 1, 2)]


def test_maf_to_sites_can_keep_snp_adjacent_to_insertion(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "AC-GT", "s1", "ATCGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "AC-GT", "s2", "AC-GT")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "chr1",
            "--out-prefix",
            str(out_prefix),
            "--samples",
            "s1",
            "s2",
            "--max-missing-count",
            "0",
            "--keep-indel-adjacent-snps",
        ],
        cwd=Path.cwd(),
    )

    all_sites = Path(str(out_prefix) + ".all_sites.vcf")
    variants = Path(str(out_prefix) + ".vcf")
    masked = Path(str(out_prefix) + ".mask.bed")

    assert [record[1] for record in _read_vcf_records(all_sites)] == ["1", "2", "3", "4"]
    assert [record[1] for record in _read_vcf_records(variants)] == ["2"]
    assert _read_bed(masked) == []


def test_maf_to_sites_can_keep_snp_adjacent_to_indel(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t8\t6\t8\t9\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGTACGT", "s1", "ACGTACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "ACGTACGT", "s2", "AT-TCCGN")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "chr1",
            "--out-prefix",
            str(out_prefix),
            "--samples",
            "s1",
            "s2",
            "--max-missing-count",
            "0",
            "--keep-indel-adjacent-snps",
        ],
        cwd=Path.cwd(),
    )

    all_sites = Path(str(out_prefix) + ".all_sites.vcf")
    variants = Path(str(out_prefix) + ".vcf")
    masked = Path(str(out_prefix) + ".mask.bed")

    assert [record[1] for record in _read_vcf_records(all_sites)] == ["1", "2", "4", "5", "6", "7"]
    assert [record[1] for record in _read_vcf_records(variants)] == ["2", "5"]
    assert _read_bed(masked) == [("chr1", 2, 3), ("chr1", 7, 8)]


def test_maf_to_sites_counts_all_deleted_sites_as_missingness(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nAAAA\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "AAAA", "s1", "A---")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "AAAA", "s2", "A---")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "chr1",
            "--out-prefix",
            str(out_prefix),
            "--samples",
            "s1",
            "s2",
            "--max-missing-count",
            "0",
        ],
        cwd=Path.cwd(),
    )

    summary = Path(str(out_prefix) + ".site_summary.tsv").read_text(encoding="utf-8")
    assert "masked_missingness\t3\n" in summary
    assert "masked_no_alignment\t0\n" in summary
    assert "masked_indel_adjacent\t0\n" in summary


def test_maf_to_sites_counts_missing_alignment_blocks_as_no_alignment(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "A", "s1", "A")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "ACGT", "s2", "ACGT")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "chr1",
            "--out-prefix",
            str(out_prefix),
            "--samples",
            "s1",
            "s2",
            "--max-missing-count",
            "0",
        ],
        cwd=Path.cwd(),
    )

    summary = Path(str(out_prefix) + ".site_summary.tsv").read_text(encoding="utf-8")
    assert "masked_no_alignment\t3\n" in summary
    assert "masked_missingness\t0\n" in summary


def test_maf_to_sites_matches_normalized_contig_names(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("1\t4\t3\t4\t5\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGT", "s1", "ACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "ACGT", "s2", "ATGT")

    out_prefix = tmp_path / "results" / "combined.1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "1",
            "--out-prefix",
            str(out_prefix),
            "--samples",
            "s1",
            "s2",
            "--max-missing-count",
            "0",
        ],
        cwd=Path.cwd(),
    )

    all_records = _read_vcf_records(Path(str(out_prefix) + ".all_sites.vcf"))
    variant_records = _read_vcf_records(Path(str(out_prefix) + ".vcf"))
    summary = Path(str(out_prefix) + ".site_summary.tsv").read_text(encoding="utf-8")

    assert [record[1] for record in all_records] == ["1", "2", "3", "4"]
    assert [record[1] for record in variant_records] == ["2"]
    assert "masked_no_alignment\t0\n" in summary


def test_summarize_site_and_mask_coverage_rejects_overlap() -> None:
    with pytest.raises(ValueError, match="overlap detected"):
        summarize_site_and_mask_coverage(
            "chr1",
            4,
            [0, 1],
            [1, 2],
        )


def test_summarize_site_and_mask_coverage_requires_full_contig_coverage() -> None:
    with pytest.raises(ValueError, match="coverage mismatch"):
        summarize_site_and_mask_coverage(
            "chr1",
            4,
            [0, 1],
            [3],
        )


def test_discover_samples_handles_maf_and_maf_gz(tmp_path: Path) -> None:
    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    (maf_dir / "sample_b.maf").write_text("", encoding="utf-8")
    (maf_dir / "sample_a.maf.gz").write_text("", encoding="utf-8")
    (maf_dir / "ignore.txt").write_text("", encoding="utf-8")

    assert discover_samples(maf_dir) == ["sample_a", "sample_b"]


def test_maf_path_for_sample_prefers_plain_then_gz(tmp_path: Path) -> None:
    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    plain = maf_dir / "sample1.maf"
    gz = maf_dir / "sample1.maf.gz"
    plain.write_text("", encoding="utf-8")
    gz.write_text("", encoding="utf-8")

    assert maf_path_for_sample(maf_dir, "sample1") == plain

    plain.unlink()
    assert maf_path_for_sample(maf_dir, "sample1") == gz


def test_maf_path_for_sample_raises_for_missing_sample(tmp_path: Path) -> None:
    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()

    with pytest.raises(FileNotFoundError, match="Missing MAF for sample 'sample1'"):
        maf_path_for_sample(maf_dir, "sample1")


def test_read_contig_sequence_raises_when_contig_is_missing(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Contig 'chr2' not found"):
        read_contig_sequence(ref, "chr2")


@pytest.mark.parametrize("fraction", [-0.1, 1.1])
def test_missing_threshold_rejects_invalid_fraction(fraction: float) -> None:
    with pytest.raises(ValueError, match="--max-missing-fraction must be between 0 and 1"):
        missing_threshold(4, None, fraction)


def test_load_quality_mask_keeps_only_below_threshold(tmp_path: Path) -> None:
    bed = tmp_path / "s.bed"
    bed.write_text(
        "# comment\n"
        "track name=quality\n"
        "chrX\t10\t20\t0.95\n"  # passes threshold -> ignored
        "chrX\t30\t40\t0.5\n"  # below threshold -> low quality
        "chrY\t0\t5\t0.1\n"
        "chrX\tbad\trow\n",  # malformed -> skipped
        encoding="utf-8",
    )
    mask = load_quality_mask(bed, 0.9)

    assert mask.is_low("chrX", 35) is True
    assert mask.is_low("chrX", 15) is False  # passing interval not stored
    assert mask.is_low("chrX", 40) is False  # half-open end excluded
    assert mask.is_low("chrY", 3) is True
    assert mask.is_low("chrZ", 0) is False


def test_quality_mask_treats_plus_strand_bases_as_missing(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGT", "s1", "ACGT")
    # s2 carries a SNP (G) at ref idx 1; without masking this is a variant site.
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "ACGT", "s2", "AGGT")

    quality_dir = tmp_path / "quality"
    quality_dir.mkdir()
    # s2's own forward coordinate of the second aligned base is 1 (start 0, + strand).
    (quality_dir / "s2.bed").write_text("s2\t1\t2\t0.5\n", encoding="utf-8")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir", str(maf_dir),
            "--reference-fasta", str(ref),
            "--contig", "chr1",
            "--out-prefix", str(out_prefix),
            "--samples", "s1", "s2",
            "--max-missing-count", "1",
            "--quality-bed-dir", str(quality_dir),
            "--quality-min", "0.9",
        ],
        cwd=Path.cwd(),
    )

    all_records = _read_vcf_records(Path(str(out_prefix) + ".all_sites.vcf"))
    variant_records = _read_vcf_records(Path(str(out_prefix) + ".vcf"))
    s2_mask = _read_bed(Path(str(out_prefix) + ".s2.missing.bed"))

    # The masked base makes ref idx 1 invariant, so no variant is emitted there.
    assert [record[1] for record in variant_records] == []
    pos2 = next(record for record in all_records if record[1] == "2")
    assert pos2[9:] == ["0", "."]  # s2 now missing at the masked position
    assert s2_mask == [("chr1", 1, 2, "s2")]


def test_quality_mask_handles_minus_strand_coordinates(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGT", "s1", "ACGT")
    # s2 aligns on the minus strand (start 0, src_size 4). The second aligned
    # base (ref idx 1) has forward coordinate 4 - 0 - 1 - 1 = 2.
    _write_pairwise_maf_strand(
        maf_dir / "s2.maf", "chr1", "ACGT", "s2", "AGGT", "-", 0, 4
    )

    quality_dir = tmp_path / "quality"
    quality_dir.mkdir()
    (quality_dir / "s2.bed").write_text("s2\t2\t3\t0.5\n", encoding="utf-8")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir", str(maf_dir),
            "--reference-fasta", str(ref),
            "--contig", "chr1",
            "--out-prefix", str(out_prefix),
            "--samples", "s1", "s2",
            "--max-missing-count", "1",
            "--quality-bed-dir", str(quality_dir),
            "--quality-min", "0.9",
        ],
        cwd=Path.cwd(),
    )

    all_records = _read_vcf_records(Path(str(out_prefix) + ".all_sites.vcf"))
    variant_records = _read_vcf_records(Path(str(out_prefix) + ".vcf"))

    # Minus-strand coordinate 2 maps to ref idx 1, masking the SNP there.
    assert [record[1] for record in variant_records] == []
    pos2 = next(record for record in all_records if record[1] == "2")
    assert pos2[9:] == ["0", "."]


def test_quality_mask_requires_quality_min(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nAC\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t2\t6\t2\t3\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "AC", "s1", "AC")

    quality_dir = tmp_path / "quality"
    quality_dir.mkdir()

    proc = subprocess.run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir", str(maf_dir),
            "--reference-fasta", str(ref),
            "--contig", "chr1",
            "--out-prefix", str(tmp_path / "results" / "combined.chr1"),
            "--samples", "s1",
            "--quality-bed-dir", str(quality_dir),
        ],
        cwd=Path.cwd(),
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert "--quality-bed-dir requires --quality-min" in proc.stderr
