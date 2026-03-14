import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


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


def _read_vcf_records(path: Path) -> list[list[str]]:
    records = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            records.append(line.rstrip("\n").split("\t"))
    return records


def _read_bed(path: Path) -> list[tuple[str, int, int]]:
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            chrom, start, end = line.rstrip("\n").split("\t")
            rows.append((chrom, int(start), int(end)))
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
            "--mask-indels",
            "--treat-n-as-missing",
        ],
        cwd=Path.cwd(),
    )

    all_sites = Path(str(out_prefix) + ".all_sites.vcf")
    variants = Path(str(out_prefix) + ".variants.vcf")
    masked = Path(str(out_prefix) + ".masked.bed")

    all_records = _read_vcf_records(all_sites)
    variant_records = _read_vcf_records(variants)
    bed = _read_bed(masked)

    assert [record[1] for record in all_records] == ["1", "4", "5", "6", "7"]
    assert [record[1] for record in variant_records] == ["5"]
    assert bed == [("chr1", 1, 3), ("chr1", 7, 8)]


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
    variants = Path(str(out_prefix) + ".variants.vcf")
    masked = Path(str(out_prefix) + ".masked.bed")

    all_records = _read_vcf_records(all_sites)
    variant_records = _read_vcf_records(variants)

    assert len(all_records) == 1
    assert len(variant_records) == 1
    assert all_records[0][3] == "A"
    assert all_records[0][4] == "C,G"
    assert all_records[0][9:] == ["0", "1", "2"]
    assert variant_records == all_records
    assert _read_bed(masked) == []


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
            "--mask-indels",
        ],
        cwd=Path.cwd(),
    )

    all_sites = Path(str(out_prefix) + ".all_sites.vcf")
    variants = Path(str(out_prefix) + ".variants.vcf")
    masked = Path(str(out_prefix) + ".masked.bed")

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
            "--mask-indels",
            "--keep-indel-adjacent-snps",
        ],
        cwd=Path.cwd(),
    )

    all_sites = Path(str(out_prefix) + ".all_sites.vcf")
    variants = Path(str(out_prefix) + ".variants.vcf")
    masked = Path(str(out_prefix) + ".masked.bed")

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
            "--mask-indels",
            "--keep-indel-adjacent-snps",
        ],
        cwd=Path.cwd(),
    )

    all_sites = Path(str(out_prefix) + ".all_sites.vcf")
    variants = Path(str(out_prefix) + ".variants.vcf")
    masked = Path(str(out_prefix) + ".masked.bed")

    assert [record[1] for record in _read_vcf_records(all_sites)] == ["1", "2", "4", "5", "6", "7"]
    assert [record[1] for record in _read_vcf_records(variants)] == ["2", "5"]
    assert _read_bed(masked) == [("chr1", 2, 3), ("chr1", 7, 8)]
