import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.check_split_coverage import summarize_site_and_mask_coverage


def _run(cmd, cwd):
    subprocess.run(cmd, cwd=cwd, check=True)


def test_all_sites_and_mask_cover_full_contig(tmp_path: Path):
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("1\t10\t0\t0\t0\n", encoding="utf-8")

    all_sites = tmp_path / "all_sites.vcf"
    all_sites.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t1\t.\tA\t.\t.\t.\tDP=1\n"
        "1\t2\t.\tC\tT\t.\t.\tDP=1\n"
        "1\t5\t.\tG\t.\t.\t.\tDP=1\n"
        "1\t6\t.\tT\tC\t.\t.\tDP=1\n"
        "1\t7\t.\tA\t.\t.\t.\tDP=1\n"
        "1\t10\t.\tC\t.\t.\t.\tDP=1\n",
        encoding="utf-8",
    )

    mask_bed = tmp_path / "masked.bed"
    mask_bed.write_text(
        "1\t2\t4\n"
        "1\t7\t9\n",
        encoding="utf-8",
    )

    summary = summarize_site_and_mask_coverage(all_sites, mask_bed, fai)
    assert summary == {
        "chrom": "1",
        "site_bp": 6,
        "mask_bed_bp": 4,
        "total_bp": 10,
        "chrom_len": 10,
    }


def test_check_split_coverage_overlap_fails(tmp_path: Path):
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("1\t5\t0\t0\t0\n", encoding="utf-8")

    all_sites = tmp_path / "all_sites.vcf"
    all_sites.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t1\t.\tA\tC\t.\t.\tDP=1\n",
        encoding="utf-8",
    )

    mask_bed = tmp_path / "masked.bed"
    mask_bed.write_text("1\t0\t2\n", encoding="utf-8")

    try:
        _run(
            [
                sys.executable,
                str(Path("scripts") / "check_split_coverage.py"),
                "--site-vcf",
                str(all_sites),
                "--mask-bed",
                str(mask_bed),
                "--fai",
                str(fai),
            ],
            cwd=Path.cwd(),
        )
    except subprocess.CalledProcessError:
        return
    raise AssertionError("check_split_coverage should have failed on overlap")
