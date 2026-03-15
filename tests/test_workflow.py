import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest


def _write_pairwise_maf(path: Path, contig: str, ref_seq: str, sample: str, sample_seq: str) -> None:
    path.write_text(
        "##maf version=1\n"
        "a score=0\n"
        f"s {contig} 0 {len(ref_seq.replace('-', ''))} + {len(ref_seq.replace('-', ''))} {ref_seq}\n"
        f"s {sample} 0 {len(sample_seq.replace('-', ''))} + {len(sample_seq.replace('-', ''))} {sample_seq}\n",
        encoding="utf-8",
    )


def _run_snakemake(tmp_path: Path, config: Path | None, *targets: str) -> subprocess.CompletedProcess[str]:
    repo_root = Path(__file__).resolve().parents[1]
    cmd = [
        sys.executable,
        "-m",
        "snakemake",
        "--snakefile",
        str(repo_root / "Snakefile"),
        "--directory",
        str(tmp_path),
        "--cores",
        "1",
    ]
    if config is not None:
        cmd.extend(["--configfile", str(config)])
    cmd.extend(targets)
    return subprocess.run(
        cmd,
        cwd=repo_root,
        text=True,
        capture_output=True,
        check=False,
    )


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_remaps_requested_contigs_to_reference_names(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGT", "s1", "ACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "ACGT", "s2", "ATGT")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'contigs: ["chr01"]',
                'samples: ["s1", "s2"]',
                "max_missing_count: 0",
                "mask_indels: false",
                "mask_indel_adjacent_snps: false",
                "treat_n_as_missing: true",
                "allow_multiallelic_snps: true",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = _run_snakemake(tmp_path, config, str(tmp_path / "results" / "summary.html"))
    assert result.returncode == 0, result.stderr

    all_sites = tmp_path / "results" / "sites" / "combined.1.all_sites.vcf"
    assert all_sites.exists()
    assert not (tmp_path / "results" / "sites" / "combined.chr01.all_sites.vcf").exists()

    records = [
        line.rstrip("\n").split("\t")
        for line in all_sites.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert [record[0] for record in records] == ["1", "1", "1", "1"]
    assert [record[1] for record in records] == ["1", "2", "3", "4"]


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_rejects_ambiguous_contig_remap(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n>chr1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGT", "s1", "ACGT")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'contigs: ["chr01"]',
                'samples: ["s1"]',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = _run_snakemake(tmp_path, config, str(tmp_path / "results" / "summary.html"))
    assert result.returncode != 0
    assert "None of the configured contigs are present in reference .fai" in (
        result.stderr + result.stdout
    )


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_requires_explicit_configfile(tmp_path: Path) -> None:
    result = _run_snakemake(tmp_path, None, str(tmp_path / "results" / "summary.html"))
    assert result.returncode != 0
    assert "A config file is required" in (result.stderr + result.stdout)


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
@pytest.mark.parametrize("missing_key", ["maf_dir", "reference_fasta"])
def test_workflow_requires_expected_config_keys(tmp_path: Path, missing_key: str) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "ACGT", "s1", "ACGT")

    config_lines = [
        "maf_dir: maf",
        "reference_fasta: ref.fa",
        "results_dir: results",
        'samples: ["s1"]',
    ]
    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(line for line in config_lines if not line.startswith(f"{missing_key}:"))
        + "\n",
        encoding="utf-8",
    )

    result = _run_snakemake(tmp_path, config, str(tmp_path / "results" / "summary.html"))
    assert result.returncode != 0
    assert f"Missing required config keys: {missing_key}" in (result.stderr + result.stdout)
