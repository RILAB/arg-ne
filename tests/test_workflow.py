import importlib.util
import os
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
    cmd.extend(["--", *targets])
    # Ensure tools (samtools, etc.) from the active conda env are on PATH
    env = os.environ.copy()
    env_bin = str(Path(sys.executable).parent)
    env["PATH"] = env_bin + os.pathsep + env.get("PATH", "")
    env["XDG_CACHE_HOME"] = str(tmp_path / ".cache")
    return subprocess.run(
        cmd,
        cwd=repo_root,
        text=True,
        capture_output=True,
        check=False,
        env=env,
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
                "mask_indel_adjacent_snps: false",
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
    assert (tmp_path / "results" / "maf_by_contig" / "s1" / "1.maf.gz").exists()
    assert (tmp_path / "results" / "maf_by_contig" / "s2" / "1.maf.gz").exists()
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
def test_workflow_default_contig_discovery_normalizes_maf_aliases(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGT", "s1", "ACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "1", "ACGT", "s2", "ATGT")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'samples: ["s1", "s2"]',
                "max_missing_count: 0",
                "mask_indel_adjacent_snps: false",
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
    assert "unmatched or ambiguous: chr01" in (
        result.stderr + result.stdout
    )


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_rejects_any_unmatched_explicit_contig(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "ACGT", "s1", "ACGT")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'contigs: ["1", "missing"]',
                'samples: ["s1"]',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = _run_snakemake(tmp_path, config, str(tmp_path / "results" / "summary.html"))
    assert result.returncode != 0
    assert "unmatched or ambiguous: missing" in (result.stderr + result.stdout)


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_default_discovery_warns_and_skips_unmatched_contigs(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    maf_text = (
        "##maf version=1\n"
        "a score=0\n"
        "s chr1 0 4 + 4 ACGT\n"
        "s s1 0 4 + 4 ACGT\n\n"
        "a score=0\n"
        "s chr2 0 4 + 4 ACGT\n"
        "s s1 0 4 + 4 ACGT\n"
    )
    (maf_dir / "s1.maf").write_text(maf_text, encoding="utf-8")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'samples: ["s1"]',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    all_sites = tmp_path / "results" / "sites" / "combined.1.all_sites.vcf"
    result = _run_snakemake(tmp_path, config, str(all_sites))
    assert result.returncode == 0, result.stderr
    assert "Skipped contigs (no unambiguous match in reference .fai): 2" in (
        result.stderr + result.stdout
    )
    assert all_sites.exists()


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
@pytest.mark.parametrize(
    ("allow_multiallelic", "expected_records", "expected_mask"),
    [(True, 1, ""), (False, 0, "1\t0\t1")],
)
def test_workflow_propagates_multiallelic_policy(
    tmp_path: Path,
    allow_multiallelic: bool,
    expected_records: int,
    expected_mask: str,
) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nA\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "A", "s1", "A")
    _write_pairwise_maf(maf_dir / "s2.maf", "1", "A", "s2", "C")
    _write_pairwise_maf(maf_dir / "s3.maf", "1", "A", "s3", "G")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'samples: ["s1", "s2", "s3"]',
                "max_missing_count: 0",
                f"allow_multiallelic_snps: {str(allow_multiallelic).lower()}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    all_sites = tmp_path / "results" / "sites" / "combined.1.all_sites.vcf"
    result = _run_snakemake(tmp_path, config, str(all_sites))
    assert result.returncode == 0, result.stderr

    records = [
        line for line in all_sites.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert len(records) == expected_records
    mask = tmp_path / "results" / "sites" / "combined.1.mask.bed"
    assert mask.read_text(encoding="utf-8").strip() == expected_mask


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_add_ref_appends_reference_sample_to_final_vcfs(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nAC\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "AC", "s1", "AC")
    _write_pairwise_maf(maf_dir / "s2.maf", "1", "AC", "s2", "AT")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'samples: ["s1", "s2"]',
                "max_missing_count: 0",
                "mask_indel_adjacent_snps: false",
                "allow_multiallelic_snps: true",
                "add_ref: true",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = _run_snakemake(tmp_path, config, str(tmp_path / "results" / "summary.html"))
    assert result.returncode == 0, result.stderr

    all_sites = tmp_path / "results" / "sites" / "combined.1.all_sites.vcf"
    variants = tmp_path / "results" / "sites" / "combined.1.vcf"

    all_lines = all_sites.read_text(encoding="utf-8").splitlines()
    variant_lines = variants.read_text(encoding="utf-8").splitlines()

    assert next(line for line in all_lines if line.startswith("#CHROM\t")).endswith("\ts1\ts2\tREF")
    assert next(line for line in variant_lines if line.startswith("#CHROM\t")).endswith("\ts1\ts2\tREF")

    all_records = [line.split("\t") for line in all_lines if line and not line.startswith("#")]
    variant_records = [line.split("\t") for line in variant_lines if line and not line.startswith("#")]

    assert [record[1] for record in all_records] == ["1", "2"]
    assert [record[1] for record in variant_records] == ["2"]
    assert all_records[0][9:] == ["0", "0", "0"]
    assert all_records[1][9:] == ["0", "1", "0"]
    assert variant_records[0][9:] == ["0", "1", "0"]


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_emit_argweaver_sites(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "ACGT", "s1", "ACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "1", "ACGT", "s2", "AGGT")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'samples: ["s1", "s2"]',
                "max_missing_count: 0",
                "emit_argweaver_sites: true",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = _run_snakemake(tmp_path, config, str(tmp_path / "results" / "summary.html"))
    assert result.returncode == 0, result.stderr

    sites = tmp_path / "results" / "sites" / "combined.1.sites"
    assert sites.exists()
    lines = [ln for ln in sites.read_text(encoding="utf-8").splitlines() if ln]
    assert lines[0] == "NAMES\ts1\ts2"
    assert lines[1] == "REGION\t1\t1\t4"
    site_positions = [ln.split("\t")[0] for ln in lines[2:]]

    variants = tmp_path / "results" / "sites" / "combined.1.vcf"
    variant_positions = [
        ln.split("\t")[1]
        for ln in variants.read_text(encoding="utf-8").splitlines()
        if ln and not ln.startswith("#")
    ]
    assert site_positions == variant_positions
    assert lines[2:] == ["2\tCG"]


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_omits_argweaver_sites_by_default(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "ACGT", "s1", "ACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "1", "ACGT", "s2", "AGGT")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'samples: ["s1", "s2"]',
                "max_missing_count: 0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = _run_snakemake(tmp_path, config, str(tmp_path / "results" / "summary.html"))
    assert result.returncode == 0, result.stderr
    assert not (tmp_path / "results" / "sites" / "combined.1.sites").exists()


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_disabling_argweaver_sites_leaves_stale_output(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "ACGT", "s1", "ACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "1", "ACGT", "s2", "AGGT")

    config = tmp_path / "config.yaml"
    base_config = [
        "maf_dir: maf",
        "reference_fasta: ref.fa",
        "results_dir: results",
        'samples: ["s1", "s2"]',
        "max_missing_count: 0",
    ]
    config.write_text("\n".join(base_config + ["emit_argweaver_sites: true"]) + "\n", encoding="utf-8")

    target = str(tmp_path / "results" / "summary.html")
    result = _run_snakemake(tmp_path, config, target)
    assert result.returncode == 0, result.stderr
    sites = tmp_path / "results" / "sites" / "combined.1.sites"
    assert sites.exists()

    config.write_text("\n".join(base_config + ["emit_argweaver_sites: false"]) + "\n", encoding="utf-8")
    result = _run_snakemake(tmp_path, config, target)
    assert result.returncode == 0, result.stderr
    assert sites.exists()


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_enabling_argweaver_sites_creates_output_when_targeted(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "ACGT", "s1", "ACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "1", "ACGT", "s2", "AGGT")

    config = tmp_path / "config.yaml"
    base_config = [
        "maf_dir: maf",
        "reference_fasta: ref.fa",
        "results_dir: results",
        'samples: ["s1", "s2"]',
        "max_missing_count: 0",
    ]
    config.write_text("\n".join(base_config + ["emit_argweaver_sites: false"]) + "\n", encoding="utf-8")

    target = str(tmp_path / "results" / "summary.html")
    result = _run_snakemake(tmp_path, config, target)
    assert result.returncode == 0, result.stderr
    sites = tmp_path / "results" / "sites" / "combined.1.sites"
    assert not sites.exists()

    config.write_text("\n".join(base_config + ["emit_argweaver_sites: true"]) + "\n", encoding="utf-8")
    result = _run_snakemake(tmp_path, config, str(sites))
    assert result.returncode == 0, result.stderr
    assert sites.exists()


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_handles_sample_names_with_spaces(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nAC\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "sample one.maf", "1", "AC", "sample_one", "AT")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'samples: ["sample one"]',
                "max_missing_count: 0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = _run_snakemake(tmp_path, config, str(tmp_path / "results" / "summary.html"))
    assert result.returncode == 0, result.stderr

    all_sites = tmp_path / "results" / "sites" / "combined.1.all_sites.vcf"
    assert next(
        line for line in all_sites.read_text(encoding="utf-8").splitlines()
        if line.startswith("#CHROM\t")
    ).endswith("\tsample one")


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
