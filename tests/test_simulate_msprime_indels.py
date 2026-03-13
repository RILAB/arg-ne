import subprocess
import sys
from pathlib import Path

import msprime
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.simulate_msprime_indels import (
    IndelEvent,
    apply_indel_events,
    align_sample_to_reference,
    sample_indel_events_on_ts,
    sample_lineage_event,
    summarize_reference_overlaps,
)


def test_apply_indel_events_mixes_insertions_and_deletions():
    sequence = "ACGTAC"
    events = [
        IndelEvent(sample="s1", shared_event_id=1, event_index=1, event_type="ins", position_1based=3, size=2, sequence="TT"),
        IndelEvent(sample="s1", shared_event_id=2, event_index=2, event_type="del", position_1based=6, size=2, sequence=""),
    ]

    new_sequence, applied = apply_indel_events(sequence, events)

    assert new_sequence == "ACTTGTA"
    assert [event.event_type for event in applied] == ["ins", "del"]
    assert applied[1].sequence == "C"


def test_sample_lineage_event_duplicates_shared_event_for_descendants():
    events = sample_lineage_event(
        shared_event_id=9,
        descendants=["sample1", "sample3"],
        position_1based=5,
        event_type="del",
        size=2,
        sequence="",
    )

    assert [event.sample for event in events] == ["sample1", "sample3"]
    assert {event.shared_event_id for event in events} == {9}
    assert [event.event_index for event in events] == [1, 2]


def test_sample_indel_events_on_ts_shares_events_across_descendants():
    ts = msprime.sim_ancestry(
        samples=3,
        ploidy=1,
        sequence_length=20,
        recombination_rate=0,
        population_size=1.0,
        discrete_genome=True,
        random_seed=4,
    )
    events = sample_indel_events_on_ts(
        ts=ts,
        sample_names=["sample1", "sample2", "sample3"],
        indel_rate=2.0,
        indel_lambda=1.0,
        rng=np.random.default_rng(5),
    )

    grouped = {}
    for event in events:
        grouped.setdefault(event.shared_event_id, []).append(event.sample)

    assert grouped
    assert any(len(samples) > 1 for samples in grouped.values())


def test_align_sample_to_reference_emits_expected_gaps():
    blocks = align_sample_to_reference(
        reference="ACGT",
        sample_name="sample1",
        sample_sequence="ATTTT",
        events=[
            IndelEvent(sample="sample1", shared_event_id=1, event_index=1, event_type="del", position_1based=2, size=1, sequence="C"),
            IndelEvent(sample="sample1", shared_event_id=2, event_index=1, event_type="ins", position_1based=3, size=2, sequence="TT"),
        ],
    )

    assert len(blocks) == 1
    assert blocks[0].ancestral_seq == "AC--GT"
    assert blocks[0].sample_seq == "A-TTTT"


def test_summarize_reference_overlaps_counts_deleted_bp_and_clean_snps():
    indel_bp, snps = summarize_reference_overlaps(
        reference="ACGT",
        haplotypes=["ACGT", "ATGT"],
        all_events=[
            IndelEvent(sample="sample1", shared_event_id=1, event_index=1, event_type="del", position_1based=1, size=1, sequence="A")
        ],
    )

    assert indel_bp == 1
    assert snps == 1


def test_simulate_msprime_indels_cli_smoke(tmp_path: Path):
    out_prefix = tmp_path / "sim"
    subprocess.run(
        [
            sys.executable,
            str(Path("scripts") / "simulate_msprime_indels.py"),
            "--sequence-length",
            "50",
            "--num-samples",
            "4",
            "--theta",
            "0.01",
            "--rho",
            "0.01",
            "--indel-rate",
            "0.05",
            "--indel-lambda",
            "1.0",
            "--seed",
            "7",
            "--out-prefix",
            str(out_prefix),
        ],
        cwd=Path.cwd(),
        check=True,
    )

    ref_fa = Path(str(out_prefix) + ".reference.fa")
    samples_fa = Path(str(out_prefix) + ".samples.fa")
    indels_tsv = Path(str(out_prefix) + ".indels.tsv")
    summary_tsv = Path(str(out_prefix) + ".summary.tsv")
    maf_dir = Path(str(out_prefix) + ".maf")

    assert ref_fa.exists()
    assert samples_fa.exists()
    assert indels_tsv.exists()
    assert summary_tsv.exists()
    assert maf_dir.exists()
    assert (maf_dir / "sample1.maf").exists()

    sample_headers = [line.strip() for line in samples_fa.read_text(encoding="utf-8").splitlines() if line.startswith(">")]
    assert sample_headers == [">sample1", ">sample2", ">sample3", ">sample4"]

    lines = indels_tsv.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "sample\tshared_event_id\tevent_index\ttype\tposition_1based\tsize\tsequence"
    summary_lines = summary_tsv.read_text(encoding="utf-8").splitlines()
    assert summary_lines[0] == "metric\tvalue"
