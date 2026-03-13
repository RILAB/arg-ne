#!/usr/bin/env python3
"""
Simulate haploid sequences with msprime SNP variation and branch-based indels.

Outputs:
  - <prefix>.reference.fa : ungapped ancestral/reference sequence
  - <prefix>.samples.fa   : simulated sample sequences after SNPs and indels
  - <prefix>.indels.tsv   : per-sample indel event table
  - <prefix>.summary.tsv  : summary counts for indel-affected ancestral bp and SNPs
  - <prefix>.maf/         : pairwise MAFs for each sample vs ancestral sequence

Notes:
  - `theta` and `rho` are interpreted as scaled values with `4Ne=4`,
    i.e. per-site mutation/recombination rates are `theta/4` and `rho/4`
    with `population_size=1`.
  - Indels are simulated on tree-sequence branch segments with rate
    `indel_rate * branch_length * genomic_span`, then projected to all
    descendant samples. This means descendant samples share indels according
    to the local genealogy.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np

DNA = np.array(list("ACGT"))


@dataclass(frozen=True)
class IndelEvent:
    sample: str
    shared_event_id: int
    event_index: int
    event_type: str
    position_1based: int
    size: int
    sequence: str


@dataclass(frozen=True)
class MafBlock:
    contig: str
    start0: int
    ancestral_size: int
    ancestral_seq: str
    sample_name: str
    sample_size: int
    sample_seq: str


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sequence-length", type=int, required=True)
    ap.add_argument("--num-samples", type=int, required=True)
    ap.add_argument("--theta", type=float, required=True, help="Scaled mutation parameter; uses theta/4 as the per-site rate")
    ap.add_argument("--rho", type=float, required=True, help="Scaled recombination parameter; uses rho/4 as the per-site rate")
    ap.add_argument(
        "--indel-rate",
        type=float,
        required=True,
        help="Indel rate per base per generation on tree-sequence branches",
    )
    ap.add_argument(
        "--indel-lambda",
        type=float,
        required=True,
        help="Rate parameter for an exponential indel size distribution; mean size is 1/lambda",
    )
    ap.add_argument("--out-prefix", required=True, help="Output prefix")
    ap.add_argument("--seed", type=int, default=1)
    return ap.parse_args()


def random_dna(length: int, rng: np.random.Generator) -> str:
    return "".join(rng.choice(DNA, size=length))


def write_fasta(path: Path, records: list[tuple[str, str]], width: int = 60) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), width):
                handle.write(sequence[start : start + width] + "\n")


def write_indel_table(path: Path, events: list[IndelEvent]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("sample\tshared_event_id\tevent_index\ttype\tposition_1based\tsize\tsequence\n")
        for event in events:
            handle.write(
                f"{event.sample}\t{event.shared_event_id}\t{event.event_index}\t{event.event_type}\t"
                f"{event.position_1based}\t{event.size}\t{event.sequence}\n"
            )


def write_summary(
    path: Path,
    *,
    sequence_length: int,
    indel_affected_bp: int,
    snps_without_indels: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"sequence_length\t{sequence_length}\n")
        handle.write(f"ancestral_bp_with_indel_in_ge1_sample\t{indel_affected_bp}\n")
        handle.write(f"snps_without_overlapping_indel\t{snps_without_indels}\n")


def write_maf(path: Path, blocks: list[MafBlock]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("##maf version=1\n")
        for block in blocks:
            handle.write("\n")
            handle.write("a score=0\n")
            handle.write(
                f"s {block.contig} {block.start0} {block.ancestral_size} + {block.ancestral_size} {block.ancestral_seq}\n"
            )
            handle.write(
                f"s {block.sample_name} {block.start0} {block.sample_size} + {block.ancestral_size} {block.sample_seq}\n"
            )


def simulate_snp_haplotypes(
    sequence_length: int,
    num_samples: int,
    theta: float,
    rho: float,
    seed: int,
) -> tuple[str, list[str], object]:
    import msprime

    if sequence_length <= 0:
        raise ValueError("--sequence-length must be positive")
    if num_samples <= 0:
        raise ValueError("--num-samples must be positive")
    if theta < 0 or rho < 0:
        raise ValueError("--theta and --rho must be non-negative")

    rng = np.random.default_rng(seed)
    reference = list(random_dna(sequence_length, rng))

    ts = msprime.sim_ancestry(
        samples=num_samples,
        ploidy=1,
        sequence_length=sequence_length,
        recombination_rate=rho / 4.0,
        population_size=1.0,
        discrete_genome=True,
        random_seed=seed,
    )
    mts = msprime.sim_mutations(
        ts,
        rate=theta / 4.0,
        model=msprime.JC69(),
        random_seed=seed + 1,
    )

    # Align the reference to the simulated ancestral states where mutations occur.
    for site in mts.sites():
        pos = int(site.position)
        if 0 <= pos < sequence_length and site.ancestral_state in {"A", "C", "G", "T"}:
            reference[pos] = site.ancestral_state

    haplotypes = [reference.copy() for _ in range(num_samples)]
    for variant in mts.variants():
        pos = int(variant.site.position)
        if not (0 <= pos < sequence_length):
            continue
        for sample_index, genotype in enumerate(variant.genotypes):
            if genotype < 0:
                continue
            allele = variant.alleles[genotype]
            if allele in {"A", "C", "G", "T"}:
                haplotypes[sample_index][pos] = allele

    return "".join(reference), ["".join(h) for h in haplotypes], mts


def sample_lineage_event(
    *,
    shared_event_id: int,
    descendants: list[str],
    position_1based: int,
    event_type: str,
    size: int,
    sequence: str,
) -> list[IndelEvent]:
    events: list[IndelEvent] = []
    for event_index, sample_name in enumerate(descendants, start=1):
        events.append(
            IndelEvent(
                sample=sample_name,
                shared_event_id=shared_event_id,
                event_index=event_index,
                event_type=event_type,
                position_1based=position_1based,
                size=size,
                sequence=sequence,
            )
        )
    return events


def sample_indel_events_on_ts(
    ts,
    sample_names: list[str],
    indel_rate: float,
    indel_lambda: float,
    rng: np.random.Generator,
) -> list[IndelEvent]:
    if indel_rate < 0:
        raise ValueError("--indel-rate must be non-negative")
    if indel_lambda <= 0:
        raise ValueError("--indel-lambda must be positive")

    events: list[IndelEvent] = []
    shared_event_id = 0
    sample_nodes = list(ts.samples())
    node_to_sample = {node_id: sample_names[i] for i, node_id in enumerate(sample_nodes)}

    for tree in ts.trees():
        left = int(tree.interval.left)
        right = int(tree.interval.right)
        span = max(0, right - left)
        if span == 0:
            continue
        for node in tree.nodes():
            parent = tree.parent(node)
            if parent == -1:
                continue
            branch_length = ts.node(parent).time - ts.node(node).time
            if branch_length <= 0:
                continue
            descendants = [node_to_sample[s] for s in tree.samples(node) if s in node_to_sample]
            if not descendants:
                continue

            num_events = int(rng.poisson(indel_rate * branch_length * span))
            for _ in range(num_events):
                shared_event_id += 1
                event_type = "ins" if rng.random() < 0.5 else "del"
                size = max(1, int(np.ceil(rng.exponential(scale=1.0 / indel_lambda))))
                if event_type == "ins":
                    pos0 = int(rng.integers(left, right + 1))
                    seq = random_dna(size, rng)
                else:
                    pos0 = int(rng.integers(left, right))
                    seq = ""
                events.extend(
                    sample_lineage_event(
                        shared_event_id=shared_event_id,
                        descendants=descendants,
                        position_1based=pos0 + 1,
                        event_type=event_type,
                        size=size,
                        sequence=seq,
                    )
                )
    return events


def apply_indel_events(sequence: str, events: list[IndelEvent]) -> tuple[str, list[IndelEvent]]:
    seq = list(sequence)
    applied: list[IndelEvent] = []

    # Apply from high to low coordinates so reference-based positions remain stable.
    for event in sorted(events, key=lambda e: (e.position_1based, e.event_index), reverse=True):
        pos0 = max(0, min(event.position_1based - 1, len(seq)))
        if event.event_type == "ins":
            seq[pos0:pos0] = list(event.sequence)
            applied.append(event)
            continue

        deleted = "".join(seq[pos0 : pos0 + event.size])
        if not deleted:
            continue
        del seq[pos0 : pos0 + len(deleted)]
        applied.append(
            IndelEvent(
                sample=event.sample,
                shared_event_id=event.shared_event_id,
                event_index=event.event_index,
                event_type=event.event_type,
                position_1based=event.position_1based,
                size=len(deleted),
                sequence=deleted,
            )
        )

    applied.reverse()
    return "".join(seq), applied


def align_sample_to_reference(
    reference: str,
    sample_name: str,
    sample_sequence: str,
    events: list[IndelEvent],
) -> list[MafBlock]:
    insertions_by_anchor: dict[int, list[str]] = defaultdict(list)
    deletions_by_start: dict[int, int] = {}
    for event in events:
        pos0 = event.position_1based - 1
        if event.event_type == "ins":
            insertions_by_anchor[pos0].append(event.sequence)
        else:
            deletions_by_start[pos0] = max(deletions_by_start.get(pos0, 0), event.size)

    ref_aligned: list[str] = []
    sample_aligned: list[str] = []
    sample_idx = 0
    ref_idx = 0

    while ref_idx < len(reference):
        ref_base = reference[ref_idx]
        for ins_seq in insertions_by_anchor.get(ref_idx, []):
            ref_aligned.extend("-" * len(ins_seq))
            sample_aligned.extend(ins_seq)
            sample_idx += len(ins_seq)

        del_size = deletions_by_start.get(ref_idx, 0)
        if del_size > 0:
            for del_ref_idx in range(ref_idx, min(ref_idx + del_size, len(reference))):
                ref_aligned.append(reference[del_ref_idx])
                sample_aligned.append("-")
            ref_idx += del_size
            continue

        ref_aligned.append(ref_base)
        if sample_idx >= len(sample_sequence):
            sample_aligned.append("-")
        else:
            sample_aligned.append(sample_sequence[sample_idx])
            sample_idx += 1
        ref_idx += 1

    for ins_seq in insertions_by_anchor.get(len(reference), []):
        ref_aligned.extend("-" * len(ins_seq))
        sample_aligned.extend(ins_seq)

    aligned_ref = "".join(ref_aligned)
    aligned_sample = "".join(sample_aligned)
    return [
        MafBlock(
            contig="ancestral",
            start0=0,
            ancestral_size=len(reference),
            ancestral_seq=aligned_ref,
            sample_name=sample_name,
            sample_size=len(sample_sequence),
            sample_seq=aligned_sample,
        )
    ]


def summarize_reference_overlaps(
    reference: str,
    haplotypes: list[str],
    all_events: list[IndelEvent],
) -> tuple[int, int]:
    indel_positions: set[int] = set()
    for event in all_events:
        if event.event_type != "del":
            continue
        start0 = event.position_1based - 1
        for pos0 in range(start0, min(start0 + event.size, len(reference))):
            indel_positions.add(pos0)

    snp_positions: set[int] = set()
    for pos0, ref_base in enumerate(reference):
        if pos0 in indel_positions:
            continue
        if any(sample[pos0] != ref_base for sample in haplotypes):
            snp_positions.add(pos0)

    return len(indel_positions), len(snp_positions)


def main() -> None:
    args = parse_args()
    out_prefix = Path(args.out_prefix)
    ref_out = Path(str(out_prefix) + ".reference.fa")
    samples_out = Path(str(out_prefix) + ".samples.fa")
    indels_out = Path(str(out_prefix) + ".indels.tsv")
    summary_out = Path(str(out_prefix) + ".summary.tsv")
    maf_dir = Path(str(out_prefix) + ".maf")

    reference, haplotypes, ts = simulate_snp_haplotypes(
        sequence_length=args.sequence_length,
        num_samples=args.num_samples,
        theta=args.theta,
        rho=args.rho,
        seed=args.seed,
    )

    rng = np.random.default_rng(args.seed + 2)
    sample_names = [f"sample{i}" for i in range(1, args.num_samples + 1)]
    proposed_events = sample_indel_events_on_ts(
        ts=ts,
        sample_names=sample_names,
        indel_rate=args.indel_rate,
        indel_lambda=args.indel_lambda,
        rng=rng,
    )
    sample_records: list[tuple[str, str]] = []
    all_events: list[IndelEvent] = []
    per_sample_events: dict[str, list[IndelEvent]] = defaultdict(list)
    pre_indel_haplotypes = list(haplotypes)
    for sample_name, haplotype in zip(sample_names, haplotypes, strict=True):
        sample_events = [event for event in proposed_events if event.sample == sample_name]
        sequence_with_indels, applied_events = apply_indel_events(haplotype, sample_events)
        sample_records.append((sample_name, sequence_with_indels))
        all_events.extend(applied_events)
        per_sample_events[sample_name].extend(applied_events)

    write_fasta(ref_out, [("reference", reference)])
    write_fasta(samples_out, sample_records)
    write_indel_table(indels_out, all_events)
    indel_bp, snps_without_indels = summarize_reference_overlaps(
        reference=reference,
        haplotypes=pre_indel_haplotypes,
        all_events=all_events,
    )
    write_summary(
        summary_out,
        sequence_length=len(reference),
        indel_affected_bp=indel_bp,
        snps_without_indels=snps_without_indels,
    )
    for sample_name, sequence in sample_records:
        maf_blocks = align_sample_to_reference(
            reference=reference,
            sample_name=sample_name,
            sample_sequence=sequence,
            events=per_sample_events[sample_name],
        )
        write_maf(maf_dir / f"{sample_name}.maf", maf_blocks)


if __name__ == "__main__":
    main()
