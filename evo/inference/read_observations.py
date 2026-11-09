"""Load prepared chromosome pileups and project reads into mutation windows."""
from __future__ import annotations

import gzip
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from joint_clone_model import ReadData


@dataclass(frozen=True)
class Mutation:
    ref: str
    alt: str


def parse_mpileup_bases(raw: str) -> List[str]:
    """Return one symbol per alignment after removing mpileup control syntax."""
    symbols: List[str] = []
    index = 0
    while index < len(raw):
        symbol = raw[index]
        if symbol == "^":
            index += 2
        elif symbol == "$":
            index += 1
        elif symbol in "+-":
            index += 1
            start = index
            while index < len(raw) and raw[index].isdigit():
                index += 1
            if start == index:
                raise ValueError(f"malformed mpileup indel in {raw!r}")
            index += int(raw[start:index])
        else:
            symbols.append(symbol)
            index += 1
    return symbols


@dataclass(frozen=True)
class Candidate:
    chrom: str
    left: int
    right: int
    alt_alt: int
    alt_ref: int
    ref_alt: int
    ref_ref: int
    total: int

    @property
    def informative_reads(self) -> int:
        return self.alt_alt + self.alt_ref + self.ref_alt


def passes_inference_gate(
    candidate: Candidate,
    minimum_informative_support: int,
    minimum_alt_alt_support: int = 2,
) -> bool:
    """Keep pairs with enough ALT-bearing reads or enough ALT/ALT reads."""
    return (
        candidate.informative_reads >= minimum_informative_support
        or candidate.alt_alt >= minimum_alt_alt_support
    )


def candidate_key(candidate: Candidate) -> Tuple[str, int, int]:
    left, right = sorted((candidate.left, candidate.right))
    return candidate.chrom, left, right


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def read_mutations(
    vcf: Path, wanted: Optional[Sequence[Tuple[str, int]]] = None,
) -> Dict[Tuple[str, int], Mutation]:
    wanted_set = set(wanted) if wanted is not None else None
    mutations = {}
    with open_text(vcf) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            key = (fields[0], int(fields[1]))
            if wanted_set is not None and key not in wanted_set:
                continue
            alts = fields[4].split(",")
            if len(fields[3]) != 1 or not alts or len(alts[0]) != 1:
                raise ValueError(f"expected a first single-base ALT at {key}")
            mutations[key] = Mutation(fields[3].upper(), alts[0].upper())
    missing = wanted_set - mutations.keys() if wanted_set is not None else set()
    if missing:
        raise ValueError(f"VCF is missing {len(missing)} candidate positions")
    return mutations


def read_observations(
    pileup_dir: Path,
    mutations: Mapping[Tuple[str, int], Mutation],
    bam_count: int,
) -> Tuple[
    Dict[str, Dict[Tuple[str, int], Tuple[int, float]]],
    Dict[Tuple[str, int], set],
    Dict[str, int],
]:
    """Load one independent observation per BAM alignment record."""
    observations: Dict[str, Dict[Tuple[str, int], Tuple[int, float]]] = {}
    conflicts = set()
    stats = Counter()
    for pileup_path in sorted(pileup_dir.glob("mpileup_*.out")):
        with pileup_path.open() as handle:
            for line_number, line in enumerate(handle, 1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 7:
                    continue
                key = (fields[0], int(fields[1]))
                mutation = mutations.get(key)
                if mutation is None:
                    continue
                stats["pileup_rows"] += 1
                # Evo assigns a unique QNAME to every alignment record,
                # including supplementary records, before creating pileups.
                if len(fields) != 3 + 4 * bam_count:
                    raise ValueError(
                        f"unexpected Evo pileup columns at {pileup_path}:{line_number}"
                    )
                for sample in range(bam_count):
                    offset = 3 + sample * 4
                    if int(fields[offset]) == 0:
                        continue
                    sample_bases = parse_mpileup_bases(fields[offset + 1])
                    sample_qualities = fields[offset + 2]
                    sample_names = fields[offset + 3].split(",")
                    if not (len(sample_bases) == len(sample_qualities) == len(sample_names)):
                        raise ValueError(
                            f"unaligned Evo sample at {pileup_path}:{line_number}"
                        )
                    for base, quality_char, read_name in zip(
                        sample_bases, sample_qualities, sample_names
                    ):
                        alignment_key = f"{read_name}\tbam{sample}:{read_name}"
                        if base in ".," or base.upper() == mutation.ref:
                            allele = 0
                        elif base.upper() == mutation.alt:
                            allele = 1
                        else:
                            stats["non_ref_alt_observations"] += 1
                            continue
                        quality = max(0, ord(quality_char) - 33)
                        error = min(0.49, max(0.01, 10.0 ** (-quality / 10.0)))
                        observation_key = (alignment_key, key)
                        if observation_key in conflicts:
                            continue
                        read_values = observations.setdefault(alignment_key, {})
                        previous = read_values.get(key)
                        if previous is None:
                            read_values[key] = (allele, error)
                            stats["usable_observations"] += 1
                        elif previous[0] == allele:
                            stats["duplicate_observations"] += 1
                            read_values[key] = (allele, min(error, previous[1]))
                        else:
                            stats["conflicting_duplicates"] += 1
                            conflicts.add(observation_key)
                            del read_values[key]

    reads_by_position: Dict[Tuple[str, int], set] = defaultdict(set)
    for read_name, read_values in observations.items():
        for key in read_values:
            reads_by_position[key].add(read_name)
    stats["alignment_records"] = len(observations)
    stats["physical_reads"] = len(observations)
    stats["retained_observations"] = sum(map(len, observations.values()))
    return observations, reads_by_position, dict(stats)


def candidates_from_observations(
    observations: Mapping[str, Mapping[Tuple[str, int], Tuple[int, float]]],
    minimum_informative_reads: int = 2,
) -> List[Candidate]:
    """Count each alignment record once per same-chromosome mutation pair.

    Input observations have already resolved duplicate/conflicting calls.
    The inference gate is applied to these same reconstructed counts.
    """
    counts: Dict[Tuple[str, int, int], List[int]] = {}
    for values in observations.values():
        for (left, left_call), (right, right_call) in combinations(
            sorted(values.items()), 2
        ):
            if left[0] != right[0]:
                continue
            key = left[0], left[1], right[1]
            if key not in counts:
                counts[key] = [0, 0, 0, 0]
            category = (1 - int(left_call[0])) * 2 + (1 - int(right_call[0]))
            counts[key][category] += 1
    return [
        Candidate(chrom, left, right, *values, sum(values))
        for (chrom, left, right), values in sorted(counts.items())
        if sum(values[:3]) >= minimum_informative_reads
    ]


def project_reads(
    positions: Sequence[Tuple[str, int]],
    observations: Mapping[str, Mapping[Tuple[str, int], Tuple[int, float]]],
    reads_by_position: Mapping[Tuple[str, int], set],
) -> List[ReadData]:
    read_names = set()
    for position in positions:
        read_names.update(reads_by_position.get(position, set()))
    position_to_index = {position: index for index, position in enumerate(positions)}
    reads = []
    for read_name in sorted(read_names):
        values = observations[read_name]
        items = [
            (position_to_index[position], values[position])
            for position in positions
            if position in values
        ]
        reads.append(
            ReadData(
                np.asarray([item[0] for item in items], dtype=np.int32),
                np.asarray([item[1][0] for item in items], dtype=np.bool_),
                np.asarray([item[1][1] for item in items], dtype=np.float64),
            )
        )
    return reads
