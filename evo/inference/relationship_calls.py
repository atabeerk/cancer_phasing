"""Convert fitted pair rows into supported relationship proposals."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Mapping, Tuple

from shared_tree_model import RELATIONS

PairKey = Tuple[str, int, int]


@dataclass(frozen=True)
class PairCall:
    chrom: str
    left: int
    right: int
    best_relation: str
    alt_alt: int
    alt_ref: int
    ref_alt: int
    ref_ref: int
    informative_reads: int
    spanning_reads: int
    bics: Mapping[str, float]

    @property
    def key(self) -> PairKey:
        return self.chrom, self.left, self.right

    @property
    def bic_margin(self) -> float:
        """Half the gap between the two lowest relationship BIC scores."""
        ordered = sorted(self.bics.values())
        return 0.5 * max(0.0, ordered[1] - ordered[0])


def calls_from_rows(rows: Iterable[Mapping[str, object]]) -> Dict[PairKey, PairCall]:
    """Convert saved TSV rows or in-memory inference rows into pair calls."""
    calls: Dict[PairKey, PairCall] = {}
    for row in rows:
        left, right = int(row["left_pos"]), int(row["right_pos"])
        bics = {relation: float(row[f"bic_{relation}"]) for relation in RELATIONS}
        key = row["chrom"], left, right
        best = row["best_pair_relationship"]
        if best not in RELATIONS:
            raise ValueError(f"unknown best relationship: {best}")
        if bics[best] > min(bics.values()) + 1e-7:
            raise ValueError(
                f"declared best relationship is not BIC-optimal: {key}"
            )
        call = PairCall(
            chrom=row["chrom"],
            left=left,
            right=right,
            best_relation=best,
            alt_alt=int(row["alt_alt"]),
            alt_ref=int(row["alt_ref"]),
            ref_alt=int(row["ref_alt"]),
            ref_ref=int(row["ref_ref"]),
            informative_reads=int(row["informative_reads"]),
            spanning_reads=int(row["spanning_reads"]),
            bics=bics,
        )
        if call.key in calls:
            raise ValueError(f"duplicate pair call: {call.key}")
        calls[call.key] = call
    return calls


def passes_relationship_support(call: PairCall) -> bool:
    """Require at least one observation of every defining allele pattern."""
    if call.best_relation == "cooccurring":
        return call.alt_alt >= 1
    if call.best_relation == "left_before_right":
        return call.alt_alt >= 1 and call.alt_ref >= 1
    if call.best_relation == "right_before_left":
        return call.alt_alt >= 1 and call.ref_alt >= 1
    if call.best_relation == "divergent":
        return call.alt_ref >= 1 and call.ref_alt >= 1
    raise ValueError(f"unknown relationship: {call.best_relation}")


def greedy_order(calls: Iterable[PairCall]):
    """Prioritize larger BIC margins, breaking ties deterministically."""
    return sorted(calls, key=lambda call: (
        -call.bic_margin, -call.informative_reads,
        call.chrom, call.left, call.right,
    ))
