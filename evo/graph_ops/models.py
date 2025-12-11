# graph_ops/models.py

from dataclasses import dataclass


@dataclass
class Edge:
    """
    Minimal edge representation for SNV relationships.

    u, v: integer positions (same chromosome; chunk is per-chrom)
    relation: "cooccurring" | "timing" | "divergent"
    loss: True if this is a *_loss relation, otherwise False
    reliability: numeric reliability score from C++ output
    source_file: which file this edge came from (for debugging / logging)
    """
    u: int
    v: int
    relation: str   # "cooccurring" | "timing" | "divergent"
    loss: bool
    reliability: float
    source_file: str = ""
