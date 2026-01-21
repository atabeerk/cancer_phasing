# graph_ops/models.py

from dataclasses import dataclass


@dataclass
class Edge:
    """
    Minimal edge representation for SNV relationships.
    """

    chrom: str            # chromosome string (from C++ file)

    u: int                # pos1 in oriented graph (source)
    v: int                # pos2 in oriented graph (target)

    relation: str         # "cooccurring" | "timing" | "divergent"
    loss: bool
    reliability: float    # RELIABILITY from C++ output

    alt_alt: int = 0
    alt_ref: int = 0
    ref_alt: int = 0
    ref_ref: int = 0

    vaf_u: float = 0.0    # VAF for u
    vaf_v: float = 0.0    # VAF for v

    best_score: float = 0.0   # BEST_SCORE from C++
    margin: float = 0.0       # MARGIN from C++

    source_file: str = ""
