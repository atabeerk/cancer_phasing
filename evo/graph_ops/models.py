from dataclasses import dataclass, field
from typing import Set, Dict, Tuple, Any, List


Edge = Tuple[str, str, str]   # (chrom, source, target)


@dataclass
class GraphData:
    base: str
    nodes: Set[str] = field(default_factory=set)

    directed: Set[Edge] = field(default_factory=set)
    directed_loss: Set[Edge] = field(default_factory=set)

    cooccurring: Set[Edge] = field(default_factory=set)
    cooccurring_loss: Set[Edge] = field(default_factory=set)

    divergent: Set[Edge] = field(default_factory=set)

    edge_meta: Dict[Tuple[str, str, str], Dict[str, Any]] = field(default_factory=dict)

    # ----------------------------------------
    # CONSTRUCTOR used by graph_building.py
    # ----------------------------------------
    @classmethod
    def from_parts(
        cls,
        base: str,
        nodes: Set[str],
        directed: Set[Edge],
        directed_loss: Set[Edge],
        cooccurring: Set[Edge],
        cooccurring_loss: Set[Edge],
        divergent: Set[Edge],
        edge_meta: Dict[Tuple[str, str, str], Dict[str, Any]]
    ):
        """
        Create a GraphData object from the pieces assembled in build_graphs().
        """
        return cls(
            base=base,
            nodes=nodes,
            directed=directed,
            directed_loss=directed_loss,
            cooccurring=cooccurring,
            cooccurring_loss=cooccurring_loss,
            divergent=divergent,
            edge_meta=edge_meta
        )

    # Optional helper to check edge counts
    def summary(self) -> Dict[str, int]:
        return {
            "nodes": len(self.nodes),
            "directed": len(self.directed),
            "directed_loss": len(self.directed_loss),
            "cooccurring": len(self.cooccurring),
            "cooccurring_loss": len(self.cooccurring_loss),
            "divergent": len(self.divergent)
        }
