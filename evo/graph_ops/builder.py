from collections import defaultdict, deque
from typing import Dict, Set, List, Optional, Tuple

from models import Edge


class GraphBuilder:
    """
    Maintains:
      - cooccurring clusters (union–find)
      - a directed timing graph (DAG)
      - divergent edges

    Enforces:
      - only one relation type per unordered pair of nodes
      - cooccurring clusters may contain only cooccurring edges internally
      - timing edges form an acyclic directed graph
      - nodes in the same cooccurring cluster may not have timing/divergent edges between them
    """

    def __init__(self):
        self.nodes: Set[int] = set()
        self.edges: List[Edge] = []  # accepted edges, in insertion order

        # Canonical edge per unordered pair (independent of direction)
        self.pair_edge: Dict[Tuple[int, int], Edge] = {}

        # Undirected adjacency: for cluster checks & general lookups
        self.adj: Dict[int, Dict[int, Edge]] = defaultdict(dict)

        # Directed adjacency for timing edges only (for cycle checks)
        self.adj_directed: Dict[int, Set[int]] = defaultdict(set)

        # Union–find for cooccurring clusters
        self.parent: Dict[int, int] = {}
        self.size: Dict[int, int] = {}
        self.cluster_nodes: Dict[int, Set[int]] = {}

        # Rejected edges and reasons
        self.inconsistencies: List[Dict] = []

    # -------- Union–find helpers --------

    def _ensure_node_dsu(self, x: int) -> None:
        if x not in self.parent:
            self.parent[x] = x
            self.size[x] = 1
            self.cluster_nodes[x] = {x}

    def _find(self, x: int) -> int:
        # Find with path compression
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != x:
            nxt = self.parent[x]
            self.parent[x] = root
            x = nxt
        return root

    def _union(self, a: int, b: int) -> int:
        ra, rb = self._find(a), self._find(b)
        if ra == rb:
            return ra
        # Union by size: attach smaller to larger
        if self.size[ra] < self.size[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        self.size[ra] += self.size[rb]
        # Merge node sets
        self.cluster_nodes[ra].update(self.cluster_nodes[rb])
        del self.cluster_nodes[rb]
        return ra

    # -------- Timing cycle check --------

    def _exists_directed_path(self, start: int, target: int) -> bool:
        """Check if there is any directed path start -> ... -> target."""
        if start == target:
            return True
        visited = {start}
        dq = deque([start])
        while dq:
            u = dq.popleft()
            for v in self.adj_directed.get(u, ()):
                if v == target:
                    return True
                if v not in visited:
                    visited.add(v)
                    dq.append(v)
        return False

    # -------- Logging --------

    def _log_conflict(
        self,
        new_edge: Edge,
        conflict_edge: Optional[Edge],
        reason: str,
        details: str = "",
    ) -> None:
        rec = {
            "u": new_edge.u,
            "v": new_edge.v,
            "relation": new_edge.relation,
            "loss": new_edge.loss,
            "reliability": new_edge.reliability,
            "reason": reason,
            "conflict_u": None,
            "conflict_v": None,
            "conflict_relation": None,
            "conflict_loss": None,
            "conflict_reliability": None,
            "details": details,
        }

        if reason == "timing_cycle" and not details:
            rec["details"] = (
                f"adding timing edge {new_edge.u} -> {new_edge.v} would create a cycle: "
                f"there exists a directed path from {new_edge.v} back to {new_edge.u}"
            )

        if conflict_edge is not None:
            rec["conflict_u"] = conflict_edge.u
            rec["conflict_v"] = conflict_edge.v
            rec["conflict_relation"] = conflict_edge.relation
            rec["conflict_loss"] = conflict_edge.loss
            rec["conflict_reliability"] = conflict_edge.reliability

        self.inconsistencies.append(rec)

    # -------- Cooccurring cluster merge check --------

    def _check_cooccurring_merge_ok(self, u: int, v: int, new_edge: Edge) -> bool:
        """
        Check if merging clusters of u and v via a cooccurring edge is allowed:
        forbids any non-cooccurring edge between any two nodes that would end up
        in the merged cluster.
        """
        self._ensure_node_dsu(u)
        self._ensure_node_dsu(v)
        ru, rv = self._find(u), self._find(v)
        if ru == rv:
            # Already in same cluster; by invariant there are only cooccurring
            # edges inside this cluster.
            return True

        # Iterate over smaller cluster for efficiency
        if len(self.cluster_nodes[ru]) <= len(self.cluster_nodes[rv]):
            small_rep, large_rep = ru, rv
        else:
            small_rep, large_rep = rv, ru
        small_nodes = self.cluster_nodes[small_rep]

        for x in small_nodes:
            for y, e2 in self.adj.get(x, {}).items():
                # If y is in the other cluster, x–y would become an internal edge
                if self._find(y) == large_rep and e2.relation != "cooccurring":
                    self._log_conflict(
                        new_edge,
                        e2,
                        "cluster_merge_noncooccurring",
                        details=(
                            f"merging clusters of {u} and {v} would place nodes "
                            f"{x} and {y} in the same cluster while they are "
                            f"connected by a non-cooccurring edge"
                        ),
                    )
                    return False

        return True

    # -------- Public API --------

    def try_add_edge(self, edge: Edge) -> bool:
        """
        Attempt to add an edge, enforcing all consistency rules.
        Returns True if added, False if rejected or skipped as a duplicate.
        """
        u, v = edge.u, edge.v

        if u == v:
            # self-loops don't make biological sense here
            self._log_conflict(edge, None, "self_loop")
            return False

        self._ensure_node_dsu(u)
        self._ensure_node_dsu(v)

        # Canonical pair key (unordered)
        key = (u, v) if u <= v else (v, u)
        existing = self.pair_edge.get(key)

        # ---- Pair-level type & direction consistency ----

        if existing is not None:
            if existing.relation != edge.relation:
                # Different relation types for same pair → inconsistent
                self._log_conflict(edge, existing, "pair_type_conflict")
                return False

            if edge.relation == "timing":
                # For timing, direction matters
                if existing.u == edge.u and existing.v == edge.v:
                    # Same directed relation; since we sort edges by reliability,
                    # this is a less reliable duplicate → ignore silently.
                    return False
                else:
                    # Opposite timing direction
                    self._log_conflict(edge, existing, "pair_direction_conflict")
                    return False

            else:
                # Same non-timing relation for this pair → less reliable duplicate
                return False

        # ---- Relation-specific global rules ----

        if edge.relation == "cooccurring":
            if not self._check_cooccurring_merge_ok(u, v, edge):
                return False

        elif edge.relation == "timing":
            # No timing edges within a cooccurring cluster
            if self._find(u) == self._find(v):
                self._log_conflict(
                    edge,
                    None,
                    "cluster_internal_timing",
                    details="timing relation within a cooccurrence cluster is forbidden",
                )
                return False

            # Enforce DAG: adding u -> v must not create a cycle
            if self._exists_directed_path(v, u):
                self._log_conflict(edge, None, "timing_cycle")
                return False

        elif edge.relation == "divergent":
            # No divergent edges within a cooccurring cluster
            if self._find(u) == self._find(v):
                self._log_conflict(
                    edge,
                    None,
                    "cluster_internal_divergent",
                    details="divergent relation within a cooccurrence cluster is forbidden",
                )
                return False

        else:
            raise ValueError(f"Unknown relation: {edge.relation}")

        # ---- Accept edge ----

        self.nodes.add(u)
        self.nodes.add(v)
        self.edges.append(edge)
        self.pair_edge[key] = edge

        # Undirected adjacency
        self.adj[u][v] = edge
        self.adj[v][u] = edge

        # Timing adjacency
        if edge.relation == "timing":
            self.adj_directed[u].add(v)

        # Cooccurrence cluster merge
        if edge.relation == "cooccurring":
            self._union(u, v)

        return True
