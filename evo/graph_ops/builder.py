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
      - between any two cooccurring clusters, accepted non-cooccurring edges must
        agree on a single relation mode:
            * timing edges only, all in the same cluster direction, or
            * divergent edges only
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

        # Node-level VAFs (one VAF per SNV/node)
        self.node_vaf: Dict[int, float] = {}

        # Node-level chromosome labels (one chrom per SNV/node)
        # NOTE: Node IDs in the JSON export are still the integer positions.
        # The chromosome label is stored as an additional attribute.
        self.node_chrom: Dict[int, str] = {}
        # Node-level haplotype tags (e.g., HP1/HP2/UNKNOWN). Optional.
        self.node_haplotype: Dict[int, str] = {}
        # Node-level HP read counts (e.g., "12/3" for HP1/HP2 reads). Optional.
        self.node_hp_reads: Dict[int, str] = {}

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

    # -------- Public DSU accessors (for stats/export) --------

    def get_cluster_rep(self, x: int) -> int:
        """Return the cooccurrence-cluster representative for node x."""
        self._ensure_node_dsu(x)
        return self._find(x)

    def get_cluster_sizes(self) -> Dict[int, int]:
        """Return {rep: size} for all current cooccurrence clusters."""
        return {rep: len(nodes) for rep, nodes in self.cluster_nodes.items()}

    # -------- VAF registration --------

    def _register_node_vaf(self, node: int, vaf: float) -> None:
        """
        Record VAF for a node if not already set.

        We assume VAF is an intrinsic property of the SNV; if multiple edges
        mention the same node with slightly different VAFs, we keep the first
        one and ignore later discrepancies.
        """
        if node not in self.node_vaf:
            self.node_vaf[node] = vaf

    def _register_node_chrom(self, node: int, chrom: str, edge: Edge) -> bool:
        """Record chromosome for a node, checking for internal consistency.

        Returns True if the chromosome label is (or remains) consistent, False if
        we detect a mismatch (we keep the original label in that case).
        """
        if node not in self.node_chrom:
            self.node_chrom[node] = chrom
            return True

        if self.node_chrom[node] != chrom:
            # This should not happen for per-chromosome chunk graphs, but if it
            # does, we log it for visibility. We do NOT reject the edge here
            # because the rest of the code assumes node IDs are positions.
            self._log_conflict(edge, None, "node_chrom_mismatch")
            return False

        return True

    def _register_node_haplotype(self, node: int, haplotype: Optional[str], edge: Edge) -> bool:
        """
        Record haplotype tag for a node, checking for internal consistency.
        If conflicting non-empty tags are seen, mark node as MIXED.
        """
        if haplotype is None:
            return True
        hp = str(haplotype).strip().upper()
        if not hp:
            return True

        existing = self.node_haplotype.get(node)
        if existing is None:
            self.node_haplotype[node] = hp
            return True

        if existing == hp:
            return True

        if existing == "MIXED":
            return True

        self.node_haplotype[node] = "MIXED"
        self._log_conflict(edge, None, "node_haplotype_mismatch")
        return False

    def _register_node_hp_reads(self, node: int, hp_reads: Optional[str], edge: Edge) -> bool:
        """
        Record HP read-count string for a node, checking for consistency.
        If conflicting non-empty values are seen, mark as MIXED.
        """
        if hp_reads is None:
            return True

        s = str(hp_reads).strip()
        if not s:
            return True

        existing = self.node_hp_reads.get(node)
        if existing is None:
            self.node_hp_reads[node] = s
            return True

        if existing == s or existing == "MIXED":
            return True

        self.node_hp_reads[node] = "MIXED"
        self._log_conflict(edge, None, "node_hp_reads_mismatch")
        return False

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
    ):
        # Read counts for the skipped edge
        read_counts_str = (
            f"{new_edge.alt_alt}/{new_edge.alt_ref}/"
            f"{new_edge.ref_alt}/{new_edge.ref_ref}"
        )

        rec = {
            "u": new_edge.u,
            "v": new_edge.v,
            "relation": new_edge.relation,
            "loss": new_edge.loss,
            "reliability": new_edge.reliability,
            "read_counts": read_counts_str,   # ALT_ALT/ALT_REF/REF_ALT/REF_REF
            "reason": reason,
            "conflict_u": None,
            "conflict_v": None,
            "conflict_relation_label": None,  # e.g. "timing_loss", "cooccurring"
            "conflict_reliability": None,
        }

        if conflict_edge is not None:
            rec["conflict_u"] = conflict_edge.u
            rec["conflict_v"] = conflict_edge.v
            # Combine relation + loss into a single label
            if conflict_edge.loss:
                rec["conflict_relation_label"] = f"{conflict_edge.relation}_loss"
            else:
                rec["conflict_relation_label"] = conflict_edge.relation
            rec["conflict_reliability"] = conflict_edge.reliability

        self.inconsistencies.append(rec)


    # -------- Cooccurring cluster merge check --------
    def _check_cooccurring_merge_internal_consistency(self, u: int, v: int, new_edge: Edge) -> bool:
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
                    )
                    return False

        return True

    def _check_cooccurring_merge_boundary_consistency(self, u: int, v: int, new_edge: Edge) -> bool:
        """
        Guard against creating mixed relation modes across cluster boundaries after
        a cooccurring merge.

        Why this exists:
          - pair/cluster checks for timing/divergent run when those edges are added
          - later cooccurring merges can collapse previously distinct cluster pairs
            into one, which can accidentally combine:
              * timing + divergent, or
              * opposite timing directions
            between the same merged cluster pair
        """
        self._ensure_node_dsu(u)
        self._ensure_node_dsu(v)
        ru, rv = self._find(u), self._find(v)
        if ru == rv:
            return True

        merged_nodes = self.cluster_nodes[ru] | self.cluster_nodes[rv]
        # For each outside cluster, track one boundary relation signature:
        #   ("divergent", None), ("timing", "out"), or ("timing", "in")
        boundary_signature: Dict[int, Tuple[str, Optional[str]]] = {}
        boundary_edge: Dict[int, Edge] = {}

        for x in merged_nodes:
            for y, e2 in self.adj.get(x, {}).items():
                ry = self._find(y)
                if ry == ru or ry == rv:
                    # Internal to the to-be-merged cluster, handled elsewhere.
                    continue

                if e2.relation == "cooccurring":
                    # Cooccurring edges are not expected across cluster boundaries.
                    continue

                if e2.relation == "timing":
                    # "out" means merged cluster -> outside cluster.
                    # "in"  means outside cluster -> merged cluster.
                    direction = "out" if e2.u == x else "in"
                    sig = ("timing", direction)
                elif e2.relation == "divergent":
                    sig = ("divergent", None)
                else:
                    # Defensive fallback; unknown relation types are rejected.
                    self._log_conflict(new_edge, e2, "unknown_relation")
                    return False

                prev = boundary_signature.get(ry)
                if prev is None:
                    boundary_signature[ry] = sig
                    boundary_edge[ry] = e2
                    continue

                if prev != sig:
                    reason = "cluster_merge_pair_relation_conflict"
                    if prev[0] == "timing" and sig[0] == "timing":
                        reason = "cluster_merge_pair_timing_direction_conflict"
                    self._log_conflict(new_edge, boundary_edge[ry], reason)
                    return False

        return True

    def _check_cooccurring_merge_timing_acyclic(self, u: int, v: int, new_edge: Edge) -> bool:
        """
        Reject cooccurring merges between clusters that are already ordered by
        timing, because contracting such clusters would create a cluster-level
        directed timing cycle.
        """
        self._ensure_node_dsu(u)
        self._ensure_node_dsu(v)
        ru, rv = self._find(u), self._find(v)
        if ru == rv:
            return True

        fwd, _ = self._build_cluster_timing_graph()
        if self._cluster_reachable(ru, rv, fwd) or self._cluster_reachable(rv, ru, fwd):
            self._log_conflict(new_edge, None, "cluster_merge_timing_cycle")
            return False
        return True

    def _check_cluster_pair_relation_ok(self, u: int, v: int, new_edge: Edge) -> bool:
        """
        Enforce cluster-level consistency for non-cooccurring edges.

        For the current cluster pair (cluster(u), cluster(v)):
          - timing can coexist only with timing in the same direction
          - divergent can coexist only with divergent
        """
        self._ensure_node_dsu(u)
        self._ensure_node_dsu(v)
        ru, rv = self._find(u), self._find(v)
        if ru == rv:
            # Internal timing/divergent is checked elsewhere.
            return True

        # Iterate over smaller cluster for efficiency.
        if len(self.cluster_nodes[ru]) <= len(self.cluster_nodes[rv]):
            small_rep, large_rep = ru, rv
        else:
            small_rep, large_rep = rv, ru

        new_from_rep = self._find(new_edge.u)
        new_to_rep = self._find(new_edge.v)

        for x in self.cluster_nodes[small_rep]:
            for y, e2 in self.adj.get(x, {}).items():
                # Only inspect edges crossing this same cluster boundary.
                if self._find(y) != large_rep:
                    continue

                if new_edge.relation == "timing":
                    if e2.relation != "timing":
                        self._log_conflict(
                            new_edge,
                            e2,
                            "cluster_pair_relation_conflict",
                        )
                        return False

                    e2_from_rep = self._find(e2.u)
                    e2_to_rep = self._find(e2.v)
                    if e2_from_rep != new_from_rep or e2_to_rep != new_to_rep:
                        self._log_conflict(
                            new_edge,
                            e2,
                            "cluster_pair_timing_direction_conflict",
                        )
                        return False

                elif new_edge.relation == "divergent":
                    if e2.relation != "divergent":
                        self._log_conflict(
                            new_edge,
                            e2,
                            "cluster_pair_relation_conflict",
                        )
                        return False

        return True

    def _build_cluster_timing_graph(self) -> Tuple[Dict[int, Set[int]], Dict[int, Set[int]]]:
        """
        Build cluster-level timing adjacency from currently accepted timing edges.
        Returns (forward_adj, reverse_adj).
        """
        fwd: Dict[int, Set[int]] = defaultdict(set)
        rev: Dict[int, Set[int]] = defaultdict(set)
        for e in self.edges:
            if e.relation != "timing":
                continue
            ru = self._find(e.u)
            rv = self._find(e.v)
            if ru == rv:
                continue
            fwd[ru].add(rv)
            rev[rv].add(ru)
        return fwd, rev

    def _cluster_reachable(self, start_rep: int, target_rep: int, fwd: Dict[int, Set[int]]) -> bool:
        """Check if there is a cluster-level timing path start_rep -> ... -> target_rep."""
        if start_rep == target_rep:
            return True
        visited = {start_rep}
        dq = deque([start_rep])
        while dq:
            cur = dq.popleft()
            for nxt in fwd.get(cur, ()):
                if nxt == target_rep:
                    return True
                if nxt not in visited:
                    visited.add(nxt)
                    dq.append(nxt)
        return False

    def _collect_cluster_reach(self, start_rep: int, adj: Dict[int, Set[int]]) -> Set[int]:
        """Collect all cluster reps reachable from start_rep in given adjacency."""
        out = {start_rep}
        dq = deque([start_rep])
        while dq:
            cur = dq.popleft()
            for nxt in adj.get(cur, ()):
                if nxt not in out:
                    out.add(nxt)
                    dq.append(nxt)
        return out

    def _build_cluster_divergent_pairs(self) -> Set[Tuple[int, int]]:
        """
        Build set of cluster-level divergent boundaries as sorted rep tuples.
        """
        pairs: Set[Tuple[int, int]] = set()
        for e in self.edges:
            if e.relation != "divergent":
                continue
            ru = self._find(e.u)
            rv = self._find(e.v)
            if ru == rv:
                continue
            pairs.add((ru, rv) if ru <= rv else (rv, ru))
        return pairs

    def _check_divergent_not_ordered_by_timing(self, u: int, v: int, new_edge: Edge) -> bool:
        """
        Reject divergent edges when clusters are already timing-ordered
        (directly or transitively) in either direction.
        """
        self._ensure_node_dsu(u)
        self._ensure_node_dsu(v)
        ru, rv = self._find(u), self._find(v)
        if ru == rv:
            return True

        fwd, _ = self._build_cluster_timing_graph()
        if self._cluster_reachable(ru, rv, fwd) or self._cluster_reachable(rv, ru, fwd):
            self._log_conflict(new_edge, None, "cluster_pair_divergent_timing_transitive_conflict")
            return False
        return True

    def _check_timing_not_imply_divergent_conflict(self, u: int, v: int, new_edge: Edge) -> bool:
        """
        Reject a timing edge u->v if, after adding it at cluster level, any
        ancestor of cluster(u) can reach any descendant of cluster(v) where that
        ancestor/descendant cluster pair is already connected by a divergent edge.
        """
        self._ensure_node_dsu(u)
        self._ensure_node_dsu(v)
        ru, rv = self._find(u), self._find(v)
        if ru == rv:
            return True

        fwd, rev = self._build_cluster_timing_graph()
        fwd[ru].add(rv)
        rev[rv].add(ru)

        ancestors = self._collect_cluster_reach(ru, rev)
        descendants = self._collect_cluster_reach(rv, fwd)
        divergent_pairs = self._build_cluster_divergent_pairs()

        for a in ancestors:
            for d in descendants:
                if a == d:
                    continue
                key = (a, d) if a <= d else (d, a)
                if key in divergent_pairs:
                    self._log_conflict(new_edge, None, "cluster_pair_timing_implies_divergent_transitive_conflict")
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
            if not self._check_cooccurring_merge_internal_consistency(u, v, edge):
                return False
            if not self._check_cooccurring_merge_boundary_consistency(u, v, edge):
                return False
            if not self._check_cooccurring_merge_timing_acyclic(u, v, edge):
                return False

        elif edge.relation == "timing":
            # No timing edges within a cooccurring cluster
            if self._find(u) == self._find(v):
                self._log_conflict(
                    edge,
                    None,
                    "cluster_internal_timing",
                )
                return False

            # Between two clusters, allow only timing in one direction OR divergent.
            if not self._check_cluster_pair_relation_ok(u, v, edge):
                return False

            # Disallow timing edges that make any existing divergent cluster pair
            # become ordered by timing through transitive reachability.
            if not self._check_timing_not_imply_divergent_conflict(u, v, edge):
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
                )
                return False

            # Between two clusters, allow only timing in one direction OR divergent.
            if not self._check_cluster_pair_relation_ok(u, v, edge):
                return False

            # Disallow divergent edges between clusters that are already ordered
            # by timing via any directed path.
            if not self._check_divergent_not_ordered_by_timing(u, v, edge):
                return False

        else:
            raise ValueError(f"Unknown relation: {edge.relation}")

        # ---- Accept edge ----

        self.nodes.add(u)
        self.nodes.add(v)
        # Record per-node chromosome labels for downstream JSON export.
        self._register_node_chrom(u, edge.chrom, edge)
        self._register_node_chrom(v, edge.chrom, edge)
        self._register_node_haplotype(u, edge.hap_u, edge)
        self._register_node_haplotype(v, edge.hap_v, edge)
        self._register_node_hp_reads(u, edge.hp_reads_u, edge)
        self._register_node_hp_reads(v, edge.hp_reads_v, edge)
        self._register_node_vaf(u, edge.vaf_u)
        self._register_node_vaf(v, edge.vaf_v)

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
