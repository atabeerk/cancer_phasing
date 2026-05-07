from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple

from models import Edge


@dataclass(frozen=True)
class BoundarySignature:
    relation: str
    direction: Optional[str]
    edge: Edge

    def flipped(self) -> "BoundarySignature":
        if self.relation != "timing":
            return self
        if self.direction == "out":
            return BoundarySignature(self.relation, "in", self.edge)
        if self.direction == "in":
            return BoundarySignature(self.relation, "out", self.edge)
        return self


class FastGraphBuilder:
    """
    Indexed implementation of GraphBuilder.

    This class preserves the legacy edge acceptance order and output-facing
    attributes, but keeps cluster boundary and reachability indexes incrementally
    instead of rebuilding cluster timing/divergent graphs for every candidate.
    """

    def __init__(self):
        self.nodes: Set[int] = set()
        self.edges: List[Edge] = []
        self.pair_edge: Dict[Tuple[int, int], Edge] = {}
        self.adj: Dict[int, Dict[int, Edge]] = defaultdict(dict)
        self.adj_directed: Dict[int, Set[int]] = defaultdict(set)

        self.parent: Dict[int, int] = {}
        self.size: Dict[int, int] = {}
        self.cluster_nodes: Dict[int, Set[int]] = {}

        self.node_vaf: Dict[int, float] = {}
        self.node_chrom: Dict[int, str] = {}
        self.node_haplotype: Dict[int, str] = {}
        self.node_hp_reads: Dict[int, str] = {}
        self.inconsistencies: List[Dict] = []

        self.node_to_bit: Dict[int, int] = {}
        self.bit_to_node: List[int] = []
        self.bit_owner: List[int] = []
        self.member_bits: Dict[int, int] = {}

        # boundary[a][b] is expressed from cluster a's perspective.
        self.boundary: Dict[int, Dict[int, BoundarySignature]] = defaultdict(dict)
        self.divergent_neighbors: Dict[int, int] = defaultdict(int)

        # Transitive closure bitsets over original node bits, interpreted through
        # bit_owner/member_bits for current cooccurrence clusters.
        self.cluster_reach: Dict[int, int] = defaultdict(int)
        self.cluster_rev_reach: Dict[int, int] = defaultdict(int)
        self.node_reach: Dict[int, int] = defaultdict(int)
        self.node_rev_reach: Dict[int, int] = defaultdict(int)

    # -------- Bit helpers --------

    def _ensure_bit(self, node: int) -> int:
        bit = self.node_to_bit.get(node)
        if bit is not None:
            return bit
        bit = len(self.bit_to_node)
        self.node_to_bit[node] = bit
        self.bit_to_node.append(node)
        self.bit_owner.append(node)
        return bit

    @staticmethod
    def _bit_mask(bit: int) -> int:
        return 1 << bit

    def _bits_to_active_reps(self, bits: int) -> Set[int]:
        reps: Set[int] = set()
        while bits:
            lsb = bits & -bits
            bit = lsb.bit_length() - 1
            reps.add(self._find(self.bit_owner[bit]))
            bits ^= lsb
        return reps

    # -------- Union-find helpers --------

    def _ensure_node_dsu(self, x: int) -> None:
        if x in self.parent:
            return
        bit = self._ensure_bit(x)
        self.parent[x] = x
        self.size[x] = 1
        self.cluster_nodes[x] = {x}
        self.member_bits[x] = self._bit_mask(bit)
        self.bit_owner[bit] = x

    def _find(self, x: int) -> int:
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != x:
            nxt = self.parent[x]
            self.parent[x] = root
            x = nxt
        return root

    def get_cluster_rep(self, x: int) -> int:
        self._ensure_node_dsu(x)
        return self._find(x)

    def get_cluster_sizes(self) -> Dict[int, int]:
        return {rep: len(nodes) for rep, nodes in self.cluster_nodes.items()}

    # -------- Public/debug-compatible helpers --------

    def _build_cluster_timing_graph(self):
        fwd: Dict[int, Set[int]] = defaultdict(set)
        rev: Dict[int, Set[int]] = defaultdict(set)
        for a, nbrs in self.boundary.items():
            if self._find(a) != a:
                continue
            for b, sig in nbrs.items():
                if sig.relation == "timing" and sig.direction == "out" and self._find(b) == b:
                    fwd[a].add(b)
                    rev[b].add(a)
        return fwd, rev

    def _build_cluster_divergent_pairs(self) -> Set[Tuple[int, int]]:
        pairs: Set[Tuple[int, int]] = set()
        for a, nbrs in self.boundary.items():
            if self._find(a) != a:
                continue
            for b, sig in nbrs.items():
                if a < b and sig.relation == "divergent" and self._find(b) == b:
                    pairs.add((a, b))
        return pairs

    def _cluster_reachable(self, start_rep: int, target_rep: int, _fwd=None) -> bool:
        start_rep = self._find(start_rep)
        target_rep = self._find(target_rep)
        if start_rep == target_rep:
            return True
        return bool(self.cluster_reach[start_rep] & self.member_bits[target_rep])

    def _collect_cluster_reach(self, start_rep: int, adj=None) -> Set[int]:
        start_rep = self._find(start_rep)
        if adj is None:
            bits = self.cluster_reach[start_rep]
        else:
            # Compatibility fallback for callers expecting legacy behavior.
            seen = {start_rep}
            stack = [start_rep]
            while stack:
                cur = stack.pop()
                for nxt in adj.get(cur, ()):
                    nxt = self._find(nxt)
                    if nxt not in seen:
                        seen.add(nxt)
                        stack.append(nxt)
            return seen
        return self._bits_to_active_reps(bits) | {start_rep}

    def _exists_directed_path(self, start: int, target: int) -> bool:
        self._ensure_node_dsu(start)
        self._ensure_node_dsu(target)
        if start == target:
            return True
        return bool(self.node_reach[start] & self._bit_mask(self.node_to_bit[target]))

    # -------- Logging and node metadata --------

    def _log_conflict(self, new_edge: Edge, conflict_edge: Optional[Edge], reason: str):
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
            "read_counts": read_counts_str,
            "reason": reason,
            "conflict_u": None,
            "conflict_v": None,
            "conflict_relation_label": None,
            "conflict_reliability": None,
        }
        if conflict_edge is not None:
            rec["conflict_u"] = conflict_edge.u
            rec["conflict_v"] = conflict_edge.v
            rec["conflict_relation_label"] = (
                f"{conflict_edge.relation}_loss" if conflict_edge.loss else conflict_edge.relation
            )
            rec["conflict_reliability"] = conflict_edge.reliability
        self.inconsistencies.append(rec)

    def _register_node_vaf(self, node: int, vaf: float) -> None:
        if node not in self.node_vaf:
            self.node_vaf[node] = vaf

    def _register_node_chrom(self, node: int, chrom: str, edge: Edge) -> bool:
        existing = self.node_chrom.get(node)
        if existing is None:
            self.node_chrom[node] = chrom
            return True
        if existing != chrom:
            self._log_conflict(edge, None, "node_chrom_mismatch")
            return False
        return True

    def _register_node_haplotype(self, node: int, haplotype: Optional[str], edge: Edge) -> bool:
        if haplotype is None:
            return True
        hp = str(haplotype).strip().upper()
        if not hp:
            return True
        existing = self.node_haplotype.get(node)
        if existing is None:
            self.node_haplotype[node] = hp
            return True
        if existing == hp or existing == "MIXED":
            return True
        self.node_haplotype[node] = "MIXED"
        self._log_conflict(edge, None, "node_haplotype_mismatch")
        return False

    def _register_node_hp_reads(self, node: int, hp_reads: Optional[str], edge: Edge) -> bool:
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

    # -------- Indexed consistency checks --------

    def _boundary_between(self, a: int, b: int) -> Optional[BoundarySignature]:
        return self.boundary.get(a, {}).get(b)

    def _cluster_reachable_by_rep(self, a: int, b: int) -> bool:
        a = self._find(a)
        b = self._find(b)
        if a == b:
            return True
        return bool(self.cluster_reach[a] & self.member_bits[b])

    def _check_cooccurring_merge_internal_consistency(self, u: int, v: int, new_edge: Edge) -> bool:
        self._ensure_node_dsu(u)
        self._ensure_node_dsu(v)
        ru, rv = self._find(u), self._find(v)
        if ru == rv:
            return True
        sig = self._boundary_between(ru, rv)
        if sig is not None:
            self._log_conflict(new_edge, sig.edge, "cluster_merge_noncooccurring")
            return False
        return True

    def _check_cooccurring_merge_boundary_consistency(self, u: int, v: int, new_edge: Edge) -> bool:
        self._ensure_node_dsu(u)
        self._ensure_node_dsu(v)
        ru, rv = self._find(u), self._find(v)
        if ru == rv:
            return True

        seen: Dict[int, BoundarySignature] = {}
        for rep in (ru, rv):
            for other, sig in self.boundary.get(rep, {}).items():
                other = self._find(other)
                if other == ru or other == rv:
                    continue
                prev = seen.get(other)
                if prev is None:
                    seen[other] = sig
                    continue
                if prev.relation != sig.relation or prev.direction != sig.direction:
                    reason = "cluster_merge_pair_relation_conflict"
                    if prev.relation == "timing" and sig.relation == "timing":
                        reason = "cluster_merge_pair_timing_direction_conflict"
                    self._log_conflict(new_edge, prev.edge, reason)
                    return False
        return True

    def _check_cooccurring_merge_timing_acyclic(self, u: int, v: int, new_edge: Edge) -> bool:
        self._ensure_node_dsu(u)
        self._ensure_node_dsu(v)
        ru, rv = self._find(u), self._find(v)
        if ru == rv:
            return True
        if self._cluster_reachable_by_rep(ru, rv) or self._cluster_reachable_by_rep(rv, ru):
            self._log_conflict(new_edge, None, "cluster_merge_timing_cycle")
            return False
        return True

    def _check_cluster_pair_relation_ok(self, u: int, v: int, new_edge: Edge) -> bool:
        self._ensure_node_dsu(u)
        self._ensure_node_dsu(v)
        ru, rv = self._find(u), self._find(v)
        if ru == rv:
            return True
        sig = self._boundary_between(ru, rv)
        if sig is None:
            return True
        if new_edge.relation == "timing":
            if sig.relation != "timing":
                self._log_conflict(new_edge, sig.edge, "cluster_pair_relation_conflict")
                return False
            if sig.direction != "out":
                self._log_conflict(new_edge, sig.edge, "cluster_pair_timing_direction_conflict")
                return False
        elif new_edge.relation == "divergent":
            if sig.relation != "divergent":
                self._log_conflict(new_edge, sig.edge, "cluster_pair_relation_conflict")
                return False
        return True

    def _check_divergent_not_ordered_by_timing(self, u: int, v: int, new_edge: Edge) -> bool:
        self._ensure_node_dsu(u)
        self._ensure_node_dsu(v)
        ru, rv = self._find(u), self._find(v)
        if ru == rv:
            return True
        if self._cluster_reachable_by_rep(ru, rv) or self._cluster_reachable_by_rep(rv, ru):
            self._log_conflict(new_edge, None, "cluster_pair_divergent_timing_transitive_conflict")
            return False
        return True

    def _check_timing_not_imply_divergent_conflict(self, u: int, v: int, new_edge: Edge) -> bool:
        self._ensure_node_dsu(u)
        self._ensure_node_dsu(v)
        ru, rv = self._find(u), self._find(v)
        if ru == rv:
            return True

        anc_bits = self.cluster_rev_reach[ru] | self.member_bits[ru]
        desc_bits = self.cluster_reach[rv] | self.member_bits[rv]
        for ancestor in self._bits_to_active_reps(anc_bits):
            if self.divergent_neighbors[ancestor] & desc_bits:
                self._log_conflict(new_edge, None, "cluster_pair_timing_implies_divergent_transitive_conflict")
                return False
        return True

    # -------- Index updates --------

    def _add_boundary(self, a: int, b: int, sig_from_a: BoundarySignature) -> None:
        self.boundary[a][b] = sig_from_a
        self.boundary[b][a] = sig_from_a.flipped()

    def _add_divergent_boundary(self, a: int, b: int, edge: Edge) -> None:
        sig = BoundarySignature("divergent", None, edge)
        self._add_boundary(a, b, sig)
        self.divergent_neighbors[a] |= self.member_bits[b]
        self.divergent_neighbors[b] |= self.member_bits[a]

    def _iter_reps_from_bits(self, bits: int) -> Iterable[int]:
        return self._bits_to_active_reps(bits)

    def _add_timing_boundary(self, a: int, b: int, edge: Edge) -> None:
        sig = BoundarySignature("timing", "out", edge)
        self._add_boundary(a, b, sig)

        anc_bits = self.cluster_rev_reach[a] | self.member_bits[a]
        desc_bits = self.cluster_reach[b] | self.member_bits[b]
        for ancestor in self._iter_reps_from_bits(anc_bits):
            self.cluster_reach[ancestor] |= desc_bits
        for descendant in self._iter_reps_from_bits(desc_bits):
            self.cluster_rev_reach[descendant] |= anc_bits

    def _add_node_timing(self, u: int, v: int) -> None:
        ub = self.node_to_bit[u]
        vb = self.node_to_bit[v]
        anc_bits = self.node_rev_reach[u] | self._bit_mask(ub)
        desc_bits = self.node_reach[v] | self._bit_mask(vb)
        for bit in self._iter_set_bits(anc_bits):
            node = self.bit_to_node[bit]
            self.node_reach[node] |= desc_bits
        for bit in self._iter_set_bits(desc_bits):
            node = self.bit_to_node[bit]
            self.node_rev_reach[node] |= anc_bits

    @staticmethod
    def _iter_set_bits(bits: int) -> Iterable[int]:
        while bits:
            lsb = bits & -bits
            bit = lsb.bit_length() - 1
            yield bit
            bits ^= lsb

    def _merge_clusters(self, u: int, v: int) -> int:
        ru, rv = self._find(u), self._find(v)
        if ru == rv:
            return ru

        keep, drop = ru, rv
        if self.size[keep] < self.size[drop]:
            keep, drop = drop, keep

        keep_bits_old = self.member_bits[keep]
        drop_bits = self.member_bits[drop]
        merged_bits = keep_bits_old | drop_bits
        desc_union = (self.cluster_reach[keep] | self.cluster_reach[drop]) & ~merged_bits
        anc_union = (self.cluster_rev_reach[keep] | self.cluster_rev_reach[drop]) & ~merged_bits

        # Update reachability before remapping bit owners, while ancestor/descendant
        # sets still resolve through the old active reps.
        for ancestor in self._bits_to_active_reps(anc_union):
            self.cluster_reach[ancestor] |= merged_bits | desc_union
        for descendant in self._bits_to_active_reps(desc_union):
            self.cluster_rev_reach[descendant] |= merged_bits | anc_union

        self.parent[drop] = keep
        self.size[keep] += self.size[drop]
        self.cluster_nodes[keep].update(self.cluster_nodes[drop])
        del self.cluster_nodes[drop]
        del self.size[drop]

        for bit in self._iter_set_bits(drop_bits):
            self.bit_owner[bit] = keep

        self.member_bits[keep] = merged_bits
        del self.member_bits[drop]
        self.cluster_reach[keep] = desc_union
        self.cluster_rev_reach[keep] = anc_union
        self.cluster_reach.pop(drop, None)
        self.cluster_rev_reach.pop(drop, None)

        self._merge_boundaries(keep, drop, merged_bits)
        return keep

    def _merge_boundaries(self, keep: int, drop: int, merged_bits: int) -> None:
        keep_boundary = self.boundary.get(keep, {})
        drop_boundary = self.boundary.get(drop, {})
        affected_neighbors = set(keep_boundary) | set(drop_boundary)
        affected_neighbors.discard(keep)
        affected_neighbors.discard(drop)

        for other, sig in list(drop_boundary.items()):
            other = self._find(other)
            if other == keep:
                continue
            self.boundary[other].pop(drop, None)
            if other not in self.boundary[keep]:
                self.boundary[keep][other] = sig
                self.boundary[other][keep] = sig.flipped()
            else:
                # Keep the representative edge/signature that was already present
                # on the kept cluster; pre-merge checks have already verified it is compatible.
                self.boundary[other][keep] = self.boundary[keep][other].flipped()

        self.boundary[keep].pop(drop, None)
        self.boundary.pop(drop, None)

        for other in affected_neighbors:
            other = self._find(other)
            if other == keep:
                continue
            sig = self.boundary.get(keep, {}).get(other)
            if sig is not None and sig.relation == "divergent":
                self.divergent_neighbors[other] &= ~merged_bits
                self.divergent_neighbors[other] |= self.member_bits[keep]
            else:
                self.divergent_neighbors[other] &= ~merged_bits

        div_bits = 0
        for other, sig in self.boundary.get(keep, {}).items():
            if sig.relation == "divergent":
                div_bits |= self.member_bits[self._find(other)]
        self.divergent_neighbors[keep] = div_bits & ~self.member_bits[keep]
        self.divergent_neighbors.pop(drop, None)

    # -------- Public API --------

    def try_add_edge(self, edge: Edge) -> bool:
        u, v = edge.u, edge.v
        if u == v:
            self._log_conflict(edge, None, "self_loop")
            return False

        self._ensure_node_dsu(u)
        self._ensure_node_dsu(v)

        key = (u, v) if u <= v else (v, u)
        existing = self.pair_edge.get(key)
        if existing is not None:
            if existing.relation != edge.relation:
                self._log_conflict(edge, existing, "pair_type_conflict")
                return False
            if edge.relation == "timing":
                if existing.u == edge.u and existing.v == edge.v:
                    return False
                self._log_conflict(edge, existing, "pair_direction_conflict")
                return False
            return False

        ru, rv = self._find(u), self._find(v)

        if edge.relation == "cooccurring":
            if not self._check_cooccurring_merge_internal_consistency(u, v, edge):
                return False
            if not self._check_cooccurring_merge_boundary_consistency(u, v, edge):
                return False
            if not self._check_cooccurring_merge_timing_acyclic(u, v, edge):
                return False
        elif edge.relation == "timing":
            if ru == rv:
                self._log_conflict(edge, None, "cluster_internal_timing")
                return False
            if not self._check_cluster_pair_relation_ok(u, v, edge):
                return False
            if not self._check_timing_not_imply_divergent_conflict(u, v, edge):
                return False
            if self._exists_directed_path(v, u):
                self._log_conflict(edge, None, "timing_cycle")
                return False
        elif edge.relation == "divergent":
            if ru == rv:
                self._log_conflict(edge, None, "cluster_internal_divergent")
                return False
            if not self._check_cluster_pair_relation_ok(u, v, edge):
                return False
            if not self._check_divergent_not_ordered_by_timing(u, v, edge):
                return False
        else:
            raise ValueError(f"Unknown relation: {edge.relation}")

        self._accept_edge(edge)
        return True

    def _accept_edge(self, edge: Edge) -> None:
        u, v = edge.u, edge.v
        ru, rv = self._find(u), self._find(v)
        self.nodes.add(u)
        self.nodes.add(v)
        self._register_node_chrom(u, edge.chrom, edge)
        self._register_node_chrom(v, edge.chrom, edge)
        self._register_node_haplotype(u, edge.hap_u, edge)
        self._register_node_haplotype(v, edge.hap_v, edge)
        self._register_node_hp_reads(u, edge.hp_reads_u, edge)
        self._register_node_hp_reads(v, edge.hp_reads_v, edge)
        self._register_node_vaf(u, edge.vaf_u)
        self._register_node_vaf(v, edge.vaf_v)

        self.edges.append(edge)
        self.pair_edge[(u, v) if u <= v else (v, u)] = edge
        self.adj[u][v] = edge
        self.adj[v][u] = edge

        if edge.relation == "timing":
            self.adj_directed[u].add(v)
            self._add_node_timing(u, v)
            self._add_timing_boundary(ru, rv, edge)
        elif edge.relation == "divergent":
            self._add_divergent_boundary(ru, rv, edge)
        elif edge.relation == "cooccurring":
            self._merge_clusters(u, v)
