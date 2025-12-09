import networkx as nx
from pathlib import Path
import json
import csv
from typing import Set, Tuple, Dict, List

Edge = Tuple[str, str, str]  # (chrom, source, target)


def conflict_free_clusters(
        nodes: Set[str],
        cooccurring_edges: Set[Edge],
        cooccurring_loss_edges: Set[Edge],
        directed_edges: Set[Edge],
        directed_loss_edges: Set[Edge],
        divergent_edges: Set[Edge]
    ) -> List[Set[str]]:
    """
    Build conflict-free clusters of nodes using cooccurrence as connectivity
    and removing nodes that create directed/divergent conflicts.

    IMPORTANT:
    - When a node is removed from a cluster, the cluster may fracture.
      Each fractured subcluster is re-evaluated independently.
    - Removed nodes are NOT discarded; after processing all clusters they
      are used to build new clusters (or singleton clusters).
    """

    # Build full cooccurrence graph
    G_co = nx.Graph()
    G_co.add_nodes_from(nodes)
    G_co.add_edges_from([(s, t) for _, s, t in cooccurring_edges])
    G_co.add_edges_from([(s, t) for _, s, t in cooccurring_loss_edges])

    # Precompute all conflict edges
    conflict_edges = set()
    for edge_set in (directed_edges, directed_loss_edges, divergent_edges):
        for _, s, t in edge_set:
            conflict_edges.add((s, t))

    # Utility: compute internal conflicts in a cluster
    def get_internal_conflicts(cluster: Set[str]):
        return [(s, t) for (s, t) in conflict_edges if s in cluster and t in cluster]

    # Utility: choose a node to remove
    def choose_node_to_remove(cluster: Set[str], internal_conflicts: List[Tuple[str, str]]):
        conflict_degree = {n: 0 for n in cluster}
        coocc_degree = {n: 0 for n in cluster}

        for s, t in internal_conflicts:
            if s in conflict_degree: conflict_degree[s] += 1
            if t in conflict_degree: conflict_degree[t] += 1

        for _, s, t in cooccurring_edges | cooccurring_loss_edges:
            if s in coocc_degree and t in coocc_degree:
                coocc_degree[s] += 1
                coocc_degree[t] += 1

        # highest conflict degree
        max_c = max(conflict_degree.values())
        candidates = [n for n in cluster if conflict_degree[n] == max_c]

        # tie-break: lowest cooccurrence degree
        if len(candidates) > 1:
            min_co = min(coocc_degree[n] for n in candidates)
            candidates = [n for n in candidates if coocc_degree[n] == min_co]

        # deterministic choice
        return sorted(candidates)[0]

    # Step 1: initial clusters = CCs of cooccurrence graph
    initial_clusters = [set(c) for c in nx.connected_components(G_co)]
    print(f"\n[INFO] Found {len(initial_clusters)} initial cooccurrence clusters")

    queue = initial_clusters[:]  # clusters needing evaluation
    final_clusters = []
    removed_nodes = set()

    # Step 2: process clusters until all are conflict-free
    while queue:
        cluster = queue.pop()
        if len(cluster) <= 1:
            final_clusters.append(cluster)
            continue

        internal_conflicts = get_internal_conflicts(cluster)

        if not internal_conflicts:
            # conflict-free cluster
            final_clusters.append(cluster)
            continue

        print(f"  - Internal conflicts detected: {len(internal_conflicts)} edges")
        for s, t in internal_conflicts:
            print(f"     * Conflict edge: {s} -> {t}")

        # remove one problematic node
        node_to_remove = choose_node_to_remove(cluster, internal_conflicts)
        print(f"  - Removing node: {node_to_remove}")
        cluster.remove(node_to_remove)
        removed_nodes.add(node_to_remove)

        # cluster may fracture — recompute connected components
        if cluster:
            subclusters = [set(c) for c in nx.connected_components(G_co.subgraph(cluster))]
            queue.extend(subclusters)

    # Step 3: process removed nodes into new clusters
    # Build cooccurrence graph restricted to removed nodes
    G_removed = nx.Graph()
    G_removed.add_nodes_from(removed_nodes)

    for _, s, t in cooccurring_edges | cooccurring_loss_edges:
        if s in removed_nodes and t in removed_nodes:
            G_removed.add_edge(s, t)

    # connected components among removed nodes
    leftover_clusters = [set(c) for c in nx.connected_components(G_removed)]

    # nodes not in any edges → singleton clusters
    nodes_in_removed_edges = set().union(*leftover_clusters) if leftover_clusters else set()
    singletons = removed_nodes - nodes_in_removed_edges
    leftover_clusters.extend([{n} for n in singletons])

    # add all leftover clusters to final result
    final_clusters.extend(leftover_clusters)

    return final_clusters



def condense_graph(
        nodes: Set[str],
        directed_edges: Set[Edge],
        directed_loss_edges: Set[Edge],
        cooccurring_edges: Set[Edge],
        cooccurring_loss_edges: Set[Edge],
        divergent_edges: Set[Edge],
    ) -> Tuple[List[Dict], Dict[str, Set[Tuple[str, str]]]]:
    """
    Collapse cooccurrence clusters into CNs using conflict-free clustering.
    Returns:
        - cn_list: list of CN dicts with members and span
        - cn_edges: dict of CN-level edges by type
    """
    # Get conflict-free clusters
    clusters = conflict_free_clusters(
        nodes,
        cooccurring_edges,
        cooccurring_loss_edges,
        directed_edges,
        directed_loss_edges,
        divergent_edges
    )

    cn_list, node_to_cn = [], {}
    for cn_id, cluster in enumerate(clusters):
        cn_name = f"CN{cn_id}"
        cn_list.append({
            "id": cn_name,
            "members": list(cluster),
            "n_nodes": len(cluster),
            "span": max(map(int, cluster)) - min(map(int, cluster)) if all(n.isdigit() for n in cluster) else None,
            "mixed_edges": False
        })
        for n in cluster:
            node_to_cn[n] = cn_name

    # Build CN-level edges (only between CNs)
    cn_edges = {"directed": set(), "directed_loss": set(), "divergent": set()}

    def add_cn_edge(edge_type, s, t, chrom):
        src, dst = node_to_cn.get(s), node_to_cn.get(t)
        if not src or not dst or src == dst:
            return
        cn_edges[edge_type].add((chrom, src, dst))

    for chrom, s, t in directed_edges:
        add_cn_edge("directed", s, t, chrom)
    for chrom, s, t in directed_loss_edges:
        add_cn_edge("directed_loss", s, t, chrom)
    for chrom, s, t in divergent_edges:
        add_cn_edge("divergent", s, t, chrom)

    return cn_list, cn_edges


def write_condensed_json_and_tsv(cn_nodes: Dict[str, Set[str]],
                                 cn_edges: List[Tuple[str, str, str]],
                                 cn_flags: Dict[str, bool],
                                 out_base: Path):
    """
    Write condensed CN graph to JSON and TSV.
    """
    elements = {"nodes": [], "edges": []}

    for cn_id, members in cn_nodes.items():
        elements["nodes"].append({
            "data": {"id": cn_id, "label": cn_id,
                     "size": len(members), "mixed_edges": cn_flags.get(cn_id, False)}
        })

    for src, tgt, edge_type in cn_edges:
        elements["edges"].append({
            "data": {"source": src, "target": tgt, "edge_type": edge_type,
                     "directed": edge_type in ("directed", "directed_loss")}
        })

    with open(f"{out_base}_condensed.json", "w") as jf:
        json.dump({"elements": elements}, jf, indent=2)

    with open(f"{out_base}_condensed.stats.tsv", "w", newline="") as tf:
        writer = csv.writer(tf, delimiter="\t")
        writer.writerow(["CN_id", "num_nodes", "span_bp", "nodes_list", "mixed_edges"])
        for cn_id, members in cn_nodes.items():
            positions = sorted(int(n) for n in members if n.isdigit())
            span = max(positions) - min(positions) if positions else 0
            writer.writerow([cn_id, len(members), span, ",".join(sorted(members)), cn_flags.get(cn_id, False)])

    print(f"Condensed JSON/TSV written: {out_base}_condensed.*")
