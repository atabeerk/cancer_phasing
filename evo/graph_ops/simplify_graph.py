from __future__ import annotations

import copy
import glob
import hashlib
import json
import os
from collections import defaultdict


KEEP_RELATIONS = {"timing", "divergent"}
CROSS_HAPLOTYPE_PAIR = {"HP1", "HP2"}

NODE_HEIGHT = 35.0
VAF_BIN_WIDTH = 0.1
BIN_HEIGHT = 2 * NODE_HEIGHT
N_VAF_BINS = int(round(1.0 / VAF_BIN_WIDTH))
VAF_BAND_HEIGHT = N_VAF_BINS * BIN_HEIGHT
HAPLOTYPE_BAND_GAP = 4 * NODE_HEIGHT
HAPLOTYPE_BAND_INDEX = {"HP1": 0, "HP2": 2}

TOP_NODE_X_GAP = 80.0
TOP_CHROMOSOME_GAP = 220.0


class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, value):
        self.parent.setdefault(value, value)
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left, right):
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            self.parent[left_root] = right_root


def chromosome_sort_key(chrom):
    value = str(chrom or "")
    token = value[3:] if value.lower().startswith("chr") else value
    upper = token.upper()
    if upper.isdigit():
        return (0, int(upper), "")
    standard = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if upper in standard:
        return (0, standard[upper], "")
    return (1, 0, upper)


def _position(node):
    return int(node["data"]["position"])


def _node_sort_key(node):
    data = node["data"]
    return (
        chromosome_sort_key(data.get("chrom")),
        _position(node),
        str(data["id"]),
    )


def vaf_bin_y(node_id, vaf):
    bin_index = min(int(float(vaf) / VAF_BIN_WIDTH), N_VAF_BINS - 1)
    bin_index = (N_VAF_BINS - 1) - bin_index
    digest = int(hashlib.md5(str(node_id).encode()).hexdigest(), 16)
    jitter = ((digest % 1000) / 1000.0 - 0.5) * NODE_HEIGHT
    return (bin_index + 0.5) * BIN_HEIGHT + jitter


def node_y(node_id, vaf, haplotype):
    band_index = HAPLOTYPE_BAND_INDEX.get(haplotype, 1)
    return (
        band_index * (VAF_BAND_HEIGHT + HAPLOTYPE_BAND_GAP)
        + vaf_bin_y(node_id, vaf)
    )


def _index_nodes(nodes):
    indexed = {}
    for node in nodes:
        node_id = str(node["data"]["id"])
        if node_id in indexed:
            raise ValueError(f"Duplicate node id in chromosome graph: {node_id}")
        node["data"]["id"] = node_id
        indexed[node_id] = node
    return indexed


def _normalize_edges(edges):
    for edge in edges:
        data = edge["data"]
        data["id"] = str(data["id"])
        data["source"] = str(data["source"])
        data["target"] = str(data["target"])


def remove_cross_haplotype_edges(nodes_by_id, edges):
    kept = []
    removed = 0
    for edge in edges:
        data = edge["data"]
        source = nodes_by_id.get(data["source"])
        target = nodes_by_id.get(data["target"])
        if source is None or target is None:
            raise ValueError(
                f"Edge {data['id']} references a node absent from the graph"
            )
        source_hp = source["data"].get("haplotype")
        target_hp = target["data"].get("haplotype")
        if {source_hp, target_hp} == CROSS_HAPLOTYPE_PAIR:
            removed += 1
            continue
        kept.append(edge)
    return kept, removed


def connected_components(node_ids, edges):
    union_find = UnionFind()
    for node_id in node_ids:
        union_find.find(node_id)
    for edge in edges:
        data = edge["data"]
        source = data["source"]
        target = data["target"]
        if source in union_find.parent and target in union_find.parent:
            union_find.union(source, target)

    component_nodes = defaultdict(list)
    for node_id in node_ids:
        component_nodes[union_find.find(node_id)].append(node_id)

    component_edges = defaultdict(list)
    for edge in edges:
        source = edge["data"]["source"]
        if source in union_find.parent:
            component_edges[union_find.find(source)].append(edge)

    return component_nodes, component_edges


def cooccurring_clusters(node_ids, edges):
    cooccurring = [
        edge for edge in edges
        if edge["data"].get("relation") == "cooccurring"
    ]
    union_find = UnionFind()
    for node_id in node_ids:
        union_find.find(node_id)
    for edge in cooccurring:
        union_find.union(edge["data"]["source"], edge["data"]["target"])

    node_to_cluster = {
        node_id: union_find.find(node_id)
        for node_id in node_ids
    }
    cluster_members = defaultdict(list)
    for node_id in node_ids:
        cluster_members[node_to_cluster[node_id]].append(node_id)
    return node_to_cluster, cluster_members


def simplify_cooccurring_edges(nodes_by_id, edges, cluster_members):
    cooccurring = [
        edge for edge in edges
        if edge["data"].get("relation") == "cooccurring"
    ]
    other = [
        edge for edge in edges
        if edge["data"].get("relation") != "cooccurring"
    ]

    adjacency = defaultdict(dict)
    for edge in cooccurring:
        source = edge["data"]["source"]
        target = edge["data"]["target"]
        adjacency[source][target] = edge
        adjacency[target][source] = edge

    positions = {
        node_id: _position(node)
        for node_id, node in nodes_by_id.items()
    }
    kept = {}

    for members in cluster_members.values():
        for node_id in members:
            node_position = positions[node_id]
            right = [
                (positions[other_id], other_id)
                for other_id in adjacency[node_id]
                if positions[other_id] > node_position
            ]
            left = [
                (positions[other_id], other_id)
                for other_id in adjacency[node_id]
                if positions[other_id] < node_position
            ]
            if right:
                _position_value, nearest = min(right, key=lambda item: item[0])
                edge = adjacency[node_id][nearest]
                kept[edge["data"]["id"]] = edge
            if left:
                _position_value, nearest = max(left, key=lambda item: item[0])
                edge = adjacency[node_id][nearest]
                kept[edge["data"]["id"]] = edge

    return list(kept.values()) + other


def fragmented_cooccurring_clusters(
    node_ids,
    original_cluster_members,
    simplified_cooccurring_edges,
):
    union_find = UnionFind()
    for node_id in node_ids:
        union_find.find(node_id)
    for edge in simplified_cooccurring_edges:
        union_find.union(edge["data"]["source"], edge["data"]["target"])

    fragmented = []
    for root, members in original_cluster_members.items():
        if len(members) < 2:
            continue
        pieces = defaultdict(list)
        for node_id in members:
            pieces[union_find.find(node_id)].append(node_id)
        if len(pieces) > 1:
            fragmented.append((root, list(pieces.values())))
    return fragmented


def informative_read_count(edge):
    counts = str(edge["data"].get("read_counts", "")).split("/")
    if len(counts) != 4:
        raise ValueError(
            f"Edge {edge['data'].get('id')} has invalid read_counts: "
            f"{edge['data'].get('read_counts')!r}"
        )
    return sum(int(value) for value in counts[:3])


def deduplicate_relation_edges(node_ids, edges, node_to_cluster, relation):
    relation_edges = [
        edge for edge in edges
        if edge["data"].get("relation") == relation
    ]
    other_edges = [
        edge for edge in edges
        if edge["data"].get("relation") != relation
    ]

    incident = defaultdict(list)
    for edge in relation_edges:
        source = edge["data"]["source"]
        target = edge["data"]["target"]
        incident[source].append((node_to_cluster[target], edge))
        incident[target].append((node_to_cluster[source], edge))

    drop_ids = set()
    for node_id in node_ids:
        by_cluster = defaultdict(list)
        for cluster_id, edge in incident[node_id]:
            by_cluster[cluster_id].append(edge)
        for group in by_cluster.values():
            if len(group) < 2:
                continue
            group.sort(
                key=lambda edge: (
                    informative_read_count(edge),
                    float(edge["data"].get("reliability", 0.0)),
                    str(edge["data"]["id"]),
                ),
                reverse=True,
            )
            for edge in group[1:]:
                drop_ids.add(edge["data"]["id"])

    survivors = [
        edge for edge in relation_edges
        if edge["data"]["id"] not in drop_ids
    ]
    return survivors + other_edges, len(drop_ids)


def build_simplified_graph(graph, source_name, name=None):
    nodes = copy.deepcopy(graph.get("elements", {}).get("nodes", []))
    edges = copy.deepcopy(graph.get("elements", {}).get("edges", []))
    nodes_by_id = _index_nodes(nodes)
    _normalize_edges(edges)

    for node_id, node in nodes_by_id.items():
        data = node["data"]
        data["source_chunk"] = source_name
        gene_name = data.get("gene_name")
        data["label"] = (
            f"{data['position']}\n{gene_name}"
            if gene_name
            else str(data["position"])
        )
        if "position" in node:
            node["position"]["y"] = node_y(
                node_id,
                float(data.get("vaf", 0.0) or 0.0),
                data.get("haplotype"),
            )

    edges, removed_cross_haplotype = remove_cross_haplotype_edges(
        nodes_by_id,
        edges,
    )

    node_ids = list(nodes_by_id)
    component_nodes, _component_edges = connected_components(node_ids, edges)
    for component_id, members in component_nodes.items():
        for node_id in members:
            nodes_by_id[node_id]["data"]["component_id"] = component_id

    node_to_cluster, cluster_members = cooccurring_clusters(node_ids, edges)
    simplified_edges = simplify_cooccurring_edges(
        nodes_by_id,
        edges,
        cluster_members,
    )

    fragmented = fragmented_cooccurring_clusters(
        node_ids,
        cluster_members,
        [
            edge for edge in simplified_edges
            if edge["data"].get("relation") == "cooccurring"
        ],
    )
    for _cluster_root, pieces in fragmented:
        pieces.sort(
            key=lambda members: min(
                _position(nodes_by_id[node_id])
                for node_id in members
            )
        )
        for index in range(len(pieces) - 1):
            source = max(
                pieces[index],
                key=lambda node_id: _position(nodes_by_id[node_id]),
            )
            target = min(
                pieces[index + 1],
                key=lambda node_id: _position(nodes_by_id[node_id]),
            )
            simplified_edges.append({
                "data": {
                    "id": f"{source}-{target}-inferred",
                    "source": source,
                    "target": target,
                    "relation": "cooccurring",
                    "loss": False,
                    "reliability": 0.0,
                    "directed": False,
                    "read_counts": "0/0/0/0",
                    "label": "inferred bridge (no direct evidence)",
                    "inferred": True,
                }
            })

    simplified_edges, timing_dropped = deduplicate_relation_edges(
        node_ids,
        simplified_edges,
        node_to_cluster,
        "timing",
    )
    simplified_edges, divergent_dropped = deduplicate_relation_edges(
        node_ids,
        simplified_edges,
        node_to_cluster,
        "divergent",
    )

    for node_id in node_ids:
        nodes_by_id[node_id]["data"]["cooc_cluster_id"] = (
            node_to_cluster[node_id]
        )

    nodes.sort(key=_node_sort_key)
    output_name = name or os.path.splitext(source_name)[0] + "_simplified"
    simplified = {
        "data": {"name": output_name},
        "elements": {
            "nodes": nodes,
            "edges": simplified_edges,
        },
    }
    stats = {
        "removed_cross_haplotype_edges": removed_cross_haplotype,
        "cooccurring_edges_before": sum(
            edge["data"].get("relation") == "cooccurring"
            for edge in edges
        ),
        "cooccurring_edges_after": sum(
            edge["data"].get("relation") == "cooccurring"
            for edge in simplified_edges
        ),
        "timing_edges_dropped": timing_dropped,
        "divergent_edges_dropped": divergent_dropped,
        "fragmented_cooccurring_clusters": len(fragmented),
    }
    return simplified, stats


def write_simplified_graph(input_path, output_path):
    with open(input_path, "r", encoding="utf-8") as handle:
        graph = json.load(handle)

    output_name = os.path.splitext(os.path.basename(output_path))[0]
    simplified, stats = build_simplified_graph(
        graph,
        source_name=os.path.basename(input_path),
        name=output_name,
    )
    with open(output_path, "w", encoding="utf-8") as handle:
        json.dump(simplified, handle, indent=2)
    return stats


def _component_records(graph, fallback_chrom):
    nodes = graph.get("elements", {}).get("nodes", [])
    edges = graph.get("elements", {}).get("edges", [])
    nodes_by_id = {
        str(node["data"]["id"]): node
        for node in nodes
    }
    _normalize_edges(edges)
    component_nodes, component_edges = connected_components(
        list(nodes_by_id),
        edges,
    )

    records = []
    for component_id, node_ids in component_nodes.items():
        component_edge_list = component_edges.get(component_id, [])
        relations = [
            edge["data"].get("relation")
            for edge in component_edge_list
        ]
        chrom = fallback_chrom
        if node_ids:
            chrom = nodes_by_id[node_ids[0]]["data"].get(
                "chrom",
                fallback_chrom,
            )
        records.append({
            "chrom": chrom,
            "component_id": component_id,
            "nodes": [nodes_by_id[node_id] for node_id in node_ids],
            "edges": component_edge_list,
            "timing_edges": sum(
                relation == "timing" for relation in relations
            ),
            "interesting_edges": sum(
                relation in KEEP_RELATIONS for relation in relations
            ),
            "total_edges": len(component_edge_list),
            "kept": bool(set(relations) & KEEP_RELATIONS),
            "min_position": min(
                _position(nodes_by_id[node_id])
                for node_id in node_ids
            ),
        })
    return records


def _namespace_component(component, rank):
    chrom = str(component["chrom"])
    node_id_map = {}
    output_nodes = []

    for node in component["nodes"]:
        copied = copy.deepcopy(node)
        data = copied["data"]
        old_id = str(data["id"])
        new_id = f"{chrom}:{old_id}"
        node_id_map[old_id] = new_id
        data["id"] = new_id
        data["chrom"] = chrom
        data["component_id"] = (
            f"{chrom}:{component['component_id']}"
        )
        if "cooc_cluster_id" in data:
            data["cooc_cluster_id"] = (
                f"{chrom}:{data['cooc_cluster_id']}"
            )
        data["component_rank"] = rank
        data["component_timing_edges"] = component["timing_edges"]
        data["component_total_edges"] = component["total_edges"]
        output_nodes.append(copied)

    output_edges = []
    for edge in component["edges"]:
        copied = copy.deepcopy(edge)
        data = copied["data"]
        data["id"] = f"{chrom}:{data['id']}"
        data["source"] = node_id_map[str(data["source"])]
        data["target"] = node_id_map[str(data["target"])]
        output_edges.append(copied)

    return output_nodes, output_edges


def assemble_top_components(components, name):
    nodes = []
    edges = []
    for rank, component in enumerate(components, start=1):
        component_nodes, component_edges = _namespace_component(
            component,
            rank,
        )
        nodes.extend(component_nodes)
        edges.extend(component_edges)

    nodes.sort(key=_node_sort_key)
    previous_chrom = None
    next_x = 0.0
    node_order = {}
    for order, node in enumerate(nodes):
        chrom = node["data"].get("chrom")
        if previous_chrom is not None and chrom != previous_chrom:
            next_x += TOP_CHROMOSOME_GAP
        node.setdefault("position", {})["x"] = next_x
        next_x += TOP_NODE_X_GAP
        previous_chrom = chrom
        node_order[node["data"]["id"]] = order

    edges.sort(
        key=lambda edge: (
            min(
                node_order[edge["data"]["source"]],
                node_order[edge["data"]["target"]],
            ),
            max(
                node_order[edge["data"]["source"]],
                node_order[edge["data"]["target"]],
            ),
            str(edge["data"]["id"]),
        )
    )
    return {
        "data": {"name": name},
        "elements": {"nodes": nodes, "edges": edges},
    }


def _ranking_key(component, metric):
    return (
        -component[metric],
        chromosome_sort_key(component["chrom"]),
        component["min_position"],
        str(component["component_id"]),
    )


def write_topn_graphs(main_output_dir, top_n, chromosomes):
    if top_n < 1:
        raise ValueError("top_n must be >= 1")

    components = []
    for chromosome in chromosomes:
        chrom = chromosome[0] if isinstance(chromosome, tuple) else chromosome
        input_path = os.path.join(
            main_output_dir,
            f"{chrom}_out",
            "graphs",
            f"{chrom}_simplified.json",
        )
        if not os.path.isfile(input_path):
            continue
        with open(input_path, "r", encoding="utf-8") as handle:
            graph = json.load(handle)
        components.extend(_component_records(graph, chrom))

    eligible = [component for component in components if component["kept"]]
    output_dir = os.path.join(main_output_dir, "topN_graphs")
    os.makedirs(output_dir, exist_ok=True)
    for stale_path in glob.glob(
        os.path.join(output_dir, "top*_components_by_*.json")
    ):
        os.remove(stale_path)

    rankings = (
        ("timing_edges", "timing_edges"),
        ("interesting_edges", "timing_divergent_edges"),
        ("total_edges", "total_edges"),
    )
    output_paths = []
    for metric, file_stub in rankings:
        selected = sorted(
            eligible,
            key=lambda component: _ranking_key(component, metric),
        )[:top_n]
        name = f"top{top_n}_components_by_{file_stub}"
        graph = assemble_top_components(selected, name)
        output_path = os.path.join(output_dir, name + ".json")
        with open(output_path, "w", encoding="utf-8") as handle:
            json.dump(graph, handle, indent=2)
        output_paths.append(output_path)

    return output_paths
