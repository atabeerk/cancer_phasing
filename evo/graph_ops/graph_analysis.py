from typing import Tuple, Dict, List
import networkx as nx
from models import GraphData


def compute_component_statistics(graph: GraphData) -> Tuple[List[Dict], Dict[str, int]]:
    """Compute components, spans and haplotype counts.

    Returns (component_stats_list, node_to_component_map).
    """
    nodes = graph.nodes
    directed_edges = graph.directed
    directed_loss_edges = graph.directed_loss
    cooccurring_edges = graph.cooccurring
    cooccurring_loss_edges = graph.cooccurring_loss
    divergent_edges = graph.divergent

    G_all = nx.Graph()
    G_all.add_nodes_from(nodes)
    G_all.add_edges_from([(s, t) for _, s, t in directed_edges])
    G_all.add_edges_from([(s, t) for _, s, t in directed_loss_edges])
    G_all.add_edges_from([(s, t) for _, s, t in cooccurring_edges])
    G_all.add_edges_from([(s, t) for _, s, t in cooccurring_loss_edges])
    G_all.add_edges_from([(s, t) for _, s, t in divergent_edges])

    component_stats = []
    node_to_component: Dict[str, int] = {}
    for comp_id, comp_nodes in enumerate(nx.connected_components(G_all)):
        comp_nodes = list(comp_nodes)
        positions = sorted(int(n) for n in comp_nodes if n.isdigit())
        span = max(positions) - min(positions) if positions else 0

        # Build cooccurring subgraph limited to this component
        G_co = nx.Graph()
        for _, s, t in cooccurring_edges | cooccurring_loss_edges:
            if s in comp_nodes and t in comp_nodes:
                G_co.add_edge(s, t)

        multi_node_haplotypes = nx.number_connected_components(G_co)
        # ensure all nodes present (singleton nodes)
        for n in comp_nodes:
            if n not in G_co:
                G_co.add_node(n)

        haplotypes = nx.number_connected_components(G_co)

        for n in comp_nodes:
            node_to_component[n] = comp_id

        component_stats.append({
            "component_id": comp_id,
            "nodes": comp_nodes,
            "span_bp": span,
            "haplotypes": haplotypes,  # includes singletons
            "multi_node_haplotypes": multi_node_haplotypes  # only multi-node clusters
        })

    return component_stats, node_to_component
