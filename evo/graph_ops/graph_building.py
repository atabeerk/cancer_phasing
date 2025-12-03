# your_package/graph_building.py
from pathlib import Path
from collections import defaultdict
from typing import Dict
from io_utils import read_edges_from_file
from models import GraphData


def build_graphs(chunk_dir: str) -> Dict[str, GraphData]:
    """Group files by chunk and type; build GraphData objects."""
    chunk_dir = Path(chunk_dir)
    graphs: Dict[str, GraphData] = {}
    grouped = defaultdict(dict)

    for f in chunk_dir.glob("chunk_*_*.txt"):
        name = f.stem
        if "_snp1_before_snp2_loss" in name:
            grouped[name.replace("_snp1_before_snp2_loss", "")]["forward_loss"] = f
        elif "_snp2_before_snp1_loss" in name:
            grouped[name.replace("_snp2_before_snp1_loss", "")]["reverse_loss"] = f
        elif "_snp1_before_snp2" in name:
            grouped[name.replace("_snp1_before_snp2", "")]["forward"] = f
        elif "_snp2_before_snp1" in name:
            grouped[name.replace("_snp2_before_snp1", "")]["reverse"] = f
        elif "_cooccurring_loss" in name:
            grouped[name.replace("_cooccurring_loss", "")]["cooccurring_loss"] = f
        elif "_cooccurring" in name:
            grouped[name.replace("_cooccurring", "")]["cooccurring"] = f
        elif "_divergent" in name:
            grouped[name.replace("_divergent", "")]["divergent"] = f

    for base, files in grouped.items():
        nodes, directed_edges, directed_loss_edges = set(), set(), set()
        cooccurring_edges, cooccurring_loss_edges, divergent_edges = set(), set(), set()
        edge_meta = {}

        # Directed edges
        for key, rev in [("forward", False), ("reverse", True)]:
            if key in files:
                n, e, meta = read_edges_from_file(files[key], reverse=rev)
                nodes |= n
                directed_edges |= e
                edge_meta.update(meta)

        # Directed loss
        for key, rev in [("forward_loss", False), ("reverse_loss", True)]:
            if key in files:
                n, e, meta = read_edges_from_file(files[key], reverse=rev)
                nodes |= n
                directed_loss_edges |= e
                edge_meta.update(meta)

        # Cooccurring
        if "cooccurring" in files:
            n, e, meta = read_edges_from_file(files["cooccurring"])
            nodes |= n
            cooccurring_edges |= e
            edge_meta.update(meta)

        # Cooccurring loss
        if "cooccurring_loss" in files:
            n, e, meta = read_edges_from_file(files["cooccurring_loss"])
            nodes |= n
            cooccurring_loss_edges |= e
            edge_meta.update(meta)

        # Divergent (only if both nodes exist)
        if "divergent" in files:
            n, e, meta = read_edges_from_file(files["divergent"])
            for chrom, s, t in e:
                if s in nodes and t in nodes:
                    divergent_edges.add((chrom, s, t))
                    edge_meta[(chrom, s, t)] = meta.get((chrom, s, t), {})

        graphs[base] = GraphData.from_parts(
            base=base,
            nodes=nodes,
            directed=directed_edges,
            directed_loss=directed_loss_edges,
            cooccurring=cooccurring_edges,
            cooccurring_loss=cooccurring_loss_edges,
            divergent=divergent_edges,
            edge_meta=edge_meta
        )


    return graphs
