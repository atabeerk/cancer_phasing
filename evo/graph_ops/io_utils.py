from pathlib import Path
from typing import Tuple, Set, Dict, Any

Edge = Tuple[str, str, str]


def read_edges_from_file(file_path: Path, reverse: bool = False) -> Tuple[Set[str], Set[Edge], Dict[Edge, Dict[str, Any]]]:
    """Read edges (chrom, pos1, pos2) with TOTAL and ORIGINAL metadata.
    Returns (nodes, edges, edge_meta).
    """
    nodes = set()
    edges = set()
    edge_meta: Dict[Edge, Dict[str, Any]] = {}

    with open(file_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            chrom, pos1, pos2 = parts[0], parts[1], parts[2]
            if reverse:
                pos1, pos2 = pos2, pos1
            total = None
            original = None
            for token in parts[3:]:
                if token.startswith("TOTAL="):
                    try:
                        total = int(token.split("=", 1)[1])
                    except ValueError:
                        total = None
                elif token.startswith("ORIGINAL="):
                    original = token.split("=", 1)[1]
            edges.add((chrom, pos1, pos2))
            nodes.update([pos1, pos2])
            edge_meta[(chrom, pos1, pos2)] = {"total": total, "original": original}

    return nodes, edges, edge_meta


def parse_pileup_depth(pileup_file: Path) -> Dict[str, int]:
    """Parse a pileup-like file to map position -> depth.
    Returns a dict with string keys for positions to match how nodes are stored.
    If pileup_file doesn't exist, returns an empty dict.
    """
    depths: Dict[str, int] = {}
    if not pileup_file or not pileup_file.exists():
        return depths

    with open(pileup_file) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            try:
                pos = int(parts[1])
                depth = int(parts[3])
            except ValueError:
                continue
            depths[str(pos)] = depth
    return depths


def write_component_statistics(out_dir: Path, base: str, component_stats: list) -> None:
    """Append component stats to a single 'component_statistics.txt' file in parent of out_dir."""
    out_dir = Path(out_dir)
    stats_path = out_dir.parent.resolve() / "component_statistics.txt"

    write_header = not stats_path.exists()

    with open(stats_path, "a") as stats_out:
        if write_header:
            stats_out.write(
                "base\tcomponent_id\tnum_nodes\tspan_bp\thaplotypes\tmulti_node_haplotypes\tnodes\n"
            )

        for comp in component_stats:
            stats_out.write(
                f"{base}\t{comp['component_id']}\t{len(comp['nodes'])}\t"
                f"{comp['span_bp']}\t{comp['haplotypes']}\t{comp['multi_node_haplotypes']}\t"
                f"{','.join(sorted(comp['nodes']))}\n"
            )
