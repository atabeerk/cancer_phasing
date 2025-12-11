import os
from typing import List, Tuple, Dict

from models import Edge


def parse_snv_line(line: str) -> Tuple[int, int, float, Dict[str, str]]:
    """
    Parse a line of the per-relationship C++ output, e.g.:

      chr1 123 456 VAF1=... VAF2=... ... RELIABILITY=0.82 BEST_SCORE=... MARGIN=...

    Returns (pos1, pos2, reliability, kv_dict).
    """
    parts = line.strip().split()
    if len(parts) < 3:
        raise ValueError(f"Line too short: {line!r}")

    # First 3 columns: chr, pos1, pos2
    # We ignore chr; the chunk is per-chromosome.
    _chr = parts[0]
    pos1 = int(parts[1])
    pos2 = int(parts[2])

    kv: Dict[str, str] = {}
    for token in parts[3:]:
        if "=" not in token:
            continue
        key, val = token.split("=", 1)
        kv[key] = val

    reliability = float(kv.get("RELIABILITY", "0.0"))
    return pos1, pos2, reliability, kv


def load_edges_from_base(base_path: str) -> List[Edge]:
    """
    Load edges for a single chunk, given the 'base' used in C++:

        base = original_chunk_filename_without_extension

    This function looks for:

        base + "_cooccurring.txt"
        base + "_cooccurring_loss.txt"
        base + "_divergent.txt"
        base + "_snp1_before_snp2.txt"
        base + "_snp2_before_snp1.txt"
        base + "_snp1_before_snp2_loss.txt"
        base + "_snp2_before_snp1_loss.txt"

    and ignores base + "_errors.txt".
    """
    suffixes = {
        "_cooccurring.txt": ("cooccurring", False, "undirected"),
        "_cooccurring_loss.txt": ("cooccurring", True, "undirected"),
        "_divergent.txt": ("divergent", False, "undirected"),
        "_snp1_before_snp2.txt": ("timing", False, "1_before_2"),
        "_snp2_before_snp1.txt": ("timing", False, "2_before_1"),
        "_snp1_before_snp2_loss.txt": ("timing", True, "1_before_2"),
        "_snp2_before_snp1_loss.txt": ("timing", True, "2_before_1"),
        # "_errors.txt" is intentionally not used
    }

    edges: List[Edge] = []

    for suffix, (relation, loss, orientation) in suffixes.items():
        path = base_path + suffix
        if not os.path.exists(path):
            continue

        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                pos1, pos2, reliability, _kv = parse_snv_line(line)

                if relation == "timing":
                    # Encode direction as u -> v (u is earlier in time)
                    if orientation == "1_before_2":
                        u, v = pos1, pos2
                    elif orientation == "2_before_1":
                        u, v = pos2, pos1
                    else:
                        raise ValueError(f"Unexpected timing orientation: {orientation}")
                else:
                    # For undirected relationships, orientation doesn't matter
                    u, v = pos1, pos2

                edges.append(
                    Edge(
                        u=u,
                        v=v,
                        relation=relation,
                        loss=loss,
                        reliability=reliability,
                        source_file=path,
                    )
                )

    # Sort by reliability descending so the most reliable edges are considered first
    edges.sort(key=lambda e: e.reliability, reverse=True)
    return edges


def find_chunk_bases(chunk_dir: str) -> List[str]:
    """
    Discover chunk 'base' paths in chunk_dir, i.e. the original

        chunk_chr_start_end.txt

    file names (without the .txt). We skip the per-relation files like
    *_cooccurring.txt, *_divergent.txt, etc.
    """
    # These suffixes belong to per-relation outputs; we don't want them as bases.
    relation_suffixes = [
        "_cooccurring.txt",
        "_cooccurring_loss.txt",
        "_divergent.txt",
        "_snp1_before_snp2.txt",
        "_snp2_before_snp1.txt",
        "_snp1_before_snp2_loss.txt",
        "_snp2_before_snp1_loss.txt",
        "_errors.txt",
    ]

    bases: List[str] = []
    for fname in os.listdir(chunk_dir):
        if not fname.endswith(".txt"):
            continue
        full = os.path.join(chunk_dir, fname)

        if any(fname.endswith(suf) for suf in relation_suffixes):
            # it's one of the per-relation outputs, skip
            continue

        # This should be the original chunk comparison file
        base = full[: -len(".txt")]
        bases.append(base)

    bases.sort()
    return bases
