# graph_ops/parsing.py

import os
from typing import List, Tuple, Dict, Optional

from models import Edge


def _first_present(kv: Dict[str, str], keys: List[str]) -> Optional[str]:
    for k in keys:
        if k in kv:
            return kv[k]
    return None


def _normalize_haplotype_tag(raw: Optional[str]) -> Optional[str]:
    if raw is None:
        return None
    s = str(raw).strip().upper()
    if not s:
        return None

    alias = {
        "1": "HP1",
        "HP1": "HP1",
        "H1": "HP1",
        "HAP1": "HP1",
        "HAPLOTYPE1": "HP1",
        "2": "HP2",
        "HP2": "HP2",
        "H2": "HP2",
        "HAP2": "HP2",
        "HAPLOTYPE2": "HP2",
        "UNK": "UNKNOWN",
        "UNKNOWN": "UNKNOWN",
        "NA": "UNKNOWN",
        ".": "UNKNOWN",
    }
    return alias.get(s, s)


def parse_haplotype_tags(kv: Dict[str, str]) -> Tuple[Optional[str], Optional[str]]:
    """
    Parse optional haplotype tags for the two SNVs in a row.
    Supports several key naming conventions for forward compatibility.
    """
    hap1_raw = _first_present(
        kv,
        ["HP1", "HAP1", "HAPLOTYPE1", "SNP1_HP", "SNP1_HAP", "SNP1_HAPLOTYPE"],
    )
    hap2_raw = _first_present(
        kv,
        ["HP2", "HAP2", "HAPLOTYPE2", "SNP2_HP", "SNP2_HAP", "SNP2_HAPLOTYPE"],
    )
    return _normalize_haplotype_tag(hap1_raw), _normalize_haplotype_tag(hap2_raw)


def _normalize_hp_reads(raw: Optional[str]) -> Optional[str]:
    if raw is None:
        return None
    s = str(raw).strip()
    if not s:
        return None
    su = s.upper()
    if su in {"UNK", "UNKNOWN", "NA", "."}:
        return "UNKNOWN"
    if "/" in s:
        a, b = s.split("/", 1)
        try:
            ai = int(a)
            bi = int(b)
            if ai >= 0 and bi >= 0:
                return f"{ai}/{bi}"
        except ValueError:
            pass
    return s


def parse_hp_reads(kv: Dict[str, str]) -> Tuple[Optional[str], Optional[str]]:
    hp1_raw = _first_present(
        kv,
        ["HP_READS1", "SNP1_HP_READS", "HPREADS1", "SNP1_HPREADS"],
    )
    hp2_raw = _first_present(
        kv,
        ["HP_READS2", "SNP2_HP_READS", "HPREADS2", "SNP2_HPREADS"],
    )
    return _normalize_hp_reads(hp1_raw), _normalize_hp_reads(hp2_raw)


def parse_snv_line(line: str) -> Tuple[str, int, int, Dict[str, str]]:
    """
    Parse a line like:

      chr7 119138007 119143860 VAF1=... VAF2=... ALT_ALT=... ALT_REF=...
      REF_ALT=... REF_REF=... TOTAL=... RELIABILITY=... BEST_SCORE=... MARGIN=...

    Returns (chrom, pos1, pos2, kv_dict) where kv_dict contains the key=value tokens.
    """
    parts = line.strip().split()
    if len(parts) < 3:
        raise ValueError(f"Line too short: {line!r}")

    chrom = parts[0]
    pos1 = int(parts[1])
    pos2 = int(parts[2])

    kv: Dict[str, str] = {}
    for token in parts[3:]:
        if "=" not in token:
            continue
        key, val = token.split("=", 1)
        kv[key] = val

    return chrom, pos1, pos2, kv


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

    # Map file suffix → (relation, loss_flag, orientation)
    # orientation is only meaningful for timing; for undirected, use "undirected".
    suffixes: Dict[str, Tuple[str, bool, str]] = {
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

                chrom, pos1, pos2, kv = parse_snv_line(line)

                # Read counts (default to 0 if missing for any reason)
                alt_alt = int(kv.get("ALT_ALT", "0"))
                alt_ref = int(kv.get("ALT_REF", "0"))
                ref_alt = int(kv.get("REF_ALT", "0"))
                ref_ref = int(kv.get("REF_REF", "0"))

                # VAFs
                vaf1 = float(kv.get("VAF1", "0.0"))
                vaf2 = float(kv.get("VAF2", "0.0"))

                # Reliability and scoring info from C++
                reliability = float(kv.get("RELIABILITY", "0.0"))
                best_score = float(kv.get("BEST_SCORE", "0.0"))
                margin = float(kv.get("MARGIN", "0.0"))
                hap1, hap2 = parse_haplotype_tags(kv)
                hp_reads1, hp_reads2 = parse_hp_reads(kv)

                # For timing relations, orient u->v using snp1_before_snp2 vs snp2_before_snp1
                if relation == "timing":
                    if orientation == "1_before_2":
                        u, v = pos1, pos2
                        vaf_u, vaf_v = vaf1, vaf2
                        hap_u, hap_v = hap1, hap2
                        hp_reads_u, hp_reads_v = hp_reads1, hp_reads2
                    elif orientation == "2_before_1":
                        u, v = pos2, pos1
                        vaf_u, vaf_v = vaf2, vaf1
                        alt_ref, ref_alt = ref_alt, alt_ref
                        hap_u, hap_v = hap2, hap1
                        hp_reads_u, hp_reads_v = hp_reads2, hp_reads1
                    else:
                        raise ValueError(f"Unexpected timing orientation: {orientation}")
                else:
                    # Undirected relations: just keep the file order
                    u, v = pos1, pos2
                    vaf_u, vaf_v = vaf1, vaf2
                    hap_u, hap_v = hap1, hap2
                    hp_reads_u, hp_reads_v = hp_reads1, hp_reads2

                edges.append(
                    Edge(
                        chrom=chrom,
                        u=u,
                        v=v,
                        relation=relation,
                        loss=loss,
                        reliability=reliability,
                        alt_alt=alt_alt,
                        alt_ref=alt_ref,
                        ref_alt=ref_alt,
                        ref_ref=ref_ref,
                        vaf_u=vaf_u,
                        vaf_v=vaf_v,
                        best_score=best_score,
                        margin=margin,
                        source_file=path,
                        hap_u=hap_u,
                        hap_v=hap_v,
                        hp_reads_u=hp_reads_u,
                        hp_reads_v=hp_reads_v,
                    )
                )

    # Sort edges by reliability descending so the most reliable are considered first.
    edges.sort(key=lambda e: e.reliability, reverse=True)
    return edges


def find_chunk_bases(chunk_dir: str) -> List[str]:
    """
    Discover chunk 'base' paths in chunk_dir, i.e. the original

        chunk_chr_start_end.txt

    file names (without the .txt). We skip the per-relation files like
    *_cooccurring.txt, *_divergent.txt, etc.
    """
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
        # Only process canonical chunk comparison files produced by the C++ pipeline.
        # This intentionally excludes sidecar files (e.g. macOS ._ AppleDouble files)
        # and any other .txt artifacts that do not start with the required prefix.
        if not fname.startswith("chunk_"):
            continue
        if not fname.endswith(".txt"):
            continue
        full = os.path.join(chunk_dir, fname)

        if any(fname.endswith(suf) for suf in relation_suffixes):
            continue

        # This should be the original chunk comparison file; strip '.txt' to get base.
        base = full[: -len(".txt")]
        bases.append(base)

    bases.sort()
    return bases


def parse_chunk_base(base_path: str) -> Tuple[str, int, int]:
    """Return ``(chromosome, core_start, core_end)`` for a chunk base path."""
    basename = os.path.basename(base_path)
    if not basename.startswith("chunk_"):
        raise ValueError(f"Not a canonical chunk base: {base_path}")

    payload = basename[len("chunk_"):]
    try:
        chrom, start_text, end_text = payload.rsplit("_", 2)
        start = int(start_text)
        end = int(end_text)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"Cannot parse chunk base: {base_path}") from exc

    if not chrom or start > end:
        raise ValueError(f"Invalid chunk base: {base_path}")
    return chrom, start, end


def group_chunk_bases_by_chromosome(
    bases: List[str],
) -> Dict[str, List[Tuple[str, int, int]]]:
    """Group canonical chunk bases by chromosome and genomic core order."""
    groups: Dict[str, List[Tuple[str, int, int]]] = {}
    for base in bases:
        chrom, start, end = parse_chunk_base(base)
        groups.setdefault(chrom, []).append((base, start, end))

    for records in groups.values():
        records.sort(key=lambda rec: (rec[1], rec[2], rec[0]))
    return dict(sorted(groups.items()))
