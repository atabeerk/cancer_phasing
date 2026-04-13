#!/usr/bin/env python3
import argparse
import csv
import gzip
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Set, Tuple, List, Optional, Union


def canonical_chrom(chrom: str) -> str:
    """
    Canonicalize chromosome names so 'chr7' and '7' match.
    Strips leading 'chr' (case-insensitive) and uppercases.
    Normalizes mitochondria naming to 'MT'.
    """
    c = str(chrom).strip()
    c_low = c.lower()
    if c_low.startswith("chr"):
        c = c[3:].strip()

    c_up = c.upper()
    if c_up in {"M", "MT", "MITO"}:
        return "MT"
    return c_up


def vcf_label(vcf_path: Path) -> str:
    """
    Return source VCF label = filename prefix before first dot.
    Examples:
      N10.filtered.sorted.vcf.gz -> N10
      sample.vcf -> sample
    """
    name = vcf_path.name
    if name.endswith(".vcf.gz"):
        name = name[:-7]
    elif name.endswith(".vcf"):
        name = name[:-4]
    return name.split(".", 1)[0]


def open_text_maybe_gzip(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


@dataclass(frozen=True)
class CnSegment:
    chrom: str
    start: int
    end: int
    coverage: float
    copy_number_state: float
    confidence: float


def _to_float_or_nan(x: str) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def load_cn_segments_by_chrom(cn_bed_path: Path) -> Tuple[Dict[str, List[CnSegment]], int, int]:
    """
    Load CN segments from a tab-delimited BED-like file.
    Expected columns:
      chr, start, end, coverage, copynumber_state, confidence, ...
    Header/comment lines starting with '#' are skipped.
    """
    by_chrom: Dict[str, List[CnSegment]] = {}
    total_rows = 0
    kept_rows = 0

    with open(cn_bed_path, "rt", encoding="utf-8", errors="replace") as f:
        for line_no, line in enumerate(f, start=1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            total_rows += 1
            parts = s.split("\t")
            if len(parts) < 6:
                print(f"[cn_bed] Skipping malformed row {line_no}: expected >=6 columns, got {len(parts)}")
                continue

            chrom_raw = parts[0].strip()
            if not chrom_raw:
                print(f"[cn_bed] Skipping row {line_no}: empty chromosome")
                continue
            chrom = canonical_chrom(chrom_raw)

            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                print(f"[cn_bed] Skipping row {line_no}: non-integer start/end")
                continue
            if end < start:
                print(f"[cn_bed] Skipping row {line_no}: end < start ({start}, {end})")
                continue

            seg = CnSegment(
                chrom=chrom,
                start=start,
                end=end,
                coverage=_to_float_or_nan(parts[3]),
                copy_number_state=_to_float_or_nan(parts[4]),
                confidence=_to_float_or_nan(parts[5]),
            )
            by_chrom.setdefault(chrom, []).append(seg)
            kept_rows += 1

    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda seg: (seg.start, seg.end))

    return by_chrom, total_rows, kept_rows


def overlap_bp_closed(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    """Closed-interval overlap length in base pairs."""
    left = max(a_start, b_start)
    right = min(a_end, b_end)
    if right < left:
        return 0
    return right - left + 1


def best_cn_segment_for_span(
    segments_by_chrom: Dict[str, List[CnSegment]],
    chrom: str,
    span_start: int,
    span_end: int,
) -> Tuple[Optional[CnSegment], int]:
    """
    Select overlapping segment using deterministic tie-breakers:
      1) maximum overlap bp
      2) higher confidence
      3) lower segment start
      4) lower segment end
    """
    cc = canonical_chrom(chrom)
    segs = segments_by_chrom.get(cc, [])
    if not segs:
        return None, 0

    best_seg: Optional[CnSegment] = None
    best_key: Optional[Tuple[int, float, int, int]] = None
    best_ov = 0
    for seg in segs:
        if seg.start > span_end:
            break
        if seg.end < span_start:
            continue
        ov = overlap_bp_closed(span_start, span_end, seg.start, seg.end)
        if ov <= 0:
            continue
        conf = seg.confidence
        conf_key = conf if not math.isnan(conf) else float("-inf")
        key = (ov, conf_key, -seg.start, -seg.end)
        if best_key is None or key > best_key:
            best_key = key
            best_seg = seg
            best_ov = ov

    return best_seg, best_ov


def normalize_haplotype_tag(hp: Optional[str]) -> str:
    if hp is None:
        return "UNKNOWN"
    s = str(hp).strip().upper()
    if s in {"1", "HP1", "H1", "HAP1", "HAPLOTYPE1"}:
        return "HP1"
    if s in {"2", "HP2", "H2", "HAP2", "HAPLOTYPE2"}:
        return "HP2"
    if s in {"MIXED"}:
        return "MIXED"
    return "UNKNOWN"


def _is_finite_number(x: Optional[float]) -> bool:
    return x is not None and not math.isnan(x)


def _sum_copy_number_states(cn1: Optional[float], cn2: Optional[float]) -> Optional[float]:
    if _is_finite_number(cn1) and _is_finite_number(cn2):
        return float(cn1) + float(cn2)
    if _is_finite_number(cn1):
        return float(cn1)
    if _is_finite_number(cn2):
        return float(cn2)
    return None


def build_dual_haplotype_cn_annotation(
    chrom: str,
    span_start: int,
    span_end: int,
    haplotype_hint: Optional[str],
    cn_segments_hp1: Dict[str, List[CnSegment]],
    cn_segments_hp2: Dict[str, List[CnSegment]],
) -> Optional[Dict[str, Optional[Union[str, float, int]]]]:
    """
    Annotate a genomic span with CN information from BOTH haplotypes.

    The two BEDs are treated independently (they do not need matching intervals):
    we pick the best-overlap segment from HP1 and HP2 separately, then combine.
    """
    seg1, ov1 = best_cn_segment_for_span(cn_segments_hp1, chrom, span_start, span_end)
    seg2, ov2 = best_cn_segment_for_span(cn_segments_hp2, chrom, span_start, span_end)

    if seg1 is None and seg2 is None:
        return None

    total_cn = _sum_copy_number_states(
        seg1.copy_number_state if seg1 is not None else None,
        seg2.copy_number_state if seg2 is not None else None,
    )

    return {
        "cn_hp1": seg1.copy_number_state if seg1 is not None else None,
        "cn_hp2": seg2.copy_number_state if seg2 is not None else None,
        "cn_total": total_cn,
    }


def write_merged_cn_segments_tsv(
    out_path: Path,
    cn_segments_hp1: Dict[str, List[CnSegment]],
    cn_segments_hp2: Dict[str, List[CnSegment]],
) -> int:
    """
    Write a merged CN table by splitting at all breakpoints from both haplotype BEDs.
    This handles imperfectly aligned intervals by keeping HP1 and HP2 assignments
    separately on each merged interval.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "chrom",
        "start",
        "end",
        "cn_hp1",
        "cn_hp2",
        "cn_total",
    ]
    rows_written = 0
    chroms = sorted(set(cn_segments_hp1.keys()) | set(cn_segments_hp2.keys()))
    with open(out_path, "wt", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for chrom in chroms:
            breakpoints: Set[int] = set()
            for seg in cn_segments_hp1.get(chrom, []):
                breakpoints.add(seg.start)
                breakpoints.add(seg.end + 1)
            for seg in cn_segments_hp2.get(chrom, []):
                breakpoints.add(seg.start)
                breakpoints.add(seg.end + 1)
            points = sorted(breakpoints)
            if len(points) < 2:
                continue

            for i in range(len(points) - 1):
                start = points[i]
                end = points[i + 1] - 1
                if end < start:
                    continue
                ann = build_dual_haplotype_cn_annotation(
                    chrom=chrom,
                    span_start=start,
                    span_end=end,
                    haplotype_hint=None,
                    cn_segments_hp1=cn_segments_hp1,
                    cn_segments_hp2=cn_segments_hp2,
                )
                if ann is None:
                    continue
                writer.writerow(
                    {
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "cn_hp1": ann.get("cn_hp1"),
                        "cn_hp2": ann.get("cn_hp2"),
                        "cn_total": ann.get("cn_total"),
                    }
                )
                rows_written += 1
    return rows_written


def iter_vcf_records(vcf_path: Path):
    """
    Yield (chrom, pos) from a VCF/VCF.GZ. Only uses CHROM and POS columns.
    Skips header lines.
    """
    with open_text_maybe_gzip(vcf_path) as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            yield chrom, pos


def find_graph_files(main_out_dir: Path) -> List[Path]:
    """Recursively find only uncondensed chunk_*.json graphs; skip *_condensed.json."""
    out: List[Path] = []
    for p in main_out_dir.rglob("chunk_*.json"):
        if p.name.endswith("_condensed.json"):
            continue
        out.append(p)
    return sorted(out)


def find_condensed_graph_files(main_out_dir: Path) -> List[Path]:
    """Recursively find condensed graphs: chunk_*_condensed.json"""
    out: List[Path] = []
    for p in main_out_dir.rglob("chunk_*_condensed.json"):
        out.append(p)
    return sorted(out)


def find_vcf_files(vcf_dir: Path, recursive: bool) -> List[Path]:
    pat = "**/*" if recursive else "*"
    files: List[Path] = []
    for p in vcf_dir.glob(pat):
        if not p.is_file():
            continue
        if p.name.endswith(".vcf") or p.name.endswith(".vcf.gz"):
            files.append(p)
    return sorted(files)


def load_needed_sites(graph_files: List[Path]) -> Tuple[Dict[str, Set[int]], int]:
    """
    Return needed[canon_chrom] = set(pos) for all nodes across all graphs.
    Also returns total node count observed (nodes with chrom+position).
    """
    needed: Dict[str, Set[int]] = {}
    total_nodes = 0

    for gf in graph_files:
        with open(gf, "rt", encoding="utf-8") as f:
            g = json.load(f)

        elements = g.get("elements", {})
        nodes = elements.get("nodes", [])
        for n in nodes:
            data = n.get("data", {})
            pos = data.get("position")
            chrom = data.get("chrom")
            if pos is None or chrom is None:
                continue
            try:
                pos_i = int(pos)
            except (TypeError, ValueError):
                continue

            cc = canonical_chrom(chrom)
            needed.setdefault(cc, set()).add(pos_i)
            total_nodes += 1

    return needed, total_nodes


def parse_pos_deltas(s: str) -> List[int]:
    """Parse comma-separated list like '0,-1,1' into [0, -1, 1]"""
    deltas: List[int] = []
    for tok in s.split(","):
        tok = tok.strip()
        if not tok:
            continue
        deltas.append(int(tok))
    if not deltas:
        raise ValueError("pos_deltas parsed to empty list")
    return deltas


def build_site_to_vcfs(
    needed: Dict[str, Set[int]],
    vcf_files: List[Path],
    pos_deltas: List[int],
) -> Tuple[Dict[Tuple[str, int], Set[str]], Dict[int, int]]:
    """
    Stream through VCFs and collect mapping (canon_chrom, graph_pos) -> {vcf_label,...}
    ONLY for sites present in `needed`, allowing graph_pos = vcf_pos + delta.

    Returns:
      mapping: (canon_chrom, graph_pos) -> set(vcf labels)
      delta_hits: delta -> number of matched records
    """
    mapping: Dict[Tuple[str, int], Set[str]] = {}
    needed_chroms = set(needed.keys())
    delta_hits: Dict[int, int] = {d: 0 for d in pos_deltas}

    for vf in vcf_files:
        label = vcf_label(vf)
        hits_this_vcf = 0

        for chrom, vcf_pos in iter_vcf_records(vf):
            cc = canonical_chrom(chrom)
            if cc not in needed_chroms:
                continue

            for d in pos_deltas:
                graph_pos = vcf_pos + d
                if graph_pos in needed[cc]:
                    key = (cc, graph_pos)
                    mapping.setdefault(key, set()).add(label)
                    delta_hits[d] += 1
                    hits_this_vcf += 1

        print(f"[vcf] {vf}: matched {hits_this_vcf} record(s) to graph sites using deltas={pos_deltas}")

    return mapping, delta_hits


def _as_vcf_list(x: Union[str, List[str], None]) -> List[str]:
    if x is None:
        return []
    if isinstance(x, list):
        return [str(v) for v in x]
    return [str(x)]


def annotate_uncondensed_graph_in_memory(
    graph_obj: dict,
    mapping: Dict[Tuple[str, int], Set[str]],
) -> Tuple[int, int]:
    """
    Annotate uncondensed nodes with data['source_vcf'].
    Returns (nodes_seen, nodes_annotated).
    """
    elements = graph_obj.get("elements", {})
    nodes = elements.get("nodes", [])

    seen = 0
    annotated = 0

    for n in nodes:
        data = n.get("data", {})
        pos = data.get("position")
        chrom = data.get("chrom")
        if pos is None or chrom is None:
            continue
        try:
            pos_i = int(pos)
        except (TypeError, ValueError):
            continue

        cc = canonical_chrom(chrom)
        key = (cc, pos_i)

        seen += 1
        if key not in mapping:
            continue

        vcfs = sorted(mapping[key])
        data["source_vcf"] = vcfs[0] if len(vcfs) == 1 else vcfs
        annotated += 1

    return seen, annotated


def annotate_uncondensed_graph_copy_number_in_memory(
    graph_obj: dict,
    cn_segments_hp1: Dict[str, List[CnSegment]],
    cn_segments_hp2: Dict[str, List[CnSegment]],
) -> Tuple[int, int]:
    """
    Annotate uncondensed graph nodes with CN values using node position.
    For each node, pick the segment that overlaps [position, position]
    using the same deterministic tie-breakers as condensed annotation.
    """
    elements = graph_obj.get("elements", {})
    nodes = elements.get("nodes", [])

    seen = 0
    annotated = 0

    for n in nodes:
        data = n.get("data", {})
        pos = data.get("position")
        chrom = data.get("chrom")
        if pos is None or chrom is None:
            continue
        try:
            pos_i = int(pos)
        except (TypeError, ValueError):
            continue
        seen += 1

        ann = build_dual_haplotype_cn_annotation(
            chrom=str(chrom),
            span_start=pos_i,
            span_end=pos_i,
            haplotype_hint=data.get("haplotype"),
            cn_segments_hp1=cn_segments_hp1,
            cn_segments_hp2=cn_segments_hp2,
        )
        if ann is None:
            continue

        # Graph JSON nodes only carry per-haplotype CN; cn_total is omitted.
        for k in ("cn_hp1", "cn_hp2"):
            data[k] = ann.get(k)
        annotated += 1

    return seen, annotated


def annotate_condensed_graph_in_memory(
    graph_obj: dict,
    mapping: Dict[Tuple[str, int], Set[str]],
) -> Tuple[int, int]:
    """
    For each condensed node, set data['source_vcf'] to the union of its members' source VCFs.
    Uses node.data['chrom'] + each member position to look up mapping.

    IMPORTANT: For condensed nodes, source_vcf is ALWAYS a string:
      - if one source: "N1"
      - if multiple:  "N1, N3"
    """
    elements = graph_obj.get("elements", {})
    nodes = elements.get("nodes", [])

    seen = 0
    annotated = 0

    for n in nodes:
        data = n.get("data", {})
        members = data.get("members")
        chrom = data.get("chrom")
        if not members or chrom is None:
            continue

        cc = canonical_chrom(chrom)
        seen += 1

        vcfs_union: Set[str] = set()
        for m in members:
            try:
                pos_i = int(m)
            except (TypeError, ValueError):
                continue
            key = (cc, pos_i)
            if key in mapping:
                vcfs_union.update(mapping[key])

        if not vcfs_union:
            continue

        vcfs = sorted(vcfs_union)
        data["source_vcf"] = ", ".join(vcfs)  # ALWAYS a string for condensed nodes
        annotated += 1

    return seen, annotated


def annotate_condensed_graph_copy_number_in_memory(
    graph_obj: dict,
    cn_segments_hp1: Dict[str, List[CnSegment]],
    cn_segments_hp2: Dict[str, List[CnSegment]],
) -> Tuple[int, int]:
    """
    Annotate condensed graph cluster nodes with CN values using maximal segment overlap.
    Cluster span is computed from member positions as [min(members), max(members)].
    """
    elements = graph_obj.get("elements", {})
    nodes = elements.get("nodes", [])

    seen = 0
    annotated = 0

    for n in nodes:
        data = n.get("data", {})
        members = data.get("members")
        chrom = data.get("chrom")
        if not members or chrom is None:
            continue

        member_positions: List[int] = []
        for m in members:
            try:
                member_positions.append(int(m))
            except (TypeError, ValueError):
                continue
        if not member_positions:
            continue

        span_start = min(member_positions)
        span_end = max(member_positions)
        seen += 1

        ann = build_dual_haplotype_cn_annotation(
            chrom=str(chrom),
            span_start=span_start,
            span_end=span_end,
            haplotype_hint=data.get("haplotype"),
            cn_segments_hp1=cn_segments_hp1,
            cn_segments_hp2=cn_segments_hp2,
        )
        if ann is None:
            continue

        # Graph JSON nodes only carry per-haplotype CN; cn_total is omitted.
        for k in ("cn_hp1", "cn_hp2"):
            data[k] = ann.get(k)
        annotated += 1

    return seen, annotated


def load_chunk_node_haplotypes(graphs_dir: Path) -> Dict[str, Dict[int, str]]:
    """
    Build mapping:
      chunk_base -> {node_position -> haplotype}
    from uncondensed graph JSON files in graphs_dir.
    """
    out: Dict[str, Dict[int, str]] = {}
    if not graphs_dir.is_dir():
        return out

    for jp in sorted(graphs_dir.glob("chunk_*.json")):
        if jp.name.endswith("_condensed.json"):
            continue
        chunk_base = jp.stem
        try:
            with open(jp, "rt", encoding="utf-8") as f:
                obj = json.load(f)
        except Exception:
            continue

        nodes = ((obj.get("elements") or {}).get("nodes") or [])
        pos_to_hp: Dict[int, str] = {}
        for n in nodes:
            data = n.get("data", {}) if isinstance(n, dict) else {}
            pos = data.get("position")
            if pos is None:
                continue
            try:
                pos_i = int(pos)
            except (TypeError, ValueError):
                continue
            pos_to_hp[pos_i] = normalize_haplotype_tag(data.get("haplotype"))
        out[chunk_base] = pos_to_hp
    return out


def majority_haplotype_from_nodes(nodes_csv: str, pos_to_hp: Dict[int, str]) -> str:
    hp1 = 0
    hp2 = 0
    if not nodes_csv:
        return "UNKNOWN"
    for tok in str(nodes_csv).split(","):
        t = tok.strip()
        if not t:
            continue
        try:
            pos = int(t)
        except ValueError:
            continue
        hp = normalize_haplotype_tag(pos_to_hp.get(pos))
        if hp == "HP1":
            hp1 += 1
        elif hp == "HP2":
            hp2 += 1
    if hp1 == 0 and hp2 == 0:
        return "UNKNOWN"
    if hp1 == hp2:
        return "MIXED"
    return "HP1" if hp1 > hp2 else "HP2"


def annotate_component_statistics_file(
    in_path: Path,
    out_path: Path,
    cn_segments_hp1: Dict[str, List[CnSegment]],
    cn_segments_hp2: Dict[str, List[CnSegment]],
    chunk_to_node_haplotype: Dict[str, Dict[int, str]],
) -> Tuple[int, int]:
    """
    Annotate component_statistics.tsv rows using chromosome + [min_pos, max_pos].
    Writes output to out_path (can be equal to in_path for in-place updates).
    """
    with open(in_path, "rt", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames:
            return 0, 0
        rows = list(reader)
        fieldnames = list(reader.fieldnames)

    legacy_cn_fields = [
        "cn_copy_number_state",
        "cn_segment_start",
        "cn_segment_end",
        "cn_overlap_bp",
        "cn_segment_coverage",
        "cn_segment_confidence",
        "cn_haplotype_source",
    ]
    for fld in legacy_cn_fields:
        if fld in fieldnames:
            fieldnames.remove(fld)

    cn_fields = ["cn_hp1", "cn_hp2", "cn_total"]
    for fld in cn_fields:
        if fld not in fieldnames:
            fieldnames.append(fld)

    seen = 0
    annotated = 0
    for row in rows:
        for fld in legacy_cn_fields:
            row.pop(fld, None)
        for fld in cn_fields:
            row.setdefault(fld, "NA")

        chrom = row.get("chromosome")
        min_pos = row.get("min_pos")
        max_pos = row.get("max_pos")
        if chrom is None or min_pos is None or max_pos is None:
            continue

        try:
            span_start = int(min_pos)
            span_end = int(max_pos)
        except (TypeError, ValueError):
            continue
        seen += 1

        chunk_base = str(row.get("chunk_base", ""))
        nodes_csv = str(row.get("nodes", ""))
        pos_to_hp = chunk_to_node_haplotype.get(chunk_base, {})
        hp_label = majority_haplotype_from_nodes(nodes_csv, pos_to_hp)
        ann = build_dual_haplotype_cn_annotation(
            chrom=str(chrom),
            span_start=span_start,
            span_end=span_end,
            haplotype_hint=hp_label,
            cn_segments_hp1=cn_segments_hp1,
            cn_segments_hp2=cn_segments_hp2,
        )
        if ann is None:
            continue

        for fld, val in ann.items():
            row[fld] = "NA" if val is None else str(val)
        annotated += 1

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "wt", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    return seen, annotated


def update_timing_stats_from_uncondensed_graph(
    graph_obj: dict,
    timing_stats: Dict[str, Dict[str, int]],
) -> int:
    """
    Update timing_stats in-place from a single (already annotated) UNCONDENSED graph JSON object.

    Stats definition (uncondensed):
      For each timing edge source->target:
        - timing_out increments for source node's source_vcf(s)
        - timing_in increments for target node's source_vcf(s)

    Returns number of timing edges processed.
    """
    elements = graph_obj.get("elements", {})
    nodes = elements.get("nodes", [])
    edges = elements.get("edges", [])

    node_to_vcfs: Dict[int, List[str]] = {}
    for n in nodes:
        data = n.get("data", {})
        pos = data.get("position")
        if pos is None:
            continue
        try:
            pos_i = int(pos)
        except (TypeError, ValueError):
            continue
        vcfs = _as_vcf_list(data.get("source_vcf"))
        if vcfs:
            node_to_vcfs[pos_i] = vcfs

    timing_edges = 0

    for e in edges:
        ed = e.get("data", {})
        if ed.get("relation") != "timing":
            continue

        src = ed.get("source")
        tgt = ed.get("target")
        if src is None or tgt is None:
            continue
        try:
            src_i = int(src)
            tgt_i = int(tgt)
        except (TypeError, ValueError):
            continue

        for vcf in node_to_vcfs.get(src_i, []):
            timing_stats.setdefault(vcf, {"timing_out": 0, "timing_in": 0})
            timing_stats[vcf]["timing_out"] += 1

        for vcf in node_to_vcfs.get(tgt_i, []):
            timing_stats.setdefault(vcf, {"timing_out": 0, "timing_in": 0})
            timing_stats[vcf]["timing_in"] += 1

        timing_edges += 1

    return timing_edges


def write_json(graph_obj: dict, out_path: Path, indent: int) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "wt", encoding="utf-8") as f:
        json.dump(graph_obj, f, indent=indent, ensure_ascii=False)
        f.write("\n")


def write_stats_tsv(stats: Dict[str, Dict[str, int]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "wt", encoding="utf-8") as f:
        f.write("source_vcf\ttiming_out\ttiming_in\ttiming_total\n")
        for vcf in sorted(stats.keys()):
            out_ = int(stats[vcf].get("timing_out", 0))
            in_ = int(stats[vcf].get("timing_in", 0))
            f.write(f"{vcf}\t{out_}\t{in_}\t{out_ + in_}\n")


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Annotate graphs under a main output directory with source_vcf (via VCF matching) "
            "and/or cluster/component copy number (via CN BED overlap)."
        )
    )
    ap.add_argument(
        "--main_out",
        required=True,
        type=Path,
        help="Main output directory (contains per-chromosome *_out folders, etc.). Script searches recursively.",
    )
    ap.add_argument(
        "--vcfs",
        required=False,
        type=Path,
        help="Directory containing .vcf / .vcf.gz files",
    )
    ap.add_argument(
        "--cn-bed-hp1",
        type=Path,
        default=None,
        help="Copy-number BED segments file for haplotype 1.",
    )
    ap.add_argument(
        "--cn-bed-hp2",
        type=Path,
        default=None,
        help="Copy-number BED segments file for haplotype 2.",
    )
    ap.add_argument(
        "--vcf_recursive",
        action="store_true",
        help="Recursively search for VCFs under --vcfs",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output root directory (default: overwrite in place). Preserves relative paths under --main_out.",
    )
    ap.add_argument(
        "--indent",
        type=int,
        default=2,
        help="JSON indentation level for output graphs (default: 2).",
    )
    ap.add_argument(
        "--pos_deltas",
        type=str,
        default="0,-1,1",
        help=(
            "Comma-separated list of position deltas to try when matching VCF POS to graph node position. "
            "Graph_pos = VCF_POS + delta. Default tries 0,-1,1."
        ),
    )
    ap.add_argument(
        "--stats_name",
        type=str,
        default="timing_edge_stats_by_source_vcf.tsv",
        help="Filename for the stats TSV (written under --out if provided else under --main_out).",
    )
    args = ap.parse_args()

    main_out: Path = args.main_out
    vcf_dir: Optional[Path] = args.vcfs
    out_root: Optional[Path] = args.out
    do_vcf = vcf_dir is not None
    if (args.cn_bed_hp1 is None) != (args.cn_bed_hp2 is None):
        raise SystemExit("Provide both --cn-bed-hp1 and --cn-bed-hp2 together.")
    do_cn = args.cn_bed_hp1 is not None and args.cn_bed_hp2 is not None
    if not do_vcf and not do_cn:
        raise SystemExit("Provide at least one annotation source: --vcfs and/or both --cn-bed-hp1/--cn-bed-hp2.")

    uncondensed_files = find_graph_files(main_out)
    condensed_files = find_condensed_graph_files(main_out)
    component_stats_files = sorted(main_out.glob("*_out/component_statistics.txt"))

    print(f"[graphs] Found {len(uncondensed_files)} uncondensed graph(s) under {main_out}")
    print(f"[condensed] Found {len(condensed_files)} condensed graph(s) under {main_out}")
    print(f"[component_stats] Found {len(component_stats_files)} component_statistics.txt file(s)")

    mapping: Dict[Tuple[str, int], Set[str]] = {}
    if do_vcf:
        pos_deltas = parse_pos_deltas(args.pos_deltas)
        if not uncondensed_files:
            raise SystemExit(
                f"No uncondensed chunk_*.json graphs found under {main_out} "
                f"(note: *_condensed.json are intentionally skipped)."
            )
        vcf_files = find_vcf_files(vcf_dir, args.vcf_recursive)
        if not vcf_files:
            raise SystemExit(f"No .vcf or .vcf.gz files found in {vcf_dir} (recursive={args.vcf_recursive})")

        print(f"[vcf] Found {len(vcf_files)} VCF file(s) under {vcf_dir}")

        needed, total_nodes_counted = load_needed_sites(uncondensed_files)
        total_needed_sites = sum(len(s) for s in needed.values())
        print(f"[sites] Collected {total_needed_sites} unique (chrom,pos) site(s) from uncondensed graphs")
        print(f"[nodes] Counted {total_nodes_counted} node(s) with chrom+position across uncondensed graphs")

        if total_needed_sites == 0:
            raise SystemExit(
                "Found 0 graph sites with (chrom, position) in uncondensed graphs. "
                "This usually means node.data['chrom'] is missing.\n"
                "Re-generate graphs after adding 'chrom' to nodes, or confirm you're pointing to the right directory."
            )

        mapping, delta_hits = build_site_to_vcfs(needed, vcf_files, pos_deltas)
        print(f"[map] Built mapping for {len(mapping)} unique graph site(s) present in both graphs and VCFs")
        print("[delta_hits] Matches by delta (graph_pos = vcf_pos + delta):")
        for d in pos_deltas:
            print(f"  delta {d:+d}: {delta_hits.get(d, 0)} matched record(s)")

    cn_segments_hp1: Dict[str, List[CnSegment]] = {}
    cn_segments_hp2: Dict[str, List[CnSegment]] = {}
    if do_cn:
        cn_bed_hp1_path = args.cn_bed_hp1.resolve()
        cn_bed_hp2_path = args.cn_bed_hp2.resolve()
        cn_segments_hp1, cn_rows_total_hp1, cn_rows_kept_hp1 = load_cn_segments_by_chrom(
            cn_bed_hp1_path
        )
        cn_segments_hp2, cn_rows_total_hp2, cn_rows_kept_hp2 = load_cn_segments_by_chrom(
            cn_bed_hp2_path
        )
        print(
            f"[cn_bed_hp1] Loaded {cn_rows_kept_hp1}/{cn_rows_total_hp1} segment row(s) "
            f"across {len(cn_segments_hp1)} chromosome(s) from {cn_bed_hp1_path}"
        )
        print(
            f"[cn_bed_hp2] Loaded {cn_rows_kept_hp2}/{cn_rows_total_hp2} segment row(s) "
            f"across {len(cn_segments_hp2)} chromosome(s) from {cn_bed_hp2_path}"
        )
        merged_cn_out_root = out_root if out_root is not None else main_out
        merged_cn_path = merged_cn_out_root / "merged_haplotype_cn_segments.tsv"
        merged_rows = write_merged_cn_segments_tsv(
            out_path=merged_cn_path,
            cn_segments_hp1=cn_segments_hp1,
            cn_segments_hp2=cn_segments_hp2,
        )
        print(
            f"[cn_bed_merged] Wrote {merged_rows} merged interval row(s) with HP1/HP2 CN values to "
            f"{merged_cn_path}"
        )
        print(
            "[cn_bed] Using dual-haplotype CN assignment (HP1 + HP2) with robust overlap matching; "
            "graph CN fields are: cn_hp1, cn_hp2."
        )

    timing_stats: Dict[str, Dict[str, int]] = {}
    total_seen = 0
    total_annotated = 0
    total_ucn_seen = 0
    total_ucn_annotated = 0
    total_timing_edges = 0

    # 1) UNCONDENSED: source_vcf annotation and/or CN annotation + write
    if do_vcf or do_cn:
        for gf in uncondensed_files:
            with open(gf, "rt", encoding="utf-8") as f:
                g = json.load(f)

            seen = 0
            annotated = 0
            cn_seen = 0
            cn_annot = 0
            timing_edges_here = 0

            if do_vcf:
                seen, annotated = annotate_uncondensed_graph_in_memory(g, mapping)
                total_seen += seen
                total_annotated += annotated

                timing_edges_here = update_timing_stats_from_uncondensed_graph(g, timing_stats)
                total_timing_edges += timing_edges_here

            if do_cn:
                cn_seen, cn_annot = annotate_uncondensed_graph_copy_number_in_memory(
                    g,
                    cn_segments_hp1=cn_segments_hp1,
                    cn_segments_hp2=cn_segments_hp2,
                )
                total_ucn_seen += cn_seen
                total_ucn_annotated += cn_annot

            out_path = gf if out_root is None else (out_root / gf.relative_to(main_out))
            write_json(g, out_path, indent=args.indent)
            print(
                f"[write] {gf}: source_vcf={annotated}/{seen} "
                f"cn={cn_annot}/{cn_seen} node(s); timing_edges={timing_edges_here} -> {out_path}"
            )

        if do_vcf:
            # Timing stats TSV
            stats_root = out_root if out_root is not None else main_out
            stats_path = stats_root / args.stats_name
            write_stats_tsv(timing_stats, stats_path)

            print(f"[stats] Wrote: {stats_path}")
            print(
                f"[done] Uncondensed: source_vcf annotated {total_annotated}/{total_seen} "
                f"node(s) across {len(uncondensed_files)} graph(s)"
            )
            print(f"[done] Uncondensed: processed {total_timing_edges} timing edge(s) total")
        if do_cn:
            print(
                f"[done] Uncondensed: CN annotated {total_ucn_annotated}/{total_ucn_seen} "
                f"node(s) across {len(uncondensed_files)} graph(s)"
            )

    # 2) CONDENSED: annotate source_vcf and/or CN on cluster nodes + write
    total_c_seen = 0
    total_c_annotated = 0
    total_cn_seen = 0
    total_cn_annotated = 0

    for cf in condensed_files:
        with open(cf, "rt", encoding="utf-8") as f:
            cg = json.load(f)

        c_seen = 0
        c_annot = 0
        cn_seen = 0
        cn_annot = 0

        if do_vcf:
            c_seen, c_annot = annotate_condensed_graph_in_memory(cg, mapping)
        if do_cn:
            cn_seen, cn_annot = annotate_condensed_graph_copy_number_in_memory(
                cg,
                cn_segments_hp1=cn_segments_hp1,
                cn_segments_hp2=cn_segments_hp2,
            )

        total_c_seen += c_seen
        total_c_annotated += c_annot
        total_cn_seen += cn_seen
        total_cn_annotated += cn_annot

        out_path = cf if out_root is None else (out_root / cf.relative_to(main_out))
        write_json(cg, out_path, indent=args.indent)
        print(
            f"[write_condensed] {cf}: source_vcf={c_annot}/{c_seen} "
            f"cn={cn_annot}/{cn_seen} cluster node(s) -> {out_path}"
        )

    if do_vcf:
        print(
            f"[condensed_done] source_vcf annotated {total_c_annotated}/{total_c_seen} "
            f"condensed cluster node(s) across {len(condensed_files)} graph(s)"
        )
    if do_cn:
        print(
            f"[condensed_done] CN annotated {total_cn_annotated}/{total_cn_seen} "
            f"condensed cluster node(s) across {len(condensed_files)} graph(s)"
        )

    # 3) COMPONENT STATS: CN annotation by component span overlap
    if do_cn:
        comp_seen = 0
        comp_annot = 0
        for comp_path in component_stats_files:
            chrom_out_dir = comp_path.parent
            graphs_dir = chrom_out_dir / "graphs"
            chunk_to_node_haplotype = load_chunk_node_haplotypes(graphs_dir)
            out_path = comp_path if out_root is None else (out_root / comp_path.relative_to(main_out))
            seen_here, annot_here = annotate_component_statistics_file(
                in_path=comp_path,
                out_path=out_path,
                cn_segments_hp1=cn_segments_hp1,
                cn_segments_hp2=cn_segments_hp2,
                chunk_to_node_haplotype=chunk_to_node_haplotype,
            )
            comp_seen += seen_here
            comp_annot += annot_here
            print(f"[write_component_stats] {comp_path}: annotated {annot_here}/{seen_here} row(s) -> {out_path}")
        print(
            f"[component_stats_done] CN annotated {comp_annot}/{comp_seen} rows "
            f"across {len(component_stats_files)} file(s)"
        )


if __name__ == "__main__":
    main()
