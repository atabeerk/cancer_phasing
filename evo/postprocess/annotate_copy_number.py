"""Annotate final graphs and component spans with haplotype copy numbers."""

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union


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


@dataclass(frozen=True)
class CnSegment:
    start: int
    end: int
    copy_number_state: float
    confidence: float


def _to_float_or_nan(x: str) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def load_dual_haplotype_cn_segments_by_chrom(
    cn_bed_path: Path,
) -> Tuple[
    Dict[str, List[CnSegment]],
    Dict[str, List[CnSegment]],
]:
    """
    Load both haplotype CN tracks from a Wakhan integer_profile.bed file.
    Expected columns:
      chr, start, end,
      hp1_coverage, hp1_copynumber_state, hp1_confidence,
      hp2_coverage, hp2_copynumber_state, hp2_confidence, ...
    Header/comment lines starting with '#' are skipped.
    """
    hp1_by_chrom: Dict[str, List[CnSegment]] = {}
    hp2_by_chrom: Dict[str, List[CnSegment]] = {}

    with open(cn_bed_path, "rt", encoding="utf-8", errors="replace") as f:
        for line_no, line in enumerate(f, start=1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split("\t")
            if len(parts) < 9:
                print(f"[cn_bed] Skipping malformed row {line_no}: expected >=9 columns, got {len(parts)}")
                continue

            chrom_raw = parts[0].strip()
            if not chrom_raw:
                print(f"[cn_bed] Skipping row {line_no}: empty chromosome")
                continue
            chrom = canonical_chrom(chrom_raw)

            try:
                bed_start = int(parts[1])
                bed_end = int(parts[2])
            except ValueError:
                print(f"[cn_bed] Skipping row {line_no}: non-integer start/end")
                continue
            if bed_end <= bed_start:
                print(
                    f"[cn_bed] Skipping row {line_no}: end <= start "
                    f"({bed_start}, {bed_end})"
                )
                continue

            # BED is zero-based and end-exclusive; graph/VCF positions are
            # one-based. Convert [bed_start, bed_end) to [start, end].
            start = bed_start + 1
            end = bed_end

            hp1_seg = CnSegment(
                start=start,
                end=end,
                copy_number_state=_to_float_or_nan(parts[4]),
                confidence=_to_float_or_nan(parts[5]),
            )
            hp2_seg = CnSegment(
                start=start,
                end=end,
                copy_number_state=_to_float_or_nan(parts[7]),
                confidence=_to_float_or_nan(parts[8]),
            )
            hp1_by_chrom.setdefault(chrom, []).append(hp1_seg)
            hp2_by_chrom.setdefault(chrom, []).append(hp2_seg)

    for by_chrom in (hp1_by_chrom, hp2_by_chrom):
        for chrom in by_chrom:
            by_chrom[chrom].sort(key=lambda seg: (seg.start, seg.end))

    return hp1_by_chrom, hp2_by_chrom


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
    cn_segments_hp1: Dict[str, List[CnSegment]],
    cn_segments_hp2: Dict[str, List[CnSegment]],
) -> Optional[Dict[str, Optional[Union[str, float, int]]]]:
    """
    Annotate a genomic span with CN information from BOTH haplotypes.

    We pick the best-overlap segment from the HP1 and HP2 tracks, then combine.
    """
    seg1, _ = best_cn_segment_for_span(cn_segments_hp1, chrom, span_start, span_end)
    seg2, _ = best_cn_segment_for_span(cn_segments_hp2, chrom, span_start, span_end)

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
) -> None:
    """
    Write a merged CN table by splitting at all breakpoints from both haplotype tracks.
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


def annotate_graph_copy_number(graph, hp1_segments, hp2_segments, condensed=False):
    """Add CN labels using mutation positions or condensed cluster spans."""
    annotated = 0
    for node in graph["elements"]["nodes"]:
        data = node["data"]
        positions = data["members"] if condensed else [data["position"]]
        positions = [int(position) for position in positions]
        annotation = build_dual_haplotype_cn_annotation(
            data["chrom"], min(positions), max(positions), hp1_segments, hp2_segments
        )
        if annotation is not None:
            data.update({key: annotation[key] for key in ("cn_hp1", "cn_hp2")})
            annotated += 1
    return annotated


def annotate_component_statistics(path, hp1_segments, hp2_segments):
    """Annotate component spans without recalculating any graph statistics."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        rows = list(reader)
    if fields is None:
        raise ValueError(f"Missing component-statistics header: {path}")
    annotated = 0
    for row in rows:
        annotation = build_dual_haplotype_cn_annotation(
            row["chromosome"], int(row["min_pos"]), int(row["max_pos"]), hp1_segments, hp2_segments
        )
        if annotation is not None:
            row.update({key: "NA" if value is None else str(value) for key, value in annotation.items()})
            annotated += 1
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return annotated


def annotate_copy_number(outdir, chromosomes, cn_bed):
    """Annotate this run's graphs; preserve read-derived VAFs and all edges."""
    outdir = Path(outdir)
    hp1, hp2 = load_dual_haplotype_cn_segments_by_chrom(Path(cn_bed))
    selected = {canonical_chrom(chrom) for chrom, _length in chromosomes}
    hp1 = {chrom: segments for chrom, segments in hp1.items() if chrom in selected}
    hp2 = {chrom: segments for chrom, segments in hp2.items() if chrom in selected}
    write_merged_cn_segments_tsv(outdir / "merged_haplotype_cn_segments.tsv", hp1, hp2)
    graph_paths = []
    annotated_components = 0
    for chrom, _length in chromosomes:
        chrom_out = outdir / f"{chrom}_out"
        graph_paths.extend((chrom_out / "graphs" / f"{chrom}{suffix}.json", suffix == "_condensed")
                           for suffix in ("", "_condensed", "_simplified"))
        stats = chrom_out / "component_statistics.txt"
        if stats.exists():
            annotated_components += annotate_component_statistics(stats, hp1, hp2)
    graph_paths.extend((path, False) for path in sorted((outdir / "topN_graphs").glob("*.json")))
    annotated_nodes = 0
    for path, condensed in graph_paths:
        if not path.exists():
            continue
        with path.open(encoding="utf-8") as handle:
            graph = json.load(handle)
        annotated_nodes += annotate_graph_copy_number(graph, hp1, hp2, condensed)
        with path.open("w", encoding="utf-8") as handle:
            json.dump(graph, handle, indent=2)
    return annotated_nodes, annotated_components
