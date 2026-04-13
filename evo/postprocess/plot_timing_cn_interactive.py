#!/usr/bin/env python3
"""
Interactive timing plot:
  - With CN BEDs:
      Top panel shows CN tracks, timing chains, rank markers, chain>maxCN alerts.
      Bottom panel shows circle-only distribution by CN state.
  - Without CN BEDs:
      Top panel shows timing chains/rank markers only.
      Bottom panel shows circle-only distributions over chain length or rank.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import plotly.graph_objects as go


@dataclass(frozen=True)
class CnSegment:
    chrom: str
    start: int
    end: int
    cn: float


@dataclass(frozen=True)
class TimingNode:
    chrom: str
    start: int
    end: int
    mut_count: int
    chain_index: int
    node_rank: int
    node_id: str
    file_rel: str
    chain_len: int
    component_id: Optional[int]


@dataclass(frozen=True)
class TimingChain:
    chrom: str
    start: int
    end: int
    chain_len: int
    chain_index: int
    file_rel: str
    total_mutations: int
    component_id: Optional[int]
    nodes: Tuple[TimingNode, ...]


@dataclass(frozen=True)
class RankedNode:
    chrom: str
    start: int
    end: int
    member_count: int
    rank: int
    node_id: str
    file_rel: str
    component_id: Optional[int]
    member_positions: Tuple[int, ...]


@dataclass(frozen=True)
class MutationVafPoint:
    vaf: float
    chrom: str
    position: int
    node_id: str
    file_rel: str


def canonical_chrom(chrom: str) -> str:
    c = str(chrom).strip()
    if c.lower().startswith("chr"):
        c = c[3:].strip()
    up = c.upper()
    if up in {"M", "MT", "MITO"}:
        return "MT"
    return up


def chrom_sort_key(chrom: str) -> Tuple[int, int, str]:
    cc = canonical_chrom(chrom)
    if cc.isdigit():
        return (0, int(cc), "")
    if cc == "X":
        return (1, 23, "")
    if cc == "Y":
        return (1, 24, "")
    if cc == "MT":
        return (1, 25, "")
    return (2, 999, cc)


def parse_float_or_nan(s: str) -> float:
    try:
        return float(s)
    except (TypeError, ValueError):
        return float("nan")


def parse_members(raw: object) -> List[int]:
    if raw is None:
        return []
    if isinstance(raw, list):
        out: List[int] = []
        for x in raw:
            try:
                out.append(int(x))
            except (TypeError, ValueError):
                continue
        return out
    text = str(raw).strip()
    if not text:
        return []
    out: List[int] = []
    for tok in text.split(","):
        tok = tok.strip()
        if not tok:
            continue
        try:
            out.append(int(tok))
        except ValueError:
            continue
    return out


def load_cn_segments(cn_bed: Path) -> Dict[str, List[CnSegment]]:
    by_chrom: Dict[str, List[CnSegment]] = {}
    with cn_bed.open("rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split("\t")
            if len(parts) < 5:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            if end < start:
                continue
            chrom = canonical_chrom(parts[0])
            cn = parse_float_or_nan(parts[4])
            by_chrom.setdefault(chrom, []).append(CnSegment(chrom, start, end, cn))
    for c in by_chrom:
        by_chrom[c].sort(key=lambda x: (x.start, x.end))
    return by_chrom


def _sum_cn(cn1: float, cn2: float) -> Optional[float]:
    has1 = not math.isnan(cn1)
    has2 = not math.isnan(cn2)
    if has1 and has2:
        return cn1 + cn2
    if has1:
        return cn1
    if has2:
        return cn2
    return None


def _max_cn(cn1: float, cn2: float) -> Optional[float]:
    has1 = not math.isnan(cn1)
    has2 = not math.isnan(cn2)
    if has1 and has2:
        return max(cn1, cn2)
    if has1:
        return cn1
    if has2:
        return cn2
    return None


def _build_cn_for_chrom(hp1: List[CnSegment], hp2: List[CnSegment], chrom: str, reducer: str) -> List[CnSegment]:
    bps = set()
    for s in hp1:
        bps.add(s.start)
        bps.add(s.end + 1)
    for s in hp2:
        bps.add(s.start)
        bps.add(s.end + 1)
    if len(bps) < 2:
        return []

    cuts = sorted(bps)
    out: List[CnSegment] = []
    i1 = 0
    i2 = 0
    for i in range(len(cuts) - 1):
        start = cuts[i]
        end = cuts[i + 1] - 1
        if end < start:
            continue
        while i1 < len(hp1) and hp1[i1].end < start:
            i1 += 1
        while i2 < len(hp2) and hp2[i2].end < start:
            i2 += 1
        cn1 = float("nan")
        cn2 = float("nan")
        if i1 < len(hp1) and hp1[i1].start <= start <= hp1[i1].end:
            cn1 = hp1[i1].cn
        if i2 < len(hp2) and hp2[i2].start <= start <= hp2[i2].end:
            cn2 = hp2[i2].cn

        val = _sum_cn(cn1, cn2) if reducer == "sum" else _max_cn(cn1, cn2)
        if val is None:
            continue
        out.append(CnSegment(chrom, start, end, val))

    if not out:
        return out
    merged = [out[0]]
    for seg in out[1:]:
        prev = merged[-1]
        if prev.end + 1 == seg.start and prev.cn == seg.cn:
            merged[-1] = CnSegment(prev.chrom, prev.start, seg.end, prev.cn)
        else:
            merged.append(seg)
    return merged


def build_total_cn_segments(cn_hp1: Dict[str, List[CnSegment]], cn_hp2: Dict[str, List[CnSegment]]) -> Dict[str, List[CnSegment]]:
    out: Dict[str, List[CnSegment]] = {}
    for c in sorted(set(cn_hp1) | set(cn_hp2), key=chrom_sort_key):
        segs = _build_cn_for_chrom(cn_hp1.get(c, []), cn_hp2.get(c, []), c, reducer="sum")
        if segs:
            out[c] = segs
    return out


def build_max_cn_segments(cn_hp1: Dict[str, List[CnSegment]], cn_hp2: Dict[str, List[CnSegment]]) -> Dict[str, List[CnSegment]]:
    out: Dict[str, List[CnSegment]] = {}
    for c in sorted(set(cn_hp1) | set(cn_hp2), key=chrom_sort_key):
        segs = _build_cn_for_chrom(cn_hp1.get(c, []), cn_hp2.get(c, []), c, reducer="max")
        if segs:
            out[c] = segs
    return out


def load_condensed_node_meta(outdir: Path, rel_file: str) -> Dict[str, Tuple[Optional[str], Optional[int], Optional[int], Optional[int]]]:
    graph_path = outdir / rel_file
    if not graph_path.exists():
        return {}
    with graph_path.open("rt", encoding="utf-8") as f:
        obj = json.load(f)
    if isinstance(obj, dict) and isinstance(obj.get("elements"), dict):
        nodes = obj["elements"].get("nodes", []) or []
    else:
        nodes = obj.get("nodes", []) or []
    out: Dict[str, Tuple[Optional[str], Optional[int], Optional[int], Optional[int]]] = {}
    for n in nodes:
        data = n.get("data", {}) or {}
        nid = data.get("id")
        if nid is None:
            continue
        members = parse_members(data.get("members"))
        chrom_raw = data.get("chrom")
        chrom = canonical_chrom(chrom_raw) if chrom_raw is not None else None
        if members:
            out[str(nid)] = (chrom, min(members), max(members), len(members))
        else:
            out[str(nid)] = (chrom, None, None, None)
    return out


def load_condensed_component_ids(outdir: Path, rel_file: str) -> Dict[str, int]:
    graph_path = outdir / rel_file
    if not graph_path.exists():
        return {}
    with graph_path.open("rt", encoding="utf-8") as f:
        obj = json.load(f)
    if isinstance(obj, dict) and isinstance(obj.get("elements"), dict):
        nodes = obj["elements"].get("nodes", []) or []
        edges = obj["elements"].get("edges", []) or []
    else:
        nodes = obj.get("nodes", []) or []
        edges = obj.get("edges", []) or []
    node_ids = []
    for n in nodes:
        nid = (n.get("data", {}) or {}).get("id")
        if nid is not None:
            node_ids.append(str(nid))
    adj: Dict[str, List[str]] = {n: [] for n in node_ids}
    for e in edges:
        d = e.get("data", {}) or {}
        s = d.get("source")
        t = d.get("target")
        if s is None or t is None:
            continue
        s = str(s)
        t = str(t)
        if s in adj and t in adj:
            adj[s].append(t)
            adj[t].append(s)

    comp: Dict[str, int] = {}
    cid = 0
    for n in node_ids:
        if n in comp:
            continue
        stack = [n]
        comp[n] = cid
        while stack:
            cur = stack.pop()
            for nb in adj.get(cur, []):
                if nb in comp:
                    continue
                comp[nb] = cid
                stack.append(nb)
        cid += 1
    return comp


def load_timing_chains(outdir: Path, timing_tsv: Path, min_chain_len: int) -> List[TimingChain]:
    if not timing_tsv.exists():
        raise FileNotFoundError(f"Timing TSV not found: {timing_tsv}")
    rows: List[dict] = []
    with timing_tsv.open("rt", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        req = {"file", "chain_index", "node_rank", "node_id", "members"}
        missing = req.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{timing_tsv} missing columns: {', '.join(sorted(missing))}")
        for r in reader:
            rows.append(r)

    chain_sizes: Dict[Tuple[str, int], int] = {}
    for r in rows:
        file_rel = (r.get("file") or "").strip()
        try:
            idx = int((r.get("chain_index") or "").strip())
        except ValueError:
            continue
        k = (file_rel, idx)
        chain_sizes[k] = chain_sizes.get(k, 0) + 1

    per_file_meta: Dict[str, Dict[str, Tuple[Optional[str], Optional[int], Optional[int], Optional[int]]]] = {}
    per_file_comp: Dict[str, Dict[str, int]] = {}
    grouped: Dict[Tuple[str, int], List[TimingNode]] = {}
    for r in rows:
        file_rel = (r.get("file") or "").strip()
        if not file_rel:
            continue
        try:
            chain_index = int((r.get("chain_index") or "").strip())
            node_rank = int((r.get("node_rank") or "").strip())
        except ValueError:
            continue
        chain_len = chain_sizes.get((file_rel, chain_index), 0)
        if chain_len < min_chain_len:
            continue

        node_id = str(r.get("node_id") or "").strip()
        if not node_id:
            continue
        if file_rel not in per_file_meta:
            per_file_meta[file_rel] = load_condensed_node_meta(outdir, file_rel)
        if file_rel not in per_file_comp:
            per_file_comp[file_rel] = load_condensed_component_ids(outdir, file_rel)
        meta = per_file_meta[file_rel].get(node_id)
        comp_id = per_file_comp[file_rel].get(node_id)

        chrom = None
        start = None
        end = None
        mut_count = None
        if meta is not None:
            chrom, start, end, mut_count = meta
        if start is None or end is None or mut_count is None:
            members = parse_members(r.get("members"))
            if members:
                start = min(members)
                end = max(members)
                mut_count = len(members)
        if chrom is None or start is None or end is None or mut_count is None:
            continue
        if end < start:
            continue
        grouped.setdefault((file_rel, chain_index), []).append(
            TimingNode(
                chrom=chrom,
                start=start,
                end=end,
                mut_count=mut_count,
                chain_index=chain_index,
                node_rank=node_rank,
                node_id=node_id,
                file_rel=file_rel,
                chain_len=chain_len,
                component_id=comp_id,
            )
        )

    out: List[TimingChain] = []
    for (file_rel, chain_index), nodes in grouped.items():
        nodes_sorted = sorted(nodes, key=lambda n: (n.node_rank, n.start, n.end, n.node_id))
        chrom = nodes_sorted[0].chrom
        chain_start = min(n.start for n in nodes_sorted)
        chain_end = max(n.end for n in nodes_sorted)
        total_mut = sum(n.mut_count for n in nodes_sorted)
        comps = [n.component_id for n in nodes_sorted if n.component_id is not None]
        comp_id = comps[0] if comps else None
        out.append(
            TimingChain(
                chrom=chrom,
                start=chain_start,
                end=chain_end,
                chain_len=nodes_sorted[0].chain_len,
                chain_index=chain_index,
                file_rel=file_rel,
                total_mutations=total_mut,
                component_id=comp_id,
                nodes=tuple(nodes_sorted),
            )
        )
    return out


def load_node_ranks(outdir: Path, ranks_tsv: Path) -> List[RankedNode]:
    if not ranks_tsv.exists():
        return []
    out: List[RankedNode] = []
    per_file_comp: Dict[str, Dict[str, int]] = {}
    with ranks_tsv.open("rt", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            rank_raw = str(r.get("rank") or "").strip()
            if not rank_raw or rank_raw == "N/A":
                continue
            try:
                rank = int(rank_raw)
                start = int(str(r.get("start") or "0"))
                end = int(str(r.get("end") or "0"))
                member_count = int(str(r.get("member_count") or "1"))
            except ValueError:
                continue
            if end < start:
                continue
            file_rel = str(r.get("file") or "").strip()
            node_id = str(r.get("node_id") or "").strip()
            chrom_raw = str(r.get("chrom") or "").strip()
            if not file_rel or not node_id or not chrom_raw:
                continue
            if file_rel not in per_file_comp:
                per_file_comp[file_rel] = load_condensed_component_ids(outdir, file_rel)
            comp_id = per_file_comp[file_rel].get(node_id)
            members = parse_members(r.get("members"))
            if not members:
                members = [start] if start == end else [start, end]
            out.append(
                RankedNode(
                    chrom=canonical_chrom(chrom_raw),
                    start=start,
                    end=end,
                    member_count=member_count,
                    rank=rank,
                    node_id=node_id,
                    file_rel=file_rel,
                    component_id=comp_id,
                    member_positions=tuple(sorted(members)),
                )
            )
    return out


def load_uncondensed_mutation_vaf_points(outdir: Path) -> List[MutationVafPoint]:
    out: List[MutationVafPoint] = []
    for jp in sorted(outdir.rglob("chunk_*.json")):
        if jp.name.endswith("_condensed.json"):
            continue
        try:
            with jp.open("rt", encoding="utf-8") as f:
                obj = json.load(f)
        except Exception:
            continue
        if isinstance(obj, dict) and isinstance(obj.get("elements"), dict):
            nodes = obj["elements"].get("nodes", []) or []
        else:
            nodes = obj.get("nodes", []) or []
        for n in nodes:
            data = n.get("data", {}) if isinstance(n, dict) else {}
            vaf_raw = data.get("vaf")
            pos_raw = data.get("position")
            nid = data.get("id")
            chrom_raw = data.get("chrom")
            if vaf_raw is None or pos_raw is None or nid is None or chrom_raw is None:
                continue
            try:
                vaf = float(vaf_raw)
                pos = int(pos_raw)
            except (TypeError, ValueError):
                continue
            if not math.isfinite(vaf):
                continue
            if vaf < 0.0 or vaf > 1.0:
                continue
            try:
                file_rel = str(jp.relative_to(outdir))
            except ValueError:
                file_rel = str(jp)
            out.append(
                MutationVafPoint(
                    vaf=vaf,
                    chrom=canonical_chrom(str(chrom_raw)),
                    position=pos,
                    node_id=str(nid),
                    file_rel=file_rel,
                )
            )
    return out


def add_vaf_distribution_trace(fig: go.Figure, mutation_vaf_points: Optional[List[MutationVafPoint]]) -> None:
    if not mutation_vaf_points:
        fig.update_layout(meta={**(fig.layout.meta or {}), "dist_has_vaf": False})
        return
    bin_w = 0.02
    n_bins = 50
    counts = [0] * n_bins
    examples: Dict[int, MutationVafPoint] = {}
    for p in mutation_vaf_points:
        v = max(0.0, min(1.0, float(p.vaf)))
        idx = min(n_bins - 1, int(v / bin_w))
        counts[idx] += 1
        if idx not in examples:
            examples[idx] = p
    xs: List[float] = []
    ys: List[int] = []
    hover: List[str] = []
    text: List[str] = []
    for i in range(n_bins):
        c = counts[i]
        if c <= 0:
            continue
        lo = i * bin_w
        hi = lo + bin_w
        center = lo + (bin_w / 2.0)
        xs.append(center)
        ys.append(c)
        text.append(str(c))
        ex = examples.get(i)
        extra = ""
        if ex is not None:
            extra = (
                f"<br>example: chr={ex.chrom} pos={ex.position:,} node_id={ex.node_id}"
            )
        hover.append(
            f"VAF bin=[{lo:.2f}, {hi:.2f})<br>mutations={c}{extra}"
        )
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=ys,
            mode="markers+text",
            text=text,
            textposition="top center",
            name=f"Mutations by VAF bins ({len(mutation_vaf_points):,})",
            marker=dict(size=[8 + 4 * math.log2(max(c, 0) + 1) for c in ys], color="#17becf", opacity=0.8),
            customdata=hover,
            hovertemplate="%{customdata}<extra></extra>",
            visible=False,
            meta={"trace_role": "dist_vaf"},
        )
    )
    fig.update_layout(meta={**(fig.layout.meta or {}), "dist_has_vaf": True})


def select_non_overlapping_chains(chains: List[TimingChain]) -> List[TimingChain]:
    grouped: Dict[Tuple[str, str], List[TimingChain]] = {}
    for c in chains:
        key = str(c.component_id) if c.component_id is not None else f"na:{c.chrom}"
        grouped.setdefault((c.file_rel, key), []).append(c)
    selected: List[TimingChain] = []
    for _, clist in grouped.items():
        used: set[str] = set()
        ordered = sorted(clist, key=lambda c: (-c.chain_len, -(c.end - c.start + 1), c.chain_index))
        for c in ordered:
            ids = {n.node_id for n in c.nodes}
            if ids & used:
                continue
            selected.append(c)
            used |= ids
    return selected


def select_longest_chain_per_cluster(chains: List[TimingChain]) -> List[TimingChain]:
    grouped: Dict[Tuple[str, str], List[TimingChain]] = {}
    for c in chains:
        key = str(c.component_id) if c.component_id is not None else f"na:{c.chrom}"
        grouped.setdefault((c.file_rel, key), []).append(c)
    out: List[TimingChain] = []
    for _, clist in grouped.items():
        out.append(max(clist, key=lambda c: (c.chain_len, (c.end - c.start + 1), -c.chain_index)))
    return out


def get_chrom_lengths(
    cn_hp1: Dict[str, List[CnSegment]],
    cn_hp2: Dict[str, List[CnSegment]],
    cn_total: Dict[str, List[CnSegment]],
    timing_chains: Iterable[TimingChain],
    ranked_nodes: Optional[List[RankedNode]] = None,
) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    for by_chrom in (cn_hp1, cn_hp2, cn_total):
        for chrom, segs in by_chrom.items():
            if segs:
                lengths[chrom] = max(lengths.get(chrom, 0), max(s.end for s in segs))
    for c in timing_chains:
        lengths[c.chrom] = max(lengths.get(c.chrom, 0), c.end)
    for n in (ranked_nodes or []):
        lengths[n.chrom] = max(lengths.get(n.chrom, 0), n.end)
    return lengths


def chrom_offsets(chrom_lengths: Dict[str, int], pad_bp: int = 1_000_000) -> Dict[str, int]:
    out: Dict[str, int] = {}
    cur = 0
    for chrom in sorted(chrom_lengths, key=chrom_sort_key):
        out[chrom] = cur
        cur += chrom_lengths[chrom] + pad_bp
    return out


def add_cn_trace(fig: go.Figure, name: str, by_chrom: Dict[str, List[CnSegment]], offsets: Dict[str, int], color: str) -> None:
    x: List[Optional[int]] = []
    y: List[Optional[float]] = []
    txt: List[Optional[str]] = []
    for chrom in sorted(by_chrom, key=chrom_sort_key):
        off = offsets.get(chrom, 0)
        for seg in by_chrom[chrom]:
            if math.isnan(seg.cn):
                continue
            x0 = off + seg.start
            x1 = off + seg.end
            hover = f"{name}<br>chr={chrom}<br>start={seg.start:,}<br>end={seg.end:,}<br>CN={seg.cn:g}"
            x.extend([x0, x1, None])
            y.extend([seg.cn, seg.cn, None])
            txt.extend([hover, hover, None])
    fig.add_trace(
        go.Scattergl(
            x=x,
            y=y,
            mode="lines",
            name=name,
            line=dict(color=color, width=2),
            text=txt,
            hovertemplate="%{text}<extra></extra>",
            yaxis="y",
        )
    )


def add_timing_traces(fig: go.Figure, chains: List[TimingChain], offsets: Dict[str, int], chain_mode: str, visible: bool) -> None:
    grouped: Dict[int, List[TimingChain]] = {}
    for c in chains:
        grouped.setdefault(c.chain_len, []).append(c)
    for chain_len in sorted(grouped):
        xs: List[Optional[int]] = []
        ys: List[Optional[float]] = []
        txt: List[Optional[str]] = []
        cds: List[Optional[Dict[str, Any]]] = []
        for c in grouped[chain_len]:
            off = offsets.get(c.chrom, 0)
            x0 = off + c.start
            x1 = off + c.end
            node_ids_text = " ".join(str(n.node_id) for n in c.nodes)
            hover = (
                f"Timing chain<br>chr={c.chrom}<br>start={c.start:,}<br>end={c.end:,}"
                f"<br>chain_levels={c.chain_len}<br>chain_index={c.chain_index}"
                f"<br>nodes_in_chain={len(c.nodes)}<br>total_mutations={c.total_mutations}"
            )
            xs.extend([x0, x1, None])
            ys.extend([float(c.chain_len), float(c.chain_len), None])
            txt.extend([hover, hover, None])
            cds.extend(
                [
                    {
                        "x0": x0,
                        "x1": x1,
                        "start_local": c.start,
                        "chain_len": c.chain_len,
                        "y": float(c.chain_len),
                        "file_rel": c.file_rel,
                        "chain_index": c.chain_index,
                        "node_ids_text": node_ids_text,
                    },
                    {
                        "x0": x0,
                        "x1": x1,
                        "start_local": c.start,
                        "chain_len": c.chain_len,
                        "y": float(c.chain_len),
                        "file_rel": c.file_rel,
                        "chain_index": c.chain_index,
                        "node_ids_text": node_ids_text,
                    },
                    None,
                ]
            )

        fig.add_trace(
            go.Scattergl(
                x=xs,
                y=ys,
                mode="lines",
                name=f"Timing chains (levels={chain_len})",
                line=dict(width=3),
                opacity=0.82,
                text=txt,
                customdata=cds,
                hovertemplate="%{text}<extra></extra>",
                visible=visible,
                meta={"trace_role": "timing_chain_line", "chain_mode": chain_mode},
                yaxis="y",
            )
        )


def add_rank_traces(fig: go.Figure, ranked_nodes: List[RankedNode], offsets: Dict[str, int]) -> None:
    if not ranked_nodes:
        return
    x: List[int] = []
    y: List[int] = []
    txt: List[str] = []
    total = 0
    for rn in ranked_nodes:
        off = offsets.get(rn.chrom, 0)
        for pos in rn.member_positions:
            x.append(off + pos)
            y.append(rn.rank)
            txt.append(
                f"Rank {rn.rank}<br>chr={rn.chrom}<br>pos={pos:,}<br>node_id={rn.node_id}<br>cluster_size={rn.member_count}"
            )
            total += 1
    fig.add_trace(
        go.Scattergl(
            x=x,
            y=y,
            mode="markers",
            name=f"Mutation ranks ({total:,} mutations)",
            marker=dict(size=4, symbol="diamond", color="#e377c2", opacity=0.6),
            text=txt,
            hovertemplate="%{text}<extra></extra>",
            meta={"trace_role": "rank_trace"},
            yaxis="y",
        )
    )


def representative_cn_for_span(chrom: str, start: int, end: int, by_chrom: Dict[str, List[CnSegment]]) -> Optional[float]:
    segs = by_chrom.get(chrom) or []
    if not segs or end < start:
        return None
    covered: Dict[float, int] = {}
    for s in segs:
        if math.isnan(s.cn):
            continue
        lo = max(start, s.start)
        hi = min(end, s.end)
        if hi < lo:
            continue
        covered[s.cn] = covered.get(s.cn, 0) + (hi - lo + 1)
    if not covered:
        return None
    return max(covered.items(), key=lambda kv: (kv[1], -kv[0]))[0]


def representative_chain_cn(chain: TimingChain, by_chrom: Dict[str, List[CnSegment]]) -> Optional[float]:
    return representative_cn_for_span(chain.chrom, chain.start, chain.end, by_chrom)


def add_chain_cn_violation_trace(
    fig: go.Figure,
    timing_chains: List[TimingChain],
    cn_max_haplotype: Dict[str, List[CnSegment]],
    offsets: Dict[str, int],
    chain_mode: str,
    visible: bool,
    showlegend: bool,
) -> None:
    xs: List[Optional[int]] = []
    ys: List[Optional[float]] = []
    texts: List[Optional[str]] = []
    n_viol = 0
    for c in timing_chains:
        rep = representative_chain_cn(c, cn_max_haplotype)
        if rep is None:
            continue
        if c.chain_len <= rep:
            continue
        off = offsets.get(c.chrom, 0)
        x0 = off + c.start
        x1 = off + c.end
        hover = (
            f"Chain > max CN<br>chr={c.chrom}<br>start={c.start:,}<br>end={c.end:,}"
            f"<br>chain_levels={c.chain_len}<br>rep_max_cn={rep:g}<br>chain_index={c.chain_index}"
        )
        xs.extend([x0, x1, None])
        ys.extend([0.0, 0.0, None])
        texts.extend([hover, hover, None])
        n_viol += 1
    fig.add_trace(
        go.Scattergl(
            x=xs,
            y=ys,
            mode="lines",
            name=f"Chain > max CN ({n_viol})",
            legendgroup="chain_cn_violation",
            showlegend=showlegend,
            line=dict(width=4, color="rgba(214,39,40,0.95)"),
            text=texts,
            hovertemplate="%{text}<extra></extra>",
            visible=visible,
            opacity=0.9,
            meta={"trace_role": "chain_cn_violation", "chain_mode": chain_mode},
            yaxis="y",
        )
    )


def add_chrom_guides(fig: go.Figure, chrom_lengths: Dict[str, int], offsets: Dict[str, int]) -> Tuple[List[int], List[str]]:
    tickvals: List[int] = []
    ticktext: List[str] = []
    for chrom in sorted(chrom_lengths, key=chrom_sort_key):
        off = offsets[chrom]
        end = off + chrom_lengths[chrom]
        mid = off + (chrom_lengths[chrom] // 2)
        tickvals.append(mid)
        ticktext.append(chrom)
        fig.add_vline(x=end, line_width=1, line_dash="dot", line_color="rgba(120,120,120,0.35)")
    return tickvals, ticktext


def _format_cn_state(cn: float) -> str:
    if not math.isfinite(cn):
        return ""
    if abs(cn - round(cn)) < 1e-9:
        return str(int(round(cn)))
    return f"{cn:.2f}".rstrip("0").rstrip(".")


def _format_mb_label(total_bp: int) -> str:
    mb = total_bp / 1_000_000.0
    if abs(mb - round(mb)) < 1e-9:
        return f"{int(round(mb))}mb"
    if mb >= 100:
        return f"{mb:.1f}mb"
    return f"{mb:.2f}mb"


def _sorted_cn_labels(labels: Iterable[str]) -> List[str]:
    vals: List[Tuple[float, str]] = []
    for lbl in labels:
        try:
            vals.append((float(lbl), lbl))
        except ValueError:
            continue
    vals.sort(key=lambda x: x[0])
    return [v for _, v in vals]


def compute_cn_state_total_bp(by_chrom: Dict[str, List[CnSegment]]) -> Dict[str, int]:
    totals: Dict[str, int] = {}
    for segs in by_chrom.values():
        for s in segs:
            if math.isnan(s.cn) or s.end < s.start:
                continue
            key = _format_cn_state(s.cn)
            totals[key] = totals.get(key, 0) + (s.end - s.start + 1)
    return totals


def build_chain_length_cn_distribution_figure(
    cn_hp1: Dict[str, List[CnSegment]],
    cn_hp2: Dict[str, List[CnSegment]],
    cn_total: Dict[str, List[CnSegment]],
    cn_max_haplotype: Dict[str, List[CnSegment]],
    timing_chains_all: List[TimingChain],
    timing_chains_non_overlap: List[TimingChain],
    timing_chains_longest_only: List[TimingChain],
    ranked_nodes: Optional[List[RankedNode]] = None,
    mutation_vaf_points: Optional[List[MutationVafPoint]] = None,
) -> go.Figure:
    mode_sources: Dict[str, Dict[str, List[CnSegment]]] = {
        "total": cn_total,
        "hp1": cn_hp1,
        "hp2": cn_hp2,
        "max": cn_max_haplotype,
    }
    mode_titles = {
        "total": "Total CN (HP1+HP2)",
        "hp1": "HP1 CN",
        "hp2": "HP2 CN",
        "max": "Max CN (max(HP1,HP2))",
    }
    mode_colors = {
        "total": "#2ca02c",
        "hp1": "#1f77b4",
        "hp2": "#d62728",
        "max": "#9467bd",
    }
    mode_state_bp: Dict[str, Dict[str, int]] = {m: compute_cn_state_total_bp(src) for m, src in mode_sources.items()}

    # row tuple: (cn_label, y_value, cluster_key, mutation_weight)
    # - chain modes use mutation_weight=1 per chain
    # - rank mode uses mutation_weight=member_count (uncondensed mutations in condensed node)
    chain_mode_rows: Dict[str, Dict[str, List[Tuple[str, int, str, int]]]] = {
        "all": {k: [] for k in mode_sources},
        "non_overlapping": {k: [] for k in mode_sources},
        "longest_only": {k: [] for k in mode_sources},
        "rank": {k: [] for k in mode_sources},
    }
    chain_mode_labels: Dict[str, Dict[str, List[str]]] = {"all": {}, "non_overlapping": {}, "longest_only": {}, "rank": {}}

    chain_mode_chains: Dict[str, List[TimingChain]] = {
        "all": timing_chains_all,
        "non_overlapping": timing_chains_non_overlap,
        "longest_only": timing_chains_longest_only,
    }
    for cm, chains in chain_mode_chains.items():
        for mode, by_chrom in mode_sources.items():
            rows: List[Tuple[str, int, str, int]] = []
            for c in chains:
                cn_val = representative_chain_cn(c, by_chrom)
                if cn_val is None:
                    continue
                ck = f"{c.file_rel}|cc:{c.component_id}" if c.component_id is not None else f"{c.file_rel}|na:{c.chrom}"
                rows.append((_format_cn_state(cn_val), c.chain_len, ck, 1))
            chain_mode_rows[cm][mode] = rows
            raw_labels = _sorted_cn_labels({x for x, _, _, _ in rows})
            chain_mode_labels[cm][mode] = [f"{x} ({_format_mb_label(mode_state_bp[mode].get(x, 0))})" for x in raw_labels]

    for mode, by_chrom in mode_sources.items():
        rows: List[Tuple[str, int, str, int]] = []
        for rn in (ranked_nodes or []):
            cn_val = representative_cn_for_span(rn.chrom, rn.start, rn.end, by_chrom)
            if cn_val is None:
                continue
            ck = f"{rn.file_rel}|cc:{rn.component_id}" if rn.component_id is not None else f"{rn.file_rel}|na:{rn.chrom}"
            rows.append((_format_cn_state(cn_val), rn.rank, ck, max(1, int(rn.member_count))))
        chain_mode_rows["rank"][mode] = rows
        raw_labels = _sorted_cn_labels({x for x, _, _, _ in rows})
        chain_mode_labels["rank"][mode] = [f"{x} ({_format_mb_label(mode_state_bp[mode].get(x, 0))})" for x in raw_labels]

    fig = go.Figure()
    cn_modes = ["total", "hp1", "hp2", "max"]
    chain_modes = ["all", "non_overlapping", "longest_only", "rank"]
    chain_mode_titles = {
        "all": "All chains",
        "non_overlapping": "Non-overlapping chains",
        "longest_only": "Longest chain per cluster",
        "rank": "Rank",
    }
    default_cn_mode = "total"
    default_chain_mode = "all"

    visible: List[bool] = []
    for cm in chain_modes:
        for cn_mode in cn_modes:
            rows = chain_mode_rows[cm][cn_mode]
            mode_bp = mode_state_bp[cn_mode]
            component_counts: Dict[Tuple[str, int], int] = {}
            mutation_counts: Dict[Tuple[str, int], int] = {}
            cluster_sets: Dict[Tuple[str, int], set[str]] = {}
            cluster_mutation_totals: Dict[Tuple[str, int], Dict[str, int]] = {}
            for cn_label, yv, cluster_key, mut_w in rows:
                k = (cn_label, yv)
                component_counts[k] = component_counts.get(k, 0) + 1
                mutation_counts[k] = mutation_counts.get(k, 0) + int(mut_w)
                cluster_sets.setdefault(k, set()).add(cluster_key)
                cluster_mutation_totals.setdefault(k, {})
                cluster_mutation_totals[k][cluster_key] = cluster_mutation_totals[k].get(cluster_key, 0) + int(mut_w)

            keys = list(component_counts.keys())
            count_x = [f"{k[0]} ({_format_mb_label(mode_bp.get(k[0], 0))})" for k in keys]
            count_y = [k[1] for k in keys]
            component_n = [component_counts[k] for k in keys]
            mutation_n = [mutation_counts[k] for k in keys]
            cluster_n = [len(cluster_sets.get(k, set())) for k in keys]
            cn_10mb = [(mode_bp.get(k[0], 0) / 10_000_000.0) for k in keys]
            component_n_norm = [(n / u if u > 0 else 0.0) for n, u in zip(component_n, cn_10mb)]
            mutation_n_norm = [(n / u if u > 0 else 0.0) for n, u in zip(mutation_n, cn_10mb)]
            cluster_n_norm = [(n / u if u > 0 else 0.0) for n, u in zip(cluster_n, cn_10mb)]
            y_axis_label = "Rank" if cm == "rank" else "Chain length"
            y_counts_label = "rank counts" if cm == "rank" else "chain counts"

            hover_raw = [
                (
                    f"CN={x}<br>{y_axis_label}={y}<br>"
                    f"Raw components={n_comp}<br>Raw clusters={n_clu}<br>"
                    f"Components/10Mb={n_comp_n:.4g}<br>Clusters/10Mb={n_clu_n:.4g}"
                )
                for x, y, n_comp, n_clu, n_comp_n, n_clu_n in zip(
                    count_x, count_y, component_n, cluster_n, component_n_norm, cluster_n_norm
                )
            ]
            hover_norm = [
                (
                    f"CN={x}<br>{y_axis_label}={y}<br>"
                    f"Components/10Mb={n_comp_n:.4g}<br>Clusters/10Mb={n_clu_n:.4g}<br>"
                    f"Raw components={n_comp}<br>Raw clusters={n_clu}"
                )
                for x, y, n_comp, n_clu, n_comp_n, n_clu_n in zip(
                    count_x, count_y, component_n, cluster_n, component_n_norm, cluster_n_norm
                )
            ]
            sizes_raw = [8 + 4 * math.log2(max(n, 0) + 1) for n in component_n]
            sizes_norm = [8 + 4 * math.log2(max(n, 0.0) + 1) for n in component_n_norm]
            mutation_sizes_raw = [8 + 4 * math.log2(max(n, 0) + 1) for n in mutation_n]
            mutation_sizes_norm = [8 + 4 * math.log2(max(n, 0.0) + 1) for n in mutation_n_norm]
            cluster_sizes_raw = [8 + 4 * math.log2(max(n, 0) + 1) for n in cluster_n]
            cluster_sizes_norm = [8 + 4 * math.log2(max(n, 0.0) + 1) for n in cluster_n_norm]
            rank_cluster_threshold_points = []
            if cm == "rank":
                for k, x, y, bp, comp_raw, comp_nrm in zip(keys, count_x, count_y, cn_10mb, component_n, component_n_norm):
                    per_cluster = cluster_mutation_totals.get(k, {})
                    rank_cluster_threshold_points.append(
                        {
                            "x": x,
                            "y": y,
                            "cn_bp_over_10mb": bp,
                            "cluster_mut_counts": sorted(int(v) for v in per_cluster.values()),
                            "component_raw": int(comp_raw),
                            "component_norm": float(comp_nrm),
                        }
                    )

            is_default = cm == default_chain_mode and cn_mode == default_cn_mode
            # chains raw
            fig.add_trace(
                go.Scatter(
                    x=count_x,
                    y=count_y,
                    mode="markers+text",
                    text=[str(n) for n in component_n],
                    textposition="middle right",
                    name=f"{mode_titles[cn_mode]} {y_counts_label}",
                    marker=dict(size=sizes_raw, color=mode_colors[cn_mode], line=dict(width=1, color="rgba(0,0,0,0.45)"), opacity=0.88),
                    hovertemplate="%{text}<br>%{customdata}<extra></extra>",
                    customdata=hover_raw,
                    visible=is_default,
                    meta={
                        "trace_role": "dist_count",
                        "dist_count_metric": "chains",
                        "dist_count_norm": "raw",
                        "dist_cn_mode": cn_mode,
                        "dist_chain_mode": cm,
                    },
                )
            )
            visible.append(is_default)
            # chains per 10mb
            fig.add_trace(
                go.Scatter(
                    x=count_x,
                    y=count_y,
                    mode="markers+text",
                    text=[f"{n:.3g}" for n in component_n_norm],
                    textposition="middle right",
                    name=f"{mode_titles[cn_mode]} {y_counts_label} (per 10Mb CN)",
                    marker=dict(size=sizes_norm, color=mode_colors[cn_mode], line=dict(width=1, color="rgba(0,0,0,0.45)"), opacity=0.88),
                    hovertemplate="%{text}<br>%{customdata}<extra></extra>",
                    customdata=hover_norm,
                    visible=False,
                    meta={
                        "trace_role": "dist_count",
                        "dist_count_metric": "chains",
                        "dist_count_norm": "per_10mb",
                        "dist_cn_mode": cn_mode,
                        "dist_chain_mode": cm,
                    },
                )
            )
            visible.append(False)
            # mutations raw (rank mode only)
            if cm == "rank":
                hover_mut_raw = [
                    (
                        f"CN={x}<br>{y_axis_label}={y}<br>"
                        f"Raw mutations={n_mut}<br>Mutations/10Mb={n_mut_n:.4g}"
                    )
                    for x, y, n_mut, n_mut_n in zip(count_x, count_y, mutation_n, mutation_n_norm)
                ]
                fig.add_trace(
                    go.Scatter(
                        x=count_x,
                        y=count_y,
                        mode="markers+text",
                        text=[str(n) for n in mutation_n],
                        textposition="middle right",
                        name=f"{mode_titles[cn_mode]} mutation counts",
                        marker=dict(size=mutation_sizes_raw, color=mode_colors[cn_mode], line=dict(width=1, color="rgba(0,0,0,0.45)"), opacity=0.88),
                        hovertemplate="%{text}<br>%{customdata}<extra></extra>",
                        customdata=hover_mut_raw,
                        visible=False,
                        meta={
                            "trace_role": "dist_count",
                            "dist_count_metric": "mutations",
                            "dist_count_norm": "raw",
                            "dist_cn_mode": cn_mode,
                            "dist_chain_mode": cm,
                        },
                    )
                )
                visible.append(False)
                hover_mut_norm = [
                    (
                        f"CN={x}<br>{y_axis_label}={y}<br>"
                        f"Mutations/10Mb={n_mut_n:.4g}<br>Raw mutations={n_mut}"
                    )
                    for x, y, n_mut, n_mut_n in zip(count_x, count_y, mutation_n, mutation_n_norm)
                ]
                fig.add_trace(
                    go.Scatter(
                        x=count_x,
                        y=count_y,
                        mode="markers+text",
                        text=[f"{n:.3g}" for n in mutation_n_norm],
                        textposition="middle right",
                        name=f"{mode_titles[cn_mode]} mutation counts (per 10Mb CN)",
                        marker=dict(size=mutation_sizes_norm, color=mode_colors[cn_mode], line=dict(width=1, color="rgba(0,0,0,0.45)"), opacity=0.88),
                        hovertemplate="%{text}<br>%{customdata}<extra></extra>",
                        customdata=hover_mut_norm,
                        visible=False,
                        meta={
                            "trace_role": "dist_count",
                            "dist_count_metric": "mutations",
                            "dist_count_norm": "per_10mb",
                            "dist_cn_mode": cn_mode,
                            "dist_chain_mode": cm,
                        },
                    )
                )
                visible.append(False)
            # clusters raw
            fig.add_trace(
                go.Scatter(
                    x=count_x,
                    y=count_y,
                    mode="markers+text",
                    text=[str(n) for n in cluster_n],
                    textposition="middle right",
                    name=f"{mode_titles[cn_mode]} cluster counts",
                    marker=dict(size=cluster_sizes_raw, color=mode_colors[cn_mode], line=dict(width=1, color="rgba(0,0,0,0.45)"), opacity=0.88),
                    hovertemplate="%{text}<br>%{customdata}<extra></extra>",
                    customdata=hover_raw,
                    visible=False,
                    meta={
                        "trace_role": "dist_count",
                        "dist_count_metric": "clusters",
                        "dist_count_norm": "raw",
                        "dist_cn_mode": cn_mode,
                        "dist_chain_mode": cm,
                        "rank_cluster_threshold_points": rank_cluster_threshold_points if cm == "rank" else [],
                    },
                )
            )
            visible.append(False)
            # clusters per 10mb
            fig.add_trace(
                go.Scatter(
                    x=count_x,
                    y=count_y,
                    mode="markers+text",
                    text=[f"{n:.3g}" for n in cluster_n_norm],
                    textposition="middle right",
                    name=f"{mode_titles[cn_mode]} cluster counts (per 10Mb CN)",
                    marker=dict(size=cluster_sizes_norm, color=mode_colors[cn_mode], line=dict(width=1, color="rgba(0,0,0,0.45)"), opacity=0.88),
                    hovertemplate="%{text}<br>%{customdata}<extra></extra>",
                    customdata=hover_norm,
                    visible=False,
                    meta={
                        "trace_role": "dist_count",
                        "dist_count_metric": "clusters",
                        "dist_count_norm": "per_10mb",
                        "dist_cn_mode": cn_mode,
                        "dist_chain_mode": cm,
                        "rank_cluster_threshold_points": rank_cluster_threshold_points if cm == "rank" else [],
                    },
                )
            )
            visible.append(False)

    fig.update_layout(
        title=(
            "Chain length distribution by CN state "
            f"({mode_titles[default_cn_mode]}, {chain_mode_titles[default_chain_mode]})"
        ),
        template="plotly_white",
        hovermode="closest",
        legend=dict(orientation="h", yanchor="bottom", y=1.01, xanchor="left", x=0.01),
        margin=dict(l=70, r=40, t=90, b=70),
        xaxis=dict(
            title="Copy-number state",
            type="category",
            categoryorder="array",
            categoryarray=chain_mode_labels[default_chain_mode][default_cn_mode],
            showgrid=False,
        ),
        yaxis=dict(title="Chain length (number of levels)", rangemode="tozero"),
        meta={
            "dist_cn_mode_titles": mode_titles,
            "dist_chain_mode_titles": chain_mode_titles,
            "dist_mode_labels": chain_mode_labels,
            "dist_default_cn_mode": default_cn_mode,
            "dist_default_chain_mode": default_chain_mode,
            "dist_default_count_metric": "chains",
            "dist_default_count_norm": "raw",
        },
    )
    for i, v in enumerate(visible):
        fig.data[i].visible = v
    add_vaf_distribution_trace(fig, mutation_vaf_points)
    return fig


def build_chain_length_distribution_no_cn_figure(
    timing_chains_all: List[TimingChain],
    timing_chains_non_overlap: List[TimingChain],
    timing_chains_longest_only: List[TimingChain],
    ranked_nodes: Optional[List[RankedNode]] = None,
    mutation_vaf_points: Optional[List[MutationVafPoint]] = None,
) -> go.Figure:
    chain_mode_rows: Dict[str, List[Tuple[int, str, int]]] = {
        "all": [],
        "non_overlapping": [],
        "longest_only": [],
        "rank": [],
    }
    chain_mode_titles = {
        "all": "All chains",
        "non_overlapping": "Non-overlapping chains",
        "longest_only": "Longest chain per cluster",
        "rank": "Rank",
    }
    chain_mode_x_titles = {
        "all": "Chain length",
        "non_overlapping": "Chain length",
        "longest_only": "Chain length",
        "rank": "Rank",
    }
    chain_mode_y_titles = {
        "chains": "# chains",
        "clusters": "# clusters",
        "mutations": "# mutations",
    }

    chain_mode_chains: Dict[str, List[TimingChain]] = {
        "all": timing_chains_all,
        "non_overlapping": timing_chains_non_overlap,
        "longest_only": timing_chains_longest_only,
    }
    for cm, chains in chain_mode_chains.items():
        rows: List[Tuple[int, str, int]] = []
        for c in chains:
            ck = f"{c.file_rel}|cc:{c.component_id}" if c.component_id is not None else f"{c.file_rel}|na:{c.chrom}"
            rows.append((c.chain_len, ck, 1))
        chain_mode_rows[cm] = rows

    rank_rows: List[Tuple[int, str, int]] = []
    for rn in (ranked_nodes or []):
        ck = f"{rn.file_rel}|cc:{rn.component_id}" if rn.component_id is not None else f"{rn.file_rel}|na:{rn.chrom}"
        rank_rows.append((rn.rank, ck, max(1, int(rn.member_count))))
    chain_mode_rows["rank"] = rank_rows

    fig = go.Figure()
    chain_modes = ["all", "non_overlapping", "longest_only", "rank"]
    default_chain_mode = "all"
    default_metric = "chains"
    color_by_mode = {
        "all": "#2ca02c",
        "non_overlapping": "#1f77b4",
        "longest_only": "#ff7f0e",
        "rank": "#9467bd",
    }

    visible: List[bool] = []
    for cm in chain_modes:
        rows = chain_mode_rows[cm]
        component_counts: Dict[int, int] = {}
        mutation_counts: Dict[int, int] = {}
        cluster_sets: Dict[int, set[str]] = {}
        for xv, cluster_key, mut_w in rows:
            component_counts[xv] = component_counts.get(xv, 0) + 1
            mutation_counts[xv] = mutation_counts.get(xv, 0) + int(mut_w)
            cluster_sets.setdefault(xv, set()).add(cluster_key)
        keys = sorted(component_counts.keys())
        count_x = keys
        component_n = [component_counts[k] for k in keys]
        mutation_n = [mutation_counts[k] for k in keys]
        cluster_n = [len(cluster_sets.get(k, set())) for k in keys]

        hover_chain = [f"{chain_mode_x_titles[cm]}={x}<br># chains={n}" for x, n in zip(count_x, component_n)]
        hover_cluster = [f"{chain_mode_x_titles[cm]}={x}<br># clusters={n}" for x, n in zip(count_x, cluster_n)]
        hover_mut = [f"{chain_mode_x_titles[cm]}={x}<br># mutations={n}" for x, n in zip(count_x, mutation_n)]
        chain_sizes = [8 + 4 * math.log2(max(n, 0) + 1) for n in component_n]
        cluster_sizes = [8 + 4 * math.log2(max(n, 0) + 1) for n in cluster_n]
        mutation_sizes = [8 + 4 * math.log2(max(n, 0) + 1) for n in mutation_n]

        is_default = cm == default_chain_mode
        fig.add_trace(
            go.Scatter(
                x=count_x,
                y=component_n,
                mode="markers+text",
                text=[str(n) for n in component_n],
                textposition="top center",
                name=f"{chain_mode_titles[cm]} # chains",
                marker=dict(
                    size=chain_sizes,
                    color=color_by_mode[cm],
                    line=dict(width=1, color="rgba(0,0,0,0.45)"),
                    opacity=0.88,
                ),
                hovertemplate="%{customdata}<extra></extra>",
                customdata=hover_chain,
                visible=is_default,
                meta={
                    "trace_role": "dist_count",
                    "dist_count_metric": "chains",
                    "dist_count_norm": "raw",
                    "dist_cn_mode": "none",
                    "dist_chain_mode": cm,
                },
            )
        )
        visible.append(is_default)

        fig.add_trace(
            go.Scatter(
                x=count_x,
                y=cluster_n,
                mode="markers+text",
                text=[str(n) for n in cluster_n],
                textposition="top center",
                name=f"{chain_mode_titles[cm]} # clusters",
                marker=dict(
                    size=cluster_sizes,
                    color=color_by_mode[cm],
                    line=dict(width=1, color="rgba(0,0,0,0.45)"),
                    opacity=0.88,
                ),
                hovertemplate="%{customdata}<extra></extra>",
                customdata=hover_cluster,
                visible=False,
                meta={
                    "trace_role": "dist_count",
                    "dist_count_metric": "clusters",
                    "dist_count_norm": "raw",
                    "dist_cn_mode": "none",
                    "dist_chain_mode": cm,
                },
            )
        )
        visible.append(False)

        if cm == "rank":
            fig.add_trace(
                go.Scatter(
                    x=count_x,
                    y=mutation_n,
                    mode="markers+text",
                    text=[str(n) for n in mutation_n],
                    textposition="top center",
                    name="Rank # mutations",
                    marker=dict(
                        size=mutation_sizes,
                        color=color_by_mode[cm],
                        line=dict(width=1, color="rgba(0,0,0,0.45)"),
                        opacity=0.88,
                    ),
                    hovertemplate="%{customdata}<extra></extra>",
                    customdata=hover_mut,
                    visible=False,
                    meta={
                        "trace_role": "dist_count",
                        "dist_count_metric": "mutations",
                        "dist_count_norm": "raw",
                        "dist_cn_mode": "none",
                        "dist_chain_mode": cm,
                    },
                )
            )
            visible.append(False)

    fig.update_layout(
        title="Chain-length distribution (All chains, #chains)",
        template="plotly_white",
        hovermode="closest",
        legend=dict(orientation="h", yanchor="bottom", y=1.01, xanchor="left", x=0.01),
        margin=dict(l=70, r=40, t=90, b=70),
        xaxis=dict(title="Chain length", showgrid=False),
        yaxis=dict(title="# chains", rangemode="tozero"),
        meta={
            "dist_cn_mode_titles": {},
            "dist_chain_mode_titles": chain_mode_titles,
            "dist_mode_labels": {},
            "dist_default_cn_mode": "none",
            "dist_default_chain_mode": default_chain_mode,
            "dist_default_count_metric": default_metric,
            "dist_default_count_norm": "raw",
            "dist_has_cn": False,
            "dist_chain_mode_x_titles": chain_mode_x_titles,
            "dist_count_metric_y_titles": chain_mode_y_titles,
        },
    )
    for i, v in enumerate(visible):
        fig.data[i].visible = v
    add_vaf_distribution_trace(fig, mutation_vaf_points)
    return fig


def build_figure(
    cn_hp1: Dict[str, List[CnSegment]],
    cn_hp2: Dict[str, List[CnSegment]],
    cn_total: Dict[str, List[CnSegment]],
    cn_max_haplotype: Dict[str, List[CnSegment]],
    timing_chains: List[TimingChain],
    title: str,
    ranked_nodes: Optional[List[RankedNode]] = None,
    show_cn: bool = True,
) -> go.Figure:
    chrom_lengths = get_chrom_lengths(cn_hp1, cn_hp2, cn_total, timing_chains, ranked_nodes)
    if not chrom_lengths:
        raise ValueError("No CN or timing segment data available to plot.")
    offsets = chrom_offsets(chrom_lengths)
    fig = go.Figure()

    if show_cn:
        add_cn_trace(fig, "HP1 copy number", cn_hp1, offsets, color="#1f77b4")
        add_cn_trace(fig, "HP2 copy number", cn_hp2, offsets, color="#d62728")
        add_cn_trace(fig, "Total copy number (HP1+HP2)", cn_total, offsets, color="#2ca02c")
        add_cn_trace(fig, "Max haplotype copy number (max(HP1,HP2))", cn_max_haplotype, offsets, color="#9467bd")

    chains_non = select_non_overlapping_chains(timing_chains)
    chains_long = select_longest_chain_per_cluster(timing_chains)
    add_timing_traces(fig, timing_chains, offsets, chain_mode="all", visible=True)
    add_timing_traces(fig, chains_non, offsets, chain_mode="non_overlapping", visible=False)
    add_timing_traces(fig, chains_long, offsets, chain_mode="longest_only", visible=False)

    if show_cn:
        add_chain_cn_violation_trace(
            fig,
            timing_chains,
            cn_max_haplotype,
            offsets,
            chain_mode="all",
            visible=True,
            showlegend=True,
        )
        add_chain_cn_violation_trace(
            fig,
            chains_non,
            cn_max_haplotype,
            offsets,
            chain_mode="non_overlapping",
            visible=False,
            showlegend=False,
        )
        add_chain_cn_violation_trace(
            fig,
            chains_long,
            cn_max_haplotype,
            offsets,
            chain_mode="longest_only",
            visible=False,
            showlegend=False,
        )

    add_rank_traces(fig, ranked_nodes or [], offsets)

    fig.add_trace(
        go.Scattergl(
            x=[None],
            y=[None],
            mode="markers",
            name="Timing markers",
            marker=dict(symbol="diamond-open", size=7, line=dict(width=2.4, color="rgba(80,80,80,0.95)"), color="rgba(80,80,80,0.2)"),
            hoverinfo="skip",
            meta={"trace_role": "marker_toggle_legend"},
            yaxis="y",
        )
    )

    tickvals, ticktext = add_chrom_guides(fig, chrom_lengths, offsets)
    genome_x_max = max((offsets[c] + chrom_lengths[c] for c in chrom_lengths), default=0)
    fig.update_layout(
        title=title,
        template="plotly_white",
        hovermode="closest",
        dragmode="zoom",
        legend=dict(orientation="h", yanchor="bottom", y=0.985, xanchor="left", x=0.01, groupclick="togglegroup"),
        margin=dict(l=70, r=70, t=110, b=70),
        xaxis=dict(
            title="Genome coordinate (chromosomes concatenated)",
            tickmode="array",
            tickvals=tickvals,
            ticktext=ticktext,
            showgrid=False,
            range=[0, genome_x_max],
        ),
        yaxis=dict(
            title="Shared scale: copy number and timing chain levels"
            if show_cn
            else "Timing chain levels / rank",
            rangemode="tozero",
        ),
    )
    return fig


def write_interactive_html(
    fig: go.Figure,
    fig_dist: go.Figure,
    out_html: Path,
    outdir: Path,
    chrom_lengths: Dict[str, int],
    offsets: Dict[str, int],
    has_cn: bool,
) -> None:
    div_id = "timing-cn-plot"
    dist_div_id = "timing-chain-dist-plot"
    plot_html = fig.to_html(include_plotlyjs="cdn", full_html=False, div_id=div_id)
    dist_html = fig_dist.to_html(include_plotlyjs=False, full_html=False, div_id=dist_div_id)
    region_lookup_json = json.dumps({"lengths": chrom_lengths, "offsets": offsets})
    outdir_json = json.dumps(str(outdir))

    bottom_title = (
        "Chain length distribution by CN state" if has_cn else "Chain/rank distribution"
    )
    if has_cn:
        dist_xaxis_controls = """
    <label for="dist-xaxis-mode-select">x-axis:</label>
    <select id="dist-xaxis-mode-select" onchange="setDistributionMode()">
      <option value="cn" selected>CN state</option>
      <option value="vaf">VAF (mutations)</option>
    </select>
"""
        dist_cn_controls = """
    <span id="dist-cn-controls">
    <label for="dist-cn-mode-select">x-axis CN:</label>
    <select id="dist-cn-mode-select" onchange="setDistributionMode()">
      <option value="total" selected>Total CN (HP1+HP2)</option>
      <option value="hp1">HP1 CN</option>
      <option value="hp2">HP2 CN</option>
      <option value="max">Max CN (max(HP1,HP2))</option>
    </select>
    </span>
"""
        dist_norm_controls = """
    <span id="dist-norm-controls">
    <label for="dist-count-normalize-select">count scale:</label>
    <select id="dist-count-normalize-select" onchange="setDistributionMode()">
      <option value="raw" selected>Raw</option>
      <option value="per_10mb">Normalized per 10Mb CN</option>
    </select>
    </span>
"""
    else:
        dist_xaxis_controls = """
    <label for="dist-xaxis-mode-select">x-axis:</label>
    <select id="dist-xaxis-mode-select" onchange="setDistributionMode()">
      <option value="mode" selected>Chain length / rank</option>
      <option value="vaf">VAF (mutations)</option>
    </select>
"""
        dist_cn_controls = ""
        dist_norm_controls = ""

    page = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Timing/CN interactive plot</title>
  <style>
    body {{ margin: 0; font-family: Arial, sans-serif; }}
    .controls {{
      display: flex; gap: 8px; align-items: center; flex-wrap: wrap;
      padding: 10px 12px; border-bottom: 1px solid #ddd; background: #fafafa;
    }}
    .controls input {{ min-width: 240px; padding: 6px 8px; font-size: 14px; }}
    .controls button {{ padding: 6px 10px; font-size: 14px; cursor: pointer; }}
    #search-status, #selected-chain-status {{ font-size: 13px; color: #444; }}
  </style>
</head>
<body>
  <div class="controls">
    <label for="node-id-input">node_id:</label>
    <input id="node-id-input" type="text" placeholder="Enter node_id and click Zoom" />
    <button type="button" onclick="zoomToNodeId()">Zoom</button>
    <label for="region-input">region:</label>
    <input id="region-input" type="text" placeholder="chrX:100000-200000" />
    <button type="button" onclick="zoomToRegion()">Zoom region</button>
    <label for="chain-mode-select">chain mode:</label>
    <select id="chain-mode-select" onchange="setChainMode(this.value)">
      <option value="all" selected>All chains</option>
      <option value="non_overlapping">Non-overlapping chains</option>
      <option value="longest_only">Longest chain per cluster</option>
    </select>
    <button type="button" onclick="openSelectedChainInCytoscape()">Open selected chain in Cytoscape</button>
    <button type="button" onclick="resetAxes()">Reset view</button>
    <span id="selected-chain-status">No chain selected.</span>
    <span id="search-status"></span>
  </div>
  {plot_html}
  <div style="padding: 10px 12px 0 12px; font-size: 15px; font-weight: 600;">{bottom_title}</div>
  <div class="controls" style="border-top: 1px solid #eee; border-bottom: 0;">
    {dist_xaxis_controls}
    {dist_cn_controls}
    <span id="dist-chain-mode-controls">
      <label id="dist-chain-mode-label" for="dist-chain-mode-select">count:</label>
      <select id="dist-chain-mode-select" onchange="setDistributionMode()">
        <option value="all" selected>All chains</option>
        <option value="non_overlapping">Non-overlapping chains</option>
        <option value="longest_only">Longest chain per cluster</option>
        <option value="rank">Rank</option>
      </select>
    </span>
    <span id="dist-count-metric-controls">
      <label for="dist-count-mode-select">circles:</label>
      <select id="dist-count-mode-select" onchange="setDistributionMode()">
        <option value="chains" selected># chains</option>
        <option value="clusters"># clusters</option>
        <option id="dist-count-mode-mutations" value="mutations" hidden># mutations</option>
      </select>
    </span>
    <span id="dist-rank-threshold-controls">
      <label id="dist-rank-cluster-min-muts-label" for="dist-rank-cluster-min-muts-input" style="display:none;">cluster min mutations:</label>
      <input
        id="dist-rank-cluster-min-muts-input"
        type="number"
        min="1"
        step="1"
        value="1"
        style="display:none; width: 90px; min-width: 90px;"
        oninput="setDistributionMode()"
      />
    </span>
    {dist_norm_controls}
  </div>
  <div style="padding: 0 8px 8px 8px;">
    {dist_html}
  </div>
  <script>
    let markersVisible = true;
    let violationsVisible = true;
    let activeChainMode = "all";
    let selectedChain = null;
    const CYTO_BRIDGE_OPEN_URL = "http://127.0.0.1:8765/open_chain";
    const regionLookup = {region_lookup_json};
    const plotOutdir = {outdir_json};

    function setStatus(msg) {{
      const el = document.getElementById("search-status");
      if (el) el.textContent = msg;
    }}
    function setSelectedChainStatus(msg) {{
      const el = document.getElementById("selected-chain-status");
      if (el) el.textContent = msg;
    }}

    function canonicalChromLabel(raw) {{
      const s = String(raw || "").trim();
      if (!s) return "";
      let c = s.replace(/^chr/i, "").toUpperCase();
      if (c === "M" || c === "MITO") c = "MT";
      return c;
    }}
    function parseIntLoose(text) {{
      const cleaned = String(text || "").replace(/,/g, "").trim();
      if (!/^-?\\d+$/.test(cleaned)) return null;
      const n = Number(cleaned);
      return Number.isFinite(n) ? n : null;
    }}
    function parseRegionSpec(spec) {{
      const raw = String(spec || "").trim();
      const m = raw.match(/^([^:]+):([^\\-]+)-([^\\-]+)$/);
      if (!m) return null;
      const chrom = canonicalChromLabel(m[1]);
      const c1 = parseIntLoose(m[2]);
      const c2 = parseIntLoose(m[3]);
      if (!chrom || c1 === null || c2 === null) return null;
      return {{ chrom, start: Math.min(c1, c2), end: Math.max(c1, c2) }};
    }}

    function setMarkersVisible(gd, showMarkers) {{
      if (!gd || !Array.isArray(gd.data)) return;
      const update = {{ visible: [] }};
      const traceIdx = [];
      for (let i = 0; i < gd.data.length; i++) {{
        const tr = gd.data[i];
        if (!tr || !tr.meta || tr.meta.trace_role !== "timing_marker") continue;
        const mode = tr.meta.chain_mode;
        if (mode && mode !== activeChainMode) continue;
        traceIdx.push(i);
        update.visible.push(showMarkers ? true : false);
      }}
      if (traceIdx.length) Plotly.restyle(gd, update, traceIdx);
    }}

    function setViolationsVisible(gd, showViolations) {{
      if (!gd || !Array.isArray(gd.data)) return;
      const update = {{ visible: [] }};
      const traceIdx = [];
      for (let i = 0; i < gd.data.length; i++) {{
        const tr = gd.data[i];
        if (!tr || !tr.meta || tr.meta.trace_role !== "chain_cn_violation") continue;
        const mode = tr.meta.chain_mode;
        if (mode && mode !== activeChainMode) continue;
        traceIdx.push(i);
        update.visible.push(showViolations ? true : false);
      }}
      if (traceIdx.length) Plotly.restyle(gd, update, traceIdx);
    }}

    function setChainMode(mode) {{
      const gd = document.getElementById("{div_id}");
      if (!gd || !Array.isArray(gd.data)) return;
      if (mode !== "all" && mode !== "non_overlapping" && mode !== "longest_only") return;
      activeChainMode = mode;
      const update = {{ visible: [] }};
      const traceIdx = [];
      for (let i = 0; i < gd.data.length; i++) {{
        const tr = gd.data[i];
        if (!tr || !tr.meta || !tr.meta.chain_mode) continue;
        const trMode = tr.meta.chain_mode;
        let show = trMode === mode;
        if (show && tr.meta.trace_role === "timing_marker" && !markersVisible) show = false;
        if (show && tr.meta.trace_role === "chain_cn_violation" && !violationsVisible) show = false;
        traceIdx.push(i);
        update.visible.push(show);
      }}
      if (traceIdx.length) Plotly.restyle(gd, update, traceIdx);
      // Ensure VAF-only traces are hidden in non-VAF mode.
      const vafIdx = [];
      for (let i = 0; i < gd.data.length; i++) {{
        const tr = gd.data[i];
        if (tr && tr.meta && tr.meta.trace_role === "dist_vaf") vafIdx.push(i);
      }}
      if (vafIdx.length) {{
        Plotly.restyle(gd, {{ visible: Array(vafIdx.length).fill(false) }}, vafIdx);
      }}
      if (mode === "all") setStatus("Showing all chains.");
      else if (mode === "non_overlapping") setStatus("Showing non-overlapping chains.");
      else setStatus("Showing longest chain per cluster.");
    }}

    const hasCn = {str(has_cn).lower()};

    function setDistributionMode() {{
      const gd = document.getElementById("{dist_div_id}");
      const cnSel = document.getElementById("dist-cn-mode-select");
      const xAxisModeSel = document.getElementById("dist-xaxis-mode-select");
      const chainSel = document.getElementById("dist-chain-mode-select");
      const chainSelLabel = document.getElementById("dist-chain-mode-label");
      const chainSelControls = document.getElementById("dist-chain-mode-controls");
      const countMetricControls = document.getElementById("dist-count-metric-controls");
      const rankThresholdControls = document.getElementById("dist-rank-threshold-controls");
      const cnControls = document.getElementById("dist-cn-controls");
      const normControls = document.getElementById("dist-norm-controls");
      const countSel = document.getElementById("dist-count-mode-select");
      const countMutOpt = document.getElementById("dist-count-mode-mutations");
      const clusterMinLbl = document.getElementById("dist-rank-cluster-min-muts-label");
      const clusterMinInput = document.getElementById("dist-rank-cluster-min-muts-input");
      const countNormSel = document.getElementById("dist-count-normalize-select");
      if (!gd || !chainSel || !countSel || !Array.isArray(gd.data)) return;
      const xAxisMode = String((xAxisModeSel && xAxisModeSel.value) || (hasCn ? "cn" : "mode"));
      const isVafMode = xAxisMode === "vaf";
      const cnMode = hasCn ? String((cnSel && cnSel.value) || "total") : "none";
      const chainMode = String(chainSel.value || "all");
      if (chainSelLabel) {{
        chainSelLabel.textContent = "count:";
      }}
      if (cnControls) cnControls.style.display = (hasCn && !isVafMode) ? "" : "none";
      if (normControls) normControls.style.display = (hasCn && !isVafMode) ? "" : "none";
      if (chainSelControls) chainSelControls.style.display = isVafMode ? "none" : "";
      if (countMetricControls) countMetricControls.style.display = isVafMode ? "none" : "";
      if (rankThresholdControls) rankThresholdControls.style.display = isVafMode ? "none" : "";
      if (isVafMode) {{
        const updateVaf = {{ visible: [] }};
        const traceIdxVaf = [];
        for (let i = 0; i < gd.data.length; i++) {{
          const tr = gd.data[i];
          if (!tr || !tr.meta) continue;
          if (tr.meta.trace_role !== "dist_vaf" && tr.meta.trace_role !== "dist_count") continue;
          traceIdxVaf.push(i);
          updateVaf.visible.push(tr.meta.trace_role === "dist_vaf");
        }}
        if (traceIdxVaf.length) Plotly.restyle(gd, updateVaf, traceIdxVaf);
        Plotly.relayout(gd, {{
          "title": "Mutation count distribution by VAF (uncondensed nodes; bin width=0.02)",
          "xaxis.title.text": "Variant allele frequency (VAF)",
          "yaxis.title.text": "Mutation count",
          "xaxis.type": "linear",
          "xaxis.range": [0, 1],
          "xaxis.tickmode": "linear",
          "xaxis.dtick": 0.1,
          "xaxis.tickformat": ".1f",
          "xaxis.categoryorder": null,
          "xaxis.categoryarray": null,
          "yaxis.rangemode": "tozero"
        }});
        return;
      }}

      // Non-VAF mode should never display VAF traces.
      const vafTraceIdx = [];
      for (let i = 0; i < gd.data.length; i++) {{
        const tr = gd.data[i];
        if (tr && tr.meta && tr.meta.trace_role === "dist_vaf") vafTraceIdx.push(i);
      }}
      if (vafTraceIdx.length) Plotly.restyle(gd, {{ visible: Array(vafTraceIdx.length).fill(false) }}, vafTraceIdx);

      if (countMutOpt) {{
        countMutOpt.hidden = chainMode !== "rank";
      }}
      if (chainMode !== "rank" && String(countSel.value || "") === "mutations") {{
        countSel.value = "chains";
      }}
      const countMetric = String(countSel.value || "chains");
      const countNorm = hasCn ? String((countNormSel && countNormSel.value) || "raw") : "raw";
      const showClusterRankThreshold = chainMode === "rank" && countMetric === "clusters";
      if (clusterMinLbl) clusterMinLbl.style.display = showClusterRankThreshold ? "" : "none";
      if (clusterMinInput) clusterMinInput.style.display = showClusterRankThreshold ? "" : "none";
      let clusterRankMinMut = 1;
      if (clusterMinInput) {{
        const parsed = Number.parseInt(String(clusterMinInput.value || "1"), 10);
        clusterRankMinMut = Number.isFinite(parsed) && parsed > 0 ? parsed : 1;
        if (String(clusterRankMinMut) !== String(clusterMinInput.value || "")) {{
          clusterMinInput.value = String(clusterRankMinMut);
        }}
      }}
      const meta = (gd.layout && gd.layout.meta) ? gd.layout.meta : {{}};
      const cnTitles = meta.dist_cn_mode_titles || {{}};
      const chainTitles = meta.dist_chain_mode_titles || {{}};
      const modeLabels = meta.dist_mode_labels || {{}};

      const update = {{ visible: [] }};
      const traceIdx = [];
      for (let i = 0; i < gd.data.length; i++) {{
        const tr = gd.data[i];
        if (!tr || !tr.meta) continue;
        const trCn = tr.meta.dist_cn_mode;
        const trChain = tr.meta.dist_chain_mode;
        const trMetric = tr.meta.dist_count_metric;
        const trNorm = tr.meta.dist_count_norm;
        if (!trCn || !trChain) continue;
        traceIdx.push(i);
        let show = trCn === cnMode && trChain === chainMode;
        if (show && tr.meta.trace_role === "dist_count") {{
          show = trMetric === countMetric && trNorm === countNorm;
        }}
        update.visible.push(show);
      }}
      if (traceIdx.length) Plotly.restyle(gd, update, traceIdx);

      // Rank + clusters supports an uncondensed-mutation threshold per cluster.
      if (showClusterRankThreshold) {{
        for (let i = 0; i < gd.data.length; i++) {{
          const tr = gd.data[i];
          if (!tr || !tr.meta) continue;
          if (tr.meta.trace_role !== "dist_count") continue;
          if (tr.meta.dist_chain_mode !== "rank") continue;
          if (tr.meta.dist_count_metric !== "clusters") continue;
          if (tr.meta.dist_cn_mode !== cnMode) continue;
          if (tr.meta.dist_count_norm !== countNorm) continue;
          if (tr.visible === false) continue;

          const pts = Array.isArray(tr.meta.rank_cluster_threshold_points) ? tr.meta.rank_cluster_threshold_points : [];
          const xVals = [];
          const yVals = [];
          const txtVals = [];
          const hoverVals = [];
          const sizeVals = [];
          for (const p of pts) {{
            if (!p) continue;
            const mutCounts = Array.isArray(p.cluster_mut_counts) ? p.cluster_mut_counts : [];
            let nRaw = 0;
            for (const v of mutCounts) {{
              const n = Number(v);
              if (Number.isFinite(n) && n >= clusterRankMinMut) nRaw += 1;
            }}
            if (nRaw <= 0) continue;
            const denom = Number(p.cn_bp_over_10mb || 0);
            const nShown = (countNorm === "per_10mb") ? (denom > 0 ? (nRaw / denom) : 0.0) : nRaw;
            xVals.push(p.x);
            yVals.push(p.y);
            txtVals.push(countNorm === "per_10mb" ? Number(nShown).toPrecision(3).replace(/\\.0+$/, "") : String(nRaw));
            hoverVals.push(
              "CN=" + String(p.x) + "<br>" +
              "Rank=" + String(p.y) + "<br>" +
              "Cluster min mutations=" + String(clusterRankMinMut) + "<br>" +
              "Clusters passing threshold=" + String(nRaw) + "<br>" +
              "Clusters/10Mb=" + (denom > 0 ? (nRaw / denom).toPrecision(4) : "0")
            );
            sizeVals.push(8 + 4 * Math.log2(Math.max(Number(nShown), 0) + 1));
          }}
          Plotly.restyle(
            gd,
            {{
              x: [xVals],
              y: [yVals],
              text: [txtVals],
              customdata: [hoverVals],
              "marker.size": [sizeVals],
            }},
            [i]
          );
        }}
      }}

      const cnTitle = cnTitles[cnMode] || cnMode;
      const chainTitle = chainTitles[chainMode] || chainMode;
      const metricTitle = countMetric === "clusters" ? "#clusters" : (countMetric === "mutations" ? "#mutations" : "#chains");
      const normTitle = countNorm === "per_10mb" ? "normalized per 10Mb CN" : "raw counts";
      const labels = (modeLabels[chainMode] || {{}})[cnMode] || [];
      const chainModeXTitles = meta.dist_chain_mode_x_titles || {{}};
      const metricYTitles = meta.dist_count_metric_y_titles || {{}};
      const yTitle = hasCn
        ? (chainMode === "rank" ? "Mutation rank" : "Chain length (number of levels)")
        : (metricYTitles[countMetric] || "# counts");
      const xTitle = hasCn
        ? "Copy-number state"
        : (chainModeXTitles[chainMode] || (chainMode === "rank" ? "Rank" : "Chain length"));
      const plotTitle = hasCn
        ? (chainMode === "rank"
        ? "Rank distribution by CN state (" + cnTitle + ", " + metricTitle + ", " + normTitle + (
            showClusterRankThreshold ? (", min mutations/cluster=" + String(clusterRankMinMut)) : ""
          ) + ")"
        : "Chain length distribution by CN state (" + cnTitle + ", " + chainTitle + ", " + metricTitle + ", " + normTitle + ")")
        : ((chainMode === "rank" ? "Rank" : chainTitle) + " distribution (" + metricTitle + ")");
      const relayoutUpdate = {{
        "title": plotTitle,
        "xaxis.title.text": xTitle,
        "yaxis.title.text": yTitle,
      }};
      if (hasCn) {{
        relayoutUpdate["xaxis.type"] = "category";
        relayoutUpdate["xaxis.categoryorder"] = "array";
        relayoutUpdate["xaxis.categoryarray"] = labels;
        relayoutUpdate["xaxis.tickmode"] = "auto";
        relayoutUpdate["xaxis.dtick"] = null;
        relayoutUpdate["xaxis.tickformat"] = null;
        relayoutUpdate["xaxis.range"] = null;
        relayoutUpdate["xaxis.autorange"] = true;
      }} else {{
        relayoutUpdate["xaxis.type"] = "linear";
        relayoutUpdate["xaxis.categoryorder"] = null;
        relayoutUpdate["xaxis.categoryarray"] = null;
        relayoutUpdate["xaxis.tickmode"] = "linear";
        relayoutUpdate["xaxis.dtick"] = 1;
        relayoutUpdate["xaxis.tickformat"] = "d";
        relayoutUpdate["xaxis.range"] = null;
        relayoutUpdate["xaxis.autorange"] = true;
      }}
      Plotly.relayout(gd, relayoutUpdate);
    }}

    async function copyTextToClipboard(text) {{
      const value = String(text);
      try {{
        if (navigator.clipboard && window.isSecureContext) {{
          await navigator.clipboard.writeText(value);
          return true;
        }}
      }} catch (_err) {{}}
      try {{
        const ta = document.createElement("textarea");
        ta.value = value;
        ta.style.position = "fixed";
        ta.style.left = "-9999px";
        document.body.appendChild(ta);
        ta.focus();
        ta.select();
        const ok = document.execCommand("copy");
        document.body.removeChild(ta);
        return !!ok;
      }} catch (_err) {{
        return false;
      }}
    }}

    async function openSelectedChainInCytoscape() {{
      if (!selectedChain || !selectedChain.file_rel) {{
        setStatus("Select a chain first by clicking a timing-chain segment.");
        return;
      }}
      const payload = {{
        outdir: plotOutdir,
        file_rel: selectedChain.file_rel,
        open_condensed: true,
        open_uncondensed: true
      }};
      try {{
        const resp = await fetch(CYTO_BRIDGE_OPEN_URL, {{
          method: "POST",
          headers: {{ "Content-Type": "application/json" }},
          body: JSON.stringify(payload)
        }});
        const data = await resp.json();
        if (!resp.ok) {{
          const msg = (data && data.error) ? data.error : ("HTTP " + resp.status);
          setStatus("Cytoscape bridge error: " + msg);
          return;
        }}
        const opened = Array.isArray(data.opened) ? data.opened.length : 0;
        setStatus("Opened " + opened + " network(s) in Cytoscape.");
      }} catch (_err) {{
        setStatus("Could not reach Cytoscape bridge at " + CYTO_BRIDGE_OPEN_URL);
      }}
    }}

    function collectNodeMatches(gd, nodeId) {{
      const out = [];
      for (const trace of gd.data || []) {{
        if (!trace || !Array.isArray(trace.customdata)) continue;
        if (trace.visible === false) continue;
        const xs = trace.x || [];
        const ys = trace.y || [];
        const cds = trace.customdata || [];
        const n = Math.min(xs.length, ys.length, cds.length);
        for (let i = 0; i < n; i++) {{
          const cd = cds[i];
          if (!cd || cd.node_id === undefined || cd.node_id === null) continue;
          if (String(cd.node_id) !== nodeId) continue;
          const x = xs[i], y = ys[i];
          if (x === null || x === undefined || y === null || y === undefined) continue;
          const x0 = (cd.x0 !== undefined && cd.x0 !== null) ? cd.x0 : x;
          const x1 = (cd.x1 !== undefined && cd.x1 !== null) ? cd.x1 : x;
          out.push({{ x0: Number(x0), x1: Number(x1), y: Number(y) }});
        }}
      }}
      return out;
    }}

    function zoomToNodeId() {{
      const gd = document.getElementById("{div_id}");
      const input = document.getElementById("node-id-input");
      if (!gd || !input) return;
      const nodeId = (input.value || "").trim();
      if (!nodeId) {{ setStatus("Enter a node_id."); return; }}
      const matches = collectNodeMatches(gd, nodeId);
      if (!matches.length) {{ setStatus("No cluster found for node_id=" + nodeId); return; }}
      let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
      for (const m of matches) {{
        if (!Number.isFinite(m.x0) || !Number.isFinite(m.x1) || !Number.isFinite(m.y)) continue;
        minX = Math.min(minX, m.x0, m.x1); maxX = Math.max(maxX, m.x0, m.x1);
        minY = Math.min(minY, m.y); maxY = Math.max(maxY, m.y);
      }}
      if (!Number.isFinite(minX) || !Number.isFinite(maxX) || !Number.isFinite(minY) || !Number.isFinite(maxY)) {{
        setStatus("Unable to zoom for node_id=" + nodeId); return;
      }}
      const xPad = Math.max((maxX - minX) * 0.1, 10000);
      const yPad = Math.max((maxY - minY) * 0.2, 0.5);
      Plotly.relayout(gd, {{
        "xaxis.range": [minX - xPad, maxX + xPad],
        "yaxis.range": [Math.max(0, minY - yPad), maxY + yPad]
      }});
      setStatus("Zoomed to node_id=" + nodeId + " (" + matches.length + " matches).");
    }}

    function resetAxes() {{
      const gd = document.getElementById("{div_id}");
      if (!gd) return;
      Plotly.relayout(gd, {{"xaxis.autorange": true, "yaxis.autorange": true}});
      setStatus("Reset to full view.");
    }}

    function zoomToRegion() {{
      const gd = document.getElementById("{div_id}");
      const input = document.getElementById("region-input");
      if (!gd || !input) return;
      const parsed = parseRegionSpec(input.value || "");
      if (!parsed) {{ setStatus("Use region format: chrX:100000-200000"); return; }}
      const lengths = (regionLookup && regionLookup.lengths) ? regionLookup.lengths : {{}};
      const offsets = (regionLookup && regionLookup.offsets) ? regionLookup.offsets : {{}};
      if (!(parsed.chrom in lengths) || !(parsed.chrom in offsets)) {{
        setStatus("Chromosome not found in plot: " + parsed.chrom); return;
      }}
      const chromLen = Number(lengths[parsed.chrom]);
      const off = Number(offsets[parsed.chrom]);
      if (!Number.isFinite(chromLen) || !Number.isFinite(off)) {{
        setStatus("Chromosome metadata unavailable for " + parsed.chrom); return;
      }}
      const startLocal = Math.max(0, Math.min(parsed.start, chromLen));
      const endLocal = Math.max(0, Math.min(parsed.end, chromLen));
      const minX = off + startLocal, maxX = off + endLocal;
      const xPad = Math.max((maxX - minX) * 0.05, 1000);
      Plotly.relayout(gd, {{"xaxis.range": [minX - xPad, maxX + xPad], "yaxis.autorange": true}});
      setStatus("Zoomed to " + parsed.chrom + ":" + startLocal + "-" + endLocal);
    }}

    (function attachHandlers() {{
      const gd = document.getElementById("{div_id}");
      if (!gd) return;
      gd.on("plotly_click", async function(ev) {{
        if (!ev || !Array.isArray(ev.points) || ev.points.length === 0) return;
        const p = ev.points[0];
        const data = p && p.data ? p.data : null;
        const cd = p ? p.customdata : null;
        if (!data || typeof data.name !== "string" || !data.name.startsWith("Timing chains")) return;
        if (!cd || cd.node_ids_text === undefined || cd.node_ids_text === null) return;
        if (!cd.file_rel) {{ setStatus("Selected chain is missing file metadata."); return; }}
        selectedChain = {{
          file_rel: String(cd.file_rel),
          chain_index: (cd.chain_index !== undefined && cd.chain_index !== null) ? Number(cd.chain_index) : null,
          node_ids_text: String(cd.node_ids_text || "")
        }};
        const chainLabel = (selectedChain.chain_index === null || Number.isNaN(selectedChain.chain_index))
          ? selectedChain.file_rel : (selectedChain.file_rel + " | chain_index=" + selectedChain.chain_index);
        setSelectedChainStatus("Selected chain: " + chainLabel);
        const nodeIdsText = String(cd.node_ids_text).trim();
        if (!nodeIdsText) return;
        const copied = await copyTextToClipboard(nodeIdsText);
        if (copied) setStatus("Chain selected. Copied condensed node IDs for selected chain.");
        else setStatus("Chain selected. Could not copy node IDs to clipboard.");
      }});
      gd.on("plotly_legendclick", function(ev) {{
        const tr = gd.data && gd.data[ev.curveNumber];
        if (tr && tr.meta && tr.meta.trace_role === "marker_toggle_legend") {{
          markersVisible = !markersVisible;
          setMarkersVisible(gd, markersVisible);
          setStatus(markersVisible ? "Timing markers shown." : "Timing markers hidden.");
          return false;
        }}
        if (tr && tr.meta && tr.meta.trace_role === "chain_cn_violation") {{
          violationsVisible = !violationsVisible;
          setViolationsVisible(gd, violationsVisible);
          setStatus(violationsVisible ? "Chain>maxCN markers shown." : "Chain>maxCN markers hidden.");
          return false;
        }}
      }});
      const nodeInput = document.getElementById("node-id-input");
      if (nodeInput) nodeInput.addEventListener("keydown", ev => {{ if (ev.key === "Enter") {{ ev.preventDefault(); zoomToNodeId(); }} }});
      const regionInput = document.getElementById("region-input");
      if (regionInput) regionInput.addEventListener("keydown", ev => {{ if (ev.key === "Enter") {{ ev.preventDefault(); zoomToRegion(); }} }});
      const chainModeSelect = document.getElementById("chain-mode-select");
      if (chainModeSelect) chainModeSelect.value = "all";
      setChainMode("all");
      const distCnSelect = document.getElementById("dist-cn-mode-select");
      if (distCnSelect) distCnSelect.value = hasCn ? "total" : "none";
      const distXAxisModeSelect = document.getElementById("dist-xaxis-mode-select");
      if (distXAxisModeSelect) distXAxisModeSelect.value = hasCn ? "cn" : "mode";
      const distChainSelect = document.getElementById("dist-chain-mode-select");
      if (distChainSelect) distChainSelect.value = "all";
      const distCountSelect = document.getElementById("dist-count-mode-select");
      if (distCountSelect) distCountSelect.value = "chains";
      const distRankClusterMinMutsInput = document.getElementById("dist-rank-cluster-min-muts-input");
      if (distRankClusterMinMutsInput) distRankClusterMinMutsInput.value = "1";
      const distCountNormSelect = document.getElementById("dist-count-normalize-select");
      if (distCountNormSelect) distCountNormSelect.value = "raw";
      setDistributionMode();
    }})();
  </script>
</body>
</html>
"""
    out_html.write_text(page, encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Create interactive genome-wide timing plot (with optional CN tracks)."
    )
    ap.add_argument("--outdir", required=True, type=Path, help="Main output directory.")
    ap.add_argument("--cn-bed-hp1", type=Path, default=None, help="HP1 CN BED (optional).")
    ap.add_argument("--cn-bed-hp2", type=Path, default=None, help="HP2 CN BED (optional).")
    ap.add_argument("--timing-tsv", type=Path, default=None, help="Timing chains TSV (default: <outdir>/timing_chains.tsv).")
    ap.add_argument("--min-chain-len", type=int, default=3, help="Minimum chain length to include.")
    ap.add_argument("--title", type=str, default="Genome-wide copy number and timing chains", help="Plot title.")
    args = ap.parse_args()

    outdir = args.outdir.resolve()
    if (args.cn_bed_hp1 is None) != (args.cn_bed_hp2 is None):
        raise SystemExit("Provide both --cn-bed-hp1 and --cn-bed-hp2 together, or neither.")
    has_cn = args.cn_bed_hp1 is not None and args.cn_bed_hp2 is not None
    timing_tsv = args.timing_tsv.resolve() if args.timing_tsv else (outdir / "timing_chains.tsv")
    ranks_tsv = outdir / "node_ranks.tsv"
    out_html = outdir / "timing_cn_interactive.html"

    if has_cn:
        cn_hp1 = load_cn_segments(args.cn_bed_hp1.resolve())
        cn_hp2 = load_cn_segments(args.cn_bed_hp2.resolve())
        cn_total = build_total_cn_segments(cn_hp1, cn_hp2)
        cn_max = build_max_cn_segments(cn_hp1, cn_hp2)
    else:
        cn_hp1 = {}
        cn_hp2 = {}
        cn_total = {}
        cn_max = {}
    timing_chains = load_timing_chains(outdir, timing_tsv, min_chain_len=max(1, args.min_chain_len))
    ranked_nodes = load_node_ranks(outdir, ranks_tsv)
    mutation_vaf_points = load_uncondensed_mutation_vaf_points(outdir)

    fig = build_figure(
        cn_hp1=cn_hp1,
        cn_hp2=cn_hp2,
        cn_total=cn_total,
        cn_max_haplotype=cn_max,
        timing_chains=timing_chains,
        title=args.title,
        ranked_nodes=ranked_nodes,
        show_cn=has_cn,
    )
    if has_cn:
        fig_dist = build_chain_length_cn_distribution_figure(
            cn_hp1=cn_hp1,
            cn_hp2=cn_hp2,
            cn_total=cn_total,
            cn_max_haplotype=cn_max,
            timing_chains_all=timing_chains,
            timing_chains_non_overlap=select_non_overlapping_chains(timing_chains),
            timing_chains_longest_only=select_longest_chain_per_cluster(timing_chains),
            ranked_nodes=ranked_nodes,
            mutation_vaf_points=mutation_vaf_points,
        )
    else:
        fig_dist = build_chain_length_distribution_no_cn_figure(
            timing_chains_all=timing_chains,
            timing_chains_non_overlap=select_non_overlapping_chains(timing_chains),
            timing_chains_longest_only=select_longest_chain_per_cluster(timing_chains),
            ranked_nodes=ranked_nodes,
            mutation_vaf_points=mutation_vaf_points,
        )
    chrom_lengths = get_chrom_lengths(cn_hp1, cn_hp2, cn_total, timing_chains, ranked_nodes)
    offsets = chrom_offsets(chrom_lengths)
    out_html.parent.mkdir(parents=True, exist_ok=True)
    write_interactive_html(fig, fig_dist, out_html, outdir, chrom_lengths, offsets, has_cn=has_cn)
    print(f"Wrote interactive plot: {out_html}")
    print(f"Timing chains plotted: {len(timing_chains)}")
    print(f"Ranked nodes plotted: {len(ranked_nodes)}")
    print(f"Mutation VAF points plotted: {len(mutation_vaf_points)}")


if __name__ == "__main__":
    main()

