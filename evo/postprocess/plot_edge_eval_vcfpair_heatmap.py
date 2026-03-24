#!/usr/bin/env python3
"""
Create interactive VCF-label-pair heatmap from edge evaluation counts.

Input:
  - edge_eval_vcfpair_relation_counts.tsv

Output:
  - edge_eval_vcfpair_heatmap.html
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict, deque
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def label_sort_key(label: str) -> Tuple[int, int, str]:
    s = str(label).strip()
    if s.upper() == "UNKNOWN":
        return (2, 10**9, s)
    if len(s) >= 2 and s[0].upper() == "N" and s[1:].isdigit():
        return (0, int(s[1:]), "")
    return (1, 10**9, s)


def load_rows(path: Path) -> List[Dict[str, int | str]]:
    out: List[Dict[str, int | str]] = []
    with path.open("rt", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {
            "vcf_label_a",
            "vcf_label_b",
            "known_total",
            "unknown",
            "consistent_timing",
            "consistent_cooccurrence",
            "consistent_divergent",
            "inconsistent_timing",
            "inconsistent_cooccurrence",
            "inconsistent_divergent",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} missing required columns: {', '.join(sorted(missing))}")

        for r in reader:
            row = {
                "vcf_label_a": str(r["vcf_label_a"]).strip() or "UNKNOWN",
                "vcf_label_b": str(r["vcf_label_b"]).strip() or "UNKNOWN",
            }
            for k in (
                "known_total",
                "unknown",
                "consistent_timing",
                "consistent_cooccurrence",
                "consistent_divergent",
                "inconsistent_timing",
                "inconsistent_cooccurrence",
                "inconsistent_divergent",
            ):
                try:
                    row[k] = int(str(r.get(k, "0")).strip() or "0")
                except ValueError:
                    row[k] = 0
            out.append(row)
    return out


def load_tree_edges(tree_path: Optional[Path]) -> List[Tuple[str, str]]:
    if tree_path is None or not tree_path.exists():
        return []
    edges: List[Tuple[str, str]] = []
    with tree_path.open("rt", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split("\t")
            if len(parts) != 2:
                continue
            parent = parts[0].strip()
            child = parts[1].strip()
            if not parent or not child:
                continue
            edges.append((parent, child))
    return edges


def build_tree_payload(edges: List[Tuple[str, str]]) -> Dict[str, object]:
    if not edges:
        return {"nodes": [], "edges": []}

    children: Dict[str, List[str]] = defaultdict(list)
    indeg: Dict[str, int] = defaultdict(int)
    nodes: set[str] = set()
    for p, c in edges:
        children[p].append(c)
        indeg[c] += 1
        nodes.add(p)
        nodes.add(c)
    for n in nodes:
        indeg.setdefault(n, 0)
        children.setdefault(n, [])
    for n in children:
        children[n].sort(key=label_sort_key)

    roots = sorted([n for n in nodes if indeg[n] == 0], key=label_sort_key)
    if not roots:
        roots = sorted(nodes, key=label_sort_key)

    depth: Dict[str, int] = {}
    q: deque[str] = deque()
    for r in roots:
        depth[r] = 0
        q.append(r)
    while q:
        cur = q.popleft()
        for ch in children.get(cur, []):
            cand = depth[cur] + 1
            if ch not in depth or cand < depth[ch]:
                depth[ch] = cand
                q.append(ch)
    for n in nodes:
        depth.setdefault(n, 0)

    by_depth: Dict[int, List[str]] = defaultdict(list)
    for n in nodes:
        by_depth[depth[n]].append(n)
    for d in by_depth:
        by_depth[d].sort(key=label_sort_key)

    node_rows: List[Dict[str, object]] = []
    for d in sorted(by_depth):
        row = by_depth[d]
        center = (len(row) - 1) / 2.0
        for i, nid in enumerate(row):
            node_rows.append({"id": nid, "x": float(i - center), "y": float(-d), "depth": d})
    return {"nodes": node_rows, "edges": [[p, c] for p, c in edges]}


def write_html(
    rows: List[Dict[str, int | str]],
    tree_payload: Dict[str, object],
    out_html: Path,
    title: str,
) -> None:
    labels = sorted(
        {str(r["vcf_label_a"]) for r in rows} | {str(r["vcf_label_b"]) for r in rows},
        key=label_sort_key,
    )
    payload = {"rows": rows, "labels": labels, "tree": tree_payload}
    payload_json = json.dumps(payload)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{title}</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    body {{ margin: 0; font-family: Arial, sans-serif; background: #fff; }}
    .controls {{
      display: flex; gap: 8px; align-items: center; flex-wrap: wrap;
      padding: 10px 12px; border-bottom: 1px solid #ddd; background: #fafafa;
    }}
    .controls select {{ padding: 6px 8px; font-size: 14px; }}
    #status {{ font-size: 13px; color: #444; }}
    .content {{
      display: flex;
      gap: 12px;
      padding: 10px 12px;
      box-sizing: border-box;
      align-items: flex-start;
    }}
    .left-panel {{
      flex: 1 1 auto;
      min-width: 0;
    }}
    .right-panel {{
      flex: 0 0 390px;
      border-left: 1px solid #e3e3e3;
      padding-left: 12px;
      display: flex;
      flex-direction: column;
      min-width: 320px;
      max-width: 450px;
    }}
    #plot {{
      width: 100%;
      height: 72vh;
      min-height: 460px;
    }}
    #selected-info {{
      margin-top: 8px;
      padding: 4px 2px 0 2px;
      font-size: 15px;
      line-height: 1.4;
      color: #222;
    }}
    #selected-info .title {{
      font-size: 16px;
      font-weight: 700;
      margin-bottom: 6px;
    }}
    #selected-info .muted {{
      color: #666;
      font-size: 14px;
    }}
    #selected-info table {{
      border-collapse: collapse;
      margin-top: 6px;
      font-size: 14px;
      min-width: 360px;
    }}
    #selected-info th, #selected-info td {{
      border: 1px solid #d9d9d9;
      padding: 5px 8px;
      text-align: right;
    }}
    #selected-info th:first-child, #selected-info td:first-child {{
      text-align: left;
      font-weight: 600;
    }}
    #selected-info th.consistent {{
      color: #2ca02c;
      font-weight: 700;
    }}
    #selected-info th.inconsistent {{
      color: #d62728;
      font-weight: 700;
    }}
    #tree-title {{
      font-size: 16px;
      font-weight: 700;
      margin: 4px 0 6px 0;
    }}
    #tree-controls {{
      margin: 0 0 8px 0;
      display: flex;
      gap: 8px;
      align-items: center;
      font-size: 14px;
    }}
    #tree-controls select {{
      padding: 4px 6px;
      font-size: 14px;
    }}
    #tree-plot {{
      width: 100%;
      height: 72vh;
      min-height: 420px;
      border: 1px solid #ececec;
      border-radius: 8px;
    }}
  </style>
</head>
<body>
  <div class="controls">
    <label for="status-mode">status:</label>
    <select id="status-mode" onchange="renderHeatmap()">
      <option value="consistent" selected>consistent</option>
      <option value="inconsistent">inconsistent</option>
      <option value="unknown">unknown</option>
    </select>

    <label for="relation-mode">relation:</label>
    <select id="relation-mode" onchange="renderHeatmap()">
      <option value="all">all</option>
      <option value="timing" selected>timing</option>
      <option value="cooccurrence">cooccurrence</option>
      <option value="divergent">divergent</option>
    </select>

    <label for="value-mode">% known option:</label>
    <select id="value-mode" onchange="renderHeatmap()">
      <option value="count" selected>count</option>
      <option value="pct_known">% of known relations</option>
    </select>

    <span id="status"></span>
  </div>
  <div class="content">
    <div class="left-panel">
      <div id="plot"></div>
      <div id="selected-info">
        <div class="title">Selected pair details</div>
        <div class="muted">Click a heatmap cell to show details here.</div>
      </div>
    </div>
    <div class="right-panel">
      <div id="tree-title">Tree (from --tree)</div>
      <div id="tree-controls">
        <label for="tree-view-mode">view:</label>
        <select id="tree-view-mode" onchange="renderTree()">
          <option value="vcf_only" selected>VCF nodes only</option>
          <option value="full">Full tree</option>
        </select>
      </div>
      <div id="tree-plot"></div>
    </div>
  </div>

  <script>
    const payload = {payload_json};
    const rows = Array.isArray(payload.rows) ? payload.rows : [];
    const labels = Array.isArray(payload.labels) ? payload.labels : [];
    const treePayload = payload.tree || {{ nodes: [], edges: [] }};
    const statusEl = document.getElementById("status");
    const selectedInfoEl = document.getElementById("selected-info");

    function getPairMap() {{
      const out = new Map();
      for (const r of rows) {{
        const a = String(r.vcf_label_a || "UNKNOWN");
        const b = String(r.vcf_label_b || "UNKNOWN");
        out.set(a + "|" + b, r);
      }}
      return out;
    }}

    function rowOrZero(r) {{
      return {{
        known_total: r ? Number(r.known_total || 0) : 0,
        unknown: r ? Number(r.unknown || 0) : 0,
        consistent_timing: r ? Number(r.consistent_timing || 0) : 0,
        consistent_cooccurrence: r ? Number(r.consistent_cooccurrence || 0) : 0,
        consistent_divergent: r ? Number(r.consistent_divergent || 0) : 0,
        inconsistent_timing: r ? Number(r.inconsistent_timing || 0) : 0,
        inconsistent_cooccurrence: r ? Number(r.inconsistent_cooccurrence || 0) : 0,
        inconsistent_divergent: r ? Number(r.inconsistent_divergent || 0) : 0
      }};
    }}

    function mergeRows(r1, r2) {{
      return {{
        known_total: Number(r1.known_total || 0) + Number(r2.known_total || 0),
        unknown: Number(r1.unknown || 0) + Number(r2.unknown || 0),
        consistent_timing: Number(r1.consistent_timing || 0) + Number(r2.consistent_timing || 0),
        consistent_cooccurrence: Number(r1.consistent_cooccurrence || 0) + Number(r2.consistent_cooccurrence || 0),
        consistent_divergent: Number(r1.consistent_divergent || 0) + Number(r2.consistent_divergent || 0),
        inconsistent_timing: Number(r1.inconsistent_timing || 0) + Number(r2.inconsistent_timing || 0),
        inconsistent_cooccurrence: Number(r1.inconsistent_cooccurrence || 0) + Number(r2.inconsistent_cooccurrence || 0),
        inconsistent_divergent: Number(r1.inconsistent_divergent || 0) + Number(r2.inconsistent_divergent || 0)
      }};
    }}

    function getCountFromRows(rDirected, rUnordered, statusMode, relationMode) {{
      if (statusMode === "unknown") return Number(rUnordered.unknown || 0);
      if (relationMode === "all") {{
        return (
          Number(rDirected[statusMode + "_timing"] || 0)
          + Number(rUnordered[statusMode + "_cooccurrence"] || 0)
          + Number(rUnordered[statusMode + "_divergent"] || 0)
        );
      }}
      if (relationMode === "timing") return Number(rDirected[statusMode + "_timing"] || 0);
      if (relationMode === "cooccurrence") return Number(rUnordered[statusMode + "_cooccurrence"] || 0);
      return Number(rUnordered[statusMode + "_divergent"] || 0);
    }}

    function renderHeatmap() {{
      const statusMode = String(document.getElementById("status-mode").value || "consistent");
      const relationMode = String(document.getElementById("relation-mode").value || "timing");
      const valueMode = String(document.getElementById("value-mode").value || "count");
      const relationSel = document.getElementById("relation-mode");
      relationSel.disabled = (statusMode === "unknown");

      const pairMap = getPairMap();
      const directionalMode = (relationMode === "timing");
      const z = [];
      const zText = [];
      const custom = [];
      for (const a of labels) {{
        const rowZ = [];
        const rowText = [];
        const rowCustom = [];
        for (const b of labels) {{
          const r_ab = rowOrZero(pairMap.get(String(a) + "|" + String(b)) || null);
          const r_ba = rowOrZero(pairMap.get(String(b) + "|" + String(a)) || null);
          const r_unordered = String(a) === String(b) ? r_ab : mergeRows(r_ab, r_ba);
          const known = Number(r_unordered.known_total || 0);
          const unknown = Number(r_unordered.unknown || 0);
          // Timing remains directional (a->b); others are symmetric (unordered).
          const c_t = Number(r_ab.consistent_timing || 0);
          const i_t = Number(r_ab.inconsistent_timing || 0);
          const c_c = Number(r_unordered.consistent_cooccurrence || 0);
          const c_d = Number(r_unordered.consistent_divergent || 0);
          const i_c = Number(r_unordered.inconsistent_cooccurrence || 0);
          const i_d = Number(r_unordered.inconsistent_divergent || 0);

          const rawCount = getCountFromRows(r_ab, r_unordered, statusMode, relationMode);
          let value = rawCount;
          if (valueMode === "pct_known") {{
            value = known > 0 ? (100.0 * rawCount / known) : 0.0;
          }}
          rowZ.push(value);
          rowText.push(
            Math.abs(Number(value) || 0) < 1e-12
              ? ""
              : (
                  valueMode === "pct_known"
                    ? (Number.isFinite(value) ? value.toFixed(1) + "%" : "NA")
                    : String(Math.trunc(rawCount))
                )
          );
          rowCustom.push([
            a, b, known, unknown,
            c_t, c_c, c_d,
            i_t, i_c, i_d,
            rawCount,
            (known > 0 ? (100.0 * rawCount / known) : NaN),
            directionalMode ? "timing directed (a->b); others unordered" : "unordered (a<->b)"
          ]);
        }}
        z.push(rowZ);
        zText.push(rowText);
        custom.push(rowCustom);
      }}

      const metricLabel = (statusMode === "unknown")
        ? "unknown"
        : (statusMode + " " + relationMode);
      const title = (valueMode === "pct_known")
        ? ("VCF-pair heatmap: " + metricLabel + " (% of known relations)")
        : ("VCF-pair heatmap: " + metricLabel + " (count)");

      const hovertemplate =
        "vcf pair: %{{customdata[0]}} - %{{customdata[1]}}<br>" +
        "known_total: %{{customdata[2]}}<br>" +
        "unknown: %{{customdata[3]}}<br>" +
        "consistent_timing: %{{customdata[4]}}<br>" +
        "consistent_cooccurrence: %{{customdata[5]}}<br>" +
        "consistent_divergent: %{{customdata[6]}}<br>" +
        "inconsistent_timing: %{{customdata[7]}}<br>" +
        "inconsistent_cooccurrence: %{{customdata[8]}}<br>" +
        "inconsistent_divergent: %{{customdata[9]}}<br>" +
        "selected_count: %{{customdata[10]}}<br>" +
        "selected_%_known: %{{customdata[11]:.2f}}%<br>" +
        "mode: %{{customdata[12]}}<extra></extra>";

      const trace = {{
        type: "heatmap",
        x: labels,
        y: labels,
        z: z,
        text: zText,
        texttemplate: "%{{text}}",
        textfont: {{ size: 11 }},
        customdata: custom,
        colorscale: "Viridis",
        colorbar: {{ title: valueMode === "pct_known" ? "% known" : "count" }},
        hovertemplate: hovertemplate
      }};

      const layout = {{
        template: "plotly_white",
        margin: {{ l: 90, r: 30, t: 40, b: 120 }},
        xaxis: {{ title: "VCF label", tickangle: -45, side: "top" }},
        yaxis: {{ title: "VCF label", autorange: "reversed" }},
        annotations: [{{
          text: title,
          xref: "paper",
          yref: "paper",
          x: 0.5,
          y: -0.22,
          showarrow: false,
          xanchor: "center",
          yanchor: "top",
          font: {{ size: 16 }},
        }}],
      }};

      Plotly.newPlot("plot", [trace], layout, {{responsive: true}});
      statusEl.textContent = "Hover a cell for full pair details.";
      const gd = document.getElementById("plot");
      gd.on("plotly_click", function(ev) {{
        if (!ev || !Array.isArray(ev.points) || ev.points.length === 0) return;
        const p = ev.points[0];
        const cd = p.customdata || [];
        const a = String(cd[0] || "");
        const b = String(cd[1] || "");
        const known = Number(cd[2] || 0);
        const unknown = Number(cd[3] || 0);
        const c_t = Number(cd[4] || 0);
        const c_c = Number(cd[5] || 0);
        const c_d = Number(cd[6] || 0);
        const i_t = Number(cd[7] || 0);
        const i_c = Number(cd[8] || 0);
        const i_d = Number(cd[9] || 0);
        const selectedCount = Number(cd[10] || 0);
        const selectedPct = Number(cd[11]);
        const pctText = Number.isFinite(selectedPct) ? (selectedPct.toFixed(2) + "%") : "N/A";
        const metricLabel = (statusMode === "unknown") ? "unknown" : (statusMode + " " + relationMode);
        selectedInfoEl.innerHTML =
          "<div class='title'>Selected pair: <b>" + a + "</b> -> <b>" + b + "</b></div>" +
          "<div><b>Selected metric</b>: " + metricLabel + " | " +
          (valueMode === "pct_known" ? "% of known relations" : "count") +
          "</div>" +
          "<div><b>Selected value</b>: " +
          (valueMode === "pct_known" ? pctText : selectedCount.toLocaleString()) +
          "</div>" +
          "<div style='margin-top:4px;'><b>Known total</b>: " + known.toLocaleString() +
          " &nbsp;|&nbsp; <b>Unknown</b>: " + unknown.toLocaleString() + "</div>" +
          "<table>" +
          "<thead><tr><th>Relation</th><th class='consistent'>Consistent</th><th class='inconsistent'>Inconsistent</th></tr></thead>" +
          "<tbody>" +
          "<tr><td>Timing</td><td>" + c_t.toLocaleString() + "</td><td>" + i_t.toLocaleString() + "</td></tr>" +
          "<tr><td>Cooccurrence</td><td>" + c_c.toLocaleString() + "</td><td>" + i_c.toLocaleString() + "</td></tr>" +
          "<tr><td>Divergent</td><td>" + c_d.toLocaleString() + "</td><td>" + i_d.toLocaleString() + "</td></tr>" +
          "</tbody></table>" +
          "<div class='muted' style='margin-top:6px;'>Detected-edge evaluation only.</div>";
      }});
    }}

    function renderTree() {{
      const treeDiv = document.getElementById("tree-plot");
      const treeViewSel = document.getElementById("tree-view-mode");
      const nodes = Array.isArray(treePayload.nodes) ? treePayload.nodes : [];
      const edges = Array.isArray(treePayload.edges) ? treePayload.edges : [];
      const viewMode = treeViewSel ? String(treeViewSel.value || "vcf_only") : "vcf_only";
      if (!treeDiv) return;
      if (!nodes.length || !edges.length) {{
        treeDiv.innerHTML = "<div style='padding:12px;color:#666;'>No tree data available.</div>";
        return;
      }}
      const axisLabels = new Set(labels.filter(v => String(v).toUpperCase() !== "UNKNOWN").map(v => String(v)));
      let activeNodes = nodes;
      let activeEdges = edges;
      if (viewMode === "vcf_only") {{
        activeNodes = nodes.filter(n => axisLabels.has(String(n.id)));
        const keep = new Set(activeNodes.map(n => String(n.id)));
        activeEdges = edges.filter(e => Array.isArray(e) && e.length === 2 && keep.has(String(e[0])) && keep.has(String(e[1])));
      }}
      if (!activeNodes.length) {{
        treeDiv.innerHTML = "<div style='padding:12px;color:#666;'>No tree nodes match VCF labels on the heatmap axes.</div>";
        return;
      }}
      const byId = new Map();
      for (const n of activeNodes) {{
        byId.set(String(n.id), n);
      }}
      const ex = [];
      const ey = [];
      for (const e of activeEdges) {{
        if (!Array.isArray(e) || e.length !== 2) continue;
        const p = byId.get(String(e[0]));
        const c = byId.get(String(e[1]));
        if (!p || !c) continue;
        ex.push(Number(p.x), Number(c.x), null);
        ey.push(Number(p.y), Number(c.y), null);
      }}
      const nx = [];
      const ny = [];
      const nt = [];
      for (const n of activeNodes) {{
        nx.push(Number(n.x));
        ny.push(Number(n.y));
        nt.push(String(n.id));
      }}
      const traces = [
        {{
          type: "scatter",
          mode: "lines",
          x: ex,
          y: ey,
          line: {{ width: 1.2, color: "#9aa0a6" }},
          hoverinfo: "skip",
          showlegend: false
        }},
        {{
          type: "scatter",
          mode: "markers+text",
          x: nx,
          y: ny,
          text: nt,
          textposition: "top center",
          marker: {{ size: 10, color: "#1f77b4" }},
          hovertemplate: "node=%{{text}}<extra></extra>",
          showlegend: false
        }}
      ];
      const layout = {{
        template: "plotly_white",
        margin: {{ l: 30, r: 20, t: 20, b: 20 }},
        xaxis: {{ visible: false, showgrid: false, zeroline: false }},
        yaxis: {{ visible: false, showgrid: false, zeroline: false }},
      }};
      Plotly.newPlot(treeDiv, traces, layout, {{ responsive: true, displayModeBar: false, staticPlot: true }});
    }}

    renderHeatmap();
    renderTree();
  </script>
</body>
</html>
"""
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html, encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description="Plot VCF-pair edge evaluation heatmap.")
    ap.add_argument(
        "--counts-tsv",
        type=Path,
        default=None,
        help="Path to edge_eval_vcfpair_relation_counts.tsv (default: <outdir>/edge_eval_vcfpair_relation_counts.tsv).",
    )
    ap.add_argument("--outdir", type=Path, required=True, help="Main output directory.")
    ap.add_argument("--tree", type=Path, default=None, help="Parent-child TSV used for static side tree panel.")
    ap.add_argument(
        "--output-name",
        type=str,
        default="edge_eval_vcfpair_heatmap.html",
        help="Output HTML filename under --outdir.",
    )
    ap.add_argument(
        "--title",
        type=str,
        default="VCF-pair edge evaluation heatmap",
        help="HTML page title.",
    )
    args = ap.parse_args()

    outdir = args.outdir.resolve()
    counts_tsv = args.counts_tsv.resolve() if args.counts_tsv else (outdir / "edge_eval_vcfpair_relation_counts.tsv")
    if not counts_tsv.exists():
        raise FileNotFoundError(f"Counts TSV not found: {counts_tsv}")
    rows = load_rows(counts_tsv)
    tree_edges = load_tree_edges(args.tree.resolve() if args.tree else None)
    tree_payload = build_tree_payload(tree_edges)
    out_html = outdir / args.output_name
    write_html(rows, tree_payload, out_html, args.title)
    print(f"Wrote: {out_html}")


if __name__ == "__main__":
    main()

