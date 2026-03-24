#!/usr/bin/env python3
"""
Bridge service from interactive plot -> Cytoscape (cyREST).

The bridge expects outdir in POST payload from the interactive HTML.
Endpoint:
  POST /open_chain
Body:
  {
    "outdir": "/path/to/main_out",
    "file_rel": "chr1_out/graphs/chunk_chr1_..._condensed.json",
    "open_condensed": true,
    "open_uncondensed": true
  }
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
import urllib.error
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Optional


def _json_response(handler: BaseHTTPRequestHandler, code: int, payload: dict[str, Any]) -> None:
    body = json.dumps(payload).encode("utf-8")
    handler.send_response(code)
    handler.send_header("Content-Type", "application/json; charset=utf-8")
    handler.send_header("Content-Length", str(len(body)))
    handler.send_header("Access-Control-Allow-Origin", "*")
    handler.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
    handler.send_header("Access-Control-Allow-Headers", "Content-Type")
    handler.end_headers()
    handler.wfile.write(body)


def cyrest_get_version(base: str, timeout_s: float = 2.0) -> Optional[dict[str, Any]]:
    url = base.rstrip("/") + "/v1/version"
    req = urllib.request.Request(url, method="GET")
    try:
        with urllib.request.urlopen(req, timeout=timeout_s) as resp:
            return json.loads(resp.read().decode("utf-8", errors="replace"))
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, ValueError):
        return None


def maybe_launch_cytoscape() -> None:
    if sys.platform == "darwin":
        subprocess.Popen(["open", "-a", "Cytoscape"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    elif sys.platform.startswith("linux"):
        subprocess.Popen(["cytoscape"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def wait_for_cyrest(base: str, timeout_s: float = 30.0) -> bool:
    deadline = time.time() + timeout_s
    while time.time() < deadline:
        if cyrest_get_version(base) is not None:
            return True
        time.sleep(0.75)
    return False


def import_network(base: str, graph_path: Path) -> dict[str, Any]:
    url = base.rstrip("/") + "/v1/networks?format=cyjs"
    req = urllib.request.Request(
        url,
        data=graph_path.read_bytes(),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=20.0) as resp:
        raw = resp.read().decode("utf-8", errors="replace")
    try:
        return json.loads(raw) if raw else {}
    except json.JSONDecodeError:
        return {"raw_response": raw}


def _resolve_paths(outdir: Path, file_rel: str, open_condensed: bool, open_uncondensed: bool) -> list[Path]:
    condensed = (outdir / file_rel).resolve()
    targets: list[Path] = []
    if open_condensed:
        targets.append(condensed)
    if open_uncondensed and condensed.name.endswith("_condensed.json"):
        uncondensed = condensed.with_name(condensed.name[: -len("_condensed.json")] + ".json")
        targets.append(uncondensed)
    return targets


def _candidate_outdirs(outdir_raw: str) -> list[Path]:
    """
    Build candidate local paths for a possibly remote path.
    Examples:
      /vf/users/KolmogorovLab/... -> /Volumes/KolmogorovLab/...
      /data/KolmogorovLab/...     -> /Volumes/KolmogorovLab/...
    """
    raw = str(outdir_raw).strip()
    candidates: list[Path] = []
    seen: set[str] = set()

    def add(path_str: str) -> None:
        p = str(Path(path_str).resolve())
        if p not in seen:
            seen.add(p)
            candidates.append(Path(p))

    add(raw)
    if raw.startswith("/vf/users/"):
        add("/Volumes/" + raw[len("/vf/users/") :])
    if raw.startswith("/data/"):
        add("/Volumes/" + raw[len("/data/") :])
    return candidates


def make_handler(cyrest_base: str, launch_if_needed: bool):
    class Handler(BaseHTTPRequestHandler):
        def do_OPTIONS(self) -> None:  # noqa: N802
            _json_response(self, 200, {"ok": True})

        def do_GET(self) -> None:  # noqa: N802
            if self.path.startswith("/health"):
                _json_response(self, 200, {"ok": True})
                return
            if self.path.startswith("/cyrest"):
                ver = cyrest_get_version(cyrest_base)
                _json_response(self, 200 if ver else 503, {"cyrest": ver})
                return
            _json_response(self, 404, {"error": "not_found"})

        def do_POST(self) -> None:  # noqa: N802
            if not self.path.startswith("/open_chain"):
                _json_response(self, 404, {"error": "not_found"})
                return

            try:
                n = int(self.headers.get("Content-Length", "0"))
            except ValueError:
                _json_response(self, 400, {"error": "invalid Content-Length"})
                return

            raw = self.rfile.read(n) if n > 0 else b"{}"
            try:
                payload = json.loads(raw.decode("utf-8", errors="replace"))
            except json.JSONDecodeError:
                _json_response(self, 400, {"error": "invalid JSON body"})
                return

            file_rel = str(payload.get("file_rel") or "").strip()
            if not file_rel:
                _json_response(self, 400, {"error": "missing file_rel"})
                return

            outdir_raw = payload.get("outdir")
            if not outdir_raw:
                _json_response(
                    self,
                    400,
                    {"error": "missing outdir in request payload (use the generated interactive HTML which includes it)"},
                )
                return
            candidates = _candidate_outdirs(str(outdir_raw))
            outdir = next((p for p in candidates if p.exists()), None)
            if outdir is None:
                _json_response(
                    self,
                    404,
                    {
                        "error": f"outdir does not exist: {outdir_raw}",
                        "tried": [str(p) for p in candidates],
                    },
                )
                return

            open_condensed = bool(payload.get("open_condensed", True))
            open_uncondensed = bool(payload.get("open_uncondensed", True))
            targets = _resolve_paths(outdir, file_rel, open_condensed, open_uncondensed)
            missing = [str(p) for p in targets if not p.exists()]
            if missing:
                _json_response(self, 404, {"error": "graph file(s) not found", "missing": missing})
                return

            if cyrest_get_version(cyrest_base) is None and launch_if_needed:
                maybe_launch_cytoscape()
                wait_for_cyrest(cyrest_base, timeout_s=30.0)

            ver = cyrest_get_version(cyrest_base)
            if ver is None:
                _json_response(
                    self,
                    503,
                    {"error": f"cyREST unavailable at {cyrest_base}. Start Cytoscape and enable cyREST."},
                )
                return

            opened: list[dict[str, Any]] = []
            for p in targets:
                try:
                    resp = import_network(cyrest_base, p)
                except urllib.error.HTTPError as e:
                    _json_response(self, 502, {"error": f"cyREST HTTP error for {p}: {e.code} {e.reason}"})
                    return
                except (urllib.error.URLError, OSError, ValueError, json.JSONDecodeError) as e:
                    _json_response(self, 500, {"error": f"failed opening {p}: {e}"})
                    return
                opened.append({"path": str(p), "cyrest_response": resp})

            _json_response(self, 200, {"ok": True, "opened": opened})

        def log_message(self, _fmt: str, *_args: Any) -> None:
            return

    return Handler


def main() -> None:
    ap = argparse.ArgumentParser(description="Open selected chain graph files in Cytoscape via cyREST.")
    ap.add_argument("--host", default="127.0.0.1", help="Bridge host (default: 127.0.0.1)")
    ap.add_argument("--port", type=int, default=8765, help="Bridge port (default: 8765)")
    ap.add_argument("--cyrest-base", default="http://127.0.0.1:1234", help="cyREST base URL.")
    ap.add_argument("--no-launch-cytoscape", action="store_true", help="Do not auto-launch Cytoscape.")
    args = ap.parse_args()

    handler = make_handler(args.cyrest_base, launch_if_needed=(not args.no_launch_cytoscape))
    server = ThreadingHTTPServer((args.host, args.port), handler)
    print(f"Cytoscape bridge listening on http://{args.host}:{args.port}")
    print(f"cyREST base: {args.cyrest_base}")
    print("outdir source: request payload from interactive HTML")
    print("health: GET /health")
    print("cyrest: GET /cyrest")
    print("open:   POST /open_chain")
    server.serve_forever()


if __name__ == "__main__":
    main()

