#!/usr/bin/env python3
"""dashboard.py — a local web dashboard for runs in progress: what is running, how far along,
how long is left, and every output downloadable the moment it exists.

    python3 scripts/dashboard.py                  # then open http://127.0.0.1:8765
    python3 scripts/dashboard.py --port 9000 --root ~/hmm_runs --root /mnt/c/Users/me/Downloads

WHY STDLIB ONLY. The pipeline already asks a lot of a machine (conda, HMMER, IQ-TREE, a 100 GB
cache). A monitoring UI that needs its own web framework is one more thing to install and break
before you can see whether your run is alive. This is `http.server` and nothing else, so it
works anywhere the pipeline does.

WHAT IT READS. Nothing bespoke — the same artefacts the pipeline already writes:
  _progress.json          streaming catalogue scans (contigs done, hits, batches)
  pipeline.log            the discovery run's own log
  collection_hits.tsv     batched collection scans
  *_summary.json/.csv     stage, coverage and scan summaries
  PACKAGE/                the deliverable
So a run started from the terminal appears here with no cooperation from the run itself, and
the dashboard can be started, stopped or restarted at any time without touching it.

SAFETY. Binds 127.0.0.1 only, serves GET only, and every served path must resolve inside one of
the --root directories — so a crafted path cannot walk out of them.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import mimetypes
import os
import re
import shutil
import socketserver
import subprocess
import sys
import threading
import time
import urllib.parse
from datetime import datetime, timedelta
from http.server import SimpleHTTPRequestHandler
from pathlib import Path

ROOTS: list = []
PROC_PAT = ("hmm_finder.py", "run_all_database_benchmark", "stream_scan_catalogue",
            "scan_genome_collection", "scan_host_genera", "scan_missed_seeds",
            "scan_full_coverage", "scan_genome.py", "preload_databases")

# Expected totals for the two catalogues, so a percentage and an ETA can be shown rather than
# a bare count that says nothing about how much is left.
# gpd is exact (a completed scan counted 142,809). gvd is approximate — the AVrC catalogue
# holds roughly 447k contigs, not the 305k first guessed, which made the percentage and the ETA
# read optimistically. An estimate that is too small is worse than none, because it turns into
# a confident "nearly finished".
CATALOGUE_TOTALS = {"gpd": 142809, "gvd": 447000}
APPROX_TOTALS = {"gvd"}          # shown with a ~ so nobody reads the % as exact


def _fmt_dur(sec: float) -> str:
    if sec is None or sec < 0 or sec != sec:
        return "—"
    sec = int(sec)
    if sec < 60:
        return f"{sec}s"
    if sec < 3600:
        return f"{sec // 60}m {sec % 60:02d}s"
    return f"{sec // 3600}h {(sec % 3600) // 60:02d}m"


def _fmt_size(n: int) -> str:
    x = float(n)
    for u in ("B", "KB", "MB", "GB", "TB"):
        if x < 1024 or u == "TB":
            return f"{x:.0f} {u}" if u == "B" else f"{x:.1f} {u}"
        x /= 1024
    return f"{x:.1f} TB"


def running_processes() -> list:
    """Pipeline processes currently alive (empty if ps is unavailable)."""
    out = []
    try:
        ps = subprocess.run(["ps", "-eo", "pid,etimes,args"], capture_output=True,
                            text=True, timeout=10)
        for ln in ps.stdout.splitlines()[1:]:
            parts = ln.strip().split(None, 2)
            if len(parts) < 3:
                continue
            pid, et, cmd = parts
            if "ps -eo" in cmd or "dashboard.py" in cmd:
                continue
            if any(p in cmd for p in PROC_PAT):
                name = next((p for p in PROC_PAT if p in cmd), "pipeline")
                out.append({"pid": pid, "elapsed": int(et) if et.isdigit() else 0,
                            "what": name, "cmd": cmd[:220]})
    except Exception:
        pass
    return out


def _scan_progress(d: Path) -> dict | None:
    """Progress of a streaming catalogue scan, with a rate-based ETA."""
    p = d / "_progress.json"
    if not p.exists():
        return None
    try:
        st = json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return None
    done = int(st.get("contigs_done", 0))
    key = "gpd" if "gpd" in d.name.lower() else ("gvd" if "gvd" in d.name.lower() else "")
    total = CATALOGUE_TOTALS.get(key, 0)
    age = time.time() - p.stat().st_mtime
    # rate from the file's own lifetime; good enough for an ETA and needs no extra state
    try:
        started = min(f.stat().st_mtime for f in d.iterdir() if f.is_file())
        el = max(1.0, p.stat().st_mtime - started)
    except Exception:
        el = 1.0
    rate = done / el if el else 0
    left = (total - done) / rate if (total and rate and total > done) else None
    return {"done": done, "total": total, "hits": int(st.get("hits", 0)),
            "batches": int(st.get("batches", 0)),
            "pct": round(100.0 * done / total, 1) if total else None,
            "approx": key in APPROX_TOTALS,
            # An end-of-scan summary is only written when a scan finishes, so its absence means
            # this scan did not complete — the same witness coverage_report now requires.
            "complete": bool(st.get("complete")) and (d / "stream_scan_summary.json").exists(),
            "eta": left, "stale": age > 900, "rate_per_min": round(rate * 60)}


def _tsv_rows(p: Path) -> int:
    try:
        return max(0, sum(1 for _ in open(p, encoding="utf-8", errors="replace")) - 1)
    except OSError:
        return 0


def discover() -> list:
    """Everything that looks like a run or a scan under the configured roots."""
    items = []
    for root in ROOTS:
        if not root.is_dir():
            continue
        try:
            entries = sorted(root.iterdir())
        except OSError:
            continue
        for d in entries:
            try:
                if not d.is_dir() or d.name.startswith("."):
                    continue
            except OSError:
                continue
            kind, detail = None, {}
            log = d / "pipeline.log"
            if log.exists() or (d / "PACKAGE").is_dir():
                kind = "discovery run"
                # Read only the tail: a completed run's log is megabytes, and re-reading all of
                # it on a 15-second refresh is pointless work.
                tail = ""
                try:
                    with log.open("rb") as fh:
                        fh.seek(max(0, log.stat().st_size - 8192))
                        tail = fh.read().decode("utf-8", "replace")
                except OSError:
                    pass
                detail["done"] = "=== DONE" in tail
                detail["tail"] = tail.splitlines()[-1][:160] if tail.strip() else ""
            else:
                try:
                    subnames = {x.name for x in d.iterdir() if x.is_dir()}
                except OSError:
                    subnames = set()
                if any(n.startswith("catalogue_") for n in subnames) or \
                        subnames & {"seed_sources", "host_genera", "results"}:
                    kind = "coverage scan"
                elif (d / "_progress.json").exists() or (d / "collection_hits.tsv").exists():
                    kind = "collection scan"
            if not kind:
                continue
            subs = []
            for sd in sorted([d] + [x for x in d.iterdir() if x.is_dir()]):
                pr = _scan_progress(sd)
                if pr:
                    subs.append({"name": sd.name if sd != d else ".", **pr})
            items.append({"path": str(d), "name": d.name, "kind": kind,
                          "mtime": d.stat().st_mtime, "progress": subs, **detail})
    items.sort(key=lambda x: -x["mtime"])
    return items


INTERESTING = (
    ("report.html", "the run report — open this first"),
    ("coverage_summary.csv", "what was and was not searched"),
    ("family_census.csv", "did the family grow"),
    ("paper_main_table.csv", "main homolog table"),
    ("hits_deduplicated.csv", "homologs, deduplicated"),
    ("overprinted_loci.csv", "overprinting, with host genes"),
    ("missed_seed_scan.csv", "per-seed verdicts"),
    ("collection_hits.tsv", "hits from this scan"),
    ("pipeline_stage_summary.csv", "every stage in one table"),
    ("stream_scan_summary.json", "catalogue scan summary"),
)


# Where an output can actually live, relative to a run directory. Checked EXPLICITLY rather
# than with rglob: a finished run holds hundreds of thousands of files (155 GenBanks, per-run
# sequence dumps, clinker HTML, a six-frame cache), and rglob over that — times ten filenames,
# times every run, on every page load — does not return. The layout is known, so look it up.
SUBDIRS = ("", "PACKAGE", "PACKAGE/01_summary_tables", "PACKAGE/09_controls",
           "PACKAGE/10_overprinting", "seed_qc", "results", "controls",
           "catalogue_gpd", "catalogue_gvd", "seed_sources", "seed_sources/results",
           "host_genera", "host_genera/results", "full_coverage", "downstream/overprinting")


def outputs_for(d: Path) -> list:
    """Named outputs that exist, looked up at known locations (never a recursive walk)."""
    found, seen = [], set()
    for name, why in INTERESTING:
        for sub in SUBDIRS:
            p = (d / sub / name) if sub else (d / name)
            try:
                if not p.is_file():
                    continue
                st = p.stat()
            except OSError:
                continue
            rp = str(p)
            if rp in seen:
                continue
            seen.add(rp)
            found.append({"name": name, "why": why, "path": rp,
                          "size": st.st_size, "mtime": st.st_mtime,
                          "rel": str(p.relative_to(d))})
    return found


def _safe(path_str: str) -> Path | None:
    """Resolve a requested path, refusing anything outside the configured roots."""
    try:
        p = Path(urllib.parse.unquote(path_str)).resolve()
    except Exception:
        return None
    for root in ROOTS:
        try:
            p.relative_to(root.resolve())
            return p if p.exists() else None
        except ValueError:
            continue
    return None


CSS = """
:root{--bg:#0f1216;--card:#171b21;--fg:#e6e9ee;--dim:#98a2b3;--line:#242a33;--ok:#3fb950;
--warn:#d29922;--bad:#f85149;--accent:#4c8dff}
@media (prefers-color-scheme:light){:root{--bg:#f6f7f9;--card:#fff;--fg:#11151a;--dim:#5b6572;
--line:#e3e7ec}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);
font:14px/1.5 ui-sans-serif,-apple-system,Segoe UI,Roboto,sans-serif}
header{padding:18px 22px;border-bottom:1px solid var(--line);display:flex;
align-items:baseline;gap:14px;flex-wrap:wrap}
h1{font-size:17px;margin:0;font-weight:650}
.dim{color:var(--dim)}.wrap{padding:18px 22px;max-width:1180px}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px 16px;
margin:0 0 14px}
.row{display:flex;justify-content:space-between;gap:12px;align-items:center;flex-wrap:wrap}
.badge{font-size:11px;padding:2px 8px;border-radius:99px;border:1px solid var(--line);
color:var(--dim)}
.run{background:rgba(63,185,80,.15);color:var(--ok);border-color:transparent}
.done{background:rgba(76,141,255,.15);color:var(--accent);border-color:transparent}
.stale{background:rgba(210,153,34,.15);color:var(--warn);border-color:transparent}
.bar{height:7px;background:var(--line);border-radius:99px;overflow:hidden;margin:7px 0 3px}
.bar>i{display:block;height:100%;background:var(--accent)}
table{border-collapse:collapse;width:100%;margin-top:8px}
td,th{padding:5px 8px;border-bottom:1px solid var(--line);text-align:left;font-size:13px}
th{color:var(--dim);font-weight:500}
a{color:var(--accent);text-decoration:none}a:hover{text-decoration:underline}
code{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:12px;color:var(--dim)}
.mono{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:12px}
.hits{color:var(--ok);font-weight:600}
.fgrid{display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));gap:10px}
.fgrid label{display:flex;flex-direction:column;gap:4px;font-size:12px;color:var(--dim)}
input[type=text],input:not([type]),input[type=number]{background:var(--bg);color:var(--fg);
border:1px solid var(--line);border-radius:6px;padding:7px 9px;font-size:13px;width:100%}
.cb{display:inline-flex;align-items:center;gap:6px;margin:10px 14px 0 0;font-size:13px;
color:var(--fg)}
button{margin-top:12px;background:var(--accent);color:#fff;border:0;border-radius:7px;
padding:9px 18px;font-size:13px;font-weight:600;cursor:pointer}
button:hover{filter:brightness(1.1)}
.dbgrid{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:4px}
.dbrow{display:flex;align-items:flex-start;gap:7px;font-size:12.5px;padding:5px 7px;
border:1px solid var(--line);border-radius:6px}
.dbrow.heavy{border-color:rgba(210,153,34,.5)}
.warn{color:var(--warn)}
"""


def page() -> str:
    procs = running_processes()
    items = discover()
    now = datetime.now().strftime("%H:%M:%S")
    p = [f"<!doctype html><meta charset=utf-8><title>HMM Homologue Finder — runs</title>",
         "<meta http-equiv=refresh content=15>",
         f"<style>{CSS}</style>",
         "<header><h1>HMM Homologue Finder</h1>",
         f"<span class=dim>live run monitor · refreshed {now} · auto every 15s</span></header>",
         "<div class=wrap>"]

    # Prefill the seed field with a FASTA already on disk, so the commonest case is one click.
    # Look in the served roots AND in Downloads, which is where seed files usually arrive.
    _seed = ""
    _look = list(ROOTS) + [Path.home() / "Downloads", Path("/mnt/c/Users") / os.environ.get(
        "WINUSER", "") / "Downloads"]
    for r in _look:
        try:
            if not r.is_dir():
                continue
            cands = sorted(list(r.glob("*.fasta")) + list(r.glob("*.faa"))
                           + list(r.glob("*.fa")), key=lambda x: -x.stat().st_mtime)
        except OSError:
            continue
        if cands:
            _seed = str(cands[0])
            break
    p.append(LAUNCH_FORM.format(fasta=html.escape(_seed),
                                outroot=html.escape(str(Path.home() / "hmm_runs")),
                                dbs=_db_checkboxes()))

    p.append("<div class=card><div class=row><b>Active processes</b>"
             f"<span class='badge {"run" if procs else ""}'>{len(procs)} running</span></div>")
    if procs:
        p.append("<table><tr><th>pid</th><th>what</th><th>elapsed</th><th>command</th></tr>")
        for pr in procs:
            p.append(f"<tr><td class=mono>{pr['pid']}</td><td>{html.escape(pr['what'])}</td>"
                     f"<td>{_fmt_dur(pr['elapsed'])}</td>"
                     f"<td><code>{html.escape(pr['cmd'])}</code></td></tr>")
        p.append("</table>")
    else:
        p.append("<div class=dim>Nothing running. Finished outputs are still listed below.</div>")
    p.append("</div>")

    for it in items:
        d = Path(it["path"])
        badge = "done" if it.get("done") else ("run" if procs else "")
        label = "complete" if it.get("done") else ("in progress" if procs else "idle")
        p.append(f"<div class=card><div class=row><b>{html.escape(it['name'])}</b>"
                 f"<span class='badge {badge}'>{it['kind']} · {label}</span></div>"
                 f"<div class=dim><code>{html.escape(it['path'])}</code></div>")
        if it.get("tail"):
            p.append(f"<div class=dim style='margin-top:6px'>{html.escape(it['tail'])}</div>")
        for pr in it["progress"]:
            pct = pr["pct"]
            # "finished" must come from the scan saying so, never from a percentage — a scan
            # that dies at 99 % is not finished, and calling it so is how a partial result
            # gets published as complete.
            eta = ("complete" if pr.get("complete") else
                   (f"~{_fmt_dur(pr['eta'])} left" if pr["eta"] else
                    ("incomplete — not finished" if pct and pct >= 99 else "estimating…")))
            width = min(100, pct or 0)
            p.append(f"<div style='margin-top:10px'><div class=row>"
                     f"<span>{html.escape(pr['name'])} · "
                     f"{pr['done']:,}{'/' + format(pr['total'], ',') if pr['total'] else ''} contigs"
                     f"{' · ' + str(pr['rate_per_min']) + '/min' if pr['rate_per_min'] else ''}</span>"
                     f"<span class='{'hits' if pr['hits'] else 'dim'}'>{pr['hits']} hits</span></div>"
                     f"<div class=bar><i style='width:{width}%'></i></div>"
                     f"<div class=dim>{'~' if pr.get('approx') else ''}"
                     f"{pct if pct is not None else '—'}% · {eta}"
                     f"{' · <span class=badge style=\"color:var(--warn)\">stale</span>' if pr['stale'] else ''}"
                     "</div></div>")
        outs = outputs_for(d)
        if outs:
            p.append("<table><tr><th>output</th><th>what it is</th><th>size</th>"
                     "<th>updated</th><th></th></tr>")
            for o in outs:
                q = urllib.parse.quote(o["path"])
                when = datetime.fromtimestamp(o["mtime"]).strftime("%d %b %H:%M")
                view = (f"<a href='/view?path={q}'>view</a> · " if o["name"].endswith(
                    (".html", ".csv", ".tsv", ".json")) else "")
                p.append(f"<tr><td class=mono>{html.escape(o['rel'])}</td>"
                         f"<td class=dim>{html.escape(o['why'])}</td>"
                         f"<td class=dim>{_fmt_size(o['size'])}</td>"
                         f"<td class=dim>{when}</td>"
                         f"<td>{view}<a href='/download?path={q}'>download</a></td></tr>")
            p.append("</table>")
        else:
            p.append("<div class=dim style='margin-top:8px'>No named outputs yet.</div>")
        p.append("</div>")

    if not items:
        p.append("<div class=card class=dim>No runs found under: "
                 + ", ".join(f"<code>{html.escape(str(r))}</code>" for r in ROOTS) + "</div>")
    p.append("</div>")
    return "".join(p)


# The status snapshot is built on a timer in the background and served from memory.
#
# Building it on demand does not work once a root lives on /mnt/c: every stat crosses WSL's
# filesystem bridge, and a page that renders in 0.15 s against native Linux paths takes 16.4 s
# against a Windows directory — longer than the refresh interval, so the browser would sit
# permanently loading. Decoupling means a page load is a dictionary lookup no matter how slow
# the filesystem is, and a slow scan makes the data a few seconds stale rather than unusable.
_CACHE = {"html": None, "json": None, "ts": 0.0, "secs": 0.0, "error": None}
_LOCK = threading.Lock()


def _rebuild() -> None:
    t0 = time.time()
    try:
        h, j = page(), {"processes": running_processes(), "runs": discover()}
        with _LOCK:
            _CACHE.update(html=h, json=j, ts=time.time(),
                          secs=round(time.time() - t0, 2), error=None)
    except Exception as e:                       # never let the refresher die
        with _LOCK:
            _CACHE["error"] = f"{type(e).__name__}: {e}"


def _refresher(interval: float) -> None:
    while True:
        _rebuild()
        time.sleep(interval)


_WAIT = ("<!doctype html><meta charset=utf-8><title>collecting…</title>"
         "<meta http-equiv=refresh content=3>"
         f"<style>{CSS}</style><div class=wrap><div class=card>"
         "<b>Collecting run status…</b><div class=dim>First scan of the run directories. "
         "Windows paths under /mnt/c are slow to stat, so this can take a few seconds.</div>"
         "</div></div>")


def _render_table(p: Path, limit: int = 400) -> str:
    """Render a CSV/TSV as an HTML table. Viewing a results table as raw text is close to
    useless for the tables this pipeline produces — 23-column homolog rows with embedded
    sequences — so the point of a viewer is to make them readable."""
    delim = "\t" if p.suffix.lower() in (".tsv", ".tbl") else ","
    try:
        with p.open("r", encoding="utf-8", errors="replace", newline="") as fh:
            rows = list(csv.reader(fh, delimiter=delim))
    except OSError as e:
        return f"<p>could not read: {html.escape(str(e))}</p>"
    if not rows:
        return "<p class=dim>empty file</p>"
    head, body = rows[0], rows[1:]
    out = [f"<h1>{html.escape(p.name)}</h1>",
           f"<p class=dim>{len(body):,} row(s) × {len(head)} column(s)"
           + (f" — showing the first {limit}" if len(body) > limit else "") + " · "
           f"<a href='/download?path={urllib.parse.quote(str(p))}'>download</a></p>",
           "<div style='overflow-x:auto'><table><tr>"]
    out += [f"<th>{html.escape(c)}</th>" for c in head]
    out.append("</tr>")
    for r in body[:limit]:
        out.append("<tr>" + "".join(
            # long cells are almost always a sequence; truncate so the table stays readable
            f"<td>{html.escape(c[:120] + ('…' if len(c) > 120 else ''))}</td>" for c in r) + "</tr>")
    out.append("</table></div>")
    return "".join(out)


# The set the reference and poster-trial runs used, and the one whose numbers are validated.
# Pre-ticking it matters: a run launched with NO --databases falls back to the generic default
# of three, silently dropping "RefSeq viral genomes" — which contributed 39 of the 155 hits.
# A GUI that quietly searches less than the run you validated is worse than no GUI.
VALIDATED_SET = ("INPHARED genomes", "RefSeq viral genomes", "INPHARED proteins",
                 "SwissProt", "RefSeq viral proteins", "VOGDB VFAM (annotation)")
# Flagged in the UI rather than hidden: these are real but not laptop jobs.
HEAVY = {"RefSeq bacterial genomes": "~21 DAYS on 4 cores — server-scale",
         "RefSeq bacterial proteins": "~80 GB download",
         "Pfam (sequences)": "~6 GB",
         "GVD-AVrC": "~5 h", "Gut Phage Database (GPD)": "~2 h"}


def catalogue() -> list:
    """The database catalog as [(name, type, size_hint, default_on, warning)]."""
    try:
        eng = str(Path(__file__).resolve().parents[1] / "engine")
        if eng not in sys.path:
            sys.path.insert(0, eng)
        from databases.builtin import BUILTIN_DATABASES
    except Exception:
        return []
    out = []
    for d in BUILTIN_DATABASES:
        n = str(d.get("name", ""))
        out.append((n, str(d.get("db_type", "") or ""), str(d.get("size_hint", "") or ""),
                    n in VALIDATED_SET, HEAVY.get(n, "")))
    return out


def _db_checkboxes() -> str:
    rows = catalogue()
    if not rows:
        return "<div class=dim>database catalog unavailable</div>"
    out = ["<div class=dbgrid>"]
    for name, typ, size, on, warn in rows:
        cid = "db_" + re.sub(r"[^A-Za-z0-9]", "_", name)
        cls = " heavy" if warn else ""
        out.append(
            f"<label class='dbrow{cls}'><input type=checkbox name=db value=\"{html.escape(name)}\""
            f"{' checked' if on else ''} id={cid}>"
            f"<span><b>{html.escape(name)}</b>"
            f"<span class=dim> · {html.escape(size)}</span>"
            + (f"<span class=warn> · {html.escape(warn)}</span>" if warn else "")
            + "</span></label>")
    out.append("</div>")
    return "".join(out)


LAUNCH_FORM = """
<div class=card>
  <div class=row><b>Start a discovery run</b>
    <span class=badge>runs detached — closing this page will not stop it</span></div>
  <form method=post action=/launch style='margin-top:10px'>
    <div class=fgrid>
      <label>Seed FASTA (full path)<input name=fasta required
        placeholder="/mnt/c/Users/you/Downloads/seeds.fasta" value="{fasta}"></label>
      <label>Run name<input name=name required placeholder="my_family" value=""></label>
      <label>Output folder<input name=outdir value="{outroot}"></label>
      <label>NCBI e-mail <span class=dim>(blank = offline)</span><input name=email
        placeholder="you@inst.edu"></label>
      <label>Iterations<input name=iterations type=number min=1 max=6 value=3></label>
      <label>CPU<input name=cpu type=number min=1 max=32 value=4></label>
    </div>
    <label class=cb><input type=checkbox name=find_interrupted checked> find interrupted /
      overprinted homologs</label>
    <label class=cb><input type=checkbox name=clear_cache> clear the database cache when the
      run finishes</label>
    <div style='margin-top:14px'><b style='font-size:13px'>Databases to search</b>
      <div class=dim style='margin:2px 0 8px'>Pre-ticked is the validated set the reference and
        poster runs used. Leaving every box unticked falls back to the built-in default of
        three, which drops <i>RefSeq viral genomes</i> — 39 of 155 hits came from it.</div>
      {dbs}
    </div>
    <button type=submit>Start run</button>
  </form>
  <div class=dim style='margin-top:8px'>Writes to a Linux path by default: more space, and no
    /mnt/c filesystem bridge. The run appears below within a few seconds.</div>
</div>
"""


def launch_run(form: dict) -> tuple:
    """Start a discovery run detached. Returns (ok, message).

    Every value goes into an explicit argv list — never a shell string — so nothing a browser
    field contains can become a command. The seed path must exist before anything is started.
    """
    fasta = (form.get("fasta") or "").strip()
    name = re.sub(r"[^A-Za-z0-9._-]", "_", (form.get("name") or "").strip())
    if not fasta or not Path(fasta).is_file():
        return False, f"Seed FASTA not found: {fasta or '(blank)'}"
    if not name:
        return False, "Run name is required"
    outroot = Path((form.get("outdir") or str(Path.home() / "hmm_runs")).strip()).expanduser()
    out = outroot / name
    try:
        outroot.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        return False, f"Cannot create {outroot}: {e}"
    # Refuse to start a run that will die on the engine's own disk check 5 minutes in.
    try:
        free_gb = shutil.disk_usage(outroot).free / 1e9
        if free_gb < 25:
            return False, (f"Only {free_gb:.1f} GB free at {outroot} — the engine requires 20 GB "
                           f"and will abort. Choose a location with more space.")
    except OSError:
        pass

    here = Path(__file__).resolve().parent
    argv = [sys.executable, str(here / "hmm_finder.py"),
            "--fasta", fasta, "--out-dir", str(out),
            "--iterations", str(int(form.get("iterations") or 3)),
            "--cpu", str(int(form.get("cpu") or 4)), "--no-overwrite"]
    # Databases: pass the tick-box selection through explicitly. Without --databases the engine
    # uses its built-in default set, which is SMALLER than the validated one — a run launched
    # here would then quietly search less than the run these results were validated against.
    dbs = [d for d in (form.get("db") or []) if d]
    known = {n for n, *_ in catalogue()}
    unknown = [d for d in dbs if known and d not in known]
    if unknown:
        return False, f"Unknown database(s): {', '.join(unknown)}"
    if dbs:
        argv += ["--databases", ",".join(dbs)]
    if (form.get("email") or "").strip():
        argv += ["--email", form["email"].strip()]
    if form.get("find_interrupted"):
        argv.append("--find-interrupted")
    if form.get("clear_cache"):
        argv.append("--clear-cache")

    log = outroot / f"{name}.log"
    try:
        with log.open("wb") as fh:
            subprocess.Popen(argv, stdout=fh, stderr=subprocess.STDOUT,
                             stdin=subprocess.DEVNULL,
                             start_new_session=True,      # survives this server exiting
                             cwd=str(here.parent))
    except Exception as e:
        return False, f"Could not start: {e}"
    return True, f"Started “{name}” → {out} (log: {log})"


class Handler(SimpleHTTPRequestHandler):
    def log_message(self, *a):        # keep the console readable
        pass

    def _send(self, body: bytes, ctype="text/html; charset=utf-8", extra=None):
        self.send_response(200)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        u = urllib.parse.urlparse(self.path)
        q = urllib.parse.parse_qs(u.query)
        if u.path == "/":
            with _LOCK:
                h, ts, secs, err = (_CACHE["html"], _CACHE["ts"], _CACHE["secs"],
                                    _CACHE["error"])
            if not h:
                return self._send(_WAIT.encode("utf-8"))
            age = int(time.time() - ts)
            note = (f"<div class='wrap dim' style='padding-top:0'>snapshot {age}s old · "
                    f"built in {secs}s"
                    + (f" · <span style='color:var(--bad)'>{html.escape(err)}</span>" if err else "")
                    + "</div>")
            return self._send((h + note).encode("utf-8"))
        if u.path == "/api/status":
            with _LOCK:
                j, ts = _CACHE["json"], _CACHE["ts"]
            if j is None:
                return self._send(b'{"status":"collecting"}', "application/json")
            return self._send(json.dumps({**j, "snapshot_age_s": round(time.time() - ts)},
                                         default=str).encode(), "application/json")
        if u.path in ("/download", "/view"):
            p = _safe((q.get("path") or [""])[0])
            if not p or not p.is_file():
                self.send_error(404, "not found, or outside the served roots")
                return
            ctype = mimetypes.guess_type(p.name)[0] or "application/octet-stream"
            if u.path == "/view":
                sfx = p.suffix.lower()
                # Render, don't dump. A 23-column homolog table or a genome map is unreadable
                # as raw bytes, and being able to LOOK at an output is the point of a viewer.
                if sfx in (".csv", ".tsv", ".tbl"):
                    page_ = (f"<!doctype html><meta charset=utf-8><title>{html.escape(p.name)}"
                             f"</title><style>{CSS}</style><div class=wrap>"
                             f"<p><a href='/'>&larr; back</a></p>{_render_table(p)}</div>")
                    return self._send(page_.encode("utf-8"))
                if sfx in (".png", ".svg", ".pdf", ".jpg", ".jpeg"):
                    return self._send(p.read_bytes(), ctype)     # browser renders these
                if sfx in (".html", ".htm"):
                    return self._send(p.read_bytes(), "text/html; charset=utf-8")
                if sfx in (".json", ".txt", ".log", ".md", ".sto", ".faa", ".fna", ".treefile"):
                    ctype = "text/plain; charset=utf-8"
                return self._send(p.read_bytes(), ctype)
            return self._send(p.read_bytes(), ctype,
                              {"Content-Disposition": f'attachment; filename="{p.name}"'})
        self.send_error(404)

    def do_POST(self):
        """Only one action exists: start a discovery run. There is no generic command endpoint."""
        u = urllib.parse.urlparse(self.path)
        if u.path != "/launch":
            return self.send_error(404)
        try:
            n = int(self.headers.get("Content-Length") or 0)
            body = self.rfile.read(min(n, 65536)).decode("utf-8", "replace")
            # `db` is a multi-value field (one per ticked database), so it must keep its list.
            parsed = urllib.parse.parse_qs(body)
            form = {k: (v if k == "db" else v[0]) for k, v in parsed.items()}
        except Exception:
            form = {}
        ok, msg = launch_run(form)
        if ok:
            _rebuild()            # so the new run appears immediately rather than in 12s
        colour = "var(--ok)" if ok else "var(--bad)"
        body = (f"<!doctype html><meta charset=utf-8><title>{'started' if ok else 'not started'}"
                f"</title><meta http-equiv=refresh content='4;url=/'>"
                f"<style>{CSS}</style><div class=wrap><div class=card>"
                f"<b style='color:{colour}'>{'Run started' if ok else 'Could not start'}</b>"
                f"<div style='margin-top:8px'>{html.escape(msg)}</div>"
                f"<div class=dim style='margin-top:10px'>returning to the dashboard…</div>"
                f"</div></div>")
        self._send(body.encode("utf-8"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--interval", type=float, default=12.0,
                    help="seconds between background status rebuilds (default 12)")
    ap.add_argument("--root", action="append", default=[],
                    help="directory to scan for runs (repeatable)")
    a = ap.parse_args()
    # Defaults are deliberately shallow. /mnt/c/Users would put every Windows profile directory
    # in the scan; the runs live in Downloads, so name that instead.
    roots = a.root or [str(Path.home() / "hmm_runs"),
                       "/mnt/c/Users/" + os.environ.get("WINUSER", "") if os.environ.get("WINUSER")
                       else str(Path.home() / "hmm_runs"),
                       str(Path.cwd())]
    roots = [r for r in dict.fromkeys(roots) if r]
    global ROOTS
    ROOTS = [Path(r).expanduser() for r in roots if Path(r).expanduser().is_dir()]
    if not ROOTS:
        print("no readable --root directories", file=sys.stderr)
        return 2
    print(f"serving {len(ROOTS)} root(s):")
    for r in ROOTS:
        print(f"  {r}")
    print(f"\n  ->  http://127.0.0.1:{a.port}\n")
    print("  (localhost only; GET only; paths confined to the roots above)")
    th = threading.Thread(target=_refresher, args=(a.interval,), daemon=True)
    th.start()
    socketserver.ThreadingTCPServer.allow_reuse_address = True
    with socketserver.ThreadingTCPServer(("127.0.0.1", a.port), Handler) as srv:
        try:
            srv.serve_forever()
        except KeyboardInterrupt:
            print("\nstopped")
    return 0


if __name__ == "__main__":
    sys.exit(main())
