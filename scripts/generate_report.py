#!/usr/bin/env python3
"""
generate_report.py — one-page, self-contained HTML summary for a discovery run.

Reads a <name>_discovery directory and writes report.html at its root: headline
stats, file links, per-run summary, the phylogenetic tree (embedded), the top
homologs, and tool versions. The tree image is base64-embedded so the single file
is portable. Never raises on a missing input.
"""
from __future__ import annotations

import argparse
import base64
import csv
import html
import json
from pathlib import Path

from package_layout import DIRS  # PACKAGE/ folder names (single source of truth)


def _read_csv(p: Path) -> list[dict]:
    try:
        with open(p, newline="") as f:
            return list(csv.DictReader(f))
    except Exception:
        return []


def _read_tsv(p: Path) -> list[dict]:
    try:
        with open(p, newline="") as f:
            return list(csv.DictReader(f, delimiter="\t"))
    except Exception:
        return []


def _b64(p: Path) -> str:
    try:
        return base64.b64encode(p.read_bytes()).decode()
    except Exception:
        return ""


_CSS = (
    "body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:0;"
    "background:#f6f7f9;color:#1a2230}.wrap{max-width:980px;margin:0 auto;padding:28px 22px 70px}"
    "h1{font-size:24px;margin:0 0 4px}.sub{color:#5f6b7a;margin:0 0 20px;font-size:14px}"
    ".cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:0 0 24px}"
    ".card{background:#fff;border:1px solid #e3e7ec;border-radius:10px;padding:12px 14px}"
    ".card .l{font-size:11px;color:#6b7888;text-transform:uppercase;letter-spacing:.04em}"
    ".card .v{font-size:18px;font-weight:600;margin-top:3px;word-break:break-word}"
    "h2{font-size:17px;margin:28px 0 10px}"
    "table{width:100%;border-collapse:collapse;background:#fff;border:1px solid #e3e7ec;border-radius:8px}"
    "th,td{text-align:left;padding:7px 10px;border-bottom:1px solid #eef1f4;font-size:13px}"
    "th{background:#f0f3f6;color:#33404f}"
    "a.btn{display:inline-block;margin:4px 8px 4px 0;padding:8px 14px;background:#fff;"
    "border:1px solid #d4dae1;border-radius:8px;text-decoration:none;color:#1558b0;font-size:13px}"
    "img{max-width:100%;border:1px solid #e3e7ec;border-radius:8px;background:#fff}"
    ".muted{color:#6b7888;font-size:13px}"
    ".msa{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:11px;line-height:1.3;"
    "white-space:pre;overflow-x:auto;border:1px solid #e3e7ec;border-radius:8px;background:#fff;padding:8px}"
    ".msa .lbl{color:#6b7888;display:inline-block;width:16ch;overflow:hidden;vertical-align:bottom}"
    ".msa .r{padding:0 .5px}"
    ".hydro{background:#bcd6ff}.pos{background:#ffc9c9}.neg{background:#efbdff}"
    ".polar{background:#c5efc5}.gly{background:#ffe0ad}.pro{background:#fff0a6}"
    ".his{background:#cdddff}.tyr{background:#cdf3e8}"
)


def _read_fasta(p: Path) -> list:
    """Minimal (name, sequence) FASTA reader — no Bio dependency in the report."""
    recs, name, seq = [], None, []
    try:
        with open(p) as f:
            for ln in f:
                if ln.startswith(">"):
                    if name is not None:
                        recs.append((name, "".join(seq)))
                    name, seq = ln[1:].rstrip("\n"), []
                elif name is not None:
                    seq.append(ln.strip())
        if name is not None:
            recs.append((name, "".join(seq)))
    except Exception:
        return []
    return recs


# Clustal-ish residue classes -> CSS class (for the inline coloured MSA)
_AA_CLASS = {a: g for g, aas in {
    "hydro": "AILMFWVC", "pos": "KR", "neg": "DE", "polar": "STNQ",
    "gly": "G", "pro": "P", "his": "H", "tyr": "Y"}.items() for a in aas}


def generate(discovery: Path) -> Path:
    discovery = Path(discovery)
    e = html.escape

    manifest = {}
    if (discovery / "run_manifest.json").exists():
        try:
            manifest = json.loads((discovery / "run_manifest.json").read_text())
        except Exception:
            manifest = {}
    params = manifest.get("parameters", {})
    label = params.get("label", discovery.name)
    summary = _read_csv(discovery / "hit_summary.csv")
    paper = _read_csv(discovery / "paper_main_table.csv")

    tree_png = discovery / "downstream" / "tree" / "hits_tree.png"
    if not tree_png.exists():
        tree_png = discovery / "PACKAGE" / DIRS["phylo"] / "hits_tree.png"
    tree_b64 = _b64(tree_png) if tree_png.exists() else ""

    aln_png = discovery / "downstream" / "tree" / "alignment_figure.png"
    if not aln_png.exists():
        aln_png = discovery / "PACKAGE" / DIRS["phylo"] / "alignment_figure.png"
    aln_b64 = _b64(aln_png) if aln_png.exists() else ""
    aln_stats = {}
    for sp in (discovery / "downstream" / "tree" / "hits.aln.stats.json",
               discovery / "PACKAGE" / DIRS["phylo"] / "hits.aln.stats.json"):
        if sp.exists():
            try:
                aln_stats = json.loads(sp.read_text())
            except Exception:
                aln_stats = {}
            break

    last = summary[-1] if summary else {}
    cards = [
        ("Family", label),
        ("Iterations", params.get("iterations", "")),
        ("Databases", params.get("databases", "")),
        ("Hits (final run)", last.get("total_hits", "")),
        ("Unique homologs", last.get("unique_sequences", "")),
        ("Organisms", last.get("unique_organisms", "")),
    ]

    p: list[str] = ["<!DOCTYPE html><html><head><meta charset='utf-8'>",
                    "<meta name='viewport' content='width=device-width, initial-scale=1'>",
                    f"<title>{e(label)} — discovery report</title><style>{_CSS}</style>",
                    "</head><body><div class='wrap'>",
                    f"<h1>{e(label)} — homologue discovery report</h1>",
                    f"<p class='sub'>{e(str(manifest.get('started_at','')))} → "
                    f"{e(str(manifest.get('finished_at','')))} · commit "
                    f"{e(str(manifest.get('code_git_commit','?')))}</p>",
                    "<div class='cards'>"]
    for l, v in cards:
        p.append(f"<div class='card'><div class='l'>{e(str(l))}</div>"
                 f"<div class='v'>{e(str(v))}</div></div>")
    p.append("</div>")

    cal = manifest.get("threshold_calibration", {}) or {}
    stop = manifest.get("iteration_stop_reason", "")
    if cal or stop:
        p.append("<h2>Calibration &amp; convergence</h2>")
        if stop:
            p.append(f"<p class='muted'>Iteration stopping criterion: {e(str(stop))}</p>")
        if cal:
            p.append("<div class='cards'>")
            _roc = cal.get("roc") or {}
            for lbl, k in [("Sensitivity (seeds)", "sensitivity"),
                           ("Specificity (controls)", "specificity"),
                           ("False-positive rate", "false_positive_rate")]:
                if k in cal:
                    p.append(f"<div class='card'><div class='l'>{e(lbl)}</div>"
                             f"<div class='v'>{e(str(cal.get(k)))}</div></div>")
            if _roc.get("auc") is not None:
                p.append(f"<div class='card'><div class='l'>ROC AUC</div>"
                         f"<div class='v'>{e(str(_roc.get('auc')))}</div></div>")
                p.append(f"<div class='card'><div class='l'>Optimal bit (Youden)</div>"
                         f"<div class='v'>{e(str(_roc.get('optimal_threshold')))}</div></div>")
            p.append("</div>")
            if _roc.get("auc") is not None:
                p.append(f"<p class='muted'>ROC over {e(str(_roc.get('n_positive','?')))} positive vs "
                         f"{e(str(_roc.get('n_negative','?')))} negative control sequences (advisory; the "
                         f"fixed strict threshold is kept for tiering). Curve: "
                         f"<code>controls/roc_curve.svg</code>.</p>")
            p.append(f"<p class='muted'>At the strict bit-score threshold: "
                     f"{e(str(cal.get('true_positives','?')))}/{e(str(cal.get('total_positives','?')))} "
                     f"seed sequences recovered; {e(str(cal.get('false_positives','?')))}/"
                     f"{e(str(cal.get('total_negatives','?')))} negative-control sequences scored "
                     f"above threshold. Detail: <code>controls/control_report.json</code>.</p>")

    p.append("<h2>Files</h2><div>")
    for txt, href in [("Main table (CSV)", "paper_main_table.csv"),
                      ("All hits — supplementary (CSV)", "all_runs_hits.csv"),
                      ("Per-run summary (CSV)", "hit_summary.csv"),
                      ("Database provenance (CSV)", "database_summary.csv"),
                      ("Methods", "METHODS.md"),
                      ("Interrupted/overprinted homologs (TSV)", "interrupted_homologs.tsv"),
                      ("Interrupted homolog proteins — domain (FASTA)", "interrupted_homologs_domain_aa.faa"),
                      ("Interrupted homolog proteins — full ORF (FASTA)", "interrupted_homologs_full_orf_aa.faa"),
                      ("Threshold calibration (JSON)", "controls/control_report.json"),
                      ("Alignment figure (SVG)", "downstream/tree/alignment_figure.svg"),
                      ("Alignment (FASTA)", "downstream/tree/hits.aln.faa"),
                      ("Per-hit HMM alignment (Stockholm)", "downstream/tree/hits_hmmalign.sto"),
                      ("Per-hit HMM alignment (A2M)", "downstream/tree/hits_hmmalign.a2m"),
                      ("Publication synteny figures", "downstream/synteny/index.html"),
                      ("Synteny figures (clinker, interactive)", "downstream/clinker/index.html")]:
        if (discovery / href).exists():
            p.append(f"<a class='btn' href='{e(href)}'>{e(txt)}</a>")
    p.append("</div>")

    if summary:
        cols = list(summary[0].keys())
        p.append("<h2>Per-run summary</h2><table><tr>")
        p += [f"<th>{e(c)}</th>" for c in cols]
        p.append("</tr>")
        for r in summary:
            p.append("<tr>" + "".join(f"<td>{e(str(r.get(c, '')))}</td>" for c in cols) + "</tr>")
        p.append("</table>")

    if aln_b64 or aln_stats:
        p.append("<h2>Multiple sequence alignment (unique homologs)</h2>")
        if aln_stats:
            p.append(
                f"<p class='muted'>{e(str(aln_stats.get('n_sequences','?')))} sequences × "
                f"{e(str(aln_stats.get('aln_length','?')))} columns · "
                f"{e(str(aln_stats.get('conserved_columns','?')))} conserved columns · "
                f"mean pairwise identity {e(str(aln_stats.get('avg_pairwise_id','?')))}% · "
                f"gaps {e(str(aln_stats.get('gap_pct','?')))}% "
                f"(MAFFT accuracy mode → trimAl; full MSA in <code>hits.aln.faa</code>).</p>")
        if aln_b64:
            p.append(f"<img alt='coloured multiple sequence alignment' "
                     f"src='data:image/png;base64,{aln_b64}'>")

    # Inline, residue-coloured MSA of the final hits (text, beyond the static figure).
    aln_fa = discovery / "downstream" / "tree" / "hits.aln.faa"
    if not aln_fa.exists():
        aln_fa = discovery / "PACKAGE" / DIRS["phylo"] / "hits.aln.faa"
    msa = _read_fasta(aln_fa)
    if msa:
        MAXSEQ, MAXCOL = 30, 180
        shown = msa[:MAXSEQ]
        ncol = max((len(s) for _, s in shown), default=0)
        col0 = min(ncol, MAXCOL)
        rows_html = []
        for name, seq in shown:
            lbl = e(name.split()[0][:16])
            cells = []
            for ch in seq[:col0]:
                cls = _AA_CLASS.get(ch.upper())
                cells.append(f"<span class='r {cls}'>{e(ch)}</span>" if cls
                             else f"<span class='r'>{e(ch)}</span>")
            rows_html.append(f"<div><span class='lbl' title='{e(name)}'>{lbl}</span> "
                             + "".join(cells) + "</div>")
        p.append("<h2>Coloured alignment (final hits)</h2>")
        note = (f"Showing {len(shown)} of {len(msa)} hit sequences"
                + (f" × first {col0} of {ncol} columns" if ncol > col0 else f" × {ncol} columns")
                + " — residues coloured by class. Full MSA in <code>hits.aln.faa</code>; "
                  "each hit aligned to the family HMM in <code>hits_hmmalign.a2m</code>.")
        p.append(f"<p class='muted'>{note}</p>")
        p.append("<div class='msa'>" + "".join(rows_html) + "</div>")

    if tree_b64:
        p.append("<h2>Phylogeny (unique homologs)</h2>")
        p.append(f"<img alt='ML tree of homologs' src='data:image/png;base64,{tree_b64}'>")

    if paper:
        show = [c for c in ("rank", "representative_organism", "accession", "copies",
                            "domain_aa_len", "best_evalue", "best_bit_score",
                            "confidence_tier") if c in paper[0]]
        p.append("<h2>Top homologs</h2><table><tr>")
        p += [f"<th>{e(c)}</th>" for c in show]
        p.append("</tr>")
        for r in paper[:25]:
            p.append("<tr>" + "".join(f"<td>{e(str(r.get(c, '')))}</td>" for c in show) + "</tr>")
        p.append("</table>")
        if len(paper) > 25:
            p.append(f"<p class='muted'>… and {len(paper) - 25} more in paper_main_table.csv</p>")

    # Stop-interrupted / overprinted homologs (only present with --find-interrupted).
    interrupted = manifest.get("interrupted_homologs", {}) or {}
    irows = _read_tsv(discovery / "interrupted_homologs.tsv")
    if interrupted.get("interrupted_candidates") or irows:
        n = interrupted.get("interrupted_candidates", len(irows))
        p.append("<h2>Stop-interrupted / overprinted homologs</h2>")
        note = (f"Read-through translation of the searched nucleotide databases (stop codons "
                f"retained, not broken on) scanned with the family HMM: <b>{e(str(n))}</b> "
                f"match(es) carry &ge;1 internal stop within the domain — candidate homologs "
                f"interrupted by a premature stop (e.g. overprinted genes whose nonsense "
                f"mutation is silent in an overlapping reading frame). The normal stop-to-stop "
                f"search cannot recover these.")
        if interrupted.get("threshold_basis"):
            note += (f" Reporting threshold {e(str(interrupted.get('min_bit')))} bits "
                     f"[{e(str(interrupted.get('threshold_basis')))}].")
        if (discovery / "interrupted_homologs_domain_aa.faa").exists():
            note += (" Protein sequences (internal stop shown as <code>*</code>) are in "
                     "<code>interrupted_homologs_domain_aa.faa</code> (domain) and "
                     "<code>interrupted_homologs_full_orf_aa.faa</code> (full ORF).")
        p.append(f"<p class='muted'>{note}</p>")
        if irows:
            show = [c for c in ("contig", "strand", "frame", "domain_nt_start",
                                "domain_nt_end", "domain_aa_len", "internal_stops",
                                "stop_nt_positions", "domain_bit_score", "i_evalue")
                    if c in irows[0]]
            p.append("<table><tr>")
            p += [f"<th>{e(c)}</th>" for c in show]
            p.append("</tr>")
            for r in irows[:25]:
                p.append("<tr>" + "".join(f"<td>{e(str(r.get(c, '')))}</td>" for c in show) + "</tr>")
            p.append("</table>")
            if len(irows) > 25:
                p.append(f"<p class='muted'>… and {len(irows) - 25} more in "
                         f"interrupted_homologs.tsv</p>")

    tv = manifest.get("tool_versions", {})
    if tv:
        p.append("<h2>Tool versions</h2><table><tr><th>Tool</th><th>Version</th></tr>")
        for t, v in sorted(tv.items()):
            ver = v.get("version", "") if isinstance(v, dict) else str(v)
            p.append(f"<tr><td>{e(t)}</td><td>{e(str(ver))}</td></tr>")
        p.append("</table>")

    p.append("</div></body></html>")
    out = discovery / "report.html"
    out.write_text("\n".join(p), encoding="utf-8")
    # mirror into PACKAGE if present
    pkg = discovery / "PACKAGE"
    if pkg.exists():
        try:
            (pkg / "report.html").write_text("\n".join(p), encoding="utf-8")
        except Exception:
            pass
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--discovery-dir", type=Path, required=True)
    args = ap.parse_args()
    print(f"  wrote {generate(args.discovery_dir)}")


if __name__ == "__main__":
    main()
