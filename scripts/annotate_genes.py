#!/usr/bin/env python3
"""
annotate_genes.py — functional annotation of proteins via VOGDB VFAM (hmmscan).

Used to put real function names ("major capsid protein", "terminase large
subunit", …) on neighbourhood genes for the publication synteny figures. The
VOGDB VFAM HMM library is downloaded once into the shared db-cache
(<cache>/annotation/vogdb) and reused across runs.
"""
from __future__ import annotations

import json
import re
import subprocess
import tempfile
from pathlib import Path

from env_paths import ensure_env_on_path  # noqa: E402
ensure_env_on_path()

VFAM_HMM_URL = "https://fileshare.csb.univie.ac.at/vog/vog230/vfam.hmm.tar.gz"
VFAM_ANN_URL = "https://fileshare.csb.univie.ac.at/vog/vog230/vfam.annotations.tsv.gz"


def vogdb_dir(cache: Path) -> Path:
    return Path(cache).expanduser() / "annotation" / "vogdb"


def is_ready(cache: Path) -> bool:
    d = vogdb_dir(cache)
    return (d / "vfam_all.hmm.h3m").exists() and (d / "vfam.annotations.tsv").exists()


def setup(cache: Path) -> bool:
    """Download + index VOGDB VFAM into the cache (idempotent). Returns is_ready."""
    d = vogdb_dir(cache)
    d.mkdir(parents=True, exist_ok=True)
    if is_ready(cache):
        return True
    try:
        ann_gz = d / "vfam.annotations.tsv.gz"
        if not (d / "vfam.annotations.tsv").exists():
            subprocess.run(["curl", "-sSL", "-o", str(ann_gz), VFAM_ANN_URL], check=True)
            subprocess.run(["gunzip", "-kf", str(ann_gz)], check=True)
        if not (d / "vfam_all.hmm").exists():
            tgz = d / "vfam.hmm.tar.gz"
            subprocess.run(["curl", "-sSL", "-o", str(tgz), VFAM_HMM_URL], check=True)
            subprocess.run(["tar", "-xzf", str(tgz), "-C", str(d)], check=True)
            allhmm = d / "vfam_all.hmm"
            with allhmm.open("wb") as out:
                for hmm in sorted(d.rglob("*.hmm")):
                    if hmm.name == "vfam_all.hmm":
                        continue
                    out.write(hmm.read_bytes())
        subprocess.run(["hmmpress", "-f", str(d / "vfam_all.hmm")],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        # reclaim space: per-VFAM files + tarballs are redundant once pressed
        import shutil as _sh
        for junk in (d / "vfam.hmm.tar.gz", d / "vfam.annotations.tsv.gz"):
            junk.unlink(missing_ok=True)
        if (d / "hmm").is_dir():
            _sh.rmtree(d / "hmm", ignore_errors=True)
        # record provenance so the figure's functional annotation is citable
        (d / "provenance.json").write_text(json.dumps({
            "database": "VOGDB VFAM",
            "release": "vog230 (VOGDB release 230)",
            "hmm_url": VFAM_HMM_URL,
            "annotations_url": VFAM_ANN_URL,
            "tool": "HMMER hmmscan",
        }, indent=2))
    except Exception as e:
        print(f"  (VOGDB setup failed: {e})")
        return False
    return is_ready(cache)


def _clean_desc(s: str) -> str:
    s = (s or "").strip()
    m = re.match(r"^(?:sp|tr)\|[^|]+\|\S+\s+(.*)$", s)
    if m:
        s = m.group(1)
    s = re.sub(r"^(Putative|Predicted|Probable|Uncharacterized)\s+", "", s, flags=re.I).strip()
    return s or "hypothetical protein"


def load_annotations(cache: Path) -> dict:
    f = vogdb_dir(cache) / "vfam.annotations.tsv"
    out = {}
    if not f.exists():
        return out
    with f.open() as fh:
        next(fh, None)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5:
                out[p[0]] = {"function": _clean_desc(p[4]), "category": p[3]}
    return out


def annotate(proteins: dict, cache: Path, cpu: int = 4, evalue: float = 1e-3) -> dict:
    """proteins: {id: aa}. Returns {id: {'vfam','function','category'}} for hits."""
    if not proteins or not is_ready(cache):
        return {}
    ann = load_annotations(cache)
    hmm = vogdb_dir(cache) / "vfam_all.hmm"
    res = {}
    with tempfile.TemporaryDirectory() as td:
        faa, tbl = Path(td) / "q.faa", Path(td) / "out.tbl"
        with faa.open("w") as fh:
            for k, aa in proteins.items():
                if aa:
                    fh.write(f">{k}\n{aa}\n")
        try:
            subprocess.run(["hmmscan", "--tblout", str(tbl), "-E", str(evalue),
                            "--cpu", str(cpu), "--noali", str(hmm), str(faa)],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        except Exception:
            return {}
        best = {}
        for line in tbl.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            if len(p) < 5:
                continue
            vfam, q = p[0], p[2]
            try:
                ev = float(p[4])
            except ValueError:
                continue
            if q not in best or ev < best[q][1]:
                best[q] = (vfam, ev)
        for q, (vfam, _) in best.items():
            a = ann.get(vfam, {})
            res[q] = {"vfam": vfam,
                      "function": a.get("function", "hypothetical protein"),
                      "category": a.get("category", "")}
    return res
