#!/usr/bin/env python3
"""
hmm_finder.py — one-command, end-to-end homologue discovery
=================================================================
Give it a single seed FASTA and it does EVERYTHING, unattended:

    HMM build  ->  six-frame search of 10 public databases  ->  ORF-validated
    sequence extraction (NT + AA + per-hit TSV)  ->  use the new sequences as
    seeds and repeat (default 3 iterations)  ->  CD-HIT clustering  ->  clinker
    synteny figures  ->  ML tree of the homologs  ->  a labelled, publication-
    ready output package.

NO interactive input is required. The ONLY required argument is --fasta.

WHY EACH STEP
-------------
* Six-frame translation lets the search find homologs encoded by genes that
  standard annotation never predicted (the homologues are such genes).
* The extractor (extract_validated_hits.py) captures the EXACT ORF the HMM
  matched, frame-correctly, and validates it is a genuine ORF (no internal
  stops, sits in a real coding locus). This is the corrected core of the
  workflow — earlier versions stored the wrong (overlapping) protein.
* Iteration tests convergence: if successive rounds stop finding new homologs,
  the family is fully captured.

REQUIREMENTS
------------
* The HMM-Discovery deployable repo (provides run_all_database_benchmark.py and
  the pipeline/ package).
* The conda env `hmm-discovery` (HMMER, MAFFT, trimAl, Prodigal, seqkit,
  CD-HIT, IQ-TREE, MEME/FIMO, clinker). The script puts the env's bin on PATH
  automatically.
* The helper scripts in the same directory: extract_validated_hits.py,
  cluster_and_clinker_corrected.py, build_tree_of_hits.py.

USAGE
-----
    python3 hmm_finder.py --fasta my_seeds.faa
    python3 hmm_finder.py --fasta my_seeds.faa --out-dir results/ --iterations 3 --cpu 8

OUTPUT (under --out-dir, default: <fasta>_discovery/)
----------------------------------------------------------
    run1/ run2/ run3/        each: validated/{hits.tsv, hits_aa.faa,
                                    hits_nt.fna, hits_unique_aa.faa, ...}
    downstream/              clusters, clinker figures, homolog tree
    PACKAGE/                 labelled, self-contained deliverable
    pipeline.log             full run log
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

# Provenance records (run_manifest.json / METHODS.md) are meant to be shared, so
# they must never embed the user's NCBI e-mail. The address is not needed to
# reproduce a run (each user supplies their own), so we redact it at the source.
_EMAIL_RE = re.compile(r"[A-Za-z0-9._%+\-]+@[A-Za-z0-9.\-]+\.[A-Za-z]{2,}")


def _redact_emails(text):
    """Replace any e-mail address in a string with a placeholder (idempotent)."""
    return _EMAIL_RE.sub("<redacted-email>", text) if isinstance(text, str) else text


def _ignore_scratch(_dir, names):
    """copytree() ignore-fn for the 08_scripts reproducibility copy: drop bytecode
    and scratch/run helpers (gitignored as scripts/_*) — they can embed local
    paths or the user's NCBI e-mail — while keeping dunder files (e.g. __init__.py)."""
    drop = {n for n in names if n == "__pycache__" or n.endswith(".pyc")}
    drop |= {n for n in names if n.startswith("_") and not n.startswith("__")}
    return drop

# --- locate the deployable repo and the conda env, put tools on PATH ---------
HOME = Path.home()
HERE = Path(__file__).resolve().parent          # this scripts/ folder
# The search engine travels bundled inside the tool ( ../engine ). Fall back to
# the development repo if the bundle is absent, so this works on any machine.
_bundled_engine = HERE.parent / "engine"
DEPLOY = (_bundled_engine
          if (_bundled_engine / "scripts" / "run_all_database_benchmark.py").exists()
          else HOME / "Documents" / "HMM-Discovery-Deployable-20260602")
from env_paths import ensure_env_on_path  # noqa: E402  (sibling helper in scripts/)
from db_catalog import list_databases, pick_databases  # noqa: E402
ensure_env_on_path()

# Reuse the engine's validated convergence rule (hit count <5 % AND ΔLENG <3)
# rather than re-implementing it; the engine package imports with DEPLOY on path.
try:
    sys.path.insert(0, str(DEPLOY))
    from pipeline.iterative import convergence_check  # noqa: E402
except Exception:
    convergence_check = None
# Threshold calibration via built-in positive/negative controls (shuffled seeds
# always available; UniProt taxon negatives optional). Wired in below.
try:
    from pipeline.controls import run_all_controls, download_control_sequences  # noqa: E402
except Exception:
    run_all_controls = None
    download_control_sequences = None

BENCHMARK = DEPLOY / "scripts" / "run_all_database_benchmark.py"
EXTRACTOR = HERE / "extract_validated_hits.py"
CLUSTER = HERE / "cluster_and_clinker_corrected.py"
SYNTENY = HERE / "synteny_figure.py"
TREE = HERE / "build_tree_of_hits.py"
GENBANK = HERE / "build_real_genbanks.py"
ANNOTATE = HERE / "annotate_organism.py"


def write_gff3(hits_tsv, gff_path) -> None:
    """Write one CDS feature per validated hit (loadable in IGV/JBrowse/Artemis)."""
    import csv
    rows = list(csv.DictReader(open(hits_tsv), delimiter="\t"))
    with open(gff_path, "w") as f:
        f.write("##gff-version 3\n# one CDS feature per validated hit\n")
        for r in rows:
            org = r.get("organism", "").replace(";", ",")
            attrs = (f"ID={r['hit_id']};Name=family_homolog;organism={org};db={r['db_name']};"
                     f"evalue={r['evalue']};bit_score={r['bit_score']};"
                     f"domain_coverage={r['domain_coverage']};in_coding_locus={r['in_coding_locus']}")
            f.write(f"{r['genome_id']}\tHMM-Discovery\tCDS\t{r['nt_start']}\t{r['nt_end']}\t"
                    f"{r['bit_score']}\t{r['strand']}\t0\t{attrs}\n")

DATABASES = (
    "INPHARED genomes,INPHARED proteins,SwissProt,RefSeq viral proteins,"
    "RefSeq viral genomes,Gut Phage Database (GPD),GVD-AVrC,"
    "Pfam (sequences),Pfam (domain scan),VOGDB VFAM (annotation)"
)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    try:
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError:
        return ""


def _annotation_provenance(db_cache: Path) -> dict:
    """VOGDB annotation-DB provenance (version/URLs) if it was set up, else {}."""
    p = Path(db_cache).expanduser() / "annotation" / "vogdb" / "provenance.json"
    try:
        return json.loads(p.read_text()) if p.exists() else {}
    except Exception:
        return {}


def _git_commit(repo_dir: Path) -> str:
    try:
        r = subprocess.run(["git", "-C", str(repo_dir), "rev-parse", "--short", "HEAD"],
                           capture_output=True, text=True, timeout=10)
        return r.stdout.strip() if r.returncode == 0 else ""
    except Exception:
        return ""


def _tool_versions() -> dict:
    """Best-effort version strings by querying each tool directly (the engine's
    reproducibility.json leaves these blank). Picks the first version-looking line.
    """
    import re as _re
    probes = {
        "hmmsearch": ["hmmsearch", "-h"], "hmmbuild": ["hmmbuild", "-h"],
        "mafft": ["mafft", "--version"], "trimal": ["trimal", "--version"],
        "prodigal": ["prodigal", "-v"], "seqkit": ["seqkit", "version"],
        "cd-hit": ["cd-hit", "-h"], "iqtree": ["iqtree", "--version"],
        "meme": ["meme", "-version"], "fimo": ["fimo", "--version"],
        "mmseqs": ["mmseqs", "version"], "clinker": ["clinker", "--version"],
    }
    versions: dict = {}
    for name, cmd in probes.items():
        ver = ""
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
            text = (r.stdout or "") + "\n" + (r.stderr or "")
            lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
            # Prefer a line with a version-number token; else a line naming the tool.
            ver = (next((ln for ln in lines if _re.search(r"\d+\.\d+", ln)), "")
                   or next((ln for ln in lines if name.lower() in ln.lower()), "")
                   or (lines[0] if lines else ""))
        except Exception:
            ver = ""
        versions[name] = ver
    return versions


# Citations for the external tools and databases the pipeline uses. Keyed so the
# per-run METHODS lists only what that run actually invoked / searched.
_TOOL_CITATIONS = {
    "hmmsearch": "HMMER 3.4 — Eddy SR (2011) Accelerated profile HMM searches. PLoS Comput Biol 7:e1002195.",
    "hmmbuild":  "HMMER 3.4 — Eddy SR (2011) Accelerated profile HMM searches. PLoS Comput Biol 7:e1002195.",
    "hmmscan":   "HMMER 3.4 — Eddy SR (2011) Accelerated profile HMM searches. PLoS Comput Biol 7:e1002195.",
    "mafft":     "MAFFT v7 — Katoh K & Standley DM (2013) Mol Biol Evol 30:772-780.",
    "trimal":    "trimAl — Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T (2009) Bioinformatics 25:1972-1973.",
    "prodigal":  "Prodigal — Hyatt D et al. (2010) BMC Bioinformatics 11:119.",
    "seqkit":    "SeqKit — Shen W et al. (2016) PLoS One 11:e0163962.",
    "cd-hit":    "CD-HIT — Fu L et al. (2012) Bioinformatics 28:3150-3152.",
    "iqtree":    "IQ-TREE 2 — Minh BQ et al. (2020) Mol Biol Evol 37:1530-1534; ModelFinder — Kalyaanamoorthy S et al. (2017) Nat Methods 14:587-589; UFBoot2 — Hoang DT et al. (2018) Mol Biol Evol 35:518-522.",
    "meme":      "MEME Suite — Bailey TL et al. (2009) Nucleic Acids Res 37:W202-W208.",
    "fimo":      "FIMO — Grant CE, Bailey TL, Noble WS (2011) Bioinformatics 27:1017-1018.",
    "clinker":   "clinker — Gilchrist CLM & Chooi YH (2021) Bioinformatics 37:2473-2475.",
    "mmseqs":    "MMseqs2 — Steinegger M & Söding J (2017) Nat Biotechnol 35:1026-1028.",
}
_BIOPYTHON_CITE = "Biopython — Cock PJA et al. (2009) Bioinformatics 25:1422-1423."
_DB_CITATIONS = {
    "INPHARED":            "INPHARED — Cook R et al. (2021) PHAGE 2:214-223.",
    "RefSeq":              "NCBI RefSeq — O'Leary NA et al. (2016) Nucleic Acids Res 44:D733-D745.",
    "SwissProt":           "UniProtKB/Swiss-Prot — The UniProt Consortium (2023) Nucleic Acids Res 51:D523-D531.",
    "Gut Phage Database":  "Gut Phage Database (GPD) — Camarillo-Guerrero LF et al. (2021) Cell 184:1098-1109.",
    "GPD":                 "Gut Phage Database (GPD) — Camarillo-Guerrero LF et al. (2021) Cell 184:1098-1109.",
    "GVD":                 "Gut Virome Database (GVD) — Gregory AC et al. (2020) Cell Host Microbe 28:724-740.",
    "AVrC":                "Gut Virome Database (GVD) — Gregory AC et al. (2020) Cell Host Microbe 28:724-740.",
    "Pfam":                "Pfam — Mistry J et al. (2021) Nucleic Acids Res 49:D412-D419.",
    "VOGDB":               "VOGDB — Virus Orthologous Groups database (https://vogdb.org).",
    "PHROG":               "PHROGs — Terzian P et al. (2021) NAR Genom Bioinform 3:lqab067.",
}


def _citation_lines(tool_versions: dict, selected_dbs: str) -> list:
    """Markdown citation lines for the tools used + databases searched (de-duplicated)."""
    seen, refs = set(), []
    for t in (tool_versions or {}):
        c = _TOOL_CITATIONS.get(t)
        if c and c not in seen:
            seen.add(c); refs.append(c)
    if refs and _BIOPYTHON_CITE not in seen:
        refs.append(_BIOPYTHON_CITE)
    for key, c in _DB_CITATIONS.items():
        if key.lower() in (selected_dbs or "").lower() and c not in seen:
            seen.add(c); refs.append(c)
    if not refs:
        return []
    return (["", "## Citations",
             "Please cite this tool (see `CITATION.cff`) and the methods/databases used:"]
            + [f"- {r}" for r in refs])


def write_methods_log(out: Path, args, fasta: Path, label: str, selected_dbs: str,
                      iter_hits: list, started_at: str, log, stop_reason: str = "",
                      control_summary: dict | None = None,
                      seed_recovery: dict | None = None,
                      interrupted: dict | None = None) -> None:
    """Write a consolidated methodology record at the run root: run_manifest.json
    (machine-readable) + METHODS.md (human-readable). Aggregates tool versions and
    per-database provenance (source URLs, SHA256, access dates) from each
    iteration's engine-generated reproducibility.json. Never raises.
    """
    try:
        finished_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        tool_versions = _tool_versions()
        db_provenance: list = []
        for i, _ in iter_hits:
            repro = out / f"run{i}" / "benchmark" / "reports" / "reproducibility.json"
            if not repro.exists():
                continue
            try:
                data = json.loads(repro.read_text())
            except Exception:
                continue
            for entry in data.get("database_provenance", []):
                db_provenance.append({"run": i, **entry})

        import shlex
        try:
            _nseed = sum(1 for ln in Path(args.fasta).read_text().splitlines()
                         if ln.startswith(">")) if getattr(args, "fasta", None) else None
        except Exception:
            _nseed = None
        manifest = {
            "tool": "hmm-homologue-finder",
            "code_git_commit": _git_commit(HERE),
            "started_at": started_at,
            "finished_at": finished_at,
            # shlex-quoted so paths containing spaces can be pasted and re-run verbatim;
            # the e-mail is redacted (provenance is shared; not needed to reproduce).
            "command_line": _redact_emails(" ".join(shlex.quote(a) for a in sys.argv)),
            "n_input_seeds": _nseed,
            "conda_env": os.environ.get("CONDA_DEFAULT_ENV", ""),
            "conda_prefix": os.environ.get("CONDA_PREFIX", ""),
            "python": sys.version.split()[0],
            "parameters": {
                "label": label, "iterations": args.iterations, "cpu": args.cpu,
                "databases": selected_dbs, "prodigal_gate": bool(args.prodigal_gate),
                "min_recovery": "0.70", "max_synteny_genomes": "200",
                # record only WHETHER an e-mail (online annotation) was used, never the address
                "email": ("<redacted>" if args.email else ""),
                "db_cache": str(args.db_cache), "out_dir": str(out),
                "input_type": args.input_type, "trans_table": args.trans_table,
                "no_annotate": bool(args.no_annotate),
            },
            "annotation_database": _annotation_provenance(args.db_cache),
            "input": {"fasta": str(fasta), "sha256": _sha256(fasta)},
            "per_iteration_unique_seeds": [{"run": i, "unique_validated_seeds": n} for i, n in iter_hits],
            "iteration_stop_reason": stop_reason,
            "threshold_calibration": control_summary or {},
            "seed_recovery_qc": seed_recovery or {},
            "interrupted_homologs": interrupted or {},
            "tool_versions": tool_versions,
            "database_provenance": db_provenance,
        }
        (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2))

        L = [f"# Methods — {label}\n",
             f"- **Tool:** hmm-homologue-finder (commit `{manifest['code_git_commit'] or 'unknown'}`)",
             f"- **Run:** {started_at} → {finished_at}",
             f"- **Command:** `{manifest['command_line']}`",
             f"- **Conda env:** {manifest['conda_env']}  |  Python {manifest['python']}",
             f"- **Input FASTA:** `{fasta}`  (sha256 `{manifest['input']['sha256'][:16]}…`)",
             "", "## Parameters"]
        L += [f"- {k}: {v}" for k, v in manifest["parameters"].items()]
        L += ["", "## Iterations"]
        L += [f"- run{i}: {n} unique validated seeds" for i, n in iter_hits]
        if stop_reason:
            L += [f"- **Stopping criterion:** {stop_reason}"]
        if control_summary:
            cs = control_summary
            L += ["", "## Threshold calibration (built-in controls)",
                  f"- At the strict bit-score threshold (45): sensitivity "
                  f"{cs.get('sensitivity')} ({cs.get('true_positives')}/{cs.get('total_positives')} "
                  f"seed sequences recovered), specificity {cs.get('specificity')}, "
                  f"false-positive rate {cs.get('false_positive_rate')} "
                  f"({cs.get('false_positives')}/{cs.get('total_negatives')} negative-control "
                  f"sequences scored ≥45).",
                  f"- Negative controls: composition-matched shuffled seeds"
                  + (" + unrelated-proteome sets" if (cs.get('n_negative_controls') or 0) > 1 else "")
                  + f" ({cs.get('n_negative_controls')} set(s)). Full detail: `controls/control_report.json`."]
            _roc = cs.get("roc") or {}
            if _roc:
                L += [f"- **ROC calibration (advisory):** AUC {_roc.get('auc')} over "
                      f"{_roc.get('n_positive')} positive vs {_roc.get('n_negative')} negative "
                      f"control sequences. The Youden's-J optimal bit-score cutoff is "
                      f"{_roc.get('optimal_threshold')} (sensitivity {_roc.get('sensitivity_at_optimum')}, "
                      f"specificity {_roc.get('specificity_at_optimum')} there); the pipeline's fixed "
                      f"strict threshold (45) is retained for tiering. See `controls/roc_curve.svg`."]
        if seed_recovery:
            sr = seed_recovery
            miss = sr.get("not_recovered_after") or []
            L += ["", "## Seed-recovery QC (per seed; `seed_qc/seed_recovery.csv`)",
                  f"- Of {sr.get('n_seeds')} input seeds, {sr.get('recovered_before')} are recovered by "
                  f"the initial model and {sr.get('recovered_after')} by the final refined model "
                  f"(strict bit≥{sr.get('strict_bits')}).",
                  ("- Every input seed is recovered by the final model."
                   if not miss else
                   f"- Not recovered by the final model ({len(miss)}): "
                   + ", ".join(miss[:12]) + (" …" if len(miss) > 12 else "")
                   + " — likely divergent outliers; see the CSV for per-seed scores.")]
        if interrupted and interrupted.get("interrupted_candidates") is not None:
            ic = interrupted
            L += ["", "## Stop-interrupted / overprinted homologs (`interrupted_homologs.tsv`)",
                  f"- A read-through translation of the searched nucleotide databases (stop codons "
                  f"retained, not broken on) was scanned with the family HMM. {ic.get('interrupted_candidates')} "
                  f"match(es) carry ≥1 internal stop within the domain — candidate homologs interrupted by "
                  f"a premature stop (e.g. overprinted genes whose nonsense mutation is silent in an "
                  f"overlapping reading frame). The stop-to-stop search cannot recover these."]
            if ic.get("threshold_basis"):
                L += [f"- Reporting threshold: {ic.get('min_bit')} bits "
                      f"[{ic.get('threshold_basis')}] — never looser than the run's own "
                      f"evidence bar, raised to the family-calibrated noise floor."]
            if ic.get("domain_faa"):
                L += ["- Protein sequences emitted with the internal stop(s) shown as `*`: "
                      "`interrupted_homologs_domain_aa.faa` (matched domain) and "
                      "`interrupted_homologs_full_orf_aa.faa` (full read-through ORF, "
                      "premature stops kept, terminal `*` = the natural gene end).",
                      "- Full read-through ORF **nucleotide** in `interrupted_homologs_full_orf_nt.fna` "
                      "(coding 5'→3', ending in the actual stop codon triplet; translates back to "
                      "`full_orf_aa`). TSV columns `orf_nt_start/end`, `natural_stop_nt` (genome "
                      "coordinate of the actual stop codon), and `full_orf_nt` carry the same."]
            if ic.get("overprinting_strong") is not None:
                L += ["- **Overprinting (silent-stop) test.** For each premature stop the scan picks "
                      "the antisense frame with the fewest stops across the domain (the candidate "
                      "overprinted ORF) and tests whether the stop is **synonymous** in that frame — "
                      "i.e. whether reverting the nonsense mutation in this gene would leave the "
                      "antisense protein unchanged. TSV columns: `overprinting_support` "
                      "(strong/partial/none), `antisense_open_frame`, `antisense_open_stops` "
                      "(0 = a fully open overlapping ORF), `stop_silent_antisense` (per-stop 1/0). "
                      f"Result: {ic.get('overprinting_strong')} strong (stop synonymous in an open "
                      f"antisense ORF) and {ic.get('overprinting_partial')} partial. The discriminating "
                      "signal is the antisense frame being OPEN across the whole domain (improbable by "
                      "chance for a long domain); a single stop being synonymous in that frame is, alone, "
                      "expected ~85-100% of the time from genetic-code geometry, so it is weak by itself. "
                      "`strong` is a necessary-but-not-sufficient SEQUENCE signature of antisense "
                      "overprinting (open frame + synonymy), NOT proof the antisense ORF is expressed."]
        if tool_versions:
            L += ["", "## Tool versions"]
            for t, info in sorted(tool_versions.items()):
                ver = info.get("version", "") if isinstance(info, dict) else str(info)
                L.append(f"- {t}: {ver}")
        if db_provenance:
            L += ["", "## Databases searched (provenance)"]
            for e in db_provenance:
                L.append(f"- **{e.get('database','?')}** (run{e.get('run','?')}, {e.get('type','?')}, "
                         f"status={e.get('status','?')}, hits={e.get('hit_count','?')})")
                if e.get("source_accessed_first"):
                    L.append(f"    - accessed: {e.get('source_accessed_first')}")
                for u in (e.get("source_urls") or [])[:8]:
                    L.append(f"    - url: {u}")
                for s in (e.get("source_sha256s") or [])[:8]:
                    L.append(f"    - sha256: {s}")
        L += _citation_lines(tool_versions, selected_dbs)
        L += ["", "> Full machine-readable provenance: `run_manifest.json` (this folder) and "
              "`run*/benchmark/reports/reproducibility.json`."]
        (out / "METHODS.md").write_text("\n".join(L) + "\n")
        log(f"Methodology log written: {out / 'METHODS.md'}  +  {out / 'run_manifest.json'}")
    except Exception as e:
        log(f"  (methods log skipped: {e})")


def _hmm_leng(hmm_path) -> int:
    """Read the LENG (match-state count) from an HMMER profile; 0 if unavailable."""
    try:
        for ln in Path(hmm_path).read_text(errors="replace").splitlines():
            if ln.startswith("LENG"):
                return int(ln.split()[1])
            if ln.startswith("HMM "):  # header reached; no LENG before it
                break
    except Exception:
        pass
    return 0


def _best_run_index(out: Path, iterations: int) -> int:
    """Pick the most complete run (most validated hit rows) as the canonical set
    for figures + the main paper table. After convergence this is the refined
    final round; ties resolve to the earliest run for determinism. Matches
    export_csv.py's `best_run = max(run_frames, key=len)`. Falls back to 1."""
    best_i, best_n = 1, -1
    for i in range(1, iterations + 1):
        tsv = out / f"run{i}" / "benchmark" / "validated" / "hits.tsv"
        if not tsv.exists():
            continue
        try:
            n = max(0, sum(1 for _ in tsv.open("r", errors="replace")) - 1)  # rows minus header
        except OSError:
            continue
        if n > best_n:
            best_i, best_n = i, n
    return best_i


def _looks_like_nucleotide(fasta: Path) -> bool:
    """Heuristic: >=90% of sequence characters are ACGTUN -> nucleotide FASTA."""
    nt = total = 0
    with fasta.open("r", errors="replace") as fh:
        for ln in fh:
            if ln.startswith(">") or not ln.strip():
                continue
            for ch in ln.strip().upper():
                total += 1
                if ch in "ACGTUN":
                    nt += 1
            if total >= 4000:
                break
    return total > 0 and nt / total >= 0.9


def translate_seed(fasta: Path, table: int, out_dir: Path, log) -> Path:
    """Translate a nucleotide CDS FASTA into a protein seed using genetic code
    `table` (frame 0; each record assumed to be a coding sequence). Returns the
    protein FASTA path. Flags any internal stop codons (wrong code / not a CDS)."""
    import warnings
    from Bio import BiopythonWarning, SeqIO
    from Bio.Seq import Seq
    warnings.simplefilter("ignore", BiopythonWarning)
    out_faa = out_dir / f"{fasta.stem}_translated_tt{table}.faa"
    n = 0
    flagged = []
    with out_faa.open("w") as oh:
        for rec in SeqIO.parse(str(fasta), "fasta"):
            s = str(rec.seq).upper().replace("U", "T")
            usable = len(s) - (len(s) % 3)
            if usable < 3:
                continue
            aa_full = str(Seq(s[:usable]).translate(table=table))
            internal = aa_full[:-1].count("*") if aa_full.endswith("*") else aa_full.count("*")
            aa = aa_full.rstrip("*")
            if internal:
                flagged.append(rec.id)
            oh.write(f">{rec.id} {rec.description}\n{aa}\n")
            n += 1
    log(f"Translated {n} nucleotide sequence(s) with genetic code {table} -> {out_faa.name}")
    if flagged:
        log(f"  WARNING: {len(flagged)} sequence(s) had internal stop codons — likely the "
            f"wrong genetic code or not a clean CDS: {', '.join(flagged[:5])}"
            + (" …" if len(flagged) > 5 else ""))
    return out_faa


def write_csv_exports(out: Path, log) -> None:
    """Write Excel-friendly CSV copies + merged tables. Never raises."""
    try:
        from export_csv import export as _export_csv
        _export_csv(out)
        log(f"CSV exports written: {out / 'all_runs_hits.csv'}, "
            f"{out / 'hit_summary.csv'}, {out / 'database_summary.csv'}")
    except Exception as e:
        log(f"  (CSV export skipped: {e})")


def write_report(out: Path, log) -> None:
    """Write the one-page HTML summary report. Never raises."""
    try:
        from generate_report import generate as _gen_report
        p = _gen_report(out)
        log(f"Summary report written: {p}")
    except Exception as e:
        log(f"  (report skipped: {e})")


def run_controls(hmm_path: Path, seed_faa: Path, out: Path, mode: str,
                 strict: float, moderate: float, cpu, log) -> dict:
    """Calibrate the bit-score thresholds with built-in controls: sensitivity on
    the seed self-test (positive) and false-positive rate on composition-matched
    shuffled seeds plus any available unrelated-proteome negatives. Writes
    controls/control_report.json + controls_summary.csv and returns the summary
    dict (empty on failure). Never raises."""
    if run_all_controls is None:
        log("  (controls skipped: controls module unavailable)")
        return {}
    if not Path(hmm_path).exists():
        log(f"  (controls skipped: HMM not found at {hmm_path})")
        return {}
    try:
        cdir = out / "controls"
        rep = run_all_controls(hmm_path=Path(hmm_path), seed_faa=Path(seed_faa),
                               out_dir=cdir, mode=mode, strict_threshold=strict,
                               moderate_threshold=moderate, cpu=int(cpu))
        summary = rep.summary()
        cdir.mkdir(parents=True, exist_ok=True)
        (cdir / "control_report.json").write_text(rep.to_json())
        try:
            rep.to_dataframe().to_csv(cdir / "controls_summary.csv", index=False)
        except Exception:
            pass
        try:
            rep.plot_roc(cdir)  # roc_curve.{png,svg,pdf} — advisory calibration figure
        except Exception:
            pass
        log(f"  Controls (strict bit≥{strict}): sensitivity {summary.get('sensitivity')}, "
            f"specificity {summary.get('specificity')}, FPR {summary.get('false_positive_rate')} "
            f"({summary.get('true_positives',0)}/{summary.get('total_positives',0)} seeds recovered; "
            f"{summary.get('false_positives',0)}/{summary.get('total_negatives',0)} negatives passed)")
        _roc = summary.get("roc") or {}
        if _roc:
            log(f"  Calibration ROC: AUC {_roc.get('auc')}; Youden-optimal bit "
                f"{_roc.get('optimal_threshold')} (strict threshold in use: {strict:g})")
        return summary
    except Exception as e:
        log(f"  (controls skipped: {e})")
        return {}


def _interrupted_min_bit(control_summary: dict | None, floor: float = 30.0):
    """Reporting threshold for the read-through interrupted scan: never below `floor`
    (the run's lenient evidence bar), raised to the family ROC-Youden optimal cutoff
    when controls were run. Returns (min_bit, basis_string). The read-through scan
    covers a much larger, noisier space than the stop-to-stop search, so the bar only
    ever tightens; with no ROC (--no-controls) the bare floor is used."""
    roc = (control_summary or {}).get("roc") or {}
    try:
        opt = float(roc.get("optimal_threshold"))
        return max(floor, opt), f"max(floor {floor:g}, ROC-Youden {opt:g})"
    except (TypeError, ValueError):
        return floor, f"floor {floor:g} bits (no ROC calibration available)"


def run_find_interrupted(out: Path, hmm: Path, db_cache: Path, databases: str,
                         cpu, log, control_summary: dict | None = None) -> dict:
    """Read-through scan of the searched NUCLEOTIDE databases for homologs that are
    interrupted by a premature stop codon (e.g. overprinted genes). Writes
    out/interrupted_homologs.tsv and returns a summary. Never raises.

    Reporting threshold is family-calibrated, never looser than the run's own
    evidence bar: the read-through scan covers a much larger, noisier space (every
    frame, stops retained, windowed) than the stop-to-stop search, so we floor the
    bit score at the validation lenient bound (30) and RAISE it to the family's ROC
    Youden-optimal cutoff when controls were run — only ever tightening. With no
    calibration (--no-controls) it falls back to the floor."""
    try:
        from find_interrupted import _run as _fi  # noqa: E402  (sibling)
    except Exception as e:
        log(f"  (find-interrupted skipped: {e})")
        return {}
    if not Path(hmm).exists():
        log("  (find-interrupted skipped: HMM not found)")
        return {}
    cache = Path(db_cache).expanduser() / "cache"
    if not cache.exists():
        log("  (find-interrupted skipped: no database cache)")
        return {}
    dbl = (databases or "").lower()
    targets: list = []
    for sub in sorted(p for p in cache.iterdir() if p.is_dir()):
        if sub.name.replace("_", " ").lower() not in dbl:
            continue
        if not any(sub.glob("*.sixframe.*")):   # only nucleotide DBs were six-frame-translated
            continue
        for pat in ("*.fa.gz", "*.fna.gz", "*.fasta.gz", "*.fa", "*.fna", "*.fasta"):
            targets += [f for f in sorted(sub.glob(pat)) if ".sixframe." not in f.name]
    if not targets:
        log("  (find-interrupted: no cached nucleotide DBs from this run to scan)")
        return {}
    # Family-calibrated reporting floor (see docstring): max(30, ROC-Youden).
    min_bit, thr_basis = _interrupted_min_bit(control_summary)
    log(f"Read-through scan for interrupted/overprinted homologs "
        f"({len(targets)} nucleotide DB file(s)); reporting threshold {min_bit:g} bits "
        f"[{thr_basis}]…")
    import csv as _csv
    out_tsv = out / "interrupted_homologs.tsv"
    all_rows, scored = [], 0
    for fa in targets:
        tmp = out / f".interrupted_{fa.stem}.part.tsv"
        try:
            # emit_fasta=False: per-DB temp runs only produce partial TSVs; the
            # protein FASTAs are written once below, from the aggregated rows.
            s = _fi(fa, Path(hmm), tmp, min_bit, int(cpu), log, emit_fasta=False)
            scored += s.get("matches_scored", 0)
            if tmp.exists():
                all_rows += list(_csv.DictReader(open(tmp), delimiter="\t"))
                tmp.unlink()
        except Exception as e:
            log(f"  (find-interrupted: {fa.name} skipped: {e})")
    domain_faa = orf_faa = orf_fna = ""
    if all_rows:
        all_rows.sort(key=lambda r: -float(r.get("domain_bit_score", 0) or 0))
        with open(out_tsv, "w", newline="") as f:
            w = _csv.DictWriter(f, fieldnames=list(all_rows[0].keys()), delimiter="\t")
            w.writeheader(); w.writerows(all_rows)
        # Protein (domain + full ORF, internal stop '*') and the full-ORF nucleotide
        # (incl. the actual stop codon triplet).
        from find_interrupted import write_aa_fastas as _wfa, write_orf_nt_fasta as _wnt  # noqa: E402
        _dom, _orf = _wfa(all_rows, out_tsv)
        _nt = _wnt(all_rows, out_tsv)
        domain_faa, orf_faa, orf_fna = str(_dom), str(_orf), str(_nt)
        log(f"  interrupted-homolog sequences -> {_dom.name}, {_orf.name}, {_nt.name}")
    # Overprinting (silent-stop) evidence tallied across the candidates.
    strong = sum(1 for r in all_rows if r.get("overprinting_support") == "strong")
    partial = sum(1 for r in all_rows if r.get("overprinting_support") == "partial")
    summary = {"matches_scored": scored, "interrupted_candidates": len(all_rows),
               "min_bit": min_bit, "threshold_basis": thr_basis,
               "tsv": str(out_tsv) if all_rows else "",
               "domain_faa": domain_faa, "orf_faa": orf_faa, "orf_fna": orf_fna,
               "overprinting_strong": strong, "overprinting_partial": partial}
    log(f"  interrupted-homolog scan: {len(all_rows)} candidate(s) carrying an internal stop"
        + (f" -> {out_tsv.name}" if all_rows else " (none found)"))
    if all_rows:
        log(f"  overprinting (silent-stop) evidence: {strong} strong "
            f"(stop synonymous in an open antisense ORF), {partial} partial")
    return summary


def run_perhit_hmm_alignment(hmm: Path, faa: Path, out_dir: Path, log) -> dict:
    """Align every unique homolog to the family HMM (hmmalign) so each hit's match
    to the model is explicit — match states vs insertions, not just the all-vs-all
    MSA. Writes hits_hmmalign.sto (Stockholm) + hits_hmmalign.a2m (A2M). Non-fatal."""
    out_dir = Path(out_dir)
    if not (Path(hmm).exists() and Path(faa).exists()):
        return {}
    out_dir.mkdir(parents=True, exist_ok=True)
    sto = out_dir / "hits_hmmalign.sto"
    a2m = out_dir / "hits_hmmalign.a2m"
    try:
        subprocess.run(["hmmalign", "--amino", "--trim", "-o", str(sto), str(hmm), str(faa)],
                       check=True, capture_output=True, text=True)
        try:
            subprocess.run(["esl-reformat", "-o", str(a2m), "a2m", str(sto)],
                           check=True, capture_output=True, text=True)
        except Exception as e:
            log(f"  (per-hit HMM alignment: A2M reformat skipped: {e})")
            a2m = None
        log(f"  per-hit HMM alignment -> {sto.name}" + (f", {a2m.name}" if a2m else ""))
        return {"sto": str(sto), "a2m": str(a2m) if a2m else ""}
    except Exception as e:
        log(f"  (per-hit HMM alignment skipped: {e})")
        return {}


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", type=Path, default=None, help="seed protein FASTA (if omitted, you'll be prompted)")
    ap.add_argument("--out-dir", type=Path, default=None, help="output root (default: <fasta>_discovery)")
    ap.add_argument("--iterations", type=int, default=3, help="number of search iterations (default 3)")
    ap.add_argument("--cpu", default="8")
    ap.add_argument("--email", default=None,
                    help="NCBI Entrez email for organism/sequence lookups. If omitted, "
                         "the NCBI_EMAIL environment variable is used; if neither is set "
                         "the run proceeds fully offline (equivalent to --no-annotate). "
                         "Never hardcoded — NCBI requires a real address.")
    ap.add_argument("--input-type", choices=("auto", "protein", "nucleotide"), default="auto",
                    help="seed FASTA type; 'auto' detects nucleotide vs protein (default auto)")
    ap.add_argument("--trans-table", type=int, default=11,
                    help="genetic code to translate a NUCLEOTIDE seed (default 11 = "
                         "bacterial/archaeal/phage; e.g. 4 = Mycoplasma/Spiroplasma)")
    ap.add_argument("--no-annotate", action="store_true",
                    help="skip the NCBI organism-name lookup (for fully offline / "
                         "fastest unattended runs); hit tables still build")
    ap.add_argument("--name", default=None,
                    help="label for the output folder (default: derived from the FASTA name)")
    ap.add_argument("--databases", default=None,
                    help="comma-separated databases to search (must match catalog names; see "
                         "--list-databases). If omitted: an interactive terminal PROMPTS you to "
                         "choose; a non-interactive/unattended run uses the full default set. "
                         "Use --all-databases to take the full set without being asked.")
    ap.add_argument("--all-databases", action="store_true",
                    help="search the full default database set without the interactive prompt "
                         "(use in scripts where you want everything but still have a TTY)")
    ap.add_argument("--no-seed-tree", action="store_true",
                    help="skip the pre-run seed QC tree + seed alignment (a quick "
                         "phylogeny/alignment of just the input seeds, for sanity-checking "
                         "the seed set before the discovery runs). The final homolog tree "
                         "always includes the seeds, marked.")
    ap.add_argument("--synteny-gene-labels", action="store_true",
                    help="label neighbour genes with their functional annotation in the synteny "
                         "figures (off by default; an interactive run will ask). Can clutter "
                         "dense neighbourhoods — colours + legend already convey function.")
    ap.add_argument("--color-by", choices=("function", "conservation", "both"), default="both",
                    help="how to colour the synteny neighbourhood genes: by functional category "
                         "(VOGDB), by cross-locus conservation, or both (default both)")
    ap.add_argument("--no-controls", action="store_true",
                    help="skip the threshold-calibration controls (sensitivity on the "
                         "seed self-test + false-positive rate on shuffled/unrelated "
                         "negatives). Controls add a few seconds and need no network "
                         "for the shuffled negative.")
    ap.add_argument("--biology-mode", choices=("generic", "phage", "bacterial"),
                    default="phage",
                    help="control panel to calibrate against (default phage, matching the "
                         "default phage/viral databases)")
    ap.add_argument("--download-controls", action="store_true",
                    help="one-time: fetch the UniProt unrelated-proteome negative controls "
                         "(fungi/mammalian/archaea) into the bundle, then continue. Without "
                         "this only the always-available shuffled-seeds negative is used.")
    ap.add_argument("--smoke", action="store_true",
                    help="fast self-test: 1 iteration against a single small database")
    ap.add_argument("--skip-tool-check", action="store_true", help="skip the startup software check")
    ap.add_argument("--prodigal-gate", action="store_true",
                    help="require six-frame hits to overlap a Prodigal-predicted gene to "
                         "pass (stricter / higher specificity). Default off: keeps genuine "
                         "antisense/alternate-frame homologs the tool is designed to find.")
    ap.add_argument("--find-interrupted", action="store_true",
                    help="ALSO scan the searched nucleotide databases with READ-THROUGH "
                         "translation (stops kept, not broken on) to find homologs interrupted "
                         "by a premature stop codon — e.g. overprinted genes where a nonsense "
                         "mutation in this gene is silent in an overlapping gene. The normal "
                         "stop-to-stop search misses these. Writes interrupted_homologs.tsv. "
                         "Opt-in (extra scan of each nucleotide DB).")
    ap.add_argument("--list-databases", action="store_true",
                    help="print the available search databases (with sizes/times) and exit")
    ap.add_argument("--pick-databases", action="store_true",
                    help="interactively choose which databases to search")
    ap.add_argument("--db-cache", type=Path,
                    default=Path.home() / ".cache" / "hmm-homologue-finder",
                    help="persistent cache dir for downloaded databases, shared across ALL "
                         "runs (default: ~/.cache/hmm-homologue-finder) — each DB downloads "
                         "once, ever, which makes repeat runs much faster")
    args = ap.parse_args()

    # Clamp thread count to available cores. Over-allocating is a real autonomy
    # hazard: IQ-TREE aborts with "more threads than CPU cores available" (exit 2)
    # and MEME/cd-hit/hmmsearch oversubscribe. This only ever REDUCES toward the
    # core count — high-core machines (e.g. an M3 Pro, or a 32-core server) keep
    # the full requested --cpu; the default 8 is untouched on any ≥8-core host.
    try:
        _avail = os.cpu_count() or 1
        _req = max(1, int(args.cpu))
        if _req > _avail:
            print(f"(--cpu {_req} exceeds {_avail} available cores; clamping to {_avail})")
        args.cpu = str(min(_req, _avail))
    except (TypeError, ValueError):
        args.cpu = str(os.cpu_count() or 4)

    # --list-databases: show the catalog and exit (no FASTA needed).
    if args.list_databases:
        list_databases(DEPLOY)
        return

    started_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Fail fast with a clear message if the search engine is missing, rather than a
    # cryptic subprocess error mid-run. Happens on a partial clone (no engine/) when
    # the dated dev-repo fallback path also isn't present.
    if not BENCHMARK.exists():
        sys.exit(f"Search engine not found at {BENCHMARK}.\n"
                 "The bundled engine/ folder is missing from this install — re-clone or "
                 "restore the full repository (it must include engine/).")

    # --smoke: minutes-long install/sanity check on a single fast database.
    if args.smoke:
        args.iterations = 1
        args.databases = "INPHARED proteins"

    # Resolve which databases to search. Precedence:
    #   explicit --databases  >  --all-databases / --pick-databases  >  default.
    # DEFAULT BEHAVIOUR: on an interactive terminal, PROMPT for the selection;
    # unattended (cron/nohup/pipe) silently uses the full default set so it never
    # blocks. --all-databases forces the full set even on a TTY.
    if not args.smoke and args.databases is None:
        if args.all_databases:
            args.databases = DATABASES
        elif sys.stdin.isatty():
            args.databases = pick_databases(DEPLOY, DATABASES.split(","))
            print(f"Selected databases: {args.databases}")
        else:
            args.databases = DATABASES
            print("(no --databases given and not interactive; searching the full default set)")
    elif not args.smoke and args.pick_databases and sys.stdin.isatty():
        # explicit --pick-databases still forces the prompt even if --databases was set
        args.databases = pick_databases(DEPLOY, DATABASES.split(","))
        print(f"Selected databases: {args.databases}")

    # (Synteny gene-name labels are opt-in via --synteny-gene-labels only — off by
    # default and not prompted, since they overlap on dense neighbourhoods; the
    # functional colours + legend already convey gene function.)

    # Preflight: refuse to start a multi-hour run if required software is missing.
    if not args.skip_tool_check:
        try:
            import check_tools
            if not check_tools.ensure(install=False):
                print("\nRequired software is missing. Install it with:")
                print(f"    bash {Path(__file__).resolve().parent.parent / 'setup.sh'}")
                print("…or re-run with --skip-tool-check to override.")
                sys.exit(1)
        except ImportError:
            pass  # check_tools not alongside; continue (tools may still be present)

    # Seed FASTA. Prompt ONLY when attached to a terminal; in an unattended run
    # (cron, nohup, background) a missing --fasta must fail fast, never block on input().
    if args.fasta is None:
        if not sys.stdin.isatty():
            sys.exit("No --fasta provided and not running interactively. "
                     "Pass --fasta <seed.faa|seed.fna> for unattended runs.")
        print("\n=== HMM-based homolog discovery pipeline ===")
        print("Drag your seed FASTA here (or type its path) and press Enter:")
        raw = input("  seed FASTA > ").strip().strip("'\"").strip()
        if not raw:
            sys.exit("No FASTA provided.")
        args.fasta = Path(raw)

    fasta = args.fasta.expanduser().resolve()
    if not fasta.exists():
        sys.exit(f"Seed FASTA not found: {fasta}")
    # Validate it really is a non-empty FASTA before committing to a long run:
    # the first non-blank line must start with '>' and there must be ≥1 record.
    try:
        n_records = 0
        first_content_line = None
        with fasta.open("r", errors="replace") as _fh:
            for _ln in _fh:
                if not _ln.strip():
                    continue
                if first_content_line is None:
                    first_content_line = _ln
                if _ln.startswith(">"):
                    n_records += 1
    except OSError as e:
        sys.exit(f"Could not read seed FASTA {fasta}: {e}")
    if first_content_line is None:
        sys.exit(f"Seed FASTA is empty: {fasta}")
    if not first_content_line.startswith(">") or n_records == 0:
        sys.exit(f"Seed FASTA does not look like FASTA (no '>' header found): {fasta}")
    label = args.name or fasta.stem
    out = (args.out_dir or fasta.parent / f"{label}_discovery").resolve()
    # The benchmark refuses to write inside the deployable repo. If the chosen
    # output would land there (e.g. the seed FASTA lives in the repo), redirect
    # to ~/Documents so the run can proceed.
    if DEPLOY == out or DEPLOY in out.parents:
        out = (Path.home() / "Documents" / f"{label}_discovery").resolve()
        print(f"(output redirected outside the deployable repo: {out})")
    out.mkdir(parents=True, exist_ok=True)
    log_path = out / "pipeline.log"

    def log(msg: str) -> None:
        line = f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}"
        print(line, flush=True)
        with log_path.open("a") as fh:
            fh.write(line + "\n")

    def sh(cmd: list[str], **kw) -> None:
        log("  $ " + " ".join(str(c) for c in cmd))
        subprocess.run(cmd, check=True, **kw)

    # Resolve the NCBI Entrez email WITHOUT ever assuming one (NCBI policy + project
    # rule). Precedence: --email > $NCBI_EMAIL > interactive prompt (TTY only) >
    # offline. Resolved here — once a real run is committed — so --list-databases and
    # arg errors never prompt. No email => force --no-annotate so the run never sends
    # a fake address or hangs on input(); local six-frame hits are unaffected, only
    # NCBI lookups (organism names, protein-DB hit sequences) are skipped.
    if not args.email:
        args.email = (os.environ.get("NCBI_EMAIL") or "").strip() or None
    if not args.email and sys.stdin.isatty():
        try:
            raw = input("NCBI Entrez email (press Enter to run offline, no organism lookups): ").strip()
        except EOFError:
            raw = ""
        args.email = raw or None
    if not args.email:
        if not args.no_annotate:
            log("No NCBI email provided (--email / $NCBI_EMAIL / prompt all empty); "
                "running offline — organism & protein-DB-sequence lookups skipped.")
        args.no_annotate = True
    else:
        log(f"NCBI Entrez email: {args.email}")

    # If the seed FASTA is nucleotide, translate it into a protein seed first.
    is_nt = args.input_type == "nucleotide" or (
        args.input_type == "auto" and _looks_like_nucleotide(fasta))
    if is_nt:
        log(f"Seed looks like nucleotide DNA; translating to protein with genetic code "
            f"{args.trans_table} (use --trans-table to change).")
        fasta = translate_seed(fasta, args.trans_table, out, log)
        if not fasta.exists() or fasta.stat().st_size == 0:
            sys.exit("Translation produced no protein sequences — check the input or --trans-table.")

    # Persistent global cache shared across ALL runs: download each DB once, ever.
    # (Override per-run with --db-cache; defaults to ~/.cache/hmm-homologue-finder.)
    shared = args.db_cache.expanduser().resolve()
    (shared / "cache").mkdir(parents=True, exist_ok=True)
    (shared / "db_setup").mkdir(parents=True, exist_ok=True)

    log(f"=== family pipeline: {args.iterations} iterations from {fasta.name} ===")
    log(f"Shared database cache: {shared}")
    log(f"Databases: {args.databases}")

    # Pre-run seed QC: a quick phylogeny + alignment of just the input seeds, so a
    # mis-curated/outlier seed is visible before committing to the discovery runs.
    # Skipped in smoke mode and with --no-seed-tree. Non-fatal.
    if not args.no_seed_tree and not args.smoke:
        log("Pre-run: seed QC tree + alignment (input seeds only)")
        try:
            sh(["python3", str(TREE), "--faa", str(fasta),
                "--out-dir", str(out / "seed_qc"), "--cpu", args.cpu,
                "--mafft-mode", "accurate"])
        except Exception as e:
            log(f"  (seed QC tree skipped: {e})")

    iter_hits: list = []
    seed = fasta
    prev_n = None          # previous round's unique-validated-seed count
    prev_leng = None       # previous round's HMM length (match states)
    stop_reason = f"ran the requested {args.iterations} iteration(s)"
    hmm_of = lambda b: b / "hmm" / "benchmark_profile.hmm"
    for i in range(1, args.iterations + 1):
        run_dir = out / f"run{i}"
        bench = run_dir / "benchmark"
        validated = bench / "validated"
        unique = validated / "hits_unique_aa.faa"

        if unique.exists() and unique.stat().st_size > 0:
            n_done = sum(1 for ln in unique.read_text().splitlines() if ln.startswith(">"))
            iter_hits.append((i, n_done))
            log(f"RUN {i}: already complete; skipping.")
            prev_n, prev_leng = n_done, _hmm_leng(hmm_of(bench))
            seed = unique
            continue

        log(f"RUN {i}: searching all databases (seed = {seed.name})")
        bench.mkdir(parents=True, exist_ok=True)
        # reuse shared cache + db_setup via symlinks
        for name in ("cache", "db_setup"):
            link = bench / name
            if not link.exists():
                link.symlink_to(shared / name)

        # 1. search (builds HMM internally, six-frame translates nucleotide DBs)
        try:
            sh(["python3", str(BENCHMARK), "--fasta", str(seed), "--out", str(bench),
                "--databases", args.databases, "--cpu", args.cpu, "--keep-cache",
                "--max-synteny-genomes", "200", "--min-recovery", "0.70", "--skip-tree"],
               cwd=str(DEPLOY))
        except subprocess.CalledProcessError as e:
            # Read the engine's fatal reason (most often HMM self-recovery).
            reason = ""
            man = bench / "benchmark_manifest.json"
            if man.exists():
                try:
                    reason = (json.loads(man.read_text()).get("fatal_error", "") or "")
                except Exception:
                    reason = ""
            is_recovery = "recovery" in reason.lower()

            # i > 1: a LATER round failing must NOT discard the earlier successful
            # rounds. Iterative refinement naturally broadens the family, so the
            # expanded seed set can dip below the self-recovery gate — that just
            # means we've gone as far as we usefully can. Stop iterating, drop the
            # partial failed round, and proceed to downstream with the best run so far.
            if i > 1:
                stop_reason = (f"stopped at run{i}: the expanded seed set failed HMM "
                               f"self-validation ({reason or 'search failed'}); kept the "
                               f"run{i - 1} results")
                log(f"RUN {i}: {reason or 'search failed'}")
                log(f"RUN {i}: not fatal — keeping run{i - 1}'s validated hits and proceeding "
                    f"to downstream (iterative refinement reached its useful limit).")
                shutil.rmtree(run_dir, ignore_errors=True)  # discard the partial round
                break

            # i == 1: nothing to fall back to → fail cleanly (no traceback) with guidance.
            if is_recovery:
                log(f"RUN 1: HMM self-validation failed — {reason}")
                sys.exit(
                    "\nThe profile HMM failed self-validation: too few seed sequences are "
                    "recovered by the model built from them (this is the --min-recovery gate, "
                    "default 0.70). It almost always means the seeds are NOT one coherent "
                    "homologous family — e.g. a broad functional label (\"endolysin\") spanning "
                    "unrelated proteins. Fixes: use a tighter single-family seed set, drop "
                    "outlier seeds, or — only for a deliberately divergent superfamily — relax "
                    "the recovery requirement. (Run aborted before wasting time on the search.)")
            sys.exit(f"\nRUN 1 search failed (exit {e.returncode}); reason: "
                     f"{reason or 'see the log above'}. Output dir: {bench}")

        # 2. ORF-validated extraction (NT + AA + TSV + next-seed FASTA).
        #    --email (only when set) lets it retrieve protein-database hit sequences;
        #    offline, those hits are kept without the NCBI-fetched AA.
        extract_cmd = ["python3", str(EXTRACTOR), "--results-dir", str(bench / "results"),
            "--hmm", str(bench / "hmm" / "benchmark_profile.hmm"),
            "--run-label", str(i), "--out-dir", str(validated),
            "--cpu", args.cpu]
        if args.email:
            extract_cmd += ["--email", args.email]
        if args.prodigal_gate:
            extract_cmd.append("--prodigal-gate")
        sh(extract_cmd)

        # 2b. add the source-organism (phage name) column to the hit table
        if args.no_annotate:
            log("  (organism annotation skipped: --no-annotate)")
        else:
            try:
                sh(["python3", str(ANNOTATE), "--hits-tsv", str(validated / "hits.tsv"),
                    "--email", args.email])
            except Exception as e:
                log(f"  (organism annotation skipped: {e})")

        n = (sum(1 for ln in unique.read_text().splitlines() if ln.startswith(">"))
             if unique.exists() else 0)
        curr_leng = _hmm_leng(hmm_of(bench))
        iter_hits.append((i, n))
        log(f"RUN {i}: complete -> {n} unique validated seeds (HMM length {curr_leng}).")
        if n == 0:
            stop_reason = f"run{i} found no validated hits (nothing to seed the next round)"
            log(f"RUN {i}: no validated hits — stopping iterations.")
            break
        # Early-stop on convergence: when the hit set AND the model both stabilise,
        # further rounds add nothing (matches METHODOLOGY.md §6). Reuses the engine's
        # convergence_check (count change <5 % AND ΔLENG <3).
        if (convergence_check is not None and prev_n is not None
                and convergence_check(prev_n, n, prev_leng or 0, curr_leng)):
            stop_reason = (f"converged at run{i} (hits {prev_n}→{n}, "
                           f"HMM length {prev_leng}→{curr_leng})")
            log(f"RUN {i}: {stop_reason} — stopping early.")
            seed = unique
            break
        prev_n, prev_leng = n, curr_leng
        seed = unique

    # Canonical run for figures + calibration = the most complete run (== the
    # paper table's best_run). After convergence this is the refined final round,
    # so figures and tables describe the SAME hit set (was figures=run1, table=best).
    best_i = _best_run_index(out, args.iterations)
    rbest = out / f"run{best_i}" / "benchmark"
    down = out / "downstream"

    # Threshold calibration on the canonical run's HMM: sensitivity (seed
    # self-test) + false-positive rate (shuffled/unrelated negatives).
    control_summary: dict = {}
    if not args.no_controls:
        if args.download_controls and download_control_sequences is not None:
            try:
                log("Downloading UniProt unrelated-proteome negative controls (one-time)…")
                download_control_sequences()
            except Exception as e:
                log(f"  (control download skipped: {e})")
        control_summary = run_controls(
            rbest / "hmm" / "benchmark_profile.hmm", fasta, out,
            args.biology_mode, 45.0, 30.0, args.cpu, log)

    # Per-seed recovery QC (named, before vs after): which of the original input
    # seeds does the model actually recover — against the INITIAL model (run1) and
    # the FINAL refined model (best run)? Surfaces divergent-outlier seeds by name
    # (the controls only give an aggregate count). Cheap; skipped in smoke.
    seed_recovery_summary: dict = {}
    if not args.smoke:
        try:
            from seed_recovery import seed_recovery_report  # noqa: E402  (sibling)
            seed_recovery_summary = seed_recovery_report(
                fasta,
                out / "run1" / "benchmark" / "hmm" / "benchmark_profile.hmm",   # before
                rbest / "hmm" / "benchmark_profile.hmm",                        # after
                out / "seed_qc", args.cpu, log)
        except Exception as e:
            log(f"  (seed-recovery QC skipped: {e})")

    # Optional: read-through scan for stop-interrupted / overprinted homologs that
    # the stop-to-stop search cannot see (opt-in; scans the nucleotide DBs again).
    interrupted_summary: dict = {}
    if getattr(args, "find_interrupted", False) and not args.smoke:
        interrupted_summary = run_find_interrupted(
            out, rbest / "hmm" / "benchmark_profile.hmm", args.db_cache,
            args.databases, args.cpu, log, control_summary)

    # Consolidated methodology / reproducibility record at the run root.
    write_methods_log(out, args, fasta, label, args.databases, iter_hits,
                      started_at, log, stop_reason, control_summary,
                      seed_recovery_summary, interrupted_summary)
    write_csv_exports(out, log)

    if args.smoke:
        # Smoke test = prove search + extraction work on new input; skip the
        # heavy downstream (tree/synteny need many hits + genomic coordinates).
        log("Smoke mode: skipping clustering/clinker/tree/GenBank downstream.")
        for i in range(1, args.iterations + 1):
            tsv = out / f"run{i}" / "benchmark" / "validated" / "hits.tsv"
            if tsv.exists():
                write_gff3(tsv, out / f"run{i}" / "benchmark" / "validated" / "hits.gff3")
        log(f"=== SMOKE TEST DONE. Hits/sequences in {rbest / 'validated'} ===")
        return

    # --- downstream analyses on the most complete run (primary discovery) ----
    log(f"Downstream: clustering + clinker + tree (on run{best_i}, the most complete run)")
    cluster_cmd = ["python3", str(CLUSTER),
        "--validated-dir", str(rbest / "validated"),
        "--cache-dir", str(rbest / "results" / "synteny_context_cache"),
        "--out-dir", str(down / "clinker")]
    if not args.no_annotate:  # resolve protein-DB hits to neighbourhoods via coded_by
        cluster_cmd += ["--email", args.email]
    sh(cluster_cmd)
    # publication-quality static synteny panels (anchored, orthogroup-coloured),
    # built from clinker's GenBank neighbourhoods. Non-fatal if it fails.
    try:
        synteny_cmd = ["python3", str(SYNTENY),
            "--clinker-dir", str(down / "clinker"),
            "--out-dir", str(down / "synteny"),
            "--annotation-cache", str(args.db_cache), "--cpu", args.cpu,
            "--color-by", args.color_by]
        if args.synteny_gene_labels:
            synteny_cmd.append("--gene-labels")
        sh(synteny_cmd)
    except Exception as e:
        log(f"  (synteny figures skipped: {e})")
    sh(["python3", str(TREE),
        "--faa", str(rbest / "validated" / "hits_unique_aa.faa"),
        "--out-dir", str(down / "tree"), "--cpu", args.cpu,
        "--hits-tsv", str(rbest / "validated" / "hits.tsv"),
        "--seeds", str(fasta)])   # place the (marked) seeds within the homolog tree
    # Per-hit alignment to the family HMM (each hit mapped onto the model), beside
    # the all-vs-all MSA the tree is built from.
    run_perhit_hmm_alignment(rbest / "hmm" / "benchmark_profile.hmm",
                             rbest / "validated" / "hits_unique_aa.faa", down / "tree", log)

    # --- per-run GFF3 (genome-browser tracks) --------------------------------
    log("Writing GFF3 tracks per run")
    for i in range(1, args.iterations + 1):
        tsv = out / f"run{i}" / "benchmark" / "validated" / "hits.tsv"
        if tsv.exists():
            write_gff3(tsv, out / f"run{i}" / "benchmark" / "validated" / "hits.gff3")

    # --- real-sequence GenBank neighbourhoods (open in Artemis/Geneious) -----
    log(f"Building real-sequence GenBank neighbourhoods (run{best_i})")
    try:
        genbank_cmd = ["python3", str(GENBANK),
            "--hits-tsv", str(rbest / "validated" / "hits.tsv"),
            "--out-dir", str(down / "genbank_with_sequence")]
        if args.email:
            genbank_cmd += ["--email", args.email]
        sh(genbank_cmd)
    except Exception as e:
        log(f"  (GenBank build skipped: {e})")

    # --- assemble a labelled package -----------------------------------------
    log("Assembling labelled output package")
    assemble_package(out, args.iterations, log, best_i)
    write_csv_exports(out, log)  # re-run now that PACKAGE exists, to mirror CSVs into it
    write_report(out, log)       # one-page HTML summary (links tables, tree, clinker)
    # README.txt in the package root + every subfolder, describing each file's
    # purpose. Done last so it reflects everything mirrored into PACKAGE/.
    try:
        from package_layout import write_readmes
        write_readmes(out / "PACKAGE", log)
    except Exception as e:
        log(f"  (package READMEs skipped: {e})")

    log(f"=== DONE. Package: {out / 'PACKAGE'} ===")


def assemble_package(out: Path, iterations: int, log, best_i: int = 1) -> None:
    """Copy the important outputs into a labelled, self-contained PACKAGE/."""
    pkg = out / "PACKAGE"
    pkg.mkdir(exist_ok=True)

    def cp(src, dst):
        src, dst = Path(src), Path(dst)
        if not src.exists():
            return
        dst.parent.mkdir(parents=True, exist_ok=True)
        if src.is_dir():
            shutil.rmtree(dst, ignore_errors=True)
            shutil.copytree(src, dst, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        else:
            shutil.copy2(src, dst)

    from package_layout import DIRS, PER_RUN  # single source of truth for the layout

    # Provenance + the human-facing report copied in so the package is self-contained.
    for f in ("METHODS.md", "run_manifest.json"):
        cp(out / f, pkg / f)
    # Stop-interrupted/overprinted homolog table + protein FASTAs (only present
    # with --find-interrupted). The .faa keep '*' at each internal stop.
    cp(out / "interrupted_homologs.tsv", pkg / DIRS["tables"] / "interrupted_homologs.tsv")
    for f in ("interrupted_homologs_domain_aa.faa", "interrupted_homologs_full_orf_aa.faa",
              "interrupted_homologs_full_orf_nt.fna"):
        cp(out / f, pkg / DIRS["sequences"] / f)
    # Publish the most refined HMM (from the canonical/most-complete run), which
    # is what the figures + paper table describe — not necessarily run1's model.
    cp(out / f"run{best_i}" / "benchmark" / "hmm" / "benchmark_profile.hmm",
       pkg / DIRS["hmm"] / "profile.hmm")
    for i in range(1, iterations + 1):
        v = out / f"run{i}" / "benchmark" / "validated"
        for f in ["hits.tsv", "hits.gff3", "hits_aa.faa", "hits_nt.fna",
                  "hits_unique_aa.faa", "orfs_aa.faa", "orfs_nt.fna"]:
            cp(v / f, pkg / DIRS["sequences"] / PER_RUN / f"run{i}" / f)
        cp(out / f"run{i}" / "benchmark" / "results" / "all_database_summary.tsv",
           pkg / DIRS["dbsum"] / f"run{i}_summary.tsv")
    cp(out / "downstream" / "clinker", pkg / DIRS["synteny"] / "clinker")
    cp(out / "downstream" / "synteny", pkg / DIRS["synteny"] / "publication_figures")
    cp(out / "downstream" / "genbank_with_sequence", pkg / DIRS["synteny"] / "genbank_with_sequence")
    cp(out / "seed_qc", pkg / DIRS["seedqc"])
    cp(out / "downstream" / "tree", pkg / DIRS["phylo"])
    # The reproducibility copy of scripts/ must EXCLUDE scratch/run helpers
    # (gitignored as scripts/_*): they can embed local paths or the user's NCBI
    # e-mail and must never be shipped. Keep dunder files (e.g. __init__.py).
    _sdst = pkg / DIRS["scripts"]
    shutil.rmtree(_sdst, ignore_errors=True)
    _sdst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(Path(__file__).resolve().parent, _sdst, ignore=_ignore_scratch)
    log(f"  package assembled at {pkg}")


if __name__ == "__main__":
    main()
