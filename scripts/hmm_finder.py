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
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

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


def write_methods_log(out: Path, args, fasta: Path, label: str, selected_dbs: str,
                      iter_hits: list, started_at: str, log, stop_reason: str = "",
                      control_summary: dict | None = None) -> None:
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

        manifest = {
            "tool": "hmm-homologue-finder",
            "code_git_commit": _git_commit(HERE),
            "started_at": started_at,
            "finished_at": finished_at,
            "command_line": " ".join(sys.argv),
            "conda_env": os.environ.get("CONDA_DEFAULT_ENV", ""),
            "conda_prefix": os.environ.get("CONDA_PREFIX", ""),
            "python": sys.version.split()[0],
            "parameters": {
                "label": label, "iterations": args.iterations, "cpu": args.cpu,
                "databases": selected_dbs, "prodigal_gate": bool(args.prodigal_gate),
                "min_recovery": "0.70", "max_synteny_genomes": "200",
                "email": args.email, "db_cache": str(args.db_cache), "out_dir": str(out),
                "input_type": args.input_type, "trans_table": args.trans_table,
                "no_annotate": bool(args.no_annotate),
            },
            "annotation_database": _annotation_provenance(args.db_cache),
            "input": {"fasta": str(fasta), "sha256": _sha256(fasta)},
            "per_iteration_unique_seeds": [{"run": i, "unique_validated_seeds": n} for i, n in iter_hits],
            "iteration_stop_reason": stop_reason,
            "threshold_calibration": control_summary or {},
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
        log(f"  Controls (strict bit≥{strict}): sensitivity {summary.get('sensitivity')}, "
            f"specificity {summary.get('specificity')}, FPR {summary.get('false_positive_rate')} "
            f"({summary.get('true_positives',0)}/{summary.get('total_positives',0)} seeds recovered; "
            f"{summary.get('false_positives',0)}/{summary.get('total_negatives',0)} negatives passed)")
        return summary
    except Exception as e:
        log(f"  (controls skipped: {e})")
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
            # Fail cleanly (no Python traceback) with an actionable diagnosis — an
            # unattended run shouldn't crash with a stack trace. The most common
            # run-1 failure is HMM self-recovery: the seeds don't form one coherent
            # homologous family, so the model can't even recover its own seeds.
            reason = ""
            man = bench / "benchmark_manifest.json"
            if man.exists():
                try:
                    reason = (json.loads(man.read_text()).get("fatal_error", "") or "")
                except Exception:
                    reason = ""
            if "recovery" in reason.lower():
                log(f"RUN {i}: HMM self-validation failed — {reason}")
                sys.exit(
                    "\nThe profile HMM failed self-validation: too few seed sequences are "
                    "recovered by the model built from them (this is the --min-recovery gate, "
                    "default 0.70). It almost always means the seeds are NOT one coherent "
                    "homologous family — e.g. a broad functional label (\"endolysin\") spanning "
                    "unrelated proteins. Fixes: use a tighter single-family seed set, drop "
                    "outlier seeds, or — only for a deliberately divergent superfamily — relax "
                    "the recovery requirement. (Run aborted before wasting time on the search.)")
            sys.exit(f"\nRUN {i} search failed (exit {e.returncode}); reason: "
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

    # Consolidated methodology / reproducibility record at the run root.
    write_methods_log(out, args, fasta, label, args.databases, iter_hits,
                      started_at, log, stop_reason, control_summary)
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
        sh(["python3", str(SYNTENY),
            "--clinker-dir", str(down / "clinker"),
            "--out-dir", str(down / "synteny"),
            "--annotation-cache", str(args.db_cache), "--cpu", args.cpu,
            "--color-by", "both"])
    except Exception as e:
        log(f"  (synteny figures skipped: {e})")
    sh(["python3", str(TREE),
        "--faa", str(rbest / "validated" / "hits_unique_aa.faa"),
        "--out-dir", str(down / "tree"), "--cpu", args.cpu,
        "--hits-tsv", str(rbest / "validated" / "hits.tsv"),
        "--seeds", str(fasta)])   # place the (marked) seeds within the homolog tree

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
            shutil.copytree(src, dst)
        else:
            shutil.copy2(src, dst)

    # Publish the most refined HMM (from the canonical/most-complete run), which
    # is what the figures + paper table describe — not necessarily run1's model.
    cp(out / f"run{best_i}" / "benchmark" / "hmm" / "benchmark_profile.hmm",
       pkg / "01_hmm_profile" / "profile.hmm")
    for i in range(1, iterations + 1):
        v = out / f"run{i}" / "benchmark" / "validated"
        for f in ["hits.tsv", "hits.gff3", "hits_aa.faa", "hits_nt.fna",
                  "hits_unique_aa.faa", "orfs_aa.faa", "orfs_nt.fna"]:
            cp(v / f, pkg / "02_sequences_per_run" / f"run{i}" / f)
        cp(out / f"run{i}" / "benchmark" / "results" / "all_database_summary.tsv",
           pkg / "03_database_summaries" / f"run{i}_summary.tsv")
    cp(out / "downstream" / "clinker", pkg / "04_synteny_clinker")
    cp(out / "downstream" / "synteny", pkg / "04_synteny_clinker" / "publication_figures")
    cp(out / "downstream" / "genbank_with_sequence", pkg / "04_synteny_clinker" / "genbank_with_sequence")
    cp(out / "seed_qc", pkg / "00_seed_qc")
    cp(out / "downstream" / "tree", pkg / "05_phylogeny")
    cp(Path(__file__).resolve().parent, pkg / "06_scripts")
    log(f"  package assembled at {pkg}")


if __name__ == "__main__":
    main()
