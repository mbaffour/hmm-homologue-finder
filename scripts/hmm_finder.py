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
_cand = ([Path(os.environ["CONDA_PREFIX"]) / "bin"] if os.environ.get("CONDA_PREFIX") else []) \
    + [Path.home() / _n / "envs" / "hmm-discovery" / "bin"
       for _n in ("miniforge3", "mambaforge", "miniconda3", "anaconda3")]
for _b in _cand:
    if _b.is_dir():
        os.environ["PATH"] = f"{_b}{os.pathsep}{os.environ.get('PATH', '')}"
        break

BENCHMARK = DEPLOY / "scripts" / "run_all_database_benchmark.py"
EXTRACTOR = HERE / "extract_validated_hits.py"
CLUSTER = HERE / "cluster_and_clinker_corrected.py"
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


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", type=Path, default=None, help="seed protein FASTA (if omitted, you'll be prompted)")
    ap.add_argument("--out-dir", type=Path, default=None, help="output root (default: <fasta>_discovery)")
    ap.add_argument("--iterations", type=int, default=3, help="number of search iterations (default 3)")
    ap.add_argument("--cpu", default="8")
    ap.add_argument("--email", default="researcher@example.com", help="NCBI Entrez email for genome fetch")
    ap.add_argument("--name", default=None,
                    help="label for the output folder (default: derived from the FASTA name)")
    ap.add_argument("--databases", default=DATABASES,
                    help="comma-separated databases to search (must exist in the deployable "
                         "repo's databases.json; default = the 10 phage/viral databases)")
    ap.add_argument("--smoke", action="store_true",
                    help="fast self-test: 1 iteration against a single small database")
    ap.add_argument("--skip-tool-check", action="store_true", help="skip the startup software check")
    args = ap.parse_args()

    # --smoke: minutes-long install/sanity check on a single fast database.
    if args.smoke:
        args.iterations = 1
        args.databases = "INPHARED proteins"

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

    # Interactive: ask for the seed FASTA when none was given on the command line.
    if args.fasta is None:
        print("\n=== HMM-based homolog discovery pipeline ===")
        print("Drag your seed protein FASTA here (or type its path) and press Enter:")
        raw = input("  seed FASTA > ").strip().strip("'\"").strip()
        if not raw:
            sys.exit("No FASTA provided.")
        args.fasta = Path(raw)

    fasta = args.fasta.expanduser().resolve()
    if not fasta.exists():
        sys.exit(f"Seed FASTA not found: {fasta}")
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

    # Shared cache + DB indexes so each iteration reuses downloads (no re-fetch).
    shared = out / "_shared"
    (shared / "cache").mkdir(parents=True, exist_ok=True)
    (shared / "db_setup").mkdir(parents=True, exist_ok=True)

    log(f"=== family pipeline: {args.iterations} iterations from {fasta.name} ===")

    seed = fasta
    for i in range(1, args.iterations + 1):
        run_dir = out / f"run{i}"
        bench = run_dir / "benchmark"
        validated = bench / "validated"
        unique = validated / "hits_unique_aa.faa"

        if unique.exists() and unique.stat().st_size > 0:
            log(f"RUN {i}: already complete; skipping.")
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
        sh(["python3", str(BENCHMARK), "--fasta", str(seed), "--out", str(bench),
            "--databases", args.databases, "--cpu", args.cpu, "--keep-cache",
            "--max-synteny-genomes", "200", "--min-recovery", "0.70"],
           cwd=str(DEPLOY))

        # 2. ORF-validated extraction (NT + AA + TSV + next-seed FASTA).
        #    --email lets it retrieve protein-database hit sequences.
        sh(["python3", str(EXTRACTOR), "--results-dir", str(bench / "results"),
            "--hmm", str(bench / "hmm" / "benchmark_profile.hmm"),
            "--run-label", str(i), "--out-dir", str(validated), "--email", args.email])

        # 2b. add the source-organism (phage name) column to the hit table
        try:
            sh(["python3", str(ANNOTATE), "--hits-tsv", str(validated / "hits.tsv"),
                "--email", args.email])
        except Exception as e:
            log(f"  (organism annotation skipped: {e})")

        n = (sum(1 for ln in unique.read_text().splitlines() if ln.startswith(">"))
             if unique.exists() else 0)
        log(f"RUN {i}: complete -> {n} unique validated seeds.")
        if n == 0:
            log(f"RUN {i}: no validated hits — stopping iterations "
                "(nothing to seed the next round).")
            break
        seed = unique

    r1 = out / "run1" / "benchmark"
    down = out / "downstream"

    if args.smoke:
        # Smoke test = prove search + extraction work on new input; skip the
        # heavy downstream (tree/synteny need many hits + genomic coordinates).
        log("Smoke mode: skipping clustering/clinker/tree/GenBank downstream.")
        for i in range(1, args.iterations + 1):
            tsv = out / f"run{i}" / "benchmark" / "validated" / "hits.tsv"
            if tsv.exists():
                write_gff3(tsv, out / f"run{i}" / "benchmark" / "validated" / "hits.gff3")
        log(f"=== SMOKE TEST DONE. Hits/sequences in {r1 / 'validated'} ===")
        return

    # --- downstream analyses on the FIRST run (primary discovery) ------------
    log("Downstream: clustering + clinker + tree (on run1)")
    sh(["python3", str(CLUSTER),
        "--validated-dir", str(r1 / "validated"),
        "--cache-dir", str(r1 / "results" / "synteny_context_cache"),
        "--out-dir", str(down / "clinker")])
    sh(["python3", str(TREE),
        "--faa", str(r1 / "validated" / "hits_unique_aa.faa"),
        "--out-dir", str(down / "tree"), "--cpu", args.cpu])

    # --- per-run GFF3 (genome-browser tracks) --------------------------------
    log("Writing GFF3 tracks per run")
    for i in range(1, args.iterations + 1):
        tsv = out / f"run{i}" / "benchmark" / "validated" / "hits.tsv"
        if tsv.exists():
            write_gff3(tsv, out / f"run{i}" / "benchmark" / "validated" / "hits.gff3")

    # --- real-sequence GenBank neighbourhoods (open in Artemis/Geneious) -----
    log("Building real-sequence GenBank neighbourhoods (run1)")
    try:
        sh(["python3", str(GENBANK),
            "--hits-tsv", str(r1 / "validated" / "hits.tsv"),
            "--out-dir", str(down / "genbank_with_sequence"),
            "--email", args.email])
    except Exception as e:
        log(f"  (GenBank build skipped: {e})")

    # --- assemble a labelled package -----------------------------------------
    log("Assembling labelled output package")
    assemble_package(out, args.iterations, log)

    log(f"=== DONE. Package: {out / 'PACKAGE'} ===")


def assemble_package(out: Path, iterations: int, log) -> None:
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

    cp(out / "run1" / "benchmark" / "hmm" / "benchmark_profile.hmm",
       pkg / "01_hmm_profile" / "profile.hmm")
    for i in range(1, iterations + 1):
        v = out / f"run{i}" / "benchmark" / "validated"
        for f in ["hits.tsv", "hits.gff3", "hits_aa.faa", "hits_nt.fna",
                  "hits_unique_aa.faa", "orfs_aa.faa", "orfs_nt.fna"]:
            cp(v / f, pkg / "02_sequences_per_run" / f"run{i}" / f)
        cp(out / f"run{i}" / "benchmark" / "results" / "all_database_summary.tsv",
           pkg / "03_database_summaries" / f"run{i}_summary.tsv")
    cp(out / "downstream" / "clinker", pkg / "04_synteny_clinker")
    cp(out / "downstream" / "genbank_with_sequence", pkg / "04_synteny_clinker" / "genbank_with_sequence")
    cp(out / "downstream" / "tree", pkg / "05_phylogeny")
    cp(Path(__file__).resolve().parent, pkg / "06_scripts")
    log(f"  package assembled at {pkg}")


if __name__ == "__main__":
    main()
