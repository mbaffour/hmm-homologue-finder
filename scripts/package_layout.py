#!/usr/bin/env python3
"""package_layout.py — the single source of truth for the PACKAGE/ folder layout.

Both the assembler (hmm_finder.assemble_package) and the table exporter
(export_csv) import the folder names from here so the structure stays consistent,
and `write_readmes()` drops a plain-text README.txt into every package folder
describing each file and its purpose (and one at the package root).
"""
from __future__ import annotations

import fnmatch
from pathlib import Path

# Canonical PACKAGE/ subfolders, in reading order. Names live ONLY here.
DIRS = {
    "tables":    "01_summary_tables",
    "sequences": "02_sequences",
    "hmm":       "03_hmm_profile",
    "phylo":     "04_alignment_phylogeny",
    "synteny":   "05_synteny",
    "dbsum":     "06_database_summaries",
    "seedqc":    "07_seed_qc",
    "scripts":   "08_scripts",
}
PER_RUN = "per_run"  # iteration subfolder under the sequences folder

# Per-folder purpose + per-file descriptions. Files are matched by exact name,
# then by glob pattern; anything unmatched gets a generic by-extension label.
_REGISTRY = {
    DIRS["tables"]: {
        "purpose": "The headline result tables and the per-database hit chart. Start here.",
        "files": {
            "paper_main_table.csv": "MAIN RESULT — one row per unique homolog (best per sequence): organism, accession, database, copies, #organisms, domain length/coverage, best E-value/bit score, confidence tier.",
            "hits_deduplicated.csv": "One row per unique homolog sequence, collapsed across databases/iterations, with 'found in N databases / N organisms' provenance.",
            "hit_summary.csv": "Per-iteration totals: hits, passed-filter, six-frame vs protein-DB, unique sequences, unique organisms, databases.",
            "database_hit_summary.csv": "Every database searched (including 0-hit ones), with unique sequence/organism counts.",
            "database_hits.png": "Bar chart of hits per database (grey = searched, 0 hits).",
            "database_hits.svg": "Bar chart, editable vector (Inkscape/Illustrator).",
            "database_hits.pdf": "Bar chart, print-ready vector.",
            "genome_metadata.csv": "Supplementary S1 — one row per source genome/contig (organism, host, databases, #hits).",
            "homolog_stats.csv": "Supplementary S3 — per-hit homology statistics (E-value, bit, domain length/coverage, confidence tier).",
            "all_runs_hits.csv": "The complete, un-collapsed hit table across all iterations (every hit, every column).",
            "database_summary.csv": "Raw engine per-database summary across all iterations.",
            "interrupted_homologs.tsv": "(--find-interrupted) Homologs interrupted by a premature stop codon, found by read-through translation — candidate overprinted/pseudogenized genes the stop-to-stop search misses. Columns: contig, strand, frame, genome coordinates (domain_nt_start/end), internal_stops, stop positions in the genome (stop_nt_positions) and in aa, bit score, residues before/after the stop, the domain DNA (domain_nt) and AA (with '*'), and the full read-through ORF (full_orf_aa).",
        },
    },
    DIRS["sequences"]: {
        "purpose": "All discovered sequences — combined multi-FASTAs plus per-iteration files.",
        "files": {
            "all_hits_aa.faa": "Every validated hit, amino-acid (all iterations).",
            "all_hits_nt.fna": "Every validated hit, nucleotide (all iterations).",
            "unique_homologs_aa.faa": "One sequence per unique homolog, rich headers (organism, #organisms, #databases).",
            PER_RUN: "Per-iteration sequences and the full per-hit evidence table (see per_run/README.txt).",
        },
    },
    PER_RUN: {
        "purpose": "Outputs for a single search iteration.",
        "files": {
            "hits.tsv": "Full per-hit evidence table (column schema in docs/OUTPUTS.md).",
            "hits.csv": "Same table, comma-separated for Excel.",
            "hits.gff3": "Genome-browser track of every hit (IGV / JBrowse / Artemis).",
            "hits_aa.faa": "Homolog domain, amino-acid.",
            "hits_nt.fna": "Homolog domain, nucleotide.",
            "orfs_aa.faa": "Full surrounding ORF, amino-acid.",
            "orfs_nt.fna": "Full surrounding ORF, nucleotide.",
            "hits_unique_aa.faa": "Deduplicated domains used to seed the next iteration.",
        },
    },
    DIRS["hmm"]: {
        "purpose": "The calibrated profile HMM — the model the whole run is built on.",
        "files": {
            "profile.hmm": "Profile HMM (most-complete run). Search with HMMER; submit to Pfam / NCBI CDD / VOGDB.",
        },
    },
    DIRS["phylo"]: {
        "purpose": "Multiple-sequence alignment of the homologs and the maximum-likelihood tree.",
        "files": {
            "hits.treefile": "ML tree (Newick) of the homologs — seeds included and marked SEED_*. Open in iTOL/FigTree.",
            "hits.contree": "Bootstrap consensus tree (Newick).",
            "hits_tree.png": "Rendered tree figure (300 dpi).",
            "hits_tree.svg": "Rendered tree, editable vector.",
            "hits_tree.pdf": "Rendered tree, print-ready vector.",
            "hits_tree_homologs_only.*": "Tree with the seed tips pruned (homologs only).",
            "hits.aln.faa": "Multiple sequence alignment of the homologs (MAFFT).",
            "hits.aln.trim.faa": "Trimmed alignment (trimAl) — the tree input.",
            "hits.aln.stats.json": "Alignment quality: length, gap %, conserved columns, mean pairwise identity.",
            "alignment_figure.*": "ClustalX-coloured alignment figure (PNG/SVG/PDF).",
            "hits.iqtree": "IQ-TREE report (model selection, parameters, support).",
            "hits.log": "IQ-TREE run log.",
            "tree_input.faa": "Exact sequences fed to the aligner (organism-labelled).",
        },
    },
    DIRS["synteny"]: {
        "purpose": "Gene-neighbourhood (synteny) comparisons and real-sequence GenBank files.",
        "files": {
            "clinker": "Interactive clinker comparison — open the .html in any browser ('Save SVG' for figures).",
            "publication_figures": "Publication synteny panels (PNG/SVG/PDF) + neighbour_gene_annotations.csv.",
            "genbank_with_sequence": "Real-sequence GenBank neighbourhoods, named by phage (Artemis/Geneious/UGENE).",
        },
    },
    DIRS["dbsum"]: {
        "purpose": "The engine's per-database summary for each iteration (status, hit counts, provenance).",
        "files": {
            "run*_summary.tsv": "Per-database results for one iteration: database, status, hit/strict counts, runtime, provenance.",
        },
    },
    DIRS["seedqc"]: {
        "purpose": "Quality control of the INPUT seeds (before committing to the search).",
        "files": {
            "seed_recovery.csv": "Per-seed QC — each input seed's best score vs the INITIAL and FINAL model + recovered flag (status: recovered / lost_after_refinement / gained_after_refinement / never_recovered).",
            "tree_input.faa": "The seed sequences as aligned (organism-labelled).",
            "hits.aln.faa": "Alignment of just the input seeds (sanity-check the set).",
            "hits.treefile": "Seed-only QC tree (Newick).",
            "hits_tree.*": "Seed-only QC tree figure.",
        },
    },
    DIRS["scripts"]: {
        "purpose": "A copy of the exact scripts that produced this run (for reproducibility).",
        "files": {},
    },
}

_EXT_DESC = {
    ".csv": "table (comma-separated)", ".tsv": "table (tab-separated)",
    ".faa": "protein FASTA", ".fna": "nucleotide FASTA", ".fasta": "FASTA",
    ".png": "figure (raster, 300 dpi)", ".svg": "figure (editable vector)",
    ".pdf": "figure (print vector)", ".hmm": "HMM profile (HMMER)",
    ".gbk": "GenBank record", ".gff3": "genome-browser track",
    ".json": "machine-readable data", ".treefile": "phylogenetic tree (Newick)",
    ".contree": "bootstrap consensus tree (Newick)", ".nwk": "phylogenetic tree (Newick)",
    ".md": "documentation", ".html": "interactive web page", ".log": "log file",
    ".txt": "text", ".py": "Python script",
}


def _describe(name: str, is_dir: bool, files_map: dict) -> str:
    if name in files_map:
        return files_map[name]
    for pat, desc in files_map.items():
        if ("*" in pat or "?" in pat) and fnmatch.fnmatch(name, pat):
            return desc
    if is_dir:
        return "folder"
    return _EXT_DESC.get(Path(name).suffix.lower(), "output file")


def _folder_readme(folder: Path, key: str) -> str:
    info = _REGISTRY.get(key, {})
    purpose = info.get("purpose", "Run outputs.")
    files_map = info.get("files", {})
    lines = [f"{folder.name}/  —  {purpose}", "=" * 70, "",
             "Contents (file : what it is):", ""]
    entries = sorted([p for p in folder.iterdir() if p.name != "README.txt"],
                     key=lambda p: (not p.is_dir(), p.name.lower()))
    if not entries:
        lines.append("  (empty)")
    for p in entries:
        tag = p.name + ("/" if p.is_dir() else "")
        lines.append(f"  {tag:<34} {_describe(p.name, p.is_dir(), files_map)}")
    lines += ["", "Generated by hmm-homologue-finder. See ../README.txt for the whole package."]
    return "\n".join(lines) + "\n"


def _top_readme(pkg: Path) -> str:
    order = [("README.txt", "this file"),
             ("METHODS.md", "how this run was produced (methods + citations)"),
             ("run_manifest.json", "machine-readable provenance (parameters, tool versions, calibration, seed recovery)")]
    lines = ["PACKAGE/  —  self-contained, shareable results of an hmm-homologue-finder run",
             "=" * 74, "",
             "Open  ../report.html  for a one-page visual summary of this run.",
             "", "At this level:", ""]
    for n, d in order:
        if (pkg / n).exists() or n == "README.txt":
            lines.append(f"  {n:<22} {d}")
    lines += ["", "Folders (in reading order):", ""]
    for key in ("tables", "sequences", "hmm", "phylo", "synteny", "dbsum", "seedqc", "scripts"):
        d = pkg / DIRS[key]
        if d.exists():
            lines.append(f"  {DIRS[key] + '/':<26} {_REGISTRY[DIRS[key]]['purpose']}")
    lines += ["", "Every folder has its own README.txt listing each file's purpose."]
    return "\n".join(lines) + "\n"


def write_readmes(pkg: Path, log=None) -> None:
    """Write README.txt into the package root and every known subfolder. Never raises."""
    pkg = Path(pkg)
    try:
        if not pkg.exists():
            return
        # subfolder READMEs (incl. the nested sequences/per_run/runN folders)
        for key in DIRS.values():
            folder = pkg / key
            if folder.is_dir():
                (folder / "README.txt").write_text(_folder_readme(folder, key))
        seqdir = pkg / DIRS["sequences"] / PER_RUN
        if seqdir.is_dir():
            (seqdir / "README.txt").write_text(_folder_readme(seqdir, PER_RUN))
            for run in sorted(p for p in seqdir.iterdir() if p.is_dir()):
                (run / "README.txt").write_text(_folder_readme(run, PER_RUN))
        (pkg / "README.txt").write_text(_top_readme(pkg))
        if log:
            log(f"  wrote README.txt in PACKAGE/ and {sum(1 for _ in pkg.rglob('README.txt'))-1} subfolder(s)")
    except Exception as e:
        if log:
            log(f"  (package READMEs skipped: {e})")
