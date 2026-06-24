# HMM Homologue Finder

A reproducible, one-command pipeline that finds distant homologues of a protein
family across public phage/viral sequence databases — **including homologues
encoded by genes that standard annotation misses** — and returns validated
sequences (DNA + protein), evidence tables, gene-neighbourhood (synteny)
figures, a phylogenetic tree, motifs, and a publication-ready output package.

Give it a seed FASTA; it does the rest — and installs the software it needs on
first run.

> A general tool: point it at any protein family's seed sequences. (It was
> originally developed for phage protein-family discovery.)

---

## Contents
- [Quick start](#quick-start)
- [Platforms](#platforms)
- [What it produces](#what-it-produces)
- [How it works](#how-it-works)
- [Documentation](#documentation)
- [Repository layout](#repository-layout)
- [Requirements](#requirements)
- [Citation & license](#citation--license)

---

## Quick start

```bash
# macOS / Linux / Windows-WSL2
git clone https://github.com/mbaffour/hmm-homologue-finder.git
cd hmm-homologue-finder
bash setup.sh                 # one-time: creates the conda env, installs all tools
bash run.sh                   # interactive: prompts for your seed FASTA, runs everything
```

Prefer flags? 
```bash
conda activate hmm-discovery
python3 scripts/hmm_finder.py --fasta examples/example_seeds.fasta --smoke --name demo
```

**First time / new machine?** Run a fast self-test to confirm the install:
```bash
python3 scripts/hmm_finder.py --fasta examples/example_seeds.fasta --smoke
```

Full instructions for every case: **[docs/INSTALL.md](docs/INSTALL.md)** and
**[docs/USAGE.md](docs/USAGE.md)**. Browse **[docs/guide.html](docs/guide.html)**
for an interactive overview.

**Two modes:**
- **Discovery** (the default) — search the whole database catalog for a family, as above.
- **Single-genome scan** — *"does **this** genome carry my gene, and where?"* Build (or supply) an HMM and scan **one** genome — a local FASTA **or an NCBI accession** (fetched for you):
  ```bash
  python3 scripts/scan_genome.py --seeds gene_seeds.faa --genome my_genome.fna --out scan_out
  python3 scripts/scan_genome.py --hmm gene.hmm --accession KX098390 --email you@inst.edu --find-interrupted
  ```
  Reports present / not-detected, the hit's coordinates + sequence, ORF validation, the **flanking genes + genome-map figures** (`scan_neighbourhood.csv` + `scan_genome_map_*` — both a controllable `--flanks N` window *and* a whole-contig map, your gene in bold **gold**) — ordered relative to your gene, named from the **genome's own annotation** (gene/gp number, product, locus_tag, protein_id) or Prodigal for plain FASTA, **including any overlapping/overprint partner** (e.g. gp75's antisense RNA polymerase) — and (with `--find-interrupted`) interrupted/overprinted copies **plus the overprinting verdict**. Maps render by default with **DNA Features Viewer** (`--map-tool dfv`: clean strand arrows, overlapping genes auto-stacked so an overprint partner never hides your gene, de-overlapped labels with leader lines, and a full genome-coordinate axis); other choices are `pub` (built-in matplotlib — arrow-direction strand, packed lanes, full-length coordinate ruler; always available and the automatic fallback), `pygenomeviz`, and `easyfig`, with any unavailable renderer falling back to `pub`. A locus GenBank is **always** written too, so you can open the map in Easyfig / Artemis / clinker / pyGenomeViz yourself. Map styling options (all optional): `--palette colorblind` (colour-blind-safe Paul Tol palette, shared with the synteny figures), `--functional-labels` (tag your gene + its overprint partner with their functional category), `--module-brackets` (bracket contiguous same-category runs as modules), a legend with per-category **counts**, and automatic **multi-line wrapping** for genomes with >40 genes. Exit code 0 = found, 1 = absent (handy in loops over many genomes). Interactive: `bash start.sh` → mode 3. See [USAGE Case 11](docs/USAGE.md).

---

## Platforms

The pipeline relies on bioinformatics tools (HMMER, MAFFT, Prodigal, CD-HIT,
IQ-TREE, MEME, clinker) distributed via **bioconda, which builds only for macOS
and Linux** — there are **no native Windows packages**.

| Platform | How to run |
|----------|-----------|
| **macOS** | `bash run.sh`, or double-click `scripts/Run HMM Homologue Finder.command` |
| **Linux** | `bash run.sh` |
| **Windows** | **via WSL2** — double-click `run.bat` (detects WSL2 and runs inside it), or in Ubuntu: `bash run.sh`. One-time: `wsl --install -d Ubuntu` (Admin PowerShell). |

The search **engine is bundled** in `engine/`, so the repository is
self-contained — clone it anywhere and it runs (nothing else to download).

---

## What it produces

For each run (under `<fasta>_discovery/PACKAGE/`):

| Output | Description |
|--------|-------------|
| `hits.tsv` | One row per hit — 37 columns: organism, genomic coordinates, ORF-validation metrics, HMM statistics, and **both nucleotide & amino-acid sequence**. |
| `hits_deduplicated.csv` | One row per **unique** homolog, collapsing the same protein found across databases/iterations, with a **"found in N databases"** provenance column. |
| `hits_aa.faa` / `hits_nt.fna` | Homologue protein / DNA sequences. |
| `hits.gff3` | Genome-browser track of every hit (IGV/JBrowse/Artemis). |
| multiple sequence alignment | High-quality MSA of the homologs (MAFFT **L-INS-i** where tractable) + trimmed copy + quality stats, and a coloured **alignment figure (PNG + SVG + PDF)**. |
| GenBank neighbourhoods | Real-sequence `.gbk` per locus, named by phage (Artemis/Geneious). |
| synteny figures | Publication panels per cluster (**PNG + SVG + PDF**, editable text) + interactive clinker. |
| phylogenetic tree | IQ-TREE ML tree (fixed seed) of the homologs **with the seeds placed in it, marked** (Newick + PNG + SVG + PDF). A pre-run **seed QC tree** is also produced. |
| threshold calibration | `controls/control_report.json` — sensitivity / specificity / false-positive rate of the bit-score threshold vs positive & negative controls. |
| profile `.hmm` | The model — submit to Pfam / NCBI CDD / VOGDB. |
| per-hit HMM alignment | every homolog aligned to the family model (`hits_hmmalign.sto`/`.a2m`) — match states vs insertions. |
| **interrupted / overprinted homologs** *(opt-in, `--find-interrupted`)* | homologs broken by a *premature stop* that the stop-to-stop search misses — with the full read-through ORF (protein + DNA, to the real stop codon) **and an overprinting verdict**: whether the stop is synonymous in an open overlapping antisense frame (`overprinting_support` = strong/partial/none). |
| convergence + report | Per-round hit counts, calibration, an embedded coloured alignment + synteny panel, and a self-contained HTML results summary. |

All vector figures (synteny, alignment, tree) are written as **SVG (editable text — Inkscape) and PDF (Illustrator)** alongside a 300-dpi PNG.

See **[docs/OUTPUTS.md](docs/OUTPUTS.md)** for the full reference.

---

## How it works

1. **Build the HMM** — align seeds (MAFFT, accuracy-first **L-INS-i**), trim (trimAl), `hmmbuild`, validate self-recovery, and calibrate the score threshold against positive/negative controls.
2. **Choose & search databases** — an interactive run **prompts you to pick** which databases to search (an unattended run uses the full default set, or pass `--databases` / `--all-databases`). `hmmsearch` (E ≤ 1e-5); genome databases are six-frame translated so unannotated genes are reachable.
3. **Extract & validate** — reconstruct each hit's ORF from genomic coordinates, delimit the domain by the HMM envelope, confirm it's a genuine ORF (no internal stops; in a real coding locus). Save NT + AA + a 37-column table. Protein-database hits are captured by accession too.
4. **Iterate to convergence** — deduplicate hits, re-seed, repeat; **stop early** when the hit set and model stabilise.
5. **Characterise** — cross-database hit deduplication, CD-HIT clustering, clinker + publication synteny, an ML tree (seeds included), high-quality alignment, MEME/FIMO motifs, GFF3 tracks, named GenBank files.
6. **Package** — assemble a labelled, self-contained output folder with full provenance (`run_manifest.json`, `METHODS.md`).

### Why it works (the science, in one breath)

- **Profile HMMs reach the "twilight zone."** Position-specific scoring (a conserved catalytic residue weighted heavily, a floppy loop loosely) detects homology down to ~15–25% identity, where pairwise BLAST fails.
- **Six-frame translation finds the unfindable.** Genome databases are translated in all six frames, so homologues encoded **antisense or out-of-frame** to predicted genes — the ones standard annotation misses — are searchable. (gp75 was found *only* this way: zero hits in curated protein/domain DBs.)
- **ORF validation proves they're real genes.** Each hit is reconstructed frame-correctly, required to have **no internal stop codons**, and checked to sit in a coding locus — with the evidence in `hits.tsv`. An HMM score alone is never enough.
- **Controls + ROC quantify that hits aren't artefacts.** Composition-matched shuffled seeds and unrelated proteomes give sensitivity / specificity / FPR and an **ROC AUC** (1.0 = perfect separation), so you can rule out the "it's just composition bias" objection with numbers.
- **Convergence proves completeness; dual tree support (SH-aLRT + UFBoot) and conserved synteny corroborate it.** And everything is recorded for exact reproduction.
- **Overprinting, demonstrated — not just suspected.** `--find-interrupted` reads *through* premature stops to recover interrupted/overprinted homologs, then tests whether each stop is **synonymous in an open overlapping antisense frame** — the sequence signature of overprinting (how gp75 hides antisense to a virion RNA polymerase). Reported as a per-candidate `strong`/`partial`/`none` verdict with the genome coordinates to follow up.

**Full, educational walkthrough — how *and why* each step works:** open
**[docs/guide.html](docs/guide.html)** → *"Deep dive — why it works."*

Scientific detail: **[docs/METHODOLOGY.md](docs/METHODOLOGY.md)**.

### What it does *not* do (scope & limitations)
Briefly, so you don't over-read the output: it finds homologs by **sequence/HMM**
(not structure, e.g. Foldseek), from **assembled databases** (not raw reads).
**RNA-seq / read-based evidence is planned future work, not wired in now.** The
overprinting test is a *necessary* sequence signature — it confirms a stop is
synonymous in an open antisense frame, not that the antisense gene is expressed.
A validated hit is a **candidate** homolog; biological function still needs
experimental validation. Full list: **[METHODOLOGY.md §10](docs/METHODOLOGY.md#10-limitations--scope-what-it-does-not-do)**.

---

## Documentation

Start with the **interactive guide** — a single self-contained HTML page (works
offline) with a live command builder, the full method, outputs reference, and
tool versions:

> **[docs/guide.html](docs/guide.html)** — open in any browser.
> A tabbed quickstart is also at **[docs/HOW_TO_RUN.html](docs/HOW_TO_RUN.html)**.

Reference docs (Markdown):

| Doc | What's in it |
|-----|--------------|
| [docs/INSTALL.md](docs/INSTALL.md) | Install per platform (incl. Apple-Silicon `CONDA_SUBDIR`). |
| [docs/USAGE.md](docs/USAGE.md) | Every run case + the complete flag reference. |
| [docs/METHODOLOGY.md](docs/METHODOLOGY.md) | The scientific method, step by step, with tool versions. |
| [docs/OUTPUTS.md](docs/OUTPUTS.md) | `PACKAGE/` layout + `hits.tsv` schema + which tool opens what. |
| [docs/DATABASES.md](docs/DATABASES.md) | The database catalog + how to add your own. |
| [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) | Provenance, pinned env, determinism. |
| [docs/ENVIRONMENT.md](docs/ENVIRONMENT.md) | The `hmm-discovery` conda environment. |
| [docs/TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md) | Common errors and fixes. |

---

## Repository layout
```
hmm-homologue-finder/
├── README.md               ← you are here
├── run.sh                  ← launcher: macOS / Linux / WSL2
├── run.bat                 ← launcher: Windows → WSL2
├── setup.sh                ← one-time environment installer
├── environment.yml         ← conda environment definition
├── requirements.txt        ← Python deps (reference)
├── scripts/                ← the pipeline + helper scripts + Mac double-click launcher
├── engine/                 ← bundled search engine (run_all_database_benchmark.py + packages)
├── examples/               ← a tiny example seed FASTA for testing
└── docs/                   ← INSTALL, USAGE, OUTPUTS, DATABASES, TROUBLESHOOTING, METHODOLOGY, guide.html
```

---

## Requirements
- **conda / mamba** (install [Miniforge](https://conda-forge.org/download/) if you don't have it).
- The `hmm-discovery` conda environment — created automatically by `setup.sh`.
  Pinned tool versions (`environment.lock.yml`; each run also records the live
  versions it used in `run_manifest.json` / `METHODS.md`):

  | Tool | Version | | Tool | Version |
  |------|---------|-|------|---------|
  | HMMER | 3.4 | | CD-HIT | 4.8.1 |
  | MAFFT | v7.526 | | IQ-TREE | 3.1.2 |
  | trimAl | v1.5.1 | | MEME / FIMO | 5.5.9 |
  | Prodigal | V2.6.3 | | clinker | v0.0.32 |
  | seqkit | v2.13.0 | | Python | 3.12.13 |
  | Biopython | 1.87 | | pandas | 3.0.3 |
  | NumPy | 2.4.6 | | matplotlib | 3.11.0 |
  | playwright | 1.60.0 *(opt.)* | | | |

  *playwright + a chromium browser are optional — only for the static clinker PNG
  export; everything else runs without them.*

- Internet access for database streaming and NCBI organism lookups — **optional**:
  with no email / `--no-annotate` the pipeline runs fully offline (six-frame
  discovery is unaffected; only NCBI lookups are skipped).

---

## Citation & license
- License: **MIT** (see `LICENSE`).
- If you use this in research, please cite it (see `CITATION.cff`) and the
  underlying tools (HMMER, MAFFT, trimAl, Prodigal, CD-HIT, IQ-TREE, MEME/FIMO,
  clinker) and databases (INPHARED, RefSeq, GPD, GVD-AVrC, Pfam, VOGDB, PHROGs).
