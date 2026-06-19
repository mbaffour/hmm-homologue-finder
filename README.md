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
| `hits.tsv` | One row per hit — 35 columns: organism, genomic coordinates, ORF-validation metrics, HMM statistics, and **both nucleotide & amino-acid sequence**. |
| `hits_aa.faa` / `hits_nt.fna` | Homologue protein / DNA sequences. |
| `hits.gff3` | Genome-browser track of every hit (IGV/JBrowse/Artemis). |
| GenBank neighbourhoods | Real-sequence `.gbk` per locus, named by phage (Artemis/Geneious). |
| clinker figures | Interactive synteny comparisons per cluster. |
| phylogenetic tree | IQ-TREE maximum-likelihood tree (Newick + PNG/SVG). |
| profile `.hmm` | The model — submit to Pfam / NCBI CDD / VOGDB. |
| convergence + report | Per-round hit counts and a results summary. |

See **[docs/OUTPUTS.md](docs/OUTPUTS.md)** for the full reference.

---

## How it works

1. **Build the HMM** — align seeds (MAFFT), trim (trimAl), `hmmbuild`, validate self-recovery.
2. **Search 10 databases** — `hmmsearch` (E ≤ 1e-5); genome databases are six-frame translated so unannotated genes are reachable.
3. **Extract & validate** — reconstruct each hit's ORF from genomic coordinates, delimit the domain by the HMM envelope, confirm it's a genuine ORF (no internal stops; in a real coding locus). Save NT + AA + a 35-column table. Protein-database hits are captured by accession too.
4. **Iterate to convergence** — deduplicate hits, re-seed, repeat (default 3 rounds).
5. **Characterise** — CD-HIT clustering, clinker synteny, IQ-TREE phylogeny, MEME/FIMO motifs, GFF3 tracks, named GenBank files.
6. **Package** — assemble a labelled, self-contained output folder.

Scientific detail: **[docs/METHODOLOGY.md](docs/METHODOLOGY.md)**.

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
- The `hmm-discovery` conda environment — created automatically by `setup.sh`
  (HMMER 3.4, MAFFT, trimAl, Prodigal, seqkit, CD-HIT, IQ-TREE, MEME/FIMO,
  clinker, Biopython, pandas).
- Internet access (databases are streamed; NCBI is queried for organism names).

---

## Citation & license
- License: **MIT** (see `LICENSE`).
- If you use this in research, please cite it (see `CITATION.cff`) and the
  underlying tools (HMMER, MAFFT, trimAl, Prodigal, CD-HIT, IQ-TREE, MEME/FIMO,
  clinker) and databases (INPHARED, RefSeq, GPD, GVD-AVrC, Pfam, VOGDB, PHROGs).
