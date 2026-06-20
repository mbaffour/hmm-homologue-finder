# Reproducibility & citation

This pipeline is designed to meet journal reproducibility expectations: pinned
software and versioned databases with access dates and checksums. Everything runs
in the terminal via conda — **no Docker or other containers required.**

## 1. Reproducing the software environment (terminal only)

- **Pinned environment (most reproducible).** `environment.lock.yml` records the
  exact version of every package:
  ```bash
  conda env create -n hmm-discovery -f environment.lock.yml
  conda activate hmm-discovery
  ```

- **Standard setup.** On a fresh machine, `bash setup.sh` builds the env from
  `environment.yml`. For an exact-version rebuild, use the lock file above.

Then run entirely in the terminal:
```bash
bash start.sh                                  # guided
# or
bash run.sh --fasta seeds.faa --name MyFamily --email you@inst.edu
```

## 2. Per-run provenance (written automatically)

Every run writes, at the output root:

- `run_manifest.json` — command line, code git commit, input FASTA SHA-256, all
  parameters, tool versions (probed at run time), and per-database provenance
  (source URLs, SHA-256, first/last access timestamps) from each iteration.
- `METHODS.md` — a human-readable methods summary.
- `run*/benchmark/reports/reproducibility.json` — the engine's per-iteration record.

These pin *which* database files were used, when they were downloaded, and their
checksums — the information a reviewer needs to confirm reproducibility.

## 3. Databases (versioned, citable)

| Database | Version / release | Role |
|----------|-------------------|------|
| INPHARED | 1 Jan 2024 genomes/proteins | six-frame discovery search |
| RefSeq viral | release recorded per run in `run_manifest.json` | protein/genome search |
| VOGDB VFAM | release 230 (vog230) | functional annotation of neighbour genes (synteny figures); **optional** |

Database selection, exact source URLs, sizes, access dates and SHA-256 are
recorded per run. Functional annotation is optional — without VOGDB the synteny
figures still build (genes shown as "hypothetical").

## 4. Tools to cite

HMMER 3.4 (Eddy); MAFFT 7.526 (Katoh & Standley); trimAl 1.5 (Capella-Gutiérrez
et al.); Prodigal 2.6.3 (Hyatt et al.); seqkit 2.13 (Shen et al.); CD-HIT 4.8.1
(Fu et al.); IQ-TREE 3 (Minh et al.) with ModelFinder (Kalyaanamoorthy et al.)
and UFBoot (Hoang et al.); MEME Suite 5.5.9 (Bailey et al.); clinker (Gilchrist &
Chooi); VOGDB (vogdb.org); INPHARED (Cook et al.). Exact versions for a given run
are in that run's `run_manifest.json`.
