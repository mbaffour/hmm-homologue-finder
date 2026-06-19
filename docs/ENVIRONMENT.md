# Environment & dependencies

## Platform support (important)
The pipeline depends on bioinformatics tools (HMMER, MAFFT, Prodigal, CD-HIT,
IQ-TREE, MEME, clinker) distributed through **bioconda, which builds only for
macOS and Linux — there are no native Windows packages**.

| Platform | How to run |
|----------|-----------|
| **macOS** | `bash run.sh` (or double-click `scripts/Run HMM Homologue Finder.command`) |
| **Linux** | `bash run.sh` |
| **Windows** | **via WSL2** — double-click `run.bat` (detects WSL and runs inside it), or in Ubuntu: `bash run.sh`. Native Windows is not supported by the tools. |

One-time Windows setup (Administrator PowerShell): `wsl --install -d Ubuntu`,
reboot, open Ubuntu, then `cd /mnt/c/…/HMM_Homologue_Finder && bash run.sh`.

The search **engine travels bundled** inside this folder (`engine/`), so nothing
external needs to be cloned or path-configured — it works wherever the folder is.

## One-time setup
```bash
bash setup.sh        # macOS / Linux / WSL2
```
This:
1. Finds conda/mamba (or tells you to install Miniforge if absent).
2. Creates the **`hmm-discovery`** conda environment from the deployable repo's
   `environment.yml`.
3. Verifies every required tool and installs any that are missing.

If you don't have conda, install Miniforge first (one time):
```bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-$(uname -m).sh
bash Miniforge3-MacOSX-$(uname -m).sh
# close & reopen the terminal, then:  bash setup.sh
```

## Required software (versions used)
| Tool | Version | Used for |
|------|---------|----------|
| HMMER | 3.4 | hmmbuild / hmmsearch / hmmscan |
| MAFFT | v7.526 | multiple sequence alignment |
| trimAl | v1.5 | alignment trimming |
| Prodigal | V2.6.3 | gene prediction / ORF validation |
| seqkit | v2.13.0 | FASTA chunking/handling |
| CD-HIT | 4.8.1 | sequence clustering |
| IQ-TREE | 3.1.2 | maximum-likelihood phylogeny |
| MEME / FIMO | 5.5.9 | motif discovery / scanning |
| clinker | v0.0.32 | synteny figures |
| curl | (system) | database downloads |
| Python | ≥3.10 | + Biopython, pandas |

## Dependency on the search engine
The heavy six-frame search (`run_all_database_benchmark.py`) and the `pipeline/`
package live in the **HMM-Discovery deployable repository**, expected at:
```
~/Documents/HMM-Discovery-Deployable-20260602/
```
`hmm_finder.py` calls it automatically and resolves the conda env and that
repo by absolute path, so the tool folder can live anywhere.

## Verify at any time
```bash
python3 scripts/check_tools.py            # report
python3 scripts/check_tools.py --install  # report + auto-install missing
```
