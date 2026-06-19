# Troubleshooting

## `conda: command not found`
Conda isn't installed or not on PATH. Install Miniforge (see [INSTALL.md](INSTALL.md)),
then close and reopen the terminal. `setup.sh` also searches common install
locations (`~/miniforge3`, `~/miniconda3`, `~/anaconda3`, `~/mambaforge`).

## The pipeline says required software is missing
Run the installer / repair:
```bash
bash setup.sh
python3 scripts/check_tools.py --install
```
`check_tools.py` prints a ✓/✗ table so you can see exactly what's absent.

## Windows: "WSL2 is required"
The tools have no native Windows build. In an **Administrator** PowerShell:
```powershell
wsl --install -d Ubuntu
```
Reboot, finish Ubuntu setup, then double-click `run.bat` again, or run
`bash run.sh` inside Ubuntu (`cd /mnt/c/…/hmm-homologue-finder`).

## A database download failed (e.g. "Could not resolve host")
A transient network error. The run records that database as `fail*` (not a
zero-hit result) and continues with the others. Re-run with the same `--out-dir`
to retry just the missing pieces. Large databases (GVD-AVrC ~4.6 GB,
GPD ~1.4 GB) need a stable connection.

## Self-search recovery failed / refuses to start a run
The HMM built from your seeds didn't recover enough of them (default ≥70%).
Usually means the seed set is too divergent or contains non-homologous/garbage
sequences. Curate the seeds (full-length, genuinely related) and retry. A quick
`--smoke` run surfaces this fast.

## 0 hits
Legitimate for a truly novel protein against annotated databases, or a sign the
seeds are off. Confirm with `--smoke` and check the seed alignment quality. The
pipeline handles 0 hits gracefully (writes empty outputs, stops iterating).

## Out of disk space
A full run streams several GB and presses HMM libraries. Keep ~20 GB free.
Outputs and caches live under your `--out-dir`; delete old run folders to reclaim
space. The `--databases` flag lets you skip the largest catalogues.

## It's taking hours
Expected for a full 3-iteration run over all databases (GVD-AVrC alone is ~2 h to
stream + translate). Use `--smoke` to validate quickly, `--databases` to limit
scope, and `nohup bash run.sh &` to run unattended.

## macOS: "cannot be opened because it is from an unidentified developer"
Right-click `scripts/Run HMM Homologue Finder.command` → **Open** → **Open**, or
just run `bash run.sh` from Terminal.

## NCBI rate-limiting / slow organism names
Provide your email (`--email you@inst.edu`) and avoid running many instances in
parallel. Organism lookups are cached within a run.

## Still stuck?
Open an issue on the GitHub repository with: your OS, the command you ran, and
the last ~30 lines of the terminal output (or the `pipeline.log` in your
out-dir).
