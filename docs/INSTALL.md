# Installation

The only prerequisite you install by hand is **conda** (via Miniforge). Then
`setup.sh` builds everything else automatically. Pick your platform below.

---

## 1. Get the code
```bash
git clone https://github.com/mbaffour/hmm-homologue-finder.git
cd hmm-homologue-finder
```
(or download the ZIP from GitHub and unzip it).

---

## 2a. macOS

```bash
# Install Miniforge once if you don't have conda:
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-$(uname -m).sh"
bash Miniforge3-MacOSX-$(uname -m).sh         # accept defaults; say "yes" to init
# close & reopen Terminal, then:
bash setup.sh
```
Run it: `bash run.sh` — or double-click **`scripts/Run HMM Homologue Finder.command`** in Finder.

> On first launch macOS may warn about an unidentified developer. Right-click
> the `.command` → **Open** → **Open**, or run `bash run.sh` from Terminal.

---

## 2b. Linux

```bash
# Install Miniforge once if needed:
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh"
bash Miniforge3-Linux-$(uname -m).sh
# restart your shell, then:
bash setup.sh
bash run.sh
```

---

## 2c. Windows — via WSL2 (required)

The bioinformatics tools have **no native Windows builds**, so the pipeline runs
inside **WSL2** (Windows Subsystem for Linux). This is the standard approach for
tools like this.

**One-time WSL setup** (in an **Administrator** PowerShell):
```powershell
wsl --install -d Ubuntu
```
Reboot when prompted, and let Ubuntu finish first-time setup (it asks for a
username/password).

**Then either:**
- Double-click **`run.bat`** in the repo — it detects WSL and launches the
  pipeline inside it, **or**
- Open the **Ubuntu** app and run:
  ```bash
  # your C: drive is at /mnt/c inside WSL
  cd /mnt/c/Users/<you>/Downloads/hmm-homologue-finder
  # install Miniforge inside WSL once:
  curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
  bash Miniforge3-Linux-x86_64.sh
  exec bash
  bash setup.sh
  bash run.sh
  ```

---

## 3. Verify the install (any platform)
```bash
conda activate hmm-discovery
python3 scripts/check_tools.py            # prints a ✓/✗ table of every tool
python3 scripts/hmm_finder.py --fasta examples/example_seeds.fasta --smoke
```
The smoke test runs one round against a single small database and finishes in a
few minutes — proving the whole chain works on your machine.

---

## What `setup.sh` does
1. Finds conda/mamba (or tells you to install Miniforge).
2. Creates the `hmm-discovery` environment from `environment.yml`.
3. Verifies every required tool and installs any that are missing.

Re-running `setup.sh` is safe (idempotent). If a tool ever goes missing,
`python3 scripts/check_tools.py --install` repairs the environment.

---

## Troubleshooting install
See **[TROUBLESHOOTING.md](TROUBLESHOOTING.md)** — covers "conda not found",
WSL path issues, slow/failed database downloads, and disk-space tips.
