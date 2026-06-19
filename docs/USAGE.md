# Usage — every case

The single entry point is `scripts/hmm_finder.py`. The launchers (`run.sh`,
`run.bat`, the Mac `.command`) just activate the environment and call it.

```
python3 scripts/hmm_finder.py [--fasta FILE] [--name LABEL] [--databases "A,B,…"]
                                 [--iterations N] [--cpu N] [--email you@inst.edu]
                                 [--smoke] [--out-dir DIR] [--skip-tool-check]
```

Always `conda activate hmm-discovery` first (the launchers do this for you).

---

## Case 1 — Interactive (no flags)
Best for newcomers. It asks for the seed FASTA, then runs everything.
```bash
bash run.sh
#   seed FASTA >  <- drag your .fasta into the terminal and press Enter
```
On macOS you can instead double-click `scripts/Run HMM Homologue Finder.command`.

## Case 2 — Explicit, full run
```bash
python3 scripts/hmm_finder.py --fasta my_seeds.faa --name my_protein --cpu 8
```
Runs the default 3 iterations against all 10 databases, then clustering, synteny,
tree, motifs, GenBanks, and assembles the package.

## Case 3 — Fast self-test / sanity check (`--smoke`)
One iteration against a single small database (INPHARED proteins). Minutes, not
hours. Use it to confirm a new install or a new protein input is wired correctly.
```bash
python3 scripts/hmm_finder.py --fasta my_seeds.faa --smoke --name my_protein
```

## Case 4 — Choose databases (`--databases`)
Comma-separated names exactly as registered (see [DATABASES.md](DATABASES.md)).
```bash
# discovery only (skip the slow gut-virome catalogues):
python3 scripts/hmm_finder.py --fasta my_seeds.faa \
    --databases "INPHARED genomes,RefSeq viral genomes"

# add a host/background control + viral annotation:
python3 scripts/hmm_finder.py --fasta my_seeds.faa \
    --databases "INPHARED genomes,RefSeq viral genomes,RefSeq bacterial proteins,PHROGs (annotation)"
```
Available: INPHARED genomes, INPHARED proteins, SwissProt, RefSeq viral proteins,
RefSeq viral genomes, Gut Phage Database (GPD), GVD-AVrC, RefSeq bacterial
proteins, Pfam (sequences), Pfam (domain scan), VOGDB VFAM (annotation),
PHROGs (annotation).

## Case 5 — Control the number of iterations
```bash
python3 scripts/hmm_finder.py --fasta my_seeds.faa --iterations 1   # single search, no re-seeding
python3 scripts/hmm_finder.py --fasta my_seeds.faa --iterations 5   # iterate more before stopping
```
The pipeline also stops early on its own when a round finds no new validated
hits (converged).

## Case 6 — A different protein family
Exactly the same — point `--fasta` at that family's seed sequences and give it a
`--name`:
```bash
python3 scripts/hmm_finder.py --fasta depolymerase_seeds.faa --name depolymerase
```
Hits from genome databases (six-frame, ORF-validated) **and** from protein
databases (captured by accession, marked `source_type=annotated_protein`) are
both reported.

## Case 7 — Choose where outputs go
```bash
python3 scripts/hmm_finder.py --fasta my_seeds.faa --out-dir ~/results/run1
```
Default is `<fasta>_discovery/` next to your FASTA.

## Case 8 — Set your NCBI email (recommended for big runs)
NCBI is queried for organism names / protein-DB hit sequences. Provide a real
email to be a good API citizen:
```bash
python3 scripts/hmm_finder.py --fasta my_seeds.faa --email you@institution.edu
```

## Case 9 — Resume / re-run
Re-running with the same `--out-dir` skips iterations that already produced their
validated outputs, so an interrupted run continues where it left off.

## Case 10 — Long unattended runs
A full 3-iteration run streams several GB and takes hours. Run it detached and
let the machine stay awake:
```bash
nohup bash run.sh > run.log 2>&1 &
```
`run.sh` automatically keeps the machine awake (caffeinate on macOS,
systemd-inhibit on Linux) when available.

---

## Tips
- **Seed quality matters more than quantity** — a handful of curated, full-length
  family members makes a sharper HMM than many fragments.
- **Try `--smoke` first** on any new machine or new protein.
- **Watch the controls**: SwissProt/Pfam/VOGDB at ~0 hits indicates a specific,
  novel family; non-zero means your profile also matches well-known proteins.
- See [OUTPUTS.md](OUTPUTS.md) for what every result file means and which tool
  opens it.
