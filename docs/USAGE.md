# Usage — every case

The single entry point is `scripts/hmm_finder.py`. The launchers (`run.sh`,
`run.bat`, the Mac `.command`) just activate the environment and call it.

> **Fully autonomous from Python (no prompts, no `conda activate`):** use
> `scripts/run_pipeline.py` (it PATH-injects the env's tools, detaches stdin, and auto-injects
> `--all-databases`/`--no-annotate`/`--skip-tool-check`; flags `--scan`/`--preset`/`--dry-run`/
> `--preload`). **Warm the database caches ahead of time** with
> `scripts/preload_databases.py` (or `run_pipeline.py --preload`) so later runs spend their time
> searching, not downloading + six-frame-translating. For the **complete flag + output reference**,
> see **[REFERENCE.md](REFERENCE.md)**.

```
python3 scripts/hmm_finder.py [--fasta FILE] [--name LABEL] [--databases "A,B,…"]
                                 [--iterations N] [--cpu N] [--email you@inst.edu]
                                 [--smoke] [--out-dir DIR] [--skip-tool-check]
```

Always `conda activate hmm-discovery` first (the launchers do this for you).

> Prefer a point-and-click reference? Open **[guide.html](guide.html)** — it has a
> live command builder that writes the exact command for you.

## Complete flag reference

| Flag | Default | Meaning |
|------|---------|---------|
| `--fasta FILE` | — | **Only required arg.** Seed protein (or nucleotide) FASTA. |
| `--name LABEL` | FASTA stem | Output-folder label. |
| `--out-dir DIR` | `<fasta>_discovery/` | Output root; re-running resumes finished rounds. |
| `--no-overwrite` | off | If the output folder already holds a run, write to a fresh numbered one (`<dir>_2`, `_3`, …) instead of overwriting it. (The interactive wizard turns this on by default.) |
| `--iterations N` | 3 | Max re-seeding rounds (stops early on convergence / 0 hits). |
| `--cpu N` | 8 | Threads; auto-clamped to available cores. |
| `--email ADDR` | none | NCBI email. Precedence `--email` > `$NCBI_EMAIL` > TTY prompt > offline. Never hardcoded. |
| `--no-annotate` | off | Fully offline; skip all NCBI lookups (tables still build). |
| `--input-type {auto,protein,nucleotide}` | auto | Seed type; nucleotide seeds are translated first. |
| `--trans-table N` | 11 | NCBI genetic code — used to translate a nucleotide seed **and** for the read-through / interrupted-gene scan (e.g. 4 = Mycoplasma). |
| `--databases "A,B,…"` | full set | Exact catalog names. On a TTY without this → interactive picker. |
| `--all-databases` | off | Full default set without the picker (even on a TTY). |
| `--pick-databases` | off | Force the interactive database picker. |
| `--list-databases` | — | Print the catalog and exit (no FASTA needed). |
| `--smoke` | off | Fast self-test: 1 round vs one small DB; skips heavy downstream. |
| `--no-controls` | off | Skip threshold-calibration controls. |
| `--biology-mode {generic,phage,bacterial}` | phage | Which control panel to calibrate against. |
| `--download-controls` | off | One-time fetch of UniProt unrelated-proteome negatives. |
| `--prodigal-gate` | off | Stricter: require six-frame hits to overlap a Prodigal gene. |
| `--find-interrupted` | off | Also read-through-scan the nucleotide DBs for homologs interrupted by a premature stop (overprinted/pseudogenized genes the stop-to-stop search misses). Writes `interrupted_homologs.tsv` (incl. genome coordinates, the full read-through ORF, the actual stop codon's coordinate `natural_stop_nt`, and an **overprinting/silent-stop verdict** `overprinting_support` = strong/partial/none) plus three FASTAs: `interrupted_homologs_domain_aa.faa`, `…_full_orf_aa.faa`, `…_full_orf_nt.fna`. See [OUTPUTS.md](OUTPUTS.md). |
| `--no-seed-tree` | off | Skip the pre-run seed-only QC tree + alignment. |
| `--detach` *(run.sh)* | off | Run in the background in its **own session** so closing the terminal window (a SIGHUP) can't kill a long / walk-away run. Re-launches itself with `setsid`; console output goes to `hmm_run_console_<timestamp>.log` next to your seed. Check progress any time with `final/CHECK_RUN.bat`. A detached run has no terminal to prompt at, so pass `--email` (else it runs offline) and `--databases` (else the defaults are used). |
| `--synteny-gene-labels` | off | Label neighbour genes with function in synteny figures. |
| `--color-by {function,conservation,both}` | both | How to colour synteny neighbourhood genes. |
| `--skip-tool-check` | off | Skip the startup software preflight. |
| `--db-cache DIR` | `~/.cache/hmm-homologue-finder` | Persistent shared DB cache (download once, ever). |

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

## Case 11 — Scan ONE genome for your gene (`scan_genome.py`)
A focused mode, separate from the database-wide discovery pipeline: build (or
supply) a profile HMM and check **a single genome** for the gene — useful for
"does *this* assembly carry my gene of interest, and where?". The scan is
read-through, so a clean copy is found with 0 internal stops, and (with
`--find-interrupted`) a stop-interrupted / overprinted copy is reported too, with
the overprinting verdict.
```bash
# build the HMM from seeds and scan a local genome:
python3 scripts/scan_genome.py --seeds gene_seeds.faa --genome my_genome.fna --out scan_out

# reuse an existing profile HMM:
python3 scripts/scan_genome.py --hmm gene.hmm --genome my_genome.fna --out scan_out

# fetch the genome from NCBI by accession instead of a local file (needs --email):
python3 scripts/scan_genome.py --hmm gene.hmm --accession KX098390 --email you@inst.edu

# also catch interrupted/overprinted copies (with the overprinting verdict):
python3 scripts/scan_genome.py --seeds gene_seeds.faa --genome my_genome.fna --find-interrupted
```
The genome comes from **either** `--genome FILE` (local) **or** `--accession ACC`
(fetched from NCBI nucleotide via Entrez — comma-separate several; assembly
`GCF_`/`GCA_` ids aren't fetched directly, use their nucleotide/contig accessions).
Or, interactively: `bash start.sh` → **mode 3 (Scan ONE genome)** — it accepts a
local path or an accession. Outputs (in `--out`): `scan_report.txt` (present /
not-detected verdict + best hits), `scan_hits.tsv` (per-hit coordinates, ORF
validation, score, overprinting), `scan_hits_aa.faa` / `scan_hits_nt.fna` (hit
sequences), `scan_neighbourhood.csv` (the **flanking genes**, see below), the fetched
`<accession>.fna` when fetched, and `profile.hmm` when built from seeds. Exit code is
**0 if the gene was detected, 1 if absent** — handy in shell loops over many genomes.
Flags: `--min-bit` (default 25), `--trans-table` (for a nucleotide seed), `--cpu`.

**Flanking genes (`scan_neighbourhood.csv`) + a genome map.** By default the scan also
describes the genes around each hit, ordered relative to your gene (`pos_index` 0 = your
gene, ± up/downstream; `relationship` = upstream / downstream / **overlapping**;
`rel_start/rel_end`, `distance_to_anchor_bp`, `strand_vs_gene`, `length_bp`). It reports
the `--flanks` genes each side (default 7, contiguous — no gaps) **and every overlapping
gene**, so an **overprint partner is shown, not hidden** — e.g. for gp75 the antisense
*virion RNA polymerase* appears as `overlapping (antisense)`. **Names come from the
genome's OWN annotation** — `gene` (often the gp number), `product`, `locus_tag`,
`protein_id` — for any annotated record (a fetched `--accession`, pulled as GenBank, or a
GenBank `--genome` `.gb/.gbk/.gbff`); for a plain unannotated FASTA the neighbours are
called **de novo with Prodigal** (the database workflow's gene-caller) + optional VOGDB
VFAM (`--db-cache`); the `annotation_source` column says which. **Two genome-map figures** are drawn per hit with
**DNA Features Viewer** (the default renderer): clean strand **arrows**, **overlapping genes
automatically stacked onto their own level** (so an overprint partner never hides your gene),
label boxes **de-overlapped automatically with leader lines**, and a real **genome-coordinate
axis** spanning the whole locus. Genes are labelled from the genome's own annotation (gp
numbers / products) and **coloured by functional category** (the same scheme as the synteny
figures — structural, transcription/regulation, lysis, …), with your gene of interest in
**bold gold** and the track labelled with the **phage/organism name** (from the record). A
**smart label policy** keeps any genome legible: your gene of interest **and any gene
overlapping it are always labelled**, while other genes are labelled only when the locus is
small enough to stay clean (so a ~279-gene phage genome is not a wall of text). Turn off all
gene-name labels with `--no-gene-labels`. `scan_genome_map_<hit>.png/.svg/.pdf` is a
**window** of your gene + `--flanks N` genes each side (controllable);
`scan_genome_map_<hit>_whole.png/.svg/.pdf` is **the whole contig**, your gene marked among
all of them. The track is labelled with the phage name **over the accession**. Output is
**PNG (300 dpi) + SVG + PDF**. Choose the renderer with `--map-tool {dfv, pub, pygenomeviz,
easyfig}` (default `dfv`, DNA Features Viewer): `pub` is the built-in matplotlib diagram
(strand = arrow direction, overlapping genes packed onto separate lanes, a full-length
coordinate ruler) — it is **always available and is the automatic fallback** when a chosen
renderer is not installed; `easyfig` needs Easyfig installed + `$EASYFIG_PY` set. Any
unavailable renderer falls back to `pub`. A **locus GenBank** (`scan_genome_map_<hit>.gb`) is
**always written regardless of renderer**, so you can open the map in **Easyfig**, Artemis,
clinker, or pyGenomeViz yourself. Disable the whole step with `--no-neighbours`. *(The database
run draws a per-hit genome map too, honours the `GENOME_MAP_TOOL` environment variable as an
override, and writes the same per-hit GenBanks.)*

**Map styling (all optional, `dfv` renderer):** genes are coloured by broad functional category
with a legend that shows **per-category counts** (e.g. `structural (4)`); the gene of interest is
bold gold. `--palette {default, colorblind}` switches to a **colour-blind-safe** scheme (Paul Tol
muted; the same scheme is used by the synteny figures so a run is consistent). `--functional-labels`
tags the gene of interest **and its overprint partner** with their category (e.g. `[transcription]`)
— colour + legend still carry function for everything else, so it never re-crowds. `--module-brackets`
draws a bracket over each contiguous run of same-category genes labelled with the **module** name
(`structural module`, `replication module`, …; hypothetical genes are not bracketed). A genome with
**>40 genes wraps onto multiple lines** automatically so every gene stays legible (the whole-contig
map of a large phage no longer collapses into one dense row). *(The database run honours the same
choices via the `GENOME_MAP_PALETTE` / `GENOME_MAP_FUNCTIONAL=1` / `GENOME_MAP_BRACKETS=1` env vars.)*

*Worked example (real overprinting):* `--hmm gp75.hmm --accession KX098390
--find-interrupted` finds gp75 as **interrupted** with `overprinting_support=strong`
— its premature stop is synonymous in the open antisense RNA-polymerase frame.

---

## Tips
- **Seed quality matters more than quantity** — a handful of curated, full-length
  family members makes a sharper HMM than many fragments.
- **Try `--smoke` first** on any new machine or new protein.
- **Watch the controls**: SwissProt/Pfam/VOGDB at ~0 hits indicates a specific,
  novel family; non-zero means your profile also matches well-known proteins.
- See [OUTPUTS.md](OUTPUTS.md) for what every result file means and which tool
  opens it.
