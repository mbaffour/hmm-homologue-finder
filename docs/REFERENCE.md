# HMM Homologue Finder — Reference

This document is the definitive, code-grounded reference for the **HMM Homologue Finder**: a Python
pipeline that discovers distant homologs of a protein family across phage/viral databases by iterative
profile-HMM search over **six-frame** translations, with an optional **read-through** mode for
stop-interrupted and overprinted genes (flagship case: **gp75**, an overprinted antisense ORF inside a
virion DNA-directed RNA polymerase). It catalogs every executable entry point, the search engine,
the interrupted/overprinting machinery, the annotation layer, the downstream figures, every output file
and column, and the complete CLI flag set — citing the real function names, flags, columns, and file
names in
`hmm-homologue-finder/scripts/*.py` and `hmm-homologue-finder/engine/`. It complements the narrative docs
(`docs/METHODOLOGY.md`, `docs/USAGE.md`, `docs/OUTPUTS.md`, `docs/README.md`); this file is the complete
catalog.

Two facts recur throughout and govern how results are read, so they are stated once here:

- **Which length to cite.** `domain_aa_len` — the **HMM-matched conserved region** (the alignment
  envelope, ≈137 aa for gp75) — is the homolog length to report. `orf_aa_len` is the *surrounding*
  six-frame ORF read from a Met start within `MET_MARGIN = 60 aa` of the domain; for an overprinted domain
  on a long, stop-free antisense frame it can be large and `has_start_M = False`, which is **correct**, not
  a defect.
- **Why `--prodigal-gate` is off by default.** Prodigal overlap is recorded as *informational* evidence,
  never exclusionary, because requiring overlap with a conventionally predicted gene would discard exactly
  the antisense/alternate-frame homologs the tool exists to find.

---

## Table of contents

1. [Entry points & how to run](#entry-points--how-to-run)
2. [The search core: HMM iteration, six-frame translation, and controls](#the-search-core-hmm-iteration-six-frame-translation-and-controls)
3. [Interrupted / overprinted genes](#interrupted--overprinted-genes)
4. [Annotation & organism naming](#annotation--organism-naming)
5. [Downstream figures & trees](#downstream-figures--trees)
6. [Outputs catalog — every file and column](#outputs-catalog--every-file-and-column)
7. [Complete CLI flag reference](#complete-cli-flag-reference)

---

## Entry points & how to run

The tool has four code entry points and two shell wizards. They form a layered stack: the shell wizards (`start.sh` → `run.sh`) handle conda activation and interactive option-gathering; the Python launcher (`scripts/run_pipeline.py`) is the zero-prompt path for scripts/cron; and the two engines underneath — `scripts/hmm_finder.py` (family discovery) and `scripts/scan_genome.py` (single-genome scan) — do the actual work. Every entry point can also be invoked directly with the conda env's `python3`.

| Entry point | Role | Required input | Interactive? |
|---|---|---|---|
| `scripts/hmm_finder.py` | Iterative family-discovery pipeline | `--fasta SEED` (prompts on a TTY if omitted) | only on a TTY |
| `scripts/scan_genome.py` | "Does THIS one genome carry my gene?" | `--seeds`/`--hmm` **and** `--genome`/`--accession` | never |
| `scripts/run_pipeline.py` | Zero-prompt launcher around the two above | same as the target it routes to | never (stdin detached) |
| `run.sh` | conda activation + tool check, then exec the chosen engine | passes through | inherits target |
| `start.sh` | Guided "run designer" that builds the command for you | none (walks you through it) | yes |

The full flag list for each script is in the [CLI flag reference](#complete-cli-flag-reference).

### `scripts/hmm_finder.py` — the main discovery pipeline

**Purpose.** End-to-end, unattended homologue discovery from a single seed FASTA. The only required argument is `--fasta`. It builds a profile HMM, six-frame-searches the configured public databases, extracts and ORF-validates the exact matched ORFs, re-seeds with the new hits, and iterates to convergence — then clusters, builds synteny figures, builds an ML tree, and assembles a labelled `PACKAGE/`.

**What it runs, in order** (see `main()`):
1. **Thread clamp.** `--cpu` is reduced to `os.cpu_count()` if it exceeds the available cores (prevents IQ-TREE's "more threads than CPU cores" abort and oversubscription); high-core hosts keep the full request.
2. **Engine check.** Resolves the search engine `DEPLOY` to the bundled `../engine` (when `engine/scripts/run_all_database_benchmark.py` exists) else the dated dev repo `~/Documents/HMM-Discovery-Deployable-20260602`; exits early if `BENCHMARK` is missing.
3. **Database resolution.** Precedence: explicit `--databases` > `--all-databases`/`--pick-databases` > default. With none given: on a TTY it **prompts** (numbered multi-select via `db_catalog.pick_databases`); unattended it silently uses the full default set (`DATABASES`, the 10 phage/viral DBs). `--smoke` overrides everything to 1 iteration on `INPHARED proteins`.
4. **Tool preflight** (`check_tools.ensure`, unless `--skip-tool-check`) — refuses to start a multi-hour run if required software is missing, pointing at `setup.sh`.
5. **Seed FASTA validation** — must exist, be non-empty, first non-blank line `>`, ≥1 record. Output redirects out of the engine repo if the chosen path would land inside it.
6. **Email resolution** (`--email` > `$NCBI_EMAIL` > TTY prompt > offline). No email forces `--no-annotate`; **the email is never hardcoded**. Only NCBI lookups (organism names, protein-DB hit sequences) are affected — local six-frame hits are unaffected. See [Offline behaviour](#offline-behaviour--no-annotate-and-the-email-gate).
7. **Nucleotide-seed translation** — `auto`/`nucleotide` input is translated with `--trans-table` (default 11) via `translate_seed`, which warns on internal stops (wrong code / not a clean CDS).
8. **Shared cache** — `--db-cache` (default `~/.cache/hmm-homologue-finder`) is reused across all runs via `cache`/`db_setup` symlinks, so each DB downloads once ever.
9. **Pre-run seed QC** (unless `--no-seed-tree`/`--smoke`) — a quick tree + alignment of just the input seeds (`build_tree_of_hits.py` → `seed_qc/`), to spot a mis-curated/outlier seed before committing.
10. **Iteration loop** (`--iterations`, default 3). Each round: runs `engine/scripts/run_all_database_benchmark.py` (HMM build + six-frame search, with `--keep-cache --max-synteny-genomes 200 --min-recovery 0.70 --skip-tree`), then `extract_validated_hits.py` (NT+AA+TSV+next-seed FASTA; `--prodigal-gate` appended only when set; `--email` only when set), then `annotate_organism.py` (adds the source-organism/phage-name column, unless `--no-annotate`). Completed rounds are skipped on re-run. **Convergence** uses the engine's `convergence_check` (hit-count change <5% AND ΔHMM-length <3); a later round failing HMM self-validation is **non-fatal** — it keeps the prior round and stops (only a run-1 self-recovery failure aborts, with a clear `--min-recovery` explanation).
11. **Canonical run selection** — `scripts/run_selection.py` is the single source of truth, called by **both** `hmm_finder._best_run_index` and `export_csv.py`. It picks the run recovering the most distinct homolog **loci**; a tie resolves to the **later** (converged) round, which is the model worth depositing. All downstream figures, the published HMM, the controls and the paper table therefore describe this same run.

    *Why one shared function:* the two modules previously used different rules — most hit **rows**, ties → earliest, in `hmm_finder`, versus most unique **sequences**, ties → latest, in `export_csv`. When two iterations tie on row count these disagree, and the package then shipped a `profile.hmm`, controls, tree and figures from one round beside tables and a report headline from another, so the quoted sensitivity/specificity did not describe the model that produced the reported hits.

12. **Homolog identity = genomic locus** — `run_selection.locus_ids` groups hits by canonical organism + strand + overlapping ORF interval, *not* by the amino-acid string. `aa_sequence` is the HMM **envelope** slice, so a refined model re-trims the same gene (one Erwinia gene appeared as a 63 aa, a 67 aa **and** a 118 aa "unique homolog") and the same physical genome appears under both a GenBank accession and its RefSeq `NC_` mirror. Both collapse to one locus. `hits_deduplicated.csv` then groups those loci by protein, so `n_loci` counts gene copies and `n_organisms` counts distinct phages.

    `n_databases` counts database **records** and is **not** evidence of corroboration — INPHARED redistributes RefSeq records, so one genome routinely appears in both. Use `n_organisms` for breadth claims. `max_domain_aa_len_any_round` records the longest envelope any round called, so a gene an earlier round called longer is never silently published short.
12. **Threshold-calibration controls** (unless `--no-controls`) — `run_controls` on the canonical HMM at strict=45 / moderate=30: positive control is the **seed self-recovery** sensitivity; negatives are composition-matched shuffled seeds (always available) plus optional UniProt unrelated-proteome sets (`--download-controls`, calibrated to `--biology-mode`, default `phage`). Writes `controls/control_report.json`, `controls_summary.csv`, `roc_curve.{png,svg,pdf}`. **Controls measure specificity vs unrelated proteomes, not a six-frame FDR.**
13. **Per-seed recovery QC** (`seed_recovery.seed_recovery_report`) — which named input seeds the initial (run1) vs final model recovers; surfaces divergent outliers by name → `seed_qc/seed_recovery.csv`.
14. **Read-through interrupted scan** (only with `--find-interrupted`, not in smoke) — `run_find_interrupted` re-scans the run's cached **nucleotide** DBs with stop-codons retained. Reporting threshold is `max(floor 30, ROC-Youden optimal)` (`_interrupted_min_bit`) — only ever tighter than the run's evidence bar. Writes `interrupted_homologs.tsv`, `..._domain_aa.faa`, `..._full_orf_aa.faa` (internal stops shown as `*`), `..._full_orf_nt.fna`, and tallies `overprinting_support` strong/partial. See [Interrupted / overprinted genes](#interrupted--overprinted-genes).
15. **Provenance + exports** — `write_methods_log` writes `run_manifest.json` + `METHODS.md` (tool versions, per-DB provenance, citations, the overprinting caveat). `write_csv_exports` (`export_csv.py`) writes `all_runs_hits.csv`, `hit_summary.csv`, `database_summary.csv`.
16. **Downstream** (skipped in smoke) — `cluster_and_clinker_corrected.py` (CD-HIT + clinker), `synteny_figure.py` (publication panels, `--color-by` function/conservation/both), `build_tree_of_hits.py` (ML tree of unique hits with the seeds marked), `run_perhit_hmm_alignment` (`hmmalign` → `hits_hmmalign.sto`/`.a2m`), per-run GFF3 (`write_gff3`), and `build_real_genbanks.py` (real-sequence GenBank neighbourhoods). Then `assemble_package` builds `PACKAGE/`, re-runs the CSV export into it, writes the one-page HTML report (`generate_report.py`), and `package_layout.write_readmes`. See [Downstream figures & trees](#downstream-figures--trees).

**Key flags.** `--fasta`, `--out-dir` (default `<fasta>_discovery/`), `--iterations` (3), `--cpu` (8), `--email`, `--input-type {auto,protein,nucleotide}`, `--trans-table` (11), `--no-annotate`, `--name`, `--databases`/`--all-databases`/`--pick-databases`/`--list-databases`, `--no-seed-tree`, `--synteny-gene-labels`, `--color-by {function,conservation,both}`, `--no-controls`, `--biology-mode {generic,phage,bacterial}`, `--download-controls`, `--smoke`, `--skip-tool-check`, **`--prodigal-gate` (OFF by default** — Prodigal overlap is informational, not exclusionary), `--find-interrupted`, `--db-cache`. Full table in the [CLI reference](#hmm_finderpy--familyhomolog-discovery).

**Outputs.** `run1/ run2/ …` (each `benchmark/validated/{hits.tsv,hits_aa.faa,hits_nt.fna,hits_unique_aa.faa,…}`), `downstream/` (clinker, synteny, tree, genbank), `controls/`, `seed_qc/`, `interrupted_homologs.tsv` (with `--find-interrupted`), `METHODS.md`, `run_manifest.json`, the three CSV exports, `pipeline.log`, and the self-contained `PACKAGE/`. Every output is detailed in the [Outputs catalog](#outputs-catalog--every-file-and-column).

**Worked command** (converged phage/overprinting run, the gp75 flagship pattern):
```
~/miniforge3/envs/hmm-discovery/bin/python3 scripts/hmm_finder.py \
    --fasta gp75_seed.faa --out-dir ~/hmm_runs/gp75_full \
    --all-databases --find-interrupted --iterations 2 --cpu 4 \
    --email you@inst.edu
```

### `scripts/scan_genome.py` — single-genome scan

**Purpose.** The focused counterpart to discovery: build (or accept) a profile HMM for one gene and answer "does **this one** genome contain it?" Reports a per-genome verdict — `GENE PRESENT` / `PRESENT but INTERRUPTED` / `GENE NOT DETECTED` — with each hit's genome coordinates, ORF validation, HMM score, and the hit protein + DNA. The scan is **read-through** (stop codons kept, not broken on), so a clean copy (0 internal stops) is found, and — only with `--find-interrupted` — a stop-interrupted/overprinted copy too, using the same silent-stop overprinting test as discovery.

**Inputs.** One HMM source (mutually exclusive, required): `--seeds FASTA` (protein or nucleotide CDS; built into an HMM via MAFFT L-INS-i ≤500 seqs else `--auto`, then `hmmbuild`) **or** `--hmm FILE` (existing profile). One genome source (mutually exclusive, required): `--genome PATH` (a `.fna/.fa[.gz]` FASTA, **or** an annotated GenBank `.gb/.gbk/.gbff` whose own gene names are used for the neighbours) **or** `--accession ACC` (NCBI nucleotide accession(s), comma-separated; needs `--email` or `$NCBI_EMAIL` — fetched as GenBank-with-parts so the genome's own annotation is available, falling back to FASTA-only).

**How it works** (`scan` → `_parse`): every contig is translated in all six frames over read-through windows (`find_interrupted._frames`/`_windows`); `hmmsearch --domtblout` scores the windows; each domain ≥`--min-bit` (default 25) is mapped back to genome coordinates (`FI.aa_to_nt`), the surrounding ORF is extended (`FI.extend_orf`), internal stops are counted, and interrupted hits are kept only with `--find-interrupted`. Overlapping-window duplicates and window-boundary truncations are de-duplicated (containment dedup keeps the full, higher-scoring hit). For interrupted hits it runs the overprinting test (`FI.analyze_overprinting`). The read-through machinery is shared verbatim with discovery — see [`find_interrupted.py`](#find_interruptedpy--read-through-detection-of-stop-interrupted-homologs).

**Reported columns** (`ROW_COLS`, written to `scan_hits.tsv`): `contig, strand, frame, nt_start, nt_end, domain_aa_len, internal_stops, status, domain_bit_score, i_evalue, orf_nt_start, orf_nt_end, orf_aa_len, has_start_M, ends_at_stop, overprinting_support, antisense_open_stops, stop_nt_positions, domain_aa, orf_aa, orf_nt`. Note `domain_aa_len` is the HMM-matched conserved region (the length to cite, ~137 for gp75); `orf_aa_len` is the surrounding Met-anchored six-frame ORF and **excludes the terminal stop** — for an overprinted domain on a stop-free antisense frame it can be long and `has_start_M=0` is correct.

**Neighbourhood + maps** (on by default; `--no-neighbours` to skip): `write_neighbourhoods` writes `scan_neighbourhood.csv` — the genes around each hit, ordered relative to your gene (`pos_index 0` = your gene; `relationship` flags upstream/downstream/overlapping). The `--flanks` genes each side are shown contiguously **and every overlapping gene is kept**, so an overprint partner (e.g. gp75's antisense RNA polymerase) is shown, not hidden. Names come from the genome's own annotation when available, else Prodigal (+ optional VOGDB VFAM via `--db-cache`). Genome-map figures `scan_genome_map_<hit>.{png,svg,…}` (windowed) and `..._whole.*` are drawn with `--map-tool` (default `dfv` = DNA Features Viewer); a locus GenBank is always written so you can open the map in Easyfig/Artemis/clinker yourself. Map options: `--no-gene-labels`, `--palette {default,colorblind}`, `--functional-labels`, `--module-brackets`.

**Other flags.** `--out` (default `genome_scan`), `--trans-table` (11, for a nucleotide seed), `--cpu` (4). **Exit code:** 0 if the gene was detected (clean or interrupted), 1 if absent — handy in shell scripts.

**Outputs.** `scan_hits.tsv`, `scan_hits_aa.faa`, `scan_hits_nt.fna`, `scan_report.txt`, `scan_neighbourhood.csv`, and the genome-map figures/GenBanks (see [single-genome scan outputs](#single-genome-scan-outputs-scan_genomepy)).

**Worked commands:**
```
# build the HMM from seeds and scan a local genome, including interrupted copies:
python3 scripts/scan_genome.py --seeds gp75_seeds.faa --genome phage.fna \
    --find-interrupted --out scan_out

# fetch the genome from NCBI by accession and scan with an existing HMM:
python3 scripts/scan_genome.py --hmm gp75.hmm --accession KX098390 \
    --email you@inst.edu --out scan_KX098390
```

### `scripts/run_pipeline.py` — the zero-prompt autonomous launcher

**Purpose.** A thin, no-prompt wrapper that runs either engine **fully autonomously, with zero interactive input** — the path to use from a script, cron, or `nohup`. It is not a re-implementation; it resolves a command and execs the real engine.

**Why it needs the env's python.** It must be invoked with the conda env's interpreter (`<env>/bin/python3 scripts/run_pipeline.py …`) because step 1 of `main()` prepends `os.path.dirname(sys.executable)` to `PATH`, making the env's tools (`hmmsearch`, `mafft`, `prodigal`, `iqtree`, …) discoverable **without `conda activate`**. The interpreter is its own pointer to the toolchain.

**Three things it does for autonomy:**
1. **PATH-injects** the running interpreter's bin dir (the tools come for free).
2. Runs the child with **`stdin=subprocess.DEVNULL`** so no step can ever block on `input()`.
3. **Injects no-prompt defaults** only when the user didn't already choose them, and forwards every other flag straight through.

**Launcher flags (consumed here, not forwarded):**
- `--scan` — route to `scan_genome.py` (single-genome mode) instead of `hmm_finder.py`.
- `--preset NAME` — apply a named default bundle (preset flags go first, so explicit user flags still win). Presets: `phage-discovery` (`--all-databases --find-interrupted`), `discovery` (`--all-databases`), `offline` (`--all-databases --no-annotate`), `smoke` (`--smoke`).
- `--list-presets` — print the presets and exit.
- `--dry-run` — print the exact resolved command (with injected defaults) and exit, running nothing.
- `-h`/`--help` — print the module docstring and exit (also printed when called with no args).

**Injected discovery defaults** (each only if absent): `--all-databases` (full set, no prompt — skipped if any of `--databases/--all-databases/--pick-databases/--list-databases/--smoke` is present); `--no-annotate` **unless** `--email` is given (the "offline-unless-email" rule — it stays fully offline and never sends a placeholder address to NCBI); and `--skip-tool-check`. In `--scan` mode nothing is injected (it validates that `--genome`/`--accession` and `--hmm`/`--seeds` are present; the scan already exits cleanly on a missing email).

**Validation.** Discovery requires `--fasta` and light-checks the file exists and begins with `>` (`_validate_fasta`); otherwise it exits with a clear `[run_pipeline]` message.

**`run_command.txt`.** On a real (non-dry) run, before exec, it writes the exact shlex-joined resolved command to `<out-dir>/run_command.txt` (`--out-dir`/`--out` is parsed by `_out_dir`) for provenance, creating the dir if needed. The launcher's exit code is the child's (0 on success, non-zero on failure).

**Worked commands:**
```
# autonomous family discovery (offline, full DB set, interrupted scan):
~/miniforge3/envs/hmm-discovery/bin/python3 scripts/run_pipeline.py \
    --fasta seed.faa --find-interrupted --out-dir ~/hmm_runs/myrun \
    --iterations 2 --cpu 4

# the same via a preset, and see exactly what it would run first:
~/.../bin/python3 scripts/run_pipeline.py --preset phage-discovery \
    --fasta seed.faa --out-dir ~/hmm_runs/myrun --dry-run

# single-genome scan, autonomous:
~/.../bin/python3 scripts/run_pipeline.py --scan --hmm gene.hmm \
    --accession NC_008720 --email you@inst.edu --out ~/scan1
```

### `run.sh` — conda activation + tool check, then exec the engine

**Purpose.** The terminal launcher for macOS/Linux/WSL2. It sources a conda **shell hook** (`etc/profile.d/conda.sh`, deliberately not `bin/activate`, which would inherit the launcher's `"$@"` and break `conda activate`), searching the common Miniforge/Mamba/(Ana)conda install locations; activates the `hmm-discovery` env, running `setup.sh` on first run if activation or `scripts/check_tools.py` fails. It then optionally wraps the run in `caffeinate -i` (macOS) or `systemd-inhibit --what=idle` (Linux; skipped under WSL where it's denied) to keep the machine awake during long runs.

**Routing.** `run.sh --scan …` execs `scripts/scan_genome.py "$@"`; otherwise it execs `scripts/hmm_finder.py "$@"`. It uses `exec` so the tool's exit code is not masked. All other flags pass straight through to the chosen engine.

**Worked command:**
```
bash run.sh --fasta seed.faa --name myrun --email you@inst.edu --iterations 3
bash run.sh --scan --seeds gp75_seeds.faa --genome phage.fna --find-interrupted
```

### `start.sh` — the guided "run designer"

**Purpose.** A friendly interactive wizard that builds the command for you and then execs `run.sh` (which handles conda activation). Best for first-time/exploratory use; it never runs anything you haven't confirmed.

**What it asks.** Step 1 — the **seed FASTA** (drag-and-drop quotes stripped; Enter accepts the bundled `examples/example_seeds.fasta`) and an optional genetic-code table for a nucleotide seed (`--trans-table`). Step 2 — the **run mode**: (1) smoke test, (2) full discovery run, or (3) **scan one genome**.

- **Mode 3 (scan)** asks for a local genome path or an NCBI accession (auto-detected: a real file → `--genome`, otherwise treated as an accession and an email is requested → `--accession --email`), an output dir, and whether to also report interrupted/overprinted copies (`--find-interrupted`), then execs `bash run.sh --scan --seeds … [--genome|--accession] …`.
- **Mode 2 (full run)** asks for output label/`--name`, NCBI email (blank = offline; a placeholder is never sent), output location, `--iterations`, `--cpu`, database choice (default set / `--pick-databases` / `--list-databases` first), strict `--prodigal-gate` (default N), `--find-interrupted` (default N), the seed QC tree (default Y; N → `--no-seed-tree`), synteny gene labels (`--synteny-gene-labels`), and `--color-by`.
- **Mode 1 (smoke)** adds `--smoke --cpu <ncpu>`.

It prints the assembled `bash run.sh …` command, a tip about `manage_cache.py`, and asks for a final confirmation before `exec bash run.sh "${ARGS[@]}"`. Default CPU is the host core count (`nproc`/`sysctl`).

**Worked command:**
```
bash start.sh        # then answer the prompts
```

---

## The search core: HMM iteration, six-frame translation, and controls

This section documents the discovery engine — the path from a seed FASTA to a converged set of profile-HMM homologs across protein and nucleotide databases, and the control machinery that calibrates the bit-score thresholds used to classify those hits. The code lives in `engine/scripts/run_all_database_benchmark.py` (the resumable benchmark driver) and the `engine/pipeline/` modules it imports: `hmm_builder.py`, `searcher.py`, `orf_prediction.py`, `hit_classifier.py`, `confidence.py`, `iterative.py`, and `controls.py`. The outer orchestrator `scripts/hmm_finder.py` wraps the engine in the iterate→validate→re-seed loop (see [the main pipeline](#scriptshmm_finderpy--the-main-discovery-pipeline)).

### Architecture: two layers

There are two layers that together implement "iterative profile-HMM search". Do not conflate them.

1. **The engine (`run_all_database_benchmark.py`)** runs *one* HMM against *all* databases. Given a seed FASTA it builds a single profile HMM (MAFFT → trimAl → `hmmbuild`), self-validates it, then for each registered database downloads it (cached), translates it if it is nucleotide (six-frame or Prodigal), runs `hmmsearch`/`hmmscan`, and aggregates hits. It is one pass — it does not re-seed itself.
2. **The orchestrator (`scripts/hmm_finder.py`)** drives the *iteration*: it calls the engine, ORF-validates the hits into a new seed set (`extract_validated_hits.py`), and feeds that expanded seed FASTA back into another engine run, repeating up to `--iterations` (default 3) or until convergence. This is the jackhmmer-style "rebuild and re-search" outer loop. The convergence rule and the candidate-promotion logic it relies on live in `engine/pipeline/iterative.py`.

### Step 1 — Seed FASTA to profile HMM (`Benchmark.build_core`, `hmm_builder.run_hmmbuild`)

`build_core()` (in `run_all_database_benchmark.py`) is the per-run HMM construction step:

- **Input summary** — `input_summary(copied)` parses the seed FASTA; `seq_count <= 0` aborts.
- **Alignment** — `accuracy_flags(seq_count)` picks the MAFFT strategy: L-INS-i (accurate, gold-standard for homologous domains) for tractable seed counts, FFT-NS-2/`--auto` only for very large sets. `run_mafft()` writes `alignments/seed.mafft.faa`.
- **Trimming** — `run_trimal()` produces `alignments/seed.mafft.trimmed.faa`; if trimAl fails the raw MAFFT alignment is used. `alignment_quality()` records columns/occupancy.
- **Build** — `run_hmmbuild(trimmed, hmm_dir/"benchmark_profile.hmm", "benchmark_profile")` shells out to `hmmbuild -n <name>`. It returns `{leng, alph, nseq, cksum, hmmbuild_version, name}`, cross-checking `hmmbuild` stdout against the `.hmm` header (`parse_hmm_file` reads NAME/LENG/ALPH/NSEQ/DATE/CKSUM). **`LENG`** — the number of match states — is the model length that drives both HMM-coverage scoring and the convergence test.

`hmm_builder.py` also exposes `parse_hmm_file()`, `logo_data()` (per-position match-emission probabilities, converting HMMER's −ln(p) log-odds back to linear probabilities via `_emissions_to_probs`, top-5 amino acids per node for sequence-logo rendering), and `self_search_recovery()` (below).

### Step 2 — Self-validation gate (`hmm_builder.self_search_recovery`)

Before any database is touched, `build_core()` calls `self_search_recovery(hmm, seed_faa)`. This runs `hmmsearch --tblout --noali` of the freshly built HMM **against its own seed sequences** and counts how many seeds score ≥ `strict_bits` (default 45.0). It returns `{recovered, total, recovery_rate, min_score, max_score}`.

The run **aborts** if `recovery_rate < --min-recovery`. The engine default is `0.95`; the orchestrator calls the engine with `--min-recovery 0.70` because iterative refinement deliberately broadens the family and an expanded seed set legitimately recovers a lower fraction. This is the gate whose failure message reads "The profile HMM failed self-validation… the seeds are NOT one coherent homologous family." Scientifically it is a sanity check that the seed set is one family before spending compute searching for its relatives. Note this is the **positive control** — the same seed-self-recovery concept that `controls.py` formalizes (see [Step 5](#step-5--threshold-calibration-and-controls-controlspy)).

### Step 3 — Searching protein databases (`searcher.run_hmmsearch_protein`, engine `run_protein_hmmsearch`)

For a protein database the engine streams the (optionally gzipped) FASTA into `hmmsearch --tblout -E <evalue> --cpu N --noali <hmm> -` (`Benchmark.run_protein_hmmsearch`). The standalone library equivalent is `searcher.run_hmmsearch_protein()`, which writes `<db>.tblout`, `<db>.domtblout`, and `<db>.out`, then returns `{tblout, domtblout, out, hit_count, strict_count}` where `strict_count` is the number of hits with `bit_score >= 45.0`.

`hmmscan`-mode databases (Pfam, VOGDB VFAM, PHROGs) take a different path: `run_pfam_hmmscan()` downloads/extracts/concatenates the HMM library, `hmmpress`es it, and scans the run's *own protein hits* against it for annotation — these are annotation layers, not discovery databases, and are excluded from the hit collection in `collect_hits()`.

**E-value** — the engine default is `--evalue 1e-5`; this is the `hmmsearch` inclusion threshold (lenient enough to expose divergent hits, which the downstream bit-score tiers then triage).

### Step 4 — Six-frame translation of nucleotide databases (`orf_prediction.predict_orfs_sixframe`, engine `run_nucleotide_hmmsearch` / `_write_sixframe_orfs`)

`hmmsearch` is a protein search, so nucleotide databases must be translated first. The mode is chosen by `--nt-orf-mode {sixframe,prodigal}`, **default `sixframe`**. Six-frame is the scientifically load-bearing choice: it is exhaustive and frame-agnostic, so it recovers genes that conventional annotation (Prodigal) misses — antisense, alternate-frame, and overprinted genes (the flagship gp75 case, an overprinted antisense ORF inside a virion RNA polymerase).

**Stop-to-stop ORFs.** The engine's `_write_sixframe_orfs` (the production translator) walks all 6 frames (`+`/`−` strands × 3 frames) codon by codon, splitting each frame at stop codons (`TAA/TAG/TGA`) into stop-bounded segments. A segment is emitted only if it is at least `min_nt = max(1, --min-orf-aa) * 3` long — i.e. ≥ `--min-orf-aa` amino acids (**default 30**, matching `MIN_AA` in `find_interrupted.py`). Each emitted ORF (`_emit_sixframe_orf`) is translated with `translate(to_stop=False)` and any internal `*` replaced with `X` (HMMER-friendly), and carries a header recording its provenance: `coords=<seqid>:<start>-<end>(<strand>) frame=<±N> nt_start= nt_end= aa_len=`. Those `coords=…` descriptions are the *only* place hit genomic coordinates are preserved (the all-ORF GFF is deliberately not retained — on large phage databases it can reach tens of GB), and they are re-parsed downstream by `collect_hits()` (`six = description.str.extract(r"coords=([^:]+):(\d+)-(\d+)\(([+-])\)")`) to populate `genome_id`, `seq_from`, `seq_to`, `strand`.

The library fallback `orf_prediction.predict_orfs_sixframe()` (used by `searcher.run_hmmsearch_nucleotide`) prefers an external `04_translate_sixframe.py` script but falls back to an inline Biopython translator that emits the same stop-split peptides with IDs `<recid>_s<strand>_f<frame>_o<n>`. Both honour `min_aa`.

**Cached translation.** Six-frame translation depends only on the DB file and the min-ORF length — **not** on the seed/HMM — so the engine persists ORFs next to the cached DB as `<dbfile>.sixframe.min<N>.faa` and reuses them across iterations and runs when `--keep-cache` is set and the cache is newer than the DB (`use_orf_cache`). This is what makes the orchestrator's 3-iteration loop affordable: the costly translation happens once and every later iteration just re-`hmmsearch`es the cached ORFs. To keep memory and disk bounded, large DBs are split with `seqkit split2 -p <cpu>` and translated/searched chunk-by-chunk, appending hits to the merged tblout. `preserve_synteny_context_records()` saves compact FASTA records for hit contigs before the raw DB is deleted, so synteny/genome-map context survives cache cleanup.

`--prodigal-gate` is **off by default**: Prodigal overlap is recorded as informational evidence (`prodigal_concordant`, `in_coding_locus` in the validated-hits table) but is **not** exclusionary, precisely because requiring overlap with a conventionally predicted gene would discard the antisense/alternate-frame homologs that are the whole point of the search.

### Building the hits table and confidence tiers (`hit_classifier.build_main_hits_table`, `confidence.classify_hits`)

`searcher.parse_tblout()` / `parse_domtblout()` parse HMMER output into DataFrames; the per-sequence tblout supplies `bit_score`, `evalue`, `bias_score`, and the per-domain domtblout supplies coordinates (`hmm_from/to`, `ali_from/to`, `env_from/to`) and `domain_bit_score`. `build_main_hits_table()` merges them — keeping only the best-scoring domain per protein — and calls `classify_hits()`.

**Confidence tiers (`confidence.score_hit`)** are assigned from bit score, HMM coverage, and bias, using the strict/moderate thresholds (defaults 45 / 30):

- `high_confidence` — `bit ≥ strict` AND `hmm_cov ≥ 0.60`. (Reciprocal validation and a domain match each *confirm* it but neither is required; if reciprocal validation ran and *failed* and there is no domain match, the hit is demoted to `putative`.)
- `putative` — `bit ≥ strict` but `0.30 ≤ hmm_cov < 0.60`.
- `divergent` — `bit ≥ strict` with `hmm_cov < 0.30`, or `moderate ≤ bit < strict`. These are the distant homologs the tool exists to surface.
- `likely_fp` — `bit < moderate`, or bias dominates the score (`bias > 0.80·bit` while `bit < strict`).

HMM coverage is `hmm_coverage_pct` from domtblout when present, else computed from `hmm_from/to` over `LENG`, else assumed 100% (so tblout-only protein searches are not penalized for missing domain coordinates). `add_qc_flags()` adds a pipe-separated `qc_flags` string: `high_bias` (`bias_score > 5.0`), `short_alignment` (`hmm_cov < 0.50`), `low_complexity` (one residue > 40% of the sequence), and `contig_edge` (a short < 80 aa protein butting a sequence end, or a hit within 30 nt of a real contig boundary when `contig_length` is known — a metagenome-fragment guard).

`hit_classifier.py` also provides `reciprocal_validate_hmmsearch()` (search hits back against the seed set, requiring `bit ≥ strict`), `extract_hit_sequences()`, `best_hit_per_genome()`, and `merge_all_databases()` (cross-database dedup keeping the highest-scoring copy of each `protein_id`).

### Iterative refinement and the convergence rule (`iterative.py`, orchestrator loop in `hmm_finder.py`)

The outer loop in `hmm_finder.py` (`for i in range(1, args.iterations + 1)`) does: engine search → `extract_validated_hits.py` writes `validated/hits_unique_aa.faa` (the next seed set) → that becomes the seed for run *i+1*. The engine rebuilds the HMM from the expanded seeds each round, so the model grows as the family is sampled more deeply — this is the "rebuild and re-search" mechanism.

`iterative.iteration_candidates(hits_df, seeds_faa, strict)` selects which hits are eligible to enter the next seed alignment: `confidence_tier == "high_confidence"` AND `bit_score ≥ strict`, excluding IDs already present in the current seed FASTA, sorted by bit score. `append_to_seeds()` concatenates the approved candidates onto the existing seeds.

**Convergence (`iterative.convergence_check`).** Refinement stops early when **both** conditions hold between consecutive rounds:

- hit-count change `< 5%` relative to the previous count (`abs(curr−prev)/max(prev,1) < 0.05`), and
- HMM-length change `< 3` match states (`abs(curr_leng − prev_leng) < 3`).

Both stabilizing means the family has been fully sampled — more rounds add neither new homologs nor new conserved columns. The orchestrator wires this in directly (`convergence_check(prev_n, n, prev_leng, curr_leng)`), counting `n` = unique validated seeds and reading `LENG` from each round's `benchmark_profile.hmm`; on convergence it records a `stop_reason` like `converged at run2 (hits 40→41, HMM length 137→138)`. The loop also stops if a round yields zero validated hits, and a *later* round whose expanded seeds fail the self-recovery gate is treated as the useful limit (keep the prior round, stop) rather than a fatal error. `convergence_data()` reshapes the per-iteration history (`iteration`, `hit_count`, `hmm_leng`, `diversity`) for the Plotly convergence chart.

### Step 5 — Threshold calibration and controls (`controls.py`)

The controls answer a specific question: **at the chosen bit-score threshold, how specific is this HMM against unrelated proteomes, and does the HMM recover its own seeds?** They are a *specificity* assay against biologically unrelated sequences — **not** a six-frame false-discovery-rate estimate (they do not model the multiple-testing burden of translating whole genomes in six frames).

**Control catalogue (`BUILTIN_CONTROLS`, `available_controls`).** Controls are mode-specific (`generic`/`phage`/`bacterial`):

- **Positive — `seeds_self_test`** (role `positive`): the seed sequences themselves, ~100% recovery expected. This is the seed self-recovery check — the positive control is recovering the input family, exactly mirroring the `self_search_recovery` gate in Step 2.
- **Negatives — unrelated proteomes:** in phage mode, `fungi_500.faa`, `mammalian_200.faa` (human/mouse housekeeping), `archaea_300.faa`; in bacterial mode, `euk_viral_200.faa`, `plant_400.faa`. These are curated Swiss-Prot sets fetched by `download_control_sequences()` from the UniProt `/search` endpoint (capped per taxon; `/stream` is deliberately avoided because it would return entire proteomes).
- **Universal negative — `shuffled_seeds`:** `generate_shuffled_control()` permutes each seed's amino acids (seeded RNG, default 42), preserving composition but destroying primary structure. A hit here means the score reflects composition/low-complexity, not homology.

**Running controls (`run_control_search`, `run_all_controls`).** Each control FASTA is searched with `hmmsearch --tblout --noali -E <evalue>` at a **lenient e-value of 1.0** so the *full score distribution* (not just passing hits) is captured. `_parse_tbl_simple` extracts full-sequence bit scores; `run_all_controls` tags each result with `strict_hits` (scores ≥ `strict_threshold`, default 45) and `moderate_hits` (≥ `moderate_threshold`, default 30) and returns a `ControlReport`.

**`ControlReport` metrics.**

- `sensitivity(t)` = fraction of *positive* sequences scoring ≥ t (divided by `n_seqs`, not by detected count).
- `specificity(t)` = fraction of *negative* sequences correctly rejected = `1 − FPR`; `false_positive_rate(t)` = `1 − specificity`.
- `summary()` reports `sensitivity`/`specificity`/`false_positive_rate` plus counts, at strict and moderate tiers, for the UI and `step_05_validate.py`.

**ROC / Youden calibration (`ControlReport.roc`).** This is **advisory only** — it reports whether the fixed strict tier (45 bits) sits near a data-driven optimum; it does **not** change the pipeline's strict/moderate tiers. It computes:

- **Exact AUC** via the Mann–Whitney U statistic over the positive-vs-negative bit-score vectors. `_score_vectors()` fills each control's undetected sequences (those below the lenient e-value, hence absent from tblout) with a `-inf` sentinel so they rank below every real detection and the AUC matches how sensitivity/specificity divide by `n_seqs`.
- **Youden's J optimal cutoff** over candidate thresholds built from the detected (finite) scores plus a finite 0-bit noise floor, midpoints, and bounds. When several cutoffs tie for max J (a clean separation gap), it picks the most robust — the maximum-margin threshold (mid-gap), ties resolving toward the lower value.
- Returned keys: `auc`, `optimal_threshold`, `youden_j`, `sensitivity_at_optimum`, `specificity_at_optimum`, `separable`, `strict_threshold`, and crucially **`optimal_threshold_defined`**.

**`optimal_threshold_defined`** is `True` only when at least one negative sequence scored above the noise floor (`n_negative_detected > 0`). When **no** negative is detected — the common, good outcome of a clean, specific HMM — the FPR is 0 at every finite cutoff, so a whole range of thresholds ties at J = 1 and the reported "optimum" is just the arbitrary midpoint between the 0-bit floor and the weakest positive. It carries **no information from the negative distribution**, so the flag is set `False` and `threshold_method` annotates it as "optimum undefined — no negative detected… advisory only." This prevents an undefined optimum from being mistaken for a calibrated one. The operative reporting filter downstream is `max(30, optimal_threshold)` (in `extract_validated_hits.py`), which is unaffected by an undefined optimum.

`plot_roc()` renders `roc_curve.{png,svg,pdf}` (house style, strict-45 and Youden-optimum points marked); `to_dataframe()` / `to_json()` / `summary()` serialize the report. The on-disk products are catalogued under [`controls/`](#controls--threshold-calibration-at-the-out-dir-root).

### `domain_aa_len` vs `orf_aa_len` — what length to cite

These two lengths describe different things and are produced for every six-frame hit (in `find_interrupted.py` and `extract_validated_hits.py`). Citing the wrong one materially misstates the homolog.

- **`domain_aa_len`** — the length of the **HMM-matched conserved region** (the alignment envelope, `env_to − env_from + 1`). This is the reliable, profile-anchored homolog length and the number to cite (≈137 aa for gp75). It is independent of how far the surrounding ORF extends and of whether a start codon is present.
- **`orf_aa_len`** — the length of the **surrounding gene ORF**, reported from its Met start codon (`met_anchor`), not from the upstream stop. The Met is the one closest to and within `MET_MARGIN = 60 aa` upstream of the domain envelope (and after the preceding in-frame stop); falling back to the domain start when no Met is near. In `find_interrupted.py`, `orf_aa_len` excludes the terminal stop codon (`(orf_to − orf_from + 1) − int(terminal_stop)`), since a `*` is not a residue — matching `scan_genome.py` and `extract_validated_hits.py`.

For an **overprinted** domain sitting on a long, stop-free antisense frame, `orf_aa_len` can be very large and **`has_start_M = False` is correct** — there is no in-frame Met within `MET_MARGIN` of the domain because the conserved region is embedded in an open reading stretch of the host gene's opposite strand, not a conventionally Met-initiated ORF. `has_start_M = False` here is a true property of an overprinted gene, not a defect. Hence `domain_aa_len` (the conserved envelope) is the meaningful length, while `orf_aa_len`/`nt_start`/`nt_end` describe the stop-bounded search ORF used for synteny and genome maps.

### Key files and flags (search core)

- `engine/scripts/run_all_database_benchmark.py` — the resumable engine. CLI: `--fasta`, `--out`, `--preset {all,partial,smoke}`, `--databases`, `--dry-run`, `--force`, `--keep-cache`, `--skip-tree`, `--cpu`, `--nt-orf-mode {sixframe,prodigal}` (default sixframe), `--min-orf-aa` (default 30), `--evalue` (default 1e-5), `--strict-bits` (default 45.0), `--min-recovery` (default 0.95), `--strict-databases` / `--stop-on-optional-failure` (a single DB failure is non-fatal by default — skip and continue, aborting only if every DB fails).
- `engine/pipeline/hmm_builder.py` — `run_hmmbuild`, `parse_hmm_file`, `logo_data`, `self_search_recovery`.
- `engine/pipeline/searcher.py` — `run_hmmsearch_protein`, `run_hmmsearch_nucleotide`, `parse_tblout`, `parse_domtblout`.
- `engine/pipeline/orf_prediction.py` — `predict_orfs_sixframe`, `predict_orfs_prodigal`, `choose_and_predict`.
- `engine/pipeline/confidence.py` — `score_hit`, `classify_hits`, `add_qc_flags`.
- `engine/pipeline/hit_classifier.py` — `build_main_hits_table`, `reciprocal_validate_hmmsearch`, `merge_all_databases`.
- `engine/pipeline/iterative.py` — `iteration_candidates`, `convergence_check`, `convergence_data`, `append_to_seeds`.
- `engine/pipeline/controls.py` — `BUILTIN_CONTROLS`, `available_controls`, `generate_shuffled_control`, `run_control_search`, `run_all_controls`, `ControlReport` (`.roc`, `.sensitivity`, `.specificity`, `.summary`, `.plot_roc`), `download_control_sequences`.
- `scripts/hmm_finder.py` — the iterate→validate→re-seed orchestrator (`--iterations` default 3; reuses `convergence_check`).

---

## Interrupted / overprinted genes

The standard search builds **stop-to-stop** six-frame ORFs, so any homolog whose gene carries an internal stop codon has its ORF truncated *at* the stop and is split or missed. This is exactly the failure mode of an **overprinted** gene: when a small gene is encoded antisense to (or in a different frame from) a large essential gene — the flagship case is **gp75 overprinted on a virion DNA-directed RNA polymerase** — a point mutation can be a **premature stop in the small-gene frame yet silent (synonymous) in the polymerase frame**. Selection keeps the polymerase intact and tolerates the truncation, so the small gene is interrupted but invisible to a stop-to-stop search.

Two scripts implement the machinery:

- `scripts/find_interrupted.py` — **read-through translation** that keeps stops, searches them with the family HMM, locates the internal stop(s), reconstructs the full ORF, maps everything to genome coordinates, and runs the overprinting (silent-stop) test. It is both a standalone CLI and the function library imported (as `FI`) by `scan_genome.py` and called (as `_run`) by the benchmark wrapper `run_find_interrupted` in `hmm_finder.py`.
- `scripts/extract_validated_hits.py` — the authoritative, ORF-validated **extractor** for a completed benchmark run. It reconstructs the exact six-frame ORF the HMM matched, anchors the gene at its Met start, slices the domain, runs a 3-pass validation, and decides `passes_orf_filter` (with the optional `--prodigal-gate`).

### `find_interrupted.py` — read-through detection of stop-interrupted homologs

CLI: `python3 find_interrupted.py --genomes genomes.fa[.gz] --hmm profile.hmm --out interrupted_homologs.tsv [--min-bit 25] [--cpu 8]`. `main()` validates that both `--genomes` and `--hmm` exist, then calls `_run`. Module-level constants: `STOP_SEARCH = "X"` (HMMER-friendly placeholder for a stop in the SEARCH sequence), `MIN_AA = 30` (frames shorter than this are ignored — matches the six-frame floor), and `MET_MARGIN = 60` (a start codon may sit at most 60 aa upstream of the domain envelope).

#### Read-through translation: `read_through_aa(nt, strand, frame)`
Translates one reading frame **keeping stops**. It uppercases, replaces `U`→`T`, reverse-complements for the `-` strand, slices off the frame offset, trims to whole codons, and translates. It returns a `(search_aa, marker_aa)` tuple from the *same* translation:
- `marker_aa` keeps `*` at every stop position — used afterward to **locate** stops and to extend the ORF.
- `search_aa` is `marker_aa` with `*`→`X` — fed to `hmmsearch` so the profile matches **across** the stop instead of stopping at it.

`_frames(seq)` yields `(strand, frame, search_aa, marker_aa)` for all six frames, skipping frames shorter than `MIN_AA`.

#### Windowing long frames
HMMER aborts (SIGABRT) on genome-length sequences (10k–100k+ aa), so read-through frames are windowed. `_windows(n)` yields `(offset, length)` windows of `WIN_AA = 5000` aa stepping by `WIN_STEP = 4400`. The overlap (`WIN_AA - WIN_STEP = 600` aa) is far larger than any domain, so every domain lies fully inside at least one window. Window offsets are carried in the FASTA name (`contig__SF__windowOffset`) so frame-relative coordinates can be made absolute later. Because overlapping windows can report the same domain twice, `_run` deduplicates by `(contig, strand, frame, stop_aa_positions)`, keeping the highest `domain_bit_score`.

#### Counting internal stops: `count_envelope_stops(marker_aa, env_from, env_to)`
Given the 1-based inclusive HMM envelope `[env_from, env_to]`, it slices `marker_aa[env_from-1:env_to]`, drops a trailing `*` (a stop that merely *ends* the envelope is not "internal"), and returns `(count, [1-based positions])` of the remaining `*`. **A match is reported as interrupted only if this count is ≥ 1.** This is the gate that distinguishes an interrupted homolog from a clean one.

#### Reconstructing the full ORF: `extend_orf(marker_aa, env_from, env_to)`
Extends the domain envelope to the flanking stops *in the read-through frame* so the sequence **after** the premature stop (the rest of the gene, down to its natural stop) is captured. Returns `(orf_from, orf_to, full_aa)`, 1-based inclusive, with `*` kept at every stop. Key behavior:
- **Start anchoring (MET_MARGIN).** It does **not** start the ORF at the upstream in-frame stop — that would prepend the residues between the upstream stop and the real ATG (a 145-aa stop-to-stop segment for a 138-aa gene). Instead it takes `lo = max(after_stop, env_from-1-MET_MARGIN)` and finds the Met **closest to and within 60 aa upstream** of the envelope start (`marker_aa.rfind("M", ...)`). Critically, it does **not** take the *earliest* Met — on a long stop-free antisense/overprint frame that Met can be hundreds of residues upstream and is not the gene start. If no Met is near, it falls back to the domain start (`orf_from = env_from`, no upstream extension).
- **3' end.** `right = marker_aa.find("*", env_to)` finds the first stop at/after the domain; `orf_to` includes that **natural terminal stop** (or runs to the window edge if none).

#### Coordinate mapping back to the genome
- `aa_to_nt(contig_seq, strand, frame, aa_from, aa_to)` maps a frame-relative aa span (1-based inclusive) back to **forward-strand** contig coordinates, returning `(fwd_start, fwd_end, coding_nt)`. For `-` hits the returned `coding_nt` is reverse-complemented to the coding 5'→3' sequence. Used for both the domain (`domain_nt_*`, `domain_nt`) and the full ORF (`orf_nt_*`, `full_orf_nt`).
- `stop_nt(contig_seq, strand, frame, aa_pos)` returns the forward-strand 1-based coordinate of the **first base** of the stop codon at a given frame-relative aa position — used to populate `stop_nt_positions` and `natural_stop_nt`.

#### Overprinting / silent-stop analysis — the *proof*, not just the location
Locating a premature stop only tells you the gene is interrupted. The overprinting verdict tests whether that stop is the **signature of an overprinted gene**: a nonsense mutation in the small gene that is **silent (synonymous)** in an **open** overlapping antisense frame, so selection on the antisense gene can tolerate (or favor) it. The genetic code (table 11, bacterial/phage) is used throughout. Helper layout:

- `_aa(codon, table=11)` — cached single-codon translation; ambiguous/short codons → `X`. `_revcomp`, `_comp` use a translation table for speed.
- `_codon_covering(contig, strand, frame, fwd_i0)` — returns the codon (read 5'→3' in the gene) and the 0-based index within it that cover a given forward 0-based position, working on a local 3-base slice (no whole-contig reverse-complement). Returns `None` if the position falls outside a complete codon.
- `_frame_stop_count(contig, strand, frame, lo1, hi1)` — counts stop codons of one reading frame whose codon lies **fully** within forward window `[lo1, hi1]`. A count of **0 means the frame is OPEN** across the domain — i.e. a candidate overlapping ORF.
- `_stop_silent_in_frame(contig, small_strand, stop_fwd1, anti, frame, table=11)` — the core synonymy test for one premature stop. It confirms the codon is a stop in the small frame, then enumerates **every single-base substitution** that *removes* that stop in the small gene; for each such stop-removing change it checks whether the corresponding base change leaves the **antisense** frame's amino acid unchanged (and that the antisense codon actually encodes a residue, not a stop). Returns `True` if **any** stop-removing substitution is synonymous antisense — i.e. the nonsense mutation can be reverted to sense in the small gene without changing the antisense protein.
- `analyze_overprinting(contig, small_strand, dom_lo1, dom_hi1, stop_fwds, table=11)` — orchestrates the verdict per interrupted homolog. It picks the antisense frame (`anti = opposite strand`, frames 0/1/2) with the **fewest** stops across the domain via `_frame_stop_count` (the candidate overprinted ORF), then runs `_stop_silent_in_frame` for each premature stop. Returns `{open_frame, open_stops, per_stop_silent, support}`.

**Verdict semantics (`overprinting_support`):**
- `strong` — the chosen antisense frame is **fully open** (`open_stops == 0`) across the whole domain **AND every** premature stop is silent in it. This is the discriminating signal: a single silent stop is ~85–100% expected by genetic-code geometry and weak alone, but an antisense frame staying **OPEN across an entire ~137-aa domain occurs only ~0.7% by chance**. `strong` is therefore a **necessary SEQUENCE signature** of an overprinted gene — **not proof of expression**.
- `partial` — at least one premature stop is silent antisense (`n_sil > 0`) but the criteria for `strong` are not all met.
- `none` — no silent-antisense evidence.

#### Output schema (`interrupted_homologs.tsv`, `ROW_COLS`)
One row per dedup'd interrupted match, sorted by descending `domain_bit_score`. Columns are catalogued in full under [`interrupted_homologs.tsv` in the outputs catalog](#interrupted_homologstsv-and-its-fastas---find-interrupted). The salient fields: **`domain_aa_len`** = `env_to - env_from + 1`, the HMM-matched conserved region (~137 for gp75) — **the homolog length to cite**; the overprinting block `overprinting_support`, `antisense_open_frame`, `antisense_open_stops`, `stop_silent_antisense` (per-stop `1`/`0`); and **`orf_aa_len`** = `(orf_to - orf_from + 1) - terminal_stop` — the surrounding Met-anchored read-through ORF, which **EXCLUDES the terminal stop** and can be long for an overprinted domain on a stop-free antisense frame.

#### Sequence FASTAs emitted beside the TSV
- `write_aa_fastas(rows, out_tsv)` → `<stem>_domain_aa.faa` (the HMM-matched domain, `*` at every internal stop) and `<stem>_full_orf_aa.faa` (the full read-through ORF; premature stops kept `*`, terminal `*` = the natural gene end).
- `write_orf_nt_fasta(rows, out_tsv)` → `<stem>_full_orf_nt.fna` (coding 5'→3', ending in the actual stop-codon triplet; translates back to `full_orf_aa`).

All three are written even when empty so report/PACKAGE links never dangle. `_fasta_header(r)` builds a unique header from contig, strand/frame, domain coordinates, stop count/positions, bit score and i-evalue.

#### Batched execution: `_run(...)` and `_search_batch(...)`
`_run` streams the nucleotide DB in `BATCH_CONTIGS = 500`-contig chunks so a huge DB never produces one giant temp file (which aborts HMMER); temp lives next to the output, not in `/tmp`. Each batch writes a `batch.faa`, runs `hmmsearch --noali --cpu N --domtblout`, and `_search_batch` parses the domtbl: it reads `dom_bits` (col 14, p[13]), `i_eval` (col 13, p[12]), and `env_from`/`env_to` (cols 20/21, p[19]/p[20]); applies the `--min-bit` floor; calls `count_envelope_stops` and **skips any match with 0 internal stops** (clean copies belong to the standard search). For the survivors it builds the full row, computes `aa_before`/`aa_after`, maps domain and ORF coordinates, and runs `analyze_overprinting`. Returns a summary dict (`matches_scored`, `interrupted_candidates`, output paths).

#### Pipeline integration (threshold tightening)
The benchmark wrapper `run_find_interrupted` in `hmm_finder.py` (opt-in via `--find-interrupted`) re-scans only this run's cached nucleotide DBs. Its reporting floor comes from `_interrupted_min_bit(control_summary, floor=30.0)`: the read-through scan covers a much larger, noisier space (every frame, stops retained, windowed) than the stop-to-stop search, so the bar **only ever tightens** — it is `max(30, ROC-Youden optimal_threshold)` when controls were run, and the bare 30-bit floor with `--no-controls`. The summary tallies `overprinting_strong`/`overprinting_partial`, feeds `METHODS.md`/`run_manifest.json` (`interrupted_homologs` block), and the TSV + three FASTAs are copied into the PACKAGE (`tables/` and `sequences/`). `scan_genome.py` reuses the same library (`import find_interrupted as FI`) for single-genome scans, emitting a `status` column of `clean`/`interrupted` and (for interrupted rows) `overprinting_support` and `antisense_open_stops`.

### `extract_validated_hits.py` — ORF-validated extraction from a completed run

This script is the **fix** for a sequence-capture defect: the benchmark found homologs by six-frame search but deleted the translated FASTA, and a later synteny step stored a *Prodigal-predicted* gene (a different gene that merely overlaps the locus) as the "hit protein," with an additional 1-based/0-based slicing bug that shifted frames and inserted stops. This extractor reconstructs the **exact** ORF the HMM matched, frame-correctly, from the recorded genomic coordinates.

CLI: `--results-dir` (must contain `hits_main.tsv` and `synteny_context_cache/<genome>.fna`), `--hmm`, `--run-label A|B|C`, `--out-dir`, optional `--email` (NCBI Entrez, for protein-DB hits — never assumed), `--cpu`, and `--prodigal-gate`. `ensure_env_on_path()` (from `env_paths`) puts the conda env's `hmmsearch`/`prodigal` on PATH even when invoked outside `conda activate`.

#### ORF reconstruction: `reconstruct_orf(genome_seq, start_1, end_1, strand)`
Mirrors the pipeline's `_emit_sixframe_orf` convention **exactly**: coordinates are always forward-strand; `sub = genome_seq[start_1-1:end_1]` (the corrected 1-based→0-based slice that the earlier code got wrong), reverse-complemented for `-` hits, trimmed to whole codons, translated, and the **terminal stop only** is stripped (`rstrip("*")`). Returns `(nt, aa)` — stop-free by construction for a genuine six-frame ORF.

#### Met start anchoring: `met_anchor(orf_aa, env_from, margin=MET_MARGIN)`
`MET_MARGIN = 60`. A six-frame ORF is stop-bounded, so its first residues lie between the upstream in-frame stop and the real ATG; quoting the whole stop-to-stop segment over-states the gene length. `met_anchor` finds the Met **closest to and within 60 aa upstream** of the domain envelope start (`orf_aa.rfind("M", max(0, env_from-margin), env_from)`) and returns `(met_offset, gene_aa)`. **Critical for overprinted families:** when a small domain rides a long stop-free antisense frame, the nearest Met before the envelope can be hundreds of residues upstream and is **not** the gene start, so the Met is *required within margin*; with none found it anchors at the domain start (no upstream extension). In that case `has_start_M` is **False, which is CORRECT**, and `domain_aa_len` (the HMM envelope) is the reliable conserved length to cite.

#### The 3-pass validation (in `main()`)
- **Pass 1 — reconstruct every ORF.** For each row in `hits_main.tsv`, parse `coords=contig:start-end(±)` from the `description`, load the cached `<contig>.fna` (or `<genome_id>.fna`), call `reconstruct_orf`, drop ORFs < 20 aa, write `orfs_aa.faa`/`orfs_nt.fna`, and record provenance (`source_url`, `source_sha256`, `accessed_at`), coordinates, frame, and `orf_aa_len`. `source_type = "six_frame_orf"`.
- **Pass 1b — protein-database hits** (rows with no `coords=`). These are already-annotated proteins (SwissProt, RefSeq/INPHARED proteins, Pfam): the matched target **is** the sequence, so there is no ORF to reconstruct. `fetch_protein_seqs(accessions, email)` retrieves them by accession via Entrez (batches of 40, retries, `clean_accession` handles bare and `sp|ACC|NAME`/`tr|...` ids). With no `--email` the hits are kept but **not** fetched (NCBI policy forbids a placeholder address); any hit whose sequence can't be fetched is written to `protein_hits_unfetched.tsv` rather than silently dropped. `source_type = "annotated_protein"`.
- **Pass 2 — domain envelopes.** `domain_envelopes(orf_aa_path, hmm, cpu)` runs `hmmsearch --domtblout -E 1e-3 --domE 1e-3` on the reconstructed ORF set and maps each ORF id → `(env_from, env_to)` of its best (lowest-evalue) domain.
- **Pass 3 — per-hit validation + domain slicing.** Using the envelope `(a, b)`: `dom_aa = orf_aa[a-1:b]`, `dom_nt = orf_nt[(a-1)*3:b*3]`; `internal_stops = dom_aa[:-1].count("*")`; `domain_aa_len = len(dom_aa)`; `domain_coverage = round(len(dom_aa)/len(orf_aa), 3)`. For **six-frame** hits it then calls `met_anchor` to report the Met→stop gene (`orf_aa_len`, `has_start_M`, adjusted `orf_nt_start`/`orf_nt_end`), sets `ends_at_stop = True` (stop-bounded by construction), and records Prodigal overlap.

#### Prodigal overlap (informational; gate only with `--prodigal-gate`)
`prodigal_genes(fna)` runs `prodigal -p meta` (cached per genome) and `prodigal_overlap(genes, start, end, strand)` returns `(same_strand_pct, any_strand_pct)`:
- `prodigal_same_strand_pct` — best overlap with a gene on the **same** strand; typically ~0 for this family because the homologs are antisense/alternate-frame to predicted genes. That **discordance is the novelty signal** (they were missed by standard annotation).
- `prodigal_any_strand_pct` — best overlap with **any** predicted gene; typically high, confirming the hit sits in a real coding-dense locus (not six-frame noise).

These set `prodigal_concordant` (`same ≥ 0.50`) and `in_coding_locus` (`any ≥ 0.50`), both **informational**.

#### `passes_orf_filter`
- Annotated-protein hits: pass iff `internal_stops == 0` (genuine proteins by definition; Prodigal/ORF checks are `NA`).
- Six-frame hits: `passes = (internal_stops == 0)`. **`--prodigal-gate` is OFF by default** — Prodigal overlap is informational, not exclusionary, because requiring overlap with a predicted gene would discard exactly the genuine antisense/alternate-frame homologs this tool is built to find. With `--prodigal-gate`, the filter additionally requires `in_coding_locus` (`prodigal_any_strand_pct ≥ 50%`) for a stricter, higher-specificity set. `domain_coverage` is likewise reported but never exclusionary.

#### Outputs
`hits.tsv` (one fully-annotated row per hit; column order fixed by `COLS` — catalogued [below](#hitstsv--per-run-hitscsv--the-per-hit-evidence-table-36-base-columns)), `hits_aa.faa`/`hits_nt.fna` (domain AA/NT per hit), `orfs_aa.faa`/`orfs_nt.fna` (full ORF context), and `hits_unique_aa.faa` — deduplicated, ORF-validated (`passes_orf_filter`) domain AA sequences that **seed the next iteration** of the search. The summary line reports total hits (six-frame vs protein-DB), `pass_filter`, `in_coding_locus`, internal-stop count, and unique seeds.

### When read-through earns its keep, and what the verdict means

`find_interrupted.py` is worth its extra full-DB scan precisely when the target is a **distant homolog carrying a mid-gene premature stop** — most importantly an **overprinted gene** where the nonsense mutation is silent in an overlapping frame. The normal stop-to-stop search structurally cannot recover these: it breaks the ORF at the stop. Read-through keeps the stop (as `X` for the HMM, `*` for the marker), so the conserved domain is searchable end-to-end, the truncation is located, and the overprinting test can ask the decisive question.

The verdict is graded and deliberately conservative. `none`/`partial` reflect weak or absent antisense synonymy. `strong` is the necessary **sequence** signature of overprinting — a fully open antisense frame across the whole domain (rare by chance) with every premature stop synonymous there — but it is evidence of a *constraint*, **not proof the small gene is expressed**. The two extractors are complementary: `extract_validated_hits.py` curates the clean, validated, ORF-grounded homolog set (and the next seed), while `find_interrupted.py` recovers and characterizes the interrupted/overprinted copies that the clean search, by design, leaves behind.

---

## Annotation & organism naming

This section covers the components that put *human-readable biology* onto raw hits: functional names on neighbourhood genes (`annotate_genes.py`, VOGDB VFAM), source-organism (phage) names on each hit's genome (`annotate_organism.py`, NCBI Entrez), and the single rule (`canonical.py`) that collapses the same phage appearing under different accessions so organism counts are honest. It also documents the offline path (`--no-annotate`) and the VOGDB-category → broad-category mapping consumed by the figure scripts.

All paths below are under `hmm-homologue-finder/scripts/`.

### `annotate_genes.py` — functional annotation via VOGDB VFAM (hmmscan)

**Purpose.** Assigns real functional names ("major capsid protein", "terminase large subunit", "DNA-directed RNA polymerase", …) to the *neighbourhood* genes drawn around the gene of interest in the publication synteny figures and genome maps. It is **not** applied to the homolog hits themselves — it annotates the flanking CDS in each locus so the synteny figure can colour genes by function and bracket conserved modules. It is consumed by `synteny_figure.py` (imported at module top as `import annotate_genes`) and, transitively, by `genome_map.py` (which reuses the synteny colour scheme and `categorize`).

**Database.** VOGDB VFAM, release **vog230** (VOGDB release 230). Two artefacts are fetched from the Vienna fileshare:
- `VFAM_HMM_URL` = `https://fileshare.csb.univie.ac.at/vog/vog230/vfam.hmm.tar.gz` — the per-VFAM profile-HMM library.
- `VFAM_ANN_URL` = `https://fileshare.csb.univie.ac.at/vog/vog230/vfam.annotations.tsv.gz` — the consensus function + functional category per VFAM.

**Cache layout & setup (one-time, idempotent).** `vogdb_dir(cache)` resolves to `<cache>/annotation/vogdb`. `setup(cache)` (called by the figure scripts before annotating):
1. Downloads + `gunzip -kf` the annotations TSV.
2. Downloads + `tar -xzf` the HMM tarball, then concatenates every per-VFAM `*.hmm` (excluding `vfam_all.hmm` itself) into a single `vfam_all.hmm`.
3. `hmmpress -f vfam_all.hmm` to build the binary index (`.h3m/.h3i/.h3f/.h3p`).
4. **Reclaims space**: deletes `vfam.hmm.tar.gz`, `vfam.annotations.tsv.gz`, and the per-VFAM `hmm/` subdirectory once pressed.
5. Writes `provenance.json` recording `database`, `release` (vog230), both source URLs, and `tool: HMMER hmmscan` — so the figure's functional annotation is citable.

`is_ready(cache)` is the readiness gate: it returns `True` only when **both** `vfam_all.hmm.h3m` and `vfam.annotations.tsv` exist in the cache. `setup()` returns the value of `is_ready()` after running, and prints `(VOGDB setup failed: …)` (returning `False`) on any download/index error — a failed or absent VOGDB never aborts the run; genes simply stay "hypothetical protein".

**Annotation call.** `annotate(proteins, cache, cpu=4, evalue=1e-3)` takes `proteins` as `{id: aa_sequence}` and returns `{id: {'vfam', 'function', 'category'}}` only for proteins that hit. Mechanics:
- Short-circuits to `{}` if `proteins` is empty or `is_ready(cache)` is False.
- Loads the annotation table via `load_annotations(cache)` (below).
- Writes the query proteins to a temp FASTA and runs `hmmscan --tblout <out.tbl> -E 1e-3 --cpu <cpu> --noali vfam_all.hmm q.faa`. Any hmmscan failure returns `{}` (best-effort; never raises into the figure).
- Parses the `--tblout`, keeps the **best (lowest) e-value VFAM per query protein** (`best[q]`), then joins to the annotation table. Missing-annotation VFAMs default to `function="hypothetical protein"`, `category=""`.

**`load_annotations(cache)`** parses `vfam.annotations.tsv` into `{VFAM_id: {"function": <cleaned col 5>, "category": <col 4>}}`. Column 4 is the VOGDB **FunctionalCategory** code; column 5 is the consensus description. `_clean_desc()` normalises the description: strips a leading `sp|…|…` / `tr|…|…` UniProt prefix, removes leading qualifiers `Putative|Predicted|Probable|Uncharacterized`, and substitutes `"hypothetical protein"` for an empty result.

**VOGDB category → broad-category mapping (where it lands).** The raw VOGDB FunctionalCategory single-letter code travels with each gene as `vog_cat` and is used **only as a fallback** by `synteny_figure.categorize(func, vog_cat="")`. The mapping `_VOG_CAT` in `synteny_figure.py` is deliberately conservative:

| VOGDB code | Broad category |
|---|---|
| `Xs` | structural |
| `Xr` | replication / nucleotide metabolism |
| `Xh` | host / metabolism / defense |

`Xp`/`Xu` (poorly characterised / unknown) and ambiguous multi-letter combos (e.g. `XhXpXrXs`) are intentionally **not** mapped and stay `"hypothetical / unknown"` (`HYPO_CAT`). Resolution order in `categorize()`: first match against the keyword rules in `CATEGORY_RULES` on the (lower-cased) function description — unless the description itself contains `hypothetical`/`uncharacterized`/`unknown`/`duf`, in which case keyword matching is skipped — and only if no keyword matches does it fall back to `_VOG_CAT`. The eight broad categories (`structural`, `packaging`, `replication / nucleotide metabolism`, `transcription / regulation`, `lysis`, `integration / mobile`, `host / metabolism / defense`, `RNA gene`, plus `hypothetical / unknown`) have fixed colours in `CATEGORY_COLORS` (and a Paul-Tol colour-blind palette `CATEGORY_COLORS_CB`, selectable per run with `--palette colorblind`), shared between synteny figures and genome maps for consistency. Note the deliberate ordering of `CATEGORY_RULES`: `transcription / regulation` (matching `rna polymerase`) is checked **before** `replication / nucleotide metabolism`, so a virion **DNA-directed RNA polymerase** is binned as transcription, not replication — directly relevant to the gp75 flagship, whose host gene is exactly that virion RNAP.

### `annotate_organism.py` — source-organism (phage) name column

**Purpose.** Adds an `organism` column to a `hits.tsv`, inserted **immediately after `genome_id`**, naming the source phage for each hit's genome (e.g. `Escherichia phage vB_EcoP_G7C`). Invoked by the iterative driver in `hmm_finder.py` once per round, after ORF-validated extraction:
```
python3 annotate_organism.py --hits-tsv <validated>/hits.tsv --email <email>
```
(only when `--no-annotate` is not set). CLI: `--hits-tsv` accepts one or more TSVs (`nargs="+"`); `--email` resolves from `--email` → `$NCBI_EMAIL` → `None`. Names are cached within a run (the same accession is queried once), and the script is safe to re-run (idempotent — it rebuilds the column).

**Accession routing (protein vs nuccore).** `fetch_names(accessions, email)` splits accessions and routes each to the correct Entrez database so a protein id can never poison a nucleotide batch:
- **Protein accessions** — prefixes in `_PROTEIN_PREFIXES` = `NP_`, `YP_`, `WP_`, `XP_`, `AP_`, `QNP`, `ADV`, `AGT` (detected by `_is_protein_acc`) → queried against the **`protein`** DB. The organism is parsed from the esummary `Title`, which is `<product> [<organism>]`; `_org_from_title()` extracts the trailing `[…]` bracket. This is what gives protein-DB hits a real parent organism.
- **Nucleotide accessions** — everything else → queried against **`nuccore`**. The organism is the genome `Title`, cleaned by `clean_title()` which trims trailing `", complete genome"` / `", partial …"` / `", whole …"` / `", genome assembly"` / `", DNA"` qualifiers.

`is_ncbi(s)` (regex `^[A-Z]{1,2}_?\d{5,8}`) decides whether a `genome_id` is an NCBI accession at all; only those are sent for lookup. `_esummary_into()` batches 40 ids per esummary call and, **on any batch error, retries each id individually** so one bad accession cannot drop the whole batch's names; it stores names under both the versioned `AccessionVersion`/`Caption` and the version-stripped base accession. Calls are throttled (`time.sleep(0.34)` between requests) and `socket.setdefaulttimeout(60)` guards against a stalled NCBI connection freezing an unattended run.

**Accession fallback for failed / skipped lookups.** The per-row resolver `org(row)` is deliberate about what it will and will not call a genome:
- If a name resolved (versioned or base accession) → use it.
- Else if the `genome_id` **is** an NCBI accession → fall back to the **bare accession** (`gid.split(".")[0]`). A real cultured RefSeq record (e.g. `NC_…`) is a genuine genome; it is **never** mislabelled `"uncultured virus"` merely because the lookup failed or the run is offline.
- Only a genuinely non-NCBI id (metagenomic) falls back to `"uncultured virus (<db_name>)"` (db_name defaulting to `metagenomic`) — matching the GVD-AVrC / GPD uncultured-and-unclassified case.

**Never sends a placeholder email.** Per NCBI policy and the project rule, the script will **not** contact Entrez without a real address. If there are NCBI accessions to resolve but `email` is empty, it **skips the lookup entirely**, prints `no --email/$NCBI_EMAIL — skipping NCBI organism lookup (… accession(s)); using generic labels.`, and `names` stays `{}` so every row takes the accession/metagenomic fallback above. No address is ever fabricated. The 0-row / no-`genome_id` case still produces a well-formed table: it just ensures an empty `organism` column exists and writes back.

### `canonical.py` — `canonical_organism()` (honest organism counts)

**Purpose.** The single source of truth for collapsing the **same phage that appears under different database accessions** (cross-database / cross-accession redundancy — e.g. a phage catalogued in INPHARED *and* RefSeq under different ids). Imported wherever a count of "unique organisms" or a per-organism dedup must be consistent across every table, figure, and the tree.

**Signature.** `canonical_organism(name, fallback="") -> str` returns a lower-cased canonical key:
- Strips assembly/status prefixes `UNVERIFIED:`, `MAG:`, `TPA:` / `TPA_asm:` from the name.
- **Host-genus aliasing**: extracts the token after `phage`/`virus`, so `"Enterobacteria phage N4"` and `"Escherichia phage N4"` both collapse to `n4` — the same physical phage counts once regardless of which host genus a database used.
- **Metagenomic / unnamed entries are NOT collapsed**: if the name is blank or matches `uncultured|unclassified|metagenom|environmental`, it falls back to the **genome accession** (then lower-cased) so each distinct uncultured genome still counts exactly once rather than all merging into one bucket.
- Returns `""` only when neither a usable name nor a fallback accession is available.

**Consumers (honest counts).** `export_csv.py` imports it (`from canonical import canonical_organism`) to compute `n_organisms` for both `paper_main_table.csv` and the per-cluster table by deduplicating `{canonical_organism(organism, genome_id)}` over a group's rows — distinct from the physical record/genome counts (`database_records`, `n_genomes`/`n_databases`) which are NOT collapsed. (Reminder of the column semantics this feeds: `paper_main_table.csv`'s `database_records` is DB records for the **exact** sequence, not gene copies/paralogs, while `genome_metadata.csv` collapses by base accession into `accessions`/`n_accessions`.) `build_tree_of_hits.py` uses it to label/dedup tree tips, and `synteny_figure._dedup_loci_by_organism()` uses it to draw each canonical organism **once** (keeping the most complete neighbourhood, preferring a RefSeq `NC_`/`NZ_` accession) rather than as two near-identical tracks.

### Offline behaviour — `--no-annotate` and the email gate

`--no-annotate` (defined in `hmm_finder.py`, `action="store_true"`) runs the pipeline **fully offline**, suppressing every NCBI-dependent step. When set, the iterative driver logs `(organism annotation skipped: --no-annotate)` and never invokes `annotate_organism.py`; hits keep their accession/metagenomic organism labels via the fallback path, and local six-frame evidence is entirely unaffected.

**The email gate forces offline automatically.** `hmm_finder.py` resolves the Entrez email without ever assuming one, with precedence `--email > $NCBI_EMAIL > interactive prompt (TTY only) > offline`. The resolution happens late — only once a real run is committed — so `--list-databases` and argument errors never trigger a prompt. If no email is found (all sources empty), it logs `running offline — organism & protein-DB-sequence lookups skipped.` and **sets `args.no_annotate = True`**, guaranteeing no fabricated address is ever sent and no `input()` hangs an unattended run. (`run_pipeline.py`, the zero-prompt autonomous launcher, injects `--no-annotate` when absent for exactly this reason, alongside `--all-databases`/`--skip-tool-check`.) The two NCBI-only capabilities that go dark offline are (1) organism-name lookups and (2) retrieval of protein-database hit **sequences** in `extract_validated_hits.py` (protein-DB hits are retained without their NCBI-fetched AA). Six-frame translation, the HMM search, validation, controls, and all local figures still run.

The functional gene annotation (`annotate_genes.py` / VOGDB) is a **separate, local** concern: it depends only on the cached VFAM library and `hmmscan`, not on NCBI or `--email`. If the VOGDB cache is absent (`is_ready()` False) the figures fall back to "hypothetical protein" labels and grey `hypothetical / unknown` colouring, but the run is otherwise unaffected.

---

## Downstream figures & trees

These six scripts (`scripts/`) consume a converged discovery run's tabular outputs (`hits.tsv`, the validated domain FASTAs, the run manifest) and turn them into the publication deliverables: linear genome maps, anchored synteny panels, real-sequence GenBank neighbourhoods, an ML tree of the discovered homologues, the interactive clinker figures, and the single-file HTML report. They are designed to **degrade gracefully** — any optional renderer, NCBI fetch, or browser that is unavailable is logged and skipped rather than crashing the run. All figure scripts emit raster (PNG, 300 dpi) plus editable vector (SVG with real `<text>`, PDF with embedded TrueType) so figures can be retouched in Inkscape/Illustrator.

A shared colour scheme runs through every figure: the **gene of interest is bold gold** (`FAMILY_COLOR = "#ffd400"`), neighbours are coloured by **broad functional category** (the eight categories from [annotation](#annotategenespy--functional-annotation-via-vogdb-vfam-hmmscan)), and unannotated genes are light grey (`HYPO_COLOR = "#dde2e8"`). `synteny_figure.py` is the single source of truth for that scheme; `genome_map.py` imports it so maps and synteny panels colour genes identically.

### `genome_map.py` — linear locus map marking the gene of interest

A linear genome map drawing the gene of interest among its neighbours, shared by the single-genome scan (`scan_genome.py`) and the database discovery run (`build_real_genbanks.py`). The public entry point is `draw(...)`; gene lists are assembled by `build_genes(...)`; a tool-agnostic GenBank of the locus is written by `write_locus_genbank(...)`.

**Renderers (`tool=`, `MAP_TOOLS`).** Default is **`'dfv'` = DNA Features Viewer** (Edinburgh Genome Foundry), the cleanest publication renderer. Selectable alternatives: `'pub'` (the built-in matplotlib diagram, always available — strand arrows on packed lanes plus a full-length coordinate ruler), `'pygenomeviz'`, and `'easyfig'`. `'auto'`/`'matplotlib'` are accepted aliases. Any renderer that is unavailable falls back along the chain `easyfig/pygenomeviz/dfv → pub`; `pub` is the final safety net. The renderer is selected per-run in `build_real_genbanks.py` from the env vars `GENOME_MAP_TOOL` (default `dfv`), `GENOME_MAP_PALETTE`, `GENOME_MAP_FUNCTIONAL`, `GENOME_MAP_BRACKETS`.

- **`build_genes(anchor, called, flank_keys=None, label_of=None)`** — builds the gene list in absolute genome coordinates. `anchor = (a_start, a_end, a_strand)` is the gene of interest; `called = [(s, e, strand, meta)]` are the neighbouring calls. A called gene that **reciprocally overlaps** the anchor (`>0.6` of both lengths) is dropped as the gene of interest's own call, but a **nested overprint partner** (a long gene that merely spans it) is kept so the overprinting relationship stays visible. Each gene is assigned a functional `category` via the synteny categorizer on its product, and a `role` of `anchor`/`overlap`/`flank`/`other`. Returns dicts `{start, end, strand, role, label, category}`. The name is stored for **every** gene (so the locus GenBank keeps gp/locus names), while the renderer decides which labels to *display* via its density budget.
- **`draw(genes, anchor, out_base, title, ...)`** — key parameters: `track_name` (phage label; may be two lines `name\naccession`), `labels` (toggle gene-name labels), `palette` (`'default'` or `'colorblind'`, the Paul Tol muted set), `functional_labels` (also tag the gene of interest + its overlap partner with their functional category — the `--functional-labels` feature), `module_brackets` (bracket contiguous same-category runs with the module name; dfv only — the `--module-brackets` feature), and `genbank` (the locus GenBank, required by Easyfig). Emits `<out_base>.png/.svg/.pdf`.
- **`_draw_dfv(...)`** — the default renderer. Notable engineering: overlapping genes are auto-stacked onto separate levels via a **patched row-packer** (`_dfv_tolerant_levels`, 60 bp tolerance) that collapses few-bp start/stop-codon overlaps onto the baseline (keeping dense phage genomes clean) while still giving a *genuine* overprint (overlap ≫ tolerance) its own row; the patch is carefully routed around DFV's label-box packing so it never forces nearby labels onto one row. **Multi-line wrapping** triggers for big genomes (`>40` genes, `multiline`), splitting onto multiple lines (~35 genes/line) so every gene has room; the label budget becomes unlimited when wrapped. A **density-aware label budget** (~5 labels/inch when single-line) always labels the gene of interest and its overlap partners, then labels other genes closest-to-anchor first. A post-draw fix raises the gold anchor arrow above the label boxes (`zorder=30`) so the subject gene is never occluded in the raster.
- **`_module_runs` / `_draw_module_brackets`** — find contiguous runs (≥2 genes) of one functional category (the phage's structural/replication/… "modules", skipping the anchor and hypotheticals) and draw a labelled bracket above each on whichever wrapped line it falls (genoPlotR/Phamerator modular-organisation convention).
- **`_draw_pub(...)`** — the always-available matplotlib renderer. Genes are strand arrows on greedily **packed lanes** (`_pack`, so overlapping genes never hide each other), labels are stacked into non-overlapping rows using the *actual rendered label widths* with leader lines, and a genome-coordinate ruler (`_nice_bar` tick spacing, kb labels) runs the whole length. Figure height grows with lane and label-row count.
- **`_legend_handles(...)`** — legend patches for the categories actually present, each annotated with a **count** (e.g. `structural (3)`), so the functional composition is readable without counting arrows.
- **`write_locus_genbank(genes, contig_seq, organism, accession, out_path)`** — the tool-agnostic deliverable: a GenBank of the window spanning `genes` with the real contig sequence, CDS features carrying the gp/annotation labels, and the gene of interest marked `/gene=gene_of_interest`. Opens directly in Easyfig, Artemis, clinker, pyGenomeViz, etc. Handles the 1-based-inclusive → 0-based-half-open coordinate conversion (the `+1` on the end is required or every CDS prints 1 bp short at its 3′ end).
- **`_draw_easyfig(...)`** — shells out to Easyfig (`github.com/mjsull/Easyfig`, a Python-2 standalone, *not* a pip/conda package); enabled by setting `$EASYFIG_PY` (and optionally `$EASYFIG_PYTHON` to a python2). Raises → caller falls back when unavailable. `_draw_pgv` (pyGenomeViz) and `_draw_mpl` (a relative-coordinate fallback) round out the renderer set.

### `synteny_figure.py` — anchored, function-coloured synteny panels (one per cluster)

Publication-quality static synteny panels built the way synteny is shown in papers: every locus is **anchored on the family gene** (aligned column) and **strand-normalised** so the family gene always points right; neighbourhood genes are functionally annotated (VOGDB VFAM via `hmmscan`) and coloured by function; homologous genes are joined by shaded links between adjacent loci. **Inputs:** a clinker directory (`--clinker-dir`, containing `genbank_files/` + `cluster_membership.tsv`). **Outputs:** per-cluster SVG+PNG+PDF (`cluster_<cid>_synteny.png`, plus `_conservation` variants), a neighbourhood heatmap, a per-gene annotation CSV (`neighbour_gene_annotations.csv`), a conservation CSV, and `index.html`. Needs matplotlib, Biopython, CD-HIT, and VOGDB (downloaded once into the shared cache; falls back to "hypothetical protein" labels if absent).

- **Functional categorisation (`categorize`, `CATEGORY_RULES`, `CATEGORY_COLORS`).** Genes are binned into broad categories by first-match keyword search over their VOGDB function description; `_VOG_CAT` provides an unambiguous single-letter fallback (see [the VOGDB mapping](#annotategenespy--functional-annotation-via-vogdb-vfam-hmmscan)). Rule ordering matters: transcription/regulation is tested **before** replication so "(DNA-directed) RNA polymerase" → transcription while "DNA polymerase" still falls through to replication, and the structural rule uses specific terms (not bare "virion", which over-greedily caught "virion DNA-directed RNA polymerase") — directly relevant to the gp75 flagship, an overprinted gene antisense within a virion RNA polymerase. Hypothetical/uncharacterized/DUF descriptions are never categorised.
- **Palette (`PALETTES`).** `'default'` (`CATEGORY_COLORS`) or `'colorblind'` (`CATEGORY_COLORS_CB`, the Paul Tol *muted* qualitative set, distinguishable under deuteranopia/protanopia/tritanopia); selected with `--palette`. The gene of interest stays gold and hypothetical stays grey in both. Genome maps and synteny figures share the scheme for consistency.
- **`--color-by` modes (`color_of`).** `function` (category colours), `conservation` (a `Blues` gradient from core→unique, where a gene's shade is its orthogroup's `_conservation` fraction = fraction of loci containing it), or `both` (draws both variants). The CLI default in `hmm_finder.py` is `both`. In **all** modes the homology links between adjacent rows are shaded by sequence similarity (`_similarity`, a cached `difflib` ratio): a **darker link = more similar** pair of homologous genes.
- **Orthogroups & consensus functions.** `assign_orthogroups` runs CD-HIT on all neighbourhood proteins (40% id, `-aL 0.6`) and tags each gene `g['og']`; `consensus_functions` then sets each gene's function to its orthogroup's most-common named function so a single mislabelled member doesn't fragment a colour.
- **Locus dedup (`_dedup_loci_by_organism`).** Collapses the *same phage* appearing under multiple database accessions (e.g. an INPHARED id and a RefSeq NC_ id) to one track, keyed by `canonical_organism` (see [`canonical.py`](#canonicalpy--canonical_organism-honest-organism-counts)), preferring the most complete neighbourhood then a RefSeq accession — so a phage is drawn once, not as two near-identical rows.
- **`MAX_LOCI` cap disclosure.** `MAX_LOCI = 12` (CLI `--max-loci`, `0` = all loci on one figure). When a cluster exceeds the cap it is subsampled by **even striding** to a representative set, and — critically — the full pre-cap size is recorded (`full_counts`) and **disclosed on the figure title** as `n=<shown> of <total> loci shown (subsampled)` via `n_full` in `draw_cluster`, plus a console note, rather than silently presenting the truncated count as the cluster size.
- **`neighbourhood_rows(cid, loci)` → `neighbour_gene_annotations.csv` (`NEIGHBOUR_COLS`).** The ordered gene-neighbourhood table that lets you describe or manually label the genes bordering the gene of interest. Loci must already be anchored. Columns: `cluster, genome_id, organism, pos_index` (gene **order** relative to the gene of interest: `0` = it, negative = upstream, positive = downstream), `is_anchor`, `rel_start, rel_end` (bp relative to the gene), `strand_vs_gene` (`+` when a neighbour runs the same direction as your gene), `length_bp`, `distance_to_anchor_bp`, `orthogroup`, `category`, `vfam`, `function` (the anchor reads `GENE OF INTEREST (family homologue)`). Catalogued in full [below](#05_synteny).
- **`build_heatmap`.** A genomes × functional-category-count heatmap (`neighbourhood_heatmap.{png,svg,pdf}`) plus `neighbourhood_conservation.csv` — shows whether the gene of interest consistently sits among the same kinds of genes (a conserved module) across genomes.
- **`set_row_labels`.** Row labels are the organism name, with the accession appended when the organism is missing, generic ("uncultured virus", etc.), or duplicated, so rows stay distinguishable.

### `build_real_genbanks.py` — real-sequence neighbourhood GenBanks (+ genome maps)

Rebuilds family gene-neighbourhood GenBank files containing **real nucleotide sequence** (not the `N`-placeholders that clinker uses), so they open directly in Artemis/Geneious/UGENE. **Inputs:** `--hits-tsv` (a run's `hits.tsv`), `--out-dir`, `--email`. **Outputs:** one `*.gbk` per genome (named `<PhageName>_<accession>.gbk` where a name is known) plus a `*_genome_map.{png,svg,pdf}` per genome.

Per-genome procedure: (1) obtain the sequence — NCBI accessions via Entrez `efetch`/`esummary` in batches of 40 (`fetch_ncbi`, which also returns parsed phage **names**), metagenomic contigs by **streaming** the source catalogue gzip (`fetch_catalogue`; GVD-AVrC for `GutCatV1_`, GPD for `uvig_` prefixes); (2) gene-call flanking context with **Prodigal** (`-p meta`, cached per genome); (3) cut a window of the hit ORF(s) ± `FLANKS = 7` nearest flanking genes (±100 bp pad); (4) write a `SeqRecord` whose central CDS is the genuine family ORF (`/gene="family"`, validated translation) and whose flanking CDS are Prodigal genes. A genome with >1 hit gets all its family ORFs in one record. It then calls `genome_map.draw(...)` (via `genome_map.build_genes`) to render the map marking the gene of interest among its neighbours.

Two safety behaviours matter: only **six-frame ORF hits** (which carry genomic coordinates) build neighbourhoods — protein-DB hits (e.g. RefSeq `YP_`) are skipped here as they have no genome (`source_type == "six_frame_orf"` filter). And the **NCBI email is never assumed**: if `--email`/`$NCBI_EMAIL` is unset, NCBI genomes are skipped outright (no placeholder address is ever sent); metagenomic-catalogue genomes are unaffected. `socket.setdefaulttimeout(60)` bounds NCBI calls so unattended runs can't hang.

### `build_tree_of_hits.py` — ML tree of the discovered homologues (SH-aLRT + UFBoot)

Builds a maximum-likelihood phylogeny of the **discovered** homologues (distinct from the seed-only tree the main pipeline makes). **Pipeline:** MAFFT align → trimAl (`-gt 0.5`) → IQ-TREE (ModelFinder `-m MFP`, `-B 1000` UFBoot **and** `-alrt 1000` SH-aLRT). **Input:** a FASTA of unique, ORF-validated family-domain proteins (`--faa`, e.g. `validated/hits_unique_aa.faa`). **Outputs:** `tree_input.faa`, `hits.aln.faa`, `hits.aln.trim.faa`, `hits.treefile`, `hits.iqtree`, `hits.aln.stats.json`, an `alignment_figure.*`, and rendered tree images `hits_tree.{svg,png,pdf}` (plus a seed-pruned `hits_tree_homologs_only.*`).

- **Tip labels (`_build_tree_input`).** Organism-first so every figure reads cleanly: discovered hits are `Organism_accession_xN`, where **N is the number of distinct organisms carrying that exact domain sequence**. The tree is deduplicated by sequence, so this `_xN` recovers the discovery breadth that exact-sequence dedup would otherwise hide; organisms are counted by `canonical_organism` (not raw accession) to avoid double-counting the same phage that appears in several databases under different accessions. Optional seeds (`--seeds`) are parsed from their headers (`_organism_from_desc` handles UniProt `OS=`, NCBI brackets, and genome-title forms) and marked `*_seed` so you can see where the starting sequences fall among the discovered homologues.
- **Alignment strategy (`--mafft-mode`, `_mafft_strategy`).** Default `accurate` uses **L-INS-i** (`--localpair --maxiterate 1000`, the gold standard for homologous domains in variable-length context) where tractable (≤500 seqs) and falls back to `--auto` for very large sets, via the engine's `accuracy_flags`; or force `linsi`/`ginsi`/`einsi`/`auto`/`fftns`. The alignment is a first-class deliverable: `alignment_quality` writes `hits.aln.stats.json` (sequence/column counts, conserved columns, mean pairwise identity, gap %) and `alignment_figure` exports a ClustalX-coloured MSA (PNG/SVG/PDF). Both come from `engine/pipeline/alignment.py` and degrade gracefully if the engine isn't importable.
- **IQ-TREE invocation.** `-m MFP -B 1000 -alrt 1000 -T AUTO -ntmax <cpu> -seed 12345 --prefix <...> -redo`. `-T AUTO` with `-ntmax` avoids IQ-TREE's "more threads than CPU cores" abort; `-seed 12345` fixes the stochastic ML search + UFBoot resampling so reruns on the same alignment yield an identical tree (golden-output regression). **Two support measures** are printed on each branch as `aLRT/UFBoot` because reviewers expect both — UFBoot alone can be over-optimistic.
- **Tree rendering (`_render_newick`).** Renders editable SVG + 300-dpi PNG + PDF on an **opaque white** canvas (toytree's default is transparent → black on dark viewers). Tips are coloured by **host genus** (first token of the organism label) with a legend. A node carries a **support dot** only when it is robustly supported — **SH-aLRT ≥ 80 AND UFBoot ≥ 95** (Minh et al. 2020), parsed from IQ-TREE's `"aLRT/UFBoot"` node label by `_supported` (with a single-value `≥ 80` fallback for older trees). Layout `'r'` rectangular or `'c'` circular. `_homologs_only_newick` writes a seed-pruned copy so the discovered homologues can be shown in a legible result figure separate from the seed-context tree (circular layout is left to iTOL/FigTree — toytree 3.0.10 has no circular layout).

### `cluster_and_clinker_corrected.py` — corrected clusters + clinker figures

Clusters the **correct** family domains and builds clinker gene-neighbourhood figures grouped by those clusters, with the central gene being the genuine family six-frame ORF (the original output mistakenly used the overlapping Prodigal CDS). **Inputs:** `--validated-dir` (`hits.tsv` + `hits_unique_aa.faa`), `--cache-dir` (a `synteny_context_cache/` of `<genome>.fna` files), `--out-dir`, optional `--email`. **Outputs:** `genbank_files/*.gbk`, `cluster_membership.tsv`, per-cluster interactive `clinker_figures/cluster_<cid>.html`, optional static `cluster_<cid>.png`, and `index.html`.

Steps: (1) **CD-HIT cluster** the validated domains (`cdhit`, 40% id, `-aL 0.8`); (2) map every hit (from `hits.tsv`) to a cluster by exact sequence match to a unique representative; (3) for each kept hit, **`build_genbank`** cuts a window of ±`FLANKS = 7` nearest Prodigal flanking genes around the real family ORF — the central CDS is named identically (`/gene="family_homologue"`, `/product="FAMILY HOMOLOGUE (HMM hit)"`) in **every** track so clinker colours and links it consistently and it is easy to spot; tracks are file-named by organism so they read as phage names; (4) group GenBanks by cluster and run **clinker** per cluster with ≥2 loci (`clinker ... -i 0.3 -j 4`).

- **Locus dedup (`dedup_synteny_loci`).** Collapses hits mapping to the **same genomic locus** (same parent accession, overlapping coordinates) — e.g. a six-frame ORF and the RefSeq protein of the *same* gene resolved via `coded_by` — preferring the six-frame ORF then higher bit score. Genuine paralogs (non-overlapping loci on one genome) are kept separate; hits without coordinates are left untouched. Reports `N coordinate hits → M distinct loci`.
- **Protein-DB neighbourhoods (`resolve_protein_neighborhoods`, `_parse_coded_by`).** When `--email` is supplied, `annotated_protein` hits are given a neighbourhood too: the GenPept record's `coded_by` qualifier is parsed to the parent genome accession + CDS coords, the genome is fetched into the cache, and coordinates are filled in so `build_genbank` can draw it. Degrades gracefully on any NCBI/parse failure.
- **`MAX_LOCI` cap.** `MAX_LOCI = 16` per clinker figure (a figure with 30–80 tracks is unreadable); larger clusters are **evenly sampled** to 16 loci. `index.html` discloses **loci shown / total** per cluster, plus interactive and static-PNG links.
- **Static PNG (`_html_to_png`).** clinker's plot is JS-rendered (clustermap.js), so a static export needs a headless browser; `_html_to_png` uses Playwright/chromium when installed, with a one-time note when it isn't (the `synteny_figure.py` panels remain the primary static figures regardless). `socket.setdefaulttimeout(60)` bounds the `coded_by`/genome fetches.

### `generate_report.py` — single-file HTML discovery report

Writes a self-contained, portable `report.html` at the root of a `<name>_discovery` directory (mirrored into `PACKAGE/` if present). Reads `run_manifest.json`, `hit_summary.csv`, and `paper_main_table.csv`, and base64-embeds the tree, alignment, and a representative synteny PNG so the single file is portable. **Never raises on a missing input** — every section is conditional on its data existing. Entry point `generate(discovery)`; CLI `--discovery-dir`.

**Report sections (in order):**
- **Headline cards** — Family (label), Iterations, Databases, Hits (final run), Unique homologs, Organisms (from the last `hit_summary.csv` row and manifest parameters).
- **Calibration & convergence** — the iteration stopping criterion, plus sensitivity (seeds), specificity (controls), false-positive rate, and ROC AUC / optimal Youden bit-score from `threshold_calibration`. The prose makes clear the ROC is **advisory** (the fixed strict threshold is kept for tiering) and that the controls count seed-recovery vs negative-control sequences scored above threshold; it points to `controls/roc_curve.svg` and `controls/control_report.json`.
- **Files** — buttons (rendered only if the target exists) to the main table, all-hits supplementary CSV, per-run summary, database provenance, `METHODS.md`, the interrupted/overprinted TSV + the three interrupted FASTAs (domain AA, full-ORF AA, full-ORF nt), the alignment figure/FASTA, the per-hit HMM alignment (Stockholm + A2M), and both the publication synteny `index.html` and the interactive clinker `index.html`.
- **Per-run summary** table; **MSA** (embedded `alignment_figure.png` + stats from `hits.aln.stats.json`); an **inline residue-coloured MSA** of up to 30 hits × 180 columns (Clustal-ish residue classes → CSS colours via `_AA_CLASS`); the **phylogeny** (embedded `hits_tree.png`); a **representative synteny panel** (first `cluster_*_synteny*.png`).
- **Top homologs** — up to 25 rows from `paper_main_table.csv`, showing `rank, representative_organism, accession, n_genomes, n_organisms, database_records, domain_aa_len, best_evalue, best_bit_score, confidence_tier`. Per the table's semantics, `domain_aa_len` is the HMM-matched conserved region (the reported homolog length, ~137 for gp75), `database_records` is DB records for the *exact* sequence (not gene copies/paralogs), and `n_genomes`/`n_organisms` are the converged-run breadth.
- **Stop-interrupted / overprinted homologs** (present only with `--find-interrupted`) — from `interrupted_homologs.tsv` + the manifest. The prose explains read-through translation (stop codons retained, not broken on) scanned with the family HMM, and renders the **overprinting test**: counts of `strong`/`partial` support and a careful statement that the discriminating, length-dependent signal is the **antisense reading frame being OPEN across the whole domain** (`antisense_open_stops = 0`, improbable by chance for a long domain), whereas a single stop being synonymous in the antisense frame is, alone, expected ~85–100% of the time from the genetic code — necessary but weak. It states explicitly this is a **necessary sequence signature of overprinting, NOT proof the antisense ORF is expressed.** Per-row columns surfaced: `contig, strand, frame, domain_nt_start, domain_nt_end, domain_aa_len, internal_stops, stop_nt_positions, overprinting_support, domain_bit_score, i_evalue` (plus `overprinting_support, antisense_open_frame, antisense_open_stops, stop_silent_antisense` referenced in the TSV).
- **Tool versions** — from `manifest.tool_versions`.

---

## Outputs catalog — every file and column

A run writes everything under `<out-dir>/` (default `<fasta>_discovery/`). That directory has three layers:

- **Working data** — per-iteration `run1/`, `run2/`, … (each with `benchmark/results/`, `benchmark/validated/`, `benchmark/hmm/`) and a `downstream/` tree built from the canonical run. This is the raw material the assembler copies from; you rarely open it directly.
- **Root-level merged tables** — the CSV/figure files `export_csv.py` and `generate_report.py` write directly at `<out-dir>/` (e.g. `all_runs_hits.csv`, `paper_main_table.csv`, `report.html`, `run_manifest.json`, `METHODS.md`, and, with `--find-interrupted`, `interrupted_homologs*.{tsv,faa,fna}`).
- **`PACKAGE/`** — the clean, self-contained, shareable copy. `assemble_package()` in `scripts/hmm_finder.py` mirrors the important files into eight numbered folders; `package_layout.write_readmes()` then drops a plain-text `README.txt` into the package root and every subfolder describing each file by exact name, glob, or extension.

The single source of truth for the `PACKAGE/` folder names is `scripts/package_layout.py` (`DIRS` dict); both the assembler and the table exporter import it so the structure never drifts.

### The `PACKAGE/` layout (reading order)

| Folder | Purpose | Key contents |
|--------|---------|--------------|
| `README.txt` | Guide to the whole package | written by `package_layout._top_readme` |
| `METHODS.md` | Human-readable methods + citations for this run | copied from `<out-dir>/METHODS.md` |
| `run_manifest.json` | Machine-readable provenance (params, tool versions, calibration, seed recovery) | copied from `<out-dir>/run_manifest.json` |
| `01_summary_tables/` | Headline result tables + per-database hit bar chart — **start here** | `paper_main_table.csv`, `hits_deduplicated.csv`, `hit_summary.csv`, `database_hit_summary.csv` (+ `database_hits.png/svg/pdf`), `genome_metadata.csv`, `homolog_stats.csv`, `all_runs_hits.csv`, `database_summary.csv`, and (with `--find-interrupted`) `interrupted_homologs.tsv` |
| `02_sequences/` | All discovered sequences | `all_hits_aa.faa`, `all_hits_nt.fna`, `unique_homologs_aa.faa`; interrupted FASTAs; `per_run/runN/` |
| `03_hmm_profile/` | The calibrated profile HMM (the model the whole run is built on) | `profile.hmm` |
| `04_alignment_phylogeny/` | MSA of homologs + per-hit HMM alignment + the ML tree | `hits.aln*`, `hits_hmmalign.sto/.a2m`, `hits.treefile`, `hits_tree.*`, etc. |
| `05_synteny/` | Gene-neighbourhood comparisons + real-sequence GenBanks | `clinker/`, `publication_figures/`, `genbank_with_sequence/` |
| `06_database_summaries/` | The engine's per-database summary, one file per iteration | `run{N}_summary.tsv` |
| `07_seed_qc/` | Quality control of the **input** seeds | `seed_recovery.csv` + seed-only alignment/QC tree |
| `08_scripts/` | A verbatim copy of the scripts that produced this run | from `scripts/` (excludes `__pycache__`/`*.pyc`) |

`controls/` (the threshold-calibration outputs) and `report.html` live at the **`<out-dir>` root**, not inside `PACKAGE/` (though `report.html` is also copied into `PACKAGE/`). `assemble_package()` publishes the **most-complete run's** HMM as `profile.hmm` (`best_i`), so the model, the figures, and `paper_main_table.csv` all describe the same converged hit set.

### `01_summary_tables/`

#### `all_runs_hits.csv` — the complete, un-collapsed hit table

The whole record: every hit from every iteration, concatenated by `export_csv.export()` from each `run*/benchmark/validated/hits.tsv` with their full column set. This is the same schema as the per-run `hits.tsv` below. It is the supplementary "everything" file; the collapsed views (`hits_deduplicated.csv`, `paper_main_table.csv`) are derived from it.

#### `hits.tsv` / per-run `hits.csv` — the per-hit evidence table (36 base columns)

Written by `extract_validated_hits.py` (the authoritative ORF-validated extractor; see [its description](#extract_validated_hitspy--orf-validated-extraction-from-a-completed-run)) to `run{N}/benchmark/validated/hits.tsv`. The column order is fixed by the `COLS` list in that script. After extraction, `annotate_organism.py` inserts an **`organism`** column immediately after `genome_id` (unless `--no-annotate`), so an annotated table has **37** columns. Each row is one validated hit — either a six-frame genome ORF or an annotated protein-DB target.

**Identity & provenance**
| column | meaning |
|--------|---------|
| `hit_id` | unique hit id — the six-frame ORF id (`protein_id`/`target_name`) for genome hits, or the protein accession for protein-DB hits |
| `genome_id` | source genome/contig id (genome hits) or the cleaned protein accession (protein hits) |
| `organism` | phage/organism name from NCBI (added by `annotate_organism.py`); `"uncultured virus (<db>)"` when unresolved, accession label when an NCBI lookup is skipped/offline |
| `contig` | source contig id |
| `db_name` | database the hit came from (e.g. INPHARED genomes, SwissProt, VOGDB VFAM) |
| `db_type` | `nucleotide` or `protein` |
| `run_label` | which iteration produced the row (stamped `--run-label`) |
| `source_type` | `six_frame_orf` (genome hit, ORF reconstructed from coordinates) or `annotated_protein` (protein-DB hit, sequence fetched from NCBI) |
| `source_url`, `source_sha256`, `accessed_at` | download provenance carried from the engine's hit table |

**Genomic location (six-frame hits; blank for protein-DB hits)**
| column | meaning |
|--------|---------|
| `nt_start`, `nt_end` | 1-based forward-strand bounds of the stop-bounded search ORF (used for synteny/maps) |
| `strand` | `+`/`-` |
| `frame` | reading frame tag parsed from the engine's `frame=` annotation |
| `orf_nt_start`, `orf_nt_end` | bounds of the **gene** read from its Met start; equal to `nt_start/nt_end` unless a Met anchor shifts the start (`+` strand shifts `orf_nt_start`, `-` strand shifts `orf_nt_end`) |

**ORF validation — the "is this a real gene?" evidence**
| column | meaning |
|--------|---------|
| `orf_aa_len` | length (aa) of the **surrounding ORF read from its Met start codon** (`met_anchor`, within `MET_MARGIN = 60` aa upstream of the domain envelope). For an overprinted domain riding a long stop-free antisense frame, no near Met exists, so the start anchors at the domain itself and this can be long with `has_start_M=False` — that is correct, not a defect |
| `domain_aa_len` | length (aa) of the **HMM-matched conserved domain** (the `env_from..env_to` envelope). **This is the homolog length to cite** (~137 for gp75) — the reliable conserved measure, independent of how the surrounding ORF is bounded |
| `domain_coverage` | `domain_aa_len / orf_aa_len`, rounded to 3 dp; reported but **not** exclusionary (a small domain in a long ORF is still real) |
| `has_start_M` | does the Met-anchored gene begin with `M`? `False` is **expected and correct** for overprinted domains on a stop-free antisense frame where no near-domain Met exists |
| `ends_at_stop` | `True` for six-frame ORFs (stop-bounded by construction); `NA` for annotated proteins |
| `internal_stops` | premature stops inside the domain; **must be 0** to pass — the core six-frame quality gate |
| `prodigal_concordant` | `True` if best same-strand Prodigal overlap ≥ 50% (informational; `NA` for proteins) |
| `prodigal_same_strand_pct` | best % overlap with a same-strand Prodigal gene (typically ~0 for antisense/alternate-frame homologs — that discordance is the novelty signal) |
| `in_coding_locus` | `True` if best any-strand Prodigal overlap ≥ 50% (the hit sits in a real coding-dense locus); informational unless `--prodigal-gate` |
| `prodigal_any_strand_pct` | best % overlap with any Prodigal gene |
| `passes_orf_filter` | keep/flag decision. Six-frame: `internal_stops == 0` (AND `in_coding_locus` **only if** `--prodigal-gate`, which is **off by default**). Protein-DB: `internal_stops == 0` |

**HMM statistics**
| column | meaning |
|--------|---------|
| `evalue`, `bit_score`, `bias_score` | hmmsearch domain statistics carried from the engine |
| `env_from`, `env_to` | 1-based domain-envelope bounds **on the ORF** (from `hmmsearch --domtblout`); define the sliced domain |
| `confidence_tier` | strict/moderate/weak tier from the engine's `hits_classified.tsv` (joined by target name) |
| `qc_flags` | engine QC flags, likewise joined from `hits_classified.tsv` |

**Sequences**
| column | meaning |
|--------|---------|
| `aa_sequence` | amino-acid sequence of the matched **domain** (the identity key used by all dedup/collapse logic) |
| `nt_sequence` | the matching DNA for the domain; **blank** for protein-DB hits |

`extract_validated_hits.py` also writes, beside `hits.tsv`: `hits_aa.faa`/`hits_nt.fna` (domain sequences), `orfs_aa.faa`/`orfs_nt.fna` (full surrounding ORFs), `hits_unique_aa.faa` (deduplicated passing domains that seed the next iteration), and — if any protein-DB hit could not be fetched from NCBI — `protein_hits_unfetched.tsv` (so the loss is visible and recoverable rather than silent). A zero-hit run still writes an empty `hits.tsv` with the full header so the pipeline continues cleanly.

#### `paper_main_table.csv` — the main result (one row per unique homolog)

Built by `export_csv._paper_table()` from the **canonical converged run** (the iteration recovering the most unique homologs, ties broken toward the later/converged round). Hits are collapsed by exact `aa_sequence`; the representative is the highest-bit-score copy. Sorted by `best_bit_score` descending, with a `rank` prepended.

| column | meaning |
|--------|---------|
| `rank` | 1-based rank by best bit score |
| `representative_organism` | organism of the best-scoring copy |
| `accession` | representative `genome_id` |
| `database` | representative `db_name` |
| `database_records` | number of **database records carrying this exact sequence** — the same gene catalogued under several accessions/DBs. **NOT** biological gene copies/paralogs |
| `n_genomes` | number of **physical** genomes carrying it, counted by **base accession** (`_base_acc` strips the `.<version>` suffix so `NC_023589.1` and `NC_023589` collapse to one) |
| `n_organisms` | number of unique organisms by canonical identity (host-genus aliases collapsed; metagenomic/unnamed fall back to genome accession) |
| `domain_aa_len` | the cited homolog length (HMM-matched domain) of the representative |
| `domain_coverage` | representative domain coverage |
| `best_evalue` | smallest E-value across the collapsed group (`%.2g`) |
| `best_bit_score` | largest bit score across the group (1 dp) |
| `confidence_tier` | representative confidence tier |
| `example_hit_id` | a representative `hit_id` to locate the row in `all_runs_hits.csv` |

#### `hits_deduplicated.csv` — unique homologs with cross-database provenance

Built by `export_csv._dedup_hits()` across **all** iterations, collapsed by exact `aa_sequence`, sorted by `n_organisms` (breadth — the discovery story), with an `H0001`-style `homolog_id` prepended.

| column | meaning |
|--------|---------|
| `homolog_id` | `H####` stable id |
| `representative_organism`, `representative_genome`, `representative_db` | from the best-scoring copy |
| `source_type` | `six_frame_orf` / `annotated_protein` |
| `n_organisms`, `organisms` | count and `;`-joined list of unique canonical organisms (the headline breadth metric, immune to one phage appearing under several accessions) |
| `n_databases`, `databases` | how many — and which — databases recovered it (multi-DB recovery = stronger evidence) |
| `n_genomes` | unique **base-accession** genomes |
| `n_runs`, `runs` | how many — and which — iterations recovered it |
| `n_copies` | raw number of collapsed rows |
| `domain_aa_len`, `best_evalue`, `best_bit_score`, `confidence_tier` | representative/best statistics |
| `aa_sequence` | the unique domain sequence (collapse key) |

#### `hit_summary.csv` — per-iteration totals

One row per `run_label` (`export_csv.export()`): `run`, `total_hits`, `passed_filter` (count of `passes_orf_filter == True`), `six_frame_hits`, `protein_db_hits`, `unique_sequences`, `unique_organisms`, `databases` (`;`-joined). Comparing rows across rounds is how you read convergence.

#### `database_hit_summary.csv` — every database searched (incl. 0-hit)

Built by `export_csv._db_hit_summary()` joining the engine's per-database record for the canonical run (`all_database_summary.tsv`) with validated-hit unique counts. One row per database **including those with zero hits** (their absence is informative — e.g. a gene missing from Pfam/Swiss-Prot), sorted by `hits` descending, plus a trailing `ALL (deduplicated across databases)` row.

| column | meaning |
|--------|---------|
| `database` | database name |
| `type` | `nucleotide (six-frame)` or `protein` |
| `status` | engine status for that database in the canonical run |
| `hits` | engine `hit_count` |
| `strict_hits` | engine `strict_count` |
| `unique_sequences` | unique `aa_sequence` among validated hits from that DB |
| `unique_organisms` | unique canonical organisms from that DB |
| `runtime_seconds` | engine runtime for that DB |

**`database_hits.png` / `.svg` / `.pdf`** — horizontal bar chart of `hits` per database (`_db_barplot()`), blue bars for hit-bearing DBs, **grey** for searched-but-0-hit DBs. 300-dpi PNG plus editable SVG and print-ready PDF.

#### `genome_metadata.csv` — Supplementary S1 (one row per physical genome)

`export_csv.export()` groups all hits by **base accession** so the same genome under versioned + unversioned ids (RefSeq `NC_023589.1` + INPHARED `NC_023589`) is **one** row, not two (otherwise the genome count inflates ~1.3× on gp75 via cross-database aliases).

| column | meaning |
|--------|---------|
| `genome_id` | the base accession (collapse key) |
| `accessions` | `;`-joined list of every raw accession alias seen for this genome |
| `n_accessions` | how many distinct accession aliases collapsed into this row |
| `organism` | first non-empty organism name |
| `host` | host genus parsed from `"<Genus> phage/virus …"` (`_host_from_organism`) |
| `databases` | `;`-joined databases that carried this genome |
| `source_type` | `;`-joined source types |
| `n_hits` | number of hits on this genome |

#### `homolog_stats.csv` — Supplementary S3 (per-hit homology statistics)

A projection of `all_runs_hits` to: `hit_id`, `organism`, `genome_id`, `db_name`, `source_type`, `run_label`, `evalue`, `bit_score`, `domain_aa_len`, `domain_coverage`, `confidence_tier`.

#### `database_summary.csv` — raw engine per-database summary

All iterations' `run*/benchmark/results/all_database_summary.tsv` concatenated with a leading `run` column — the engine's untransformed per-database record (database, status, hit/strict counts, `nt_orf_mode`, runtime, provenance).

### `02_sequences/`

Combined multi-FASTAs written by `export_csv._write_multifastas()`:

- **`all_hits_aa.faa`** — every validated hit, amino acid (all iterations). Header: `{hit_id} {organism} [{db_name}]`.
- **`all_hits_nt.fna`** — every validated hit, nucleotide (genome hits only; protein-DB hits have no DNA).
- **`unique_homologs_aa.faa`** — one sequence per unique homolog, rich header `{homolog_id} {representative_organism} n_organisms=… n_databases=…`.

All three strip any `*` so the files load cleanly in aligners/BLAST.

**`per_run/runN/`** — one folder per iteration (copied by `assemble_package`): `hits.tsv` (the evidence table above), `hits.csv` (comma-separated copy from `export_csv`), `hits.gff3`, `hits_aa.faa`/`hits_nt.fna` (domains), `orfs_aa.faa`/`orfs_nt.fna` (full surrounding ORFs), `hits_unique_aa.faa` (deduplicated domains that seeded the next iteration).

**`hits.gff3`** — written by `hmm_finder.write_gff3()`: one `CDS` feature per validated hit, `<genome_id>\tHMM-Discovery\tCDS\t<nt_start>\t<nt_end>\t<bit_score>\t<strand>\t0\t<attrs>`, attributes `ID`, `Name=family_homolog`, `organism`, `db`, `evalue`, `bit_score`, `domain_coverage`, `in_coding_locus`. Load in IGV/JBrowse/Artemis alongside the genome FASTA.

### `03_hmm_profile/profile.hmm`

The calibrated profile HMM from the **most-complete run** (`benchmark_profile.hmm`) — the model the figures and paper table describe. Search other databases with HMMER; submit to Pfam / NCBI CDD / VOGDB.

### `04_alignment_phylogeny/`

Copied from `downstream/tree/`: `hits.aln.faa` (MAFFT MSA of the homologs), `hits.aln.trim.faa` (trimAl-trimmed, the tree input), `hits.aln.stats.json` (length, gap %, conserved columns, mean pairwise identity), `alignment_figure.{png,svg,pdf}` (ClustalX-coloured), the per-hit HMM alignment `hits_hmmalign.sto` (Stockholm) and `hits_hmmalign.a2m` (A2M; UPPERCASE = HMM match columns, lowercase/`.` = insertions), the ML tree `hits.treefile` / bootstrap consensus `hits.contree` (Newick; **seeds included and marked `SEED_*`**), rendered `hits_tree.{png,svg,pdf}`, a seed-pruned `hits_tree_homologs_only.*`, the IQ-TREE report `hits.iqtree` and `hits.log`, and `tree_input.faa` (the exact organism-labelled sequences fed to the aligner). Built by [`build_tree_of_hits.py`](#build_tree_of_hitspy--ml-tree-of-the-discovered-homologues-sh-alrt--ufboot).

### `05_synteny/`

- **`clinker/`** — clinker gene-neighbourhood comparison (from [`cluster_and_clinker_corrected.py`](#cluster_and_clinker_correctedpy--corrected-clusters--clinker-figures)): interactive `cluster_*.html` (open in a browser; "Save SVG" for figures) and, when a headless browser is installed, a static `cluster_*.png` per cluster (the "Static PNG" column in `clinker/index.html`).
- **`publication_figures/`** — publication synteny panels (PNG/SVG/PDF) **plus `neighbour_gene_annotations.csv`**, the ordered gene-neighbourhood table written by `synteny_figure.py` (`NEIGHBOUR_COLS`), one row per neighbouring gene per locus, **anchored on your gene of interest**. Sort by `genome_id` then `pos_index` to walk each neighbourhood:

  | column | meaning |
  |--------|---------|
  | `cluster`, `genome_id`, `organism` | which synteny cluster / source locus the gene belongs to |
  | `pos_index` | gene order relative to your gene: `0` = your gene, `-1/-2…` upstream, `+1/+2…` downstream |
  | `is_anchor` | `1` for your gene of interest (the family homolog), else `0` |
  | `rel_start`, `rel_end` | gene coordinates relative to your gene (bp; your gene sits at 0), strand-normalised |
  | `strand_vs_gene` | `+` = same orientation as your gene, `-` = opposite |
  | `length_bp`, `distance_to_anchor_bp` | gene length; signed gap to your gene (− upstream, + downstream) |
  | `orthogroup`, `category`, `vfam`, `function` | cross-locus orthogroup, functional category, VOGDB VFAM, product/function |

- **`genbank_with_sequence/`** — real-sequence GenBank neighbourhoods named by phage (`<PhageName>_<accession>.gbk`, from [`build_real_genbanks.py`](#build_real_genbankspy--real-sequence-neighbourhood-genbanks--genome-maps)), each loadable in Artemis/Geneious/UGENE/clinker/pyGenomeViz, plus a `<name>_genome_map.png/.svg` marking the gene of interest (the HMM hit, gold) among its neighbours.

### `06_database_summaries/run{N}_summary.tsv`

One file per iteration (`assemble_package` copies each run's `all_database_summary.tsv`). Per-database results for that iteration: database, status, hit/strict counts, `nt_orf_mode`, runtime, and source provenance.

### `07_seed_qc/`

Quality control of the **input** seeds, copied from `<out-dir>/seed_qc/`.

**`seed_recovery.csv`** — written by `seed_recovery.seed_recovery_report()`, one row per input seed (by name, file order):

| column | meaning |
|--------|---------|
| `seed_id` | seed identifier (first token of its FASTA header) |
| `before_bit` | best full-sequence bit score vs the **initial** model (run1 HMM) |
| `before_recovered` | `before_bit ≥ STRICT_BITS` (45.0) |
| `after_bit` | best bit score vs the **final** refined model |
| `after_recovered` | `after_bit ≥ 45.0` |
| `status` | `recovered` / `lost_after_refinement` / `gained_after_refinement` / `never_recovered` |

A seed not recovered by the final model is usually a divergent outlier (consider dropping it or treating it as a sub-family). Aggregate counts go to `run_manifest.json → seed_recovery_qc` and `METHODS.md`. The folder also holds a seed-only `hits.aln.faa`, `tree_input.faa`, `hits.treefile`, and `hits_tree.*` for sanity-checking the seed set.

### `08_scripts/`

A verbatim copy of the `scripts/` directory that produced this run (excluding `__pycache__`/`*.pyc`), so the exact code is reproducible alongside the data.

### `controls/` — threshold calibration (at the `<out-dir>` root)

Written by `hmm_finder.run_controls()` via `engine/pipeline/controls.py` (skipped with `--no-controls`; mechanics in [Step 5 of the search core](#step-5--threshold-calibration-and-controls-controlspy)). Controls measure **specificity against unrelated proteomes — NOT a six-frame false-discovery rate.** The **positive** control is a seed **self-recovery** check (`seeds_self_test`, 100% recovery expected); the **negatives** are composition-matched shuffled seeds plus optional unrelated-proteome sets (fungal, mammalian housekeeping, archaeal proteins, etc., from `BUILTIN_CONTROLS`).

**`control_report.json`** — `ControlReport.to_json()` (the `summary()` dict):

| field | meaning |
|-------|---------|
| `sensitivity` | fraction of positive (seed) sequences scoring ≥ strict threshold (45) |
| `specificity` | `(total_negatives − false_positives)/total_negatives` at the strict threshold |
| `false_positive_rate` | `1 − specificity` |
| `total_positives`, `true_positives` | seed sequences scored / recovered at strict |
| `total_negatives`, `false_positives` | negative-control sequences scored / scoring ≥ strict |
| `n_controls`, `n_positive_controls`, `n_negative_controls` | counts of control sets |
| `sensitivity_strict`, `sensitivity_moderate`, `specificity_strict` | tiered aliases |
| `roc` | nested ROC calibration (below); `{}` if either side is empty |
| `results` | per-control detail (name, role, n_seqs, n_hits, score distribution) |

The nested **`roc`** object (`ControlReport.roc()`): `auc` (exact Mann-Whitney; 1.0 = perfect separation), `optimal_threshold` (Youden's-J optimal bit cutoff), **`optimal_threshold_defined`** (`True` only when ≥1 negative scores above the noise floor — **`False` when no negative is detected**, in which case the "optimum" carries no information from the negative distribution and is advisory only), **`n_negative_detected`** (count of negatives above the floor), `youden_j`, `sensitivity_at_optimum`, `specificity_at_optimum`, `threshold_method` (notes when the optimum is undefined), `n_positive`, `n_negative`, `separable` (`True` when J ≥ 0.999), `strict_threshold` (45). The ROC is **advisory** — the fixed strict threshold (45) is always retained for tiering.

**`controls_summary.csv`** — `ControlReport.to_dataframe()`: one row per control (`name`, `role`, `desc`, `n_sequences`, `n_hits`, `n_hits_strict`, `hit_rate_pct`, `min/max/mean_score`, `pass`).

**`roc_curve.{png,svg,pdf}`** — the ROC figure (`plot_roc()`), marking the Youden-optimal point and the strict threshold. A summary of all of this is mirrored into `run_manifest.json → threshold_calibration` and `METHODS.md`.

### `report.html` — one-page visual summary

`generate_report.py` writes `report.html` at the `<out-dir>` root (and copies it into `PACKAGE/`). Self-contained HTML (base64-embedded figures); section breakdown in [`generate_report.py`](#generate_reportpy--single-file-html-discovery-report).

### `interrupted_homologs.tsv` and its FASTAs (`--find-interrupted`)

Written by `find_interrupted.py` (`_run()` → `ROW_COLS`) and copied to `PACKAGE/01_summary_tables/`. This is the read-through scan: the searched nucleotide databases are six-frame translated **with stop codons retained** and searched with the family HMM, so it recovers homologs the normal stop-to-stop search cannot — domains carrying a premature internal stop (overprinted/pseudogenized genes). One row per candidate, deduplicated by `(contig, strand, frame, stop positions)` keeping the best-scoring copy, sorted by `domain_bit_score`. (Mechanics in [`find_interrupted.py`](#find_interruptedpy--read-through-detection-of-stop-interrupted-homologs).)

| column | meaning |
|--------|---------|
| `contig` | source contig accession |
| `organism` | organism / phage name for the contig (joined offline from `genome_metadata.csv`; falls back to the accession if unknown) |
| `strand`, `frame` | the read-through reading frame (strand + 0/1/2) |
| `domain_nt_start`, `domain_nt_end` | forward-strand 1-based genome coordinates of the matched domain |
| `domain_aa_len` | domain length (aa), `env_to − env_from + 1` — **the homolog length to cite** |
| `internal_stops` | number of premature internal stops in the domain |
| `stop_nt_positions`, `stop_aa_positions` | per-stop forward genome coordinate(s) and aa position(s), `;`-separated |
| `overprinting_support` | **`strong`** = the chosen antisense frame is fully open (0 stops) over the domain **AND** every premature stop is synonymous in it; **`partial`** = some stops silent; **`none`** = no evidence |
| `antisense_open_frame` | the antisense frame (0/1/2) with the fewest stops over the domain (the candidate overprinted ORF); `-1` if no DNA |
| `antisense_open_stops` | that frame's stop count (`0` = a fully open overlapping antisense ORF) |
| `stop_silent_antisense` | per-stop `1/0` — the premature stop is synonymous in that antisense frame |
| `domain_bit_score`, `i_evalue` | HMM domain score / independent E-value |
| `orf_aa_len` | full read-through ORF length, **EXCLUDING the terminal stop** (so it matches `scan_genome.py` / `extract_validated_hits.py` for the same gene) |
| `aa_before_first_stop`, `aa_after_last_stop` | intact residues before the first / after the last premature stop |
| `orf_nt_start`, `orf_nt_end` | full read-through ORF genome bounds |
| `natural_stop_nt` | forward genome coordinate of the **actual** (natural/terminal) stop codon; `0` if the ORF ran to the contig/window edge |
| `domain_nt` | matched-domain DNA |
| `domain_aa_with_stops` | matched-domain protein, each internal stop shown as `*` |
| `full_orf_aa` | full read-through ORF protein — premature stops kept as `*`, terminal `*` = the natural gene end |
| `full_orf_nt` | full read-through ORF nucleotide (coding 5'→3', ending in the actual stop-codon triplet; translates back to `full_orf_aa`) |

**The science of `overprinting_support`** (`analyze_overprinting()`): the **discriminating** signal is the antisense frame being **OPEN across the whole domain** — for a 137-aa domain that is ~0.7% likely by chance, so an open antisense ORF is strong evidence of a real overlapping gene. A single premature stop being *synonymous* in that frame is, **alone**, expected ~85–100% of the time from genetic-code geometry, so it is weak by itself. `strong` is therefore a **necessary but not sufficient SEQUENCE signature** of antisense overprinting (open frame + synonymy) — **not proof** the antisense ORF is expressed.

The three sibling FASTAs (always written so report/package links never dangle, named off the TSV stem):
- **`interrupted_homologs_domain_aa.faa`** — each matched domain protein with every internal stop as `*` (`domain_aa_with_stops`).
- **`interrupted_homologs_full_orf_aa.faa`** — the full read-through ORF protein (`full_orf_aa`); premature stops kept as `*`, terminal `*` = natural gene end.
- **`interrupted_homologs_full_orf_nt.fna`** — the full read-through ORF nucleotide (`full_orf_nt`).

A summary (candidate count, threshold basis, strong/partial counts) is mirrored into `run_manifest.json → interrupted_homologs` and `METHODS.md`.

### Single-genome scan outputs (`scan_genome.py`)

A targeted scan of one genome/accession (not the database sweep) writes into its own output dir (`ROW_COLS` in `scan_genome.py`). It reuses `find_interrupted`'s read-through machinery, so it reports both **clean** and **interrupted** copies. See [`scan_genome.py`](#scriptsscan_genomepy--single-genome-scan).

**`scan_hits.tsv`** — one row per hit. Columns: `contig`, `strand`, `frame`, `nt_start`, `nt_end`, `domain_aa_len`, `internal_stops`, **`status`** (`clean` / `interrupted`), `domain_bit_score`, `i_evalue`, `orf_nt_start`, `orf_nt_end`, `orf_aa_len` (excludes the terminal stop), `has_start_M`, `ends_at_stop`, `overprinting_support`, `antisense_open_stops`, `stop_nt_positions`, `domain_aa`, `orf_aa`, `orf_nt`.

**`scan_hits_aa.faa` / `scan_hits_nt.fna`** — the full-ORF protein/nucleotide of each hit; header `{contig}|{strand}{frame}|{nt_start}-{nt_end}|{status}|bit=…`.

**`scan_report.txt`** — a plain-text verdict (`_finish()`): **`GENE PRESENT`** (≥1 clean hit), **`PRESENT but INTERRUPTED`** (only stop-interrupted/overprinted copies), or **`GENE NOT DETECTED`**; contigs scanned, reporting threshold, and a bulleted list of the top 10 hits (domain/ORF lengths, M/stop flags, and for interrupted hits the overprinting verdict).

**`scan_neighbourhood.csv`** — the ordered neighbour table (`write_neighbourhoods()`, `SCAN_NB_COLS`), one row per neighbouring gene around each scan hit, anchored on the gene of interest. Gene names come from the genome's **own** GenBank annotation when available, else from de novo **Prodigal** calling (+ optional VOGDB VFAM). Columns: `hit`, `contig`, `pos_index`, `relationship` (upstream/downstream/overlapping — overlapping captures the overprint partner / nested genes, e.g. gp75's antisense RNA polymerase), `is_anchor`, `gene`, `product`, `locus_tag`, `protein_id`, `rel_start`, `rel_end`, `strand_vs_gene`, `length_bp`, `distance_to_anchor_bp`, `annotation_source`, `category`, `vfam`.

**`scan_genome_map_<hit>.png/.svg`** (+ `_whole` variants and `.gb` GenBank locus files) — genome-map figures (default renderer **DNA Features Viewer**, `dfv`) marking the gene of interest among its neighbours, plus a whole-contig view. When `--accession` is used, the fetched record is also saved (`<accession>.gb` GenBank with the genome's own annotation, `<accession>.fna` FASTA).

### `run_manifest.json` and `METHODS.md` (provenance, at the `<out-dir>` root, also in `PACKAGE/`)

Both written by `hmm_finder.write_methods_log()` (never raises). **`run_manifest.json`** is the machine-readable provenance record:

| field | meaning |
|-------|---------|
| `tool`, `code_git_commit` | tool name and short git commit of the code |
| `started_at`, `finished_at`, `command_line` | run window and the shlex-quoted command line (re-runnable verbatim) |
| `n_input_seeds`, `input` | seed count; input FASTA path + its SHA256 |
| `conda_env`, `conda_prefix`, `python` | environment provenance |
| `parameters` | label, iterations, cpu, databases, `prodigal_gate`, `min_recovery` (0.70), `max_synteny_genomes` (200), email, db_cache, out_dir, input_type, trans_table, no_annotate |
| `annotation_database` | provenance of the local annotation DB |
| `per_iteration_unique_seeds` | unique validated seeds recovered per run |
| `iteration_stop_reason` | why iteration stopped (convergence / no hits / max iterations) |
| `threshold_calibration` | the controls `summary()` (incl. the `roc` block above) |
| `seed_recovery_qc` | the `seed_recovery` summary (`strict_bits`, `n_seeds`, `recovered_before/after`, `not_recovered_after`, …) |
| `interrupted_homologs` | the read-through scan summary (candidate count, strong/partial overprinting counts, FASTA paths) |
| `tool_versions` | resolved versions of HMMER, MAFFT, trimAl, CD-HIT, IQ-TREE, Prodigal, clinker, etc. |
| `database_provenance` | per-database source URLs, SHA256s, and access dates aggregated from each run's `reproducibility.json` |

**`METHODS.md`** renders the same information as human-readable Markdown with citations: command/env, parameters, per-iteration recovery, stopping criterion, the calibration paragraph (sensitivity at strict 45, specificity, FPR, and the advisory ROC/Youden optimum), the seed-recovery QC, the overprinting/read-through narrative (with the open-frame-vs-silent-stop caveat), tool versions, and per-database provenance.

### Other root-level files

`run_command.txt` — written by **`run_pipeline.py`** (the zero-prompt autonomous launcher) to `<out-dir>/run_command.txt`: the fully resolved command it executed (after PATH-injecting the interpreter bin and injecting `--all-databases`/`--no-annotate`/`--skip-tool-check` when absent). Per-iteration `run{N}/benchmark/reports/reproducibility.json` (engine-generated) is the source of the `database_provenance` aggregated into the manifest.

**Key files for this catalog:**
- `scripts/extract_validated_hits.py` — `hits.tsv` 36-column writer (`COLS`), domain/ORF FASTAs, `MET_MARGIN`/`met_anchor`
- `scripts/export_csv.py` — `paper_main_table.csv`, `hits_deduplicated.csv`, `genome_metadata.csv`, `database_hit_summary.csv`, multi-FASTAs
- `scripts/find_interrupted.py` — `interrupted_homologs.tsv` (`ROW_COLS`), overprinting analysis, the three interrupted FASTAs
- `scripts/scan_genome.py` — `scan_hits.tsv`, `scan_neighbourhood.csv`, `scan_report.txt`, `scan_genome_map_*`
- `engine/pipeline/controls.py` — `control_report.json`, `roc()` (`optimal_threshold_defined`, `n_negative_detected`), `roc_curve.*`
- `scripts/seed_recovery.py` — `seed_recovery.csv`
- `scripts/hmm_finder.py` — `run_manifest.json`/`METHODS.md`, `assemble_package()`, `write_gff3()`, organism annotation step
- `scripts/package_layout.py` — `PACKAGE/` layout (`DIRS`) and README writers
- `scripts/generate_report.py` — `report.html`
- `scripts/synteny_figure.py` — `neighbour_gene_annotations.csv` (`NEIGHBOUR_COLS`)
- `scripts/annotate_organism.py` — adds the `organism` column to `hits.tsv`

---

## Complete CLI flag reference

This section is the canonical, exhaustive catalog of every command-line flag accepted by the three executable entry points, taken directly from their `argparse` blocks: `scripts/hmm_finder.py` (family discovery), `scripts/scan_genome.py` (single-genome scan), and `scripts/run_pipeline.py` (the zero-prompt autonomous launcher). Defaults shown are the literal `default=` values in the code. "Non-interactive-safe" means the flag does not, by itself, trigger a terminal prompt or block on `input()`; the launcher (`run_pipeline.py`) additionally detaches `stdin` so no step can ever block.

### `hmm_finder.py` — family/homolog discovery

Defined in `main()` of `scripts/hmm_finder.py` (argparse block at lines 666–742). `RawDescriptionHelpFormatter` is used so the module docstring prints verbatim. There are no required flags; if `--fasta` is omitted an interactive prompt asks for the seed. All flags are parsed by a single `ArgumentParser` (no subcommands).

#### Input / seed

| Flag | Type / choices | Default | Effect |
|------|----------------|---------|--------|
| `--fasta PATH` | `Path` | `None` | Seed protein FASTA. If omitted, the program prompts on a TTY (so omitting it is NOT non-interactive-safe). The launcher requires it. |
| `--input-type {auto,protein,nucleotide}` | choice | `auto` | Declares the seed FASTA type; `auto` detects nucleotide vs protein. A nucleotide seed is translated before the HMM is built. |
| `--trans-table N` | `int` | `11` | Genetic code used to translate a **nucleotide** seed. Default 11 = bacterial/archaeal/phage; e.g. 4 = Mycoplasma/Spiroplasma. No effect on a protein seed. |
| `--name LABEL` | `str` | `None` | Label for the output folder; default is derived from the FASTA filename. |
| `--out-dir PATH` | `Path` | `None` | Output root. Default is `<fasta>_discovery`. The launcher writes `run_command.txt` here. |

#### Iteration / compute

| Flag | Type | Default | Effect |
|------|------|---------|--------|
| `--iterations N` | `int` | `3` | Number of jackhmmer-style search iterations. Each iteration re-builds the profile from validated hits and re-searches, pulling in progressively more distant homologs. |
| `--cpu N` | `str` (numeric) | `"8"` | Thread count handed to hmmsearch/MAFFT/cd-hit/IQ-TREE/MEME. After parsing it is **clamped down** to `os.cpu_count()` if it exceeds available cores (prevents IQ-TREE exit-2 "more threads than cores" aborts and tool oversubscription); machines with ≥ the requested cores keep the full value, so the default 8 is untouched on any ≥8-core host. |
| `--smoke` | flag | off | Fast self-test: 1 iteration against a single small database. Plumbing check, not a real search. Non-interactive-safe. |
| `--skip-tool-check` | flag | off | Skip the startup external-software check (hmmsearch, mafft, etc.). Non-interactive-safe; the launcher injects this automatically. |

#### Databases

| Flag | Type | Default | Effect |
|------|------|---------|--------|
| `--databases NAMES` | comma-list `str` | `None` | Explicit comma-separated database set (names must match the catalog; see `--list-databases`). If omitted: a TTY **prompts** you to choose; a non-interactive run falls back to the full default set. |
| `--all-databases` | flag | off | Take the full default database set with **no** interactive prompt (use in scripts that want everything but still have a TTY). Non-interactive-safe. The launcher injects this when no DB-selection flag is present. |
| `--pick-databases` | flag | off | Interactively choose which databases to search. NOT non-interactive-safe (prompts). |
| `--list-databases` | flag | off | Print the catalog (names, sizes, rough times) and exit. No FASTA required. Non-interactive-safe. |
| `--db-cache PATH` | `Path` | `~/.cache/hmm-homologue-finder` | Persistent, cross-run cache for downloaded databases — each DB downloads once, ever, so repeat runs are much faster. Shared across all runs. |

#### Interrupted / overprinting mode

| Flag | Type | Default | Effect |
|------|------|---------|--------|
| `--find-interrupted` | flag | off | ALSO scan the searched nucleotide databases with **read-through** translation (stop codons kept, not broken on) to recover homologs interrupted by a premature stop — e.g. overprinted genes where a nonsense mutation here is silent in an overlapping gene. The normal stop-to-stop search misses these. Writes `interrupted_homologs.tsv`. Opt-in; adds an extra scan of each nucleotide DB. This is the flagship gp75 mode. |
| `--prodigal-gate` | flag | **off** | Require six-frame hits to overlap a Prodigal-predicted gene to pass (stricter, higher specificity). **Off by default** — Prodigal overlap is informational, not exclusionary — so that genuine antisense/alternate-frame homologs (the whole point of the tool) are not discarded. |

#### Controls / calibration

| Flag | Type / choices | Default | Effect |
|------|----------------|---------|--------|
| `--no-controls` | flag | off | Skip the threshold-calibration controls: the **positive control** (seed self-recovery sensitivity check) and the **specificity** test against unrelated/shuffled negatives. Controls measure specificity vs unrelated proteomes (NOT a six-frame FDR); the shuffled-seed negative needs no network. |
| `--biology-mode {generic,phage,bacterial}` | choice | `phage` | Control panel to calibrate against. Default `phage`, matching the default phage/viral databases. |
| `--download-controls` | flag | off | One-time fetch of the UniProt unrelated-proteome negatives (fungi/mammalian/archaea) into the bundle, then continue. Without it only the always-available shuffled-seeds negative is used (so `controls.roc` may emit `optimal_threshold_defined=false` when no real negative is detected). |

#### Figures / annotation

| Flag | Type / choices | Default | Effect |
|------|----------------|---------|--------|
| `--no-seed-tree` | flag | off | Skip the pre-run seed-QC tree + seed alignment (a quick phylogeny/alignment of just the input seeds, for sanity-checking the seed set). The **final** homolog tree always includes the seeds, marked. |
| `--synteny-gene-labels` | flag | off | Label neighbour genes with their functional annotation in the synteny figures. Off by default (an interactive run asks); can clutter dense neighbourhoods, where colour + legend already convey function. |
| `--color-by {function,conservation,both}` | choice | `both` | How synteny neighbourhood genes are coloured: by functional category (VOGDB), by cross-locus conservation, or both. |
| `--no-annotate` | flag | off | Skip the NCBI organism-name lookup (for fully offline / fastest unattended runs); hit tables still build. The launcher injects this when no `--email` is given. Non-interactive-safe. |
| `--email ADDR` | `str` | `None` | NCBI Entrez email for organism/sequence lookups. If omitted, the `NCBI_EMAIL` environment variable is used; if neither is set the run proceeds fully offline (equivalent to `--no-annotate`). Never hardcoded — NCBI requires a real address. |

### `scan_genome.py` — single-genome scan ("does THIS genome carry my gene?")

Defined in `main()` of `scripts/scan_genome.py` (argparse block at lines 556–607). Two **mutually exclusive, required** groups gate the call: one for the query model, one for the target genome. `--map-tool` default is `dfv` and `--trans-table` default is `11`.

#### Query model (mutually exclusive, exactly one required)

| Flag | Type | Effect |
|------|------|--------|
| `--seeds PATH` | `Path` | Seed FASTA (protein or nucleotide CDS) from which the profile HMM is **built**. |
| `--hmm PATH` | `Path` | An existing profile HMM to scan with (skip building). |

#### Target genome (mutually exclusive, exactly one required)

| Flag | Type | Effect |
|------|------|--------|
| `--genome PATH` | `Path` | A local genome to scan: nucleotide FASTA (`.fna`/`.fa[.gz]`) OR an annotated GenBank (`.gb`/`.gbk`/`.gbff`, whose gene names are used to label neighbours). |
| `--accession ACC[,ACC...]` | `str` | NCBI nucleotide accession(s) to fetch & scan, comma-separated (e.g. `KX098390` or `NC_031062`). Requires `--email`. |

#### Scan parameters

| Flag | Type / choices | Default | Effect |
|------|----------------|---------|--------|
| `--email ADDR` | `str` | `None` | NCBI email, required with `--accession` (or set `$NCBI_EMAIL`). Unused for a local `--genome`. |
| `--out PATH` | `Path` | `genome_scan` | Output directory (created on start). Holds `scan_hits.tsv`, `scan_hits_aa.faa`, `scan_hits_nt.fna`, `scan_report.txt`, the locus GenBank, and the genome map. |
| `--min-bit X` | `float` | `25.0` | Minimum domain bit score to report a hit. |
| `--trans-table N` | `int` | `11` | Genetic code for translating a **nucleotide seed** (default 11). |
| `--cpu N` | `int` | `4` | Threads for hmmsearch / Prodigal. |

#### Interrupted / overprinting

| Flag | Type | Default | Effect |
|------|------|---------|--------|
| `--find-interrupted` | flag | off | Also report stop-interrupted / overprinted copies, running the overprinting test. Note that the `find_interrupted` `orf_aa_len` it reports **excludes** the terminal stop. |

#### Neighbourhood / synteny

| Flag | Type / choices | Default | Effect |
|------|----------------|---------|--------|
| `--no-neighbours` | flag (`dest=neighbours`, `store_false`) | neighbours **on** | Skip Prodigal gene-calling of the flanking genes. When on (default, via `set_defaults(neighbours=True)`), writes `scan_neighbourhood.csv` (the ordered neighbour table). |
| `--flanks N` | `int` | `7` | Number of flanking genes to report on **each** side of the gene of interest. Overlapping genes are always included regardless of this count. |
| `--db-cache PATH` | `Path` | `~/.cache/hmm-homologue-finder` | Cache holding the VOGDB VFAM DB used to annotate neighbour genes (optional). |

#### Genome-map rendering

| Flag | Type / choices | Default | Effect |
|------|----------------|---------|--------|
| `--map-tool {dfv,pub,pygenomeviz,easyfig}` | choice | **`dfv`** | Genome-map renderer. `dfv` = DNA Features Viewer (clean strand arrows, overlapping genes auto-stacked onto their own level, auto label de-overlap, real coordinate axis; emits PNG+SVG+PDF). `pub` = the built-in matplotlib diagram (always available). `pygenomeviz` and `easyfig` are alternatives; `easyfig` needs Easyfig installed and `$EASYFIG_PY` set. A locus GenBank is always written so you can also open the map in Easyfig/Artemis/clinker yourself. |
| `--no-gene-labels` | flag (`dest=gene_labels`, `store_false`) | labels **on** | Draw the genome map without gene-name labels (just coloured arrows). |
| `--palette {default,colorblind}` | choice | `default` | Gene-colour palette: `default`, or `colorblind` (Paul Tol muted, colour-blind-safe). Applies to the genome maps. |
| `--functional-labels` | flag | off | Also tag the gene of interest + its overprint partner with their functional category (e.g. `[transcription]`); colour+legend still carry function for the rest. |
| `--module-brackets` | flag | off | Bracket contiguous runs of same-category genes with the module name (e.g. `structural module`) above the map. |

### `run_pipeline.py` — zero-prompt autonomous launcher

`scripts/run_pipeline.py` does **not** use `argparse`. It parses `sys.argv` by hand (functions `_has`, `_out_dir`, `main`), consumes a few launcher-only flags, injects no-prompt defaults, and forwards every remaining flag unchanged to either `hmm_finder.py` (default) or `scan_genome.py` (`--scan`). It runs the target with the running interpreter's `bin` dir prepended to `PATH` (so the conda env's tools resolve without `conda activate`) and with `stdin=subprocess.DEVNULL` (so no step can block on `input()`). On a real run it writes the resolved command to `<out-dir>/run_command.txt`.

#### Launcher-only flags (consumed here, never forwarded)

| Flag | Effect |
|------|--------|
| `--scan` | Route to `scan_genome.py` (single-genome mode) instead of `hmm_finder.py`. Removed from the forwarded args. |
| `--preset NAME` | Prepend a named bundle of defaults (explicit user flags still win because they come after). Must be one of the `PRESETS` keys, else it exits with the valid list. `--preset smoke` also forces discovery mode (`scan_mode=False`). |
| `--list-presets` | Print the presets and exit. |
| `--dry-run` | Print the exact resolved command (with injected defaults) and exit; run nothing. |
| `-h`, `--help` | Print the module docstring and exit. Also printed when called with no args at all. |

#### Presets (the `PRESETS` dict)

| Preset | Expands to | Purpose |
|--------|-----------|---------|
| `phage-discovery` | `--all-databases --find-interrupted` | The typical overprinting search (gp75-style). |
| `discovery` | `--all-databases` | Clean family discovery, no read-through. |
| `offline` | `--all-databases --no-annotate` | No NCBI contact at all. |
| `smoke` | `--smoke` | Fast plumbing self-test. |

#### Injected no-prompt defaults

The launcher only injects a default when the user did **not** already supply an equivalent flag (`_has` matches both `--flag` and `--flag=value` forms):

- **Discovery mode (default target `hmm_finder.py`):**
  - Requires `--fasta` (exits otherwise); validates that the seed path exists and begins with a `>` FASTA header.
  - If none of `--databases / --all-databases / --pick-databases / --list-databases / --smoke` is present → appends `--all-databases` (full set, no prompt).
  - If neither `--email` nor `--no-annotate` is present → appends `--no-annotate` (offline unless an email is given).
  - If `--skip-tool-check` is absent → appends it.
- **Scan mode (`--scan`, target `scan_genome.py`):**
  - Requires one of `--genome` / `--accession` (exits otherwise) and one of `--hmm` / `--seeds` (exits otherwise).
  - No database/annotate/tool-check defaults are injected; `scan_genome.py` already exits cleanly on a missing email for `--accession`.

#### Forwarded flags

Every flag that is not a launcher-only flag is passed straight through to the target script, so any `hmm_finder.py` flag (in default mode) or any `scan_genome.py` flag (in `--scan` mode) is valid on the `run_pipeline.py` command line. Run the target scripts with `--help` for their own lists. The launcher exits 0 on success and forwards the child's non-zero exit code on failure.

### Cross-cutting notes

- **Defaults to memorize:** `--prodigal-gate` is **off** (overlap is informational, not exclusionary); `--map-tool` defaults to **`dfv`** (DNA Features Viewer); `--color-by` defaults to **`both`**; `--trans-table` defaults to **`11`** in both `hmm_finder.py` and `scan_genome.py`; `--biology-mode` defaults to **`phage`**; `hmm_finder.py --cpu` defaults to `"8"` (auto-clamped to core count), `scan_genome.py --cpu` to `4`.
- **Email handling:** `--email` is never hardcoded; both scripts fall back to `$NCBI_EMAIL`, and absence simply forces offline behavior in discovery (`hmm_finder.py`) or a clean exit when `--accession` is used in `scan_genome.py`. See [Offline behaviour](#offline-behaviour--no-annotate-and-the-email-gate).
- **Non-interactive operation:** for fully unattended runs prefer `run_pipeline.py` (detaches stdin and injects `--all-databases`, `--no-annotate`, `--skip-tool-check`), or on the bare scripts combine `--all-databases` + `--no-annotate` + `--skip-tool-check` and supply `--fasta`/`--out-dir` explicitly. Avoid `--pick-databases` and omitting `--fasta`, which prompt.

**Relevant source files:**
- `scripts/hmm_finder.py` (argparse: lines 666–742)
- `scripts/scan_genome.py` (argparse: lines 556–607)
- `scripts/run_pipeline.py` (manual argv parsing: lines 50–147)
