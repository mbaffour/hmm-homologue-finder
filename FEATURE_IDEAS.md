# Feature ideas

A running list of proposed, **not-yet-implemented** enhancements for HMM
Homologue Finder. Each entry states the idea, why it helps, and the files it
would most likely touch. These are proposals for discussion — nothing here has
been implemented, and none of them should be built without validating the
behaviour against a real run first.

---

## 1. Resumable runs (`--resume`)

**What.** Let an interrupted discovery run pick up where it stopped instead of
restarting the multi-hour, multi-database search from scratch.

**Why.** A full run performs a six-frame search across up to ten databases,
ORF validation, iteration, clustering, tree building, and figure generation. If
the process dies late (network drop, host reboot, out-of-memory during figures)
the user currently re-runs everything. The pipeline already records per-step
completion, so most of the machinery for a safe resume already exists.

**Touch-points.**
- `engine/core/state.py` — `PipelineState` already persists step completion and
  params to `.pipeline_state.json` (`mark_complete` / `is_complete` /
  `get_params`). A resume would read this to skip completed steps.
- `scripts/hmm_finder.py` — add a `--resume` flag alongside the existing
  `--no-overwrite` (around line 780) and gate each stage on
  `state.is_complete(...)`.
- `engine/core/runner.py` — thread the resume decision through the stage driver.

---

## 2. HTML summary report

**What.** A single self-contained `summary_report.html` that a reader can open
in a browser: run settings, per-database hit counts, confidence-tier breakdown,
and links (or thumbnails) for the synteny figures and tree.

**Why.** The pipeline already produces a `reproducibility.json` audit record and
a compact run summary; a rendered HTML view makes results shareable with
collaborators who won't open a terminal, which suits a portfolio/publication
context.

**Touch-points.**
- `engine/pipeline/reporter.py` — already renders `reports/summary_report.html`
  from Jinja2 and writes `reproducibility.json`; extend the context passed to
  the template rather than adding a new renderer.
- `engine/www/templates/report.html.j2` — the existing template to enrich.
- `engine/core/run_summary.py` — source of the factual per-run fields
  (inputs, settings, database statuses, hit counts, confidence tiers).

---

## 3. Machine-readable JSON output (`--json-summary`)

**What.** Emit one stable, documented, top-level JSON summary of a run (hits per
database, confidence tiers, output-file paths, versions) intended for
downstream tooling.

**Why.** `reproducibility.json` exists as an audit record, but a small, stable,
explicitly-versioned summary schema lets other scripts, notebooks, or a CI job
diff runs and assert on hit counts without parsing HTML or scraping logs.

**Touch-points.**
- `engine/core/run_summary.py` — already assembles the summary dict from files
  the app/benchmark wrote; add a documented, versioned serialization.
- `scripts/hmm_finder.py` — add a `--json-summary PATH` flag near the other
  output flags (around lines 779–800) to choose the destination.
- `scripts/export_csv.py` — sibling exporter to keep field names consistent
  with.

---

## 4. Parallelize the per-database search (`--threads`)

**What.** Search multiple databases concurrently, in addition to the existing
per-search HMMER threading.

**Why.** `engine/pipeline/searcher.py` already passes `--cpu` to HMMER for a
single search (see the `cpu` handling around lines 78–79), but the ten
databases are searched one after another. On a multi-core host, overlapping the
independent per-database searches is the biggest wall-clock win available.

**Touch-points.**
- `scripts/hmm_finder.py` — add a `--threads N` flag (distinct from the existing
  `--cpu`, which is HMMER threads-per-search) near line 784.
- `engine/pipeline/searcher.py` — the per-database search entry point that would
  be dispatched across a bounded worker pool.
- `engine/core/runner.py` — where the database loop is driven; the natural place
  to fan out and join.

_Care: bound the pool and keep the aggregate thread count (`--threads` ×
`--cpu`) under the available cores — `hmm_finder.py` already clamps `--cpu` to
the host core count (around line 872), so `--threads` should respect the same
budget._

---

## 5. Config presets (`--preset`)

**What.** Named bundles of settings (e.g. `quick`, `standard`, `exhaustive`)
that expand to the corresponding flag values, overridable individually.

**Why.** A run has many knobs (`--iterations`, `--databases`, `--all-databases`,
`--color-by`, `--biology-mode`, control handling, …). Presets give newcomers a
sensible default and make documented runs reproducible without a long command
line — a good fit for a portfolio demo.

**Touch-points.**
- `scripts/hmm_finder.py` — add a `--preset {quick,standard,exhaustive}` flag
  and expand it into defaults before argument resolution (the argparse block
  begins around line 776); explicit flags should still win over the preset.
- `docs/USAGE.md` — document each preset and the flags it sets.
- `README.md` — mention the presets in the Quick start section.

---

_See also: the three HMM repositories (`hmm-homologue-finder`, `hmm-discovery`,
`hmm-discovery-app`) share the `pipeline/` engine. Extracting that shared engine
into a single installable package that all three consume would stop the copies
drifting apart (the same class of problem as the CI dependency drift fixed in
this branch). This is a cross-repo refactor, noted here only as future
direction — not something to attempt inside this repo alone._
