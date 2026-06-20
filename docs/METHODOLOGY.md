# Methodology

This describes the method the pipeline implements, generically. Replace the
example seed set with your own protein family; the workflow is identical.

## 1. Overview
Distant homologues of a protein family are discovered by iterative profile-HMM
search of public phage/viral sequence databases. Genome databases are translated
in all six reading frames before search, so homologues encoded by genes that
standard annotation does not predict can be recovered. For every hit, the matched
open reading frame (ORF) is reconstructed from its genomic coordinates, validated
as a genuine ORF, and exported as both nucleotide and amino-acid sequence. The
discovered, deduplicated sequences re-seed further search rounds to test for
convergence.

## 2. Seed sequences
A curated set of family protein sequences (the only manual input). Quality
matters more than quantity: full-length, genuinely related members give a sharper
profile than many fragments.

## 3. Profile HMM construction
Seeds are aligned with **MAFFT** using an accuracy-first strategy — **L-INS-i**
(`--localpair --maxiterate 1000`, the gold standard for homologous domains in
variable-length context) for tractable seed counts (≤500 sequences), falling back
to `--auto` (FFT-NS-2/PartTree) only for very large sets — then trimmed with
**trimAl** (`automated1`). A profile HMM is built with **hmmbuild** (HMMER 3.4)
and validated by self-search against the seeds; a run proceeds only if seed
recovery exceeds a threshold (default 70%).

## 4. Database search
The profile is searched (**hmmsearch**, E ≤ 1 × 10⁻⁵) against the selected
databases. Each download is integrity-checked (SHA-256) with its URL and access
date recorded. Nucleotide databases are translated in all six frames into
stop-to-stop ORFs (minimum length 30 aa) and searched with the protein profile.
Large databases are chunked with **seqkit** and processed in parallel.

## 5. ORF-validated sequence extraction
For each hit the ORF is reconstructed directly from the genomic coordinates
(1-based, strand/frame-correct) and translated; the family domain within it is
delimited by the HMM envelope (`hmmsearch --domtblout`). Recorded per hit:
full-ORF length, domain length and coverage, internal-stop count (required to be
0), and overlap with **Prodigal** gene predictions (same-strand and any-strand).
A hit passes the ORF filter if it has no internal stop codons and sits within a
genuine coding locus. Hits in annotated protein databases are captured by
accession and marked accordingly. Both nucleotide and amino-acid sequences are
exported with a per-hit evidence table.

## 6. Iterative refinement and convergence
Unique, ORF-validated domains from one round seed the next; identical databases,
parameters, and extraction are applied each round. Iteration **stops early on
convergence**: when the unique-validated-hit count changes by <5 % AND the HMM
length (match states) changes by <3 between consecutive rounds, the detectable
family is considered fully recovered and no further round is run. A round that
yields zero validated hits also stops iteration. The stopping reason is recorded
in `run_manifest.json` (`iteration_stop_reason`) and `METHODS.md`. (Identical hit
sequences recurring across many genomes are deduplicated before re-seeding, so
the unique-sequence count is the meaningful measure of family diversity.) The
most complete round (most validated hits — after convergence, the refined final
round) is the canonical set used for the figures, the published profile HMM, and
the main paper table, so tables and figures describe the same homolog set.

## 7. Downstream characterisation
- **Clustering** — CD-HIT (40% identity, 80% coverage).
- **Synteny** — Prodigal gene calls provide flanking-gene context; neighbourhoods
  are compared per cluster with clinker; real-sequence GenBank files are written.
- **Alignment of homologs** — the unique ORF-validated domains are aligned with
  the same accuracy-first MAFFT strategy (L-INS-i where tractable, else `--auto`).
  The alignment is a first-class deliverable: the full MSA (`hits.aln.faa`), a
  trimmed copy, quality statistics (`hits.aln.stats.json`: length, gap %, conserved
  columns, mean pairwise identity), and a publication-ready ClustalX-coloured
  figure (`alignment_figure.{png,svg,pdf}`) are written and embedded in the report.
- **Phylogenetics** — trimAl (`-gt 0.5`) -> IQ-TREE (ModelFinder; 1000 ultrafast
  bootstrap; fixed random seed for reproducibility). The deliverable tree is built
  once, after the runs, on the most-complete run's homologs **with the seeds
  included and marked (`SEED_*`)** so the reader sees where the starting sequences
  fall among everything discovered. A separate seed-only QC tree + alignment is
  built once *before* the runs (skippable with `--no-seed-tree`) to sanity-check
  the input set. No per-iteration trees are built (they answer no scientific
  question and waste compute).
- **Motifs** — MEME (<=3 motifs, width 6-30 aa); scanned with FIMO.

## 7b. Threshold calibration (controls)
The bit-score thresholds used to tier hits (strict 45, moderate 30) are calibrated
per run against built-in controls on the profile HMM:
- **Positive** — the seed set itself (sensitivity = fraction of seeds recovered at
  the strict threshold; expected ≈1.0).
- **Negative** — composition-matched shuffled seeds (same amino-acid composition,
  randomised order) plus, when present, curated unrelated-proteome sets
  (reviewed Swiss-Prot fungi/mammalian/archaea, fetched once with
  `--download-controls`). The false-positive rate is the fraction of negative
  sequences scoring ≥ strict.
Sensitivity, specificity, and false-positive rate are written to
`controls/control_report.json` + `controls_summary.csv` and summarised in
`run_manifest.json` (`threshold_calibration`) and `METHODS.md`. This quantifies
that hits above threshold are not an artefact of amino-acid composition or generic
similarity to unrelated proteins. (Disable with `--no-controls`.)

## 8. Reproducibility
The whole workflow runs from a single command requiring only a seed FASTA. Tool
versions, database URLs, access dates, and checksums are recorded in the run's
`reproducibility.json`.

Software: HMMER 3.4, MAFFT v7.526, trimAl v1.5, Prodigal V2.6.3, seqkit v2.13.0,
CD-HIT 4.8.1, IQ-TREE 3.1.2, MEME/FIMO 5.5.9, clinker v0.0.32; genome retrieval
via NCBI Entrez and direct catalogue streaming.

## 9. Interpreting results
- **Converged** — hit counts stop growing between rounds.
- **Novel & specific** — zero hits in reviewed-protein and domain databases
  (SwissProt, Pfam, VOGDB) across rounds, with hits found only via six-frame
  translation of genome databases.
- **Validated** — every reported hit is a real ORF (no internal stops; in a
  coding locus), with both DNA and protein sequence recorded.
