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

**Interrupted / overprinted homologs (optional, `--find-interrupted`).** Because
six-frame ORFs are stop-to-stop, a homolog whose gene carries a *premature stop*
is truncated at the stop and missed — which is precisely how an overprinted gene
behaves (a nonsense mutation in the overprinted gene can be synonymous in the
overlapping reading frame, so it is tolerated). With this flag the nucleotide
databases are additionally translated with **read-through** (stop codons retained)
and searched with the family HMM; matches whose domain envelope contains ≥1
internal stop are candidate interrupted/overprinted homologs the stop-to-stop
search cannot see. For each candidate the scan records (`interrupted_homologs.tsv`):
the genome coordinates and DNA of the matched domain; the per-stop genome
coordinates (`stop_nt_positions`); the **full read-through ORF** read *through* the
premature stop(s) to the natural gene end, as protein (`full_orf_aa`, internal
stops shown as `*`) and as nucleotide (`full_orf_nt`, coding 5'→3', ending in the
actual stop codon triplet — `natural_stop_nt` gives its genome coordinate); and
three FASTAs (`interrupted_homologs_domain_aa.faa`, `…_full_orf_aa.faa`,
`…_full_orf_nt.fna`). The **reporting threshold is family-calibrated**:
`max(30 bits, the ROC Youden-optimal cutoff)` from the same run controls — never
below the run's lenient evidence bar and raised to the family's calibrated noise
floor when controls were run (the read-through scan covers a much larger, noisier
space than the stop-to-stop search, so the bar only ever tightens); with
`--no-controls` the bare 30-bit floor is used.

*Overprinting (silent-stop) test — the proof, not just the location.* Locating a
premature stop is necessary but not sufficient to call a gene overprinted. For
each candidate the scan therefore picks the antisense frame with the **fewest
stops across the domain** (the candidate overlapping ORF; `antisense_open_frame`,
`antisense_open_stops` — 0 = a fully open overlapping reading frame) and tests
whether each premature stop is **synonymous** in that frame — i.e. whether some
single-base change that reverts the nonsense mutation in *this* gene would leave
the antisense protein unchanged (`stop_silent_antisense`, per-stop). The verdict
`overprinting_support` is **strong** when the antisense frame is fully open *and*
every premature stop is synonymous in it (direct evidence the gene is overprinted
antisense to another gene), **partial** when only some stops are silent, **none**
otherwise. This is a *necessary* signature — it confirms synonymy in an open
overlapping frame — but does not by itself prove the antisense ORF is expressed
(see Limitations).

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
  are compared per cluster with clinker (interactive `cluster_*.html`; a static
  `cluster_*.png` is also exported per cluster when a headless browser is installed)
  and rendered as anchored, orthogroup-coloured publication panels (PNG/SVG/PDF);
  real-sequence GenBank files are written.
- **Alignment of homologs** — the unique ORF-validated domains are aligned with
  the same accuracy-first MAFFT strategy (L-INS-i where tractable, else `--auto`).
  The alignment is a first-class deliverable: the full MSA (`hits.aln.faa`), a
  trimmed copy, quality statistics (`hits.aln.stats.json`: length, gap %, conserved
  columns, mean pairwise identity), and a publication-ready ClustalX-coloured
  figure (`alignment_figure.{png,svg,pdf}`) are written and embedded in the report.
  The report also embeds an inline, residue-coloured view of the first hits.
  Separately, every unique homolog is aligned **to the family HMM itself** with
  `hmmalign` (`hits_hmmalign.sto` Stockholm + `hits_hmmalign.a2m` A2M) so each hit's
  match states vs insertions relative to the model are explicit.
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

The positive and negative score distributions are also summarised as an **ROC
curve** (`controls/roc_curve.{png,svg,pdf}`): the area under the curve (AUC, exact
Mann–Whitney form) measures how cleanly the profile separates true family members
from negatives, and the **Youden's-J optimal** bit-score cutoff (the maximum-margin
threshold) is reported alongside the fixed strict threshold. The ROC is **advisory**
— it shows whether the fixed strict threshold sits within the separating gap; the
pipeline keeps the fixed strict/moderate tiers so results stay comparable across
runs. AUC and the optimal cutoff are recorded in `run_manifest.json`
(`threshold_calibration.roc`) and `METHODS.md`.

## 8. Reproducibility
The whole workflow runs from a single command requiring only a seed FASTA. Tool
versions, database URLs, access dates, and checksums are recorded in the run's
`reproducibility.json`.

Software: HMMER 3.4, MAFFT v7.526, trimAl v1.5.1, Prodigal V2.6.3, seqkit v2.13.0,
CD-HIT 4.8.1, IQ-TREE 3.1.2, MEME/FIMO 5.5.9, clinker v0.0.32; genome retrieval
via NCBI Entrez and direct catalogue streaming.

## 9. Interpreting results
- **Converged** — hit counts stop growing between rounds.
- **Novel & specific** — zero hits in reviewed-protein and domain databases
  (SwissProt, Pfam, VOGDB) across rounds, with hits found only via six-frame
  translation of genome databases.
- **Validated** — every reported hit is a real ORF (no internal stops; in a
  coding locus), with both DNA and protein sequence recorded.

## 10. Limitations / scope (what it does *not* do)
The tool is a sequence-homology discovery pipeline; reading these bounds prevents
over-interpreting its output.

- **Sequence/HMM homology only — no structural search.** Detection is profile-HMM
  based (down to the ~15–25 % "twilight zone"), not structure-based. Homologs whose
  sequence has diverged past HMM detectability but whose fold is conserved (the
  domain of tools like Foldseek/DALI) are out of scope.
- **Assembled databases, not raw reads.** Input is a seed FASTA searched against
  assembled genome/protein databases. **Read-level data is not an input.**
  *RNA-seq / read-based evidence (e.g. expression or transcript support) is planned
  future work and is **not wired in now**.*
- **Overprinting test is necessary, not sufficient.** `overprinting_support=strong`
  confirms the premature stop is synonymous in an *open* overlapping antisense
  frame — strong sequence evidence of overprinting — but it does **not** prove that
  antisense frame is a transcribed, translated, selected gene. Confirming
  expression needs orthogonal data (e.g. RNA-seq/ribo-seq, conservation of the
  antisense ORF, dN/dS), which the tool does not generate.
- **Interrupted-scan threshold is heuristic for arbitrary families.** The
  read-through reporting bar is `max(30 bits, ROC-Youden)`. The 30-bit floor is a
  fixed heuristic (the validation lenient bound), and the ROC cutoff is calibrated
  on the *stop-to-stop* control set, then transferred to the larger read-through
  space; it is not re-derived against a read-through-specific null. It is validated
  on gp75 (where the floor dominates); for a family where the ROC term would
  dominate, treat low-scoring interrupted candidates with extra caution.
- **Genetic code & target domain.** Translation defaults to code **11**
  (bacterial/phage); set `--trans-table` for others. The database catalog and
  controls are tuned for **phage/viral** discovery — usable on other families, but
  the curated databases are viral.
- **Candidate homologs, not function.** A validated ORF with HMM homology is a
  *candidate*; biological function still requires experimental validation. The tool
  does no wet-lab/primer/expression design.
- **Optional components degrade gracefully.** NCBI annotation (organism names,
  GenBank neighbourhoods) needs network + an email and is skipped offline; static
  clinker PNGs need a headless browser (`playwright install chromium`) and are
  skipped if absent — the run still completes and the static synteny panels are
  produced regardless.
