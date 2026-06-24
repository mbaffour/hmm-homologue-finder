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
A six-frame hit passes the ORF filter if it has **no internal stop codons** (a clean,
stop-bounded ORF that matched the profile). Overlap with a Prodigal-predicted gene is
**recorded** (`in_coding_locus`, `prodigal_*_pct`) but is **NOT exclusionary by default** —
the discovery targets are precisely the antisense / alternate-frame genes that standard
annotation misses, so requiring a Prodigal call would discard true positives. The optional
`--prodigal-gate` flag additionally requires that overlap for stricter specificity. Hits in
annotated protein databases are captured by accession and marked accordingly. Both nucleotide
and amino-acid sequences are exported with a per-hit evidence table. (The reported gene length
`orf_aa_len` is the Met-start→stop ORF; for a domain riding a long stop-free antisense frame
with no clear start, `domain_aa_len` — the conserved HMM envelope — is the length to cite.)

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
  and rendered as anchored, orthogroup-coloured publication panels (PNG/SVG/PDF).
  An ordered gene-neighbourhood table (`neighbour_gene_annotations.csv`) lists every
  neighbour with its order relative to the gene of interest (`pos_index`: 0 = the gene,
  ± = down/upstream), position (`rel_start/rel_end`, `distance_to_anchor_bp`), strand
  vs. the gene, and function — so the bordering genes can be described or manually
  labelled. A locus GenBank file is **always** written (so the map can also be opened
  in Easyfig / Artemis / clinker / pyGenomeViz directly), each accompanied by a
  **genome-map figure** (`<name>_genome_map.{png,svg,pdf}`; PNG at 300 dpi). The default
  renderer is **DNA Features Viewer** (`--map-tool dfv`, the Edinburgh Genome Foundry
  library): clean strand arrows, a real genome-coordinate axis spanning the whole locus,
  overlapping genes automatically stacked onto their own level (so an overprint partner
  never hides the gene of interest), and de-overlapped label boxes with leader lines.
  Genes are coloured by broad functional category (the synteny colour scheme used
  elsewhere); the gene of interest (the HMM hit) is bold **gold**. A smart label policy
  keeps any locus legible — the gene of interest and any gene overlapping it are always
  labelled, while other genes are labelled only when the locus is small enough to stay
  clean (so a ~279-gene phage genome is not a wall of text); all gene-name labels can be
  turned off with `--no-gene-labels`. Other `--map-tool` choices are `pub` (the built-in
  matplotlib diagram — strand shown as arrow direction, overlapping genes packed onto
  separate lanes, a full-length coordinate ruler; always available and used as the
  automatic fallback when a chosen renderer is not installed), `pygenomeviz`, and
  `easyfig` (needs Easyfig installed + `$EASYFIG_PY` set); any unavailable renderer falls
  back to `pub`. This is the same linear genome map the single-genome scan produces
  (`scan_genome_map_*`: a controllable window via `--flanks N` and a whole-contig view).
  The database discovery run emits a per-hit genome map the same way and honours the
  `GENOME_MAP_TOOL` environment variable as a renderer override.
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
- **Positive** — the seed set itself. Sensitivity = fraction of seeds recovered at the
  strict threshold (expected ≈1.0). **This is a self-recovery check** (the seeds built the
  model), NOT independent validation of distant-homolog detection.
- **Negative** — composition-matched shuffled seeds (same amino-acid composition,
  randomised order) plus, when present, curated unrelated-proteome sets
  (reviewed Swiss-Prot fungi/mammalian/archaea, fetched once with
  `--download-controls`). The false-positive rate is the fraction of negative
  sequences scoring ≥ strict.

**Scope of these controls (important, and what to state in a paper).** The negatives are
unrelated **protein** proteomes searched directly with `hmmsearch`; they measure **specificity
against unrelated proteins**, i.e. that hits are not an artefact of amino-acid composition or
generic similarity. They do **NOT** estimate a **false-discovery rate over the six-frame /
read-through search space**, where this pipeline's actual false positives (spurious read-through
ORFs in genomic sequence) would arise — no negative is a six-frame translation of phage genomes,
and the pipeline computes no decoy/q-value/empirical FDR. A high AUC / specificity here therefore
means "the profile cleanly separates family members from unrelated proteins", not "the genomic
six-frame discovery has FDR≈0". A genome-space FDR would require a decoy (e.g. reversed or
codon-shuffled phage genomes run through the same six-frame/read-through path) — a recommended
addition for a methods paper. Sensitivity, specificity, and FPR are written to
`controls/control_report.json` + `controls_summary.csv`, `run_manifest.json`
(`threshold_calibration`) and `METHODS.md`. (Disable with `--no-controls`.)

The positive and negative score distributions are also summarised as an **ROC
curve** (`controls/roc_curve.{png,svg,pdf}`): the area under the curve (AUC, exact
Mann–Whitney form) measures how cleanly the profile separates true family members
from negatives, and the **Youden's-J optimal** bit-score cutoff (the maximum-margin
threshold) is reported alongside the fixed strict threshold. The ROC is **advisory**
— it shows whether the fixed strict threshold sits within the separating gap; the
pipeline keeps the fixed strict/moderate tiers so results stay comparable across
runs. **Caveat:** when no negative scores above the noise floor (the usual case for a
specific family), the FPR is 0 at every finite cutoff and the "optimal" threshold is merely
the midpoint of the 0-bit floor and the weakest seed — it carries no negative-distribution
information. This is flagged in `control_report.json` (`optimal_threshold_defined: false`) and
the value should be read as advisory/upper-bounded, never as a data-driven calibrated optimum.
AUC and the optimal cutoff are recorded in `run_manifest.json`
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
- **Validated** — every reported hit is a real, stop-bounded ORF (no internal stops; Prodigal
  coding-locus overlap recorded but not required by default — see §5/`--prodigal-gate`), with
  both DNA and protein sequence recorded.

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
- **Overprinting test is necessary, not sufficient — and the signal is the OPEN frame,
  not the silent stop.** `overprinting_support=strong` requires that the overlapping antisense
  frame is **fully open across the whole domain** (`antisense_open_stops=0`) **and** that the
  premature stop is synonymous in it. The discriminating, length-dependent component is the
  open-frame condition (improbable by chance for a long domain, ~0.7 % at the 137-aa gp75
  length). A single stop being synonymous in the antisense frame is, on its own, expected
  **~85–100 % of the time from the genetic-code geometry** (the base that removes a stop tends
  to land on a synonymous antisense-wobble position), so the silent-stop clause alone carries
  little information — do not present silentness as standalone "direct evidence". Even with both
  conditions, this does **not** prove the antisense frame is a transcribed, translated, selected
  gene; confirming expression needs orthogonal data (RNA-seq/ribo-seq, conservation of the
  antisense ORF, dN/dS), which the tool does not generate.
- **Interrupted-scan threshold is heuristic for arbitrary families.** The
  read-through reporting bar is `max(30 bits, ROC-Youden)`. The 30-bit floor is a
  fixed heuristic (the validation lenient bound), and the ROC cutoff is calibrated
  on the *stop-to-stop* control set, then transferred to the larger read-through
  space; it is not re-derived against a read-through-specific null. It is validated
  on gp75 (where the floor dominates); for a family where the ROC term would
  dominate, treat low-scoring interrupted candidates with extra caution.
- **Genetic code & target domain.** Translation defaults to code **11**
  (bacterial/phage); set `--trans-table` for others. Note `--trans-table` reassigns the
  **seed** translation and the standard codon table; the read-through / six-frame stop set is
  the standard `TAA/TAG/TGA` (identical for code 11), so a stop-reassigning code (e.g. 4/25)
  would need the read-through path adjusted. The database catalog and controls are tuned for
  **phage/viral** discovery — usable on other families, but the curated databases are viral.
- **Single-genome scan threshold.** `scan_genome.py` reports six-frame/read-through matches above
  a **fixed 25-bit floor** (independent of the database run's ROC-Youden calibration); it is a
  per-genome presence/absence tool, so treat low-scoring single-genome calls with the same
  caution as low-scoring interrupted candidates.
- **Candidate homologs, not function.** A validated ORF with HMM homology is a
  *candidate*; biological function still requires experimental validation. The tool
  does no wet-lab/primer/expression design.
- **Optional components degrade gracefully.** NCBI annotation (organism names,
  GenBank neighbourhoods) needs network + an email and is skipped offline; static
  clinker PNGs need a headless browser (`playwright install chromium`) and are
  skipped if absent — the run still completes and the static synteny panels are
  produced regardless.
