# Adding databases to the HMM Homologue Finder

The search engine streams **any gzipped FASTA over https** — it downloads and
searches on the fly, nucleotide databases are six-frame translated. To make a
database selectable with `--databases "<name>,…"`, add an entry to the registry:

```
~/Documents/HMM-Discovery-Deployable-20260602/databases.json
```

## Entry schema (copy an existing one and edit)
```json
{
  "name": "My Database",
  "type": "nucleotide",            // or "protein"
  "download_url": "https://host/path/file.fasta.gz",
  "streaming": true,
  "enabled": true,
  "optional": true,
  "size_hint": "~2 GB",
  "est_time": "10-20 min",
  "notes": "what it is",
  "relevance": "why you'd search it",
  "search_mode": "hmmsearch",      // "hmmscan" for HMM libraries (needs a setup handler)
  "setup_handler": null
}
```
- **Nucleotide** databases are translated in all six frames automatically.
- **Protein** FASTA databases are searched directly; their hits are now captured
  by accession (`source_type=annotated_protein`).
- A plain gzipped-FASTA URL is all that's required — no other code changes.

## Currently available (11 in the registry)
Stream: INPHARED genomes, INPHARED proteins, SwissProt, RefSeq viral proteins,
RefSeq viral genomes, GPD, GVD-AVrC, Pfam (sequences), **RefSeq bacterial
proteins**. Download-once HMM libraries: Pfam (domain scan), VOGDB VFAM.

> **RefSeq bacterial proteins** is ready to use now — just include it:
> `--databases "INPHARED genomes,RefSeq viral genomes,RefSeq bacterial proteins"`
> (≈80 GB stream; good as a host/background specificity check.)

## PHROGs — ✅ WIRED IN & TESTED (selectable now)

PHROGs (38,880 prokaryotic-virus protein HMM families) is now a built-in
`hmmscan` annotation database. Use it in any run:
```bash
python3 scripts/hmm_finder.py --fasta my_seeds.faa \
    --databases "INPHARED genomes,RefSeq viral genomes,PHROGs (annotation)"
```
- **Source (verified, reliable):** the Pharokka v1.8.0 database bundle on Zenodo
  (656 MB) — it ships the PHROGs models as a pre-pressed `all_phrogs.h3m`:
  `https://zenodo.org/records/17110353/files/pharokka_v1.8.0_databases.tar.gz?download=1`
- **How it was wired in:** the engine's hmmscan handler
  (`run_all_database_benchmark.py`) now falls back to a pre-pressed `.h3m` when a
  bundle has no raw `.hmm` files — it extracts the largest `.h3m`, runs
  `hmmpress` to rebuild the aux files (verified: HMMER reads the binary `.h3m`
  fine), then `hmmscan`. Registered in `databases.json` as `PHROGs (annotation)`.
- **Tested:** download → extract → hmmpress (38,880 models) → hmmscan all
  succeeded; runs in seconds against a small protein set.
- *Note:* the official site `phrogs.lmge.uca.fr` has an expired SSL cert and the
  old CLIMB mirror 404s, which is why the reliable Zenodo/Pharokka bundle is used.

## Still needs your input — IMG/VR (JGI uncultivated viral genomes)
- Distributed **only through the JGI Data Portal** (login/Globus); there is **no
  public gzipped-FASTA URL** to stream. It cannot be auto-fetched.
- To use it: download a release with a free JGI account (e.g. the
  high-confidence genomes FASTA), then add a local-file entry:
  ```json
  { "name": "IMG/VR genomes", "type": "nucleotide",
    "download_url": "file:///Users/you/Downloads/IMGVR_high_confidence.fna.gz",
    "streaming": true, "enabled": true, "optional": true,
    "search_mode": "hmmsearch", "setup_handler": null }
  ```
  (the engine reads a `file://` path the same way it streams a URL).

Both are a few minutes' work once the input is in hand — the engine already
supports `hmmsearch` (FASTA) and `hmmscan` (HMM library) modes.
