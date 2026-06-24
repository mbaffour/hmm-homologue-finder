#!/usr/bin/env python3
"""scan_genome.py — does ONE genome contain my gene?

A focused, single-genome counterpart to the full discovery pipeline. Build (or
accept) a profile HMM for your gene of interest and scan **one** genome for it,
reporting a per-genome result: present / not detected, the hit's genome
coordinates, ORF validation (start M, ends at a stop, internal stops), the HMM
score, and the hit sequence (protein + DNA). Because the scan is read-through
(stop codons kept, not broken on), a clean copy is found with 0 internal stops
**and** a stop-interrupted / overprinted copy is found too (reported only with
``--find-interrupted``), with the same overprinting (silent-stop) test the
discovery pipeline uses.

    # from seed sequences (builds the HMM for you; nucleotide seeds are translated):
    python3 scan_genome.py --seeds gene_seeds.faa --genome genome.fna --out scan_out

    # from an existing profile HMM:
    python3 scan_genome.py --hmm gene.hmm --genome genome.fna --out scan_out

    # fetch the genome from NCBI by accession instead of a local file:
    python3 scan_genome.py --hmm gene.hmm --accession KX098390 --email you@inst.edu

    # also report stop-interrupted / overprinted copies (with the overprinting test):
    python3 scan_genome.py --seeds gene_seeds.faa --genome genome.fna --find-interrupted

The genes around each hit are written to scan_neighbourhood.csv (ordered relative to your
gene: pos_index 0 = your gene, ± up/downstream; a relationship column flags upstream /
downstream / overlapping) and drawn as a genome-map figure scan_genome_map_<hit>.png/svg.
The --flanks genes each side are shown contiguously AND every overlapping gene is kept,
so an overprint partner (e.g. gp75's antisense RNA polymerase) is shown, not hidden.
Names come from the genome's OWN annotation (gene / gp number, product, locus_tag,
protein_id) for an annotated record (a fetched --accession, or a .gb/.gbk/.gbff
--genome); otherwise via Prodigal + optional VOGDB VFAM. Disable with --no-neighbours.
"""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import socket
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

import find_interrupted as FI   # read-through translation, coord mapping, overprinting

ROW_COLS = ["contig", "strand", "frame", "nt_start", "nt_end", "domain_aa_len",
            "internal_stops", "status", "domain_bit_score", "i_evalue",
            "orf_nt_start", "orf_nt_end", "orf_aa_len", "has_start_M", "ends_at_stop",
            "overprinting_support", "antisense_open_stops", "stop_nt_positions",
            "domain_aa", "orf_aa", "orf_nt"]


def _efetch(ids, rettype, email, log):
    from Bio import Entrez
    Entrez.email = email
    socket.setdefaulttimeout(120)
    for attempt in (1, 2, 3):
        try:
            with Entrez.efetch(db="nucleotide", id=",".join(ids),
                               rettype=rettype, retmode="text") as h:
                data = h.read()
            if data and data.strip():
                return data
        except Exception as e:
            log(f"  (NCBI {rettype} fetch attempt {attempt} failed: {e})")
        time.sleep(3 * attempt)
    return ""


def _genbank_to_fasta(recs, fa_path: Path):
    with open(fa_path, "w") as f:
        for rec in recs:
            f.write(f">{rec.id}\n{str(rec.seq)}\n")


def fetch_genome(accessions: str, email: str | None, out_dir: Path, log):
    """Pull NCBI nucleotide accession(s) (comma-separated) as a **GenBank** record so the
    genome's OWN gene annotation (names, products, locus tags, protein IDs, gp numbers)
    is available to describe the neighbours; extract the nucleotide FASTA for scanning.
    Returns (fasta_path, genbank_path_or_None). Falls back to FASTA-only. Needs an email."""
    email = email or os.environ.get("NCBI_EMAIL")
    if not email:
        sys.exit("--accession needs an email for NCBI Entrez: pass --email you@inst.edu "
                 "(or set $NCBI_EMAIL).")
    ids = [a.strip() for a in accessions.replace(",", " ").split() if a.strip()]
    if not ids:
        sys.exit("no accession given to --accession")
    stem = ids[0].replace("/", "_") if len(ids) == 1 else "fetched_genome"
    fa_path, gb_path = out_dir / f"{stem}.fna", out_dir / f"{stem}.gb"
    log(f"Fetching {len(ids)} accession(s) from NCBI (GenBank + annotation): {', '.join(ids)} …")
    gb = _efetch(ids, "gbwithparts", email, log)
    if gb.lstrip().startswith("LOCUS"):
        gb_path.write_text(gb)
        try:
            recs = list(SeqIO.parse(str(gb_path), "genbank"))
        except Exception:
            recs = []
        if recs and any(len(r.seq) for r in recs):
            _genbank_to_fasta(recs, fa_path)
            n_cds = sum(1 for r in recs for ft in r.features if ft.type == "CDS")
            log(f"  {len(recs)} record(s), {n_cds} annotated CDS -> {fa_path.name}"
                + (f" (+ {gb_path.name}, using the genome's own gene names)" if n_cds else ""))
            return fa_path, (gb_path if n_cds else None)
    # fallback: sequence-only FASTA (no annotation) -> Prodigal will call neighbours
    log("  (no usable GenBank annotation; fetching FASTA only — neighbours via Prodigal)")
    fa = _efetch(ids, "fasta", email, log)
    if not fa.strip().startswith(">"):
        sys.exit(f"NCBI returned nothing for: {accessions} (is it a nucleotide record? "
                 "assembly GCF_/GCA_ ids aren't fetched directly — use their contig accessions).")
    fa_path.write_text(fa)
    log(f"  fetched {fa.count('>')} sequence(s) -> {fa_path.name}")
    return fa_path, None


def genbank_genes(gb_path) -> dict:
    """{contig_id: [(start_1based, end, strand, meta)]} from a GenBank file's CDS features.
    meta carries the genome's OWN annotation: gene (or locus_tag) — often the gp number —
    product, locus_tag, protein_id."""
    out = {}
    try:
        recs = list(SeqIO.parse(str(gb_path), "genbank"))
    except Exception:
        return out
    for rec in recs:
        genes = []
        for ft in rec.features:
            if ft.type != "CDS":
                continue
            try:
                s, e = int(ft.location.start) + 1, int(ft.location.end)
            except Exception:
                continue
            q = ft.qualifiers
            genes.append((s, e, -1 if ft.location.strand == -1 else 1, {
                "gene": (q.get("gene", [""])[0] or q.get("locus_tag", [""])[0]),
                "product": q.get("product", [""])[0],
                "locus_tag": q.get("locus_tag", [""])[0],
                "protein_id": q.get("protein_id", [""])[0]}))
        genes.sort()
        out[rec.id] = genes
    return out


def genbank_organisms(gb_path) -> dict:
    """{contig_id: organism/phage name} from a GenBank file (the record's /organism)."""
    out = {}
    try:
        for rec in SeqIO.parse(str(gb_path), "genbank"):
            org = rec.annotations.get("organism") or ""
            src = next((f for f in rec.features if f.type == "source"), None)
            if not org and src:
                org = src.qualifiers.get("organism", [""])[0]
            out[rec.id] = (org or rec.id).strip()
    except Exception:
        pass
    return out


def resolve_local_genome(path: Path, out_dir: Path, log):
    """Return (fasta_path, genbank_path_or_None) for a local genome. A GenBank input
    (.gb/.gbk/.gbff, or a file starting with LOCUS) is used for its annotation and its
    sequence is extracted to FASTA; a plain FASTA is used as-is (neighbours via Prodigal)."""
    path = Path(path)
    is_gb = path.suffix.lower() in (".gb", ".gbk", ".gbff", ".genbank", ".gbf")
    if not is_gb:
        try:
            with FI.open_maybe_gz(path) as fh:
                is_gb = fh.read(64).lstrip().startswith("LOCUS")
        except Exception:
            is_gb = False
    if not is_gb:
        return path, None
    try:
        recs = list(SeqIO.parse(str(path), "genbank"))
    except Exception:
        return path, None
    fa = out_dir / (path.stem + ".fna")
    _genbank_to_fasta(recs, fa)
    n_cds = sum(1 for r in recs for ft in r.features if ft.type == "CDS")
    log(f"GenBank input: {len(recs)} record(s), {n_cds} annotated CDS"
        + (" -> using the genome's own gene names for neighbours" if n_cds else ""))
    return fa, (path if n_cds else None)


def build_hmm_from_seeds(seeds: Path, table: int, out_dir: Path, cpu: int, log) -> Path:
    """Build a profile HMM from seed sequences (translating a nucleotide seed first).
    MAFFT-align (accuracy-first L-INS-i for ≤500 seqs, else --auto) then hmmbuild;
    a single seed is built directly. Same toolchain as the discovery pipeline."""
    import hmm_finder as H   # sibling: reuse the nucleotide-seed detection + translation
    seeds = Path(seeds)
    if H._looks_like_nucleotide(seeds):
        log(f"Seed looks like nucleotide CDS; translating with genetic code {table}…")
        seeds = H.translate_seed(seeds, table, out_dir, log)
    n = sum(1 for _ in SeqIO.parse(str(seeds), "fasta"))
    if n == 0:
        sys.exit("No sequences found in the seed FASTA.")
    aln = out_dir / "seeds.aln.faa"
    if n == 1:
        shutil.copy(str(seeds), str(aln))           # hmmbuild accepts a single sequence
    else:
        mode = (["--localpair", "--maxiterate", "1000"] if n <= 500 else ["--auto"])
        with open(aln, "w") as f:
            subprocess.run(["mafft", *mode, "--thread", str(cpu), str(seeds)],
                           check=True, stdout=f, stderr=subprocess.DEVNULL)
    hmm = out_dir / "profile.hmm"
    subprocess.run(["hmmbuild", "--amino", str(hmm), str(aln)],
                   check=True, capture_output=True, text=True)
    log(f"Built HMM from {n} seed(s) -> {hmm.name}")
    return hmm


def scan(genome: Path, hmm: Path, out_dir: Path, min_bit: float, find_interrupted: bool,
         cpu: int, log, neighbours: bool = False, db_cache=None, annotation_gb=None,
         flanks: int = 7, map_tool: str = "dfv", gene_labels: bool = True) -> dict:
    """Read-through six-frame scan of one genome with the HMM. Returns a summary dict
    and writes scan_hits.tsv + scan_hits_{aa.faa,nt.fna} + scan_report.txt. When
    `neighbours` is set, also describes the flanking genes (from the genome's own
    annotation `annotation_gb` when available, else Prodigal) -> scan_neighbourhood.csv."""
    out_dir.mkdir(parents=True, exist_ok=True)
    workdir = Path(tempfile.mkdtemp(prefix="scan_", dir=str(out_dir)))
    frames, markers, contig_nt = [], {}, {}
    n_contigs = 0
    try:
        with FI.open_maybe_gz(genome) as fh:
            for rec in SeqIO.parse(fh, "fasta"):
                seq = str(rec.seq).upper().replace("U", "T")
                contig_nt[rec.id] = seq
                n_contigs += 1
                for strand, frame, search, marker in FI._frames(seq):
                    for off, wlen in FI._windows(len(search)):
                        sw = search[off:off + wlen]
                        if len(sw) < FI.MIN_AA:
                            continue
                        name = f"{rec.id}__{strand}{frame}__{off}"
                        frames.append((name, sw))
                        markers[name] = marker[off:off + wlen]
        if not frames:
            log("  (genome too short to scan in any frame)")
            return _finish(out_dir, [], min_bit, find_interrupted, n_contigs, log)
        sfa = workdir / "frames.faa"
        with open(sfa, "w") as f:
            for name, search in frames:
                f.write(f">{name}\n{search}\n")
        dt = workdir / "scan.domtbl"
        log(f"hmmsearch: {len(frames)} read-through frame window(s) over {n_contigs} contig(s)…")
        subprocess.run(["hmmsearch", "--noali", "--cpu", str(cpu), "--domtblout", str(dt),
                        str(hmm), str(sfa)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        rows = _parse(dt, markers, contig_nt, min_bit, find_interrupted)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    s = _finish(out_dir, rows, min_bit, find_interrupted, n_contigs, log)
    if rows and neighbours:
        nb = write_neighbourhoods(
            rows, contig_nt, out_dir,
            db_cache or Path("~/.cache/hmm-homologue-finder"), cpu, log,
            annotation_gb=annotation_gb, flanks=flanks, map_tool=map_tool,
            gene_labels=gene_labels)
        s["neighbourhood"] = nb
        if nb:
            with open(out_dir / "scan_report.txt", "a") as f:
                f.write(f"\nNeighbouring genes (Prodigal): {Path(nb).name} — ordered by "
                        f"position relative to your gene (pos_index 0 = your gene).\n")
    return s


def _parse(dt: Path, markers: dict, contig_nt: dict, min_bit: float,
           find_interrupted: bool) -> list:
    rows = []
    for ln in dt.read_text(errors="replace").splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        p = ln.split()
        if len(p) < 23:
            continue
        try:
            name, dom_bits, i_eval = p[0], float(p[13]), float(p[12])
            env_from, env_to = int(p[19]), int(p[20])
        except (ValueError, IndexError):
            continue
        if dom_bits < min_bit:
            continue
        marker = markers.get(name, "")
        n_stops, positions = FI.count_envelope_stops(marker, env_from, env_to)
        if n_stops > 0 and not find_interrupted:
            continue   # interrupted copy, but the user asked only for clean detection
        contig, sf_tag, off_s = name.rsplit("__", 2)
        strand, frame, offset = sf_tag[0], int(sf_tag[1]), int(off_s)
        nt = contig_nt.get(contig, "")
        nt_start, nt_end, dom_nt = FI.aa_to_nt(nt, strand, frame, offset + env_from, offset + env_to)
        orf_from, orf_to, orf_aa = FI.extend_orf(marker, env_from, env_to)
        orf_nt_start, orf_nt_end, orf_nt = FI.aa_to_nt(nt, strand, frame, offset + orf_from, offset + orf_to)
        terminal_stop = orf_aa.endswith("*")
        # ORF validation: a clean gene starts with M (after the upstream stop) and ends at a stop.
        has_M = orf_aa[:1] == "M"
        status = "clean" if n_stops == 0 else "interrupted"
        overpr, anti_open = "", ""
        stop_pos = ""
        if n_stops > 0:
            stop_fwds = [FI.stop_nt(nt, strand, frame, offset + x) for x in positions]
            stop_pos = ";".join(str(s) for s in stop_fwds)
            opa = FI.analyze_overprinting(nt, strand, nt_start, nt_end, stop_fwds)
            overpr, anti_open = opa["support"], opa["open_stops"]
        rows.append({
            "contig": contig, "strand": strand, "frame": frame,
            "nt_start": nt_start, "nt_end": nt_end,
            "domain_aa_len": env_to - env_from + 1, "internal_stops": n_stops,
            "status": status, "domain_bit_score": round(dom_bits, 1),
            "i_evalue": f"{i_eval:.2g}",
            "orf_nt_start": orf_nt_start, "orf_nt_end": orf_nt_end,
            "orf_aa_len": orf_to - orf_from + 1,
            "has_start_M": int(has_M), "ends_at_stop": int(terminal_stop),
            "overprinting_support": overpr, "antisense_open_stops": anti_open,
            "stop_nt_positions": stop_pos,
            "domain_aa": marker[env_from - 1:env_to], "orf_aa": orf_aa, "orf_nt": orf_nt,
        })
    # one read-through frame can report a domain twice across overlapping windows; keep best
    best = {}
    for r in rows:
        key = (r["contig"], r["strand"], r["frame"], r["nt_start"], r["nt_end"])
        if key not in best or r["domain_bit_score"] > best[key]["domain_bit_score"]:
            best[key] = r
    out = sorted(best.values(), key=lambda r: -r["domain_bit_score"])
    return out


FLANKS = 7   # ORFs each side of the gene of interest — matches the discovery neighbourhoods
SCAN_NB_COLS = ["hit", "contig", "pos_index", "relationship", "is_anchor", "gene",
                "product", "locus_tag", "protein_id", "rel_start", "rel_end",
                "strand_vs_gene", "length_bp", "distance_to_anchor_bp",
                "annotation_source", "category", "vfam"]


def _select_neighbours(genes, a_s, a_e, flanks=FLANKS):
    """Contiguous neighbourhood of the gene of interest [a_s, a_e]: (upstream, downstream,
    overlapping). upstream/downstream = the `flanks` genes immediately flanking it (no
    gaps); overlapping = EVERY gene that intersects it — the overprint partner / nested
    genes (e.g. gp75's antisense RNA polymerase), which must not be dropped. [(s,e,st,meta)]."""
    up = sorted([g for g in genes if g[1] < a_s], key=lambda g: g[1])[-flanks:]
    down = sorted([g for g in genes if g[0] > a_e], key=lambda g: g[0])[:flanks]
    over = sorted([g for g in genes if g[0] <= a_e and g[1] >= a_s], key=lambda g: g[0])
    return up, down, over


def write_neighbourhoods(rows: list, contig_nt: dict, out_dir: Path, db_cache: Path,
                         cpu: int, log, annotation_gb=None, flanks=FLANKS,
                         map_tool="dfv", gene_labels=True) -> str:
    """Write an ordered neighbourhood table for each hit, anchored on the gene of
    interest, and a genome-map figure. Gene names come from **the genome's OWN
    annotation** (gene / gp number, product, locus_tag, protein_id) when a GenBank record
    is available; else genes are called de novo with **Prodigal** (the database workflow's
    caller) + optional VOGDB VFAM. The FLANKS genes each side are shown contiguously AND
    every overlapping gene (the overprint partner / nested genes — e.g. the antisense RNA
    polymerase) is included, labelled by `relationship`. Writes scan_neighbourhood.csv."""
    import csv as _csv
    try:
        import synteny_figure as SY           # anchor() + categorize()
    except Exception as e:
        log(f"  (neighbour table unavailable: {e})")
        return ""
    gb_genes = genbank_genes(annotation_gb) if annotation_gb else {}
    gb_org = genbank_organisms(annotation_gb) if annotation_gb else {}
    cache = Path(db_cache).expanduser()
    try:
        ann_ready = __import__("annotate_genes").is_ready(cache)
    except Exception:
        ann_ready = False
    all_rows = []
    for hi, r in enumerate(rows, 1):
        contig = r["contig"]
        a_s, a_e = int(r["nt_start"]), int(r["nt_end"])
        a_st = 1 if r["strand"] == "+" else -1
        if gb_genes.get(contig):
            source, genes_list = "genome annotation", gb_genes[contig]
        else:
            seq = contig_nt.get(contig, "")
            if not seq:
                continue
            try:
                import build_real_genbanks as BRG
                pg = BRG.prodigal_genes(f"{contig}__scan{hi}", seq)
            except Exception as e:
                log(f"  (neighbour gene-calling skipped for {contig}: {e})")
                continue
            source = "Prodigal"
            genes_list = [(s, e, st, {"gene": "", "product": "", "locus_tag": "",
                                      "protein_id": "", "_prot": p}) for (s, e, st, p) in pg]
        up, down, over = _select_neighbours(genes_list, a_s, a_e, flanks)
        selected = list(up) + list(down) + list(over)
        if source == "Prodigal" and ann_ready:
            prot_of = {f"g{i}": m.get("_prot", "") for i, (s, e, st, m) in enumerate(selected) if m.get("_prot")}
            if prot_of:
                try:
                    hits = __import__("annotate_genes").annotate(prot_of, cache, cpu=cpu)
                    for i, (s, e, st, m) in enumerate(selected):
                        h = hits.get(f"g{i}")
                        if h:
                            m["product"], m["vfam"] = h.get("function", ""), h.get("vfam", "")
                            m["category"] = SY.categorize(h.get("function", ""), h.get("category", ""))
                    source = "Prodigal + VOGDB VFAM"
                except Exception as e:
                    log(f"  (neighbour annotation skipped: {e})")
        # tag each gene with role + pos_index, then re-zero coords on the gene of interest
        anchor_g = {"s": a_s, "e": a_e, "st": a_st, "fam": True,
                    "meta": {"gene": "GENE OF INTEREST", "product": "family homologue"}}
        tagged = [(anchor_g, 0, "gene of interest")]
        for k, (s, e, st, m) in enumerate(reversed(up)):          # nearest upstream = -1
            tagged.append(({"s": s, "e": e, "st": st, "fam": False, "meta": m}, -(k + 1), "upstream"))
        for k, (s, e, st, m) in enumerate(down):
            tagged.append(({"s": s, "e": e, "st": st, "fam": False, "meta": m}, k + 1, "downstream"))
        for (s, e, st, m) in over:
            rel = "overlapping (antisense)" if st != a_st else "overlapping (same strand)"
            tagged.append(({"s": s, "e": e, "st": st, "fam": False, "meta": m}, 0, rel))
        if SY.anchor({"genes": [g for g, _, _ in tagged]}) is None:
            continue
        ae = anchor_g["e"]                                        # normalised gene-of-interest end
        for g, pidx, rel in tagged:
            if g["fam"] or rel.startswith("overlapping"):
                dist = 0
            elif g["e"] < 0:
                dist = g["e"]
            elif g["s"] > ae:
                dist = g["s"] - ae
            else:
                dist = 0
            m = g["meta"]
            all_rows.append({
                "hit": f"hit{hi}", "contig": contig, "pos_index": pidx, "relationship": rel,
                "is_anchor": int(bool(g["fam"])),
                "gene": m.get("gene", ""), "product": m.get("product", ""),
                "locus_tag": m.get("locus_tag", ""), "protein_id": m.get("protein_id", ""),
                "rel_start": g["s"], "rel_end": g["e"],
                "strand_vs_gene": "+" if g["st"] >= 0 else "-",
                "length_bp": g["e"] - g["s"], "distance_to_anchor_bp": dist,
                "annotation_source": "" if g["fam"] else source,
                "category": m.get("category", ""), "vfam": m.get("vfam", "")})
        # genome maps for this hit: a controllable ±flanks window AND the whole contig
        try:
            import genome_map as GM
            hl = f"hit{hi}"
            org = gb_org.get(contig, contig)
            tname = f"{org}\n{contig}" if org and org != contig else contig   # name over accession
            seq = contig_nt.get(contig, "")
            fk = {(g[0], g[1]) for g in up} | {(g[0], g[1]) for g in down}
            win = GM.build_genes((a_s, a_e, a_st), list(up) + list(down) + list(over), flank_keys=fk)
            whole = GM.build_genes((a_s, a_e, a_st), genes_list, flank_keys=fk)
            # tool-agnostic locus GenBank (open in Easyfig / Artemis / clinker / pyGenomeViz)
            gbw = GM.write_locus_genbank(win, seq, org, contig, out_dir / f"scan_genome_map_{hl}.gb")
            gba = GM.write_locus_genbank(whole, seq, org, contig, out_dir / f"scan_genome_map_{hl}_whole.gb")
            GM.draw(win, (a_s, a_e), out_dir / f"scan_genome_map_{hl}",
                    f"your gene + {flanks} genes each side", log,
                    track_name=tname, tool=map_tool, genbank=gbw, labels=gene_labels)
            GM.draw(whole, (a_s, a_e), out_dir / f"scan_genome_map_{hl}_whole",
                    f"whole genome — {len(genes_list)} genes; your gene in gold",
                    log, track_name=tname, tool=map_tool, genbank=gba, labels=gene_labels)
        except Exception as e:
            log(f"  (genome-map figure skipped: {e})")
    if not all_rows:
        return ""
    out = out_dir / "scan_neighbourhood.csv"
    with open(out, "w", newline="") as f:
        w = _csv.DictWriter(f, fieldnames=SCAN_NB_COLS)
        w.writeheader()
        w.writerows(all_rows)
    src = sorted({r["annotation_source"] for r in all_rows if r["annotation_source"]})
    n_over = sum(1 for r in all_rows if str(r["relationship"]).startswith("overlapping"))
    log(f"  neighbouring genes [{', '.join(src) or 'called'}]: {len(all_rows)} gene(s) "
        f"across {len(rows)} hit(s)" + (f", incl. {n_over} overlapping (overprint partner)" if n_over else "")
        + f" -> {out.name}")
    return str(out)


def _finish(out_dir: Path, rows: list, min_bit: float, find_interrupted: bool,
            n_contigs: int, log) -> dict:
    tsv = out_dir / "scan_hits.tsv"
    with open(tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ROW_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open(out_dir / "scan_hits_aa.faa", "w") as fa, \
         open(out_dir / "scan_hits_nt.fna", "w") as fn:
        for i, r in enumerate(rows, 1):
            hdr = (f"{r['contig']}|{r['strand']}{r['frame']}|{r['nt_start']}-{r['nt_end']}|"
                   f"{r['status']}|bit={r['domain_bit_score']}")
            if r["orf_aa"]:
                fa.write(f">{hdr}\n{r['orf_aa']}\n")
            if r["orf_nt"]:
                fn.write(f">{hdr}\n{r['orf_nt']}\n")
    clean = [r for r in rows if r["status"] == "clean"]
    interrupted = [r for r in rows if r["status"] == "interrupted"]
    if clean:
        verdict = f"GENE PRESENT — {len(clean)} clean hit(s)"
    elif interrupted:
        verdict = (f"PRESENT but INTERRUPTED — {len(interrupted)} stop-interrupted/"
                   f"overprinted copy(ies), no clean ORF")
    else:
        verdict = "GENE NOT DETECTED"
    lines = [f"Single-genome scan — {verdict}",
             f"  contigs scanned: {n_contigs};  reporting threshold: {min_bit:g} bits"
             + ("  (read-through: clean + interrupted)" if find_interrupted else "  (clean ORFs only)"),
             ""]
    for r in (rows[:10] or []):
        line = (f"  • {r['contig']} {r['strand']}{r['frame']} {r['nt_start']}-{r['nt_end']} "
                f"| {r['status']} | bit {r['domain_bit_score']} E {r['i_evalue']} "
                f"| ORF {r['orf_aa_len']}aa (M={r['has_start_M']}, stop={r['ends_at_stop']})")
        if r["status"] == "interrupted":
            line += f" | overprinting={r['overprinting_support']} (antisense_open_stops={r['antisense_open_stops']})"
        lines.append(line)
    if len(rows) > 10:
        lines.append(f"  … and {len(rows) - 10} more in scan_hits.tsv")
    report = "\n".join(lines) + "\n"
    (out_dir / "scan_report.txt").write_text(report)
    log("\n" + report)
    log(f"  -> {tsv.name}, scan_hits_aa.faa, scan_hits_nt.fna, scan_report.txt")
    return {"verdict": verdict, "n_clean": len(clean), "n_interrupted": len(interrupted),
            "n_contigs": n_contigs, "hits": len(rows), "tsv": str(tsv)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--seeds", type=Path, help="seed FASTA (protein or nucleotide CDS) to BUILD the HMM from")
    src.add_argument("--hmm", type=Path, help="an existing profile HMM to scan with")
    gsrc = ap.add_mutually_exclusive_group(required=True)
    gsrc.add_argument("--genome", type=Path,
                      help="a local genome to scan: nucleotide FASTA (.fna/.fa[.gz]) OR an "
                           "annotated GenBank (.gb/.gbk/.gbff — its gene names are used for neighbours)")
    gsrc.add_argument("--accession",
                      help="NCBI nucleotide accession(s) to fetch & scan, comma-separated "
                           "(e.g. KX098390 or NC_031062). Needs --email.")
    ap.add_argument("--email", default=None,
                    help="NCBI email (required with --accession; or set $NCBI_EMAIL)")
    ap.add_argument("--out", type=Path, default=Path("genome_scan"), help="output directory")
    ap.add_argument("--min-bit", type=float, default=25.0,
                    help="minimum domain bit score to report (default 25)")
    ap.add_argument("--find-interrupted", action="store_true",
                    help="also report stop-interrupted / overprinted copies (with the overprinting test)")
    ap.add_argument("--trans-table", type=int, default=11,
                    help="genetic code for translating a nucleotide SEED (default 11)")
    ap.add_argument("--no-neighbours", dest="neighbours", action="store_false",
                    help="skip Prodigal gene-calling of the flanking genes (on by default; "
                         "writes scan_neighbourhood.csv = the ordered neighbour table)")
    ap.add_argument("--db-cache", type=Path, default=Path("~/.cache/hmm-homologue-finder"),
                    help="cache holding the VOGDB VFAM DB used to annotate neighbour genes (optional)")
    ap.add_argument("--flanks", type=int, default=7,
                    help="number of flanking genes to report EACH side of your gene (default 7); "
                         "overlapping genes are always included")
    ap.add_argument("--no-gene-labels", dest="gene_labels", action="store_false",
                    help="draw the genome map WITHOUT gene-name labels (just coloured arrows)")
    ap.add_argument("--map-tool", choices=["dfv", "pub", "pygenomeviz", "easyfig"],
                    default="dfv",
                    help="genome-map renderer (default 'dfv' = DNA Features Viewer: clean strand "
                         "arrows, overlapping genes auto-stacked onto their own level, auto label "
                         "de-overlap, real coordinate axis, PNG+SVG+PDF). 'pub' = the built-in "
                         "matplotlib diagram (always available); 'pygenomeviz'; 'easyfig' needs "
                         "Easyfig installed + $EASYFIG_PY set. A locus GenBank is always written so "
                         "you can also open the map in Easyfig/Artemis/clinker yourself.")
    ap.add_argument("--cpu", type=int, default=4)
    ap.set_defaults(neighbours=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    log = print
    # Resolve the genome (+ any annotation): fetch by accession, or a local FASTA/GenBank.
    if args.accession:
        genome, annotation_gb = fetch_genome(args.accession, args.email, args.out, log)
    else:
        if not args.genome.exists():
            sys.exit(f"genome not found: {args.genome}")
        genome, annotation_gb = resolve_local_genome(args.genome, args.out, log)
    if args.hmm:
        if not args.hmm.exists():
            sys.exit(f"HMM not found: {args.hmm}")
        hmm = args.hmm
    else:
        if not args.seeds.exists():
            sys.exit(f"seeds not found: {args.seeds}")
        hmm = build_hmm_from_seeds(args.seeds, args.trans_table, args.out, args.cpu, log)
    s = scan(genome, hmm, args.out, args.min_bit, args.find_interrupted, args.cpu, log,
             neighbours=args.neighbours, db_cache=args.db_cache, annotation_gb=annotation_gb,
             flanks=args.flanks, map_tool=args.map_tool, gene_labels=args.gene_labels)
    # exit code 0 if the gene was detected (clean or interrupted), 1 if absent — handy in scripts
    sys.exit(0 if (s["n_clean"] or s["n_interrupted"]) else 1)


if __name__ == "__main__":
    main()
