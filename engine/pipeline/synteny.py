"""
pipeline/synteny.py — Synteny neighbourhood analysis.

Three modes, tried in order per hit:
  1. Local GenBank files  — parse .gb/.gbk directly (fast, offline)
  2. NCBI Entrez          — fetch by accession (requires internet)
  3. Script 18 wrapper    — delegate to 18_synteny_analysis.py subprocess

Accession extraction handles:
  - NCBI nucleotide  : NC_001604.1, NZ_CP012345.1, MK321214.1, MN123456.1
  - NCBI protein     : WP_012345678.1, YP_001234567.1  (→ nucleotide via elink)
  - INPHARED headers : MK321214|1_1|Phage name|...  →  MK321214
  - GPD / custom     : genome_acc_orf_N  →  genome_acc
"""
from __future__ import annotations

import io
import gzip
import re
import shutil
import sys
import time
import warnings
from pathlib import Path
from typing import Optional
from urllib.request import urlopen

import pandas as pd

# ---------------------------------------------------------------------------
# Colour palette
# ---------------------------------------------------------------------------

_FUNCTION_COLORS: dict[str, str] = {
    "replication":      "#2196F3",
    "lysis":            "#F44336",
    "structural":       "#4CAF50",
    "tail":             "#FF9800",
    "integrase":        "#9C27B0",
    "unknown":          "#9E9E9E",
    "hypothetical":     "#9E9E9E",
    "capsid":           "#4CAF50",
    "head":             "#4CAF50",
    "baseplate":        "#FF9800",
    "fiber":            "#FF9800",
    "portal":           "#4CAF50",
    "terminase":        "#2196F3",
    "polymerase":       "#2196F3",
    "recombinase":      "#9C27B0",
    "holin":            "#F44336",
    "endolysin":        "#F44336",
    "lysin":            "#F44336",
    "antirestriction":  "#00BCD4",
    "methyltransferase":"#00BCD4",
    "nuclease":         "#FF5722",
    "protease":         "#795548",
}

_KEYWORD_MAP: list[tuple[str, str]] = [
    ("replic",          "replication"),
    ("polymerase",      "replication"),
    ("terminase",       "replication"),
    ("dna bind",        "replication"),
    ("primase",         "replication"),
    ("helicase",        "replication"),
    ("holin",           "lysis"),
    ("endolysin",       "lysis"),
    ("lysin",           "lysis"),
    ("lytic",           "lysis"),
    ("spanin",          "lysis"),
    ("capsid",          "structural"),
    ("head",            "structural"),
    ("coat",            "structural"),
    ("portal",          "structural"),
    ("scaffold",        "structural"),
    ("major capsid",    "structural"),
    ("tail",            "tail"),
    ("baseplate",       "tail"),
    ("fiber",           "tail"),
    ("spike",           "tail"),
    ("tape measure",    "tail"),
    ("integrase",       "integrase"),
    ("recombinase",     "integrase"),
    ("transposase",     "integrase"),
    ("excisionase",     "integrase"),
    ("methyltransf",    "replication"),
    ("nuclease",        "nuclease"),
    ("protease",        "protease"),
    ("hypothetical",    "hypothetical"),
    ("unknown",         "unknown"),
    ("uncharacterized", "unknown"),
    ("putative",        "unknown"),
]

# Nucleotide accession patterns (covers RefSeq, GenBank, ENA/DDBJ, INPHARED, GPD, etc.)
_NT_ACCESSION_RE = re.compile(
    r'\b('
    r'(?:'
    # RefSeq
    r'NC|NZ|NG|NT|NW|NM|NR|'
    # GenBank
    r'AC|AP|AE|CP|AY|DQ|EF|EU|FJ|GQ|HM|HQ|JF|JN|JQ|JX|'
    r'KC|KF|KJ|KM|KP|KR|KT|KU|KX|KY|'
    r'MF|MG|MH|MK|MN|MT|MW|MZ|'
    r'OK|OL|OM|ON|OP|OQ|OR|OV|OW|OX|OY|OZ|PP|PS|'
    # ENA/DDBJ (common in INPHARED)
    r'AL|BK|BX|CR|CU|FN|FO|FP|FQ|FR|HE|HF|HG|LK|LM|LN|LO|LP|LR|LS|LT'
    r')_?[0-9]{5,9}'
    r'(?:\.[0-9]+)?'
    r')\b',
    re.IGNORECASE,
)

# Protein accession patterns
_PROT_ACCESSION_RE = re.compile(
    r'\b((?:WP|YP|NP|XP|AP)_[0-9]{6,9}(?:\.[0-9]+)?)\b',
    re.IGNORECASE,
)

_SIXFRAME_ORF_RE = re.compile(
    r"^(?P<genome>.+)_s(?P<strand>-?1)_f(?P<frame>[012])_o(?P<orf>\d+)$"
)


def _function_color(gene_name: str, function: str) -> str:
    combined = f"{gene_name} {function}".lower()
    for keyword, category in _KEYWORD_MAP:
        if keyword in combined:
            return _FUNCTION_COLORS.get(category, _FUNCTION_COLORS["unknown"])
    return _FUNCTION_COLORS["unknown"]


def _clean_label_value(value) -> str:
    """Return a display-safe label value, treating NaN/unknown placeholders as blank."""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    text = str(value or "").strip()
    if text.lower() in {"", ".", "nan", "none", "null", "na", "n/a", "unknown"}:
        return ""
    return text


def _short_label(value, max_len: int = 14) -> str:
    text = _clean_label_value(value)
    if len(text) <= max_len:
        return text
    if max_len <= 4:
        return text[:max_len]
    return text[: max_len - 1] + "…"


def _is_gp_or_gene_name(value: str) -> bool:
    """True for readable gene/locus names like gp63, LD35_gp81, rpoB, etc.
    False for cryptic protein accessions like UWJ04322.1, AHV82680.1."""
    import re as _re
    v = value.strip()
    if not v:
        return False
    # gp-style: gp63, PP768_gp063, LD35_gp81
    if _re.search(r'gp\d+', v, _re.IGNORECASE):
        return True
    # Named genes: rpoB, dnaA, terL, etc. (short, no dots, no version suffix)
    if len(v) <= 10 and '.' not in v and not _re.match(r'^[A-Z]{2,3}\d{5,}', v):
        return True
    # Locus_tag with underscore: PhAPEC7_56, ECBP1_0058
    if '_' in v and not _re.match(r'^[A-Z]{2,3}\d{5,}', v.split('_')[0]):
        return True
    return False


def _gene_display_label(row, hit_label: str = "hit gene", max_len: int = 14) -> str:
    """Build a readable gene label.

    Priority:
    1. Readable gene/locus name (gp-style, short gene symbols) from gene_name
    2. Informative function (not 'hypothetical protein')
    3. gp-style name from flank_gene_id or protein_id
    4. Any gene_name or protein_id as fallback
    5. 'hypothetical' as last resort
    """
    _get = (lambda k: _clean_label_value(row.get(k, ""))) if hasattr(row, "get") else (lambda k: "")

    # 1. Readable gene / locus name
    gene_name = _get("gene_name")
    if gene_name and _is_gp_or_gene_name(gene_name):
        return _short_label(gene_name, max_len=max_len)

    # 2. Informative function (skip hypothetical / unknown / Prodigal-predicted)
    func = _get("function")
    func_lower = func.lower() if func else ""
    is_informative_func = func and not any(
        kw in func_lower for kw in ("hypothetical", "unknown", "uncharacterized",
                                     "putative protein", "prodigal-predicted",
                                     "six-frame orf")
    )
    if is_informative_func:
        return _short_label(func, max_len=max_len)

    # 3. gp-style name from other ID fields
    for key in ("flank_gene_id", "protein_id"):
        val = _get(key)
        if val and _is_gp_or_gene_name(val):
            return _short_label(val, max_len=max_len)

    # 4. Any non-empty gene_name or protein_id
    if gene_name:
        return _short_label(gene_name, max_len=max_len)
    for key in ("flank_gene_id", "protein_id"):
        val = _get(key)
        if val:
            return _short_label(val, max_len=max_len)

    # 5. Fall back to function (even if hypothetical) or hit_label
    if func:
        return _short_label(func, max_len=max_len)
    return _short_label(hit_label, max_len=max_len)


# ---------------------------------------------------------------------------
# Coordinate-layout helpers  (coordinate-aware gene-map rendering)
# ---------------------------------------------------------------------------

def _coerce_int(value, default: int = 0) -> int:
    """Best-effort integer conversion for table values that may be NaN/strings."""
    try:
        if pd.isna(value):
            return default
    except Exception:
        pass
    try:
        return int(float(value))
    except Exception:
        return default


def _strand_symbol(value) -> str:
    text = str(value or "+").strip()
    if text in {"-", "-1", "reverse", "minus"}:
        return "-"
    return "+"


def _interval_note(row, group: pd.DataFrame) -> str:
    """Describe overlap/nesting relationships for hover text."""
    start = _coerce_int(row.get("start", 0))
    end   = _coerce_int(row.get("end",   0))
    if start <= 0 or end <= 0:
        return ""
    lo, hi = sorted((start, end))
    overlaps: list[str] = []
    nested = False
    row_index = getattr(row, "name", None)
    for other_index, other in group.iterrows():
        if row_index is not None and other_index == row_index:
            continue
        o_start = _coerce_int(other.get("start", 0))
        o_end   = _coerce_int(other.get("end",   0))
        if o_start <= 0 or o_end <= 0:
            continue
        o_lo, o_hi = sorted((o_start, o_end))
        if lo <= o_hi and hi >= o_lo:
            other_label  = _gene_display_label(other, max_len=18)
            other_strand = _strand_symbol(other.get("strand", "+"))
            overlaps.append(f"{other_label} ({other_strand})")
            if (lo >= o_lo and hi <= o_hi) or (o_lo >= lo and o_hi <= hi):
                nested = True
    if not overlaps:
        return ""
    prefix = "Nested/overlapping ORF" if nested else "Overlaps"
    return f"{prefix}: {', '.join(overlaps[:3])}{'...' if len(overlaps) > 3 else ''}"


def _coordinate_layout_for_genome(group: pd.DataFrame) -> "tuple[list[dict], bool]":
    """
    Return coordinate-aware gene layout records for one locus.

    Coordinates are centered on the hit gene. Overlapping genes are placed on
    separate vertical lanes so nested opposite-strand ORFs remain visible.

    Returns
    -------
    (items, ok)
        *items* — list of dicts with keys: row, start, end, x0, x1, x_mid,
        lane, lane_count, lane_offset, overlap_note.
        *ok* — False when no valid coordinates were found (fall back to
        relative-index mode).
    """
    if group is None or group.empty:
        return [], False

    rows: list[dict] = []
    for _, row in group.iterrows():
        start = _coerce_int(row.get("start", 0))
        end   = _coerce_int(row.get("end",   0))
        if start <= 0 or end <= 0:
            continue
        lo, hi = sorted((start, end))
        if hi <= lo:
            continue
        rows.append({"row": row, "start": lo, "end": hi})

    if not rows:
        return [], False

    hit_rows = [r for r in rows if float(r["row"].get("position_rel", 99) or 99) == 0]
    origin   = (
        (hit_rows[0]["start"] + hit_rows[0]["end"]) / 2.0
        if hit_rows else rows[0]["start"]
    )

    # Greedy lane assignment — prevent arrows from overlapping vertically
    lanes: list[int] = []
    for item in sorted(rows, key=lambda r: (r["start"], r["end"])):
        lane = 0
        while lane < len(lanes) and item["start"] <= lanes[lane]:
            lane += 1
        if lane == len(lanes):
            lanes.append(item["end"])
        else:
            lanes[lane] = item["end"]
        item["lane"] = lane

    lane_count = max(len(lanes), 1)
    out: list[dict] = []
    for item in rows:
        lane        = int(item.get("lane", 0))
        lane_offset = (lane - (lane_count - 1) / 2.0) * 0.20
        row         = item["row"]
        out.append({
            "row":          row,
            "start":        item["start"],
            "end":          item["end"],
            "x0":           (item["start"] - origin) / 1000.0,
            "x1":           (item["end"]   - origin) / 1000.0,
            "x_mid":        ((item["start"] + item["end"]) / 2.0 - origin) / 1000.0,
            "lane":         lane,
            "lane_count":   lane_count,
            "lane_offset":  lane_offset,
            "overlap_note": _interval_note(row, group),
        })
    return out, True


# ---------------------------------------------------------------------------
# Accession extraction
# ---------------------------------------------------------------------------

def extract_nucleotide_accession(protein_id: str, genome_id: str = "") -> str:
    """
    Try every known header format to extract a nucleotide accession.

    Priority order:
    1. Nucleotide accession regex in genome_id
    2. Nucleotide accession regex in protein_id
    3. INPHARED pipe-delimited format  (first field before |)
    4. Strip trailing _orf / _N / .N to recover genome id
    5. Protein accession (caller must convert via elink separately)
    """
    for candidate in [genome_id, protein_id]:
        if not candidate:
            continue
        # Direct nucleotide accession
        m = _NT_ACCESSION_RE.search(candidate)
        if m:
            return m.group(1)

    # INPHARED / GPD pipe format: MK321214|1_1|Phage name|...
    for candidate in [protein_id, genome_id]:
        if "|" in (candidate or ""):
            first_field = candidate.split("|")[0].strip()
            m = _NT_ACCESSION_RE.match(first_field)
            if m:
                return m.group(1)

    # Strip 6-frame translation suffixes: acc_s1_f2_o824, acc_s-1_f0_o1006
    for candidate in [genome_id, protein_id]:
        if not candidate:
            continue
        cleaned = re.sub(r'_s-?[12]_f[012]_o\d+$', '', candidate).strip()
        if cleaned != candidate:
            m = _NT_ACCESSION_RE.match(cleaned)
            if m:
                return m.group(1)

    # Strip common ORF suffixes: acc_orf_1, acc_1, acc.prot1
    for candidate in [genome_id, protein_id]:
        if not candidate:
            continue
        cleaned = re.sub(r'[_\.\-](?:orf|prot|CDS|gene|protein|ORF)?[_\.\-]?\d+$',
                         '', candidate, flags=re.IGNORECASE).strip()
        m = _NT_ACCESSION_RE.match(cleaned)
        if m:
            return m.group(1)

    return ""


def extract_protein_accession(protein_id: str) -> str:
    """Return a NCBI protein accession (WP_, YP_, NP_, XP_) if present."""
    m = _PROT_ACCESSION_RE.search(protein_id or "")
    return m.group(1) if m else ""


# ---------------------------------------------------------------------------
# Six-frame FASTA context rescue
# ---------------------------------------------------------------------------

def _record_key(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or ""))[:120]


def _candidate_accessions(protein_id: str, genome_id: str = "") -> list[str]:
    """Return accession candidates ordered from most to least specific."""
    vals: list[str] = []
    for val in (genome_id, extract_nucleotide_accession(protein_id, genome_id), protein_id):
        val = str(val or "").strip()
        if val and val.lower() != "nan" and val not in vals:
            vals.append(val)
            base = val.split(".")[0]
            if base and base not in vals:
                vals.append(base)
    parsed = _parse_sixframe_orf_id(protein_id)
    if parsed:
        genome = parsed["genome"]
        if genome not in vals:
            vals.append(genome)
        base = genome.split(".")[0]
        if base and base not in vals:
            vals.append(base)
    return vals


def _parse_sixframe_orf_id(protein_id: str) -> Optional[dict]:
    """Parse IDs emitted by the app's 6-frame translator: genome_s1_f0_o123."""
    m = _SIXFRAME_ORF_RE.match(str(protein_id or ""))
    if not m:
        return None
    return {
        "genome": m.group("genome"),
        "strand": int(m.group("strand")),
        "frame": int(m.group("frame")),
        "orf_index": int(m.group("orf")),
    }


_PRODIGAL_DESC_RE = re.compile(
    r"^\s*#\s*(?P<start>\d+)\s*#\s*(?P<end>\d+)\s*#\s*(?P<strand>-?1)\s*(?:#|$)"
)


def _parse_prodigal_hit(protein_id: str, description: str = "") -> Optional[dict]:
    """Parse Prodigal FASTA headers preserved by HMMER tblout descriptions."""
    desc_match = _PRODIGAL_DESC_RE.match(str(description or ""))
    if not desc_match:
        return None

    pid = str(protein_id or "").strip()
    genome = re.sub(r"_\d+$", "", pid).strip()
    if not genome:
        return None

    start = int(desc_match.group("start"))
    end = int(desc_match.group("end"))
    strand = "+" if int(desc_match.group("strand")) >= 0 else "-"

    return {
        "genome": genome,
        "start": min(start, end),
        "end": max(start, end),
        "strand": strand,
    }


def _sixframe_segment_coords(seq, strand: int, frame: int, orf_index: int) -> Optional[tuple[int, int, str, int]]:
    """
    Reconstruct nucleotide coordinates for one six-frame ORF segment.

    The search step labels ORFs by their index in ``translate().split("*")``.
    Replaying that exact indexing lets synteny recover the original genomic
    position even when the streamed database was not saved locally.
    """
    nuc = seq if strand == 1 else seq.reverse_complement()
    usable = ((len(nuc) - frame) // 3) * 3
    trans = str(nuc[frame:frame + usable].translate())
    aa_offset = 0
    seq_len = len(seq)
    for idx, aa in enumerate(trans.split("*")):
        aa_len = len(aa)
        if idx == orf_index:
            nt_start = frame + aa_offset * 3
            nt_end = min(frame + (aa_offset + aa_len) * 3, seq_len)
            if strand == 1:
                start = nt_start + 1
                end = nt_end
                out_strand = "+"
            else:
                start = seq_len - nt_end + 1
                end = seq_len - nt_start
                out_strand = "-"
            return (min(start, end), max(start, end), out_strand, aa_len)
        aa_offset += aa_len + 1
    return None


def _sixframe_orf_genes(record, min_aa: int = 30) -> list[dict]:
    """Fallback gene calls from all six-frame ORFs when Prodigal is unavailable."""
    genes: list[dict] = []
    seq = record.seq
    for strand in (1, -1):
        nuc = seq if strand == 1 else seq.reverse_complement()
        for frame in range(3):
            usable = ((len(nuc) - frame) // 3) * 3
            trans = str(nuc[frame:frame + usable].translate())
            aa_offset = 0
            for idx, aa in enumerate(trans.split("*")):
                if len(aa) >= min_aa:
                    coords = _sixframe_segment_coords(seq, strand, frame, idx)
                    if coords:
                        start, end, strand_txt, aa_len = coords
                        gid = f"{record.id}_s{strand}_f{frame}_o{idx}"
                        genes.append({
                            "gene_name": gid,
                            "protein_id": gid,
                            "start": start,
                            "end": end,
                            "strand": strand_txt,
                            "function": f"six-frame ORF ({aa_len} aa)",
                            "color": _function_color(gid, "hypothetical protein"),
                            "protein_seq": aa,
                        })
                aa_offset += len(aa) + 1
    genes.sort(key=lambda g: (g["start"], g["end"]))
    return genes


def _prodigal_genes_for_record(record, cache_dir: Path, accession: str) -> list[dict]:
    """Predict CDS with Prodigal and parse its GFF output."""
    try:
        from pipeline.utils import find_tool
        prodigal = find_tool("prodigal")
    except Exception:
        prodigal = shutil.which("prodigal")
    if not prodigal:
        return []
    try:
        from Bio import SeqIO
        import subprocess
    except Exception:
        return []

    cache_dir.mkdir(parents=True, exist_ok=True)
    stem = _record_key(accession or record.id)
    fna = cache_dir / f"{stem}.fna"
    gff = cache_dir / f"{stem}.prodigal.gff"
    if not fna.exists():
        SeqIO.write(record, str(fna), "fasta")
    if not gff.exists():
        cmd = [prodigal, "-i", str(fna), "-o", str(gff), "-f", "gff", "-p", "meta", "-q"]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            return []

    genes: list[dict] = []
    try:
        for line in gff.read_text().splitlines():
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6] if parts[6] in ("+", "-") else "+"
            attrs = {}
            for item in parts[8].split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    attrs[k] = v
            gid = attrs.get("ID") or attrs.get("locus_tag") or f"{record.id}_{start}_{end}"
            # Translate CDS from the nucleotide record so clinker can draw links
            try:
                sub = record.seq[start - 1:end]
                if strand == "-":
                    sub = sub.reverse_complement()
                trim = len(sub) - (len(sub) % 3)
                protein_seq = str(sub[:trim].translate(to_stop=True))
            except Exception:
                protein_seq = ""
            genes.append({
                "gene_name": gid,
                "protein_id": gid,
                "start": start,
                "end": end,
                "strand": strand,
                "function": "Prodigal-predicted CDS",
                "color": _function_color(gid, "hypothetical protein"),
                "protein_seq": protein_seq,
            })
    except Exception:
        return []
    genes.sort(key=lambda g: (g["start"], g["end"]))
    return genes


def _window_from_gene_calls(
    genes: list[dict],
    hit_start: int,
    hit_end: int,
    flanks: int,
    hit_protein_id: str,
    hit_name: str = "hit gene",
) -> list[dict]:
    """Return a fixed flanking window around the hit coordinates."""
    if not genes or hit_start <= 0 or hit_end <= 0:
        return []

    hit_mid = (hit_start + hit_end) / 2.0
    overlap_idx = [
        i for i, g in enumerate(genes)
        if int(g.get("start", 0)) <= hit_end and int(g.get("end", 0)) >= hit_start
    ]
    if overlap_idx:
        hit_idx = min(
            overlap_idx,
            key=lambda i: abs(((genes[i]["start"] + genes[i]["end"]) / 2.0) - hit_mid),
        )
        genes = [dict(g) for g in genes]
        genes[hit_idx]["protein_id"] = genes[hit_idx].get("protein_id") or hit_protein_id
        genes[hit_idx]["gene_name"] = genes[hit_idx].get("gene_name") or hit_name
        genes[hit_idx]["function"] = genes[hit_idx].get("function") or "hit gene"
        genes[hit_idx]["color"] = "#E91E63"
    else:
        hit_gene = {
            "gene_name": hit_name,
            "protein_id": hit_protein_id,
            "start": hit_start,
            "end": hit_end,
            "strand": "+",
            "function": "six-frame HMM hit",
            "color": "#E91E63",
        }
        insert_at = 0
        while insert_at < len(genes) and int(genes[insert_at].get("start", 0)) < hit_start:
            insert_at += 1
        genes = [dict(g) for g in genes]
        genes.insert(insert_at, hit_gene)
        hit_idx = insert_at

    lo = max(0, hit_idx - flanks)
    hi = min(len(genes) - 1, hit_idx + flanks)
    result: list[dict] = []
    for i in range(lo, hi + 1):
        gene = dict(genes[i])
        gene["position_rel"] = i - hit_idx
        if i == hit_idx:
            gene["color"] = "#E91E63"
            gene["function"] = gene.get("function") or "hit gene"
        result.append(gene)
    return result


def _read_cached_record(path: Path):
    try:
        from Bio import SeqIO
        return SeqIO.read(str(path), "fasta")
    except Exception:
        return None


def _prepare_sequence_context_cache(
    hits_df: pd.DataFrame,
    max_genomes: int,
    cache_dir: Optional[Path],
    log_callback=None,
) -> dict[str, Path]:
    """
    Cache source nucleotide records needed for streamed six-frame synteny.

    The search stage stores ``source_url`` for streamed nucleotide databases.
    This function streams each source FASTA once, extracts only the hit
    genomes, and writes tiny per-genome FASTA files used by the synteny rescue.
    """
    if hits_df is None or hits_df.empty or "source_url" not in hits_df.columns:
        return {}
    try:
        from Bio import SeqIO
    except Exception:
        return {}

    cache_root = Path(cache_dir or Path.cwd() / "synteny_context_cache")
    cache_root.mkdir(parents=True, exist_ok=True)

    needed_by_url: dict[str, set[str]] = {}
    seen: set[str] = set()
    for _, row in hits_df.iterrows():
        pid = str(row.get("protein_id", row.get("target_name", "")) or "")
        parsed = _parse_sixframe_orf_id(pid)
        acc = str(
            row.get("source_contig", "")
            or row.get("genome_id", "")
            or (parsed["genome"] if parsed else "")
        )
        try:
            has_coords = int(row.get("seq_from", 0) or 0) > 0 and int(row.get("seq_to", 0) or 0) > 0
        except Exception:
            has_coords = False
        if not parsed and not (acc and has_coords):
            continue
        if not acc or acc.lower() == "nan" or acc in seen:
            continue
        source_url = str(row.get("source_url", "") or row.get("db_path", "") or "")
        if not source_url:
            continue
        seen.add(acc)
        parsed_genome = parsed["genome"] if parsed else acc
        source_candidates = {
            acc,
            acc.split(".")[0],
            parsed_genome,
            parsed_genome.split(".")[0],
        }
        needed_by_url.setdefault(source_url, set()).update(c for c in source_candidates if c)
        if len(seen) >= max_genomes:
            break

    found: dict[str, Path] = {}

    def _log(msg: str):
        if log_callback:
            log_callback(msg)

    for url, candidates in needed_by_url.items():
        pending = set(candidates)
        for cand in list(pending):
            p = cache_root / f"{_record_key(cand)}.fna"
            if p.exists():
                found[cand] = p
                pending.discard(cand)
        if not pending:
            continue

        label = url.split("/")[-1] or url
        _log(f"  📥  Fetching source records for synteny context from {label} …")
        try:
            if re.match(r"https?://", url):
                raw = urlopen(url, timeout=60)
                handle = (
                    io.TextIOWrapper(gzip.GzipFile(fileobj=raw))
                    if url.endswith(".gz")
                    else io.TextIOWrapper(raw)
                )
            else:
                path = Path(url)
                handle = gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()
            with handle:
                for rec in SeqIO.parse(handle, "fasta"):
                    rec_ids = {rec.id, rec.id.split()[0], rec.id.split(".")[0]}
                    match = next((cand for cand in pending if cand in rec_ids or cand.split(".")[0] in rec_ids), None)
                    if not match:
                        continue
                    out = cache_root / f"{_record_key(match)}.fna"
                    SeqIO.write(rec, str(out), "fasta")
                    found[match] = out
                    for alias in list(pending):
                        if alias == match or alias.split(".")[0] == match.split(".")[0]:
                            found[alias] = out
                            pending.discard(alias)
                    _log(f"    ✅  cached {rec.id}")
                    if not pending:
                        break
        except Exception as exc:
            _log(f"    ⚠️  Could not fetch synteny context from {label}: {exc}")

    return found


def fetch_neighborhood_sixframe_context(
    hit_row,
    flanks: int,
    sequence_cache: dict[str, Path],
    gene_cache_dir: Optional[Path] = None,
) -> list[dict]:
    """Build a neighbourhood from a streamed nucleotide FASTA hit."""
    protein_id = str(hit_row.get("protein_id", hit_row.get("target_name", "")) or "")
    parsed = _parse_sixframe_orf_id(protein_id)
    genome_id = str(
        hit_row.get("source_contig", "")
        or hit_row.get("genome_id", "")
        or (parsed["genome"] if parsed else "")
    )
    try:
        row_start = int(hit_row.get("seq_from", 0) or 0)
        row_end = int(hit_row.get("seq_to", 0) or 0)
    except Exception:
        row_start = row_end = 0
    row_strand = str(hit_row.get("strand", "") or "")
    if not parsed and (not genome_id or row_start <= 0 or row_end <= 0):
        return []
    record_path = None
    for cand in _candidate_accessions(protein_id, genome_id):
        if cand in sequence_cache:
            record_path = sequence_cache[cand]
            break
    if record_path is None:
        return []
    record = _read_cached_record(record_path)
    if record is None:
        return []
    coords = None
    if parsed:
        coords = _sixframe_segment_coords(
            record.seq,
            parsed["strand"],
            parsed["frame"],
            parsed["orf_index"],
        )
    if coords:
        hit_start, hit_end, hit_strand, _aa_len = coords
    elif row_start > 0 and row_end > 0:
        hit_start, hit_end = sorted((row_start, row_end))
        hit_strand = row_strand if row_strand in {"+", "-"} else "+"
    else:
        return []
    prodigal_genes = _prodigal_genes_for_record(
        record,
        Path(gene_cache_dir or record_path.parent),
        genome_id or (parsed["genome"] if parsed else ""),
    )
    prodigal_window = _window_from_gene_calls(
        prodigal_genes,
        hit_start,
        hit_end,
        flanks,
        protein_id,
        str(hit_row.get("hit_name", "hit gene") or "hit gene"),
    )
    source = "streamed_fasta_prodigal"
    window = prodigal_window

    # Prodigal can occasionally merge artificial or very atypical contigs into
    # one giant CDS. In that case the app still needs a usable neighbourhood,
    # so prefer the six-frame reconstruction when it gives a fuller window.
    sixframe_window = []
    if len(prodigal_window) < (2 * flanks + 1):
        sixframe_window = _window_from_gene_calls(
            _sixframe_orf_genes(record),
            hit_start,
            hit_end,
            flanks,
            protein_id,
            str(hit_row.get("hit_name", "hit gene") or "hit gene"),
        )
    if len(sixframe_window) > len(prodigal_window):
        window = sixframe_window
        source = "streamed_fasta_sixframe"

    for g in window:
        g["accession"] = genome_id or (parsed["genome"] if parsed else "")
        g["source"] = source
        if g.get("protein_id") == protein_id or g.get("position_rel") == 0:
            g["strand"] = hit_strand
    return window


# ---------------------------------------------------------------------------
# Local GenBank parsing
# ---------------------------------------------------------------------------

def _find_local_genbank(accession: str, search_dirs: list[Path]) -> Optional[Path]:
    """Locate a GenBank file (.gb / .gbk / .gbff) for *accession* in search_dirs."""
    for d in search_dirs:
        if not d or not Path(d).is_dir():
            continue
        for suffix in (".gb", ".gbk", ".gbff", ".gb.gz", ".gbk.gz"):
            # exact match
            p = Path(d) / f"{accession}{suffix}"
            if p.exists():
                return p
            # glob with version stripped
            acc_base = accession.split(".")[0]
            for hit in Path(d).glob(f"{acc_base}*{suffix}"):
                return hit
    return None


def _parse_neighborhood_from_record(
    record, seq_from: int, seq_to: int, flanks: int,
    hit_strand: str = "",
    hit_protein_id: str = "",
    hit_name: str = "hit gene",
) -> list[dict]:
    """Extract flanking CDS from a BioPython SeqRecord.

    When *hit_strand* is provided, the hit gene is matched to the nearest
    CDS on the **same strand**.  If no same-strand CDS overlaps the hit
    coordinates, a synthetic "novel ORF" entry is inserted so the hit is
    clearly distinguished from any overlapping gene on the opposite strand
    (e.g., a small ORF encoded inside a large vRNAP on the other strand).
    """
    cds_list: list[dict] = []
    for feat in record.features:
        if feat.type != "CDS":
            continue
        f_start = int(feat.location.start) + 1
        f_end   = int(feat.location.end)
        qualifiers = feat.qualifiers
        gene_name  = qualifiers.get("gene",  qualifiers.get("locus_tag", [""]))[0]
        protein_id = qualifiers.get("protein_id", [""])[0]
        function   = qualifiers.get("product", ["hypothetical protein"])[0]
        strand     = "+" if feat.location.strand == 1 else "-"
        color      = _function_color(gene_name, function)
        protein_seq = qualifiers.get("translation", [""])[0]
        cds_list.append({
            "gene_name":  gene_name,
            "protein_id": protein_id,
            "start":      f_start,
            "end":        f_end,
            "strand":     strand,
            "function":   function,
            "color":      color,
            "protein_seq": protein_seq,
        })

    if not cds_list:
        return []

    cds_list.sort(key=lambda x: x["start"])
    hit_mid = (seq_from + seq_to) / 2.0

    # ── Strand-aware hit matching ──────────────────────────────────────
    # Prefer a CDS on the same strand that overlaps the hit coordinates.
    # Fall back to closest-by-midpoint on any strand if none found.
    hit_idx = None
    if hit_strand:
        # Same-strand overlapping CDS
        same_strand_overlaps = [
            i for i, g in enumerate(cds_list)
            if g["strand"] == hit_strand
            and g["start"] <= seq_to and g["end"] >= seq_from
        ]
        if same_strand_overlaps:
            hit_idx = min(same_strand_overlaps,
                          key=lambda i: abs((cds_list[i]["start"] + cds_list[i]["end"]) / 2.0 - hit_mid))
        else:
            # No same-strand CDS overlaps the hit → novel ORF inside
            # an opposite-strand gene. Insert a synthetic entry so the
            # hit is NOT confused with the overlapping gene.
            novel_hit = {
                "gene_name":  hit_name,
                "protein_id": hit_protein_id or hit_name,
                "start":      seq_from,
                "end":        seq_to,
                "strand":     hit_strand,
                "function":   "novel ORF (opposite strand from annotated gene)",
                "color":      "#E91E63",
                "is_novel":   True,
            }
            # Find insertion point by coordinate
            insert_at = 0
            while insert_at < len(cds_list) and cds_list[insert_at]["start"] < seq_from:
                insert_at += 1
            cds_list.insert(insert_at, novel_hit)
            hit_idx = insert_at

    if hit_idx is None:
        # Original behaviour: closest CDS by midpoint
        hit_idx = min(range(len(cds_list)),
                      key=lambda i: abs((cds_list[i]["start"] + cds_list[i]["end"]) / 2.0 - hit_mid))

    lo       = max(0, hit_idx - flanks)
    hi       = min(len(cds_list) - 1, hit_idx + flanks)
    window   = cds_list[lo: hi + 1]

    result = []
    for i, gene in enumerate(window):
        gene = dict(gene)
        gene["position_rel"] = i - (hit_idx - lo)
        result.append(gene)
    return result


def fetch_neighborhood_local(
    accession: str,
    seq_from: int,
    seq_to: int,
    flanks: int,
    search_dirs: list[Path],
    hit_strand: str = "",
    hit_protein_id: str = "",
    hit_name: str = "hit gene",
) -> list[dict]:
    """Get flanking genes from a local GenBank file — no internet required."""
    try:
        from Bio import SeqIO
    except ImportError:
        return []

    gb_path = _find_local_genbank(accession, search_dirs)
    if gb_path is None:
        return []

    try:
        import gzip
        opener = gzip.open if str(gb_path).endswith(".gz") else open
        with opener(gb_path, "rt") as fh:
            record = SeqIO.read(fh, "genbank")
        genes = _parse_neighborhood_from_record(
            record, seq_from, seq_to, flanks,
            hit_strand=hit_strand,
            hit_protein_id=hit_protein_id,
            hit_name=hit_name,
        )
        for g in genes:
            g["accession"] = accession
            g["source"]    = "local_genbank"
        return genes
    except Exception as exc:
        print(f"WARNING [local_genbank] {gb_path}: {exc}", file=sys.stderr)
        return []


# ---------------------------------------------------------------------------
# NCBI Entrez fetch
# ---------------------------------------------------------------------------

def fetch_neighborhood_entrez(
    accession: str,
    seq_from: int,
    seq_to: int,
    strand: str,
    flanks: int,
    email: str,
    hit_protein_id: str = "",
    hit_name: str = "hit gene",
) -> list[dict]:
    """Fetch flanking genes from NCBI Entrez (requires internet + valid accession)."""
    try:
        from Bio import Entrez, SeqIO
    except ImportError:
        return []

    Entrez.email = email
    time.sleep(0.35)   # NCBI rate limit

    # If we only have a protein accession, convert to nucleotide via elink
    if _PROT_ACCESSION_RE.match(accession):
        try:
            link_handle = Entrez.elink(dbfrom="protein", db="nuccore", id=accession)
            link_record = Entrez.read(link_handle)
            link_handle.close()
            links = link_record[0].get("LinkSetDb", [])
            if links:
                accession = links[0]["Link"][0]["Id"]
                time.sleep(0.35)
        except Exception:
            pass

    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
    except Exception as exc:
        print(f"WARNING [entrez] {accession}: {exc}", file=sys.stderr)
        return []

    genes = _parse_neighborhood_from_record(
        record, seq_from, seq_to, flanks,
        hit_strand=strand,
        hit_protein_id=hit_protein_id,
        hit_name=hit_name,
    )
    for g in genes:
        g["accession"] = accession
        g["source"]    = "ncbi_entrez"
    return genes


# ---------------------------------------------------------------------------
# Build synteny table  (multi-mode)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Missing-context classifier
# ---------------------------------------------------------------------------

def classify_missing_context(protein_id: str, genome_id: str) -> tuple:
    """
    Categorise *why* a hit cannot be placed into a synteny neighbourhood.

    Returns
    -------
    (status, reason, rescue_hint) — all strings.

    status values
    -------------
    ``"no_id"``             Both protein_id and genome_id are empty / NaN.
    ``"protein_only_ncbi"`` NCBI WP_/YP_/NP_ protein accession only — Entrez
                             elink can recover the nucleotide locus.
    ``"no_coords"``         Accession found but seq_from/seq_to are both 0 —
                             the tblout had no coordinate data.
    ``"custom_db"``         Header doesn't match any known accession pattern.
                             Requires local GenBank files.
    ``"non_ncbi_db"``       INPHARED / GPD / MGV header detected — these
                             databases are not on NCBI Entrez.
    """
    pid = str(protein_id or "").strip()
    gid = str(genome_id  or "").strip()

    if not pid and not gid:
        return (
            "no_id",
            "Both protein_id and genome_id are empty",
            "Check that the hits table was parsed correctly from the tblout file.",
        )

    # INPHARED / GPD-style headers: contain pipe-delimited fields
    for cand in (pid, gid):
        if "|" in cand:
            fields = cand.split("|")
            if _NT_ACCESSION_RE.match(fields[0].strip()):
                return (
                    "non_ncbi_db",
                    f"INPHARED/GPD-style header: {cand[:60]}",
                    "Download the matching GenBank files from the database source and "
                    "provide the folder path in 'Local GenBank folder'.",
                )

    # NCBI protein accession (WP_, YP_, NP_, XP_)
    for cand in (pid, gid):
        if _PROT_ACCESSION_RE.search(cand):
            return (
                "protein_only_ncbi",
                f"NCBI protein accession only: {_PROT_ACCESSION_RE.search(cand).group(1)}",  # type: ignore[union-attr]
                "Provide your NCBI email — the app will use Entrez elink to convert "
                "protein → nucleotide and fetch the neighbourhood automatically.",
            )

    # Nucleotide accession exists — must be coords missing (handled by caller)
    for cand in (pid, gid):
        if _NT_ACCESSION_RE.search(cand):
            return (
                "no_coords",
                f"Accession found ({_NT_ACCESSION_RE.search(cand).group(1)}) "  # type: ignore[union-attr]
                "but seq_from / seq_to are both 0",
                "Coordinates are missing from the hmmsearch tblout. Re-run the search "
                "making sure the correct parser is used, or use NCBI Entrez mode "
                "which fetches coordinates from the protein record.",
            )

    return (
        "custom_db",
        f"No recognisable accession in: pid='{pid[:40]}' | gid='{gid[:40]}'",
        "This hit is from a custom database. Provide the matching GenBank files "
        "in 'Local GenBank folder', with filenames that include part of the "
        "accession string (e.g. genome_id.gbk).",
    )


# ---------------------------------------------------------------------------
# Protein-coordinate rescue  (WP_ / YP_ with seq_from == 0)
# ---------------------------------------------------------------------------

def fetch_protein_coordinates(prot_accession: str, email: str) -> tuple:
    """
    Fetch the nucleotide locus and coordinates for an NCBI protein accession.

    Uses ``Entrez.efetch(db="protein", rettype="gp")`` to read the
    ``/coded_by`` qualifier, which gives the nucleotide accession and
    coordinate range directly (e.g. ``"NC_001604.1:12345..67890"``).

    Returns
    -------
    (nt_accession, seq_from, seq_to, strand)  on success,
    or ``("", 0, 0, "+")`` on failure.
    """
    try:
        from Bio import Entrez, SeqIO
    except ImportError:
        return ("", 0, 0, "+")

    Entrez.email = email
    time.sleep(0.35)

    try:
        handle = Entrez.efetch(db="protein", id=prot_accession,
                               rettype="gp", retmode="text")
        rec = SeqIO.read(handle, "genbank")
        handle.close()
    except Exception:
        return ("", 0, 0, "+")

    # The CDS feature on a protein record carries /coded_by
    for feat in rec.features:
        if feat.type != "CDS":
            continue
        coded_by = feat.qualifiers.get("coded_by", [""])[0]
        if not coded_by:
            continue
        # Format: "NC_001604.1:12345..67890"  or "complement(NC_001604.1:12345..67890)"
        strand = "-" if coded_by.startswith("complement") else "+"
        m = re.search(
            r'([A-Z]{1,2}_?[0-9]{5,9}(?:\.[0-9]+)?):(\d+)\.\.(\d+)',
            coded_by,
        )
        if m:
            nt_acc  = m.group(1)
            sf      = int(m.group(2))
            st      = int(m.group(3))
            return (nt_acc, sf, st, strand)

    return ("", 0, 0, "+")


# ---------------------------------------------------------------------------
# Main synteny table builder  (returns neighbourhood + placement report)
# ---------------------------------------------------------------------------

def build_synteny_table(
    hits_df: pd.DataFrame,
    email: str = "researcher@example.com",
    flanks: int = 5,
    max_genomes: int = 30,
    local_genbank_dirs: Optional[list] = None,
    sequence_cache_dir: Optional[Path] = None,
    scripts_dir: Optional[Path] = None,
    log_callback=None,
) -> tuple:
    """
    Build a long-format synteny neighbourhood table from a hits DataFrame.

    Tries in order for each hit:
      1. Local GenBank files (fast, offline)
      2. NCBI Entrez fetch  (requires internet + valid accession)
      3. Protein-coordinate rescue  (WP_/YP_ without seq coords → efetch)

    Parameters
    ----------
    hits_df : pd.DataFrame
        Must have ``protein_id``, ``genome_id``.
        Optional but helpful: ``seq_from``, ``seq_to``, ``strand``, ``bit_score``.
    email : str
        NCBI Entrez email.  Required for modes 2 and 3.
    flanks : int
        Genes each side of the hit gene.
    max_genomes : int
        Maximum number of unique loci to query.
    local_genbank_dirs : list of Path-like, optional
        Directories to search for local .gb / .gbk / .gbff files.
    sequence_cache_dir : Path-like, optional
        Cache directory for extracting source nucleotide records from streamed
        six-frame FASTA database hits.
    log_callback : callable, optional
        Receives progress strings during execution.

    Returns
    -------
    (syn_df, placement_df)

    ``syn_df`` : pd.DataFrame
        Long-format table with one row per flanking gene.
        Columns: hit_protein_id, accession, flank_gene_id, position_rel,
        gene_name, function, strand, start, end, color, source.

    ``placement_df`` : pd.DataFrame
        One row per input hit showing placement outcome.
        Columns: protein_id, genome_id, accession, status, reason,
        rescue_hint, placed.
    """
    _EMPTY_SYN = pd.DataFrame(columns=[
        "hit_protein_id", "accession", "flank_gene_id", "position_rel",
        "gene_name", "function", "strand", "start", "end", "color", "source",
        "protein_seq",
    ])
    _EMPTY_REPORT = pd.DataFrame(columns=[
        "protein_id", "genome_id", "accession",
        "status", "reason", "rescue_hint", "placed",
    ])

    if hits_df is None or hits_df.empty:
        return _EMPTY_SYN, _EMPTY_REPORT

    gb_dirs: list[Path] = [Path(d) for d in (local_genbank_dirs or []) if d]

    use_entrez = bool(email and email not in ("", "researcher@example.com"))

    # Sort best hits first
    if "bit_score" in hits_df.columns:
        df = hits_df.sort_values("bit_score", ascending=False).copy()
    else:
        df = hits_df.copy()

    seen_accessions: set[str] = set()
    rows: list[dict]          = []
    report_rows: list[dict]   = []

    def _log(msg: str):
        if log_callback:
            log_callback(msg)

    sequence_cache = _prepare_sequence_context_cache(
        df,
        max_genomes=max_genomes,
        cache_dir=sequence_cache_dir,
        log_callback=log_callback,
    )

    for _, hit_row in df.iterrows():
        protein_id = str(hit_row.get("protein_id", hit_row.get("target_name", "")) or "")
        genome_id  = str(hit_row.get("genome_id",  "") or "")
        # _coerce_int is NaN-safe: `int(x or 0)` raises on NaN because NaN is
        # truthy, so a hits table without coordinate columns would crash here.
        seq_from   = _coerce_int(hit_row.get("seq_from", 0))
        seq_to     = _coerce_int(hit_row.get("seq_to",   0))
        strand_raw = hit_row.get("strand", "+")
        strand_val = str(strand_raw) if strand_raw is not None and not (
            isinstance(strand_raw, float) and pd.isna(strand_raw)
        ) else "+"

        prodigal_hit = _parse_prodigal_hit(
            protein_id,
            str(hit_row.get("description", "") or ""),
        )
        if prodigal_hit:
            if not genome_id or genome_id == protein_id:
                genome_id = prodigal_hit["genome"]
            if seq_from <= 0 or seq_to <= 0:
                seq_from = prodigal_hit["start"]
                seq_to = prodigal_hit["end"]
                strand_val = prodigal_hit["strand"]

        # ── Already at capacity? ────────────────────────────────────────────
        if len(seen_accessions) >= max_genomes:
            report_rows.append({
                "protein_id": protein_id, "genome_id": genome_id,
                "accession": "", "status": "capped",
                "reason": f"Reached max_genomes limit ({max_genomes})",
                "rescue_hint": "Increase 'Max genomes' slider to include more loci.",
                "placed": False,
            })
            continue

        # ── Extract nucleotide accession ────────────────────────────────────
        accession = extract_nucleotide_accession(protein_id, genome_id)
        if not accession:
            accession = genome_id or protein_id   # last-resort: use raw id for local lookup

        if not accession or accession in ("nan", ""):
            status, reason, hint = classify_missing_context(protein_id, genome_id)
            report_rows.append({
                "protein_id": protein_id, "genome_id": genome_id,
                "accession": "", "status": status,
                "reason": reason, "rescue_hint": hint, "placed": False,
            })
            _log(f"  ⚠️  No accession for: {protein_id[:40]} — {status}")
            continue

        # ── Duplicate accession ─────────────────────────────────────────────
        if accession in seen_accessions:
            report_rows.append({
                "protein_id": protein_id, "genome_id": genome_id,
                "accession": accession, "status": "duplicate",
                "reason": "Same genome accession already processed",
                "rescue_hint": "", "placed": True,  # neighbourhood already in table
            })
            continue
        seen_accessions.add(accession)

        # ── Mode 1: Local GenBank ───────────────────────────────────────────
        genes = fetch_neighborhood_local(
            accession, seq_from, seq_to, flanks, gb_dirs,
            hit_strand=strand_val, hit_protein_id=protein_id,
            hit_name=str(hit_row.get("hit_name", "hit gene") or "hit gene"),
        )
        source_mode = "local_genbank" if genes else ""

        # ── Mode 2: NCBI Entrez ─────────────────────────────────────────────
        if not genes and use_entrez:
            nt_acc = accession
            if not _NT_ACCESSION_RE.match(accession):
                prot_acc = extract_protein_accession(protein_id)
                if prot_acc:
                    nt_acc = prot_acc
            genes = fetch_neighborhood_entrez(
                nt_acc, seq_from, seq_to, strand_val, flanks, email,
                hit_protein_id=protein_id,
                hit_name=str(hit_row.get("hit_name", "hit gene") or "hit gene"),
            )
            if genes:
                source_mode = "ncbi_entrez"

        # ── Mode 3: Protein-coordinate rescue (WP_/YP_ + seq_from==0) ──────
        if not genes and use_entrez and seq_from == 0 and seq_to == 0:
            prot_acc = extract_protein_accession(protein_id)
            if prot_acc:
                _log(f"  🔍  Rescuing coordinates for {prot_acc} via Entrez protein efetch …")
                nt_acc2, sf2, st2, strand2 = fetch_protein_coordinates(prot_acc, email)
                if nt_acc2 and sf2 > 0:
                    genes = fetch_neighborhood_entrez(
                        nt_acc2, sf2, st2, strand2, flanks, email,
                    )
                    if genes:
                        accession   = nt_acc2
                        seq_from    = sf2
                        seq_to      = st2
                        strand_val  = strand2
                        source_mode = "protein_coord_rescue"
                        _log(f"    ✅  Rescued: {nt_acc2}:{sf2}..{st2} ({strand2})")

        # ── Mode 4: Streamed nucleotide FASTA context rescue ───────────────
        if not genes:
            genes = fetch_neighborhood_sixframe_context(
                hit_row,
                flanks=flanks,
                sequence_cache=sequence_cache,
                gene_cache_dir=Path(sequence_cache_dir) if sequence_cache_dir else None,
            )
            if genes:
                source_mode = genes[0].get("source", "streamed_fasta_context")
                accession = genes[0].get("accession", accession)
                _log(
                    f"  ✅  Recovered {len(genes)} genes around {protein_id[:40]} "
                    f"from streamed FASTA context"
                )

        # ── Record placement outcome ────────────────────────────────────────
        if genes:
            for g in genes:
                rows.append({
                    "hit_protein_id": protein_id,
                    "accession":      accession,
                    "flank_gene_id":  g.get("protein_id", ""),
                    "position_rel":   g.get("position_rel", 0),
                    "gene_name":      g.get("gene_name", ""),
                    "function":       g.get("function", "hypothetical protein"),
                    "strand":         g.get("strand", "+"),
                    "start":          g.get("start", 0),
                    "end":            g.get("end", 0),
                    "color":          g.get("color", _FUNCTION_COLORS["unknown"]),
                    "source":         g.get("source", source_mode or "unknown"),
                    "protein_seq":    g.get("protein_seq", ""),
                })
            report_rows.append({
                "protein_id": protein_id, "genome_id": genome_id,
                "accession": accession, "status": "placed",
                "reason": f"Found {len(genes)} flanking genes via {source_mode}",
                "rescue_hint": "", "placed": True,
            })
        else:
            # Classify why nothing was returned
            status, reason, hint = classify_missing_context(protein_id, genome_id)
            if status == "no_coords":
                reason = (
                    f"Accession {accession} found but no flanking genes returned. "
                    "seq_from/seq_to are 0 — coordinate data is missing from the "
                    "tblout. Try NCBI Entrez mode with a valid email."
                )
            elif not use_entrez:
                hint = "Provide your NCBI email to enable Entrez fallback."
                reason += " (Entrez disabled — no email provided)"
            report_rows.append({
                "protein_id": protein_id, "genome_id": genome_id,
                "accession": accession, "status": status,
                "reason": reason, "rescue_hint": hint, "placed": False,
            })
            _log(f"  ⚠️  No genes for {accession[:40]} — {status}")

    syn_df     = pd.DataFrame(rows)     if rows        else _EMPTY_SYN.copy()
    report_df  = pd.DataFrame(report_rows) if report_rows else _EMPTY_REPORT.copy()

    if not rows:
        syn_df = _EMPTY_SYN.copy()
    else:
        syn_df = _normalise_synteny_df(syn_df)

    return syn_df, report_df


# ---------------------------------------------------------------------------
# Build synteny table by calling scripts 18 and 24 as subprocesses
# ---------------------------------------------------------------------------

def run_script18(
    proj_dir: Path,
    master_hits_path: Path,
    email: str,
    flanks: int,
    log_callback=None,
) -> Path:
    """
    Call 18_synteny_analysis.py as a subprocess.

    Returns path to the produced synteny_table.tsv.
    """
    import subprocess, shutil

    scripts_dir = Path(proj_dir).parent.parent / "scripts"
    script18 = scripts_dir / "18_synteny_analysis.py"
    if not script18.exists():
        # Try adjacent to app
        script18 = Path(proj_dir).parent / "18_synteny_analysis.py"
    if not script18.exists():
        raise FileNotFoundError(f"18_synteny_analysis.py not found near {proj_dir}")

    python = shutil.which("python3") or sys.executable
    cmd = [
        python, str(script18),
        "--email", email,
        "--flanks", str(flanks),
    ]
    if log_callback:
        log_callback(f"$ {' '.join(cmd)}")

    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=str(proj_dir.parent.parent),   # run from project root
    )
    if log_callback:
        for line in (proc.stdout + proc.stderr).splitlines():
            log_callback(line)

    out_path = proj_dir.parent.parent / "results" / "synteny_table.tsv"
    if out_path.exists():
        return out_path
    raise RuntimeError(f"Script 18 exited {proc.returncode}; no synteny_table.tsv produced")


def run_script24(proj_dir: Path, log_callback=None) -> tuple[Path, Path]:
    """
    Call 24_synteny_figure.py as a subprocess.

    Returns (synteny_frequency.png, synteny_genemap.png).
    """
    import subprocess, shutil

    scripts_dir = Path(proj_dir).parent.parent / "scripts"
    script24 = scripts_dir / "24_synteny_figure.py"
    if not script24.exists():
        script24 = Path(proj_dir).parent / "24_synteny_figure.py"
    if not script24.exists():
        raise FileNotFoundError(f"24_synteny_figure.py not found near {proj_dir}")

    python = shutil.which("python3") or sys.executable
    cmd = [python, str(script24)]
    if log_callback:
        log_callback(f"$ {' '.join(cmd)}")

    proc = subprocess.run(
        cmd, capture_output=True, text=True,
        cwd=str(proj_dir.parent.parent),
    )
    if log_callback:
        for line in (proc.stdout + proc.stderr).splitlines():
            log_callback(line)

    fig_dir = proj_dir.parent.parent / "results" / "figures"
    freq  = fig_dir / "synteny_frequency.png"
    genemap = fig_dir / "synteny_genemap.png"
    if not freq.exists() or not genemap.exists():
        raise RuntimeError(f"Script 24 exited {proc.returncode}; figures not produced")
    return freq, genemap


# ---------------------------------------------------------------------------
# Conservation scoring
# ---------------------------------------------------------------------------


def _normalise_synteny_df(df):
    """Ensure synteny DataFrame has expected columns, accepting common aliases."""
    import pandas as _pd
    if df is None or df.empty:
        return df
    df = df.copy()
    # hit_protein_id
    if 'hit_protein_id' not in df.columns:
        for alias in ('protein_id', 'target_name', 'query_name'):
            if alias in df.columns:
                df['hit_protein_id'] = df[alias]
                break
        else:
            df['hit_protein_id'] = 'unknown'
    # position_rel
    if 'position_rel' not in df.columns:
        df['position_rel'] = 0
    # gene_name
    if 'gene_name' not in df.columns:
        for alias in ('flank_gene_id', 'gene'):
            if alias in df.columns:
                df['gene_name'] = df[alias]
                break
        else:
            df['gene_name'] = 'unknown'
    # is_hit
    if 'is_hit' not in df.columns:
        df['is_hit'] = False
    # end (alias for stop)
    if 'end' not in df.columns and 'stop' in df.columns:
        df['end'] = df['stop']
    elif 'end' not in df.columns:
        df['end'] = df.get('start', _pd.Series([0] * len(df)))
    # start
    if 'start' not in df.columns:
        df['start'] = 0
    # function (gene function description)
    if 'function' not in df.columns:
        for alias in ('product', 'description', 'annotation'):
            if alias in df.columns:
                df['function'] = df[alias]
                break
        else:
            df['function'] = 'hypothetical protein'
    # strand
    if 'strand' not in df.columns:
        df['strand'] = '+'
    # accession
    if 'accession' not in df.columns:
        df['accession'] = df.get('genome_id', _pd.Series(['unknown'] * len(df)))
    # protein_seq — optional; present when source GenBank had /translation
    if 'protein_seq' not in df.columns:
        df['protein_seq'] = ''
    df['protein_seq'] = df['protein_seq'].fillna('').astype(str)

    # display_label: visible gene label used by static and interactive plots.
    # Fill blank/NaN gene names with useful identifiers so figures do not show
    # empty labels or literal "nan".
    df['display_label'] = df.apply(lambda row: _gene_display_label(row, max_len=16), axis=1)
    df['gene_name'] = df.apply(
        lambda row: _clean_label_value(row.get('gene_name', '')) or row.get('display_label', 'gene'),
        axis=1,
    )
    return df


def _normalise_conservation_df(df):
    """Ensure conservation DataFrame has expected columns."""
    if df is None or df.empty:
        import pandas as _pd
        return _pd.DataFrame(columns=['position_rel','gene_name','function','presence_fraction','conservation_fraction','is_core'])
    df = df.copy()
    if 'conservation_fraction' not in df.columns:
        for alias in ('conservation', 'fraction', 'score'):
            if alias in df.columns:
                df['conservation_fraction'] = df[alias]
                break
        else:
            df['conservation_fraction'] = 1.0
    if 'presence_fraction' not in df.columns:
        df['presence_fraction'] = df.get('conservation_fraction', 1.0)
    if 'is_core' not in df.columns:
        df['is_core'] = df['conservation_fraction'] >= 0.9
    if 'function' not in df.columns:
        for alias in ('product','description'):
            if alias in df.columns:
                df['function'] = df[alias]
                break
        else:
            df['function'] = 'unknown'
    if 'gene_name' not in df.columns:
        df['gene_name'] = 'unknown'
    df['gene_name'] = df['gene_name'].map(lambda v: _clean_label_value(v) or "gene")
    if 'position_rel' not in df.columns:
        df['position_rel'] = 0
    return df

def conservation_scores(synteny_df: pd.DataFrame) -> pd.DataFrame:
    """Per-position conservation statistics across genomes."""
    _EMPTY = pd.DataFrame(
        columns=["position_rel", "gene_name", "function",
                 "presence_fraction", "conservation_fraction", "is_core"]
    )
    if synteny_df is None or synteny_df.empty:
        return _EMPTY

    synteny_df = _normalise_synteny_df(synteny_df)
    n_genomes = synteny_df["hit_protein_id"].nunique()
    if n_genomes == 0:
        return _EMPTY

    records: list[dict] = []
    for pos, grp in synteny_df.groupby("position_rel"):
        func_counts = grp["function"].value_counts()
        top_func    = func_counts.index[0] if not func_counts.empty else "unknown"
        top_func_n  = func_counts.iloc[0]  if not func_counts.empty else 0
        name_counts = grp["gene_name"].value_counts()
        top_name    = name_counts.index[0] if not name_counts.empty else ""
        present_n   = grp["hit_protein_id"].nunique()
        presence    = present_n / n_genomes
        conserv     = top_func_n / present_n if present_n else 0.0
        records.append({
            "position_rel":          int(pos),
            "gene_name":             top_name,
            "function":              top_func,
            "presence_fraction":     round(presence, 3),
            "conservation_fraction": round(conserv,  3),
            "is_core":               conserv > 0.70,
        })

    return pd.DataFrame(records).sort_values("position_rel").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Matplotlib figure  (gene-arrow + conservation panel)
# ---------------------------------------------------------------------------

def synteny_figure_matplotlib(
    synteny_df: pd.DataFrame,
    conservation_df: pd.DataFrame,
    max_genomes: int = 15,
    flanks: int = 5,
    hit_label: str = "hit gene",
    fmt: str = "png",
    dpi: int = 300,
) -> bytes:
    """
    Publication-ready gene-arrow synteny figure.

    Parameters
    ----------
    fmt : str
        Output format: ``"png"`` (300 dpi raster), ``"svg"`` (vector),
        or ``"pdf"`` (vector, embeddable in manuscripts).
    dpi : int
        Raster resolution — only meaningful for ``fmt="png"``.
        Default 300 (print quality).  Use 150 for quick screen previews.

    Returns
    -------
    bytes
        Encoded figure bytes in the requested format, or ``b""`` on failure.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.backends.backend_pdf import PdfPages  # noqa: F401 — ensure PDF backend
    except ImportError:
        print("ERROR: matplotlib not installed.", file=sys.stderr)
        return b""

    fmt = fmt.lower().lstrip(".")
    if fmt not in ("png", "svg", "pdf"):
        fmt = "png"

    if synteny_df is None or synteny_df.empty:
        return b""

    synteny_df = _normalise_synteny_df(synteny_df)
    genomes   = synteny_df["hit_protein_id"].unique()[:max_genomes]
    n_genomes = len(genomes)
    conservation_df = _normalise_conservation_df(conservation_df)
    has_cons  = conservation_df is not None and not conservation_df.empty

    # ── Coordinate-aware layout (falls back to relative-index mode) ────────
    coord_layout: dict[str, list[dict]] = {}
    coord_bounds: list[float] = []
    for _gid in genomes:
        _layout, _ok = _coordinate_layout_for_genome(
            synteny_df[synteny_df["hit_protein_id"] == _gid]
        )
        if _ok:
            coord_layout[str(_gid)] = _layout
            for _item in _layout:
                coord_bounds.extend([float(_item["x0"]), float(_item["x1"])])
    coordinate_mode = bool(coord_layout)

    # ── Sizing: scale for print (A4-landscape friendly at ~180 mm wide) ──
    fig_w      = 14.0     # inches — ~178 mm; fits single column if scaled 50%
    row_h      = 0.70     # inch per genome row (generous for readability)
    fig_height = max(5.0, n_genomes * row_h + 3.0)
    h_ratios   = [1, n_genomes] if has_cons else [n_genomes]
    n_rows     = 2 if has_cons else 1

    # Use serifless fonts common in journals
    plt.rcParams.update({
        "font.family":      "sans-serif",
        "font.size":        8,
        "axes.linewidth":   0.8,
        "xtick.major.width":0.6,
        "ytick.major.width":0.6,
        "pdf.fonttype":     42,   # TrueType in PDF — avoids Type 3 font warnings
        "svg.fonttype":     "none",   # editable text in SVG (Inkscape/Illustrator)
    })

    fig, axes = plt.subplots(
        n_rows, 1,
        figsize=(fig_w, fig_height),
        gridspec_kw={"height_ratios": h_ratios},
    )
    if n_rows == 1:
        axes = [axes]
    ax_cons, ax_genes = (axes[0], axes[1]) if has_cons else (None, axes[0])

    # ── Conservation bar panel ─────────────────────────────────────────────
    if has_cons and ax_cons is not None:
        positions  = conservation_df["position_rel"].values
        fractions  = conservation_df["conservation_fraction"].values
        bar_colors = [
            "#4CAF50" if f > 0.70 else "#FFC107" if f > 0.40 else "#9E9E9E"
            for f in fractions
        ]
        ax_cons.bar(positions, fractions, color=bar_colors, edgecolor="white", linewidth=0.6,
                    zorder=3)
        ax_cons.axhline(0.70, color="#D32F2F", linestyle="--", lw=0.9,
                        label="Core threshold (70%)", zorder=4)
        ax_cons.set_ylim(0, 1.12)
        ax_cons.set_ylabel("Conservation", fontsize=8)
        ax_cons.set_title(
            f"Synteny Neighbourhood — Conservation & Gene Map  (n={n_genomes} loci)",
            fontsize=9, fontweight="bold", pad=6,
        )
        ax_cons.legend(fontsize=7, framealpha=0.85, edgecolor="#cccccc")
        ax_cons.tick_params(axis="both", labelsize=7)
        ax_cons.set_xticks(range(-flanks, flanks + 1))
        ax_cons.set_xlim(-flanks - 0.65, flanks + 0.65)
        ax_cons.spines[["top", "right"]].set_visible(False)

    # ── Gene arrow rows ────────────────────────────────────────────────────
    if coordinate_mode and coord_bounds:
        _xmin, _xmax = min(coord_bounds), max(coord_bounds)
        _pad = max(0.35, (_xmax - _xmin) * 0.05)
        ax_genes.set_xlim(_xmin - _pad, _xmax + _pad)
        ax_genes.set_xlabel("Position relative to query hit (kb)", fontsize=8)
    else:
        ax_genes.set_xlim(-flanks - 0.9, flanks + 0.9)
        ax_genes.set_xlabel("Relative gene position", fontsize=8)
    ax_genes.set_ylim(-0.6, n_genomes - 0.4)
    ax_genes.set_yticks(range(n_genomes))
    ax_genes.set_yticklabels([str(g)[:38] for g in genomes], fontsize=6.5)
    ax_genes.tick_params(axis="x", labelsize=7)
    ax_genes.grid(axis="x", linestyle=":", linewidth=0.35, alpha=0.55, zorder=0)
    if not coordinate_mode:
        ax_genes.set_xticks(range(-flanks, flanks + 1))
    ax_genes.spines[["top", "right"]].set_visible(False)

    arrow_h = 0.38
    head_len = 0.15

    for row_i, genome_id in enumerate(genomes):
        genome_genes = synteny_df[synteny_df["hit_protein_id"] == genome_id]
        if coordinate_mode and str(genome_id) in coord_layout:
            for item in coord_layout[str(genome_id)]:
                g       = item["row"]
                col     = str(g.get("color", "#9E9E9E"))
                strand  = _strand_symbol(g.get("strand", "+"))
                label   = _gene_display_label(g, hit_label=hit_label, max_len=14)
                is_hit  = float(g.get("position_rel", 99) or 99) == 0
                x0, x1  = float(item["x0"]), float(item["x1"])
                y       = row_i + float(item.get("lane_offset", 0.0))
                width   = max(abs(x1 - x0), 0.03)
                dx      = width if strand == "+" else -width
                x_start = min(x0, x1) if strand == "+" else max(x0, x1)
                h       = max(0.16, arrow_h * 0.62)
                local_head = min(max(width * 0.25, 0.025), 0.18)
                arrow = mpatches.FancyArrow(
                    x=x_start, y=y, dx=dx, dy=0,
                    width=h, head_width=h * 1.35, head_length=local_head,
                    length_includes_head=True,
                    facecolor=col,
                    edgecolor="#111111" if is_hit else col,
                    linewidth=2.0 if is_hit else 0.5,
                    alpha=1.0 if is_hit else 0.88, zorder=2,
                )
                ax_genes.add_patch(arrow)
                ax_genes.text(
                    float(item["x_mid"]), y + h * 0.85, label,
                    ha="center", va="bottom",
                    fontsize=4.8 if not is_hit else 5.5,
                    fontweight="bold" if is_hit else "normal",
                    color="#B71C1C" if is_hit else "#333333",
                    clip_on=True,
                )
        else:
            for _, g in genome_genes.iterrows():
                pos    = float(g["position_rel"])
                col    = str(g.get("color", "#9E9E9E"))
                strand = _strand_symbol(g.get("strand", "+"))
                label  = _gene_display_label(g, hit_label=hit_label, max_len=14)
                is_hit = (pos == 0)
                dx      = 0.82 if strand == "+" else -0.82
                edge_c  = "#111111" if is_hit else col
                lw      = 2.0    if is_hit else 0.5
                alpha   = 1.0    if is_hit else 0.88
                arrow = mpatches.FancyArrow(
                    x=pos - dx * 0.46, y=row_i, dx=dx * 0.82, dy=0,
                    width=arrow_h, head_width=arrow_h * 1.35, head_length=head_len,
                    length_includes_head=True,
                    facecolor=col, edgecolor=edge_c, linewidth=lw,
                    alpha=alpha, zorder=2,
                )
                ax_genes.add_patch(arrow)
                if abs(pos) <= flanks:
                    ax_genes.text(
                        pos, row_i + arrow_h * 0.82, label,
                        ha="center", va="bottom",
                        fontsize=4.8 if not is_hit else 5.5,
                        fontweight="bold" if is_hit else "normal",
                        color="#B71C1C" if is_hit else "#333333",
                        clip_on=True,
                    )

    # ── Functional colour legend ───────────────────────────────────────────
    legend_cats = [k for k in _FUNCTION_COLORS
                   if k not in ("unknown", "hypothetical", "putative")][:9]
    legend_items = [
        mpatches.Patch(facecolor=_FUNCTION_COLORS[k], edgecolor="#555",
                       linewidth=0.5, label=k.capitalize())
        for k in legend_cats
    ]
    ax_genes.legend(handles=legend_items, loc="lower right",
                    fontsize=5.5, ncol=3, framealpha=0.88,
                    edgecolor="#cccccc", handlelength=1.2)

    plt.tight_layout(pad=1.2)

    buf = io.BytesIO()
    save_kw: dict = {"format": fmt, "bbox_inches": "tight"}
    if fmt == "png":
        save_kw["dpi"] = dpi
    fig.savefig(buf, **save_kw)
    plt.close(fig)
    buf.seek(0)
    return buf.read()


# ---------------------------------------------------------------------------
# Multi-format export helper  (PNG 300 dpi  +  SVG  +  PDF)
# ---------------------------------------------------------------------------

def export_synteny_figures(
    synteny_df: pd.DataFrame,
    conservation_df: pd.DataFrame,
    out_dir: Path,
    hit_label: str = "hit gene",
    max_genomes: int = 15,
    flanks: int = 5,
    dpi: int = 300,
    log_callback=None,
) -> dict:
    """
    Save the synteny figure in three publication-ready formats.

    Files written:
      * ``synteny_map.png``  — 300 dpi raster (microscopy journals, Word)
      * ``synteny_map.svg``  — editable vector (Inkscape / Illustrator)
      * ``synteny_map.pdf``  — embeddable PDF (LaTeX / Keynote / PowerPoint)

    Returns
    -------
    dict  mapping ``"png" | "svg" | "pdf"``  →  ``Path``  for each produced file.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    produced: dict = {}

    for fmt in ("png", "svg", "pdf"):
        try:
            data = synteny_figure_matplotlib(
                synteny_df, conservation_df,
                max_genomes=max_genomes,
                flanks=flanks,
                hit_label=hit_label,
                fmt=fmt,
                dpi=dpi,
            )
            if data:
                p = out_dir / f"synteny_map.{fmt}"
                p.write_bytes(data)
                produced[fmt] = p
                if log_callback:
                    log_callback(
                        f"  ✅  synteny_map.{fmt.upper()}"
                        f"  ({len(data) // 1024} KB)"
                        + (f"  @ {dpi} dpi" if fmt == "png" else "  [vector]")
                    )
        except Exception as exc:
            if log_callback:
                log_callback(f"  ⚠️  {fmt.upper()} export failed: {exc}")

    return produced


# ---------------------------------------------------------------------------
# Plotly interactive figure
# ---------------------------------------------------------------------------

def synteny_figure_plotly(
    synteny_df: pd.DataFrame,
    conservation_df: pd.DataFrame,
    max_genomes: int = 20,
    hit_label: str = "hit gene",
):
    """Interactive plotly synteny figure with hover info."""
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError:
        print("ERROR: plotly not installed.", file=sys.stderr)
        class _E:
            def to_json(self): return "{}"
        return _E()

    if synteny_df is None or synteny_df.empty:
        return go.Figure()

    synteny_df = _normalise_synteny_df(synteny_df)
    genomes   = synteny_df["hit_protein_id"].unique()[:max_genomes]
    n_genomes = len(genomes)
    conservation_df = _normalise_conservation_df(conservation_df)
    has_cons  = conservation_df is not None and not conservation_df.empty

    # ── Coordinate-aware layout ────────────────────────────────────────────
    coord_layout: dict[str, list[dict]] = {}
    coord_bounds: list[float] = []
    for _gid in genomes:
        _layout, _ok = _coordinate_layout_for_genome(
            synteny_df[synteny_df["hit_protein_id"] == _gid]
        )
        if _ok:
            coord_layout[str(_gid)] = _layout
            for _item in _layout:
                coord_bounds.extend([float(_item["x0"]), float(_item["x1"])])
    coordinate_mode = bool(coord_layout)

    fig = make_subplots(
        rows=2 if has_cons else 1, cols=1,
        row_heights=[0.18, 0.82] if has_cons else [1.0],
        vertical_spacing=0.04,
        shared_xaxes=False,
    )

    if has_cons:
        bar_colors = [
            "#4CAF50" if f > 0.70 else "#FFC107" if f > 0.40 else "#9E9E9E"
            for f in conservation_df["conservation_fraction"]
        ]
        fig.add_trace(go.Bar(
            x=conservation_df["position_rel"],
            y=conservation_df["conservation_fraction"],
            marker_color=bar_colors,
            name="Conservation",
            hovertemplate="Position %{x} · Conservation %{y:.0%}<extra></extra>",
        ), row=1, col=1)
        fig.add_hline(y=0.70, line_dash="dash", line_color="#F44336",
                      annotation_text="Core (70%)", row=1, col=1)

    gene_row = 2 if has_cons else 1
    hover_traces: list = []

    def _arrow_path(x0: float, x1: float, y: float, strand: str, height: float = 0.28) -> str:
        left, right = min(x0, x1), max(x0, x1)
        length = max(right - left, 0.04)
        head   = min(max(length * 0.28, 0.025), 0.18)
        y0p, y1p = y - height / 2.0, y + height / 2.0
        if strand == "-":
            return (f"M {right} {y0p} L {left + head} {y0p} L {left} {y} "
                    f"L {left + head} {y1p} L {right} {y1p} Z")
        return (f"M {left} {y0p} L {right - head} {y0p} L {right} {y} "
                f"L {right - head} {y1p} L {left} {y1p} Z")

    shapes: list = []
    annotations: list = []

    for row_i, genome_id in enumerate(genomes):
        genome_genes = synteny_df[synteny_df["hit_protein_id"] == genome_id]
        if coordinate_mode and str(genome_id) in coord_layout:
            for item in coord_layout[str(genome_id)]:
                g      = item["row"]
                col    = str(g.get("color", "#9E9E9E"))
                s      = _strand_symbol(g.get("strand", "+"))
                label  = _gene_display_label(g, hit_label=hit_label, max_len=18)
                pid    = _clean_label_value(g.get("flank_gene_id", g.get("protein_id", "")))
                func   = _clean_label_value(g.get("function", ""))
                source = str(g.get("source", ""))
                pos    = float(g.get("position_rel", 0) or 0)
                is_hit = pos == 0
                y      = row_i + float(item.get("lane_offset", 0.0))
                x0, x1, xm = float(item["x0"]), float(item["x1"]), float(item["x_mid"])
                start, end  = int(item["start"]), int(item["end"])
                overlap_note = str(item.get("overlap_note", ""))

                shapes.append(dict(
                    type="path",
                    path=_arrow_path(x0, x1, y, s),
                    fillcolor=col,
                    line=dict(color="#212121" if is_hit else col,
                              width=2 if is_hit else 0.7),
                    layer="below",
                ))
                hover_traces.append(go.Scatter(
                    x=[xm], y=[y],
                    mode="markers",
                    marker=dict(size=22, opacity=0, color=col),
                    hovertemplate=(
                        f"<b>{label}</b><br>"
                        f"Protein ID: {pid or 'not available'}<br>"
                        f"Function: {func or 'not annotated'}<br>"
                        f"Relative position: {int(pos):+d}<br>"
                        f"Coordinates: {start}..{end} nt<br>"
                        f"Strand: {s} ({'forward' if s == '+' else 'reverse'})<br>"
                        f"{overlap_note + '<br>' if overlap_note else ''}"
                        f"Source: {source}<br>"
                        f"Genome: {str(genome_id)[:50]}<extra></extra>"
                    ),
                    showlegend=False,
                ))
                annotations.append(dict(
                    x=xm, y=y, text=label, showarrow=False,
                    font=dict(size=8 if not is_hit else 9,
                              color="#B71C1C" if is_hit else "#263238"),
                    yshift=16,
                ))
                if is_hit:
                    annotations.append(dict(
                        x=xm, y=y, text="▼", showarrow=False,
                        font=dict(size=10, color="#F44336"), yshift=-17,
                    ))
        else:
            for _, g in genome_genes.iterrows():
                pos    = float(g["position_rel"])
                col    = str(g.get("color", "#9E9E9E"))
                s      = _strand_symbol(g.get("strand", "+"))
                label  = _gene_display_label(g, hit_label=hit_label, max_len=18)
                pid    = _clean_label_value(g.get("flank_gene_id", g.get("protein_id", "")))
                func   = _clean_label_value(g.get("function", ""))
                source = str(g.get("source", ""))
                is_hit = pos == 0

                shapes.append(dict(
                    type="rect",
                    x0=pos - 0.42, x1=pos + 0.42,
                    y0=row_i - 0.30, y1=row_i + 0.30,
                    fillcolor=col,
                    line=dict(color="#212121" if is_hit else col,
                              width=2 if is_hit else 0.5),
                    layer="below",
                ))
                hover_traces.append(go.Scatter(
                    x=[pos], y=[row_i],
                    mode="markers",
                    marker=dict(size=20, opacity=0, color=col),
                    hovertemplate=(
                        f"<b>{label}</b><br>"
                        f"Protein ID: {pid or 'not available'}<br>"
                        f"Function: {func}<br>"
                        f"Position: {int(pos):+d}<br>"
                        f"Strand: {s}<br>"
                        f"Source: {source}<br>"
                        f"Genome: {str(genome_id)[:50]}<extra></extra>"
                    ),
                    showlegend=False,
                ))
                annotations.append(dict(
                    x=pos, y=row_i, text=label, showarrow=False,
                    font=dict(size=8 if not is_hit else 9,
                              color="#B71C1C" if is_hit else "#263238"),
                    yshift=16,
                ))
                if is_hit:
                    annotations.append(dict(
                        x=pos, y=row_i, text="▼", showarrow=False,
                        font=dict(size=10, color="#F44336"), yshift=-18,
                    ))

    for trace in hover_traces:
        fig.add_trace(trace, row=gene_row, col=1)

    x_title = "Position relative to query hit (kb)" if coordinate_mode else "Relative gene position"
    fig.update_layout(
        shapes=shapes, annotations=annotations,
        title="Synteny Neighbourhood (interactive — hover for gene details)",
        height=max(400, 80 + n_genomes * 38),
        plot_bgcolor="white", paper_bgcolor="white",
        font=dict(size=11), showlegend=False,
        xaxis=dict(title=x_title, dtick=1 if not coordinate_mode else None),
        yaxis=dict(
            tickvals=list(range(n_genomes)),
            ticktext=[str(g)[:40] for g in genomes],
            autorange="reversed",
        ),
        margin=dict(l=240, r=30, t=60, b=50),
    )
    return fig


# ---------------------------------------------------------------------------
# Neighbourhood GenBank builder  (input for clinker / pyGenomeViz / EasyFig)
# ---------------------------------------------------------------------------

def build_neighborhood_genbanks(
    synteny_df: pd.DataFrame,
    out_dir: Path,
    min_genes: int = 2,
    log_callback=None,
) -> list:
    """
    Create one synthetic GenBank file per hit neighbourhood.

    Required by clinker, pyGenomeViz, and EasyFig — all three expect GenBank
    files as primary input.

    Parameters
    ----------
    synteny_df : pd.DataFrame
        Long-format table produced by ``build_synteny_table()``.
    out_dir : Path
        Directory where .gbk files are written.
    min_genes : int
        Skip loci with fewer than this many CDS features.

    Returns
    -------
    List of Path objects for the produced .gbk files.
    """
    try:
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
    except ImportError:
        msg = "ERROR: BioPython required — pip install biopython"
        if log_callback:
            log_callback(msg)
        print(msg, file=sys.stderr)
        return []

    if synteny_df is None or synteny_df.empty:
        return []

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    produced: list = []

    synteny_df = _normalise_synteny_df(synteny_df)
    for hit_id, group in synteny_df.groupby("hit_protein_id"):
        if len(group) < min_genes:
            continue

        safe_name = re.sub(r"[^A-Za-z0-9._-]", "_", str(hit_id))[:40]
        gbk_path  = out_dir / f"{safe_name}.gbk"

        # Anchor coordinates
        starts = group["start"].dropna().astype(int)
        ends   = group["end"].dropna().astype(int)
        min_start = int(starts.min()) if not starts.empty else 0
        max_end   = int(ends.max())   if not ends.empty   else 1000
        if max_end <= min_start:
            max_end = min_start + 300 * len(group)
        seq_len = max_end - min_start + 200

        record = SeqRecord(
            Seq("N" * seq_len),
            id=safe_name,
            name=safe_name[:16],
            description=f"Neighbourhood of {hit_id}",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )

        for _, row in group.iterrows():
            g_start = max(0, int(row.get("start", min_start) or min_start) - min_start)
            g_end   = max(g_start + 10, int(row.get("end", min_start + 300) or min_start + 300) - min_start)
            strand_v = 1 if str(row.get("strand", "+")) == "+" else -1

            qualifiers: dict = {
                "product":   [str(row.get("function", "hypothetical protein"))],
                "locus_tag": [str(row.get("flank_gene_id", ""))],
            }
            if row.get("gene_name"):
                qualifiers["gene"] = [str(row["gene_name"])]
            if int(row.get("position_rel", -99)) == 0:
                qualifiers["note"] = ["HMM discovery hit"]
            prot_seq = str(row.get("protein_seq", "") or "").strip()
            if prot_seq:
                qualifiers["translation"] = [prot_seq]

            # clinker expects a gene feature alongside each CDS to avoid
            # "Could not find parent gene" warnings.
            gene_qualifiers: dict = {}
            if qualifiers.get("gene"):
                gene_qualifiers["gene"] = qualifiers["gene"]
            elif qualifiers.get("locus_tag"):
                gene_qualifiers["gene"] = qualifiers["locus_tag"]
            record.features.append(SeqFeature(
                FeatureLocation(g_start, g_end, strand=strand_v),
                type="gene",
                qualifiers=gene_qualifiers,
            ))
            feat = SeqFeature(
                FeatureLocation(g_start, g_end, strand=strand_v),
                type="CDS",
                qualifiers=qualifiers,
            )
            record.features.append(feat)

        try:
            SeqIO.write(record, str(gbk_path), "genbank")
            produced.append(gbk_path)
            if log_callback:
                log_callback(f"  Written: {gbk_path.name}")
        except Exception as exc:
            msg = f"WARNING: Could not write {gbk_path}: {exc}"
            if log_callback:
                log_callback(msg)
            print(msg, file=sys.stderr)

    return produced


# ---------------------------------------------------------------------------
# clinker wrapper
# ---------------------------------------------------------------------------

def run_clinker(
    gbk_dir: Path,
    out_dir: Path,
    log_callback=None,
    identity: float = 0.3,
) -> "Optional[Path]":
    """
    Run clinker on all .gbk files in *gbk_dir* and produce an interactive HTML.

    Requires: ``pip install clinker``

    Parameters
    ----------
    gbk_dir : Path
        Directory containing .gbk neighbourhood files.
    out_dir : Path
        Output directory; ``clinker_output.html`` is written here.
    identity : float
        Minimum protein identity threshold for linking (0–1).

    Returns
    -------
    Path to HTML, or ``None`` on failure.
    """
    import subprocess, shutil

    gbk_files = sorted(Path(gbk_dir).glob("*.gbk"))
    if len(gbk_files) < 2:
        msg = f"⚠️  clinker needs ≥2 GenBank files; found {len(gbk_files)}"
        if log_callback:
            log_callback(msg)
        return None

    try:
        from pipeline.utils import find_tool as _find_tool
        clinker_exe = _find_tool("clinker")
    except Exception:
        clinker_exe = shutil.which("clinker")
    if not clinker_exe:
        if log_callback:
            log_callback(
                "❌  clinker not on PATH.\n"
                "    Install with:  pip install clinker"
            )
        return None

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    html_out = out_dir / "clinker_output.html"

    cmd = [
        clinker_exe,
        *[str(f) for f in gbk_files],
        "--plot",     str(html_out),
        "--identity", str(identity),
        "--jobs",     "1",
    ]
    if log_callback:
        log_callback(f"$ {' '.join(cmd)}")

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if log_callback:
        for line in (proc.stdout + proc.stderr).splitlines():
            if line.strip():
                log_callback(line)

    if html_out.exists():
        if log_callback:
            log_callback(f"✅  clinker HTML → {html_out}")
        return html_out

    if log_callback:
        log_callback(f"❌  clinker exited {proc.returncode}; no HTML produced.")
    return None


# ---------------------------------------------------------------------------
# pyGenomeViz wrapper
# ---------------------------------------------------------------------------

def run_pygenomeviz(
    synteny_df: pd.DataFrame,
    out_dir: Path,
    log_callback=None,
    max_genomes: int = 15,
) -> bytes:
    """
    Build a synteny figure with pyGenomeViz (Python API, no subprocess).

    Requires: ``pip install pygenomeviz``

    Parameters
    ----------
    synteny_df : pd.DataFrame
        Long-format neighbourhood table.
    out_dir : Path
        Output directory; ``synteny_pygenomeviz.png`` is written here.
    max_genomes : int
        Cap on number of loci shown.

    Returns
    -------
    PNG bytes, or ``b""`` on failure.
    """
    try:
        from pygenomeviz import GenomeViz
    except ImportError:
        if log_callback:
            log_callback(
                "❌  pygenomeviz not installed.\n"
                "    Install with:  pip install pygenomeviz"
            )
        return b""

    try:
        import io as _io
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        if synteny_df is None or synteny_df.empty:
            if log_callback:
                log_callback("⚠️  Empty synteny table — nothing to plot.")
            return b""

        genomes = synteny_df["hit_protein_id"].unique()[:max_genomes]
        n = len(genomes)
        if n == 0:
            return b""

        if log_callback:
            log_callback(f"Building pyGenomeViz figure for {n} loci …")

        gv = GenomeViz(
            genome_track_ratio=0.15,
            fig_track_height=0.5,
            feature_track_ratio=0.25,
        )

        for gid in genomes:
            grp    = synteny_df[synteny_df["hit_protein_id"] == gid]
            starts = grp["start"].dropna().astype(int)
            ends   = grp["end"].dropna().astype(int)
            g_min  = int(starts.min()) if not starts.empty else 0
            g_max  = int(ends.max())   if not ends.empty   else 1000
            span   = max(g_max - g_min + 100, 100)

            label = str(gid)[:30]
            track = gv.add_feature_track(label, span)

            for _, row in grp.iterrows():
                s      = max(0, int(row.get("start", g_min) or g_min) - g_min)
                e      = max(s + 10, int(row.get("end",   g_min) or g_min) - g_min)
                strand = 1 if str(row.get("strand", "+")) == "+" else -1
                color  = str(row.get("color", "#9E9E9E"))
                is_hit = int(row.get("position_rel", -99)) == 0
                glabel = _gene_display_label(row, max_len=14)

                track.add_feature(
                    start=s,
                    end=e,
                    strand=strand,
                    fc=color,
                    ec="#212121" if is_hit else color,
                    lw=2.0 if is_hit else 0.5,
                    label=glabel or None,
                    label_size=6,
                )

        fig = gv.plotfig(figsize=(14, max(4, n * 1.4 + 2)))

        buf = _io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        png_bytes = buf.read()

        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "synteny_pygenomeviz.png").write_bytes(png_bytes)

        if log_callback:
            log_callback(f"✅  pyGenomeViz PNG saved ({len(png_bytes)//1024} KB)")

        return png_bytes

    except Exception as exc:
        import traceback
        if log_callback:
            log_callback(f"❌  pyGenomeViz error: {exc}")
            log_callback(traceback.format_exc())
        return b""


# ---------------------------------------------------------------------------
# EasyFig wrapper
# ---------------------------------------------------------------------------

def run_easyfig(
    gbk_dir: Path,
    out_dir: Path,
    log_callback=None,
) -> "Optional[Path]":
    """
    Run EasyFig on .gbk neighbourhood files.

    EasyFig is **not** on pip; install from https://github.com/mjsull/Easyfig
    and place ``easyfig`` (or ``Easyfig.py``) on PATH.

    Returns
    -------
    Path to produced SVG or PNG, or ``None`` on failure.
    """
    import subprocess, shutil

    gbk_files = sorted(Path(gbk_dir).glob("*.gbk"))
    if not gbk_files:
        if log_callback:
            log_callback("⚠️  No .gbk files found for EasyFig.")
        return None

    try:
        from pipeline.utils import find_tool as _find_tool
        easyfig_exe = (
            _find_tool("easyfig")
            or _find_tool("Easyfig.py")
            or _find_tool("easyfig.py")
            or _find_tool("EasyFig")
        )
    except Exception:
        easyfig_exe = (
            shutil.which("easyfig")
            or shutil.which("Easyfig.py")
            or shutil.which("easyfig.py")
            or shutil.which("EasyFig")
        )
    if not easyfig_exe:
        if log_callback:
            log_callback(
                "❌  EasyFig not found on PATH.\n"
                "    Download from https://github.com/mjsull/Easyfig\n"
                "    and add to PATH as 'easyfig' or 'Easyfig.py'."
            )
        return None

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_svg = out_dir / "synteny_easyfig.svg"

    cmd = [easyfig_exe] + [str(f) for f in gbk_files] + ["-o", str(out_svg)]
    if log_callback:
        log_callback(f"$ {' '.join(cmd)}")

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if log_callback:
        for line in (proc.stdout + proc.stderr).splitlines():
            if line.strip():
                log_callback(line)

    if out_svg.exists():
        if log_callback:
            log_callback(f"✅  EasyFig SVG → {out_svg}")
        return out_svg

    out_png = out_dir / "synteny_easyfig.png"
    if out_png.exists():
        if log_callback:
            log_callback(f"✅  EasyFig PNG → {out_png}")
        return out_png

    if log_callback:
        log_callback(f"❌  EasyFig exited {proc.returncode}; no output produced.")
    return None


# ---------------------------------------------------------------------------
# GFF3 export
# ---------------------------------------------------------------------------

def export_gff3(synteny_df: pd.DataFrame, out_path: Path) -> None:
    """Write synteny neighbourhood as GFF3 (loadable in Artemis, IGV, Geneious)."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if synteny_df is None or synteny_df.empty:
        out_path.write_text("##gff-version 3\n")
        return

    lines = ["##gff-version 3", "# Generated by HMM Discovery App — synteny.py"]
    for _, row in synteny_df.iterrows():
        seqid    = str(row.get("accession", ".")).replace(" ", "_") or "."
        start    = int(row.get("start", 1))
        end      = int(row.get("end", start))
        strand   = str(row.get("strand", "."))
        if strand not in ("+", "-"):
            strand = "."
        gene_id  = str(row.get("flank_gene_id", "")).replace(";", "%3B") or "."
        name     = str(row.get("gene_name",     "")).replace(";", "%3B") or "."
        func     = str(row.get("function",      "")).replace(";", "%3B") or "hypothetical protein"
        pos_rel  = int(row.get("position_rel", 0))
        hit_pid  = str(row.get("hit_protein_id","")).replace(";", "%3B") or "."
        attrs    = (f"ID={gene_id};Name={name};Note={func};"
                    f"position_rel={pos_rel};hit_protein_id={hit_pid}")
        lines.append("\t".join([seqid, "synteny_pipeline", "CDS",
                                 str(start), str(end), ".", strand, ".", attrs]))

    out_path.write_text("\n".join(lines) + "\n")
