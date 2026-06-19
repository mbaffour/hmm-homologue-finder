#!/usr/bin/env python3
"""Run an exhaustive, resumable all-database HMM Discovery benchmark.

This script is intentionally separate from the Shiny app so a long validation
run can continue outside a browser session. It writes all research outputs to a
benchmark directory outside the deployable repository.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import hashlib
import io
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import time
import warnings
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Iterable

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from core.env_setup import check_environment  # noqa: E402
from databases.builtin import BUILTIN_DATABASES  # noqa: E402
from databases.downloader import check_ncbi_ftp_files, format_bytes  # noqa: E402
from pipeline.alignment import alignment_quality, run_mafft, run_trimal  # noqa: E402
from pipeline.clustering import cluster_cdhit, cluster_summary  # noqa: E402
from pipeline.confidence import add_qc_flags, classify_hits  # noqa: E402
from pipeline.hmm_builder import run_hmmbuild, self_search_recovery  # noqa: E402
from pipeline.input_handler import input_summary  # noqa: E402
from pipeline.matrix import build_matrix, heatmap_png  # noqa: E402
from pipeline.motifs import run_fimo, run_meme  # noqa: E402
from pipeline.phylo import render_tree, run_iqtree  # noqa: E402
from pipeline.reporter import (  # noqa: E402
    build_report_context,
    build_reproducibility_json,
    create_export_zip,
    generate_methods_text,
    render_html_report,
)
from pipeline.searcher import parse_tblout  # noqa: E402
from pipeline.synteny import (  # noqa: E402
    build_neighborhood_genbanks,
    build_synteny_table,
    conservation_scores,
    export_gff3,
    export_synteny_figures,
)
from pipeline.taxonomy import taxonomy_table  # noqa: E402
from pipeline.utils import find_tool, get_env  # noqa: E402


DEFAULT_FASTA = ROOT / "example_data" / "demo_protein_family.fasta"
DEFAULT_OUT = Path.home() / "hmm_homologue_finder_results"
SMOKE_DBS = {"INPHARED proteins", "SwissProt"}
PARTIAL_DBS = {"INPHARED proteins", "SwissProt", "RefSeq viral proteins"}
NUCLEOTIDE_MODES = {"nucleotide"}


def safe_name(value: str, limit: int = 160) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "unknown"))[:limit].strip("_")


def q(value: os.PathLike | str) -> str:
    return shlex.quote(str(value))


def now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def count_tblout_hits(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open(errors="replace") as fh:
        return sum(1 for line in fh if line.strip() and not line.startswith("#"))


def disk_free_gb(path: Path) -> float:
    path.mkdir(parents=True, exist_ok=True)
    return shutil.disk_usage(path).free / (1024**3)


def require_free_space(path: Path, min_free_gb: float) -> None:
    free = disk_free_gb(path)
    if free < min_free_gb:
        raise RuntimeError(
            f"Low disk space at {path}: {free:.1f} GiB free, "
            f"{min_free_gb:.1f} GiB required. Free disk or resume later."
        )


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_cmd(
    cmd: str,
    *,
    desc: str,
    cwd: Path,
    logs_dir: Path,
    timeout: int | None = None,
) -> tuple[int, float, Path, Path]:
    logs_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    stem = safe_name(f"{stamp}_{desc}", 220)
    stdout_path = logs_dir / f"{stem}.stdout.txt"
    stderr_path = logs_dir / f"{stem}.stderr.txt"
    start = time.time()
    with stdout_path.open("w") as out, stderr_path.open("w") as err:
        proc = subprocess.run(
            ["bash", "-lc", cmd],
            cwd=str(cwd),
            env=get_env(),
            stdout=out,
            stderr=err,
            text=True,
            timeout=timeout,
        )
    return proc.returncode, time.time() - start, stdout_path, stderr_path


def write_tsv(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


class Benchmark:
    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.out = Path(args.out).resolve()
        self.fasta = Path(args.fasta).resolve()
        self.data = self.out / "data"
        self.alignments = self.out / "alignments"
        self.hmm_dir = self.out / "hmm"
        self.search_results = self.out / "search_results"
        self.results = self.out / "results"
        self.figures = self.out / "figures"
        self.reports = self.out / "reports"
        self.logs = self.out / "logs"
        self.cache = self.out / "cache"
        self.runtime = self.out / "runtime"
        self.manifest_path = self.out / "benchmark_manifest.json"
        self.metrics_path = self.results / "per_database_metrics.tsv"
        self.summary_path = self.results / "all_database_summary.tsv"
        self.manifest = self._load_manifest()

    def _load_manifest(self) -> dict:
        if self.manifest_path.exists():
            try:
                return json.loads(self.manifest_path.read_text())
            except Exception:
                pass
        return {
            "created_at": now_iso(),
            "updated_at": now_iso(),
            "status": "created",
            "root": str(self.out),
            "fasta": str(self.fasta),
            "databases": {},
            "core": {},
            "failures": [],
        }

    def save_manifest(self) -> None:
        self.manifest["updated_at"] = now_iso()
        self.manifest_path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.manifest_path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(self.manifest, indent=2, default=str))
        tmp.replace(self.manifest_path)

    def log(self, message: str) -> None:
        line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}"
        print(line, flush=True)
        self.logs.mkdir(parents=True, exist_ok=True)
        with (self.logs / "benchmark.log").open("a") as fh:
            fh.write(line + "\n")

    def init_dirs(self) -> None:
        for d in [
            self.data,
            self.alignments,
            self.hmm_dir,
            self.search_results,
            self.results,
            self.figures,
            self.reports,
            self.logs,
            self.cache,
            self.runtime,
        ]:
            d.mkdir(parents=True, exist_ok=True)

    def preflight(self) -> None:
        self.init_dirs()
        if ROOT in self.out.parents or self.out == ROOT:
            raise RuntimeError(
                "Benchmark output root must not be inside the deployable repository."
            )
        if not self.fasta.exists():
            raise FileNotFoundError(f"Input FASTA not found: {self.fasta}")
        require_free_space(self.out, self.args.min_free_gb)
        env = check_environment()
        self.manifest["environment"] = {
            "all_required_ok": env.get("all_required_ok"),
            "all_python_ok": env.get("all_python_ok"),
            "all_full_run_ok": env.get("all_full_run_ok"),
            "missing_full_run_tools": [
                t.get("name") for t in env.get("missing_full_run_tools", [])
            ],
            "missing_full_run_python": [
                p.get("pkg") for p in env.get("missing_full_run_python", [])
            ],
        }
        self.save_manifest()
        if not env.get("all_full_run_ok"):
            raise RuntimeError("Full-run environment check failed; see manifest.")
        for tool in ["hmmsearch", "hmmbuild", "mafft", "trimal", "curl"]:
            if not find_tool(tool):
                raise RuntimeError(f"Required tool not found: {tool}")

    def selected_databases(self) -> list[dict]:
        if self.args.preset == "smoke":
            wanted = SMOKE_DBS
        elif self.args.preset == "partial":
            wanted = PARTIAL_DBS
        elif self.args.databases:
            wanted = {d.strip() for d in self.args.databases.split(",") if d.strip()}
        else:
            wanted = {db["name"] for db in BUILTIN_DATABASES}
        return [db for db in BUILTIN_DATABASES if db["name"] in wanted]

    def resolve_urls(self, db: dict) -> list[str]:
        url = db.get("download_url") or db.get("url") or ""
        if not url:
            return []
        if "*" not in url:
            return [url]
        base = url[: url.rfind("/") + 1]
        pattern = url[url.rfind("/") + 1 :]
        return check_ncbi_ftp_files(base, pattern)

    def dry_run(self) -> None:
        rows = []
        for db in self.selected_databases():
            urls = self.resolve_urls(db)
            rows.append(
                {
                    "database": db["name"],
                    "type": db.get("type"),
                    "search_mode": db.get("search_mode"),
                    "optional": db.get("optional"),
                    "size_hint": db.get("size_hint"),
                    "relevance": db.get("relevance", ""),
                    "file_count": len(urls),
                    "first_url": urls[0] if urls else "",
                }
            )
            self.log(f"DRY {db['name']}: {len(urls)} files")
        write_tsv(rows, self.results / "dry_run_database_expansion.tsv")
        self.manifest["status"] = "dry_run_complete"
        self.manifest["dry_run"] = rows
        self.save_manifest()

    def build_core(self) -> tuple[Path, Path, dict]:
        if self.manifest.get("core", {}).get("status") == "complete" and not self.args.force:
            hmm = Path(self.manifest["core"]["hmm_path"])
            trimmed = Path(self.manifest["core"]["trimmed_alignment"])
            if hmm.exists() and trimmed.exists():
                self.log("Core HMM already complete; resuming.")
                return hmm, trimmed, self.manifest["core"]

        copied = self.data / self.fasta.name
        if not copied.exists() or self.args.force:
            shutil.copy2(self.fasta, copied)
        summary = input_summary(copied)
        if summary.get("seq_count", 0) <= 0:
            raise RuntimeError(f"No sequences parsed from {copied}")
        self.log(f"Input summary: {summary}")
        aln = run_mafft(copied, self.alignments / "seed.mafft.faa", cpu=self.args.cpu)
        if not aln or not aln.exists():
            raise RuntimeError("MAFFT failed")
        trimmed = run_trimal(aln, self.alignments / "seed.mafft.trimmed.faa")
        if not trimmed or not trimmed.exists():
            self.log("trimAl failed; using raw MAFFT alignment.")
            trimmed = aln
        quality = alignment_quality(trimmed)
        hmm = self.hmm_dir / "benchmark_profile.hmm"
        hmm_info = run_hmmbuild(trimmed, hmm, "benchmark_profile")
        if not hmm_info:
            raise RuntimeError("hmmbuild failed")
        recovery = self_search_recovery(hmm, copied)
        if recovery.get("total", 0) <= 0 or recovery.get("recovery_rate", 0) < self.args.min_recovery:
            raise RuntimeError(f"Self-search recovery failed: {recovery}")
        core = {
            "status": "complete",
            "input_path": str(copied),
            "input_summary": summary,
            "alignment_quality": quality,
            "hmm_info": hmm_info,
            "self_search": recovery,
            "hmm_path": str(hmm),
            "trimmed_alignment": str(trimmed),
        }
        self.manifest["core"] = core
        self.save_manifest()
        shutil.copy2(copied, self.results / "hits_proteins.faa")
        shutil.copy2(trimmed, self.results / "hits_aligned.faa")
        return hmm, trimmed, core

    def cache_url(self, url: str, db_name: str, file_idx: int) -> tuple[Path, int, dict]:
        require_free_space(self.out, self.args.min_free_gb)
        label = safe_name(Path(url.split("?")[0]).name or f"file_{file_idx}")
        db_cache = self.cache / safe_name(db_name)
        db_cache.mkdir(parents=True, exist_ok=True)
        dest = db_cache / f"{file_idx:04d}_{label}"
        marker = Path(str(dest) + ".complete")
        meta_path = Path(str(dest) + ".meta.json")
        before = dest.stat().st_size if dest.exists() else 0
        if marker.exists() and dest.exists() and dest.stat().st_size > 0:
            meta = self.cache_metadata(url, db_name, file_idx, dest, meta_path, 0, None)
            return dest, 0, meta
        cmd = (
            f"curl -sS -f -L -C - --retry 20 --retry-delay 10 "
            f"--retry-all-errors --connect-timeout 30 -o {q(dest)} {q(url)} "
            f"&& touch {q(marker)}"
        )
        rc, elapsed, _, stderr = run_cmd(
            cmd,
            desc=f"download_{safe_name(db_name)}_{file_idx}",
            cwd=self.out,
            logs_dir=self.logs,
            timeout=self.args.download_timeout,
        )
        if rc != 0 or not marker.exists():
            tail = stderr.read_text(errors="replace")[-1200:] if stderr.exists() else ""
            raise RuntimeError(f"Download failed for {db_name} file {file_idx}: {tail}")
        after = dest.stat().st_size
        self.log(
            f"Cached {db_name} file {file_idx + 1}: {format_bytes(after)} "
            f"in {elapsed:.1f}s"
        )
        downloaded = max(after - before, 0)
        meta = self.cache_metadata(url, db_name, file_idx, dest, meta_path, downloaded, elapsed)
        return dest, downloaded, meta

    def cache_metadata(
        self,
        url: str,
        db_name: str,
        file_idx: int,
        dest: Path,
        meta_path: Path,
        downloaded_bytes: int,
        elapsed: float | None,
    ) -> dict:
        try:
            if meta_path.exists():
                meta = json.loads(meta_path.read_text())
            else:
                meta = {}
        except Exception:
            meta = {}
        current_size = dest.stat().st_size
        if meta.get("sha256") and meta.get("size_bytes") == current_size:
            checksum = meta["sha256"]
        else:
            checksum = sha256_file(dest)
        accessed_at = meta.get("accessed_at") or now_iso()
        meta.update(
            {
                "url": url,
                "database": db_name,
                "file_index": file_idx,
                "cached_path": str(dest),
                "size_bytes": current_size,
                "sha256": checksum,
                "accessed_at": accessed_at,
                "downloaded_bytes": downloaded_bytes,
            }
        )
        if elapsed is not None:
            meta["download_elapsed_seconds"] = round(float(elapsed), 3)
        try:
            meta_path.write_text(json.dumps(meta, indent=2, default=str))
        except Exception as exc:
            self.log(f"WARNING: Could not write cache metadata for {db_name}: {exc}")
        return meta

    def cleanup_cache(self, path: Path) -> None:
        if self.args.keep_cache:
            return
        try:
            path.unlink(missing_ok=True)
            Path(str(path) + ".complete").unlink(missing_ok=True)
            Path(str(path) + ".meta.json").unlink(missing_ok=True)
        except Exception as exc:
            self.log(f"WARNING: Could not remove cache {path}: {exc}")

    def preserve_synteny_context_records(self, cache_file: Path, tbl: Path) -> int:
        """Save compact FASTA records for nucleotide hits before raw DB cleanup."""
        try:
            from Bio import SeqIO
        except Exception:
            return 0
        parsed = parse_tblout(tbl)
        if parsed.empty or "description" not in parsed.columns:
            return 0
        coords = parsed["description"].astype(str).str.extract(
            r"coords=([^:]+):(\d+)-(\d+)\(([+-])\)"
        )
        accessions = {
            str(value).strip()
            for value in coords[0].dropna().tolist()
            if str(value).strip()
        }
        if not accessions:
            return 0

        context_dir = self.results / "synteny_context_cache"
        context_dir.mkdir(parents=True, exist_ok=True)
        wanted = set(accessions)
        wanted.update(acc.split(".")[0] for acc in accessions if acc)
        saved = 0
        try:
            handle = (
                gzip.open(cache_file, "rt")
                if str(cache_file).endswith(".gz")
                else cache_file.open()
            )
            with handle:
                for rec in SeqIO.parse(handle, "fasta"):
                    rec_ids = {str(rec.id), str(rec.id).split()[0], str(rec.id).split(".")[0]}
                    matches = {acc for acc in accessions if acc in rec_ids or acc.split(".")[0] in rec_ids}
                    if not matches:
                        continue
                    for acc in matches:
                        out = context_dir / f"{safe_name(acc, 120)}.fna"
                        if not out.exists():
                            SeqIO.write(rec, str(out), "fasta")
                            saved += 1
                    if saved >= len(accessions):
                        break
        except Exception as exc:
            self.log(f"WARNING: could not preserve synteny context from {cache_file.name}: {exc}")
        if saved:
            self.log(f"Preserved {saved} compact synteny context records from {cache_file.name}")
        return saved

    def run_protein_hmmsearch(
        self, cache_file: Path, hmm: Path, tbl: Path, cpu: int
    ) -> tuple[int, float]:
        source = f"gzip -cd {q(cache_file)}" if cache_file.suffix == ".gz" else f"cat {q(cache_file)}"
        cmd = (
            "set -o pipefail; "
            f"{source} | {q(find_tool('hmmsearch') or 'hmmsearch')} "
            f"--tblout {q(tbl)} -E {self.args.evalue} --cpu {cpu} --noali {q(hmm)} -"
        )
        rc, elapsed, _, stderr = run_cmd(
            cmd,
            desc=f"search_{safe_name(tbl.stem)}",
            cwd=self.out,
            logs_dir=self.logs,
            timeout=self.args.search_timeout,
        )
        if rc != 0:
            tail = stderr.read_text(errors="replace")[-1200:] if stderr.exists() else ""
            raise RuntimeError(f"hmmsearch failed: {tail}")
        return count_tblout_hits(tbl), elapsed

    def run_nucleotide_hmmsearch(
        self, cache_file: Path, hmm: Path, tbl: Path, gff: Path, cpu: int
    ) -> tuple[int, float]:
        prodigal = find_tool("prodigal")
        seqkit = find_tool("seqkit")
        hmmsearch = find_tool("hmmsearch") or "hmmsearch"
        nt_mode = self.args.nt_orf_mode
        if nt_mode == "prodigal" and not prodigal:
            raise RuntimeError("Nucleotide search requires prodigal for this benchmark")
        work = self.runtime / f"{safe_name(tbl.stem)}_chunks"
        start = time.monotonic()
        shutil.rmtree(work, ignore_errors=True)
        chunks_dir = work / "chunks"
        plain_dir = work / "plain"
        proteins_dir = work / "proteins"
        for directory in (chunks_dir, plain_dir, proteins_dir):
            directory.mkdir(parents=True, exist_ok=True)

        def _decompress_to_plain(in_file: Path, plain: Path) -> None:
            if str(in_file).endswith(".gz"):
                with gzip.open(in_file, "rb") as src, plain.open("wb") as dst:
                    shutil.copyfileobj(src, dst)
            else:
                shutil.copy2(in_file, plain)

        def _write_sixframe_orfs(nt_fasta: Path, faa: Path, part_gff: Path | None = None) -> int:
            from Bio import BiopythonWarning, SeqIO
            from Bio.Seq import Seq

            stops = {"TAA", "TAG", "TGA"}
            min_nt = max(1, int(self.args.min_orf_aa)) * 3
            n_orfs = 0
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonWarning)
                with faa.open("w") as out_faa:
                    out_gff = part_gff.open("w") if part_gff is not None else None
                    if out_gff is not None:
                        out_gff.write("##gff-version 3\n")
                    for record in SeqIO.parse(str(nt_fasta), "fasta"):
                        seq_id = str(record.id)
                        seq = str(record.seq).upper().replace("U", "T")
                        seq_len = len(seq)
                        strands = [("+", seq), ("-", str(record.seq.reverse_complement()).upper().replace("U", "T"))]
                        for strand, scan_seq in strands:
                            for frame in range(3):
                                segment_start = frame
                                last_full = frame + ((len(scan_seq) - frame) // 3) * 3
                                for pos in range(frame, last_full, 3):
                                    codon = scan_seq[pos:pos + 3]
                                    if codon in stops:
                                        if pos - segment_start >= min_nt:
                                            n_orfs += 1
                                            self._emit_sixframe_orf(
                                                seq_id, scan_seq, seq_len, strand, frame,
                                                segment_start, pos, n_orfs, out_faa, out_gff
                                            )
                                        segment_start = pos + 3
                                if last_full - segment_start >= min_nt:
                                    n_orfs += 1
                                    self._emit_sixframe_orf(
                                        seq_id, scan_seq, seq_len, strand, frame,
                                        segment_start, last_full, n_orfs, out_faa, out_gff
                                    )
                    if out_gff is not None:
                        out_gff.close()
            return n_orfs

        def _translate_chunk(in_file: Path) -> tuple[Path, Path, str]:
            base = in_file.name
            plain = plain_dir / f"{base}.fna"
            faa = proteins_dir / f"{base}.faa"
            part_gff = proteins_dir / f"{base}.gff"
            _decompress_to_plain(in_file, plain)
            try:
                if nt_mode == "sixframe":
                    # Do not write all-ORF GFF for exhaustive mode. It can reach
                    # tens of GB on large phage databases; coordinates are kept in
                    # FASTA descriptions and copied into compact hit tables later.
                    n_orfs = _write_sixframe_orfs(plain, faa, None)
                    msg = "" if n_orfs else f"No six-frame ORFs >= {self.args.min_orf_aa} aa"
                    return faa, None, msg
                proc = subprocess.run(
                    [
                        str(prodigal),
                        "-i",
                        str(plain),
                        "-a",
                        str(faa),
                        "-o",
                        str(part_gff),
                        "-f",
                        "gff",
                        "-p",
                        "meta",
                        "-q",
                    ],
                    capture_output=True,
                    text=True,
                    timeout=self.args.search_timeout,
                )
            finally:
                plain.unlink(missing_ok=True)
            msg = ""
            if proc.returncode != 0:
                msg = (proc.stderr or proc.stdout or "").strip()
            return faa, part_gff, msg

        if seqkit:
            parts = max(cpu, 1)
            split_cmd = (
                f"{q(seqkit)} split2 -p {parts} -O {q(chunks_dir)} {q(cache_file)}"
            )
            rc, _, _, stderr = run_cmd(
                split_cmd,
                desc=f"split_{safe_name(tbl.stem)}",
                cwd=self.out,
                logs_dir=self.logs,
                timeout=self.args.search_timeout,
            )
            if rc != 0:
                tail = stderr.read_text(errors="replace")[-1600:] if stderr.exists() else ""
                raise RuntimeError(f"seqkit split failed: {tail}")
            fasta_suffixes = (".fa", ".faa", ".fasta", ".fna", ".fa.gz", ".faa.gz", ".fasta.gz", ".fna.gz")
            chunk_files = sorted(
                p for p in chunks_dir.iterdir()
                if p.is_file() and p.name.lower().endswith(fasta_suffixes)
            )
        else:
            chunk_file = chunks_dir / cache_file.name
            shutil.copy2(cache_file, chunk_file)
            chunk_files = [chunk_file]

        if not chunk_files:
            raise RuntimeError(f"No nucleotide chunks created for {cache_file}")

        failures: list[str] = []
        if nt_mode == "sixframe":
            total_hits = 0
            tbl.write_text("")
            total_chunks = len(chunk_files)
            for chunk_idx, chunk in enumerate(chunk_files, start=1):
                try:
                    self.log(
                        f"{tbl.stem} chunk {chunk_idx}/{total_chunks}: "
                        f"translating {chunk.name} with six-frame ORFs"
                    )
                    faa, _, msg = _translate_chunk(chunk)
                    if msg:
                        failures.append(f"{chunk.name}: {msg[-500:]}")
                    if not faa.exists() or faa.stat().st_size == 0:
                        self.log(
                            f"{tbl.stem} chunk {chunk_idx}/{total_chunks}: "
                            "no translated ORFs passed the length cutoff"
                        )
                        continue
                    self.log(
                        f"{tbl.stem} chunk {chunk_idx}/{total_chunks}: "
                        f"searching {format_bytes(faa.stat().st_size)} translated ORFs"
                    )
                    part_tbl = proteins_dir / f"{chunk.stem}.tblout"
                    search_cmd = (
                        f"{q(hmmsearch)} --tblout {q(part_tbl)} -E {self.args.evalue} "
                        f"--cpu {cpu} --noali {q(hmm)} {q(faa)}"
                    )
                    rc, _, _, stderr = run_cmd(
                        search_cmd,
                        desc=f"search_{safe_name(tbl.stem)}_chunk{chunk_idx:04d}",
                        cwd=self.out,
                        logs_dir=self.logs,
                        timeout=self.args.search_timeout,
                    )
                    if rc != 0:
                        tail = stderr.read_text(errors="replace")[-1600:] if stderr.exists() else ""
                        raise RuntimeError(f"nucleotide hmmsearch failed on {chunk.name}: {tail}")
                    chunk_hits = 0
                    with tbl.open("a") as out_tbl, part_tbl.open(errors="replace") as in_tbl:
                        for line in in_tbl:
                            if line.strip() and not line.startswith("#"):
                                out_tbl.write(line)
                                chunk_hits += 1
                                total_hits += 1
                    self.log(
                        f"{tbl.stem} chunk {chunk_idx}/{total_chunks}: "
                        f"{chunk_hits} hits, {total_hits} cumulative"
                    )
                    faa.unlink(missing_ok=True)
                    part_tbl.unlink(missing_ok=True)
                    chunk.unlink(missing_ok=True)
                except Exception as exc:
                    failures.append(f"{chunk.name}: {exc}")
                    raise
            if total_hits == 0 and failures:
                self.log("WARNING: six-frame nucleotide search completed with no hits; " + "; ".join(failures[-3:]))
            elapsed = time.monotonic() - start
            if not self.args.keep_translation_chunks:
                shutil.rmtree(work, ignore_errors=True)
            gff.write_text(
                "##gff-version 3\n"
                "# Exhaustive six-frame benchmark does not retain all-ORF GFF files.\n"
                "# Hit coordinates are stored in hits_main.tsv from tblout descriptions.\n"
            )
            return total_hits, elapsed

        with concurrent.futures.ThreadPoolExecutor(max_workers=max(cpu, 1)) as executor:
            future_map = {executor.submit(_translate_chunk, p): p for p in chunk_files}
            for future in concurrent.futures.as_completed(future_map):
                chunk = future_map[future]
                try:
                    _, _, msg = future.result()
                    if msg:
                        failures.append(f"{chunk.name}: {msg[-500:]}")
                except Exception as exc:
                    failures.append(f"{chunk.name}: {exc}")

        faa_files = sorted(p for p in proteins_dir.glob("*.faa") if p.stat().st_size > 0)
        if not faa_files:
            detail = "; ".join(failures[-6:]) if failures else f"{nt_mode} produced no proteins."
            raise RuntimeError(f"No proteins translated from nucleotide DB chunks. {detail}")

        with gff.open("wb") as out_gff:
            for part_gff in sorted(proteins_dir.glob("*.gff")):
                if part_gff.stat().st_size > 0:
                    out_gff.write(part_gff.read_bytes())
        combined_faa = proteins_dir / "translated_all.faa"
        with combined_faa.open("wb") as out_faa:
            for faa in faa_files:
                out_faa.write(faa.read_bytes())

        search_cmd = (
            f"{q(hmmsearch)} --tblout {q(tbl)} -E {self.args.evalue} "
            f"--cpu {cpu} --noali {q(hmm)} {q(combined_faa)}"
        )
        rc, _, _, stderr = run_cmd(
            search_cmd,
            desc=f"search_{safe_name(tbl.stem)}",
            cwd=self.out,
            logs_dir=self.logs,
            timeout=self.args.search_timeout,
        )
        if rc != 0:
            tail = stderr.read_text(errors="replace")[-1600:] if stderr.exists() else ""
            raise RuntimeError(f"nucleotide hmmsearch failed: {tail}")
        elapsed = time.monotonic() - start
        if not self.args.keep_translation_chunks:
            shutil.rmtree(work, ignore_errors=True)
        return count_tblout_hits(tbl), elapsed

    @staticmethod
    def _emit_sixframe_orf(
        seq_id: str,
        scan_seq: str,
        seq_len: int,
        strand: str,
        frame: int,
        segment_start: int,
        segment_end: int,
        n_orfs: int,
        out_faa,
        out_gff,
    ) -> None:
        from Bio.Seq import Seq

        aa = str(Seq(scan_seq[segment_start:segment_end]).translate(to_stop=False)).replace("*", "X")
        if not aa:
            return
        if strand == "+":
            start_1 = segment_start + 1
            end_1 = segment_end
            frame_label = f"+{frame + 1}"
        else:
            start_1 = seq_len - segment_end + 1
            end_1 = seq_len - segment_start
            frame_label = f"-{frame + 1}"
        orf_id = f"{safe_name(seq_id, 80)}_sixframe_orf{n_orfs:08d}"
        header = (
            f"{orf_id} coords={seq_id}:{start_1}-{end_1}({strand}) "
            f"frame={frame_label} nt_start={start_1} nt_end={end_1} aa_len={len(aa)}"
        )
        out_faa.write(f">{header}\n")
        for i in range(0, len(aa), 80):
            out_faa.write(aa[i:i + 80] + "\n")
        if out_gff is not None:
            out_gff.write(
                f"{seq_id}\tHMMDiscovery\tsixframe_ORF\t{start_1}\t{end_1}\t.\t"
                f"{strand}\t{frame}\tID={orf_id};Name={orf_id};aa_len={len(aa)}\n"
            )

    def run_pfam_hmmscan(self, db: dict, input_faa: Path) -> dict:
        urls = self.resolve_urls(db)
        if not urls:
            raise RuntimeError(f"{db['name']} hmmscan URL not configured")
        cache_file, downloaded, provenance = self.cache_url(urls[0], db["name"], 0)
        handler = db.get("setup_handler") or safe_name(db["name"])
        setup = self.out / "db_setup" / safe_name(handler)
        setup.mkdir(parents=True, exist_ok=True)
        filename = Path(urls[0].split("?")[0]).name
        is_tar_gz = filename.endswith((".tar.gz", ".tgz"))
        is_gz = filename.endswith(".gz") and not is_tar_gz
        if is_tar_gz:
            hmm = setup / ("vfam.hmm" if "vfam" in filename.lower() else f"{safe_name(db['name'])}.hmm")
        elif is_gz:
            hmm = setup / filename[:-3]
        else:
            hmm = setup / filename
        if not hmm.exists():
            if is_tar_gz:
                with tarfile.open(cache_file, "r:gz") as tf:
                    members = tf.getmembers()
                    hmm_members = [
                        m for m in members
                        if m.isfile() and m.name.lower().endswith(".hmm")
                    ]
                    if hmm_members:
                        # Raw text .hmm files (Pfam, VOGDB VFAM): concatenate them.
                        with hmm.open("wb") as out_hmm:
                            for member in hmm_members:
                                src = tf.extractfile(member)
                                if src is None:
                                    continue
                                shutil.copyfileobj(src, out_hmm)
                                out_hmm.write(b"\n")
                        self.log(f"{db['name']}: extracted/concatenated {len(hmm_members)} HMM file(s)")
                    else:
                        # Pre-pressed binary .h3m (e.g. PHROGs via the Pharokka
                        # bundle): take the largest .h3m and copy it out. hmmpress
                        # (below) re-creates the .h3i/.h3f/.h3p aux files from it.
                        h3m_members = [
                            m for m in members
                            if m.isfile() and m.name.lower().endswith(".h3m")
                        ]
                        if not h3m_members:
                            raise RuntimeError(f"No .hmm or .h3m files found in {cache_file}")
                        h3m_members.sort(key=lambda m: m.size, reverse=True)
                        src = tf.extractfile(h3m_members[0])
                        with hmm.open("wb") as out_hmm:
                            shutil.copyfileobj(src, out_hmm)
                        self.log(f"{db['name']}: extracted pre-pressed {h3m_members[0].name}; will hmmpress")
            elif is_gz:
                rc, _, _, stderr = run_cmd(
                    f"gzip -cd {q(cache_file)} > {q(hmm)}",
                    desc=f"{safe_name(db['name'])}_decompress",
                    cwd=self.out,
                    logs_dir=self.logs,
                    timeout=self.args.search_timeout,
                )
                if rc != 0:
                    raise RuntimeError(stderr.read_text(errors="replace")[-1200:])
            else:
                shutil.copy2(cache_file, hmm)
        if not self.args.keep_cache:
            self.cleanup_cache(cache_file)
        if not Path(str(hmm) + ".h3i").exists():
            rc, _, _, stderr = run_cmd(
                f"{q(find_tool('hmmpress') or 'hmmpress')} -f {q(hmm)}",
                desc=f"{safe_name(db['name'])}_hmmpress",
                cwd=self.out,
                logs_dir=self.logs,
                timeout=self.args.search_timeout,
            )
            if rc != 0:
                raise RuntimeError(stderr.read_text(errors="replace")[-1200:])
        db_key = safe_name(db["name"])
        tbl = self.search_results / f"{db_key}.tblout"
        domtbl = self.search_results / f"{db_key}.domtblout"
        rc, elapsed, _, stderr = run_cmd(
            f"{q(find_tool('hmmscan') or 'hmmscan')} --domtblout {q(domtbl)} "
            f"--tblout {q(tbl)} -E {self.args.evalue} --cpu {self.args.cpu} "
            f"--noali {q(hmm)} {q(input_faa)}",
            desc=f"{db_key}_hmmscan",
            cwd=self.out,
            logs_dir=self.logs,
            timeout=self.args.search_timeout,
        )
        if rc != 0:
            raise RuntimeError(stderr.read_text(errors="replace")[-1200:])
        hits = count_tblout_hits(tbl)
        annotation_tsv = ""
        annotation_provenance = {}
        annotation_url = db.get("annotation_url") or ""
        if db.get("setup_handler") == "vogdb_hmmscan" and annotation_url:
            annotation_cache, annotation_downloaded, annotation_provenance = self.cache_url(
                annotation_url, f"{db['name']} annotations", 1
            )
            annotation_file = setup / Path(annotation_url.split("?")[0]).name
            if not annotation_file.exists():
                shutil.copy2(annotation_cache, annotation_file)
            if not self.args.keep_cache:
                self.cleanup_cache(annotation_cache)
            annotation_tsv = str(
                self.write_vogdb_annotation_table(
                    tbl,
                    domtbl,
                    annotation_file,
                    self.results / "vogdb_vfam_annotation.tsv",
                )
            )
            downloaded += annotation_downloaded
        return {
            "database": db["name"],
            "status": "complete",
            "optional": db.get("optional", False),
            "type": db.get("type"),
            "search_mode": db.get("search_mode"),
            "setup_handler": db.get("setup_handler"),
            "release": db.get("release", ""),
            "file_count": 1,
            "downloaded_bytes": downloaded,
            "hit_count": hits,
            "strict_count": hits,
            "runtime_seconds": elapsed,
            "tblout": str(tbl),
            "domtblout": str(domtbl),
            "annotation_tsv": annotation_tsv,
            "kept_cache": str(hmm),
            "files": {
                "0000": {
                    **provenance,
                    "status": "complete",
                    "hit_count": hits,
                    "strict_count": hits,
                    "runtime_seconds": elapsed,
                    "tblout": str(tbl),
                    "domtblout": str(domtbl),
                    "annotation": annotation_provenance,
                }
            },
        }

    def write_vogdb_annotation_table(
        self,
        tblout: Path,
        domtblout: Path,
        annotation_file: Path,
        out_tsv: Path,
    ) -> Path:
        annotations = self.load_vogdb_annotations(annotation_file)
        qcov = self.parse_domtbl_query_coverage(domtblout)
        rows = []
        if tblout.exists():
            with tblout.open(encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if not line.strip() or line.startswith("#"):
                        continue
                    parts = line.split(maxsplit=18)
                    if len(parts) < 6:
                        continue
                    vfam_id = parts[0]
                    query_id = parts[2] if len(parts) > 2 else ""
                    ann = annotations.get(vfam_id, {})
                    rows.append(
                        {
                            "query_protein_id": query_id,
                            "vfam_id": vfam_id,
                            "evalue": parts[4],
                            "bit_score": parts[5],
                            "query_coverage": qcov.get((query_id, vfam_id), ""),
                            "annotation": ann.get("annotation", ""),
                            "function": ann.get("function", ""),
                            "category": ann.get("category", ""),
                        }
                    )
        write_tsv(rows, out_tsv)
        return out_tsv

    @staticmethod
    def parse_domtbl_query_coverage(path: Path) -> dict[tuple[str, str], float]:
        qcov = {}
        if not path.exists():
            return qcov
        with path.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split(maxsplit=22)
                if len(parts) < 19:
                    continue
                try:
                    target, query = parts[0], parts[3]
                    qlen = float(parts[5])
                    ali_from, ali_to = int(parts[17]), int(parts[18])
                    cov = round(max(0, ali_to - ali_from + 1) / qlen, 4) if qlen else 0.0
                except Exception:
                    continue
                key = (query, target)
                if key not in qcov or cov > qcov[key]:
                    qcov[key] = cov
        return qcov

    @staticmethod
    def load_vogdb_annotations(path: Path) -> dict[str, dict]:
        def open_text(p: Path):
            if str(p).endswith(".gz"):
                return gzip.open(p, "rt", encoding="utf-8", errors="replace")
            return p.open("r", encoding="utf-8", errors="replace")

        def pick(row: dict, header: list[str], tokens: tuple[str, ...]) -> str:
            for h in header:
                if any(token in h.lower() for token in tokens) and row.get(h):
                    return str(row[h]).strip()
            return ""

        if not path.exists():
            return {}
        with open_text(path) as fh:
            first = ""
            for line in fh:
                if line.strip() and not line.startswith("#"):
                    first = line.rstrip("\n")
                    break
            if not first:
                return {}
            raw_header = first.split("\t")
            first_lower = [value.lower().strip() for value in raw_header]
            has_header = (
                first_lower[0] in {"vfam", "vfam_id", "vog", "vog_id", "group", "group_id", "groupname", "hmm"}
                or any(value in {"function", "annotation", "category", "description"} for value in first_lower)
            )
            if has_header:
                header = raw_header
                reader = pd.read_csv(fh, sep="\t", names=header, dtype=str).fillna("")
            else:
                header = [f"column_{idx + 1}" for idx in range(len(raw_header))]
                reader = pd.read_csv(
                    io.StringIO(first + "\n" + fh.read()),
                    sep="\t",
                    names=header,
                    dtype=str,
                ).fillna("")
            lower = {h.lower().strip(): h for h in header}
            rows = {}
            for _, row_series in reader.iterrows():
                row = row_series.to_dict()
                group_id = ""
                for key in ("vfam", "vfam_id", "vog", "vog_id", "group", "group_id", "groupname", "hmm"):
                    if key in lower and row.get(lower[key]):
                        group_id = str(row[lower[key]]).strip()
                        break
                if not group_id:
                    group_id = str(row.get(header[0], "")).strip()
                if not group_id:
                    continue
                annotation = pick(row, header, ("annot", "description", "consensus", "name"))
                function = pick(row, header, ("function", "description"))
                category = pick(row, header, ("category", "class"))
                if header and header[0].startswith("column_"):
                    function = function or str(row.get("column_2", "")).strip()
                    category = category or str(row.get("column_3", "")).strip()
                    annotation = annotation or str(row.get("column_4", "")).strip() or function
                rows[group_id] = {
                    "annotation": annotation,
                    "function": function,
                    "category": category,
                }
            return rows

    def run_database(self, db: dict, hmm: Path) -> dict:
        db_name = db["name"]
        db_key = safe_name(db_name)
        existing = self.manifest["databases"].get(db_name, {})
        db_is_nt = db.get("type") in NUCLEOTIDE_MODES
        desired_nt_mode = self.args.nt_orf_mode if db_is_nt else None
        existing_mode_ok = (
            not db_is_nt
            or existing.get("nt_orf_mode") == desired_nt_mode
        )
        if existing.get("status") == "complete" and existing_mode_ok and not self.args.force:
            self.log(f"Skipping complete database: {db_name}")
            return existing
        if existing.get("status") == "complete" and db_is_nt and not existing_mode_ok:
            self.log(
                f"Re-running {db_name}: prior nucleotide ORF mode was "
                f"{existing.get('nt_orf_mode') or 'unknown'}, requested {desired_nt_mode}"
            )
        self.log(f"Starting database: {db_name}")
        self.manifest["databases"][db_name] = {
            "status": "running",
            "started_at": now_iso(),
            "optional": db.get("optional", False),
            "nt_orf_mode": desired_nt_mode,
        }
        self.save_manifest()

        if db.get("search_mode") == "hmmscan":
            result = self.run_pfam_hmmscan(db, self.results / "hits_proteins.faa")
            self.manifest["databases"][db_name] = result
            self.save_manifest()
            return result

        urls = self.resolve_urls(db)
        if not urls:
            raise RuntimeError(f"No URLs resolved for {db_name}")
        total_hits = 0
        total_strict = 0
        downloaded_bytes = 0
        runtime = 0.0
        tblouts = []
        file_metrics = []
        for idx, url in enumerate(urls):
            part_key = f"{idx:04d}"
            part_info = existing.get("files", {}).get(part_key, {})
            part_mode_ok = (
                not db_is_nt
                or part_info.get("nt_orf_mode") == desired_nt_mode
            )
            if part_info.get("status") == "complete" and part_mode_ok and not self.args.force:
                total_hits += int(part_info.get("hit_count", 0))
                total_strict += int(part_info.get("strict_count", 0))
                runtime += float(part_info.get("runtime_seconds", 0))
                downloaded_bytes += int(part_info.get("downloaded_bytes", 0) or 0)
                tblouts.append(part_info.get("tblout", ""))
                continue
            require_free_space(self.out, self.args.min_free_gb)
            cache_file, dl_bytes, provenance = self.cache_url(url, db_name, idx)
            downloaded_bytes += dl_bytes
            tbl = self.search_results / f"{db_key}_part{idx:04d}.tblout"
            try:
                if db_is_nt:
                    gff = self.search_results / f"{db_key}_part{idx:04d}.{desired_nt_mode}.gff"
                    hits, elapsed = self.run_nucleotide_hmmsearch(
                        cache_file, hmm, tbl, gff, self.args.cpu
                    )
                else:
                    hits, elapsed = self.run_protein_hmmsearch(
                        cache_file, hmm, tbl, self.args.cpu
                    )
                strict = 0
                parsed = parse_tblout(tbl)
                if not parsed.empty and "bit_score" in parsed.columns:
                    strict = int((parsed["bit_score"] >= self.args.strict_bits).sum())
                if db_is_nt and hits:
                    self.preserve_synteny_context_records(cache_file, tbl)
                total_hits += hits
                total_strict += strict
                runtime += elapsed
                tblouts.append(str(tbl))
                file_metrics.append(
                    {
                        **provenance,
                        "url": url,
                        "status": "complete",
                        "hit_count": hits,
                        "strict_count": strict,
                        "runtime_seconds": elapsed,
                        "downloaded_bytes": dl_bytes,
                        "tblout": str(tbl),
                        "nt_orf_mode": desired_nt_mode,
                    }
                )
                self.log(f"{db_name} file {idx + 1}/{len(urls)}: {hits} hits")
                self.manifest["databases"].setdefault(db_name, {}).setdefault("files", {})[
                    part_key
                ] = file_metrics[-1]
                self.save_manifest()
            finally:
                if db.get("search_mode") == "hmmsearch":
                    self.cleanup_cache(cache_file)

        merged = self.search_results / f"{db_key}.merged.tblout"
        with merged.open("w") as out_fh:
            for tbl in tblouts:
                p = Path(tbl)
                if p.exists():
                    for line in p.read_text(errors="replace").splitlines():
                        if line.strip() and not line.startswith("#"):
                            out_fh.write(line + "\n")
        result = {
            "database": db_name,
            "status": "complete",
            "optional": db.get("optional", False),
            "type": db.get("type"),
            "search_mode": db.get("search_mode"),
            "nt_orf_mode": desired_nt_mode,
            "file_count": len(urls),
            "downloaded_bytes": downloaded_bytes,
            "hit_count": total_hits,
            "strict_count": total_strict,
            "runtime_seconds": runtime,
            "tblout": str(merged),
            "files": self.manifest["databases"].get(db_name, {}).get("files", {}),
        }
        self.manifest["databases"][db_name] = result
        self.save_manifest()
        return result

    def collect_hits(self) -> pd.DataFrame:
        frames = []
        for db_name, info in self.manifest.get("databases", {}).items():
            if info.get("status") != "complete":
                continue
            if info.get("search_mode") == "hmmscan":
                continue
            db = next((d for d in BUILTIN_DATABASES if d["name"] == db_name), {})
            file_infos = list((info.get("files") or {}).values())
            if file_infos:
                for file_info in file_infos:
                    tbl = Path(file_info.get("tblout", ""))
                    if not tbl.exists():
                        continue
                    df = parse_tblout(tbl)
                    if df.empty:
                        continue
                    df["db_name"] = db_name
                    df["database_source"] = db_name
                    df["db_type"] = db.get("type", "")
                    df["source_url"] = file_info.get("url", "")
                    df["source_sha256"] = file_info.get("sha256", "")
                    df["source_accessed_at"] = file_info.get("accessed_at", "")
                    df["source_size_bytes"] = file_info.get("size_bytes", "")
                    frames.append(df)
            else:
                tbl = Path(info.get("tblout", ""))
                if not tbl.exists():
                    continue
                df = parse_tblout(tbl)
                if df.empty:
                    continue
                df["db_name"] = db_name
                df["database_source"] = db_name
                df["db_type"] = db.get("type", "")
                df["source_url"] = ""
                frames.append(df)
        if frames:
            hits = pd.concat(frames, ignore_index=True)
        else:
            hits = pd.DataFrame(
                columns=[
                    "target_name",
                    "query_name",
                    "evalue",
                    "bit_score",
                    "bias_score",
                    "description",
                    "db_name",
                    "db_type",
                    "source_url",
                ]
            )
        if not hits.empty:
            hits["protein_id"] = hits["target_name"]
            hits["hit_id"] = hits["target_name"]
            six = hits["description"].astype(str).str.extract(
                r"coords=([^:]+):(\d+)-(\d+)\(([+-])\)"
            )
            prod = hits["description"].astype(str).str.extract(
                r"^\s*#\s*(\d+)\s*#\s*(\d+)\s*#\s*(-?1)\s*#"
            )
            fallback_genome = (
                hits["target_name"]
                .astype(str)
                .str.replace(r"_sixframe_orf\d+$", "", regex=True)
                .str.replace(r"_s-?1_f\d+_o\d+$", "", regex=True)
            )
            hits["genome_id"] = six[0].where(six[0].notna() & (six[0] != ""), fallback_genome)
            hits["source_contig"] = six[0].fillna("")
            hits["seq_from"] = (
                pd.to_numeric(six[1], errors="coerce")
                .fillna(pd.to_numeric(prod[0], errors="coerce"))
                .fillna(0)
                .astype(int)
            )
            hits["seq_to"] = (
                pd.to_numeric(six[2], errors="coerce")
                .fillna(pd.to_numeric(prod[1], errors="coerce"))
                .fillna(0)
                .astype(int)
            )
            hits["strand"] = six[3].fillna(prod[2].map({"1": "+", "-1": "-"})).fillna("")
            hits = hits.sort_values(["bit_score", "evalue"], ascending=[False, True])
        hits.to_csv(self.results / "hits_main.tsv", sep="\t", index=False)
        if not hits.empty:
            best = (
                hits.sort_values(["genome_id", "bit_score"], ascending=[True, False])
                .drop_duplicates("genome_id")
                .copy()
            )
        else:
            best = hits.copy()
        best.to_csv(self.results / "hits_best_per_genome.tsv", sep="\t", index=False)
        return hits

    def downstream_analysis(self, hits: pd.DataFrame, trimmed_alignment: Path) -> None:
        if not hits.empty:
            hmm_length = self.hmm_length()
            classified = classify_hits(hits, hmm_length=hmm_length)
        else:
            classified = hits.copy()
        classified.to_csv(self.results / "hits_classified.tsv", sep="\t", index=False)

        matrix = build_matrix(hits) if not hits.empty else pd.DataFrame()
        matrix.to_csv(self.results / "presence_absence_matrix.tsv", sep="\t")
        if not matrix.empty:
            try:
                heatmap_png(matrix, self.figures / "presence_absence_heatmap.png")
            except Exception as exc:
                self.log(f"WARNING: heatmap failed: {exc}")

        tax = taxonomy_table(hits) if not hits.empty else pd.DataFrame()
        tax.to_csv(self.results / "taxonomy_table.tsv", sep="\t", index=False)

        seed = Path(self.manifest["core"]["input_path"])
        cluster = cluster_cdhit(seed, self.out / "clusters", identity=0.4, coverage=0.8, threads=self.args.cpu)
        if cluster.get("membership_df") is not None and not cluster["membership_df"].empty:
            cluster["membership_df"].to_csv(
                self.results / "cluster_membership.tsv", sep="\t", index=False
            )
            cluster_summary(cluster["membership_df"]).to_csv(
                self.results / "cluster_summary.tsv", sep="\t", index=False
            )
        if cluster.get("rep_faa") and Path(cluster["rep_faa"]).exists():
            shutil.copy2(cluster["rep_faa"], self.results / "cluster_reps.faa")

        meme = run_meme(seed, self.out / "motifs" / "meme", n_motifs=3, min_width=6, max_width=30, cpu=self.args.cpu)
        if meme.get("success"):
            fimo = run_fimo(Path(meme["meme_txt"]), seed, self.out / "motifs" / "fimo")
            if isinstance(fimo, pd.DataFrame):
                fimo.to_csv(self.results / "fimo_hits.tsv", sep="\t", index=False)

        phy = run_iqtree(trimmed_alignment, self.out / "trees", model="MFP", bootstrap=1000, cpu=self.args.cpu)
        if phy.get("success") and phy.get("treefile"):
            render_tree(Path(phy["treefile"]), hits, self.figures)

        nt_hits = hits[hits.get("db_type", pd.Series(dtype=str)).eq("nucleotide")].copy() if not hits.empty else hits
        syn_df, placement_df = build_synteny_table(
            nt_hits,
            flanks=5,
            max_genomes=self.args.max_synteny_genomes,
            sequence_cache_dir=self.results / "synteny_context_cache",
            log_callback=self.log,
        )
        syn_df.to_csv(self.results / "synteny_table.tsv", sep="\t", index=False)
        placement_df.to_csv(
            self.results / "synteny_placement_report.tsv", sep="\t", index=False
        )
        if not syn_df.empty:
            export_gff3(syn_df, self.results / "synteny_neighborhoods.gff3")
            try:
                export_synteny_figures(
                    syn_df,
                    conservation_scores(syn_df),
                    out_dir=self.figures,
                    flanks=5,
                    max_genomes=self.args.max_synteny_genomes,
                    dpi=300,
                    log_callback=self.log,
                )
            except Exception as exc:
                self.log(f"WARNING: synteny figure export failed: {exc}")
            try:
                build_neighborhood_genbanks(syn_df, self.results / "synteny_genbanks")
            except Exception as exc:
                self.log(f"WARNING: synteny GenBank export failed: {exc}")

    def hmm_length(self) -> int:
        hmm = self.hmm_dir / "benchmark_profile.hmm"
        try:
            with hmm.open(errors="replace") as fh:
                for line in fh:
                    if line.startswith("LENG"):
                        return int(line.split()[1])
        except Exception:
            pass
        return 0

    def write_reports(self, hits: pd.DataFrame, metrics: list[dict]) -> None:
        write_tsv(metrics, self.metrics_path)
        summary_rows = []
        for row in metrics:
            files = row.get("files") or {}
            file_records = files.values() if isinstance(files, dict) else files
            urls = []
            sha256s = []
            accessed = []
            source_size = 0
            for info in file_records:
                if not isinstance(info, dict):
                    continue
                if info.get("url") and info.get("url") not in urls:
                    urls.append(info.get("url"))
                if info.get("sha256") and info.get("sha256") not in sha256s:
                    sha256s.append(info.get("sha256"))
                if info.get("accessed_at"):
                    accessed.append(str(info.get("accessed_at")))
                try:
                    source_size += int(info.get("size_bytes") or 0)
                except Exception:
                    pass
            summary_rows.append(
                {
                    "database": row.get("database"),
                    "status": row.get("status"),
                    "optional": row.get("optional", ""),
                    "file_count": row.get("file_count", 0),
                    "hit_count": row.get("hit_count", 0),
                    "strict_count": row.get("strict_count", 0),
                    "nt_orf_mode": row.get("nt_orf_mode", ""),
                    "runtime_seconds": round(float(row.get("runtime_seconds", 0)), 2),
                    "downloaded_bytes": row.get("downloaded_bytes", 0),
                    "source_url_count": len(urls),
                    "source_urls": ";".join(urls[:10]),
                    "source_size_bytes": source_size,
                    "source_accessed_first": min(accessed) if accessed else "",
                    "source_accessed_last": max(accessed) if accessed else "",
                    "source_sha256_count": len(sha256s),
                    "source_sha256_prefixes": ";".join(value[:16] for value in sha256s[:10]),
                    "error": row.get("error", ""),
                }
            )
        write_tsv(summary_rows, self.summary_path)

        tools = {
            name: {"available": bool(find_tool(name)), "version": "", "description": name}
            for name in [
                "hmmbuild",
                "hmmsearch",
                "hmmscan",
                "mafft",
                "trimal",
                "prodigal",
                "seqkit",
                "mmseqs",
                "meme",
                "fimo",
                "iqtree",
            ]
        }
        state = {
            "input": {"params": self.manifest.get("core", {}).get("input_summary", {})},
            "benchmark": self.manifest,
        }
        repro = build_reproducibility_json(self.out, hits, state, tools)
        methods = generate_methods_text(self.out, repro)
        context = build_report_context(self.out, hits, repro, methods, tools)
        report = render_html_report(self.out, context)
        report2 = self.reports / "all_database_benchmark_report.html"
        shutil.copy2(report, report2)
        zip_path = create_export_zip(self.out)
        verdict = self.deploy_readiness(metrics)
        (self.reports / "DEPLOY_READINESS_VERDICT.txt").write_text(verdict + "\n")
        self.manifest["status"] = "complete"
        self.manifest["final_report"] = str(report2)
        self.manifest["export_zip"] = str(zip_path)
        self.manifest["deploy_readiness_verdict"] = verdict
        self.save_manifest()

    def deploy_readiness(self, metrics: list[dict]) -> str:
        selected = set(self.manifest.get("selected_databases", []))
        all_registered = {db["name"] for db in BUILTIN_DATABASES}
        scope = "all registered databases" if selected == all_registered else "selected databases"
        required_failures = [
            m for m in metrics if m.get("status") != "complete" and not m.get("optional")
        ]
        optional_failures = [
            m for m in metrics if m.get("status") != "complete" and m.get("optional")
        ]
        if required_failures:
            return (
                "NOT READY: required databases failed: "
                + ", ".join(m.get("database", "") for m in required_failures)
            )
        if optional_failures:
            return (
                f"READY WITH OPTIONAL ANNOTATION WARNINGS: required/core discovery {scope} passed; optional databases failed or were skipped: "
                + ", ".join(m.get("database", "") for m in optional_failures)
                + ". VOGDB VFAM is the preferred viral ortholog annotation layer."
            )
        return f"READY: {scope} completed successfully."

    def run(self) -> None:
        self.init_dirs()
        self.manifest["status"] = "running"
        self.manifest["active_command"] = {
            "preset": self.args.preset,
            "databases": self.args.databases,
            "dry_run": self.args.dry_run,
            "force": self.args.force,
            "keep_cache": self.args.keep_cache,
            "cpu": self.args.cpu,
            "min_free_gb": self.args.min_free_gb,
            "nt_orf_mode": self.args.nt_orf_mode,
            "min_orf_aa": self.args.min_orf_aa,
            "fasta": str(self.fasta),
            "out": str(self.out),
        }
        self.save_manifest()
        self.log(
            "Benchmark runner starting: "
            f"preset={self.args.preset}, out={self.out}, fasta={self.fasta}"
        )
        self.preflight()
        dbs = self.selected_databases()
        self.manifest["selected_databases"] = [db["name"] for db in dbs]
        self.save_manifest()
        self.log(f"Selected {len(dbs)} databases: {', '.join(self.manifest['selected_databases'])}")
        if self.args.dry_run:
            self.dry_run()
            return

        hmm, trimmed, _ = self.build_core()
        metrics = []
        for db in dbs:
            try:
                result = self.run_database(db, hmm)
                metrics.append(result)
                write_tsv(metrics, self.metrics_path)
            except Exception as exc:
                result = {
                    "database": db["name"],
                    "status": "failed",
                    "optional": db.get("optional", False),
                    "type": db.get("type"),
                    "search_mode": db.get("search_mode"),
                    "error": str(exc),
                    "file_count": 0,
                    "hit_count": 0,
                    "strict_count": 0,
                    "runtime_seconds": 0,
                    "downloaded_bytes": 0,
                }
                self.log(f"ERROR {db['name']}: {exc}")
                self.manifest["databases"][db["name"]] = result
                self.manifest.setdefault("failures", []).append(result)
                self.save_manifest()
                metrics.append(result)
                write_tsv(metrics, self.metrics_path)
                if not db.get("optional", False) or self.args.stop_on_optional_failure:
                    raise

        hits = self.collect_hits()
        self.downstream_analysis(hits, trimmed)
        self.write_reports(hits, metrics)
        self.log(self.manifest["deploy_readiness_verdict"])


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", default=str(DEFAULT_FASTA), help="Input seed FASTA")
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="Benchmark output root")
    parser.add_argument(
        "--preset",
        choices=["all", "partial", "smoke"],
        default="all",
        help="Database preset. Ignored if --databases is set.",
    )
    parser.add_argument("--databases", default="", help="Comma-separated database names")
    parser.add_argument("--dry-run", action="store_true", help="Expand DB URLs only")
    parser.add_argument("--force", action="store_true", help="Re-run completed steps")
    parser.add_argument("--keep-cache", action="store_true", help="Keep downloaded DB cache files")
    parser.add_argument(
        "--keep-translation-chunks",
        action="store_true",
        help="Keep seqkit/translation chunk working directories",
    )
    parser.add_argument("--cpu", type=int, default=4)
    parser.add_argument(
        "--nt-orf-mode",
        choices=["sixframe", "prodigal"],
        default="sixframe",
        help=(
            "How nucleotide databases are translated before hmmsearch. "
            "sixframe is exhaustive and best for annotation-missed genes; "
            "prodigal is a faster conventional annotation baseline."
        ),
    )
    parser.add_argument(
        "--min-orf-aa",
        type=int,
        default=30,
        help="Minimum amino-acid length for six-frame ORFs",
    )
    parser.add_argument("--evalue", default="1e-5")
    parser.add_argument("--strict-bits", type=float, default=45.0)
    parser.add_argument("--min-recovery", type=float, default=0.95)
    parser.add_argument("--min-free-gb", type=float, default=20.0)
    parser.add_argument("--max-synteny-genomes", type=int, default=24)
    parser.add_argument("--download-timeout", type=int, default=7200)
    parser.add_argument("--search-timeout", type=int, default=86400)
    parser.add_argument(
        "--stop-on-optional-failure",
        action="store_true",
        help="Treat optional DB failures as fatal",
    )
    return parser.parse_args(list(argv))


def main(argv: Iterable[str] = sys.argv[1:]) -> int:
    args = parse_args(argv)
    bench = Benchmark(args)
    try:
        bench.run()
        return 0
    except Exception as exc:
        bench.log(f"FATAL: {exc}")
        bench.manifest["status"] = "failed"
        bench.manifest["fatal_error"] = str(exc)
        bench.save_manifest()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
