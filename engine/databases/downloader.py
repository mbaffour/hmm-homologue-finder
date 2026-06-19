"""
databases/downloader.py — Resumable async database download with progress tracking.
"""

import asyncio
import gzip
import logging
import shutil
import time
import urllib.error
import urllib.request
from fnmatch import fnmatch
from html.parser import HTMLParser
from pathlib import Path
from typing import Callable, Optional

import aiofiles
import aiohttp

logger = logging.getLogger(__name__)

# ------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------

_CHUNK_SIZE = 1 << 20          # 1 MiB read/write chunks
_PROGRESS_INTERVAL = 0.5       # seconds between progress_callback calls
_MAX_RETRIES = 3
_RETRY_DELAY = 5.0             # seconds between retries
_CONNECT_TIMEOUT = 30          # seconds
_READ_TIMEOUT = 60             # seconds


# ------------------------------------------------------------------
# Public helpers
# ------------------------------------------------------------------

def format_bytes(n: int) -> str:
    """
    Return a human-readable representation of *n* bytes.

    Examples
    --------
    >>> format_bytes(1_500_000_000)
    '1.4 GB'
    >>> format_bytes(734_003_200)
    '699.9 MB'
    >>> format_bytes(512)
    '512 B'
    """
    if n < 0:
        return "? B"
    for unit, threshold in (("GB", 1 << 30), ("MB", 1 << 20), ("KB", 1 << 10)):
        if n >= threshold:
            return f"{n / threshold:.1f} {unit}"
    return f"{n} B"


def estimate_download_time(size_bytes: int, speed_mbps: float) -> str:
    """
    Return a human-readable ETA string.

    Parameters
    ----------
    size_bytes:
        Remaining bytes to download.
    speed_mbps:
        Current download speed in **megabytes per second** (MB/s, not Mbit/s).

    Returns
    -------
    str
        e.g. ``"~3 minutes"``, ``"~2 hours"``, ``"< 1 minute"``.
    """
    if speed_mbps <= 0 or size_bytes <= 0:
        return "unknown"

    seconds = size_bytes / (speed_mbps * 1_000_000)

    if seconds < 60:
        return "< 1 minute"
    if seconds < 3600:
        minutes = round(seconds / 60)
        return f"~{minutes} minute{'s' if minutes != 1 else ''}"
    hours = seconds / 3600
    if hours < 48:
        h = round(hours)
        return f"~{h} hour{'s' if h != 1 else ''}"
    days = round(hours / 24)
    return f"~{days} day{'s' if days != 1 else ''}"


# ------------------------------------------------------------------
# NCBI FTP directory listing
# ------------------------------------------------------------------

class _HrefParser(HTMLParser):
    """Minimal HTML parser that collects all href attribute values."""

    def __init__(self) -> None:
        super().__init__()
        self.hrefs: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, Optional[str]]]) -> None:
        if tag == "a":
            for attr, value in attrs:
                if attr == "href" and value:
                    self.hrefs.append(value)


def check_ncbi_ftp_files(base_url: str, pattern: str) -> list[str]:
    """
    List files at *base_url* whose names match *pattern* (glob syntax).

    Fetches the FTP/HTTP directory index at *base_url* and returns the
    full URLs of all files whose basenames match *pattern*.

    Parameters
    ----------
    base_url:
        Directory URL, e.g.
        ``"https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/"``.
    pattern:
        Glob pattern applied to the filename only, e.g.
        ``"viral.*.protein.faa.gz"``.

    Returns
    -------
    list[str]
        Sorted list of matching full URLs.  Empty on error.
    """
    if not base_url.endswith("/"):
        base_url += "/"

    try:
        with urllib.request.urlopen(base_url, timeout=30) as resp:
            html = resp.read().decode("utf-8", errors="replace")
    except (urllib.error.URLError, OSError) as exc:
        logger.error("Could not fetch FTP directory %s: %s", base_url, exc)
        return []

    parser = _HrefParser()
    parser.feed(html)

    matched: list[str] = []
    for href in parser.hrefs:
        # Strip trailing slashes and relative prefixes
        filename = href.rstrip("/").split("/")[-1]
        if fnmatch(filename, pattern):
            # Build absolute URL
            if href.startswith("http://") or href.startswith("https://"):
                full_url = href
            else:
                full_url = base_url + filename
            matched.append(full_url)

    matched.sort()
    logger.debug(
        "check_ncbi_ftp_files: found %d files matching %r at %s",
        len(matched), pattern, base_url,
    )
    return matched


# ------------------------------------------------------------------
# Core download logic
# ------------------------------------------------------------------

async def _download_single(
    url: str,
    dest_path: Path,
    progress_callback: Optional[Callable[[int, int, float], None]],
    session: aiohttp.ClientSession,
) -> bool:
    """
    Download *url* to *dest_path*, resuming if the file already exists.

    Returns ``True`` on success.
    """
    dest_path.parent.mkdir(parents=True, exist_ok=True)

    # Determine bytes already on disk (for resume)
    existing_size = dest_path.stat().st_size if dest_path.exists() else 0
    headers: dict[str, str] = {}
    if existing_size > 0:
        headers["Range"] = f"bytes={existing_size}-"
        logger.debug("Attempting resume for %s at byte %d", dest_path.name, existing_size)

    timeout = aiohttp.ClientTimeout(connect=_CONNECT_TIMEOUT, sock_read=_READ_TIMEOUT)

    try:
        async with session.get(url, headers=headers, timeout=timeout) as resp:
            # Handle resume: 206 = partial content, 200 = full content
            if resp.status == 416:
                # Range not satisfiable — file already complete
                logger.info("File already complete: %s", dest_path.name)
                return True
            if resp.status == 200 and existing_size > 0:
                # Server doesn't support Range; start over
                logger.debug("Server ignored Range header; restarting download for %s", dest_path.name)
                existing_size = 0
            if resp.status not in (200, 206):
                logger.error("HTTP %d for %s", resp.status, url)
                return False

            total_from_server = int(resp.headers.get("Content-Length", 0))
            total_bytes = (existing_size + total_from_server) if total_from_server else 0
            downloaded = existing_size

            open_mode = "ab" if existing_size > 0 else "wb"
            last_progress = time.monotonic()
            speed_window_start = last_progress
            speed_window_bytes = 0

            try:
                async with aiofiles.open(dest_path, open_mode) as fh:
                    async for chunk in resp.content.iter_chunked(_CHUNK_SIZE):
                        await fh.write(chunk)
                        downloaded += len(chunk)
                        speed_window_bytes += len(chunk)

                        now = time.monotonic()
                        elapsed_window = now - speed_window_start
                        if elapsed_window >= _PROGRESS_INTERVAL:
                            speed_mbps = (
                                speed_window_bytes / elapsed_window / 1_000_000
                            )
                            if progress_callback is not None:
                                try:
                                    progress_callback(downloaded, total_bytes, speed_mbps)
                                except Exception:
                                    pass  # never let a callback crash the download
                            speed_window_start = now
                            speed_window_bytes = 0
                            last_progress = now

            except OSError as exc:
                logger.error("Disk error writing %s: %s", dest_path, exc)
                return False

    except aiohttp.ClientError as exc:
        logger.error("Network error downloading %s: %s", url, exc)
        return False

    return True


async def download_database(
    url: str,
    dest_path: Path,
    progress_callback: Optional[Callable[[int, int, float], None]] = None,
) -> bool:
    """
    Download a database from *url* to *dest_path* with resume support.

    Handles:

    * **Resume** — if *dest_path* already exists, a ``Range`` header is
      sent so only missing bytes are fetched.  Falls back to a full
      download if the server returns HTTP 200 instead of 206.
    * **Multi-file NCBI FTP patterns** — if *url* contains ``*``,
      :func:`check_ncbi_ftp_files` is called to expand the glob, and
      each file is downloaded into the same directory as *dest_path*.
    * **Retry** — up to :data:`_MAX_RETRIES` attempts with
      :data:`_RETRY_DELAY` second gaps on connection-level failures.
    * **Disk full** — ``OSError`` during write returns ``False``.
    * **Empty URL** — returns ``False`` immediately.

    Parameters
    ----------
    url:
        Source URL.  May contain ``*`` for NCBI FTP multi-file patterns.
    dest_path:
        Local destination file (or directory when *url* contains ``*``).
    progress_callback:
        Optional callable invoked roughly every 0.5 s with
        ``(downloaded_bytes: int, total_bytes: int, speed_mbps: float)``.

    Returns
    -------
    bool
        ``True`` if all files downloaded successfully.
    """
    if not url:
        logger.warning("download_database called with empty URL")
        return False

    # Multi-file mode: URL contains a glob wildcard
    if "*" in url:
        # Split into base directory and pattern
        last_slash = url.rfind("/")
        base_url = url[: last_slash + 1]
        pattern = url[last_slash + 1 :]
        file_urls = check_ncbi_ftp_files(base_url, pattern)
        if not file_urls:
            logger.error("No files matched pattern %r at %s", pattern, base_url)
            return False

        dest_dir = dest_path if dest_path.suffix == "" else dest_path.parent
        dest_dir.mkdir(parents=True, exist_ok=True)

        all_ok = True
        for file_url in file_urls:
            filename = file_url.split("/")[-1]
            file_dest = dest_dir / filename
            ok = await download_database(file_url, file_dest, progress_callback)
            if not ok:
                all_ok = False
        return all_ok

    # Single-file mode with retry
    connector = aiohttp.TCPConnector(limit=1)
    async with aiohttp.ClientSession(connector=connector) as session:
        for attempt in range(1, _MAX_RETRIES + 1):
            try:
                ok = await _download_single(url, dest_path, progress_callback, session)
                if ok:
                    return True
                logger.warning(
                    "Download attempt %d/%d failed for %s",
                    attempt, _MAX_RETRIES, url,
                )
            except (aiohttp.ServerConnectionError, asyncio.TimeoutError) as exc:
                logger.warning(
                    "Connection error on attempt %d/%d for %s: %s",
                    attempt, _MAX_RETRIES, url, exc,
                )
            if attempt < _MAX_RETRIES:
                logger.info("Retrying in %g s…", _RETRY_DELAY)
                await asyncio.sleep(_RETRY_DELAY)

    logger.error("All %d download attempts failed for %s", _MAX_RETRIES, url)
    return False


# ------------------------------------------------------------------
# Decompression helper
# ------------------------------------------------------------------

async def gunzip_file(gz_path: Path, out_path: Path) -> bool:
    """
    Decompress *gz_path* to *out_path* asynchronously.

    Reads the compressed file and writes decompressed content in
    :data:`_CHUNK_SIZE` byte chunks using :mod:`aiofiles` for non-blocking
    I/O.  The ``gzip`` module is used for decompression.

    Parameters
    ----------
    gz_path:
        Path to the ``.gz`` input file.
    out_path:
        Destination path for the decompressed output.

    Returns
    -------
    bool
        ``True`` on success, ``False`` on any error.
    """
    if not gz_path.exists():
        logger.error("gunzip_file: source file not found: %s", gz_path)
        return False

    out_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        # Decompress entirely in memory for files that fit; for very large
        # files the caller should consider a streaming approach.  We use a
        # thread pool to keep the event loop responsive.
        loop = asyncio.get_running_loop()

        def _decompress() -> bytes:
            with gzip.open(gz_path, "rb") as f:
                return f.read()

        decompressed = await loop.run_in_executor(None, _decompress)

        async with aiofiles.open(out_path, "wb") as out:
            # Write in chunks to avoid one giant awaitable
            for offset in range(0, len(decompressed), _CHUNK_SIZE):
                await out.write(decompressed[offset : offset + _CHUNK_SIZE])

        logger.info("Decompressed %s -> %s", gz_path.name, out_path.name)
        return True

    except (OSError, gzip.BadGzipFile, EOFError) as exc:
        logger.error("gunzip_file failed for %s: %s", gz_path, exc)
        # Remove partial output
        if out_path.exists():
            try:
                out_path.unlink()
            except OSError:
                pass
        return False


async def gunzip_file_streaming(gz_path: Path, out_path: Path) -> bool:
    """
    Streaming variant of :func:`gunzip_file` for very large files.

    Reads the compressed file in :data:`_CHUNK_SIZE` chunks and
    decompresses on the fly, writing output asynchronously.

    Parameters
    ----------
    gz_path:
        Path to the ``.gz`` input file.
    out_path:
        Destination path for the decompressed output.

    Returns
    -------
    bool
        ``True`` on success, ``False`` on any error.
    """
    if not gz_path.exists():
        logger.error("gunzip_file_streaming: source not found: %s", gz_path)
        return False

    out_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        async with aiofiles.open(out_path, "wb") as out_fh:
            with gzip.open(gz_path, "rb") as gz_fh:
                while True:
                    chunk = gz_fh.read(_CHUNK_SIZE)
                    if not chunk:
                        break
                    await out_fh.write(chunk)
        logger.info("Streaming-decompressed %s -> %s", gz_path.name, out_path.name)
        return True

    except (OSError, gzip.BadGzipFile, EOFError) as exc:
        logger.error("gunzip_file_streaming failed for %s: %s", gz_path, exc)
        if out_path.exists():
            try:
                out_path.unlink()
            except OSError:
                pass
        return False
