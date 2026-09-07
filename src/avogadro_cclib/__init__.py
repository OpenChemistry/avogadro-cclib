"""
/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the 3-Clause BSD License, (see "LICENSE").
******************************************************************************/
"""

import argparse
import gzip
import io
import json
import sys
import traceback
import zipfile

try:
    import bz2
except ImportError:  # pragma: no cover - some CPython builds omit bz2
    bz2 = None

try:
    import lzma
except ImportError:  # pragma: no cover - some CPython builds omit lzma
    lzma = None

# Magic-byte signatures for the compressed formats we sniff. We deliberately
# do this ourselves rather than relying on cclib's own filename-suffix-based
# decompression (cclib/parser/logfilewrapper.py): our plugin never hands
# cclib a filename (it always reads the bytes itself and hands cclib a
# StringIO), and sniffing content instead of the extension lets Avogadro
# open e.g. a compressed file with no matching suffix at all.
_GZIP_MAGIC = b"\x1f\x8b"
_BZIP2_MAGIC = b"BZh"
_XZ_MAGIC = b"\xfd7zXZ\x00"
_ZIP_MAGIC = b"PK\x03\x04"
_ZSTD_MAGIC = b"\x28\xb5\x2f\xfd"


def _decompress(raw, debug=False):
    """Transparently decompress ``raw`` if it looks compressed.

    Returns ``raw`` unchanged when no known signature matches. This is
    content sniffing only -- no file extension is consulted or required.
    """
    if raw.startswith(_GZIP_MAGIC):
        if debug:
            print("avogadro-cclib: detected gzip-compressed input", file=sys.stderr)
        return gzip.decompress(raw)

    if raw.startswith(_BZIP2_MAGIC):
        if bz2 is None:
            raise ValueError(
                "input looks bzip2-compressed, but this Python build has no "
                "'bz2' module"
            )
        if debug:
            print("avogadro-cclib: detected bzip2-compressed input", file=sys.stderr)
        return bz2.decompress(raw)

    if raw.startswith(_XZ_MAGIC):
        if lzma is None:
            raise ValueError(
                "input looks xz-compressed, but this Python build has no "
                "'lzma' module"
            )
        if debug:
            print("avogadro-cclib: detected xz-compressed input", file=sys.stderr)
        return lzma.decompress(raw)

    if raw.startswith(_ZIP_MAGIC):
        if debug:
            print("avogadro-cclib: detected zip-compressed input", file=sys.stderr)
        with zipfile.ZipFile(io.BytesIO(raw)) as archive:
            names = archive.namelist()
            if len(names) != 1:
                raise ValueError(
                    "zip input must contain exactly one member to be read "
                    f"unambiguously; found {len(names)}: {names}"
                )
            return archive.read(names[0])

    if raw.startswith(_ZSTD_MAGIC):
        if debug:
            print("avogadro-cclib: detected zstd-compressed input", file=sys.stderr)
        return _decompress_zstd(raw)

    return raw


def _decompress_zstd(raw):
    """Decompress zstd-compressed bytes without adding a hard dependency.

    Python's stdlib gained ``compression.zstd`` in 3.14; this project
    supports 3.11+, so we try that first, then fall back to the
    third-party ``zstandard`` package if it happens to be installed.
    Neither is required by pyproject.toml -- a clear error tells the user
    what to install if they hit a zstd file without either available.
    """
    try:
        from compression import zstd

        return zstd.decompress(raw)
    except ImportError:
        pass

    try:
        import zstandard
    except ImportError as exc:
        raise ValueError(
            "input looks zstd-compressed, but no zstd decoder is available "
            "in this Python environment; install the 'zstandard' package "
            "(pip install zstandard) to read zstd-compressed files"
        ) from exc

    return zstandard.ZstdDecompressor().decompress(raw)


def main():
    # Avogadro calls the plugin as:
    #   avogadro-cclib <identifier> [--read] [--write] [--lang <locale>] [--debug]
    #
    # Two input modes are possible, and we auto-detect which one we got:
    #   - "string" mode (default, used by all released Avogadro builds): the
    #     raw bytes of the computational chemistry output file are written to
    #     stdin, with no wrapper of any kind.
    #   - "file" mode (opt-in via `input-mode = "file"` in pyproject.toml, only
    #     understood by recent Avogadro builds): stdin carries a small JSON
    #     object, e.g. {"operation": "read", "filename": "/path/to/file"}.
    #
    # In both modes, stdout must be nothing but the raw CJSON document -- no
    # envelope, no extra text -- since it is passed directly to Avogadro's
    # CjsonFormat::readString().
    parser = argparse.ArgumentParser()
    parser.add_argument("feature")
    parser.add_argument("--lang", nargs="?", default="en")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--read", action="store_true")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()

    # Read stdin as bytes -- never let the platform locale decide the
    # encoding (this is the source of mojibake / UnicodeDecodeError on
    # Windows, where stdin would otherwise be decoded as cp1252 or cp932).
    raw = sys.stdin.buffer.read()

    output = None
    try:
        # Sniff for the file-mode JSON envelope on the raw bytes decoded as
        # UTF-8 text. A compressed or otherwise binary payload will simply
        # fail this parse (or fail the dict/key check) and fall through to
        # string mode, which is what we want -- this decoded text is used
        # only for that detection, never as the actual file content.
        mode = "string"
        payload = None
        try:
            payload = json.loads(raw.decode("utf-8-sig", errors="replace"))
        except (json.JSONDecodeError, ValueError):
            payload = None

        if isinstance(payload, dict) and ("filename" in payload or "file" in payload):
            mode = "file"
            filename = payload.get("filename", payload.get("file"))
            with open(filename, "rb") as handle:
                raw_content = handle.read()
        else:
            raw_content = raw

        # Transparently decompress gzip/bzip2/xz/zip/zstd content, sniffed
        # by magic bytes rather than by file extension -- we never register
        # compression extensions with Avogadro (see pyproject.toml).
        raw_content = _decompress(raw_content, debug=args.debug)
        content = raw_content.decode("utf-8-sig", errors="replace")

        if args.debug:
            print(
                f"avogadro-cclib: feature={args.feature} mode={mode} "
                f"input_bytes={len(raw_content)}",
                file=sys.stderr,
            )

        match args.feature:
            case "cclib":
                if args.read:
                    from .reader import read

                    cjson = read(content)
                    output = json.dumps(cjson, ensure_ascii=True)
    except Exception as exc:  # noqa: BLE001 - report cleanly, never crash
        print(f"avogadro-cclib: error reading input: {exc}", file=sys.stderr)
        if args.debug:
            traceback.print_exc(file=sys.stderr)
        sys.exit(1)

    if output is None:
        print(
            f"avogadro-cclib: no output produced for feature='{args.feature}'",
            file=sys.stderr,
        )
        sys.exit(1)

    sys.stdout.buffer.write(output.encode("utf-8"))
    sys.stdout.buffer.flush()
