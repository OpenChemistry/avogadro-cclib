# SPDX-License-Identifier: BSD-3-Clause

"""Transparent-decompression tests.

avogadro-cclib sniffs magic bytes on the raw input (in both string and file
input modes) to accept gzip/bzip2/xz/zip/zstd-compressed output files
without Avogadro needing to know about compression file extensions at all
-- see avogadro_cclib._decompress. All compressed variants are built in
memory / in tmp_path from the existing tests/data/ fixtures; no new fixture
files are added.
"""

import bz2
import gzip
import io
import json
import lzma
import zipfile

import pytest

from avogadro_cclib import _decompress, _ZSTD_MAGIC

from conftest import READABLE_FIXTURE_NAMES, read_fixture_text, run_cli

_SAMPLE_NAME = READABLE_FIXTURE_NAMES[0]


def _zstd_decoder():
    """Return a zstd-compress callable, or None if no encoder is importable.

    We only need to *compress* for the test; avogadro_cclib._decompress_zstd
    is what performs the decode under test.
    """
    try:
        from compression import zstd

        return zstd.compress
    except ImportError:
        pass
    try:
        import zstandard
    except ImportError:
        return None
    return zstandard.ZstdCompressor().compress


@pytest.fixture(scope="module")
def sample_bytes():
    return read_fixture_text(_SAMPLE_NAME).encode("utf-8")


@pytest.fixture(scope="module")
def baseline_cjson_bytes(sample_bytes):
    result = run_cli(["cclib", "--read"], sample_bytes)
    assert result.returncode == 0, result.stderr.decode("utf-8", errors="replace")
    return result.stdout


# -- unit tests for _decompress directly -------------------------------------


def test_decompress_passes_through_plain_text():
    text = b"just some plain ascii text, not compressed at all\n"
    assert _decompress(text) == text


def test_decompress_gzip():
    payload = b"hello gzip world"
    assert _decompress(gzip.compress(payload)) == payload


def test_decompress_bzip2():
    payload = b"hello bzip2 world"
    assert _decompress(bz2.compress(payload)) == payload


def test_decompress_xz():
    payload = b"hello xz world"
    assert _decompress(lzma.compress(payload)) == payload


def test_decompress_zip_single_member():
    payload = b"hello zip world"
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as archive:
        archive.writestr("only-member.txt", payload)
    assert _decompress(buf.getvalue()) == payload


def test_decompress_zip_multi_member_raises_clear_error():
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as archive:
        archive.writestr("a.txt", b"a")
        archive.writestr("b.txt", b"b")
    with pytest.raises(ValueError, match="exactly one member"):
        _decompress(buf.getvalue())


def test_decompress_zstd():
    compress = _zstd_decoder()
    payload = b"hello zstd world"
    if compress is None:
        # Neither the stdlib compression.zstd (3.14+) nor the third-party
        # zstandard package is importable here -- a recognised-but-unusable
        # signature must give a clear, actionable error, never a garbled
        # parse.
        with pytest.raises(ValueError, match="zstandard"):
            _decompress(_ZSTD_MAGIC + b"\x00" * 8)
    else:
        compressed = compress(payload)
        assert compressed.startswith(_ZSTD_MAGIC)
        assert _decompress(compressed) == payload


# -- end-to-end CLI tests, parametrized over compression scheme --------------


@pytest.mark.parametrize(
    "compress",
    [gzip.compress, bz2.compress, lzma.compress],
    ids=["gzip", "bzip2", "xz"],
)
def test_string_mode_compressed_matches_uncompressed(
    compress, sample_bytes, baseline_cjson_bytes
):
    result = run_cli(["cclib", "--read"], compress(sample_bytes))
    assert result.returncode == 0, result.stderr.decode("utf-8", errors="replace")
    assert result.stdout == baseline_cjson_bytes


@pytest.mark.parametrize(
    "compress,suffix",
    [(gzip.compress, ".gz"), (bz2.compress, ".bz2"), (lzma.compress, ".xz")],
    ids=["gzip", "bzip2", "xz"],
)
def test_file_mode_compressed_matches_uncompressed(
    compress, suffix, sample_bytes, baseline_cjson_bytes, tmp_path
):
    compressed_path = tmp_path / f"sample{suffix}"
    compressed_path.write_bytes(compress(sample_bytes))

    payload = json.dumps(
        {"operation": "read", "filename": str(compressed_path)}
    ).encode("utf-8")
    result = run_cli(["cclib", "--read"], payload)
    assert result.returncode == 0, result.stderr.decode("utf-8", errors="replace")
    assert result.stdout == baseline_cjson_bytes


def test_file_mode_gzip_with_no_gz_extension_is_still_sniffed(
    sample_bytes, baseline_cjson_bytes, tmp_path
):
    # The whole point of magic-byte sniffing over cclib's own
    # extension-dispatched decompression: it works even when the filename
    # gives no hint that the content is compressed.
    compressed_path = tmp_path / "weird_name_with_no_compression_suffix"
    compressed_path.write_bytes(gzip.compress(sample_bytes))

    payload = json.dumps(
        {"operation": "read", "filename": str(compressed_path)}
    ).encode("utf-8")
    result = run_cli(["cclib", "--read"], payload)
    assert result.returncode == 0, result.stderr.decode("utf-8", errors="replace")
    assert result.stdout == baseline_cjson_bytes


def test_string_mode_single_member_zip_matches_uncompressed(
    sample_bytes, baseline_cjson_bytes
):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as archive:
        archive.writestr("output.out", sample_bytes)
    result = run_cli(["cclib", "--read"], buf.getvalue())
    assert result.returncode == 0, result.stderr.decode("utf-8", errors="replace")
    assert result.stdout == baseline_cjson_bytes


def test_string_mode_multi_member_zip_fails_cleanly(sample_bytes):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as archive:
        archive.writestr("a.out", sample_bytes)
        archive.writestr("b.out", sample_bytes)
    result = run_cli(["cclib", "--read"], buf.getvalue())
    assert result.returncode != 0
    assert result.stdout == b""
    assert b"exactly one member" in result.stderr


def test_uncompressed_input_is_unaffected(sample_bytes, baseline_cjson_bytes):
    result = run_cli(["cclib", "--read"], sample_bytes)
    assert result.returncode == 0
    assert result.stdout == baseline_cjson_bytes
