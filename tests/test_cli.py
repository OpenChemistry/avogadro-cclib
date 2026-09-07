# SPDX-License-Identifier: BSD-3-Clause

"""End-to-end CLI tests: invoke the plugin exactly the way Avogadro does.

Avogadro spawns ``avogadro-cclib <identifier> --read [--debug] --lang <locale>``
as a subprocess and writes to its stdin. We reproduce that here with
``sys.executable -c "from avogadro_cclib import main; main()" ...`` so the
tests work the same way on Windows (where the installed console-script is a
``.exe`` wrapper) as everywhere else.
"""

import json
import os

import pytest

from conftest import (
    ALL_FIXTURE_PARAMS,
    READABLE_FIXTURE_NAMES,
    read_fixture_text,
    run_cli,
    strict_json_loads,
)


@pytest.mark.parametrize("name", ALL_FIXTURE_PARAMS)
def test_string_mode_read(name):
    text = read_fixture_text(name)
    result = run_cli(["cclib", "--read"], text.encode("utf-8"))
    assert result.returncode == 0, result.stderr.decode("utf-8", errors="replace")
    cjson = strict_json_loads(result.stdout.decode("utf-8"))
    assert cjson["chemicalJson"] == 1


@pytest.mark.parametrize("name", ALL_FIXTURE_PARAMS)
def test_file_mode_matches_string_mode(name, fixture_paths):
    text = read_fixture_text(name)
    string_result = run_cli(["cclib", "--read"], text.encode("utf-8"))
    assert string_result.returncode == 0
    string_cjson = strict_json_loads(string_result.stdout.decode("utf-8"))

    payload = json.dumps(
        {"operation": "read", "filename": str(fixture_paths[name])}
    ).encode("utf-8")
    file_result = run_cli(["cclib", "--read"], payload)
    assert file_result.returncode == 0, file_result.stderr.decode(
        "utf-8", errors="replace"
    )
    file_cjson = strict_json_loads(file_result.stdout.decode("utf-8"))

    assert file_cjson == string_cjson


def test_stdout_is_bare_cjson_with_no_leaked_log_lines():
    # Guards against cclib's own logging (or anything else) leaking into
    # stdout: stdout must be exactly one JSON document, nothing more.
    name = READABLE_FIXTURE_NAMES[0]
    text = read_fixture_text(name)
    result = run_cli(["cclib", "--read"], text.encode("utf-8"))
    assert result.returncode == 0
    stdout_text = result.stdout.decode("utf-8")
    assert stdout_text.startswith("{")
    decoder = json.JSONDecoder()
    obj, end = decoder.raw_decode(stdout_text)
    assert (
        stdout_text[end:].strip() == ""
    ), "trailing content found after the JSON document"
    assert obj["chemicalJson"] == 1


def test_failure_path_on_garbage_input():
    result = run_cli(["cclib", "--read"], b"not a chemistry file\n")
    assert result.returncode != 0
    assert result.stdout == b""
    assert result.stderr.strip() != b""


def test_windows_locale_simulation_with_non_ascii_content():
    # Regression test for the Windows encoding bug: a naive `json.load(sys.stdin)`
    # or text-mode stdin read lets the platform locale (e.g. cp1252 on
    # Windows) decide the input encoding, which mangles or crashes on
    # non-ASCII bytes. Force a hostile locale via env vars and inject a
    # UTF-8 BOM plus a non-ASCII character into otherwise-valid content.
    name = READABLE_FIXTURE_NAMES[0]
    text = read_fixture_text(name)
    injected = "﻿" + text + "\n! stray non-ascii note: r = 1.0 Å\n"

    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "cp1252"
    env["PYTHONUTF8"] = "0"
    env["LC_ALL"] = "C"

    result = run_cli(["cclib", "--read"], injected.encode("utf-8"), env=env)
    assert result.returncode == 0, result.stderr.decode("utf-8", errors="replace")
    cjson = strict_json_loads(result.stdout.decode("utf-8"))
    assert cjson["chemicalJson"] == 1
