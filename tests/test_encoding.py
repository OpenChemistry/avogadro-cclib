# SPDX-License-Identifier: BSD-3-Clause

"""Encoding robustness tests for reader.read().

Real output files come from many platforms: Windows-authored files are
CRLF, some carry a UTF-8 BOM (again, typically from Windows tools), and
comments or basis-set/element labels occasionally contain non-ASCII
characters (e.g. "Å"). None of that should ever crash the reader.
"""

from conftest import READABLE_FIXTURE_NAMES, quiet_read, read_fixture_text


def test_lf_and_crlf_are_equivalent():
    name = READABLE_FIXTURE_NAMES[0]
    text = read_fixture_text(name)
    lf_text = text.replace("\r\n", "\n")
    crlf_text = lf_text.replace("\n", "\r\n")

    assert quiet_read(lf_text) == quiet_read(crlf_text)


def test_utf8_bom_does_not_change_result():
    name = READABLE_FIXTURE_NAMES[0]
    text = read_fixture_text(name)

    plain = quiet_read(text)
    with_bom = quiet_read("﻿" + text)

    assert plain == with_bom


def test_non_ascii_content_does_not_raise():
    name = READABLE_FIXTURE_NAMES[0]
    text = read_fixture_text(name)
    injected = text + "\n! stray non-ascii note: r = 1.0 Å\n"

    cjson = quiet_read(injected)
    assert cjson["chemicalJson"] == 1
