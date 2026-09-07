# SPDX-License-Identifier: BSD-3-Clause

"""Shared fixtures for the avogadro-cclib test suite.

The test corpus lives under ``tests/data/``. See ``tests/README.md`` for how
to extend it -- dropping a new (optionally gzip-compressed) real output file
in that directory is enough for it to be picked up automatically by the
corpus-wide contract tests, with no code changes required.
"""

import contextlib
import gzip
import io
import json
import shutil
import subprocess
import sys
import warnings
from pathlib import Path

import pytest

FIXTURE_DIR = Path(__file__).parent / "data"

# Avogadro spawns ``avogadro-cclib <identifier> --read [--debug] --lang
# <locale>`` as a subprocess and writes to its stdin. We reproduce that here
# with ``sys.executable -c "from avogadro_cclib import main; main()" ...`` so
# the tests work the same way on Windows (where the installed console-script
# is a ``.exe`` wrapper) as everywhere else.
MAIN_SNIPPET = "from avogadro_cclib import main; main()"

# Fixtures that cclib cannot turn into a readable molecule (e.g. no atomic
# coordinates at all), and so must be exercised by a targeted regression
# test rather than by the generic corpus contract tests.
EXPECTED_UNREADABLE = {
    # NWChem b2plyp job: cclib parses only calculation metadata -- there is
    # no atomcoords/atomnos data at all. reader.read() must raise a clean
    # ValueError instead of crashing with AttributeError. See
    # test_regressions.py::test_missing_coordinates_raises_clean_error.
    "nwchem_He_b2plyp.out",
}


def _iter_fixture_files():
    """Yield (logical_name, path) for every fixture under tests/data/."""
    for path in sorted(FIXTURE_DIR.iterdir()):
        if path.name.startswith("."):
            continue
        if not path.is_file():
            continue
        name = path.name[: -len(".gz")] if path.name.endswith(".gz") else path.name
        yield name, path


# Sorted logical names of every fixture under tests/data/.
FIXTURE_NAMES = sorted({name for name, _path in _iter_fixture_files()})
READABLE_FIXTURE_NAMES = [n for n in FIXTURE_NAMES if n not in EXPECTED_UNREADABLE]

if not READABLE_FIXTURE_NAMES:
    # Fail loudly rather than reporting a green run over an empty corpus.
    # This has bitten us once already: a global gitignore rule for "*.gz"
    # kept the fixtures out of the commit, so CI had nothing to test.
    raise RuntimeError(
        f"no readable test fixtures found in {FIXTURE_DIR} -- "
        "check that tests/data/ was committed"
    )


def fixture_params():
    """Parametrize values for the whole discovered corpus.

    Fixtures listed in EXPECTED_UNREADABLE are included (so a maintainer can
    see they exist) but marked to be skipped, rather than silently omitted
    or left to fail -- they get their own targeted test in
    test_regressions.py instead.
    """
    params = []
    for name in FIXTURE_NAMES:
        if name in EXPECTED_UNREADABLE:
            params.append(
                pytest.param(
                    name,
                    marks=pytest.mark.skip(
                        reason=(
                            "cclib cannot produce a readable molecule for this "
                            "fixture; see its dedicated test in "
                            "test_regressions.py"
                        )
                    ),
                )
            )
        else:
            params.append(pytest.param(name))
    return params


ALL_FIXTURE_PARAMS = fixture_params()


@pytest.fixture(scope="session")
def fixture_paths(tmp_path_factory):
    """Decompress every fixture once and return {name: real on-disk path}.

    Needed for tests that exercise the CLI's file-input mode, which takes a
    filename rather than file content.
    """
    tmp_dir = tmp_path_factory.mktemp("cclib-fixtures")
    paths = {}
    for name, path in _iter_fixture_files():
        dest = tmp_dir / name
        if path.name.endswith(".gz"):
            with gzip.open(path, "rb") as src, open(dest, "wb") as out:
                shutil.copyfileobj(src, out)
        else:
            shutil.copy(path, dest)
        paths[name] = dest
    return paths


def read_fixture_text(name):
    """Return the decoded text content of a fixture, by logical name."""
    for logical, path in _iter_fixture_files():
        if logical != name:
            continue
        if path.name.endswith(".gz"):
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
                return handle.read()
        return path.read_text(encoding="utf-8", errors="replace")
    raise FileNotFoundError(f"no fixture named {name!r} in {FIXTURE_DIR}")


def quiet_read(text):
    """Run reader.read(text), suppressing cclib's own logging chatter.

    Used both by read_cjson() (fixture files) and directly by tests that
    build the input text themselves (e.g. encoding variants).
    """
    from avogadro_cclib.reader import read

    buf = io.StringIO()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with contextlib.redirect_stderr(buf), contextlib.redirect_stdout(buf):
            return read(text)


def read_cjson(name):
    """Return the CJSON dict for a fixture, running the real reader.read().

    cclib's own logging chatter is suppressed so it doesn't clutter test
    output.
    """
    return quiet_read(read_fixture_text(name))


def run_cli(args, input_bytes, env=None):
    cmd = [sys.executable, "-c", MAIN_SNIPPET, *args]
    return subprocess.run(cmd, input=input_bytes, capture_output=True, env=env)


def strict_json_loads(text):
    """json.loads that rejects bare NaN/Infinity/-Infinity tokens."""

    def _reject(token):
        raise ValueError(f"non-finite constant in JSON output: {token}")

    return json.loads(text, parse_constant=_reject)
