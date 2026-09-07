# SPDX-License-Identifier: BSD-3-Clause

"""One test per bug fixed, each named for the symptom it guards against."""

import json

from avogadro_cclib.reader import read

from conftest import read_cjson, read_fixture_text, run_cli


def test_more_scf_energies_than_geometries():
    # Symptom: nwchem_He_scs-mp2 has 18 scfenergies but only 1 atomcoords.
    # The old code assumed one geometry per SCF energy and indexed
    # data.atomcoords[i] for i up to 17, raising:
    #   IndexError: index 1 is out of bounds for axis 0 with size 1
    cjson = read_cjson("nwchem_He_scs-mp2.out")

    coords = cjson["atoms"]["coords"]
    if "3dSets" in coords:
        assert "properties" in cjson and "energies" in cjson["properties"]
        assert len(coords["3dSets"]) == len(cjson["properties"]["energies"])

    assert "properties" in cjson
    assert "totalEnergy" in cjson["properties"]


def test_missing_coordinates_raises_clean_error():
    # Symptom: nwchem_He_b2plyp has no atomcoords/atomnos at all -- cclib
    # parses only metadata. The old code accessed data.atomcoords[-1]
    # unconditionally, raising:
    #   AttributeError: 'ccData_optdone_bool' object has no attribute 'atomnos'
    # reader.read() must raise a clean ValueError instead.
    text = read_fixture_text("nwchem_He_b2plyp.out")
    try:
        read(text)
        assert False, "expected ValueError, no exception was raised"
    except ValueError:
        pass
    except AttributeError as exc:  # pragma: no cover - the bug we're guarding against
        raise AssertionError(
            f"reader.read() raised AttributeError instead of a clean ValueError: {exc}"
        )

    # And the CLI must surface this as a clean failure: exit 1, a message on
    # stderr, and nothing on stdout.
    result = run_cli(["cclib", "--read"], text.encode("utf-8"))
    assert result.returncode != 0
    assert result.stdout == b""
    assert result.stderr.strip() != b""


def test_envelope_free_stdout():
    # Symptom (commit 69927f5): the CLI used to wrap the CJSON document in
    # {"cjson": ...} and print() it with the platform's text-mode codec.
    # Avogadro's C++ side passes stdout straight to CjsonFormat::readString(),
    # so the output must be the bare CJSON object, with no envelope.
    text = read_fixture_text("nwchem_tce_ccsd_t_h2o.out")
    result = run_cli(["cclib", "--read"], text.encode("utf-8"))
    assert result.returncode == 0, result.stderr.decode("utf-8", errors="replace")

    obj = json.loads(result.stdout.decode("utf-8"))
    assert "cjson" not in obj
    assert obj["chemicalJson"] == 1


def test_scalar_non_finite_value_is_dropped():
    # Symptom: json.dumps emits a bare NaN token for a non-finite scalar such
    # as properties.totalEnergy. nlohmann::json discards the *entire*
    # document on that token, so Avogadro reports only "Root is discarded"
    # and nothing loads. Scalars must be dropped like lists are.
    from avogadro_cclib.reader import _sanitize

    cleaned = _sanitize({"properties": {"totalEnergy": float("nan"), "totalCharge": 0}})
    assert "totalEnergy" not in cleaned["properties"]
    assert cleaned["properties"]["totalCharge"] == 0
    json.dumps(cleaned, allow_nan=False)


def test_unidentifiable_format_reports_cleanly():
    # Symptom: ccopen() returns None for a file it cannot attribute to any
    # program, and the old code called .parse() on it, surfacing
    # "'NoneType' object has no attribute 'parse'" to the user.
    try:
        read("not a chemistry file\n")
        assert False, "expected ValueError, no exception was raised"
    except ValueError as exc:
        assert "identify" in str(exc)
    except AttributeError as exc:  # pragma: no cover - the bug we guard against
        raise AssertionError(f"reader.read() raised AttributeError: {exc}")
