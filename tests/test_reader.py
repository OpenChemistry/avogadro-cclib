# SPDX-License-Identifier: BSD-3-Clause

"""CJSON contract tests, parametrized over the whole discovered corpus.

For every readable fixture, every invariant below that applies must hold.
An assertion is skipped when the CJSON key it covers is absent -- but if the
key is present, the invariant must hold.
"""

import json

import pytest

from conftest import ALL_FIXTURE_PARAMS, read_cjson, strict_json_loads


@pytest.mark.parametrize("name", ALL_FIXTURE_PARAMS)
def test_cjson_contract(name):
    cjson = read_cjson(name)

    assert cjson["chemicalJson"] == 1

    n_atoms = len(cjson["atoms"]["elements"]["number"])
    coords = cjson["atoms"]["coords"]
    assert len(coords["3d"]) == 3 * n_atoms

    if "3dSets" in coords:
        for frame in coords["3dSets"]:
            assert len(frame) == 3 * n_atoms
        if "properties" in cjson and "energies" in cjson["properties"]:
            assert len(coords["3dSets"]) == len(cjson["properties"]["energies"])

    if "partialCharges" in cjson:
        for values in cjson["partialCharges"].values():
            assert len(values) == n_atoms

    if "vibrations" in cjson:
        vib = cjson["vibrations"]
        if "frequencies" in vib:
            n_modes = len(vib["frequencies"])
            for key in ("eigenVectors", "intensities", "modes"):
                if key in vib:
                    assert len(vib[key]) == n_modes
            if "eigenVectors" in vib:
                for eigvec in vib["eigenVectors"]:
                    assert len(eigvec) == 3 * n_atoms
            if "ramanIntensities" in vib:
                assert len(vib["ramanIntensities"]) == n_modes
            if "symmetries" in vib:
                assert len(vib["symmetries"]) == n_modes

    if "orbitals" in cjson:
        orb = cjson["orbitals"]
        if "occupations" in orb and "energies" in orb:
            assert len(orb["occupations"]) == len(orb["energies"])
        if "moCoefficients" in orb and "energies" in orb and len(orb["energies"]) > 0:
            assert len(orb["moCoefficients"]) % len(orb["energies"]) == 0

    if "basisSet" in cjson:
        basis = cjson["basisSet"]
        assert (
            len(basis["primitivesPerShell"])
            == len(basis["shellToAtomMap"])
            == len(basis["shellTypes"])
        )
        assert (
            sum(basis["primitivesPerShell"])
            == len(basis["exponents"])
            == len(basis["coefficients"])
        )
        assert 0 <= min(basis["shellToAtomMap"])
        assert max(basis["shellToAtomMap"]) < n_atoms

    if "spectra" in cjson and "electronic" in cjson["spectra"]:
        electronic = cjson["spectra"]["electronic"]
        assert len(electronic["energies"]) == len(electronic["intensities"])
        if "rotation" in electronic:
            assert len(electronic["rotation"]) == len(electronic["energies"])

    # The whole dict must survive a strict JSON round trip that rejects
    # NaN/Infinity, and must contain no non-JSON-serializable types.
    serialized = json.dumps(cjson, allow_nan=False)
    round_tripped = strict_json_loads(serialized)
    assert round_tripped == cjson
