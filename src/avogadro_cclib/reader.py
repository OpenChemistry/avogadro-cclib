"""
/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the 3-Clause BSD License, (see "LICENSE").
******************************************************************************/
"""

import math
import sys
import warnings
from io import StringIO

warnings.simplefilter("ignore", UserWarning)

from cclib.io.ccio import ccopen

from .utils import (
    _cclib_to_cjson_basis,
    _cclib_to_cjson_vibdisps,
    _cclib_to_cjson_mocoeffs,
)

# these are in cclib conversion too
HARTREE_TO_J_MOL = 2625499.638933033
EV_TO_J_MOL = 96485.33290025658
EV_TO_KJ_MOL = EV_TO_J_MOL / 1000.0


def _has_data(data, attr):
    """True if cclib set ``attr`` on ``data`` and it is not empty."""
    return hasattr(data, attr) and len(getattr(data, attr)) > 0


def read(file_content):
    log = ccopen(StringIO(file_content))
    if log is None:
        # ccopen returns None when it cannot identify the program that
        # produced the file. Report that plainly rather than letting an
        # AttributeError on None escape to the user.
        raise ValueError("cclib could not identify the program that produced this file")
    data = log.parse()

    has_coords = _has_data(data, "atomcoords")
    has_elements = _has_data(data, "atomnos")
    if not has_coords or not has_elements:
        raise ValueError("no atomic coordinates found in the output file")

    cjson = {"chemicalJson": 1, "atoms": {}}
    cjson["atoms"]["coords"] = {}
    cjson["atoms"]["coords"]["3d"] = data.atomcoords[-1].flatten().tolist()

    cjson["atoms"]["elements"] = {}
    cjson["atoms"]["elements"]["number"] = data.atomnos.tolist()

    if hasattr(data, "charge"):
        cjson.setdefault("properties", {})["totalCharge"] = data.charge

    if hasattr(data, "mult"):
        cjson.setdefault("properties", {})["totalSpinMultiplicity"] = data.mult

    # check for geometry optimization coords or scancoords
    has_scancoords = _has_data(data, "scancoords")
    has_scanenergies = _has_data(data, "scanenergies")
    if has_scancoords and has_scanenergies:
        steps = min(len(data.scancoords), len(data.scanenergies))
        energies = []
        cjson["atoms"]["coords"]["3dSets"] = []
        for i in range(steps):
            coords = data.scancoords[i].flatten().tolist()
            cjson["atoms"]["coords"]["3dSets"].append(coords)
            energies.append(data.scanenergies[i] * EV_TO_KJ_MOL)
        cjson.setdefault("properties", {})["energies"] = energies

    # Add calculated properties
    if _has_data(data, "scfenergies"):
        energy = data.scfenergies[-1] * EV_TO_KJ_MOL
        cjson.setdefault("properties", {})["totalEnergy"] = energy
        # A geometry per SCF energy is not guaranteed (e.g. single-point or
        # multi-step jobs may report many scfenergies for one geometry, or
        # vice versa), so drive this off the number of geometries we
        # actually have and only pair up as many steps as both sequences
        # provide. properties.energies stays exactly parallel to
        # atoms.coords.3dSets.
        if len(data.atomcoords) > 1:
            steps = min(len(data.atomcoords), len(data.scfenergies))
            # first frame defaults to optimized (i.e. the final geometry)
            energies = [data.scfenergies[-1] * EV_TO_KJ_MOL]
            coords = data.atomcoords[-1].flatten().tolist()
            cjson["atoms"]["coords"]["3dSets"] = [coords]
            for i in range(steps - 1):
                coords = data.atomcoords[i].flatten().tolist()
                cjson["atoms"]["coords"]["3dSets"].append(coords)
                energies.append(data.scfenergies[i] * EV_TO_KJ_MOL)
            cjson.setdefault("properties", {})["energies"] = energies

    # atomic partial charges
    if hasattr(data, "atomcharges"):
        for set in data.atomcharges.items():
            cjson.setdefault("partialCharges", {})[set[0]] = set[1].tolist()

    if _has_data(data, "gbasis"):
        basis = _cclib_to_cjson_basis(data.gbasis)
        cjson["basisSet"] = basis

    # Convert mo coefficients
    if _has_data(data, "mocoeffs"):
        mocoeffs = _cclib_to_cjson_mocoeffs(data.mocoeffs)
        cjson.setdefault("orbitals", {})["moCoefficients"] = mocoeffs

    # Convert mo energies
    if _has_data(data, "moenergies"):
        moenergies = list(data.moenergies[-1])
        cjson.setdefault("orbitals", {})["energies"] = moenergies

    if hasattr(data, "nelectrons"):
        cjson.setdefault("orbitals", {})["electronCount"] = int(data.nelectrons)

    if hasattr(data, "homos") and hasattr(data, "nmo"):
        homos = data.homos
        nmo = data.nmo
        if len(homos) == 1:
            occupations = [2 if i <= homos[0] else 0 for i in range(nmo)]
            cjson.setdefault("orbitals", {})["occupations"] = occupations

    if _has_data(data, "mosyms"):
        alpha_syms = data.mosyms[0]
        beta_syms = data.mosyms[1] if len(data.mosyms) > 1 else alpha_syms
        cjson.setdefault("orbitals", {})["symmetries"] = [alpha_syms, beta_syms]

    if _has_data(data, "vibfreqs"):
        vibfreqs = list(data.vibfreqs)
        cjson.setdefault("vibrations", {})["frequencies"] = vibfreqs

    if _has_data(data, "vibdisps"):
        vibdisps = _cclib_to_cjson_vibdisps(data.vibdisps)
        cjson.setdefault("vibrations", {})["eigenVectors"] = vibdisps

    # electronic spectra
    if _has_data(data, "etenergies"):
        # reported as wavenumbers, convert to eV
        etenergies = list(data.etenergies / 8065.544)
        etoscs = (
            list(data.etoscs) if hasattr(data, "etoscs") else [1.0] * len(etenergies)
        )
        cjson.setdefault("spectra", {})["electronic"] = {
            "energies": etenergies,
            "intensities": etoscs,
        }
        if hasattr(data, "etrotats"):
            cjson["spectra"]["electronic"]["rotation"] = list(data.etrotats)

    # nmr spectra
    if hasattr(data, "nmrtensors"):
        nmrshifts = []
        for atom in data.nmrtensors:
            total = data.nmrtensors[atom]["total"]
            isotropic = (total[0][0] + total[1][1] + total[2][2]) / 3.0
            nmrshifts.append(isotropic)
        cjson.setdefault("spectra", {})["nmr"] = {"shifts": nmrshifts}

    # vibrational intensities and metadata
    if "vibrations" in cjson and "frequencies" in cjson["vibrations"]:
        if hasattr(data, "vibirs"):
            cjson["vibrations"]["intensities"] = list(data.vibirs)
        else:
            cjson["vibrations"]["intensities"] = [
                1 for i in range(len(cjson["vibrations"]["frequencies"]))
            ]

        if hasattr(data, "vibramans"):
            cjson["vibrations"]["ramanIntensities"] = list(data.vibramans)

        if hasattr(data, "vibsyms"):
            cjson["vibrations"]["symmetries"] = list(data.vibsyms)

        if "modes" not in cjson["vibrations"]:
            cjson["vibrations"]["modes"] = [
                i + 1 for i in range(len(cjson["vibrations"]["frequencies"]))
            ]

    # Convert calculation metadata
    if hasattr(data, "metadata"):
        metadata = data.metadata
        if "basis_set" in metadata:
            cjson.setdefault("inputParameters", {})["basis"] = metadata[
                "basis_set"
            ].lower()
        if "functional" in metadata:
            cjson.setdefault("inputParameters", {})["functional"] = metadata[
                "functional"
            ].lower()
        if "methods" in metadata and len(metadata["methods"]) > 0:
            cjson.setdefault("inputParameters", {})["theory"] = metadata["methods"][
                -1
            ].lower()

    return _sanitize(cjson)


def _sanitize_value(value):
    """Convert numpy scalars to native types, reporting any non-finite float.

    Returns ``(converted_value, saw_non_finite)`` in a single pass over
    ``value`` (recursing through nested lists), instead of walking the data
    once to detect non-finite floats and a second time to unwrap numpy
    scalars.
    """
    if isinstance(value, list):
        converted = []
        dirty = False
        for item in value:
            item_value, item_dirty = _sanitize_value(item)
            dirty = dirty or item_dirty
            converted.append(item_value)
        return converted, dirty
    if hasattr(value, "item"):
        # numpy scalar (e.g. numpy.float64, numpy.int64, numpy.bool_) --
        # unwrap to a native Python type. Note numpy.float64 is a subclass
        # of the built-in float, so this check must come before any
        # isinstance(value, float) test.
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return value, True
    return value, False


def _sanitize(obj):
    """Make ``obj`` safe to pass to ``json.dumps``.

    ``json.dumps`` happily emits bare ``NaN``/``Infinity`` tokens, which are
    not valid JSON and are rejected by Avogadro's C++ CJSON reader. Any
    key whose value (scalar or list, however deeply nested) contains a
    non-finite float is dropped entirely (with a warning to stderr) rather
    than having its data replaced -- we never fabricate scientific data.
    Numpy scalar types (e.g. numpy.float64, numpy.int64) are coerced to
    plain Python int/float so the result is always JSON-serializable.
    """
    sanitized = {}
    for key, value in obj.items():
        if isinstance(value, dict):
            # Nested dicts (e.g. spectra.electronic) are structure, not
            # leaf values -- recurse instead of testing them as a whole.
            sanitized[key] = _sanitize(value)
            continue
        converted, dirty = _sanitize_value(value)
        if dirty:
            print(
                f"avogadro-cclib: dropping non-finite values found in '{key}'",
                file=sys.stderr,
            )
            continue
        sanitized[key] = converted
    return sanitized
