"""
/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the 3-Clause BSD License, (see "LICENSE").
******************************************************************************/
"""

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


def read(file_content):
    log = ccopen(StringIO(file_content))
    data = log.parse()

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
    if hasattr(data, "scancoords"):
        steps = len(data.scanenergies)
        energies = []
        cjson["atoms"]["coords"]["3dSets"] = []
        for i in range(steps):
            coords = data.scancoords[i].flatten().tolist()
            cjson["atoms"]["coords"]["3dSets"].append(coords)
            energies.append(data.scanenergies[i] * EV_TO_KJ_MOL)
        cjson.setdefault("properties", {})["energies"] = energies

    # Add calculated properties
    if hasattr(data, "scfenergies"):
        if len(data.scfenergies) > 0:
            energy = data.scfenergies[-1] * EV_TO_KJ_MOL
            cjson.setdefault("properties", {})["totalEnergy"] = energy
        if len(data.scfenergies) > 1:  # optimization!
            steps = len(data.scfenergies)
            # first frame defaults to optimized
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

    if hasattr(data, "gbasis"):
        basis = _cclib_to_cjson_basis(data.gbasis)
        cjson["basisSet"] = basis

    # Convert mo coefficients
    if hasattr(data, "mocoeffs"):
        mocoeffs = _cclib_to_cjson_mocoeffs(data.mocoeffs)
        cjson.setdefault("orbitals", {})["moCoefficients"] = mocoeffs

    # Convert mo energies
    if hasattr(data, "moenergies"):
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

    if hasattr(data, "mosyms"):
        alpha_syms = data.mosyms[0]
        beta_syms = data.mosyms[1] if len(data.mosyms) > 1 else alpha_syms
        cjson.setdefault("orbitals", {})["symmetries"] = [alpha_syms, beta_syms]

    if hasattr(data, "vibfreqs"):
        vibfreqs = list(data.vibfreqs)
        cjson.setdefault("vibrations", {})["frequencies"] = vibfreqs

    if hasattr(data, "vibdisps"):
        vibdisps = _cclib_to_cjson_vibdisps(data.vibdisps)
        cjson.setdefault("vibrations", {})["eigenVectors"] = vibdisps

    # electronic spectra
    if hasattr(data, "etenergies"):
        # reported as wavenumbers, convert to eV
        etenergies = list(data.etenergies / 8065.544)
        etoscs = list(data.etoscs) if hasattr(data, "etoscs") else [1.0] * len(etenergies)
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
            cjson.setdefault("inputParameters", {})["basis"] = metadata["basis_set"].lower()
        if "functional" in metadata:
            cjson.setdefault("inputParameters", {})["functional"] = metadata["functional"].lower()
        if "methods" in metadata and len(metadata["methods"]) > 0:
            cjson.setdefault("inputParameters", {})["theory"] = metadata["methods"][-1].lower()

    return cjson
