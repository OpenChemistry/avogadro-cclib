"""
/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the 3-Clause BSD License, (see "LICENSE").
******************************************************************************/
"""
import argparse
import json
import sys
from io import StringIO 

import warnings
warnings.simplefilter("ignore", UserWarning)

from cclib.io.ccio import ccopen

# Code from openchemistry
# https://github.com/OpenChemistry/openchemistrypy/blob/master/openchemistry/io/
from utils import (
    _cclib_to_cjson_basis,
    _cclib_to_cjson_vibdisps,
    _cclib_to_cjson_mocoeffs,
)

# these are in cclib conversion too
HARTREE_TO_J_MOL = 2625499.638933033
EV_TO_J_MOL = 96485.33290025658
EV_TO_KJ_MOL = EV_TO_J_MOL / 1000.0


def getMetaData():
    metaData = {}
    metaData["inputFormat"] = "cjson"
    metaData["outputFormat"] = "cjson"
    metaData["operations"] = ["read"]
    metaData["identifier"] = "cclib"
    metaData["name"] = "cclib"
    metaData["description"] = (
        "The cclib script provided by the cclib repository is used to "
        + "write the CJSON format using the input file provided "
        + "to Avogadro2."
    )
    metaData["fileExtensions"] = [
        "out",
        "log",
        "adfout",
        "g09",
        "g03",
        "g98",
        "fchk",
    ]
    metaData["mimeTypes"] = [""]
    metaData["bond"] = True
    return metaData


def read():
    # Pass the standard input to ccopen:
    input = sys.stdin.read()
    log = ccopen(StringIO(input))
    data = log.parse()

    cjson = {"chemicalJson": 1, "atoms": {}}  # version number
    cjson["atoms"]["coords"] = {}
    cjson["atoms"]["coords"]["3d"] = data.atomcoords[-1].flatten().tolist()
    # 3dSets

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
        # add in the energies
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
        # iterate through the methods and charges
        for set in data.atomcharges.items():
            #if set.endswith("_sum"):
            #    continue  # adds hydrogens into the other atoms .. skip

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
        if len(data.mosyms) > 1:
            beta_syms = data.mosyms[1]
        else:
            beta_syms = alpha_syms
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
        if hasattr(data, "etoscs"):
            etoscs = list(data.etoscs)
        else:
            etoscs = [1.0] * len(etenergies)
        cjson.setdefault("spectra", {})["electronic"] = {
            "energies": etenergies,
            "intensities": etoscs,
        }
        if hasattr(data, "etrotats"):
            etrotats = list(data.etrotats)
            cjson["spectra"]["electronic"]["rotation"] = etrotats

    # nmr spectra
    if hasattr(data, "nmrtensors"):
        nmrshifts = []
        # we want the isotropic shift, so take the trace
        # from the "total" tensor
        for atom in data.nmrtensors:
            total = data.nmrtensors[atom]["total"]
            isotropic = (total[0][0] + total[1][1] + total[2][2]) / 3.0
            nmrshifts.append(isotropic)
        cjson.setdefault("spectra", {})["nmr"] = {
            "shifts": nmrshifts,
        }

    # Add an intensities array
    if "vibrations" in cjson and "frequencies" in cjson["vibrations"]:
        if hasattr(data, "vibirs"):
            cjson["vibrations"]["intensities"] = list(data.vibirs)
        else:
            # placeholder = everything is 1
            cjson["vibrations"]["intensities"] = [
                1 for i in range(len(cjson["vibrations"]["frequencies"]))
            ]

        # check for raman intensities
        if hasattr(data, "vibramans"):
            cjson["vibrations"]["ramanIntensities"] = list(data.vibramans)

        # check for symmetries
        if hasattr(data, "vibsyms"):
            cjson["vibrations"]["symmetries"] = list(data.vibsyms)

        if "modes" not in cjson["vibrations"]:
            #  count of modes ..  seems  redundant, but required
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

    return json.dumps(cjson)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Read files using cclib")
    parser.add_argument("--metadata", action="store_true")
    parser.add_argument("--read", action="store_true")
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--display-name", action="store_true")
    parser.add_argument("--lang", nargs="?", default="en")
    args = vars(parser.parse_args())

    if args["metadata"]:
        print(json.dumps(getMetaData()))
    elif args["display_name"]:
        print(getMetaData()["name"])
    elif args["read"]:
        print(read())
    elif args["write"]:
        pass
