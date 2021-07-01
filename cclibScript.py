"""
/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the 3-Clause BSD License, (see "LICENSE").
******************************************************************************/
"""
import argparse
import json
import sys

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

def getMetaData():
    metaData = {}
    metaData['inputFormat'] = 'cjson'
    metaData['outputFormat'] = 'cjson'
    metaData['operations'] = ['read']
    metaData['identifier'] = 'cclib'
    metaData['name'] = 'cclib'
    metaData['description'] = "The cclib script provided by the cclib repository is used to " +\
                              "write the CJSON format using the input file provided " +\
                              "to Avogadro2."
    metaData['fileExtensions'] = ['out', 'log', 'adfout', 'g09', '.g03', '.g98', '.fchk']
    metaData['mimeTypes'] = ['']
    return metaData


def read():
    # Pass the standard input to ccopen:
    log = ccopen(sys.stdin)
    data = log.parse()

    cjson = { 
        'chemicalJson': 1, # version number
        'atoms': {}
    }
    cjson['atoms']['coords'] = {}
    cjson['atoms']['coords']['3d'] = data.atomcoords[-1].flatten().tolist()
    # 3dSets

    cjson['atoms']['elements'] = {}
    cjson['atoms']['elements']['number'] = data.atomnos.tolist()

    # Add calculated properties
    if hasattr(data, "scfenergies"):
        if len(data.scfenergies) > 0:
            energy = data.scfenergies[-1] * EV_TO_J_MOL
            cjson.setdefault("properties", {})["totalEnergy"] = energy

    if hasattr(data, "gbasis"):
        basis = _cclib_to_cjson_basis(data.gbasis)
        cjson["basisSet"] = basis

    # Convert mo coefficients
    if hasattr(data, 'mocoeffs'):
        mocoeffs = _cclib_to_cjson_mocoeffs(data.mocoeffs)
        cjson.setdefault('orbitals', {})['moCoefficients'] = mocoeffs

    # Convert mo energies
    if hasattr(data, 'moenergies'):
        moenergies = list(data.moenergies[-1])
        cjson.setdefault('orbitals', {})['energies'] = moenergies

    if hasattr(data, 'nelectrons'):
        cjson.setdefault('orbitals', {})['electronCount'] = int(data.nelectrons)

    if hasattr(data, 'homos') and hasattr(data, 'nmo'):
        homos = data.homos
        nmo = data.nmo
        if len(homos) == 1:
            occupations = [2 if i <= homos[0] else 0 for i in range(nmo)]
            cjson.setdefault('orbitals', {})['occupations'] = occupations

    if hasattr(data, "vibfreqs"):
        vibfreqs = list(data.vibfreqs)
        cjson.setdefault("vibrations", {})["frequencies"] = vibfreqs

    if hasattr(data, "vibdisps"):
        vibdisps = _cclib_to_cjson_vibdisps(data.vibdisps)
        cjson.setdefault("vibrations", {})["eigenVectors"] = vibdisps

    # Add a placeholder intensities array
    if "vibrations" in cjson and "frequencies" in cjson["vibrations"]:
        if hasattr(data, "vibirs"):
            cjson["vibrations"]["intensities"] = list(data.vibirs)
        else:
            # placeholder = everything is 1
            cjson["vibrations"]["intensities"] = [
                1 for i in range(len(cjson["vibrations"]["frequencies"]))
            ]
        if "modes" not in cjson["vibrations"]:
            #  count of modes ..  seems  redundant, but required
            cjson["vibrations"]["modes"] = [
                i + 1 for i in range(len(cjson["vibrations"]["frequencies"]))
            ]

    # Convert calculation metadata
    if hasattr(data, 'metadata'):
        metadata = data.metadata
        if 'basis_set' in metadata:
            cjson.setdefault('metadata', {})['basisSet'] = metadata['basis_set'].lower()
        if 'functional' in metadata:
            cjson.setdefault('metadata', {})['functional'] = metadata['functional'].lower()
        if 'methods' in metadata and len(metadata['methods']) > 0:
            cjson.setdefault('metadata', {})['theory'] = metadata['methods'][0].lower()

    return json.dumps(cjson)


if __name__ == "__main__":
    parser = argparse.ArgumentParser('Read files using cclib')
    parser.add_argument('--metadata', action='store_true')
    parser.add_argument('--read', action='store_true')
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--display-name', action='store_true')
    parser.add_argument('--lang', nargs='?', default='en')
    args = vars(parser.parse_args())

    if args['metadata']:
        print(json.dumps(getMetaData()))
    elif args['display_name']:
        print(getMetaData()['name'])
    elif args['read']:
        print(read())
    elif args['write']:
        pass
