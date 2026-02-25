"""Avogadro plugin for reading computational chemistry output files via cclib."""

import argparse
import json
import sys


def main():
    # Avogadro calls the plugin as:
    #   avogadro-cclib <identifier> [--read] [--write] [--lang <locale>] [--debug]
    # with the file content (or molecule JSON) on stdin.
    parser = argparse.ArgumentParser()
    parser.add_argument("feature")
    parser.add_argument("--lang", nargs="?", default="en")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--read", action="store_true")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()

    avo_input = json.load(sys.stdin)
    output = None

    match args.feature:
        case "cclib":
            from .reader import read
            if args.read:
                output = {"cjson": read(avo_input["file"])}

    if output is not None:
        print(json.dumps(output))
