"""
clashscore.py

Open-source replacement for `phenix.clashscore`, built directly on the
mmtbx.validation.clashscore API (the same code phenix.clashscore calls).

Computes the MolProbity all-atom clashscore for a PDB/mmCIF model:
the number of serious steric clashes (>0.4 Angstrom overlap between
non-bonded atoms) per 1000 atoms.

Usage:
    python clashscore.py model.pdb
    python clashscore.py model.pdb --out clashscore.txt

As a library:
    from clashscore import run_clashscore
    result = run_clashscore("model.pdb")
    print(result.get_clashscore())
"""

import argparse
import sys

import iotbx.pdb
from mmtbx.validation import clashscore


def run_clashscore(model_path, nuclear=False, keep_hydrogens=False,
                    fast=False):
    """
    Run clashscore analysis on a PDB/mmCIF file.

    Args:
        model_path: path to a .pdb or .cif file
        nuclear: use nuclear (neutron) X-H distances instead of electron-cloud
        keep_hydrogens: if False, existing hydrogens are stripped and
            clashscore adds its own via its internal reduce call
        fast: skip the slower, more exhaustive symmetry-clash checks

    Returns:
        an mmtbx.validation.clashscore.clashscore results object
    """
    pdb_input = iotbx.pdb.input(file_name=model_path)
    pdb_hierarchy = pdb_input.construct_hierarchy()

    result = clashscore.clashscore(
        pdb_hierarchy=pdb_hierarchy,
        keep_hydrogens=keep_hydrogens,
        nuclear=nuclear,
        fast=fast,
        condensed_probe=True,
    )
    return result


def format_report(result, model_path):
    lines = []
    lines.append("clashscore analysis for {}".format(model_path))
    lines.append("-" * 60)
    lines.append("Clashscore (clashes per 1000 atoms): {:.2f}".format(
        result.get_clashscore()))
    lines.append("")
    lines.append("Individual clashes:")
    # show_old_output mirrors phenix.clashscore's per-clash table
    lines.append(result.as_text() if hasattr(result, "as_text") else "")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", help="Path to PDB/mmCIF model file")
    parser.add_argument("--out", help="Write report to this file "
                                       "(default: stdout)")
    parser.add_argument("--nuclear", action="store_true",
                         help="Use nuclear X-H distances")
    parser.add_argument("--keep-hydrogens", action="store_true",
                         help="Do not strip/re-add hydrogens before scoring")
    parser.add_argument("--fast", action="store_true",
                         help="Skip exhaustive symmetry clash checks")
    args = parser.parse_args()

    result = run_clashscore(
        args.model,
        nuclear=args.nuclear,
        keep_hydrogens=args.keep_hydrogens,
        fast=args.fast,
    )

    report = format_report(result, args.model)

    if args.out:
        with open(args.out, "w") as fh:
            fh.write(report + "\n")
    else:
        sys.stdout.write(report + "\n")


if __name__ == "__main__":
    main()
