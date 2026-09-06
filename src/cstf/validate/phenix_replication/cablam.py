"""
cablam.py

Open-source replacement for `phenix.cablam`, built on
mmtbx.validation.cablam - the CaBLAM backbone-conformation validation
(CA-based Local Analysis of Backbone), useful particularly for flagging
suspect loop geometry that Ramachandran analysis alone can miss.

Usage:
    python cablam.py model.pdb
    python cablam.py model.pdb --out cablam.out
    python cablam.py model.pdb --kinemage cablam.kin
"""

import argparse
import sys

import iotbx.pdb
from mmtbx.validation import cablam as cablam_module


def run_cablam(model_path):
    """
    Run CaBLAM analysis on a PDB/mmCIF file.

    Returns:
        an mmtbx.validation.cablam.cablamalyze results object
    """
    pdb_input = iotbx.pdb.input(file_name=model_path)
    pdb_hierarchy = pdb_input.construct_hierarchy()
    pdb_hierarchy.atoms().reset_i_seq()
    print(type(cablam_module.cablamalyze))
    result = cablam_module.cablamalyze(
        pdb_hierarchy=pdb_hierarchy,
        outliers_only=False,
        out=None,
        quiet=True,
    )
    return result


def format_report(result, model_path):
    lines = []
    lines.append("CaBLAM analysis for {}".format(model_path))
    lines.append("-" * 60)
    # as_text()/show_summary() names vary slightly by version; try the
    # common ones and fall back to iterating results directly.
    if hasattr(result, "as_text"):
        lines.append(result.as_text())
    elif hasattr(result, "show_old_output"):
        import io

        buf = io.StringIO()
        result.show_old_output(out=buf)
        lines.append(buf.getvalue())
    else:
        for r in result.results:
            lines.append(str(r))
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", help="Path to PDB/mmCIF model file")
    parser.add_argument("--out", help="Write report to this file " "(default: stdout)")
    parser.add_argument(
        "--kinemage", help="Optional path to write a " "kinemage visualization file"
    )
    args = parser.parse_args()

    result = run_cablam(args.model)
    report = format_report(result, args.model)

    if args.out:
        with open(args.out, "w") as fh:
            fh.write(report + "\n")
    else:
        sys.stdout.write(report + "\n")

    if args.kinemage and hasattr(result, "as_kinemage"):
        with open(args.kinemage, "w") as fh:
            fh.write(result.as_kinemage())


if __name__ == "__main__":
    main()
