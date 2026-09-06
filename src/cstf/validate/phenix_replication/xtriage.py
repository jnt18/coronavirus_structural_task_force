"""
xtriage.py

Open-source replacement for `phenix.xtriage`, built directly on
mmtbx.scaling.xtriage.run - the exact function phenix.xtriage's CLI
calls. Runs data-quality analysis (Wilson statistics, twinning tests,
anisotropy, translational NCS, ice rings) on reflection data, and
optionally cross-checks against a structure model.

Usage:
    python xtriage.py data.mtz
    python xtriage.py data.mtz --model model.pdb --out Xtriage_output.log
"""

import argparse
import sys

from iotbx.reflection_file_reader import any_reflection_file
from iotbx import reflection_file_reader
from mmtbx.scaling import xtriage
import iotbx.pdb


def load_miller_arrays(mtz_path):
    hkl_in = any_reflection_file(mtz_path)
    miller_arrays = hkl_in.file_content().as_miller_arrays()
    if not miller_arrays:
        raise ValueError(
            "No miller arrays found in {} - check the file is a valid "
            "reflection file with intensity/amplitude data.".format(mtz_path)
        )
    return miller_arrays


def load_xray_structure(model_path):
    pdb_input = iotbx.pdb.input(file_name=model_path)
    return pdb_input.xray_structure_simple()


def run_xtriage(mtz_path, model_path=None, out_stream=None):
    """
    Run Xtriage analysis.

    Args:
        mtz_path: path to reflection file (.mtz, .sca, etc.)
        model_path: optional path to a PDB/mmCIF model for
            model-dependent checks (e.g. twin-aware R-factor analysis)
        out_stream: file-like object to receive the human-readable
            report as it's generated (e.g. an open file or sys.stdout).
            Pass None to suppress and only get the results object back.

    Returns:
        the xtriage results object (mmtbx.scaling.xtriage.xtriage_analyses)
    """
    miller_arrays = load_miller_arrays(mtz_path)

    kwargs = dict(miller_arrays=miller_arrays, text_out=out_stream)

    if model_path:
        kwargs["xray_structure"] = load_xray_structure(model_path)

    results = xtriage.run(**kwargs)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mtz", help="Path to reflection file (.mtz, etc.)")
    parser.add_argument(
        "--model", help="Optional PDB/mmCIF model for " "model-dependent checks"
    )
    parser.add_argument(
        "--out", help="Write the report to this file " "(default: stdout)"
    )
    args = parser.parse_args()

    if args.out:
        with open(args.out, "w") as fh:
            run_xtriage(args.mtz, model_path=args.model, out_stream=fh)
        sys.stderr.write("Xtriage report written to {}\n".format(args.out))
    else:
        run_xtriage(args.mtz, model_path=args.model, out_stream=sys.stdout)


if __name__ == "__main__":
    main()
