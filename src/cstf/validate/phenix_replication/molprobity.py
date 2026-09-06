"""
molprobity.py

Open-source replacement for `phenix.molprobity`, the aggregate
multi-criterion validation report. Built on
mmtbx.validation.molprobity.molprobity, which is the same class the
phenix.molprobity CLI instantiates: it runs clashscore, Ramachandran
(ramalyze), rotamer (rotalyze), CaBLAM, and (when reflection data is
supplied) real-space/model-vs-data metrics, then bundles them into one
summary object with an overall MolProbity score.

When an --mtz is supplied, this module builds a scaled fmodel via
fmodel_builder.build_fmodel() (bulk-solvent correction + anisotropic
scaling handled by CCTBX's own loader, not reimplemented here) so the
real-space correlation and other model-vs-data metrics actually
populate instead of being silently skipped.

Usage:
    python molprobity.py model.pdb
    python molprobity.py model.pdb --mtz data.mtz --out molprobity.out
    python molprobity.py model.pdb --mtz data.mtz --twin-law "-h,-k,l"
"""

import argparse
import sys

import iotbx.pdb
from mmtbx.validation import molprobity as molprobity_module

from cstf.validate.phenix_replication.fmodel_builder import build_fmodel


def run_molprobity(model_path, mtz_path=None, keep_hydrogens=False, twin_law=None):
    """
    Run the full MolProbity multi-criterion validation suite.

    Args:
        model_path: path to a PDB/mmCIF file
        mtz_path: optional reflection file for model-vs-data metrics
            (real-space correlation, etc.) -- when given, a scaled
            fmodel is built via fmodel_builder.build_fmodel()
        keep_hydrogens: whether to keep existing hydrogens rather than
            stripping/re-adding them internally
        twin_law: optional twin law (e.g. "-h,-k,l") passed through to
            build_fmodel() if Xtriage flagged the data as twinned;
            ignored if mtz_path is not given

    Returns:
        an mmtbx.validation.molprobity.molprobity results object

    Raises:
        RuntimeError if mtz_path is given but the fmodel could not be
        built (mismatched symmetry, missing R-free flags, etc.) --
        see fmodel_builder.build_fmodel() for the underlying cause.
        Model-only validation is not silently substituted in this
        case, since a caller who asked for data-based metrics should
        know they didn't get them rather than get a partial report
        that looks complete.
    """
    pdb_input = iotbx.pdb.input(file_name=model_path)
    pdb_hierarchy = pdb_input.construct_hierarchy()

    fmodel = None
    if mtz_path:
        fmodel = build_fmodel(model_path, mtz_path, twin_law=twin_law)

    result = molprobity_module.molprobity(
        pdb_hierarchy=pdb_hierarchy,
        fmodel=fmodel,
        keep_hydrogens=keep_hydrogens,
        nuclear=False,
        save_probe_unformatted_file=None,
    )
    return result


def format_report(result, model_path):
    lines = []
    lines.append("MolProbity summary for {}".format(model_path))
    lines.append("-" * 60)

    # The molprobity object exposes a summarize() / show() style report;
    # also surface the headline numbers explicitly for scripting.
    if hasattr(result, "clashscore"):
        try:
            lines.append("Clashscore: {:.2f}".format(result.clashscore()))
        except Exception:  # noqa: BLE001
            pass
    if hasattr(result, "rama_favored"):
        try:
            lines.append("Ramachandran favored: {:.1f}%".format(result.rama_favored()))
            lines.append(
                "Ramachandran outliers: {:.1f}%".format(result.rama_outliers())
            )
        except Exception:  # noqa: BLE001
            pass
    if hasattr(result, "rota_outliers"):
        try:
            lines.append("Rotamer outliers: {:.1f}%".format(result.rota_outliers()))
        except Exception:  # noqa: BLE001
            pass
    if hasattr(result, "overall_score"):
        try:
            lines.append("MolProbity score: {:.2f}".format(result.overall_score()))
        except Exception:  # noqa: BLE001
            pass

    # Model-vs-data metrics: only present when an fmodel was passed in.
    if hasattr(result, "r_work") and hasattr(result, "r_free"):
        try:
            lines.append("R-work: {:.4f}".format(result.r_work()))
            lines.append("R-free: {:.4f}".format(result.r_free()))
        except Exception:  # noqa: BLE001
            pass
    if hasattr(result, "real_space_correlation"):
        try:
            lines.append(
                "Real-space correlation: {:.3f}".format(result.real_space_correlation())
            )
        except Exception:  # noqa: BLE001
            pass

    lines.append("")
    lines.append("Full detail:")
    if hasattr(result, "as_text"):
        lines.append(result.as_text())
    elif hasattr(result, "show"):
        import io

        buf = io.StringIO()
        result.show(out=buf)
        lines.append(buf.getvalue())

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", help="Path to PDB/mmCIF model file")
    parser.add_argument(
        "--mtz", help="Optional reflection file for " "model-vs-data metrics"
    )
    parser.add_argument(
        "--twin-law",
        help="Twin law, e.g. '-h,-k,l', "
        "if Xtriage flagged twinning "
        "(only used with --mtz)",
    )
    parser.add_argument("--out", help="Write report to this file " "(default: stdout)")
    parser.add_argument(
        "--keep-hydrogens",
        action="store_true",
        help="Do not strip/re-add hydrogens internally",
    )
    args = parser.parse_args()

    try:
        result = run_molprobity(
            args.model,
            mtz_path=args.mtz,
            keep_hydrogens=args.keep_hydrogens,
            twin_law=args.twin_law,
        )
    except RuntimeError as exc:
        sys.stderr.write("Error building fmodel: {}\n".format(exc))
        sys.exit(1)

    report = format_report(result, args.model)

    if args.out:
        with open(args.out, "w") as fh:
            fh.write(report + "\n")
    else:
        sys.stdout.write(report + "\n")


if __name__ == "__main__":
    main()
