"""
fmodel_builder.py

Builds a scaled mmtbx.f_model.manager ("fmodel") from a model +
reflection-data pair, handling bulk-solvent correction, anisotropic
scaling, and R-free flag interpretation via CCTBX's own helper
(mmtbx.command_line.load_model_and_data) rather than reimplementing
any of that math.

This is the piece `molprobity.py` needs in order to report real-space
correlation and other model-vs-data metrics rather than model-only
geometry statistics.

Usage as a library:
    from fmodel_builder import build_fmodel
    fmodel = build_fmodel("model.pdb", "data.mtz")
    print(fmodel.r_work(), fmodel.r_free())

Usage as a script (quick R-factor sanity check):
    python fmodel_builder.py model.pdb data.mtz
    python fmodel_builder.py model.pdb data.mtz --twin-law "-h,-k,l"
"""

import argparse
import sys


def build_fmodel(model_path, mtz_path, twin_law=None,
                  prefer_anomalous=False):
    """
    Build a scaled fmodel from a model file and reflection file.

    Args:
        model_path: path to a PDB/mmCIF model
        mtz_path: path to a reflection file with observed
            amplitudes/intensities and (ideally) an R-free flag column
        twin_law: optional twin law string (e.g. "-h,-k,l") if Xtriage
            flagged the data as twinned; passing this builds a
            twin-aware fmodel so R-factors match what
            phenix.molprobity/phenix.refine would report
        prefer_anomalous: whether to prefer anomalous data arrays if
            both anomalous and merged arrays are present in the file

    Returns:
        an mmtbx.f_model.manager object with scaling already applied
        (update_all_scales() has been called)

    Raises:
        RuntimeError if the model and data could not be reconciled
        (mismatched symmetry, no usable amplitude/intensity array,
        missing/unusable R-free flags, etc.) -- the underlying
        exception message from load_model_and_data is preserved since
        it's usually specific about which check failed.
    """
    from mmtbx.command_line import load_model_and_data
    import mmtbx.utils
    from iotbx import phil as iotbx_phil

    # Build a minimal PHIL scope covering the parameters
    # load_model_and_data expects; twin_law is passed through if given.
    master_phil = iotbx_phil.parse(mmtbx.utils.cmdline_input_phil_str)

    args = [model_path, mtz_path]
    if twin_law:
        args.append("twin_law={}".format(twin_law))

    try:
        cmdline = load_model_and_data(
            args=args,
            master_phil=master_phil,
            process_pdb_file=True,
            create_fmodel=True,
            prefer_anomalous=prefer_anomalous,
            out=sys.stderr,  # loader's own diagnostic messages
        )
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(
            "Could not build fmodel from {} + {}: {}\n"
            "Common causes: mismatched space group/cell between model "
            "and data, no R-free flag column, or twinned data that "
            "needs --twin-law.".format(model_path, mtz_path, exc)
        ) from exc

    fmodel = cmdline.fmodel

    # Apply bulk-solvent correction + anisotropic scaling explicitly,
    # in case the loader's defaults didn't already trigger it.
    fmodel.update_all_scales()

    return fmodel


def summarize(fmodel):
    lines = []
    lines.append("R-work: {:.4f}".format(fmodel.r_work()))
    lines.append("R-free: {:.4f}".format(fmodel.r_free()))
    try:
        lines.append("Overall CC (work): {:.4f}".format(
            fmodel.r_work_scale_k1()))
    except Exception:  # noqa: BLE001
        pass
    lines.append("Resolution range: {:.2f} - {:.2f} A".format(
        *fmodel.f_obs().d_max_min()))
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", help="Path to PDB/mmCIF model file")
    parser.add_argument("mtz", help="Path to reflection file")
    parser.add_argument("--twin-law", help="Twin law, e.g. '-h,-k,l', "
                                            "if Xtriage flagged twinning")
    parser.add_argument("--prefer-anomalous", action="store_true")
    args = parser.parse_args()

    fmodel = build_fmodel(
        args.model, args.mtz,
        twin_law=args.twin_law,
        prefer_anomalous=args.prefer_anomalous,
    )
    print(summarize(fmodel))


if __name__ == "__main__":
    main()
