"""
orchestrator.py

Runs the open-source validation pipeline (xtriage, cablam, clashscore,
reduce, molprobity) over a batch of PDB entries described by a
DataFrame indexed by PDB ID with a `path_in_repo` column pointing at a
directory containing one .pdb, .cif, and .mtz file each -- the same
layout the original run_xtriage.sh / run_cablam.sh /
run_molprobity.sh scripts walked.

Output layout mirrors the original scripts:
    {path_in_repo}/validation/Xtriage_output.log
    {path_in_repo}/validation/molprobity/cablam.out
    {path_in_repo}/validation/molprobity/clashscore.txt
    {path_in_repo}/validation/molprobity/{pdb_id}.H.pdb
    {path_in_repo}/validation/molprobity/molprobity.out

KNOWN GAP: the original run_molprobity.sh also called
`rama_chart_pdf` and `multichart` to generate Ramachandran/
multi-criterion PDF plots. No module for those was written in this
set (they're Phenix-specific plotting utilities, not just data
computed by mmtbx.validation) -- this orchestrator computes all the
same underlying statistics but does not reproduce those plot files.
If you need them, the ramalyze/rotalyze/cablam result objects here
carry the residue-level data needed to plot equivalents with
matplotlib; that would be a separate module.

Usage:
    import pandas as pd
    from orchestrator import run_pipeline

    df = pd.DataFrame({"path_in_repo": [...]}, index=["1abc", "2xyz"])
    summary = run_pipeline(df)
    summary.to_csv("validation_run_summary.csv")
"""

import logging
import sys
import traceback
from pathlib import Path

import pandas as pd

from cstf.validate.phenix_replication import xtriage as xtriage_mod
from cstf.validate.phenix_replication import cablam as cablam_mod
from cstf.validate.phenix_replication import clashscore as clashscore_mod
from cstf.validate.phenix_replication import reduce as reduce_mod
from cstf.validate.phenix_replication import molprobity as molprobity_mod

logger = logging.getLogger("xtal_validation.orchestrator")

STEP_NAMES = ["xtriage", "cablam", "clashscore", "reduce", "molprobity"]


def find_structure_files(path_in_repo):
    """
    Locate the .pdb, .cif, and .mtz files in a directory.

    Returns a dict with keys 'pdb', 'cif', 'mtz', each either a Path
    or None if not found. Matching is by extension only (case
    insensitive), not by filename, since callers may use varying
    PDB-ID casing conventions (e.g. pdb_redo's lowercase entries).
    """
    path = Path(path_in_repo)
    found = {"pdb": None, "cif": None, "mtz": None}
    if not path.is_dir():
        return found

    for f in path.iterdir():
        if not f.is_file():
            continue
        suffix = f.suffix.lower()
        if suffix == ".pdb" and found["pdb"] is None:
            found["pdb"] = f
        elif suffix == ".cif" and found["cif"] is None:
            found["cif"] = f
        elif suffix == ".mtz" and found["mtz"] is None:
            found["mtz"] = f
    return found


def _model_path(files):
    """Prefer .pdb over .cif for the model input, matching the
    original scripts' behavior of operating on the .pdb file."""
    return files["pdb"] or files["cif"]


def process_one(pdb_id, path_in_repo, skip_existing=True):
    """
    Run the full validation pipeline for a single PDB entry.

    Args:
        pdb_id: identifier, used for output filenames (e.g. the
            {pdb_id}.H.pdb hydrogenated model)
        path_in_repo: directory containing the entry's .pdb/.cif/.mtz
        skip_existing: if True (default, matching the original bash
            scripts), skip a step whose output file already exists
            rather than recomputing it

    Returns:
        a dict with one status entry per step in STEP_NAMES, each
        either 'done', 'skipped', 'no_input' (required input file
        missing), or 'error: <message>'. Never raises -- a failure in
        one step is recorded and the remaining steps for that entry
        are still attempted where feasible, and failures in one entry
        never stop processing of other entries in run_pipeline().
    """
    path = Path(path_in_repo)
    status = {step: "not_run" for step in STEP_NAMES}

    files = find_structure_files(path)
    model_path = _model_path(files)

    validation_dir = path / "validation"
    molprobity_dir = validation_dir / "molprobity"

    if model_path is None:
        for step in STEP_NAMES:
            status[step] = "no_input: no .pdb or .cif file found"
        return status

    validation_dir.mkdir(exist_ok=True)
    molprobity_dir.mkdir(exist_ok=True)

    # --- xtriage: needs reflection data ---
    xtriage_log = validation_dir / "Xtriage_output.log"
    if files["mtz"] is None:
        status["xtriage"] = "no_input: no .mtz file found"
    elif skip_existing and xtriage_log.exists():
        status["xtriage"] = "skipped"
    else:
        try:
            with open(xtriage_log, "w") as fh:
                xtriage_mod.run_xtriage(
                    str(files["mtz"]),
                    model_path=str(model_path),
                    out_stream=fh,
                )
            status["xtriage"] = "done"
        except Exception as exc:  # noqa: BLE001
            status["xtriage"] = "error: {}".format(exc)
            logger.warning("xtriage failed for %s: %s", pdb_id, exc)
            logger.debug(traceback.format_exc())

    # --- cablam ---
    cablam_out = molprobity_dir / "cablam.out"
    if skip_existing and cablam_out.exists():
        status["cablam"] = "skipped"
    else:
        try:
            result = cablam_mod.run_cablam(str(model_path))
            report = cablam_mod.format_report(result, str(model_path))
            cablam_out.write_text(report + "\n")
            status["cablam"] = "done"
        except Exception as exc:  # noqa: BLE001
            status["cablam"] = "error: {}".format(exc)
            logger.warning("cablam failed for %s: %s", pdb_id, exc)
            logger.debug(traceback.format_exc())

    # --- clashscore ---
    clashscore_txt = molprobity_dir / "clashscore.txt"
    if skip_existing and clashscore_txt.exists():
        status["clashscore"] = "skipped"
    else:
        try:
            result = clashscore_mod.run_clashscore(str(model_path))
            report = clashscore_mod.format_report(result, str(model_path))
            clashscore_txt.write_text(report + "\n")
            status["clashscore"] = "done"
        except Exception as exc:  # noqa: BLE001
            status["clashscore"] = "error: {}".format(exc)
            logger.warning("clashscore failed for %s: %s", pdb_id, exc)
            logger.debug(traceback.format_exc())

    # --- reduce (add hydrogens) ---
    h_pdb = molprobity_dir / "{}.H.pdb".format(pdb_id)
    if skip_existing and h_pdb.exists():
        status["reduce"] = "skipped"
    else:
        try:
            backend = reduce_mod.add_hydrogens(str(model_path), str(h_pdb))
            status["reduce"] = "done ({})".format(backend)
        except Exception as exc:  # noqa: BLE001
            status["reduce"] = "error: {}".format(exc)
            logger.warning("reduce failed for %s: %s", pdb_id, exc)
            logger.debug(traceback.format_exc())

    # --- molprobity aggregate report ---
    molprobity_out = molprobity_dir / "molprobity.out"
    if skip_existing and molprobity_out.exists():
        status["molprobity"] = "skipped"
    else:
        try:
            mtz_arg = str(files["mtz"]) if files["mtz"] else None
            result = molprobity_mod.run_molprobity(str(model_path), mtz_path=mtz_arg)
            report = molprobity_mod.format_report(result, str(model_path))
            molprobity_out.write_text(report + "\n")
            status["molprobity"] = "done" + (
                "" if mtz_arg else " (model-only, no .mtz found)"
            )
        except Exception as exc:  # noqa: BLE001
            status["molprobity"] = "error: {}".format(exc)
            logger.warning("molprobity failed for %s: %s", pdb_id, exc)
            logger.debug(traceback.format_exc())

    return status


def run_pipeline(
    df, repo_path, path_column="path_in_repo", skip_existing=True, verbose=True
):
    """
    Run the validation pipeline for every entry in df.

    Args:
        df: DataFrame indexed by PDB ID, with a column (default
            'path_in_repo') giving the directory containing that
            entry's .pdb/.cif/.mtz files
        path_column: name of the column holding the directory path
        skip_existing: skip steps whose output already exists,
            matching the original bash scripts' caching behavior
        verbose: log a line per entry as it's processed

    Returns:
        a DataFrame indexed the same as df, with one column per
        pipeline step (xtriage, cablam, clashscore, reduce,
        molprobity) holding a status string for that step. Never
        raises for individual-entry failures -- check the returned
        DataFrame for 'error: ...' values to find problem entries.
    """
    if path_column not in df.columns:
        raise ValueError("df must have a '{}' column".format(path_column))

    records = {}
    total = len(df)
    for i, (pdb_id, row) in enumerate(df.iterrows(), start=1):
        if verbose:
            logger.info("[%d/%d] processing %s", i, total, pdb_id)
        try:
            records[pdb_id] = process_one(
                pdb_id, Path(repo_path / row[path_column]), skip_existing=skip_existing
            )
        except Exception as exc:  # noqa: BLE001
            # process_one is designed not to raise, but guard here too
            # so one catastrophic entry (e.g. a permissions error
            # creating the validation dir) never aborts the batch.
            logger.error("Unhandled failure for %s: %s", pdb_id, exc)
            logger.debug(traceback.format_exc())
            records[pdb_id] = {
                step: "error: unhandled failure - {}".format(exc) for step in STEP_NAMES
            }

    summary = pd.DataFrame.from_dict(records, orient="index")
    summary.index.name = df.index.name or "pdb_id"
    return summary[STEP_NAMES]


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Run the validation pipeline over a CSV of "
        "pdb_id,path_in_repo rows."
    )
    parser.add_argument(
        "csv", help="CSV with pdb_id and path_in_repo " "columns (pdb_id as index)"
    )
    parser.add_argument("--path-column", default="path_in_repo")
    parser.add_argument(
        "--no-skip-existing",
        action="store_true",
        help="Recompute steps even if output exists",
    )
    parser.add_argument("--out", default="validation_run_summary.csv")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        stream=sys.stderr,
    )

    df = pd.read_csv(args.csv, index_col=0)
    summary = run_pipeline(
        df,
        path_column=args.path_column,
        skip_existing=not args.no_skip_existing,
    )
    summary.to_csv(args.out)
    logger.info("Summary written to %s", args.out)

    n_errors = summary.apply(lambda col: col.str.startswith("error").sum()).sum()
    if n_errors:
        logger.warning(
            "%d step-level errors across the batch; see %s", n_errors, args.out
        )


if __name__ == "__main__":
    main()
