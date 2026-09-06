#!/usr/bin/env python3

"""
Set up CCTBX chem_data for a uv-managed environment.

The actual chem_data directory lives in the project, rather than inside
the virtual environment.  A symlink is created in the location where the
installed CCTBX/libtbx environment searches for chem_data.

This script assumes the project contains:

    chem_data/
        cablam_data/
        geostd/
        rotarama_data/
        ...

It creates:

    .venv/lib/pythonX.Y/site-packages/chem_data
        -> <project>/chem_data
"""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CHEM_DATA = PROJECT_ROOT / "chem_data"


def main() -> int:
    print(f"Project root: {PROJECT_ROOT}")
    print(f"chem_data:   {CHEM_DATA}")

    if not CHEM_DATA.is_dir():
        print("\nERROR: chem_data directory does not exist:\n" f"  {CHEM_DATA}\n")
        print("Obtain the version-matched CCTBX chem_data package first.")
        return 1

    cablam_data = CHEM_DATA / "cablam_data"

    if not cablam_data.is_dir():
        print("\nERROR: chem_data/cablam_data does not exist:\n" f"  {cablam_data}\n")
        return 1

    required_stat_files = [
        "cablam.8000.expected.general.stat",
        "cablam.8000.expected.gly.stat",
        "cablam.8000.expected.transpro.stat",
        "cablam.8000.expected.cispro.stat",
    ]

    missing = [
        name for name in required_stat_files if not (cablam_data / name).is_file()
    ]

    if missing:
        print("\nERROR: missing CaBLAM source tables:")
        for name in missing:
            print(f"  {name}")
        return 1

    # Locate the Python installation used by this script.
    site_packages = Path(next(p for p in sys.path if p.endswith("site-packages")))

    target = site_packages / "chem_data"

    print(f"Python:      {sys.executable}")
    print(f"site-packages: {site_packages}")
    print(f"CCTBX target:  {target}")

    if target.is_symlink():
        existing = target.resolve()
        if existing == CHEM_DATA.resolve():
            print("\nchem_data symlink already correct.")
            return 0

        print("\nRemoving existing chem_data symlink:\n" f"  {target} -> {existing}")
        target.unlink()

    elif target.exists():
        print(
            "\nERROR: target already exists and is not a symlink:\n"
            f"  {target}\n\n"
            "Refusing to delete it automatically."
        )
        return 1

    target.parent.mkdir(parents=True, exist_ok=True)

    # Use an absolute symlink so it remains unambiguous.
    target.symlink_to(CHEM_DATA.resolve(), target_is_directory=True)

    print("\nCreated:")
    print(f"  {target} -> {target.resolve()}")

    # Verify using libtbx itself.
    print("\nTesting libtbx chem_data lookup...")

    import libtbx.load_env
    import libtbx

    found = libtbx.env.find_in_repositories(
        relative_path="chem_data/cablam_data",
        test=os.path.isdir,
    )

    print(f"libtbx found: {found}")

    if found is None:
        print("\nERROR: libtbx still cannot find chem_data.")
        return 1

    print("\nchem_data setup successful.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
