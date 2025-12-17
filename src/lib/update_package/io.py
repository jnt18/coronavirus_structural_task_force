"""This module provides functions for downloading structure files and managing superseded entries.

It includes utilities for:
- Downloading structure files (cif, pdb, mtz) from RCSB PDB with concurrent request limiting
  while archiving previous versions of downloaded files with timestamps.
- Removing directories of superseded PDB entries and updating their status in the dataframe.

Typical usage example:
    get ids and dataframe using update_pipeline.query.
    taxonomy = "SARS-CoV-2"
    repo_path = Path.cwd().parent / "data"

    download_files(ids, df, repo_path)
    delete_superseded(ids, df, repo_path)
"""

import asyncio
import aiohttp
import shutil
from pathlib import Path
from tqdm.asyncio import tqdm_asyncio
import pandas as pd
from typing import Iterable
from .utils import async_wrapper


@async_wrapper
async def download_files(ids: Iterable, new_df: pd.DataFrame, repo_path: str) -> None:
    """Download structure files for all PDB entries in the given DataFrame.

    Uses asyncio for concurrent downloads and aiohttp for async HTTP requests.
    When used in a jupyter notebook do
    Can be run without the wrapper using asyncio.run(download_files.__wrapper__(*args)).
    Args:
        ids: Files will be downloaded for these ids if they are in the dataframe.
        df: Make sure that the dataframe is up-to-date using functions in the query module.
        repo_path: Base directory path where downloaded files will be stored
    """

    sem = asyncio.Semaphore(20)
    connector = aiohttp.TCPConnector(limit=50)

    repo_path = Path(repo_path)  # ensure repo_path is a Path object

    async with aiohttp.ClientSession(connector=connector) as session:

        tasks = []
        for pdb_id in ids:
            directory = repo_path / Path(new_df.loc[pdb_id, "path_in_repo"])
            timestamp = new_df.loc[pdb_id, "release_date"]

            for ext in ("cif", "pdb", "mtz"):
                tasks.append(
                    download_files_handle_file(
                        session=session,
                        pdb_id=pdb_id,
                        file_dir=directory,
                        ext=ext,
                        timestamp=timestamp,
                        sem=sem,
                    )
                )

        await tqdm_asyncio.gather(*tasks)


async def download_files_handle_file(
    session: aiohttp.ClientSession,
    pdb_id: str,
    file_dir: Path,
    ext: str,
    timestamp: str,
    sem: asyncio.Semaphore,
):
    """Download a file from RCSB PDB and archive existing file with a timestamp.
    Args:
        session: An aiohttp ClientSession object for making async HTTP requests.
        pdb_id: The 4-character PDB identifier.
        file_dir: Path object specifying the directory where the file will be saved.
        ext: File extension (e.g., 'pdb', 'cif') indicating the file format to download.
        timestamp: Timestamp string appended to archived filenames for version tracking.
        sem: asyncio.Semaphore object to limit concurrent downloads.
    """
    async with sem:
        file_dir.mkdir(parents=True, exist_ok=True)

        file_path = file_dir / f"{pdb_id}.{ext}"

        # Archive old file if it exists
        if file_path.exists():
            old_dir = file_dir / "old"
            old_dir.mkdir(exist_ok=True)

            archived_name = f"{pdb_id}_{timestamp}.{ext}"
            archived_path = old_dir / archived_name

            file_path.replace(archived_path)

        # Download the new file
        url = f"https://files.rcsb.org/download/{pdb_id}.{ext}"
        async with session.get(url) as resp:
            data = await resp.read()

        file_path.write_bytes(data)


def delete_superseded(ids: Iterable, new_df: pd.DataFrame, repo_path: str) -> None:
    """Delete directories of superseded PDB entries and update their status in the dataframe inplace.
    Args:
        ids: List of PDB IDs to process.
        df: Pandas DataFrame containing PDB entry metadata with columns including "path_in_repo"
            and "superseded_by". Index should contain PDB IDs.
        repo_path: Path to the repository root directory containing PDB structure folders.
    Raises:
        FileNotFoundError: If a specified PDB file cannot be read (non-fatal, continues processing).
    """
    repo_path = Path(repo_path)

    for pdb_id in ids:
        id_path = repo_path / Path(new_df.loc[pdb_id, "path_in_repo"])
        pdb_path = id_path / f"{pdb_id}.pdb"

        if not pdb_path.exists():
            print(f"File {pdb_path} not found.")
            continue

        with pdb_path.open("r") as f:
            for line in f:
                if line.startswith("SPRSDE"):
                    # lines have the form: SPRSDE ... old_ids ...
                    old_ids = line.split()[3:]
                    for old_id in old_ids:
                        old_id = old_id.lower()
                        if old_id in new_df.index:
                            # mark superseded
                            new_df.loc[old_id, "superseded_by"] = pdb_id
                            # remove the old file folder if it exists
                            old_path = repo_path / Path(
                                new_df.loc[old_id, "path_in_repo"]
                            )
                            if old_path.exists() and old_path.is_dir():
                                shutil.rmtree(old_path)
                                print(
                                    f"deleted subdirectory for {old_id} as it was superseded by {pdb_id}"
                                )
