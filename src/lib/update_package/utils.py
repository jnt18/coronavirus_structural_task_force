"""Utility functions for updating packages in the coronavirus structural task force project.

Functions:
    async_wrapper: Decorator to run async functions synchronously
    get_time: Get the most recent Wednesday or optionally next week's Wednesday
    reset_index: Reset DataFrame index if indexed by pdb_id
    edit_dict_via_file: Interactively edit a dictionary through manual file editing"""

import asyncio
from datetime import date, timedelta
import json
from pathlib import Path
import functools
import pandas as pd

import time
import random
import functools
from typing import Callable, Any
from contextlib import contextmanager
from concurrent.futures import ThreadPoolExecutor


def retry_with_backoff(func):
    """
    Retry decorator with exponential backoff.

    Retries up to 3 times with base delay 1.0s and jitter 0.5s.
    """

    MAX_RETRIES = 3
    BASE_DELAY = 1.0
    JITTER = 0.5

    @functools.wraps(func)
    def wrapper(*args, **kwargs) -> Any:
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                return func(*args, **kwargs)

            except Exception as e:
                if attempt == MAX_RETRIES:
                    print(f"[ERROR] failed after {MAX_RETRIES} attempts: {e}")
                    raise

                delay = BASE_DELAY * (2 ** (attempt - 1))
                delay += random.uniform(0, JITTER)

                print(
                    f"[WARN] attempt {attempt} failed ({e}); "
                    f"retrying in {delay:.2f}s"
                )

                time.sleep(delay)

    return wrapper


def async_wrapper(func):
    """wrapper for async functions"""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return asyncio.run(func(*args, **kwargs))

    return wrapper


def retrieve_nested_attribute(data, rcsb_data_path: str) -> str:
    """
    Extract values from nested data structures following a dot-notation path.

    Args:
        data: The data structure to search (dict or list)
        rcsb_data_path: Dot-separated path (e.g., "entry.rcsb_accession_info.date")

    Returns:
        String of unique values joined by "__", or empty string if none found
    """
    # Remove "polymer_entities" from path segments
    path_segments = [s for s in rcsb_data_path.split(".") if s != "polymer_entities"]
    found_values = set()

    def traverse(current_obj, remaining_segments):
        """Recursively walk through nested data structure."""
        if not remaining_segments or current_obj is None:
            return

        current_key = remaining_segments[0]
        remaining = remaining_segments[1:]

        # Handle lists by recursively checking each item
        if isinstance(current_obj, list):
            for item in current_obj:
                traverse(item, remaining_segments)

        # Handle dictionaries
        elif isinstance(current_obj, dict):
            value = current_obj.get(current_key)

            # If we've reached the end of the path, collect the value
            if not remaining and value:
                # Extract first item from lists
                if isinstance(value, list):
                    value = value[0]
                found_values.add(str(value))
            else:
                # Continue traversing deeper
                traverse(value, remaining)

    traverse(data, path_segments)
    return "__".join(sorted(found_values))


def get_time(day: str = None, next_week: bool = False) -> date:
    """Get the date of the most recent Wednesday, or optionally the Wednesday of the next week.

    If the passed in day is a Wednesday it returns it.
    Args:
        day (optional): ISO format date string (YYYY-MM-DD). Defaults to today's date if not provided.
        next_week (optional): If True, returns the Wednesday of the following week instead of the current week.
                              Defaults to False.
    Returns:
        date: A date object representing Wednesday of the specified week.
              If the input day is already a Wednesday, returns that day
              (or next week's Wednesday if next_week is True).
    Examples:
        >>> get_time("2024-01-10")  # Returns Wednesday of that week
        >>> get_time(next_week=True)  # Returns next week's Wednesday
    """

    day = date.fromisoformat(day) if day else date.today()
    offset = (day.weekday() - 2) % 7  # Wednesday = 2
    if offset == 0:
        return day
    return day - timedelta(offset) + next_week * timedelta(7)


def reset_index(df: pd.DataFrame) -> pd.DataFrame:
    """Used in tcp_manual, tcp_scan_results and analyse_and_fix_dataframe
    while transitioning to having the dataframe be indexed by pdb_ids.

    Args:
        df: Usually read from data/dataframes

    Returns:
        The same dataframe with index reset to integers.
    """
    if df.index.name == "pdb_id":
        return df.reset_index()
    return df


def edit_dict_via_file(data: dict) -> dict:
    """Edit a dictionary through manual file editing.

    Creates a JSON file in the current working directory that the user can edit.
    After editing and saving the file, the function reads the updated content back
    and returns it as a dictionary. The created file is deleted.
    Args:
        data: The dictionary to be edited.
    Returns:
        dict: The updated dictionary after manual editing.

    """

    path = Path.cwd() / "assign_manually_edit_me.json"

    # Prevent overwriting
    if path.exists():
        raise FileExistsError(f"File already exists: {path}")

    # Write JSON
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=4)

    print(f"\nFile created: {path}")
    print("Edit the file, save it, then press ENTER here to continue.")
    print(
        "You can also copy the dictionary into a text editor, escape the function,\n"
        "and copy it back in when you run the function again."
    )
    input()

    # Read edited JSON
    with open(path, "r", encoding="utf-8") as f:
        updated = json.load(f)

    if path.is_file():
        path.unlink()
    else:
        raise RuntimeError(f"Refusing to delete non-file path: {path}")

    return updated
