import asyncio
from datetime import date, timedelta
import functools
import pandas as pd


def async_wrapper(func):
    """wrapper for async functions"""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return asyncio.run(func(*args, **kwargs))

    return wrapper


def get_time(day: str = None, next_week: bool = False) -> date:
    """
    Get the date of the most recent Wednesday, or optionally the Wednesday of the next week.
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
