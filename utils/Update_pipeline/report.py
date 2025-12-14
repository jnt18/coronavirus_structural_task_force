"""Module for writing weekly reports and last update report.

Each report lists newly released and revised structures,
with new structures additionally grouped by protein.
Make sure that the dataframe is up-to-date using functions in the query module.
Typical usage:
    start = '2023-03-02'
    end = '2023-04-14'
    repo_path = Path.cwd().parent / "data"
    df = pd.read_pickle(repo / "dataframes/repo_database_SARS-CoV-2_copy.pkl")
    taxonomy = "SARS-CoV-2"

    write_reports(start, end, df, taxonomy, repo_path)
"""

from pathlib import Path
from itertools import groupby
import pandas as pd
import numpy as np
from .utils import get_time


def write_reports(
    start: str, end: str, df: pd.DataFrame, taxonomy: str, repo_path: str | Path
) -> None:
    """Generates weekly reports in the date range and a report summarising the full period.

    For the final report it calls write_single_update with start, end (as date types) and
    overwrites weekly_reports/latest_update_report.
    For weekly reports it calls write_single_update with start = end = d for all Wednesdays in the date range.
    Each report lists newly released and revised structures, with new structures additionally grouped by protein.
    Args:
        start: Start date for filtering structures in the report.
        end: End date for filtering structures in the report.
        df: Make sure that the dataframe is up-to-date using functions in the query module.
        taxonomy: Taxonomy identifier for report header.
        repo_path: Reports are written to folder called "weekly_reports" in this directory.
    """
    repo_path = Path(repo_path)
    reports_path = repo_path / "weekly_reports"
    reports_path.mkdir(parents=True, exist_ok=True)

    # We start on the Wednesday after the first date
    start_dt = get_time(start, next_week=True)
    # We end on the Wednesday after the second date.
    end_dt = get_time(end)

    # Latest report
    report_path = reports_path / f"latest_update_report_{taxonomy}.txt"
    write_single_report(start_dt, end_dt, taxonomy, report_path, df, latest_report=True)

    dates = pd.date_range(
        start_dt, end_dt, periods=(end_dt - start_dt).days // 7 + 1
    ).date

    for day in dates:
        df_date = df[(df.release_date == day) | (df.last_revision == day)]
        report_path = reports_path / f"{day}_update_report_{taxonomy}.txt"
        write_single_report(day, day, taxonomy, report_path, df_date)


def write_single_report(
    start,
    end,
    taxonomy,
    report_path: Path,
    df: pd.DataFrame,
    latest_report: bool = False,
) -> None:
    """Write a single report file containing newly released and revised structures.

    Generates a formatted text report that lists revised structures, new structures,
    and new structures grouped by protein type. The report is written to the specified file path.
    Args:
        start, end: Must be date objects and Wednesdays!
        taxonomy: Taxonomy identifier for report header.
        report_path: File path where the report will be written.
        df: DataFrame containing ids that will be written into the report, no matter their release dates.
        latest_report: If True, formats header as "start until end" for latest report.
            If False, formats header as "end weekly" for weekly reports. Defaults to False.
    """
    new_ids = set(df[df.release_date.between(start, end)].index)
    revised_ids = set(df.index) - new_ids

    date_header = f"{start} until {end}" if latest_report else f"{end} weekly"

    with report_path.open("w") as doc:
        doc.write(f"{date_header} report for {taxonomy}\n")
        doc.write(f"##### {len(revised_ids)} revised structures #####\n")
        doc.write(", ".join(sorted(revised_ids)) + "\n\n")
        doc.write(f"##### {len(new_ids)} new structures #####\n")
        doc.write(", ".join(sorted(new_ids)) + "\n\n")
        doc.write("##### new structures by protein #####")
        for protein, ids in groupby(new_ids, key=lambda k: df.loc[k, "protein"]):
            doc.write(f"\n{protein}\n>" + " ".join(ids))
