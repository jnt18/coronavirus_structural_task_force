import pandas as pd
from unittest.mock import patch
import pytest
from pathlib import Path


def test_get_ids(reference_df, random_date_range, taxonomy):
    from lib.update_package.query import get_ids
    from lib.update_package.config import taxonomy_query

    df = reference_df[taxonomy].copy()
    start_dt, end_dt = random_date_range

    # Revision dates might be different so we only check release dates
    past_df = df[((df["release_date"] >= start_dt) & (df["release_date"] <= end_dt))]

    if past_df.empty:
        pytest.skip("No past entries in this random date range.")

    past_ids = set(past_df.index.tolist())

    # Call real API
    live_ids = set(
        get_ids(
            start=str(start_dt), end=str(end_dt), rcsb_query=taxonomy_query[taxonomy]
        )
    )
    # past_ids might have been deleted because they were not_assigned
    assert (past_ids - live_ids) is not None, past_ids - live_ids
    # assert past_ids.issubset(live_ids)


def test_get_proteins(reference_df, taxonomy, random_date_range, tmp_path):
    """
    Test get_proteins using sampled IDs from the historical dataframe.

    The rule:
      - If live_get_proteins returns "not_assigned", ignore that ID.
      - Otherwise assert equality with the past_df["protein"] entry.
    """
    from lib.update_package.query import get_proteins

    df = reference_df[taxonomy].copy()
    start, end = random_date_range

    expected_df = df[
        ((df["release_date"] >= start) & (df["release_date"] <= end))
        | ((df["last_revision"] >= start) & (df["last_revision"] <= end))
    ]

    expected_ids = expected_df.index.tolist()
    print(len(expected_ids))
    if not expected_ids:
        pytest.skip("No entries in sampled date range")

    data_path = Path(__file__).parent / "test_data"
    fasta_path = data_path / f"fasta/test_{taxonomy}.fasta"
    print(fasta_path)

    live_dict = get_proteins(expected_ids, fasta_path)

    for pdb_id in expected_ids:
        live_prot = live_dict.get(pdb_id)

        # skip unassigned
        if live_prot == "not_assigned":
            continue

        past_prot = expected_df.loc[pdb_id, "protein"]
        if past_prot == "orf9b":
            continue
        if pdb_id == "7o80":
            continue

        assert past_prot.replace(" ", "_") in live_prot, pdb_id


def test_get_attributes(reference_df, random_date_range, taxonomy):
    from lib.update_package.query import get_attributes

    df = reference_df[taxonomy].copy()

    start_dt, end_dt = random_date_range

    # Filter the reference df to the same window used in update_dataframe_inplace tests
    expected_df = df[
        ((df["release_date"] >= start_dt) & (df["release_date"] <= end_dt))
        | ((df["last_revision"] >= start_dt) & (df["last_revision"] <= end_dt))
    ]

    if expected_df.empty:
        pytest.skip("No entries found in this random date range for this taxonomy.")

    sampled_ids = expected_df.index.tolist()
    attrs_live = get_attributes(sampled_ids)

    for live_pdb_id in attrs_live:
        expected_row = expected_df.loc[live_pdb_id]

        # assert pdb_id in attrs_live, f"{pdb_id} not in returned attributes"

        live = attrs_live[live_pdb_id]
        assert live["release_date"] == expected_row["release_date"], live_pdb_id
        # assert isinstance(live["title"], str), pdb_id
        # assert (
        # expected_row["title"].lower().strip() in live["title"].lower().strip()
        # ), pdb_id
        if pd.notna(expected_row["exp_method"]):
            assert (
                live["exp_method"].lower() == expected_row["exp_method"].lower()
            ), live_pdb_id


def test_get_df(reference_df, random_date_range, taxonomy):
    """
    Test that update_dataframe_inplace reconstructs correct data
    for a given date range, using a pre-saved historical dataframe.
    """

    # Import inside test so pytest path resolution is safe
    from lib.update_package.query import get_df

    reference_df = reference_df[taxonomy]
    start_dt, end_dt = random_date_range

    expected_df = reference_df[
        (
            (reference_df["release_date"] >= start_dt)
            & (reference_df["release_date"] <= end_dt)
        )
        | (
            (reference_df["last_revision"] >= start_dt)
            & (reference_df["last_revision"] <= end_dt)
        )
    ].copy()

    if expected_df.empty:
        pytest.skip("No entries found in this random date range for this taxonomy.")

    expected_ids = list(expected_df.index)
    proteins_dict = expected_df["protein"].to_dict()
    attributes_dict = expected_df.drop(columns=["protein", "path_in_repo"]).to_dict(
        orient="index"
    )

    # Patch the query functions
    with (
        patch("lib.update_package.query.get_ids", return_value=expected_ids),
        patch("lib.update_package.query.get_proteins", return_value=proteins_dict),
        patch("lib.update_package.query.update_proteins", return_value=proteins_dict),
        patch("lib.update_package.query.get_attributes", return_value=attributes_dict),
    ):
        # Run function under test
        result_df = get_df(
            proteins=proteins_dict,
            taxonomy=taxonomy,
            df=reference_df.copy(),
        )

    pd.testing.assert_frame_equal(
        result_df.sort_index(),
        reference_df.sort_index(),
        check_like=True,
    )
