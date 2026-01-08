"""Module for querying and updating PDB structural data.

It contains utilities for:
- Retrieving PDB ids from an rcsb query that were released or revised in a
  specified date range.
- Retrieving Protein names for a list of pdb_ids through sequence similarity search
  with a specified FASTA file.
- Retrieving release date, revision date, version, experimental methods and title for ids
  and storing them in a dictionary.
- Updating / creating a dataframe from a dictionary containing ids and protein names.
  calls get_attributes to get data and merges with an existing dataframe if provided,
  while making sure that ids, that were previously assigned a protein manually, are
  not overwritten by 'not_assigned'.

Typical usage example:
    import pandas as pd
    from rcsbapi.search import search_attributes as attrs
    from pathlib import Path

    q1 = attrs.rcsb_entity_source_organism.taxonomy_lineage.id == "2697049"
    q2 = attrs.rcsb_entity_source_organism.taxonomy_lineage.id != "6645"
    query = q1 & q2
    start = '2023-03-02'
    end = '2023-04-14'
    taxonomy = "SARS-CoV-2"
    fasta_path = "../data/fasta/seq_SARS-CoV-2.fasta"
    repo_path = Path.cwd().parent / "data"
    df = pd.read_pickle(repo_path / "dataframes/repo_database_SARS-CoV-2_copy.pkl")

    ids = get_ids(start, end, query)
    proteins = get_proteins(ids, fasta_path)
    new_df = get_df(proteins, taxonomy, df)
"""

import re
import pandas as pd
from Bio import SeqIO
from typing import Iterable, TypedDict
from pathlib import Path
from datetime import date
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
from rcsbapi.search import search_attributes as attrs
from rcsbapi.search import SeqSimilarityQuery
from rcsbapi.search.search_query import Group
from rcsbapi.data import DataQuery
from .utils import edit_dict_via_file


def get_ids(start: str, end: str, rcsb_query: Group) -> list:
    """Retrieve PDB accession IDs within a specified date range.
    Args:
        start: Start date in ISO format (YYYY-MM-DD) for filtering.
        end: End date in ISO format (YYYY-MM-DD) for filtering.
        query: An rcsb-api Group object to filter and constrain.
        Check out rcsbapi.readthedocs.io for information on how to
        create custom queries.
    Returns:
        A list of lowercase PDB accession IDs matching the query criteria.
    A usage example is described in help(update_pipeline.query).
    """

    revised_query = (start <= attrs.rcsb_accession_info.revision_date) & (
        attrs.rcsb_accession_info.revision_date <= end
    )
    release_query = (start <= attrs.rcsb_accession_info.initial_release_date) & (
        attrs.rcsb_accession_info.initial_release_date <= end
    )

    query = rcsb_query & (revised_query | release_query)
    ids = [id.lower() for id in query()]

    return ids


def get_proteins(
    ids: Iterable, fasta_path: str | Path, repo_path=None
) -> dict[str, str]:
    """Retrieve and match protein sequences to PDB IDs using sequence similarity search.

    Args:
        ids: List of PDB IDs to search for matches.
        fasta_path: File path to a FASTA file containing protein sequences.
    Returns:
        Dictionary mapping each PDB ID to either a hyphen-separated string
        of matched protein names (sorted alphabetically) or "not_assigned"
        if no matches were found.
    """

    fasta = list(SeqIO.parse(fasta_path, "fasta"))
    result = {pdb_id: [] for pdb_id in ids}

    def search_seq(seq: str) -> set[str]:
        try:
            query = SeqSimilarityQuery(
                value=seq,
                evalue_cutoff=1,
                sequence_type="protein",
            )
        except Exception as e:
            print(f"Exception occured for {seq}: {e}")
            return set()

        return set(e[:4].lower() for e in query("polymer_entity"))

    # Run searches in parallel threads
    with ThreadPoolExecutor() as executor:
        hits_list = list(
            tqdm(
                executor.map(lambda rec: search_seq(str(rec.seq)), fasta),
                total=len(fasta),
                desc="Sequence similarity search",
            )
        )

    for rec, hits in zip(fasta, hits_list):
        prot_name = rec.description.split(" ")[1]
        brackets = re.search(r"\((.*?)\)", prot_name)
        if brackets:
            prot_name = brackets.group(1)
        for hit in hits:
            if hit in ids:
                result[hit].append(prot_name)

    return {
        pdb_id: ("-".join(sorted(names)) if names else "not_assigned")
        for pdb_id, names in result.items()
    }


def get_attributes(ids: Iterable) -> dict[str, dict]:
    """Retrieve and reformat structural attributes for PDB entries from RCSB PDB.
    Args:
        ids: A list of PDB entry IDs to query.
    Returns:
        A dictionary mapping PDB IDs (lowercase) to dictionaries containing the following keys:
            protein (str), path_in_repo (str), release_date (date), last_revision (date), version (float),
            exp_method (str or None), resolution (float) title (str),superseded_by (float): NaN.
        Returns an empty dictionary if the input list is empty or if the query yields no results.
    """

    if not ids:
        return {}
    query = DataQuery(
        input_type="entries",
        input_ids=list(ids),
        return_data_list=[
            "entries.rcsb_id",
            "entries.rcsb_accession_info.initial_release_date",
            "entries.rcsb_accession_info.revision_date",
            "entries.rcsb_accession_info.major_revision",
            "entries.rcsb_accession_info.minor_revision",
            "entries.exptl.method",
            "entries.rcsb_entry_info.resolution_combined",
            "struct.title",
        ],
    )
    query.exec()
    response = query.get_response()
    raw_entries = response.get("data", {}).get("entries", [])
    reformatted = {}

    for entry in raw_entries:
        pdb_id = entry.get("rcsb_id", "").lower()
        release_date = entry.get("rcsb_accession_info", {}).get("initial_release_date")
        revision_date = entry.get("rcsb_accession_info", {}).get("revision_date")
        major = entry.get("rcsb_accession_info", {}).get("major_revision")
        minor = entry.get("rcsb_accession_info", {}).get("minor_revision")
        version = f"{major}.{minor}" if minor else str(major)
        exp_method = None
        if isinstance(entry["exptl"], list) and entry["exptl"]:
            exp_method = "; ".join([e.get("method") for e in entry["exptl"]])
        resolution = entry.get("rcsb_entry_info", {}).get("resolution_combined")
        if isinstance(resolution, list) and resolution:
            resolution = resolution[0]
        title = entry.get("struct", {}).get("title")
        reformatted[pdb_id] = {
            "protein": float("NaN"),
            "path_in_repo": float("NaN"),
            "release_date": date.fromisoformat(release_date.split("T")[0]),
            "last_revision": date.fromisoformat(revision_date.split("T")[0]),
            "version": float(version),
            "exp_method": exp_method,
            "resolution": float(resolution) if resolution else float("NaN"),
            "title": title,
            "superseded_by": None,
        }
    return reformatted


def update_proteins(proteins: dict[str, str], df: pd.DataFrame = None) -> dict:
    """
    Update protein assignments in a dictionary based on a DataFrame and user input.

    Ids that are not assigned a protein are updated if that same id is already in
    the dataframe. If there are still not_assigned ids, it writes a JSON file called
    '"assign_manually_edit_me.json" to the current working directory. You can then edit
    the file manually, save it and the updated dict is returned. The JSON file is deleted.
    Args:
        proteins: dictionary mapping ids to proteins.
        df (optional): A DataFrame indexed by ID containing a "protein" column.
                       If None or empty, the ids are not updated with the df,
                       but a JSON file is created nonetheless.
    Returns:
        dict: The updated proteins dictionary with newly assigned protein names.
              Entries previously marked as "not_assigned" are updated where possible.
    """

    if df is not None and not df.empty:
        common_ids = list(set(proteins.keys()).intersection(df.index))
        proteins.update(
            df.loc[[i for i in common_ids if proteins[i] == "not_assigned"], "protein"]
        )

    not_assigned_dict = {
        id: protein for id, protein in proteins.items() if protein == "not_assigned"
    }

    if not_assigned_dict:
        newly_assigned = edit_dict_via_file(not_assigned_dict)
        proteins.update(newly_assigned)

    return proteins


def get_df(
    proteins: dict[str, str],
    taxonomy: str,
    df: pd.DataFrame = None,
):
    """Create a DataFrame with protein metadata and repository paths.

    Calls update_proteins to update not_assigned proteins using the dataframe (if provided)
    and user input via a JSON file called "assign_manually_edit_me.json" in the current working
    directory. Calls get_attributes to get data and merges with an existing dataframe (if provided).
    Args:
        proteins: dictionary mapping ids to proteins.
        taxonomy: Taxonomy identifier for id folder path.
        df (optional): Will be updated with the new data, if provided.
            Must be indexed by pdb_id and have the same following columns:
            protein, release_date, last_revision, version, exp_method,
            path_in_repo, resolution, title, superseded_by.
    """

    attributes = get_attributes(proteins.keys())
    proteins = update_proteins(proteins, df)

    attributes_df = pd.DataFrame.from_dict(attributes, orient="index")
    proteins_df = pd.DataFrame.from_dict(proteins, orient="index", columns=["protein"])

    new_df = attributes_df.combine_first(proteins_df)
    new_df["path_in_repo"] = (
        "pdb/" + new_df["protein"].astype(str) + f"/{taxonomy}/" + new_df.index
    )
    cols = [
        "protein",
        "release_date",
        "last_revision",
        "version",
        "exp_method",
        "path_in_repo",
        "resolution",
        "title",
        "superseded_by",
    ]
    new_df = new_df[cols]
    new_df.index.rename("pdb_id", inplace=True)
    new_df["superseded_by"] = new_df["superseded_by"].astype(object)

    if df is None or df.empty:
        return new_df

    return new_df.combine_first(df)
