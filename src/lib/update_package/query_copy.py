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
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from rcsbapi.search import search_attributes as attrs
from rcsbapi.search import SeqSimilarityQuery, SeqMotifQuery
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
    ids = {entity.lower() for entity in query("polymer_entity")}
    return ids


def get_polymer_type_and_sequences(entity_ids):
    """_summary_

    Args:
        entity_ids (_type_): _description_
    """
    query = DataQuery(
        input_type="polymer_entities",
        input_ids=list(entity_ids),
        return_data_list=[
            "entity_poly.rcsb_entity_polymer_type",
            "entity_poly.pdbx_seq_one_letter_code_can",
        ],
    )
    result = query.exec()

    polymer_types = {}
    sequences = {}
    for d in result["data"]["polymer_entities"]:
        id = d["rcsb_id"].lower()
        polymer_types[id] = d["entity_poly"]["rcsb_entity_polymer_type"].lower()
        sequences[id] = d["entity_poly"]["pdbx_seq_one_letter_code_can"]

    return polymer_types, sequences


def get_proteins(ids: Iterable, fasta_path: str | Path, workers=100) -> dict[str, str]:
    """Retrieve and match protein sequences to PDB IDs using sequence similarity search.

    Args:
        ids: List of PDB IDs to search for matches.
        fasta_path: File path to a FASTA file containing protein sequences.
    Returns:
        Dictionary mapping each PDB ID to either a hyphen-separated string
        of matched protein names (sorted alphabetically) or "not_assigned"
        if no matches were found.
    """

    # List to accomodate multiple sequences for the same protein
    protein_sequences = []
    for protein in list(SeqIO.parse(fasta_path, "fasta")):
        protein_name = protein.description.split(" ")[1]
        brackets = re.search(r"\((.*?)\)", protein_name)
        if brackets:
            protein_name = brackets.group(1)
        protein_sequences.append((protein_name, str(protein.seq)))

    def search_seq(item: tuple) -> tuple[str, set[str]]:
        protein_name, seq = item
        try:
            if len(seq) >= 25:
                query = SeqSimilarityQuery(
                    value=seq,
                    evalue_cutoff=1,
                    sequence_type="protein",
                )
            else:
                query = SeqMotifQuery(value=seq)

            return protein_name, set(id.lower() for id in query("polymer_entity"))

        except Exception as e:
            print(f"Exception occurred for {protein_name}: {e}")
            return protein_name, set()

    def search_concurently(sequences):
        similar_ids = {}
        with ThreadPoolExecutor(max_workers=workers) as executor:
            # Submit all tasks
            future_to_id = {
                executor.submit(search_seq, item): item[0] for item in sequences
            }

            # Process completed tasks with progress bar
            for future in tqdm(
                as_completed(future_to_id),
                total=len(sequences),
                desc="Sequence similarity search",
            ):
                key, hits = future.result()
                if key in similar_ids:
                    similar_ids[key].update(hits)
                else:
                    similar_ids[key] = hits
        return similar_ids

    protein_hits = search_concurently(protein_sequences)

    result = {id: "not_assigned" for id in ids}
    for protein_name, hits in protein_hits.items():
        for id in ids:
            if id in hits:
                result[id] = protein_name

    not_assigned_entities = {
        id for id, name in result.items() if name == "not_assigned"
    }

    na_polymer_types, na_sequences = get_polymer_type_and_sequences(
        not_assigned_entities
    )

    for entity, polymer_type in na_polymer_types.items():
        if polymer_type != "protein":
            result[entity] = polymer_type
            na_sequences.pop(entity)

    entity_hits = search_concurently(na_sequences.items())

    max_proportion = {entity_id: 0 for entity_id in not_assigned_entities}
    for entity_id, similar_ids in entity_hits.items():
        for protein_name, reference_ids in protein_hits.items():
            current_proportion = len(similar_ids.intersection(reference_ids)) / len(
                similar_ids
            )
            if current_proportion > max_proportion[entity_id]:
                result[entity_id] = protein_name
                max_proportion[entity_id] = current_proportion

    return {
        pdb_id: names if names else "not_assigned" for pdb_id, names in result.items()
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
