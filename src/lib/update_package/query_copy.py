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

import time
import random
from typing import Tuple, Set

import re
import pandas as pd
import numpy as np
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


def get_ids(start: str, end: str, rcsb_queries: dict) -> list:
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

    ids_by_taxonomy = {}
    for taxonomy, rcsb_query in rcsb_queries.items():
        revised_query = (start <= attrs.rcsb_accession_info.revision_date) & (
            attrs.rcsb_accession_info.revision_date <= end
        )
        release_query = (start <= attrs.rcsb_accession_info.initial_release_date) & (
            attrs.rcsb_accession_info.initial_release_date <= end
        )

        query = rcsb_query & (revised_query | release_query)
        ids = {entity.lower() for entity in query("polymer_entity")}
        for id in ids:
            if id in ids_by_taxonomy:
                ids_by_taxonomy[id].add(taxonomy)
            else:
                ids_by_taxonomy[id] = {taxonomy}
    return {id: "__".join(sorted(taxonomy)) for id, taxonomy in ids_by_taxonomy.items()}


def get_polymer_type_and_sequences(entity_ids):
    """Helper function used in get_proteins.

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


def get_proteins(
    ids_by_taxonomy: dict, fasta_path: str | Path, workers=100
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

    # List to accomodate multiple sequences for the same protein
    protein_sequences = []
    for protein in list(SeqIO.parse(fasta_path, "fasta")):
        if "OS" in protein.description:
            protein_name = protein.description[
                protein.description.find(" ") : protein.description.find("OS")
            ].strip()
        else:
            protein_name = protein.description.split(" ")[1]
            brackets = re.search(r"\((.*?)\)", protein_name)
            if brackets:
                protein_name = brackets.group(1)
        protein_name = protein_name.replace("-", "_").replace(" ", "_")
        protein_sequences.append((protein_name, str(protein.seq)))

    MAX_RETRIES = 3
    BASE_DELAY = 1.0  # seconds
    JITTER = 0.5  # add randomness to avoid thundering herd

    def _run_query_with_retry(query, protein_name: str) -> Set[str]:
        """
        Execute an RCSB query with retry and exponential backoff.
        """
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                # Execute query
                results = query("polymer_entity")
                return {rid.lower() for rid in results}

            except Exception as e:
                if attempt == MAX_RETRIES:
                    print(
                        f"[ERROR] {protein_name}: failed after {MAX_RETRIES} attempts: {e}"
                    )
                    return set()

                delay = BASE_DELAY * (2 ** (attempt - 1))
                delay += random.uniform(0, JITTER)

                print(
                    f"[WARN] {protein_name}: attempt {attempt} failed ({e}); "
                    f"retrying in {delay:.2f}s"
                )
                time.sleep(delay)

        return set()

    def search_seq(item: Tuple[str, str]) -> Tuple[str, Set[str]]:
        """
        Perform sequence similarity or motif search with retry handling.
        """
        protein_name, seq = item

        if not seq:
            return protein_name, set()

        try:
            if len(seq) >= 25:
                query = SeqSimilarityQuery(
                    value=seq,
                    evalue_cutoff=1,
                    sequence_type="protein",
                )
            else:
                query = SeqMotifQuery(value=seq)

            hits = _run_query_with_retry(query, protein_name)
            return protein_name, hits

        except Exception as e:
            # This catches construction-time errors (very rare)
            print(f"[ERROR] {protein_name}: query setup failed: {e}")
            return protein_name, set()

    def search_concurrently(sequences, workers: int = 8):
        """
        Run sequence searches concurrently with progress tracking.
        """
        similar_ids: dict[str, Set[str]] = {}

        with ThreadPoolExecutor(max_workers=workers) as executor:
            future_to_name = {
                executor.submit(search_seq, item): item[0] for item in sequences
            }

            for future in tqdm(
                as_completed(future_to_name),
                total=len(sequences),
                desc="Sequence similarity search",
            ):
                protein_name, hits = future.result()

                if protein_name in similar_ids:
                    similar_ids[protein_name].update(hits)
                else:
                    similar_ids[protein_name] = hits

        return similar_ids

    protein_hits = search_concurrently(protein_sequences, workers=workers)

    ids_by_protein = {id: "" for id in ids_by_taxonomy}

    for id in ids_by_protein:
        for protein_name, hits in protein_hits.items():
            if id in hits:
                ids_by_protein[id] = protein_name

    na_ids = {id for id, protein_name in ids_by_protein.items() if not protein_name}

    if na_ids:
        na_ids_by_polymer_types, na_ids_by_sequences = get_polymer_type_and_sequences(
            na_ids
        )

        # Assign rna if polymer type is rna
        for id, polymer_type in na_ids_by_polymer_types.items():
            if (polymer_type != "protein") and id in ids_by_protein:
                ids_by_protein[id] = polymer_type
                na_ids_by_sequences.pop(id)

        entity_hits = search_concurrently(na_ids_by_sequences.items())

        max_proportion = {entity_id: 0 for entity_id in na_ids}
        for entity_id, similar_ids in entity_hits.items():
            for protein_name, reference_ids in protein_hits.items():
                if similar_ids:
                    current_proportion = len(
                        similar_ids.intersection(reference_ids)
                    ) / len(similar_ids)
                    if current_proportion > max_proportion[entity_id]:
                        max_proportion[entity_id] = current_proportion
                        ids_by_protein[entity_id] = protein_name

    result = {
        id: {"taxonomy": ids_by_taxonomy[id], "protein": ids_by_protein[id]}
        for id in ids_by_taxonomy
    }
    return result


def retrieve_nested_attribute(data, rcsb_data_path: str):
    """Helper function used by get_df

    Args:
        data (_type_): _description_
        rcsb_data_path (str): _description_

    Returns:
        _type_: _description_
    """
    paths = rcsb_data_path.split(".")
    if "polymer_entities" in paths:
        paths.remove("polymer_entities")
    results = set()

    def retrieve_single_attribute(data, paths):
        if isinstance(data, list):
            for d in data:
                retrieve_single_attribute(d, paths.copy())

        elif isinstance(data, dict):
            if not paths:
                return

            attribute = paths[0]
            new_data = data[attribute]
            if (len(paths) == 1) and new_data:
                if isinstance(new_data, list):
                    new_data = new_data[0]
                if new_data not in results:
                    results.add(str(new_data))
            else:
                retrieve_single_attribute(new_data, paths[1:])

    retrieve_single_attribute(data, paths)
    return "__".join(sorted(results))


def get_df(
    proteins_and_taxonomy: dict,
    rcsb_data_attributes: dict = {},
    functions_to_combine_columns: dict = {},
    old_df: pd.DataFrame | None = None,
    aggregate: bool = True,
) -> dict[str, dict]:
    """Retrieve and reformat structural attributes for PDB entries from RCSB PDB.
    Args:
        ids: A list of PDB entry IDs to query.
    Returns:
        A dictionary mapping PDB IDs (lowercase) to dictionaries containing the following keys:
            protein (str), path_in_repo (str), release_date (date), last_revision (date), version (float),
            exp_method (str or None), resolution (float) title (str),superseded_by (float): NaN.
        Returns an empty dictionary if the input list is empty or if the query yields no results.
    """

    ids = proteins_and_taxonomy.keys()
    if not ids:
        return {}

    rcsb_data_attributes = {
        "entry": "entry.entry.id",
        "release_date": "entry.rcsb_accession_info.initial_release_date",
        "last_revised": "entry.rcsb_accession_info.revision_date",
        **rcsb_data_attributes,
    }

    query = DataQuery(
        input_type="polymer_entities",
        input_ids=list(ids),
        return_data_list=list(rcsb_data_attributes.values()),
    )
    query.exec()
    response = query.get_response()
    raw_entries = response["data"]["polymer_entities"]

    attrs = {}
    for entry in raw_entries:

        pdb_id = entry["rcsb_id"].lower()
        attrs[pdb_id] = {
            attribute_name: retrieve_nested_attribute(entry, rcsb_data_path)
            for attribute_name, rcsb_data_path in rcsb_data_attributes.items()
        }

    df = pd.DataFrame.from_dict(attrs, orient="index")
    df = pd.concat(
        [pd.DataFrame.from_dict(proteins_and_taxonomy, orient="index"), df], axis=1
    )

    df.entry = [id.lower() for id in df.entry]

    if aggregate:

        def aggregate_entities(s):
            entity_attribute = "-".join(
                sorted(set().union(*[set(d.split("__")) for d in s]))
            )
            if entity_attribute:
                if entity_attribute[0] == "-":
                    entity_attribute = entity_attribute[1:]
            return entity_attribute

        df = df.groupby("entry").agg(aggregate_entities)

    else:
        df.replace({"__": "-"}, regex=True)

    for column in df.columns:
        try:
            df[column] = pd.to_numeric(df[column], downcast="integer")
        except:
            ValueError

    df.release_date = pd.to_datetime(df.release_date).dt.date
    df.last_revised = pd.to_datetime(df.last_revised).dt.date

    for columns, f in functions_to_combine_columns.items():
        if "=" in columns:
            new_column, old_columns = columns.split("=")
        else:
            new_column, old_columns = columns, ""
        df[new_column.strip()] = df.apply(f, axis=1)
        if old_columns:
            old_columns = [c.strip() for c in old_columns.split("+")]
            print(f"dropping columns {', '.join(old_columns)}")
            df = df.drop(old_columns, axis=1)

    df = df.replace("", float("Nan"))

    if old_df is not None:
        df = df.combine_first(old_df)

    return df
