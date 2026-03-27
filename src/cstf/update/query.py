"""Utilities for querying PDB structural data.

This module provides helpers to:
- Retrieve PDB polymer entity IDs released or revised within a given date range
  using RCSB search queries.
- Assign protein names to PDB entities via sequence similarity searches against
  a reference FASTA file.
- Fetch metadata (e.g. release dates, revisions, experimental methods, titles)
  for PDB entities and normalize nested RCSB responses.
- Create or update a pandas DataFrame containing PDB metadata, optionally
  merging with an existing dataset while preserving manually assigned proteins.

Typical usage:
    >>> import pandas as pd
    >>> from pathlib import Path

    >>> taxonomies = ["H1N1", "H3N2", "H5N1", "H5N8"]
    >>> rcsb_queries = {k: v for k, v in config.taxonomy_query.items() if k in taxonomies}

    >>> start = "2023-03-02"
    >>> end = "2023-04-14"
    >>> fasta_path = "../data/fasta/seq_SARS-CoV-2.fasta"
    >>> repo_path = Path.cwd().parent / "data"

    >>> old_df = pd.read_pickle(repo_path / "dataframes/repo_database_SARS-CoV-2_copy.pkl")

    >>> attributes = {k: v for k, v in config.rcsb_data_attributes.items() if k in ["version_1", "version_2", "exp_method", "resolution", "title"]}

    Additional data attributes can be explored via:

    >>> from rcsbapi.data import DataSchema
    >>> DataSchema().find_field_names(string)

    >>> functions = {k: v for k, v in config.functions_to_combine_columns.items() if k in ["version=version_1+version_2", "path_in_repo", "exp_method", "supersed_by"]}

    >>> ids = get_ids(start, end, rcsb_queries)
    >>> proteins = get_proteins(ids, fasta_path)
    >>> new_df = get_df(proteins, attributes, functions, old_df)
"""

import pandas as pd
import numpy as np
from typing import Iterable
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from rcsbapi.search import search_attributes as attrs
from rcsbapi.search import SeqSimilarityQuery, SeqMotifQuery
from rcsbapi.data import DataQuery
from .utils import (
    get_protein_seq_from_fastas,
    retry_with_backoff,
    retrieve_nested_attribute,
)


def get_ids(start: str, end: str, rcsb_queries: dict) -> dict[str:str]:
    """Retrieve PDB polymer entity IDs released or revised within a date range.

    For each supplied RCSB query (typically representing a taxonomy or category),
    this function finds polymer entities whose initial release date or most recent
    revision date falls within the specified interval.

    Args:
        start: Start date (inclusive), ISO format: YYYY-MM-DD.
        end: End date (inclusive), ISO format: YYYY-MM-DD.
        rcsb_queries: Mapping of taxonomy name to an RCSB search query object.

    Returns:
        Dictionary mapping lowercase polymer entity IDs to a taxonomy string.
        If an entity matches multiple taxonomies, they are joined using "__".
    """
    # Track which taxonomies each entity ID matches
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
            ids_by_taxonomy.setdefault(id, set()).add(taxonomy)
    return {
        id: {"taxonomy": "__".join(sorted(taxonomy))}
        for id, taxonomy in ids_by_taxonomy.items()
    }


def get_attributes(ids: Iterable, attributes_to_fetch: dict):
    """Fetch metadata for PDB polymer entities using rcsbapi.data.

    This helper executes a bulk DataQuery and flattens nested JSON responses
    into a dictionary suitable for DataFrame construction.

    Args:
        ids: Iterable of polymer entity IDs (e.g. from get_ids).
        attributes_to_fetch: Mapping of output column name to RCSB schema path.
            Available paths for a given string can be explored via:
            from rcsbapi.data import DataSchema
            DataSchema().find_field_names(string)

    Returns:
        Dictionary keyed by lowercase polymer entity ID, where each value is
        a dictionary of fetched attributes. Can be converted to a DataFrame via:
            pd.DataFrame.from_dict(result, orient="index")
    """
    query = DataQuery(
        input_type="polymer_entities",
        input_ids=list(ids),
        return_data_list=list(attributes_to_fetch.values()),
    )
    query.exec()

    # Extract and structure the response data
    raw_entries = query.get_response()["data"]["polymer_entities"]

    return {
        entry["rcsb_id"].lower(): {
            attr_name: retrieve_nested_attribute(entry, path)
            for attr_name, path in attributes_to_fetch.items()
        }
        for entry in raw_entries
    }


def get_proteins(
    ids_by_taxonomy: dict,
    fasta_paths: list[str | Path] = None,
    workers: int = 10,
) -> dict[str, dict]:
    """
    Assign protein names to PDB polymer entities using sequence similarity search.

    Protein sequences from the provided FASTA file are queried against RCSB using
    sequence similarity (or motif search for short sequences). Entities matching
    a reference sequence are assigned the corresponding protein name.

    For entities that remain unassigned:
    - Non-protein polymers are labeled by their polymer type (e.g. "rna", "dna").
    - Protein entities are compared via entity-to-entity similarity, and assigned
      the protein with the highest overlap to reference hits.

    Args:
        ids_by_taxonomy: Mapping of polymer entity ID to taxonomy label
            (as returned by get_ids).
        fasta_path: Path to FASTA file containing reference protein sequences.
        workers: Maximum number of threads for concurrent RCSB queries.

    Returns:
        Dictionary mapping polymer entity ID to a dictionary with:
            - "taxonomy": taxonomy label(s)
            - "protein": assigned protein name (empty string if unresolved)
    """

    # Parse FASTA sequences
    protein_sequences = get_protein_seq_from_fastas(fasta_paths)

    @retry_with_backoff
    def run_rcsb_query(query) -> set[str]:
        """
        Execute an RCSB sequence query with retry and exponential backoff.
        """
        results = query("polymer_entity")
        return {id.lower() for id in results}

    # Sequence search function
    def search_seq(item: tuple[str, str]) -> tuple[str, set[str]]:
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
            hits = run_rcsb_query(query)
            return protein_name, hits
        except Exception as e:
            print(f"[ERROR] {protein_name}: query failed: {e}")
            return protein_name, set()

    def run_concurrent_tasks(items: Iterable, desc: str) -> dict[str, set[str]]:
        results: dict[str, set[str]] = {}
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(search_seq, item): item[0] for item in items}
            for key, hits in tqdm(
                (future.result() for future in as_completed(futures)),
                total=len(futures),
                desc=desc,
            ):
                results.setdefault(key, set()).update(hits)
        return results

    # Run sequence similarity searches for all reference proteins
    protein_hits = run_concurrent_tasks(protein_sequences, "Sequence similarity search")
    # Initialize protein assignments (empty string = unassigned)
    ids_by_protein = {id: "" for id in ids_by_taxonomy}
    for id in ids_by_protein:
        for protein_name, hits in protein_hits.items():
            if id in hits:
                ids_by_protein[id] = protein_name

    # Handle unassigned entities
    not_assigned_ids = {id for id, name in ids_by_protein.items() if not name}
    if not_assigned_ids:
        attributes_to_fetch = {
            "polymer_type": "entity_poly.rcsb_entity_polymer_type",
            "seq": "entity_poly.pdbx_seq_one_letter_code_can",
        }
        attributes = get_attributes(not_assigned_ids, attributes_to_fetch)
        not_assigned_polymer_types = {
            id: d["polymer_type"].lower() for id, d in attributes.items()
        }
        not_assigned_sequences = {id: d["seq"] for id, d in attributes.items()}
        # Assign non-protein polymer types directly
        for id, polymer_type in not_assigned_polymer_types.items():
            if polymer_type != "protein" and id in ids_by_protein:
                ids_by_protein[id] = polymer_type
                not_assigned_sequences.pop(id)
        # Entity-to-entity similarity search
        entity_hits = run_concurrent_tasks(
            list(not_assigned_sequences.items()),
            "Entity similarity search",
        )
        max_proportion = {id: 0 for id in not_assigned_ids}
        for id, similar_ids in entity_hits.items():
            for protein_name, reference_ids in protein_hits.items():
                if similar_ids:
                    proportion = len(similar_ids.intersection(reference_ids)) / len(
                        similar_ids
                    )
                    if proportion > max_proportion[id]:
                        max_proportion[id] = proportion
                        ids_by_protein[id] = protein_name

    return {
        id: {
            **ids_by_taxonomy[id],
            "protein": ids_by_protein[id],
        }
        for id in ids_by_taxonomy
    }


def get_df(
    proteins_and_taxonomy: dict[str:dict],
    rcsb_data_attributes: dict | None = None,
    functions_to_combine_columns: dict | None = None,
    old_df: pd.DataFrame | None = None,
    aggregate: bool = True,
) -> pd.DataFrame:
    """
    Build or update a DataFrame containing PDB metadata and protein annotations.

    This function fetches RCSB metadata for polymer entities, merges it with
    taxonomy and protein assignments, and optionally aggregates multiple
    entities belonging to the same PDB entry.

    Args:
        proteins_and_taxonomy: Mapping of polymer entity ID to metadata such as
            taxonomy and protein name (from get_proteins).
        rcsb_data_attributes: Optional mapping of additional RCSB attributes
            to fetch. If omitted, only release and revision dates are included.
        functions_to_combine_columns: Optional mapping of column specifications
            to transformation functions. Keys may encode source columns using
            the format "new_col=col1+col2" where col1 and col2 would be dropped.
        old_df: Existing DataFrame to merge with. Existing values take precedence.
        aggregate: Whether to aggregate multiple polymer entities per PDB entry.

    Returns:
        A pandas DataFrame containing PDB metadata, protein annotations,
        and any derived columns.
    """
    if not proteins_and_taxonomy:
        return pd.DataFrame()

    # Define default attributes with any custom overrides
    attributes_to_fetch = {
        "entry": "entry.entry.id",
        "release_date": "entry.rcsb_accession_info.initial_release_date",
        "last_revised": "entry.rcsb_accession_info.revision_date",
        **(rcsb_data_attributes or {}),
    }

    extracted_attributes = get_attributes(
        proteins_and_taxonomy.keys(), attributes_to_fetch
    )

    # Combine input data with fetched attributes
    df = pd.concat(
        [
            pd.DataFrame.from_dict(proteins_and_taxonomy, orient="index"),
            pd.DataFrame.from_dict(extracted_attributes, orient="index"),
        ],
        axis=1,
    )

    df["entry"] = df["entry"].str.lower()
    df.sort_index(inplace=True)

    # not_date_cols = df.columns.difference(["release_date", "last_revised"])
    # df[not_date_cols] = df[not_date_cols].replace({"-": "_"}, regex=True)

    # Optionally aggregate multiple entities per entry
    if aggregate:

        def aggregate_values(series):
            """Combine multiple values, removing duplicates"""
            unique_values = {}
            for value in series:
                if not value:
                    continue
                for v in value.split("__"):
                    unique_values[v] = None  # keeps insertion order
            return "-".join(unique_values)
            # return "-".join(
            #     dict.fromkeys(v for value in series if value for v in value.split("__"))
            # )

        df = df.groupby("entry", as_index=True).agg(aggregate_values)
    else:
        df.replace({"__": "-"}, regex=True, inplace=True)

    # Convert numeric columns to appropriate types
    for column in df.columns:
        try:
            df[column] = pd.to_numeric(df[column], downcast="integer")
        except (ValueError, TypeError):
            pass  # Keep as string if conversion fails

    # Convert date columns to datetime
    for date_column in ("release_date", "last_revised"):
        df[date_column] = pd.to_datetime(df[date_column], errors="coerce").dt.date

    # Apply custom column transformations
    for column_spec, transform_func in (functions_to_combine_columns or {}).items():
        # Parse "new_column = old_col1 + old_col2" format
        new_column_name, *source_columns = column_spec.split("=")
        df[new_column_name.strip()] = df.apply(transform_func, axis=1)

        # Drop source columns if specified
        if source_columns:
            columns_to_drop = [col.strip() for col in source_columns[0].split("+")]
            df.drop(columns=columns_to_drop, inplace=True)

    # Replace empty strings with NaN
    df.replace("", float("nan"), inplace=True)

    # Merge with existing DataFrame if provided
    if old_df is not None:
        df = df.combine_first(old_df)

    return df
