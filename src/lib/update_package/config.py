from rcsbapi.search import search_attributes as attrs
import pandas as pd

taxonomy_query = {
    "SARS-CoV": (
        attrs.rcsb_entity_source_organism.taxonomy_lineage.name
        == "Severe acute respiratory syndrome-related coronavirus"
    )
    & (
        attrs.rcsb_entity_source_organism.taxonomy_lineage.name
        != "Severe acute respiratory syndrome coronavirus 2"
    ),
    "SARS-CoV-2": (
        attrs.rcsb_entity_source_organism.taxonomy_lineage.name
        == "Severe acute respiratory syndrome coronavirus 2"
    )
    & (attrs.rcsb_entity_source_organism.taxonomy_lineage.name != "Octopus"),
    "H1N1": (attrs.rcsb_entity_source_organism.taxonomy_lineage.name == "H1N1"),
    "H3N2": (attrs.rcsb_entity_source_organism.taxonomy_lineage.name == "H3N2"),
    "H5N1": (attrs.rcsb_entity_source_organism.taxonomy_lineage.name == "H5N1"),
    "H5N8": (attrs.rcsb_entity_source_organism.taxonomy_lineage.name == "H5N8"),
}


rcsb_data_attributes = {
    "version_1": "entry.rcsb_accession_info.major_revision",
    "version_2": "entry.rcsb_accession_info.minor_revision",
    "exp_method": "entry.exptl.method",
    "resolution": "entry.rcsb_entry_info.resolution_combined",
    "title": "entry.struct.title",
    "gene_1": "uniprots.rcsb_uniprot_protein.gene.name.value",
    "gene_2": "rcsb_entity_source_organism.rcsb_gene_name.value",
    "common_names": "rcsb_polymer_entity.rcsb_macromolecular_names_combined.name",
    "scientific_name": "polymer_entities.rcsb_entity_source_organism.scientific_name",
}


functions_to_combine_columns = {
    "version=version_1+version_2": lambda x: float(
        f"{str(x.version_1)}.{str(0 if pd.isna(x.version_2) else int(x.version_2))}"
    ),
    "gene=gene_1+gene_2": lambda x: "-".join(
        sorted((set(x.gene_1.split("-")) | set(x.gene_2.split("-"))))
    ),
    "gene": lambda x: x.gene[1:] if x.gene.startswith("-") else x.gene,
    "path_in_repo": lambda x: f"pdb/{x.taxonomy}/{x.protein}/{x.name}",
    "exp_method": lambda x: "; ".join(
        x.exp_method.replace("X-RAY", "X_RAY").split("-")
    ),
    "superseded_by": lambda x: float("NaN"),
}
