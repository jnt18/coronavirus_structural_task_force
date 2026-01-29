from rcsbapi.search import search_attributes as attrs

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
