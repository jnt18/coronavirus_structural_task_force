from rcsbapi.search import search_attributes as attrs

taxonomy_query = {
    "SARS-CoV": (attrs.rcsb_entity_source_organism.taxonomy_lineage.id == "694009")
    & (attrs.rcsb_entity_source_organism.taxonomy_lineage.id != "2697049"),
    "SARS-CoV-2": (attrs.rcsb_entity_source_organism.taxonomy_lineage.id == "2697049")
    & (attrs.rcsb_entity_source_organism.taxonomy_lineage.id != "6645"),
}
