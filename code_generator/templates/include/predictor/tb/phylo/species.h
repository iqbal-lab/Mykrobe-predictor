{% extends 'include/predictor/common/phylo/species.h' %}

{% block extra %}
boolean is_MTBC_present(SpeciesInfo* species_info);
void update_phylo_group_presence_and_coverage_from_species(SpeciesInfo* species_info);
boolean* create_MTBC_mask();
boolean* create_NTM_mask();
Species get_best_MTBC_species(SpeciesInfo* species_info );
Species get_best_NTM_species(SpeciesInfo* species_info );
Lineage get_best_lineage(SpeciesInfo* species_info );
boolean MTBC_panels_are_present(SpeciesInfo* species_info);
boolean NTM_panels_are_present(SpeciesInfo* species_info);
boolean no_MTBC_panels_are_present(SpeciesInfo* species_info);
boolean no_NTM_panels_are_present(SpeciesInfo* species_info);
boolean no_lineage_panels_are_present(SpeciesInfo* species_info);


boolean myco_is_present(SpeciesInfo* species_info);
boolean tuberculosis_is_present(SpeciesInfo* species_info);

int get_contamination_covg(SpeciesInfo* species_info);


{% endblock %}