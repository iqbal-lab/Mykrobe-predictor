{% extends 'include/predictor/common/phylo/species.h' %}


{% block extra %}
char* get_char_name_of_species_enum(Species species);
int get_best_hit(CovgInfo* covg_info, boolean* mask);
boolean* create_staph_mask();
boolean* create_non_aureus_mask();
boolean non_aureus_panels_are_present(SpeciesInfo* species_info);
boolean no_non_aureus_panels_are_present(SpeciesInfo* species_info);
boolean staphylococcus_is_present(SpeciesInfo* species_info);
Species get_best_staph_species(SpeciesInfo* species_info );
Species get_best_non_aureus_species(SpeciesInfo* species_info );
boolean is_aureus_present(SpeciesInfo* species_info);
boolean is_non_aureus_staph_present(SpeciesInfo* species_info);
int get_contamination_covg(SpeciesInfo* species_info);
{% endblock %}