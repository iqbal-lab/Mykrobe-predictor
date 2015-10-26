{% extends 'src/predictor/common/phylo/species.c' %}

{% block extra %}

boolean* create_MTBC_mask()
{
  boolean* mask= create_mask(false);
  mask[Tuberculosis] = true;
  mask[Bovis] = true;
  mask[Caprae] = true;
  mask[Africanum] = true;
  return (mask);
}

boolean* create_NTM_mask()
{
  boolean* mask= create_mask(true);
  mask[Tuberculosis] = false;
  mask[Bovis] = false;
  mask[Caprae] = false;
  mask[Africanum] = false;
  return (mask);
}

Species get_best_MTBC_species(SpeciesInfo* species_info ){
  boolean* mask = create_MTBC_mask();
  int species_enum = get_best_hit(species_info->species_covg_info,mask);
  // Finished with the mask so we need to free it. 
  free(mask);    
  Species species = species_enum;
  return (species);
}

Species get_best_NTM_species(SpeciesInfo* species_info ){
  boolean* mask = create_NTM_mask();
  int species_enum  = get_best_hit(species_info->species_covg_info,mask);
  // Finished with the mask so we need to free it. 
  free(mask);    
  Species species = species_enum;
  return (species);
}

Lineage get_best_lineage(SpeciesInfo* species_info ){
  boolean* mask= create_mask(true);
  int lineage_enum  = get_best_hit(species_info->lineage_covg_info,mask);
  // Finished with the mask so we need to free it. 
  free(mask);    
  Lineage lineage = lineage_enum;
  return (lineage);
}

boolean MTBC_panels_are_present(SpeciesInfo* species_info){
  boolean* mask = create_MTBC_mask();
  boolean MTBC_species_panels_are_present = panels_are_present(species_info->species_covg_info,mask);
  // Finished with the mask so we need to free it. 
  free(mask);    
  return (MTBC_species_panels_are_present);
}
boolean NTM_panels_are_present(SpeciesInfo* species_info){
  boolean* mask = create_NTM_mask();
  boolean NTM_species_panels_are_present = panels_are_present(species_info->species_covg_info,mask);
  // Finished with the mask so we need to free it. 
  free(mask);    
  return (NTM_species_panels_are_present);  
}

boolean no_MTBC_panels_are_present(SpeciesInfo* species_info){
  return (!MTBC_panels_are_present(species_info));
}
boolean no_NTM_panels_are_present(SpeciesInfo* species_info){
  return (!NTM_panels_are_present(species_info));  
}

boolean no_lineage_panels_are_present(SpeciesInfo* species_info){
  if (species_info->lineage_covg_info->num_panels_present >0 ){
    return (false);
  }
  else{
    return (true);
  }
}

boolean tuberculosis_is_present(SpeciesInfo* species_info){
  return (species_info->species_covg_info->present[Tuberculosis]);
}
boolean myco_is_present(SpeciesInfo* species_info){
  boolean MTBC_is_present = species_info->complex_covg_info->present[Mtbc];
  boolean NTM_is_present = species_info->complex_covg_info->present[Ntm];
  return (MTBC_is_present || NTM_is_present);
}

int get_contamination_covg(SpeciesInfo* species_info){
	
  boolean is_contamination_present = NTM_panels_are_present(species_info);
  int contamination_covg = 0;
	if (is_contamination_present){
    contamination_covg = get_best_NTM_species(species_info);
	}
  return contamination_covg;
	
}

int get_expected_covg(SpeciesInfo* species_info){
  int expected_covg = species_info->complex_covg_info->median_coverage[Mtbc];
  return expected_covg;
}
{% endblock %}