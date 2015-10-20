{% extends 'src/predictor/common/phylo/species.c' %}

{% block extra %}

int get_best_hit(CovgInfo* covg_info, boolean* mask)
{
  int i;
  int best_perc_cov_so_far=0;
  int best_median_cov_so_far=0;
  Species unknown_enum = unknownspecies;
  int curr=(int) unknown_enum;

  for (i=0; i<covg_info->NUM_PANELS; i++)
  {
    if (mask[i]){
      if (covg_info->percentage_coverage[i] >= best_perc_cov_so_far)
      {
         best_perc_cov_so_far = covg_info->percentage_coverage[i];
        // Only update if the median coverage has also improved
        if (covg_info->median_coverage[i] > best_median_cov_so_far){
          best_median_cov_so_far = covg_info->median_coverage[i];
          curr=i;
        }          
      }      
    }
  }
  return  curr;
}

boolean* create_staph_mask(){
  boolean* mask= create_mask(true);
  return (mask);
}

boolean* create_non_aureus_mask(){
  boolean* mask= create_mask(true);
  mask[Saureus] = false;
  return (mask);
}

boolean non_aureus_panels_are_present(SpeciesInfo* species_info){
  boolean* mask = create_non_aureus_mask();
  boolean non_aureus_species_panels_are_present = panels_are_present(species_info->species_covg_info,mask);
  // Finished with the mask so we need to free it. 
  free(mask);
  return (non_aureus_species_panels_are_present);
}

boolean no_non_aureus_panels_are_present(SpeciesInfo* species_info){
  return (!non_aureus_panels_are_present(species_info));
}

boolean staphylococcus_is_present(SpeciesInfo* species_info){
  boolean* mask = create_staph_mask();
  boolean staph_species_panels_are_present = panels_are_present(species_info->species_covg_info,mask);
  // Finished with the mask so we need to free it. 
  free(mask);  
  return (staph_species_panels_are_present);
}

Species get_best_staph_species(SpeciesInfo* species_info ){
  boolean* mask = create_staph_mask();
  int species_enum = get_best_hit(species_info->species_covg_info,mask);
  Species species = species_enum;
  // Finished with the mask so we need to free it. 
  free(mask);  
  return (species);
}

Species get_best_non_aureus_species(SpeciesInfo* species_info ){
  boolean* mask = create_non_aureus_mask();
  int species_enum  = get_best_hit(species_info->species_covg_info,mask);
  Species species = species_enum;
  // Finished with the mask so we need to free it. 
  free(mask);  
  return (species);
}


boolean is_aureus_present(SpeciesInfo* species_info){
  return (species_info->species_covg_info->present[Saureus]);
}

boolean is_non_aureus_staph_present(SpeciesInfo* species_info){
  boolean is_epi_present = species_info->species_covg_info->present[Sepidermidis];
  boolean is_haem_present = species_info->species_covg_info->present[Shaemolyticus];
  boolean is_sother_present = species_info->species_covg_info->present[Sother];
  if (is_epi_present || is_haem_present  || is_sother_present){
    return (true);
  }
  else{
    return (false);
  }
}

int get_contamination_covg(SpeciesInfo* species_info){
	
	boolean is_contamination_present = is_non_aureus_staph_present(species_info);
  int contamination_covg = 0;
	if (is_contamination_present){
		Species best_non_aureus_species = get_best_non_aureus_species(species_info);
    contamination_covg = species_info->species_covg_info->median_coverage[best_non_aureus_species];
	}
  printf("contamination_covg : %i\n", contamination_covg);
  return contamination_covg;
	
}
{% endblock %}