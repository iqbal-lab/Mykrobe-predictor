{% extends 'src/predictor/common/phylo/species.c' %}

{% block extra %}

void print_json_lineage_start()
{
  printf("\t\t\"lineage\": {\n");
}


char* get_char_name_of_species_enum(Species species){
  StrBuf* species_name = strbuf_new(); 
  map_species_enum_to_str(species, species_name);
  return species_name->buff;
}


void update_phylo_group_presence_and_coverage_from_species(SpeciesInfo* species_info){
  if (non_aureus_panels_are_present(species_info)){
    if (! species_info->phylo_group_covg_info->present[Coagneg]){
      species_info->phylo_group_covg_info->present[Coagneg] = true;
      species_info->phylo_group_covg_info->num_panels_present = species_info->phylo_group_covg_info->num_panels_present + 1;
    }
    Species best_staph_species = get_best_non_aureus_species(species_info);
    species_info->phylo_group_covg_info->percentage_coverage[Coagneg] = max(species_info->phylo_group_covg_info->percentage_coverage[Coagneg] , species_info->species_covg_info->percentage_coverage[best_staph_species] );
    species_info->phylo_group_covg_info->median_coverage[Coagneg] = max(species_info->phylo_group_covg_info->median_coverage[Coagneg] , species_info->species_covg_info->median_coverage[best_staph_species] );
  }
}

void load_all_cat_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir )
{
  panel_file_paths[0] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[0], "data/staph/species/coag_neg.fasta" ); 
}

void cat_threshold(int* thresholds){
  thresholds[0] = 20;
}

SpeciesInfo* get_species_info(dBGraph *db_graph,int max_branch_len, 
                            StrBuf* install_dir,int expected_covg,
                            int ignore_first,int ignore_last)

{
  StrBuf* phylo_group_file_paths[NUM_Phylo_Group];
  load_all_phylo_group_file_paths(phylo_group_file_paths,install_dir);
  StrBuf* species_file_paths[NUM_Species];
  load_all_species_file_paths(species_file_paths,install_dir);


  CovgInfo* phylo_group_covg_info = get_coverage_info(db_graph,
                                                  phylo_group_file_paths,
                                                  max_branch_len,NUM_Phylo_Group,
                                                  ignore_first,ignore_last,
                                                  phylo_group_threshold);
  CovgInfo* species_covg_info = get_coverage_info(db_graph,
                                                  species_file_paths,
                                                  max_branch_len,NUM_Species,
                                                  ignore_first,ignore_last,
                                                  species_threshold);

  StrBuf* cat_file_paths[1];
  load_all_cat_file_paths(cat_file_paths,install_dir);    
  CovgInfo* cat_covg_info = get_coverage_info(db_graph,
                                            cat_file_paths,
                                            max_branch_len,1,
                                            ignore_first,ignore_last,
                                            cat_threshold);   



  SpeciesInfo* species_info=(SpeciesInfo *)malloc(sizeof(SpeciesInfo)); 
  species_info->phylo_group_covg_info = phylo_group_covg_info;
  species_info->species_covg_info = species_covg_info;
  species_info->other_covg_info = cat_covg_info;

  update_phylo_group_presence_and_coverage_from_species(species_info);

  return species_info;
}

void print_json_aureus(SpeciesInfo* species_info, boolean last){
    print_json_called_variant_item( get_char_name_of_species_enum (Saureus) ,species_info->species_covg_info->median_coverage[Saureus], last);
}

void print_json_best_hit_non_aureus(SpeciesInfo* species_info){
  if (no_non_aureus_panels_are_present(species_info)){
    print_json_called_variant_item( "Unknown Species", -1 , true);
  }
  else{
  Species staph_species = get_best_non_aureus_species(species_info);
  print_json_called_variant_item( get_char_name_of_species_enum(staph_species), species_info->species_covg_info->median_coverage[staph_species], true);    
  }
}

void print_json_aureus_and_best_hit_non_aureus(SpeciesInfo* species_info){
  if (is_aureus_present(species_info)){
    print_json_aureus(species_info,false);
  }
  else
  {  
    print_json_called_variant_item( "Unknown Species", -1 , false);
  }
  if (no_non_aureus_panels_are_present(species_info)){
    print_json_called_variant_item( "Unknown Species", -1 , true);
  }
  else
  {  
    Species staph_species = get_best_non_aureus_species(species_info);  
    print_json_called_variant_item( get_char_name_of_species_enum(staph_species), species_info->species_covg_info->median_coverage[staph_species], true);    
  }
}

boolean catalayse_exists_in_sample(SpeciesInfo* species_info){    
    if (species_info->other_covg_info->percentage_coverage[0] > 20)
    {
      return true;
    }else
    {
      return false;
    }
}
int get_coverage_on_catalayse(SpeciesInfo* species_info){    
  return(species_info->other_covg_info->median_coverage[0]);
}

void print_json_phylo_group(SpeciesInfo* species_info){
    CovgInfo* covg_info =species_info->phylo_group_covg_info;    
    int num_panels_present = covg_info->num_panels_present;
    print_json_phylo_group_start();
    if (num_panels_present > 0){
      print_json_indiv_phylo(covg_info,get_ith_phylo_group_name);
    }
    else
    {
      if (catalayse_exists_in_sample(species_info)){
        print_json_called_variant_item( "Coagulase-Negative Staphylococcus", get_coverage_on_catalayse(species_info), true);
      }
      else{
        print_json_called_variant_item( "Non Staphylococcus", -1, true);
      }
    }
    print_json_phylo_group_end(false);  
}

void print_json_species(SpeciesInfo* species_info){
    Species aureus_is_present = is_aureus_present(species_info);
    Species non_aureus_staph_is_present = is_non_aureus_staph_present(species_info);
    print_json_species_start();
    if (aureus_is_present && non_aureus_staph_is_present){
      print_json_aureus_and_best_hit_non_aureus(species_info);
    }
    else if (aureus_is_present){
      print_json_aureus(species_info,true);
    }
    else if (non_aureus_staph_is_present){
      print_json_best_hit_non_aureus(species_info);
    }
    else
    {
      print_json_called_variant_item( "Unknown Species", -1, true);
    }    
    print_json_phylo_group_end(false);  
}

void print_json_lineage(SpeciesInfo* species_info){
    print_json_lineage_start();
    print_json_called_variant_item( "N/A", -1, true);
    print_json_phylo_group_end(true); 
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
	
	boolean is_contamination_present = species_info->phylo_group_covg_info->present[Coagneg];
  int contamination_covg = 0;
	if (is_contamination_present){
    contamination_covg = species_info->phylo_group_covg_info->median_coverage[Coagneg];
	}
  return contamination_covg;
	
}

int get_expected_covg(SpeciesInfo* species_info){
  
  int expected_covg = species_info->species_covg_info->median_coverage[Saureus];
  return expected_covg;
  
}

{% endblock %}