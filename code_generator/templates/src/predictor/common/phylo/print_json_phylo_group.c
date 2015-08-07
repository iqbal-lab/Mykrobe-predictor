
void print_json_{{phylogroup.name}}(SpeciesInfo* species_info){
    CovgInfo* covg_info =species_info->{{phylogroup.name}}_covg_info;    
    int num_panels_present = covg_info->num_panels_present;
    print_json_{{phylogroup.name}}_start();
    if (num_panels_present > 0){
      print_json_indiv_phylo(covg_info,get_ith_{{phylogroup.name}}_name);
    }
    else
    {
      print_json_called_variant_item( "Unknown {{phylogroup.name}}", -1, true);
    }
    print_json_phylo_group_end();  
}
