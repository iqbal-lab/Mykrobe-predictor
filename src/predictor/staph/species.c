/*
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
  species.c
*/
// system headers
#include <stdlib.h>
#include <limits.h>

// third party headers
#include <string_buffer.h>

#include "build.h" 
#include "maths.h" 
#include "element.h"
#include "seq.h"
#include "open_hash/hash_table.h"
#include "dB_graph.h"
#include "species.h"
#include "gene_presence.h"
#include "genotyping_known.h"


  void map_phylo_group_enum_to_str(Phylo_Group sp, StrBuf* sbuf)
{
  
     if(sp==Cat){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "cat");
      }
  
     else if(sp==Coagneg){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "coagneg");
      }
  
  else
    {
      die("Coding error - I would expect the compiler to prevent assigning a bad enum value %i \n",sp );
    }
}
  void load_all_phylo_group_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir )
{
	
  panel_file_paths[Cat] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Cat], "data/staph/phylo/phylo_group/cat.fasta" );
  
  panel_file_paths[Coagneg] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Coagneg], "data/staph/phylo/phylo_group/coag_neg.fasta" );
  
}
  void phylo_group_threshold(int* thresholds){
	
	  thresholds[Cat] = 50;
	  thresholds[Cat] = 50;
	
	  thresholds[Coagneg] = 50;
	  thresholds[Coagneg] = 50;
	
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
      print_json_called_variant_item( "Unknown phylo_group", -1, true);
    }
    print_json_phylo_group_end();  
}
  char* get_ith_phylo_group_name(CovgInfo* covg_info, int i)
{
  Phylo_Group phylo_group;
  StrBuf* phylo_group_name = strbuf_new(); 
  phylo_group = get_ith_present_panel( covg_info, i);
  map_phylo_group_enum_to_str(phylo_group, phylo_group_name);
  return phylo_group_name->buff;
}

  void map_species_enum_to_str(Species sp, StrBuf* sbuf)
{
  
     if(sp==Saureus){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "Saureus");
      }
  
     else if(sp==Sepidermidis){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "Sepidermidis");
      }
  
     else if(sp==Shaemolyticus){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "Shaemolyticus");
      }
  
     else if(sp==Sother){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "Sother");
      }
  
  else
    {
      die("Coding error - I would expect the compiler to prevent assigning a bad enum value %i \n",sp );
    }
}
  void load_all_species_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir )
{
	
  panel_file_paths[Saureus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Saureus], "data/staph/phylo/species/Saureus.fasta" );
  
  panel_file_paths[Sepidermidis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Sepidermidis], "data/staph/phylo/species/Sepidermidis.fasta" );
  
  panel_file_paths[Shaemolyticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shaemolyticus], "data/staph/phylo/species/Shaemolyticus.fasta" );
  
  panel_file_paths[Sother] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Sother], "data/staph/phylo/species/Sother.fasta" );
  
}
  void species_threshold(int* thresholds){
	
	  thresholds[Saureus] = 90;
	  thresholds[Saureus] = 90;
	
	  thresholds[Sepidermidis] = 50;
	  thresholds[Sepidermidis] = 50;
	
	  thresholds[Shaemolyticus] = 50;
	  thresholds[Shaemolyticus] = 50;
	
	  thresholds[Sother] = 50;
	  thresholds[Sother] = 50;
	
}
  
void print_json_species(SpeciesInfo* species_info){
    CovgInfo* covg_info =species_info->species_covg_info;    
    int num_panels_present = covg_info->num_panels_present;
    print_json_species_start();
    if (num_panels_present > 0){
      print_json_indiv_phylo(covg_info,get_ith_species_name);
    }
    else
    {
      print_json_called_variant_item( "Unknown species", -1, true);
    }
    print_json_phylo_group_end();  
}
  char* get_ith_species_name(CovgInfo* covg_info, int i)
{
  Species species;
  StrBuf* species_name = strbuf_new(); 
  species = get_ith_present_panel( covg_info, i);
  map_species_enum_to_str(species, species_name);
  return species_name->buff;
}




SpeciesInfo* get_species_info(dBGraph *db_graph,int max_branch_len, 
                            StrBuf* install_dir,int expected_covg,
                            int ignore_first,int ignore_last)

{

  char* get_ith_phylo_group_name(CovgInfo* covg_info, int i);
  StrBuf* phylo_group_file_paths[NUM_Phylo_Group];
  load_all_phylo_group_file_paths(phylo_group_file_paths,install_dir);
  CovgInfo* phylo_group_covg_info = get_coverage_info(db_graph,
                                                  phylo_group_file_paths,
                                                  max_branch_len,NUM_Phylo_Group,
                                                  ignore_first,ignore_last,
                                                  phylo_group_threshold);  


  char* get_ith_species_name(CovgInfo* covg_info, int i);
  StrBuf* species_file_paths[NUM_Species];
  load_all_species_file_paths(species_file_paths,install_dir);
  CovgInfo* species_covg_info = get_coverage_info(db_graph,
                                                  species_file_paths,
                                                  max_branch_len,NUM_Species,
                                                  ignore_first,ignore_last,
                                                  species_threshold);  




  SpeciesInfo* species_info=(SpeciesInfo *)malloc(sizeof(SpeciesInfo)); 
    
  species_info->phylo_group_covg_info = phylo_group_covg_info;
    
  species_info->species_covg_info = species_covg_info;
  

  // update_phylo_group_presence_and_coverage_from_species(species_info);

  return species_info;
}




void print_json_phylo_group_start()
{
  printf("\t\t\"phylo_group\": {\n");
}


void print_json_species_start()
{
  printf("\t\t\"species\": {\n");
}


void print_json_phylogenetics(SpeciesInfo* species_info){
    print_json_phylogenetics_start();

  
    print_json_phylo_group(species_info);
  
    print_json_species(species_info);
  
    print_json_phylogenetics_end();  
}

