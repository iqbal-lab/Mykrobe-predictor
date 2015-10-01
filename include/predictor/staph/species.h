/*
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 *  species.h
*/
#include "dB_graph.h"
#include "base_species.h"



void print_json_phylo_group_start();
void print_json_phylo_group_end();
char* get_ith_phylo_group_name(CovgInfo* covg_info, int i);

	typedef enum 
	 {
	 	
	 	Cat = 0,
	 	
	 	Coagneg = 1,
	 	
    unknownphylo_group=2
	   	} Phylo_Group ;
	#define NUM_Phylo_Group 2
   	
  void map_phylo_group_enum_to_str(Phylo_Group sp, StrBuf* sbuf);
  void load_all_phylo_group_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir );
  void phylo_group_threshold(int* thresholds);


void print_json_species_start();
void print_json_species_end();
char* get_ith_species_name(CovgInfo* covg_info, int i);

	typedef enum 
	 {
	 	
	 	Saureus = 0,
	 	
	 	Sepidermidis = 1,
	 	
	 	Shaemolyticus = 2,
	 	
	 	Sother = 3,
	 	
    unknownspecies=4
	   	} Species ;
	#define NUM_Species 4
   	
  void map_species_enum_to_str(Species sp, StrBuf* sbuf);
  void load_all_species_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir );
  void species_threshold(int* thresholds);


typedef struct
{

  CovgInfo* phylo_group_covg_info;

  CovgInfo* species_covg_info;

} SpeciesInfo;

SpeciesInfo* get_species_info(dBGraph *db_graph,int max_branch_len, 
                            StrBuf* install_dir,int expected_covg,
                            int ignore_first,int ignore_last);


  void print_json_phylo_group(SpeciesInfo* species_info);

  void print_json_species(SpeciesInfo* species_info);
