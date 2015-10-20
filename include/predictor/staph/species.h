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
	 	
	 	Staphaureus = 0,
	 	
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


  CovgInfo* other_covg_info;
  
} SpeciesInfo;

void print_json_phylogenetics(SpeciesInfo* species_info);


SpeciesInfo* get_species_info(dBGraph *db_graph,int max_branch_len, 
                            StrBuf* install_dir,int expected_covg,
                            int ignore_first,int ignore_last);


  
    void print_json_phylo_group(SpeciesInfo* species_info);
  

  
    void print_json_species(SpeciesInfo* species_info);
  




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
