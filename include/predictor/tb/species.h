/*
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 *  species.h
*/
#include "dB_graph.h"
#include "base_species.h"




void print_json_complex_start();
// void print_json_complex_end();
char* get_ith_complex_name(CovgInfo* covg_info, int i);

	typedef enum 
	 {
	 	
	 	Mtbc = 0,
	 	
	 	Ntm = 1,
	 	
    unknowncomplex=2
	   	} Complex ;
	#define NUM_Complex 2
   	
  void map_complex_enum_to_str(Complex sp, StrBuf* sbuf);
  void load_all_complex_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir );
  void complex_threshold(int* thresholds);


void print_json_lineage_start();
// void print_json_lineage_end();
char* get_ith_lineage_name(CovgInfo* covg_info, int i);

	typedef enum 
	 {
	 	
	 	Lineage1 = 0,
	 	
	 	Lineage2 = 1,
	 	
	 	Lineage3 = 2,
	 	
	 	Lineage4 = 3,
	 	
	 	Lineage5 = 4,
	 	
	 	Lineage6 = 5,
	 	
	 	Animallineage = 6,
	 	
	 	Beijingsublineage = 7,
	 	
    unknownlineage=8
	   	} Lineage ;
	#define NUM_Lineage 8
   	
  void map_lineage_enum_to_str(Lineage sp, StrBuf* sbuf);
  void load_all_lineage_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir );
  void lineage_threshold(int* thresholds);


void print_json_species_start();
// void print_json_species_end();
char* get_ith_species_name(CovgInfo* covg_info, int i);

	typedef enum 
	 {
	 	
	 	Abscessus = 0,
	 	
	 	Africanum = 1,
	 	
	 	Aromaticivorans = 2,
	 	
	 	Avium = 3,
	 	
	 	Bovis = 4,
	 	
	 	Branderi = 5,
	 	
	 	Caprae = 6,
	 	
	 	Chelonae = 7,
	 	
	 	Chlorophenolicum = 8,
	 	
	 	Chubuense = 9,
	 	
	 	Colombiense = 10,
	 	
	 	Crocinum = 11,
	 	
	 	Flavescens = 12,
	 	
	 	Fluoranthenivorans = 13,
	 	
	 	Fortuitum = 14,
	 	
	 	Gilvum = 15,
	 	
	 	Gordonae = 16,
	 	
	 	Hodleri = 17,
	 	
	 	Interjectum = 18,
	 	
	 	Intracellulare = 19,
	 	
	 	Kansasii = 20,
	 	
	 	Lentiflavum = 21,
	 	
	 	Leprae = 22,
	 	
	 	Malmoense = 23,
	 	
	 	Marinum = 24,
	 	
	 	Mucogenicum = 25,
	 	
	 	Pallens = 26,
	 	
	 	Peregrinum = 27,
	 	
	 	Phage = 28,
	 	
	 	Pyrenivorans = 29,
	 	
	 	Rhodesiae = 30,
	 	
	 	Rufum = 31,
	 	
	 	Rutilum = 32,
	 	
	 	Scrofulaceum = 33,
	 	
	 	Senegalense = 34,
	 	
	 	Smegmatis = 35,
	 	
	 	Sphagni = 36,
	 	
	 	Szulgai = 37,
	 	
	 	Triplex = 38,
	 	
	 	Tuberculosis = 39,
	 	
	 	Tusciae = 40,
	 	
	 	Ulcerans = 41,
	 	
	 	Vaccae = 42,
	 	
	 	Xenopi = 43,
	 	
    unknownspecies=44
	   	} Species ;
	#define NUM_Species 44
   	
  void map_species_enum_to_str(Species sp, StrBuf* sbuf);
  void load_all_species_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir );
  void species_threshold(int* thresholds);


typedef struct
{

  CovgInfo* complex_covg_info;

  CovgInfo* lineage_covg_info;

  CovgInfo* species_covg_info;

  
} SpeciesInfo;

void print_json_phylogenetics(SpeciesInfo* species_info);


SpeciesInfo* get_species_info(dBGraph *db_graph,int max_branch_len, 
                            StrBuf* install_dir,int expected_covg,
                            int ignore_first,int ignore_last);


  
    void print_json_complex(SpeciesInfo* species_info);
  

  
    void print_json_lineage(SpeciesInfo* species_info);
  

  
    void print_json_species(SpeciesInfo* species_info);
  


int get_expected_covg(SpeciesInfo* species_info);


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


