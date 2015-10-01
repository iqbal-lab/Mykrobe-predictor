/*
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 *  species.h
*/
#include "dB_graph.h"
#include "base_species.h"



void print_json_complex_start();
void print_json_complex_end();
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
void print_json_lineage_end();
char* get_ith_lineage_name(CovgInfo* covg_info, int i);

	typedef enum 
	 {
	 	
	 	Lineage1 = 0,
	 	
	 	Lineage2 = 1,
	 	
	 	Lineage3 = 2,
	 	
	 	Lineage4 = 3,
	 	
	 	Lineage5 = 4,
	 	
	 	Lineage6 = 5,
	 	
    unknownlineage=6
	   	} Lineage ;
	#define NUM_Lineage 6
   	
  void map_lineage_enum_to_str(Lineage sp, StrBuf* sbuf);
  void load_all_lineage_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir );
  void lineage_threshold(int* thresholds);


void print_json_species_start();
void print_json_species_end();
char* get_ith_species_name(CovgInfo* covg_info, int i);

	typedef enum 
	 {
	 	
	 	Abscessus = 0,
	 	
	 	Africanum = 1,
	 	
	 	Animallineage = 2,
	 	
	 	Aromaticivorans = 3,
	 	
	 	Avium = 4,
	 	
	 	Beijingsublineage = 5,
	 	
	 	Bovis = 6,
	 	
	 	Branderi = 7,
	 	
	 	Caprae = 8,
	 	
	 	Chelonae = 9,
	 	
	 	Chlorophenolicum = 10,
	 	
	 	Chubuense = 11,
	 	
	 	Colombiense = 12,
	 	
	 	Crocinum = 13,
	 	
	 	Flavescens = 14,
	 	
	 	Fluoranthenivorans = 15,
	 	
	 	Fortuitum = 16,
	 	
	 	Gilvum = 17,
	 	
	 	Gordonae = 18,
	 	
	 	Hodleri = 19,
	 	
	 	Interjectum = 20,
	 	
	 	Intracellulare = 21,
	 	
	 	Kansasii = 22,
	 	
	 	Lentiflavum = 23,
	 	
	 	Leprae = 24,
	 	
	 	Malmoense = 25,
	 	
	 	Marinum = 26,
	 	
	 	Mucogenicum = 27,
	 	
	 	Pallens = 28,
	 	
	 	Peregrinum = 29,
	 	
	 	Phage = 30,
	 	
	 	Pyrenivorans = 31,
	 	
	 	Rhodesiae = 32,
	 	
	 	Rufum = 33,
	 	
	 	Rutilum = 34,
	 	
	 	Scrofulaceum = 35,
	 	
	 	Senegalense = 36,
	 	
	 	Smegmatis = 37,
	 	
	 	Sphagni = 38,
	 	
	 	Szulgai = 39,
	 	
	 	Triplex = 40,
	 	
	 	Tuberculosis = 41,
	 	
	 	Tusciae = 42,
	 	
	 	Ulcerans = 43,
	 	
	 	Vaccae = 44,
	 	
	 	Xenopi = 45,
	 	
    unknownspecies=46
	   	} Species ;
	#define NUM_Species 46
   	
  void map_species_enum_to_str(Species sp, StrBuf* sbuf);
  void load_all_species_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir );
  void species_threshold(int* thresholds);


typedef struct
{

  CovgInfo* complex_covg_info;

  CovgInfo* lineage_covg_info;

  CovgInfo* species_covg_info;

} SpeciesInfo;

SpeciesInfo* get_species_info(dBGraph *db_graph,int max_branch_len, 
                            StrBuf* install_dir,int expected_covg,
                            int ignore_first,int ignore_last);


  void print_json_complex(SpeciesInfo* species_info);

  void print_json_lineage(SpeciesInfo* species_info);

  void print_json_species(SpeciesInfo* species_info);
