/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 *  species.h
*/

/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 *  species.h
*/
#include "dB_graph.h"
#include "base_species.h"

void are_mtbc_and_ntm_present(int* percentage_coverage,boolean* present, 
                                  int* num_panels);
boolean sample_is_mixed(boolean NTM_is_present,boolean MTBC_is_present);
boolean sample_is_MTBC(boolean NTM_is_present,boolean MTBC_is_present);
boolean sample_is_NTM(boolean NTM_is_present,boolean MTBC_is_present);
boolean is_NTM_present(SpeciesInfo* species_info);
boolean is_MTBC_present(SpeciesInfo* species_info);

void load_all_mtbc_and_ntm_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir );
void load_all_lineage_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir );
void load_all_species_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir );


SpeciesInfo* get_species_info(dBGraph *db_graph,int max_branch_len, 
                            StrBuf* install_dir,int expected_covg,
                            int ignore_first,int ignore_last);
void find_which_panels_are_present(CovgInfo* covg_info);

void map_complex_enum_to_str(Complex sp, StrBuf* sbuf);
void map_species_enum_to_str(Species sp, StrBuf* sbuf);
void map_lineage_enum_to_str(Lineage sp, StrBuf* sbuf);


char* get_ith_complex_name(CovgInfo* covg_info, int i);
char* get_ith_species_name(CovgInfo* covg_info, int i);
char* get_ith_lineage_name(CovgInfo* covg_info, int i);


char* get_pure_lineage_name(SpeciesInfo* species_info);
void print_json_indiv_phylo(CovgInfo* covg_info,
                           char* (*get_ith_name)(CovgInfo*, int));
boolean* create_mask(boolean default_value);
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
