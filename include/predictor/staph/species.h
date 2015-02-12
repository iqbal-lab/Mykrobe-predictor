/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 *  species.h
*/
#include "dB_graph.h"
#include "base_species.h"



void load_all_catalayse_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir );

void load_all_species_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir );

void load_catalayse_threshold(int* thresholds);
void load_all_species_thresholds(int* thresholds);

boolean* create_staph_mask();
boolean* create_non_aureus_mask();

boolean staphylococcus_is_present(SpeciesInfo* species_info);

Species get_best_staph_species(SpeciesInfo* species_info );

void update_complex_presence_and_coverage_from_species(SpeciesInfo* species_info);

SpeciesInfo* get_species_info(dBGraph *db_graph,int max_branch_len, 
                            StrBuf* install_dir,int expected_covg,
                            int ignore_first,int ignore_last);


boolean is_aureus_present(SpeciesInfo* species_info);
boolean is_non_aureus_staph_present(SpeciesInfo* species_info);
void print_json_aureus(SpeciesInfo* species_info,boolean last);

Species get_best_non_aureus_species(SpeciesInfo* species_info );

boolean non_aureus_panels_are_present(SpeciesInfo* species_info);

boolean no_non_aureus_panels_are_present(SpeciesInfo* species_info);



void print_json_best_hit_non_aureus(SpeciesInfo* species_info);


void print_json_aureus_and_best_hit_non_aureus(SpeciesInfo* species_info);

void print_json_species(SpeciesInfo* species_info);

void print_json_lineage(SpeciesInfo* species_info);