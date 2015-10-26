/*
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 *  base_species.h
*/

#ifndef BASE_SPECIES_H_
#define BASE_SPECIES_H_

#include "dB_graph.h"


typedef struct
{
  boolean present[10000];
  int percentage_coverage[10000];
  int median_coverage[10000];
  int percentage_threshold[10000];  
  int num_panels_present;
  int NUM_PANELS;
} CovgInfo;


void get_coverage_on_panels(int* percentage_coverage,int* median_coverage,
                            StrBuf** panel_file_paths,
                            int max_branch_len, dBGraph *db_graph,
                            int ignore_first, int ignore_last , int NUM_PANELS);

int get_ith_coverage_panel(CovgInfo* covg_info, int i);
int get_ith_present_panel(CovgInfo* covg_info, int i);


boolean is_percentage_coverage_above_threshold(int per_cov,int threshold);


CovgInfo* alloc_and_init_covg_info();
CovgInfo* get_coverage_info(dBGraph *db_graph,
                          StrBuf** file_paths,
                          int max_branch_len,
                          int NUM_PANELS,
                          int ignore_first,
                          int ignore_last,
                          void (*load_thresholds)()
                          );

int max(int int1,int int2);
int get_best_hit(CovgInfo* covg_info,boolean* mask);
boolean* create_mask(boolean default_value);
boolean panels_are_present(CovgInfo* covg_info ,  boolean* mask);

void print_json_indiv_phylo(CovgInfo* covg_info,
                           char* (*get_ith_name)(CovgInfo*, int));




#endif

