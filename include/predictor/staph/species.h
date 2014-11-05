/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 * **********************************************************************
 *
 * This file is part of myKrobe.
 *
 * myKrobe is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * myKrobe is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with myKrobe.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */
/*
  species.h
*/
#include "dB_graph.h"

typedef enum
  {
    PureStaphAureus =0,
    MixedStaph=1,
    NonStaphylococcal = 2,
  } SampleType;


typedef enum 
{
Saureus=0,
Sepidermidis = 1,
Shaemolyticus = 2,
Sother=3
} Staph_species ;

#define NUM_SPECIES 4

// void add_file_path(int index, StrBuf* install_dir, StrBuf* file_paths,  char file_path );

void load_all_species_panel_file_paths(StrBuf** panel_file_paths ,
                                         StrBuf* install_dir );
void get_coverage_on_panels(int* percentage_coverage,int* median_coverage,
                            int*  total_kmers, StrBuf** panel_file_paths,
                            int max_branch_len, dBGraph *db_graph,
                            int ignore_first,int ignore_last);

void get_coverage_on_best_catalayse_gene(dBGraph *db_graph,int max_branch_len,
                                        StrBuf* install_dir,int ignore_first, 
                                        int ignore_last,
                                        int* percentage_cov_cat,
                                        Covg* median_cov_cat);
boolean catalayse_exists_in_sample(dBGraph *db_graph,int max_branch_len,
                                StrBuf* install_dir,int ignore_first, 
                                int ignore_last);
boolean sample_is_staph(dBGraph *db_graph,int max_branch_len,
                    StrBuf* install_dir,int ignore_first, 
                    int ignore_last);		
SampleType get_species_type(dBGraph *db_graph,int max_branch_len, 
                            StrBuf* install_dir,int expected_covg,
                            int ignore_first,int ignore_last);
// boolean sample_is_mixed();







