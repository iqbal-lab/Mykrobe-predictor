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
void load_individual_species_panel_file_path(int species_index, char file_path, StrBuf* install_dir, StrBuf* species_file_paths);
void load_all_species_panel_file_paths(StrBuf* install_dir,StrBuf* species_file_paths);


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






