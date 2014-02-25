/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */
/*
  graph_info.h
*/

#ifndef GRAPH_INFO_H_
#define GRAPH_INFO_H_

#include <stdio.h>
#include <stdint.h>

#include "global.h"
#include "element.h"

#define MAX_LEN_NAME_GRAPH_CLEANED_AGAINST 10000;
extern int MAX_LEN_SAMPLE_NAME;

// This is alloced and initialised by the GraphInfo object. I miss C++
typedef struct
{
  boolean tip_clipping;
  boolean remv_low_cov_sups;
  boolean remv_low_cov_nodes;
  int len_name_of_graph_against_which_was_cleaned;

  // These two are ints as they can be set to -1
  // DEV: convert to Covg, use remv_low_cov_sups & remv_low_cov_nodes
  //      instead of remv_low_cov_sups_thresh == -1 and
  //      remv_low_cov_nodes_thresh == -1
  int remv_low_cov_sups_thresh;
  int remv_low_cov_nodes_thresh;

  // eg cleaning a low covg sample against cleaned pool of population
  boolean cleaned_against_another_graph;
  char* name_of_graph_against_which_was_cleaned;
} ErrorCleaning;

ErrorCleaning* error_cleaning_alloc_and_init();
void error_cleaning_initialise_except_pool_cleaning(ErrorCleaning* cl);
void error_cleaning_initialise(ErrorCleaning* cl);
void error_cleaning_free(ErrorCleaning* );
void error_cleaning_assign_with_OR(ErrorCleaning* target, ErrorCleaning* src, boolean dont_set_pool_cleaning);

typedef struct{
  int           sample_id_lens[NUMBER_OF_COLOURS];
  char*          sample_ids[NUMBER_OF_COLOURS];
  uint64_t       total_sequence[NUMBER_OF_COLOURS];
  uint32_t       mean_read_length[NUMBER_OF_COLOURS];
  long double    seq_err[NUMBER_OF_COLOURS];
  ErrorCleaning* cleaning[NUMBER_OF_COLOURS];
} GraphInfo;

GraphInfo* graph_info_alloc_and_init();
void graph_info_free(GraphInfo* ginfo);
void graph_info_initialise(GraphInfo* ginfo);
void graph_info_initialise_one_colour_except_pool_cleaning(GraphInfo* ginfo, int colour);

//set total amount of sequence in a colour
void graph_info_set_seq(GraphInfo* ginfo, int colour, long long num_bp);
//increment total amoutn of sequence in a colour (eg when you merge in a new binary)
long long graph_info_increment_seq(GraphInfo* ginfo, int colour, long long num_bp);
//set mean read length in a colour
void graph_info_set_mean_readlen(GraphInfo* ginfo, int colour, int len);
void graph_info_set_tip_clipping(GraphInfo* ginfo, int colour);
void graph_info_unset_tip_clipping(GraphInfo* ginfo, int colour);
void graph_info_set_remv_low_cov_sups(GraphInfo* ginfo, int colour, int thresh);
void graph_info_unset_remv_low_cov_sups(GraphInfo* ginfo, int colour);
void graph_info_set_remv_low_cov_nodes(GraphInfo* ginfo, int colour, int thresh);
void graph_info_UNset_remv_low_cov_nodes(GraphInfo* ginfo, int colour);
void graph_info_set_seq_err(GraphInfo* ginfo, int col, long double err);
void graph_info_set_all_metadata(GraphInfo* target, GraphInfo* src, int colour, boolean dont_set_pool_cleaning);
void graph_info_set_specific_colour_to_cleaned_against_pool(GraphInfo* ginfo, int colour, char* multicol_binary, int colour_in_multicol);
void graph_info_unset_specific_colour_from_cleaned_against_pool(GraphInfo* ginfo, int colour);

//update mean read length in a colour, eg when you merge a new binary
int graph_info_update_mean_readlen(GraphInfo* ginfo, int colour, int previous_mean, long long previous_seq, int mean_readlen_in_added_data, long long added_seq);
void graph_info_update_mean_readlen_and_total_seq(GraphInfo* ginfo, int colour,
                                                  unsigned long mean_readlen_in_added_data,
				                                  unsigned long long added_seq);
void graph_info_set_sample_ids(char** sample_ids, int num_ids, GraphInfo* ginfo, int first_col);

double get_total_depth_of_coverage_across_colours(GraphInfo* ginfo, long long genome_length);
int get_mean_readlen_across_colours(GraphInfo* ginfo);
void read_estimated_seq_errors_from_file(GraphInfo* ginfo, FILE* fp);
void print_seq_err_rates_to_screen(GraphInfo* ginfo);

#endif /* GRAPH_INFO_H_ */
