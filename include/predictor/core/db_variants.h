/*
 * Copyright 2009-2013 Zamin Iqbal and Mario Caccamo
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
  db_variants.h
*/

#ifndef DB_VARIANTS_H_
#define DB_VARIANTS_H_

// system libraries                                                                                                                                                                      
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "graph_info.h"
#include "element.h"
#include "experiment.h"
#include "dB_graph.h"
#include "experiment.h"
#include "dB_graph_supernode.h"

#define MAX_VARNAME_LEN 200


//utility functions
boolean get_num_effective_reads_on_branch(Covg* array, dBNode** allele, int how_many_nodes, 
					  boolean use_median, CovgArray* working_ca, GraphInfo* ginfo, int kmer);
boolean get_num_effective_reads_on_unique_part_of_branch(Covg* array1, dBNode** allele1, int len1, 
							 Covg* array2, dBNode** allele2, int len2,
							 CovgArray* working_ca, GraphInfo* ginfo, int kmer,
							 boolean use_ref_allele_info, int ref_colour, int ref_allele);


Covg count_reads_on_allele_in_specific_colour(dBNode** allele, int len, int colour, boolean* too_short);
Covg count_reads_on_allele_in_specific_colour_given_array_of_cvgs(Covg* covgs, int len, boolean* too_short);
Covg count_reads_on_allele_in_specific_func_of_colours(
						       dBNode** allele, int len,
						       Covg (*sum_of_covgs_in_desired_colours)(const Element *),
						       boolean* too_short);

int scan_allele_and_remove_repeats_in_covgarray(dBNode** allele, int len, CovgArray* working_ca,
						 int colour, 
						int ignore_first, int ignore_last, int expected_covg,
						int* len_populated_in_covg_array);

Covg median_covg_on_allele_in_specific_colour(dBNode** allele, int len, CovgArray* working_ca,
					      int colour, boolean* too_short,
					      int ignore_first, int ignore_last,
					      boolean working_ca_pre_initialised);

Covg median_covg_ignoring_zeroes_on_allele_in_specific_colour(dBNode** allele, int len, CovgArray* working_ca,
							      int colour, boolean* too_short, 
							      int ignore_first, int ignore_last);

Covg median_of_nonzero_values_of_sorted_covg_array(int len_populated_in_array, CovgArray* working_ca,
						   int colour, 
						   int ignore_first, int ignore_last);

Covg min_covg_on_allele_in_specific_colour(dBNode** allele, int len, int colour, boolean* too_short,
					   int ignore_first, int ignore_last);

int percent_nonzero_on_allele_in_specific_colour(dBNode** allele, int len, int colour, boolean* too_short,
						 int ignore_first, int ignore_last);

int longest_gap_on_allele_in_specific_colour(dBNode** allele, 
					     int len, 
					     int colour, 
					     boolean* too_short, 
					     int ignore_first, int ignore_last);

#endif /* DB_VARIANTS_H_ */
