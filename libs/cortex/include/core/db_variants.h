/*
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@tgac.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * The MIT License (MIT)
 * Copyright (c) 2009-2014 <Z. Iqbal and M. Caccamo>
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
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

Covg median_covg_on_allele_in_specific_colour(dBNode** allele, int len, CovgArray* working_ca,
					      int colour, boolean* too_short,
					      int ignore_first, int ignore_last);

Covg median_covg_ignoring_zeroes_on_allele_in_specific_colour(dBNode** allele, int len, CovgArray* working_ca,
							      int colour, boolean* too_short, 
							      int ignore_first, int ignore_last);

Covg min_covg_on_allele_in_specific_colour(dBNode** allele, int len, int colour, boolean* too_short,
					   int ignore_first, int ignore_last);

int percent_nonzero_on_allele_in_specific_colour(dBNode** allele, int len, int colour, boolean* too_short,
						 int ignore_first, int ignore_last);

int num_gaps_on_allele_in_specific_colour(dBNode** allele, 
					  int len, 
					  int colour, 
					  boolean* too_short, 
					  int ignore_first, int ignore_last);

#endif /* DB_VARIANTS_H_ */
