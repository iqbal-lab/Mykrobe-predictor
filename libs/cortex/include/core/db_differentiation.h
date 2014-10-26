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
  db_differentiation.h
*/

#ifndef DB_DIFFERENTIATION_H_
#define DB_DIFFERENTIATION_H_

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "global.h"
#include "element.h"
#include "dB_graph.h"
#include "seq.h"
#include "file_reader.h"
#include "db_variants.h"


void align_list_of_fastaq_to_graph_and_print_coverages_in_all_colours(
  FileFormat format, char* list_of_fastaq, int max_read_length,
  int* array_of_colours, char** array_of_names_of_colours,
  int num_of_colours, dBGraph* db_graph,int fastq_ascii_offset,
  boolean is_for_testing, char** for_test_array_of_strings,
  int* for_test_index, boolean mark_nodes_for_dumping);

void print_percent_agreement_for_each_colour_for_each_read(char* fasta, int max_read_length, 
							   dBGraph* db_graph, char** list_sample_ids);
void print_supernodes_touched_by_fasta( char* fasta, int max_read_length, 
					dBGraph* db_graph, char* outfile_stub );
#endif /* DB_DIFFERENTIATION_H_ */
