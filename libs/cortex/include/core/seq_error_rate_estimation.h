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
  seq_error_rate_estimation.h
*/

#ifndef SEQ_ERR_RATE_H_
#define SEQ_ERR_RATE_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "element.h"
#include "file_reader.h"
#include "graph_info.h"
#include "db_variants.h"
#include "maths.h"



void estimate_seq_error_rate_from_snps_for_each_colour(char* colourlist_snp_alleles, 
                                                       GraphInfo* db_graph_info, 
                                                       dBGraph* db_graph, 
                                                       int ref_colour, 
                                                       //long long genome_size,
                                                       long double default_seq_err_rate, 
                                                       char* output_file);



//pass in a fasta of SNP alleles that we know from SNP-chip genotyping, this sample (colour) does not contain.
//format 
//   >allele with expected zero covg
//   GAAGAGAGAG
//   >allele which we think is homozygous in this sample
//   GAAGATAGAG
//  should be k-1 bases before and after the SNP base itself.
long double estimate_seq_error_rate_for_one_colour_from_snp_allele_fasta(char* fasta, dBGraph* dbg, int colour, 
									 Sequence* seq, KmerSlidingWindow* kmer_window,
									 int (*file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry),
									 dBNode** array_nodes, Orientation* array_or, //int* num_snps_tested,
									 int max_read_length, long double default_seq_err, FILE* fout);

#endif /* SEQ_ERR_RATE_H_ */
