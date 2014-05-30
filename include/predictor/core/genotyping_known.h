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
  genotyping_known.h - taking known resistance mutations and genotyping
*/



#ifndef GENOTYPING_KNOWN_H_
#define GENOTYPING_KNOWN_H_

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "global.h"
#include "element.h"
#include "seq.h"
#include "file_reader.h"
#include "db_variants.h"
#include "known_mutations.h"

#define MIN_PERCENT_MUT_ALLELE_PRESENT 80

typedef struct
{
  Covg median_covg;
  Covg min_covg;
  int  percent_nonzero;
} AlleleInfo;

AlleleInfo* alloc_allele_info();
void free_allele_info(AlleleInfo* ai);

int get_next_single_allele_info(FILE* fp, dBGraph* db_graph, AlleleInfo* ainfo,
				Sequence* seq, KmerSlidingWindow* kmer_window,
				int (*file_reader)(FILE * fp, 
						   Sequence * seq, 
						   int max_read_length, 
						   boolean new_entry, 
						   boolean * full_entry),
				dBNode** array_nodes, Orientation*  array_or,
				CovgArray* working_ca, int max_read_length,
				int ignore_first, int ignore_last);



typedef struct
{
  AlleleInfo susceptible_allele;
  AlleleInfo resistant_alleles[6];
  int num_resistant_alleles;
  boolean some_resistant_allele_present;
  int working_current_max_res_allele_present;
  int working_current_max_sus_allele_present;
  KnownMutation var_id;
  GeneMutationGene gene;
}ResVarInfo;

ResVarInfo* alloc_and_init_res_var_info();
void free_res_var_info(ResVarInfo* rvi);
void reset_res_var_info(ResVarInfo* rvi);
void copy_res_var_info(ResVarInfo* from_rvi, ResVarInfo* to_rvi);

void get_next_mutation_allele_info(FILE* fp, dBGraph* db_graph, ResVarInfo* rinfo,
				   Sequence* seq, KmerSlidingWindow* kmer_window,
				   int (*file_reader)(FILE * fp, 
						      Sequence * seq, 
						      int max_read_length, 
						      boolean new_entry, 
						      boolean * full_entry),
				   dBNode** array_nodes, Orientation*  array_or,
				   CovgArray* working_ca, int max_read_length,
				   StrBuf* temp_readid_buf, 
				   StrBuf* temp_mut_buf,
				   StrBuf* temp_gene_name_buf,
				   int ignore_first, int ignore_last, 
				   int expected_covg, KnownMutation* prev_mut);


boolean both_alleles_null(ResVarInfo* rvi);
Covg get_max_covg_on_any_resistant_allele(ResVarInfo* rvi);
int get_max_perc_covg_on_any_resistant_allele(ResVarInfo* rvi);
#endif
