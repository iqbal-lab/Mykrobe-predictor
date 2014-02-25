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


typedef struct
{
  Covg median_covg;
  Covg min_covg;
  int  percent_nonzero;
} AlleleInfo;

typedef struct
{
  AlleleInfo* susceptible_allele;
  AlleleInfo* resistant_allele;
  StrBuf* name_of_gene;
}ResVarInfo;

void get_next_mutation_allele_info(FILE* fp, dBGraph* db_graph, ResVarInfo* r1info,
				   Sequence* seq, KmerSlidingWindow* kmer_window,
				   (*file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry),
				   dBNode** array_nodes, Orientation*  array_or);


#endif
