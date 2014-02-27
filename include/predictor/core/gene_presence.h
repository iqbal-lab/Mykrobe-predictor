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
  gene_presence.h - build the de Bruijn graph
*/

#include "global.h" // Covg def 
#include "seq.h" // Sequence def 
#include "binary_kmer.h" // KmerSlidingWindow
#include "dB_graph.h"// dBGraph
#include "dB_graph_supernode.h"

typedef struct
{
  Covg median_covg;
  Covg min_covg;
  int  percent_nonzero;
  StrBuf* strbuf;
  StrBuf* name;
} GeneInfo;

GeneInfo* alloc_and_init_gene_info();
void free_gene_info(GeneInfo* rvi);
void reset_gene_info(GeneInfo* rvi);

void get_next_gene_info(FILE* fp, 
			dBGraph* db_graph, 
			GeneInfo* gene_info,
			Sequence* seq, 
			KmerSlidingWindow* kmer_window,
			int (*file_reader)(FILE * fp, 
				       Sequence * seq, 
				       int max_read_length, 
				       boolean new_entry, 
				       boolean * full_entry),
			dBNode** array_nodes, 
			Orientation*  array_or,
			CovgArray* working_ca, int max_read_length);
