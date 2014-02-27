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
  gene_presence.c
*/


#include "gene_presence.h"
#include "element.h"
#include "file_reader.h"
#include "build.h"
#include "db_differentiation.h"
#include "maths.h"


GeneInfo* alloc_and_init_gene_info()
{
  GeneInfo* gi = (GeneInfo*) calloc(1, sizeof(GeneInfo));
  if (gi==NULL)
    {
      die("Can't alloc gi\n");
    }
  gi->name=strbuf_new(); 
  return gi;
}
void free_gene_info(GeneInfo* gi)
{
  strbuf_free(gi->name);
  free(gi);
}

void reset_gene_info(GeneInfo* gi)
{
  gi->median_covg=0;
  gi->min_covg=0;
  gi->percent_nonzero=0;
  strbuf_reset(gi->name);
}

//assume you are going to repeatedly call this on a fasta
//file of genes
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
			CovgArray* working_ca, int max_read_length)
{

  boolean full_entry=true;
  int dummy_colour_ignored=0;
  int num_kmers = 
    align_next_read_to_graph_and_return_node_array(fp, 
						   max_read_length, 
						   array_nodes, 
						   array_or, 
						   false, 
						   &full_entry, 
						   file_reader,
						   seq, 
						   kmer_window, 
						   db_graph, 
						   dummy_colour_ignored);

  
  //collect min, median covg on gene and also percentage of kmers with any covg
  boolean too_short=false;
  gene_info->median_covg = 
    median_covg_on_allele_in_specific_colour(array_nodes,
					     num_kmers,
					     working_ca,
					     0,
					     &too_short);
  gene_info->min_covg = 
    min_covg_on_allele_in_specific_colour(array_nodes,
					  num_kmers,
					  0,
					  &too_short);
  gene_info->percent_nonzero = 
    percent_nonzero_on_allele_in_specific_colour(array_nodes,
						 num_kmers,
						 0,
						 &too_short);
  strbuf_reset(gene_info->name);
  strbuf_append_str(gene_info->name,
		    seq->seq);
  
}
