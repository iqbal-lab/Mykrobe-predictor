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
#include "string_buffer.h"

GenePresenceGene map_string_to_gene_presence_gene(StrBuf* sbuf)
{
  if (strcmp(sbuf->buff, "aacA-aphD")==0)
    {
      return aacAaphD;
    }
  else if (strcmp(sbuf->buff, "blaZ")==0)
    {
      return blaZ ;
    }
  else if (strcmp(sbuf->buff, "dfrA")==0)
    {
      return dfrA;
    }
  else if (strcmp(sbuf->buff, "dfrG")==0)
    {
      return dfrG ;
    }
  else if (strcmp(sbuf->buff, "ermA")==0)
    {
      return ermA;
    }
  else if (strcmp(sbuf->buff, "ermB")==0)
    {
      return ermB;
    }
  else if (strcmp(sbuf->buff, "ermC")==0)
    {
      return ermC;
    }
  else if (strcmp(sbuf->buff, "ermT")==0)
    {
      return ermT;
    }
  else if (strcmp(sbuf->buff, "fusB")==0)
    {
      return fusB;
    }
  else if (strcmp(sbuf->buff, "fusC")==0)
    {
      return fusC;
    }
  else if (strcmp(sbuf->buff, "vga(A)LC")==0)
    {
      return vga_A_LC;
    }
  else if (strcmp(sbuf->buff, "msrA")==0)
    {
      return msrA;
    }
  else if (strcmp(sbuf->buff, "mecA")==0)
    {
      return mecA;
    }
  else if (strcmp(sbuf->buff, "tetK")==0)
    {
      return tetK;
    }
  else if (strcmp(sbuf->buff, "tetL")==0)
    {
      return tetL;
    }
  else if (strcmp(sbuf->buff, "tetM")==0)
    {
      return tetM;
    }
  else if (strcmp(sbuf->buff, "vanA")==0)
    {
      return vanA;
    }
  else if (strcmp(sbuf->buff, "mupA")==0)
    {
      return mupA;
    }
  else if (strcmp(sbuf->buff, "mupB")==0)
    {
      return mupB;
    }
  else if (strcmp(sbuf->buff, "luk")==0)
    {
      return luk;
    }
  else 
    {
      die("Unknown gene %s\n", sbuf->buff);
    }

}


boolean map_gene_to_fasta(GenePresenceGene gene, StrBuf* fa, StrBuf* install_dir)
{
  strbuf_append_str(fa, install_dir->buff);
  strbuf_append_str(fa, "data/staph/antibiotics/");

  if (gene==aacAaphD)
    {
      strbuf_append_str(fa, "aacAaphD.fa");
    }
  else if (gene==blaZ)
    {
      strbuf_append_str(fa, "blaZ.fa");
    }
  else if (gene==dfrA)
    {
      strbuf_append_str(fa, "dfrA.fa");
    }
  else if (gene==dfrG)
    {
      strbuf_append_str(fa, "dfrG.fa");
    }
  else if (gene==ermA)
    {
      strbuf_append_str(fa, "ermA.fa");
    }
  else if (gene==ermB)
    {
      strbuf_append_str(fa, "ermB.fa");
    }
  else if (gene==ermC)
    {
      strbuf_append_str(fa, "ermC.fa");
    }
  else if (gene==ermT)
    {
      strbuf_append_str(fa, "ermT.fa");
    }
  else if (gene==fusB)
    {
      strbuf_append_str(fa, "fusB.fa");
    }
  else if (gene==fusC)
    {
      strbuf_append_str(fa, "fusC.fa");
    }
  else if (gene==vga_A_LC)
    {
      strbuf_append_str(fa, "vga_A_LC.fa");
    }
  else if (gene==msrA)
    {
      strbuf_append_str(fa, "msrA.fa");
    }
  else if (gene==mecA)
    {
      strbuf_append_str(fa, "mecA.fa");
    }
  else if (gene==tetK)
    {
      strbuf_append_str(fa, "tetK.fa");
    }
  else if (gene==tetL)
    {
      strbuf_append_str(fa, "tetL.fa");
    }
  else if (gene==tetM)
    {
      strbuf_append_str(fa, "tetM.fa");
    }
  else if (gene==vanA)
    {
      strbuf_append_str(fa, "vanA.fa");
    }
  else if (gene==mupA)
    {
      strbuf_append_str(fa, "mupA.fa");
    }
  else if (gene==mupB)
    {
      strbuf_append_str(fa, "mupB.fa");
    }
  else if (gene==luk)
    {
      strbuf_append_str(fa, "luk.fa");
    }
  else if (gene==unspecified_gpg)
    {
      strbuf_reset(fa);
      return false;
    }
  return true;
}

GeneInfo* alloc_and_init_gene_info()
{
  GeneInfo* gi = (GeneInfo*) calloc(1, sizeof(GeneInfo));
  if (gi==NULL)
    {
      die("Can't alloc gi\n");
    }
  gi->strbuf=strbuf_new(); 
  gi->name = unspecified_gpg;
  return gi;
}
void free_gene_info(GeneInfo* gi)
{
  strbuf_free(gi->strbuf);
  free(gi);
}
void copy_gene_info(GeneInfo* from_gi, GeneInfo* to_gi)
{
  to_gi->median_covg = from_gi->median_covg;
  to_gi->min_covg = from_gi->min_covg;
  to_gi->percent_nonzero = from_gi->percent_nonzero;
  to_gi->median_covg_on_nonzero_nodes = from_gi->median_covg_on_nonzero_nodes;
  to_gi->longest_gap=from_gi->longest_gap;
  to_gi->len=from_gi->len;
  strbuf_reset(to_gi->strbuf);
  strbuf_append_str(to_gi->strbuf, from_gi->strbuf->buff);
  to_gi->name = from_gi->name;
}

void reset_gene_info(GeneInfo* gi)
{
  if (gi==NULL)
    {
      return;
    }
  gi->median_covg=0;
  gi->min_covg=0;
  gi->percent_nonzero=0;
  gi->median_covg_on_nonzero_nodes=0;
  gi->len=0;
  gi->longest_gap=0;
  if (gi->strbuf!=NULL)
    {
      strbuf_reset(gi->strbuf);
    }
  gi->name = unspecified_gpg;
}

//assume you are going to repeatedly call this on a fasta
//file of genes
int get_next_gene_info(FILE* fp, 
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
		       CovgArray* working_ca, 
		       int max_read_length,
		       int expected_covg)

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

  if (num_kmers==0)
    {
      return 0;
    }
  else
    {
      gene_info->len=num_kmers;
    }
  //collect min, median covg on gene and also percentage of kmers with any covg
  boolean too_short=false;
  int ignore_first=0;
  int ignore_last=0;
  gene_info->median_covg = 
    median_covg_on_allele_in_specific_colour(array_nodes,
					     num_kmers,
					     working_ca,
					     0,
					     &too_short,
					     ignore_first, ignore_last, false);

  gene_info->min_covg = 
    min_covg_on_allele_in_specific_colour(array_nodes,
					  num_kmers,
					  0,
					  &too_short,
					  ignore_first, ignore_last);
  gene_info->percent_nonzero = 
    percent_nonzero_on_allele_in_specific_colour(array_nodes,
						 num_kmers,
						 0,
						 &too_short,
						 ignore_first, ignore_last);


  gene_info->median_covg_on_nonzero_nodes = 
    median_covg_ignoring_zeroes_on_allele_in_specific_colour(array_nodes,
							     num_kmers,
							     working_ca,
							     0,
							     &too_short,
							     ignore_first, ignore_last);

  /*  if (gene_info->percent_nonzero==0)
    {
      gene_info->longest_gap = gene_info->len;
    }
  else
    {
      gene_info->longest_gap = 
	longest_gap_on_allele_in_specific_colour(array_nodes, 
						 num_kmers,
						 0,
						 &too_short,
						 ignore_first, ignore_last);
    }

  //in the case when we have partial coverage of a gene, we will want to distinguish
  //minor infections (covg across whole gene, but patchy, and low median given by frequency)
  //from susceptible (very little covg on gene, localised, and low median given by error rate)
  //in both cases, repeat sequence are unhelpful, so get rid of them
  if ( (gene_info->percent_nonzero>0) && (gene_info->percent_nonzero<30) )
    {
      if (gene_info->median_covg_on_nonzero_nodes > 0.5*expected_covg)
	{

	  //dbNode objects not modified
	  int len_covg_array=0;
	  int removed_kmers = 
	    scan_allele_and_remove_repeats_in_covgarray(array_nodes, num_kmers, 
							working_ca, 0, 
							ignore_first, ignore_last,
							expected_covg, &len_covg_array);

	  //this only looks at the covg array, because of the true in final argument.
	  //after this func call the array is sorted.
	  gene_info->median_covg = 
	    median_covg_on_allele_in_specific_colour(array_nodes,
						     num_kmers,
						     working_ca,
						     0,
						     &too_short,
						     ignore_first, ignore_last, true);

	  
	  //leave min covg alone

	  //slightly tortuous:
	  gene_info->percent_nonzero = 
	    (int) (gene_info->percent_nonzero*num_kmers/100) - removed_kmers;
	  if (gene_info->percent_nonzero <0)
	    {
	      gene_info->percent_nonzero =0;
	    }
	  gene_info->percent_nonzero =(int) (100*gene_info->percent_nonzero /num_kmers);

	  if (gene_info->percent_nonzero==0)
	    {
	      gene_info->longest_gap = gene_info->len;
	    }


	  //can't use median on dBNodes as I dont want to modify them
	  gene_info->median_covg_on_nonzero_nodes = 
	    median_of_nonzero_values_of_sorted_covg_array(len_covg_array,
							  working_ca,
							  0,
							  ignore_first, ignore_last);
	  
	}
    }

  */


/*
  if ( (gene_info->percent_nonzero<=20)
       &&
       (gene_info->median_covg_on_nonzero_nodes> 0.5*expected_covg) )
    {
      //this is just repeats. Ignore.
      gene_info->percent_nonzero=0;
      gene_info->median_covg_on_nonzero_nodes=0;
      gene_info->median_covg=0;
      gene_info->min_covg=0;
      gene_info->longest_gap=0;
    }
*/

  strbuf_reset(gene_info->strbuf);
  strbuf_append_str(gene_info->strbuf,
		    seq->name);
  gene_info->name 
    = map_string_to_gene_presence_gene(gene_info->strbuf);

  return num_kmers;
}
