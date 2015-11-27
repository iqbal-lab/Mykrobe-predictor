/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 *
 *  gene_presence.c
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
  
    die("Unknown gene %s\n", sbuf->buff);
  

}


boolean map_gene_to_fasta(GenePresenceGene gene, StrBuf* fa, StrBuf* install_dir)
{
  
    return false;  
    
}

const char* map_enum_to_gene_name(GenePresenceGene gene)
{
   switch (gene) 
   {
    
    case unspecified_gpg  : return "unknown";
   }
   return "unknown";
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
  to_gi->num_gaps=from_gi->num_gaps;
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
  gi->num_gaps=0;
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
		       int max_read_length)

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
					     ignore_first, ignore_last);
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
  gene_info->num_gaps = 
    num_gaps_on_allele_in_specific_colour(array_nodes, 
					  num_kmers,
					  0,
					  &too_short,
					  ignore_first, ignore_last);

  gene_info->fasta_id = seq->name;
  char gene_family[LINE_MAX];
  int i = 0;
  while(i<LINE_MAX && ! (seq->name[i] == '\n' || seq->name[i] == ' ' || seq->name[i] == '\t' || seq->name[i] == '\r')){
    gene_family[i] = seq->name[i];
    i++;
  }
  // printf("%s\n", gene_family);
  strbuf_reset(gene_info->strbuf);
  strbuf_append_str(gene_info->strbuf,
		    gene_family);
  memset( gene_family, 0, sizeof(gene_family) );
  gene_info->name 
    = map_string_to_gene_presence_gene(gene_info->strbuf);

  return num_kmers;
}