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
  genotyping_known.c - taking known resistance mutations and genotyping
*/

#include <string.h>
#include "genotyping_known.h"
#include "string_buffer.h"
#include "known_mutations.h"

AlleleInfo* alloc_allele_info()
{
  AlleleInfo* ai = (AlleleInfo*) calloc(1, sizeof(AlleleInfo));
  if (ai==NULL)
    {
      die("Disaster - cant evben alloc a tiny alleleinfo");
    }
  return ai;
}

void free_allele_info(AlleleInfo* ai)
{
  free(ai);
}

ResVarInfo* alloc_and_init_res_var_info()
{
  ResVarInfo* rvi = calloc(1, sizeof(ResVarInfo));
  if (rvi==NULL)
    {
      die("Disaster - cant evben alloc a tiny resvarinfo object");
    }
  rvi->num_resistant_alleles = 0;
  rvi->var_id = NotSpecified;
  rvi->gene = Unknown;
  rvi->some_resistant_allele_present = false;
  return rvi;
}

void free_res_var_info(ResVarInfo* rvi)
{
  free(rvi);
}

void copy_res_var_info(ResVarInfo* from_rvi, ResVarInfo* to_rvi)
{
  memcpy(to_rvi, from_rvi, sizeof(ResVarInfo));
}

void reset_res_var_info(ResVarInfo* rvi)
{
  memset(rvi,0, sizeof(ResVarInfo));
  rvi->some_resistant_allele_present=false;
  /*  rvi->susceptible_allele.median_covg=0;
  rvi->susceptible_allele.min_covg=0;
  rvi->susceptible_allele.percent_nonzero=0;
  int i;
  for (i=0; i<rvi->num_resistant_alleles; i++)
    {
      rvi->resistant_alleles[i].median_covg=0;
      rvi->resistant_alleles[i].min_covg=0;
      rvi->resistant_alleles[i].percent_nonzero=0;
      } */
}



//util funcs

//finds the (single-digit_ number after the first minus sign
int find_number_resistant_alleles(StrBuf* sbuf)
{
  uint32_t i;
  for (i=0; i<sbuf->len; i++)
    {
      char c = sbuf->buff[i];
      if (c=='-')
	{
	  char d = sbuf->buff[i+1];
	  int id = d - '0';
	  return id;
	}
    }
  die("Cannot parse %s to find the number of resistant alleles in this read-id\n", sbuf->buff);
}


void find_gene_name(StrBuf* sbuf_in, StrBuf* sbuf_out)
{
  uint32_t i;
  int count=0;
  for (i=0; i<sbuf_in->len; i++)
    {
      char c = sbuf_in->buff[i];
      if (c=='-')
	{
	  count++;
	  if (count==2)
	    {
	      strbuf_reset(sbuf_out);
	      strbuf_append_str(sbuf_out, sbuf_in->buff+i+1);
	      i = sbuf_in->len;
	    }
	}
    }
  if (count==2)
    {
      return;
    }
  die("Cannot parse %s to find the gene name\n", sbuf_in->buff);
}

void find_mutation_name(StrBuf* sbuf_in, StrBuf* sbuf_out)
{
  uint32_t i;
  int count=0;
  int first=0;
  int second = 0;
  strbuf_reset(sbuf_out);
  strbuf_append_str(sbuf_out, sbuf_in->buff);

  for (i=0; i<sbuf_in->len; i++)
    {
      char c = sbuf_in->buff[i];
      if (c=='_')
	{
	  count++;
	  if (count==1)
	    {
	      first = i;
	    }
	  else if (count==2)
	    {
	      second=i;
	    }
	}
    }
  if ( (first==0) || (second==0) )
    {
      die("Cannot parse %s to find the gene name\n", sbuf_in->buff);
    }
  strbuf_delete(sbuf_out, second, sbuf_out->len-second);
  strbuf_delete(sbuf_out, 0, first-1);
  printf("Got mut name %s\n", sbuf_out->buff);
}


//assume you are going to repeatedly call this on a fasta
//file of known mutations, with the first allele being the susceptible allele, in this format
//>ref_F99I_panel_fasta_sub--3
// the "ref" says it's the susceptible allele. 
//then it has the amino acid substitution (in future will end up with more general, but right now ALL the cases from
//the literature are AA substitutions or indels in genes
//the "sub" says it is a substitution
//then after the -- it says how many resistant alleles there are
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
				   StrBuf* temp_gene_name_buf)
{


  boolean full_entry=true;
  int dummy_colour_ignored = -1;
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

  //read- is in this format

  //>ref_F99I_panel_fasta_sub--3--dfrB        so this is >ref means susceptible allele, then identifier of which reference, then sub=susbstitution
  //                                          then 3 means there are 3 possibe resistance alleles all implying F-->Y, and finally the gene name
  //CATGTTTTTATATTTGGAGGGCAAACATTATTTGAAGAAATGATTGATAAAGTGGACGAC
  //>ATA_F99I_panel_fasta_alt--1--dfrB
  //CATGTTTTTATATTTGGAGGGCAAACATTAATAGAAGAAATGATTGATAAAGTGGACGAC
  //ATG
  //>ATT_F99I_panel_fasta_alt--2--dfrB
  //CATGTTTTTATATTTGGAGGGCAAACATTAATTGAAGAAATGATTGATAAAGTGGACGAC
  //ATG
  //>ATC_F99I_panel_fasta_alt--3--dfrB
  //CATGTTTTTATATTTGGAGGGCAAACATTAATCGAAGAAATGATTGATAAAGTGGACGAC
  //ATG

  strbuf_reset(temp_readid_buf);
  strbuf_reset(temp_mut_buf);
  strbuf_append_str(temp_readid_buf, seq->name);
  rinfo->num_resistant_alleles = find_number_resistant_alleles(temp_readid_buf);
  find_gene_name(temp_readid_buf, temp_gene_name_buf);
  rinfo->gene = map_gene_name_str_to_genename(temp_gene_name_buf);
  find_mutation_name(temp_readid_buf, temp_mut_buf);
  rinfo->var_id = map_mutation_name_to_enum(temp_mut_buf ,rinfo->gene);

  //collect min, median covg on allele and also percentage of kmers with any covg
  boolean too_short=false;
  rinfo->susceptible_allele.median_covg = 
    median_covg_on_allele_in_specific_colour(array_nodes, 
					     num_kmers, 
					     working_ca, 
					     0, 
					     &too_short);
  rinfo->susceptible_allele.min_covg = 
    min_covg_on_allele_in_specific_colour(array_nodes, 
					  num_kmers, 
					  0, 
					  &too_short);
  rinfo->susceptible_allele.percent_nonzero = 
    percent_nonzero_on_allele_in_specific_colour(array_nodes, 
						 num_kmers, 
						 0, 
						 &too_short);

  int i;
  for (i=0; i<rinfo->num_resistant_alleles; i++)
    {
      num_kmers = align_next_read_to_graph_and_return_node_array(fp, 
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
      too_short=false;
      rinfo->resistant_alleles[i].median_covg = 
	median_covg_on_allele_in_specific_colour(array_nodes, 
						 num_kmers, 
						 working_ca,
						 0,
						 &too_short);
      rinfo->resistant_alleles[i].min_covg = 
	min_covg_on_allele_in_specific_colour(array_nodes,
					      num_kmers,
					      0,
					      &too_short);
      rinfo->resistant_alleles[i].percent_nonzero = 
	percent_nonzero_on_allele_in_specific_colour(array_nodes,
						     num_kmers,
						     0,
						     &too_short);
      if (rinfo->resistant_alleles[i].percent_nonzero>=MIN_PERCENT_MUT_ALLELE_PRESENT)
	{
	  rinfo->some_resistant_allele_present=true;
	}
    }
}
