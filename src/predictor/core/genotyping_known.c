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



//we want to use the same FASTA files irrespective of k.
//so we expect alleles ay have extra bases on the front and en
//ie we have fasta with 30 bases before/after a SNP, but if we are using k=21,
//then we must ignore the first 9 and last 9 bases
int get_next_single_allele_info(FILE* fp, dBGraph* db_graph, AlleleInfo* ainfo,
				boolean get_median_on_nonzero,
				Sequence* seq, KmerSlidingWindow* kmer_window,
				int (*file_reader)(FILE * fp, 
						    Sequence * seq, 
						   int max_read_length, 
						   boolean new_entry, 
						    boolean * full_entry),
				dBNode** array_nodes, Orientation*  array_or,
				CovgArray* working_ca, int max_read_length,
				int ignore_first, int ignore_last)

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

  if (num_kmers>0)
    {
      //collect min, median covg on allele and also percentage of kmers with any covg
      boolean too_short=false;
      ainfo->median_covg = 
	median_covg_on_allele_in_specific_colour(array_nodes, 
						 num_kmers, 
						 working_ca, 
						 0, 
						 &too_short,
						 ignore_first, ignore_last);
      ainfo->min_covg = 
	min_covg_on_allele_in_specific_colour(array_nodes, 
					      num_kmers, 
					      0, 
					      &too_short,
					      ignore_first, ignore_last);
      ainfo->percent_nonzero = 
	percent_nonzero_on_allele_in_specific_colour(array_nodes, 
						     num_kmers, 
						     0, 
						     &too_short,
						     ignore_first, ignore_last);
      if (get_median_on_nonzero==true)
	{
	  ainfo->median_covg_on_nonzero_nodes = 
	    median_covg_ignoring_zeroes_on_allele_in_specific_colour(array_nodes, 
								     num_kmers, 
								     working_ca,
								     0, 
								     &too_short,
								     ignore_first, ignore_last);
	}
    }
  return num_kmers;

}

CalledVariant* alloc_and_init_called_variant_array()
{
  int i;
  CalledVariant* called_variants = malloc(NUM_KNOWN_MUTATIONS * sizeof(*called_variants));
  for (i=0; i<=NUM_KNOWN_MUTATIONS; i++){  	 
	 called_variants[i].var_id = NotSpecified;
	}
  return called_variants;
}
void free_called_variant_array(CalledVariant* cva)
{
  free(cva);
}

void print_called_variants(CalledVariant* called_variants)
{
	// CalledVariant* current_called_variant;
	int i;
	// Iterate through all the variants and print the enum strings
	for (i=0; i<=NUM_KNOWN_MUTATIONS; i++){
		if (called_variants[i].var_id != NotSpecified){
			printf("%s(R:%iS%i)\n", map_enum_to_mutation_name(called_variants[i].var_id),
				called_variants[i].max_res_allele_present,called_variants[i].max_sus_allele_present);
		}
	}
}

CalledGene* alloc_and_init_called_genes_array()
{
  CalledGene* called_genes = malloc(NUM_KNOWN_GENES * sizeof(*called_genes));
  int i;
  for (i=0; i<=NUM_KNOWN_MUTATIONS; i++){
	  called_genes[i].gene = Unknown;
	}
  return called_genes;
}
void free_called_genes_array(CalledGene* cg)
{
  free(cg);
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
  rvi->working_current_max_res_allele_present=0;
  rvi->working_current_max_sus_allele_present=0;
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
  rvi->working_current_max_res_allele_present=0;
  rvi->working_current_max_sus_allele_present=0;
}



//util funcs

//if both alleles have median zero
boolean both_alleles_null(ResVarInfo* rvi)
{
  Covg c = get_max_perc_covg_on_any_resistant_allele(rvi);

  if ( (rvi->susceptible_allele.percent_nonzero==0)
       && (c==0) )
    {
      return true;
    }
  return false;
}
Covg get_max_covg_on_any_resistant_allele(ResVarInfo* rvi)
{
  int i;
  Covg max=0;
  for (i=0; i<rvi->num_resistant_alleles; i++)
    {
      Covg c = rvi->resistant_alleles[i].median_covg;
      if (c>max)
	{
	  max=c;
	}
    }
  return max;
}


int get_max_perc_covg_on_any_resistant_allele(ResVarInfo* rvi)
{
  int i;
  int max=0;
  for (i=0; i<rvi->num_resistant_alleles; i++)
    {
      int c = rvi->resistant_alleles[i].percent_nonzero;
      if (c>max)
	{
	  max=c;
	}
    }
  return max;
}


//finds the (number after the first minus sign
int find_number_resistant_alleles(StrBuf* sbuf)
{
  uint32_t i;
  for (i=0; i<sbuf->len; i++)
    {
      char c = sbuf->buff[i];
      if (c=='-')
	{
	  
	  int id;
	  char d;
	  // 4 digit max 
	    if (sbuf->buff[i+2] == '-')
	      {
	       d = sbuf->buff[i+1];
	       id = d - '0';
	      }
	    else if (sbuf->buff[i+3] == '-')
	      {
			d = sbuf->buff[i+1];
			int idM = d - '0';
			char dd = sbuf->buff[i+2];
			int idd = dd - '0';
			id = (idM*10) + idd;
	      }
	      else if (sbuf->buff[i+4] == '-'){
	      	d = sbuf->buff[i+1];
			int idM = d - '0';

			char dd = sbuf->buff[i+2];
			int idd = dd - '0';

			char ddd = sbuf->buff[i+3];
			int iddd = ddd - '0';

			id = (idM*100) + (idd*10) + iddd;

	      }
	      else if (sbuf->buff[i+5] == '-'){
	      	d = sbuf->buff[i+1];
			int idM = d - '0';

			char dd = sbuf->buff[i+2];
			int idd = dd - '0';

			char ddd = sbuf->buff[i+3];
			int iddd = ddd - '0';

			char dddd = sbuf->buff[i+4];
			int idddd = dddd - '0';

			id = (idM*1000) + (idd*100) + (iddd*10) + idddd;
	      	
	      }
	      else{
	      	 die("Too many alternates ");
	      }
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
	      //from here to the end
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
  strbuf_delete(sbuf_out, 0, first+1);
}


//assume you are going to repeatedly call this on a fasta
//file of known mutations, with the first allele being the susceptible allele, in this format
//>ref_F99I_panel_fasta_sub-3
// the "ref" says it's the susceptible allele. 
//then it has the amino acid substitution (in future will end up with more general, but right now ALL the cases from
//the literature are AA substitutions or indels in genes
//the "sub" says it is a substitution
//then after the - it says how many resistant alleles there are
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
				   int expected_covg, KnownMutation* prev_mut)
				   
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
  if (num_kmers==0)
    {
      return;
    }
   // reset resinfo

  //read- is in this format
  //>ref_F99I_panel_fasta_sub-3-dfrB        so this is >ref means susceptible allele, then identifier of which reference, then sub=susbstitution
  //                                          then 3 means there are 3 possibe resistance alleles all implying F-->Y, and finally the gene name
  //CATGTTTTTATATTTGGAGGGCAAACATTATTTGAAGAAATGATTGATAAAGTGGACGAC
  //>ATA_F99I_panel_fasta_alt-1-dfrB
  //CATGTTTTTATATTTGGAGGGCAAACATTAATAGAAGAAATGATTGATAAAGTGGACGAC
  //ATG
  //>ATT_F99I_panel_fasta_alt-2-dfrB
  //CATGTTTTTATATTTGGAGGGCAAACATTAATTGAAGAAATGATTGATAAAGTGGACGAC
  //ATG
  //>ATC_F99I_panel_fasta_alt-3-dfrB
  //CATGTTTTTATATTTGGAGGGCAAACATTAATCGAAGAAATGATTGATAAAGTGGACGAC
  //ATG

  strbuf_reset(temp_readid_buf);
  strbuf_reset(temp_mut_buf);
  strbuf_append_str(temp_readid_buf, seq->name);
  int r = find_number_resistant_alleles(temp_readid_buf);
  find_gene_name(temp_readid_buf, temp_gene_name_buf);
  GeneMutationGene g  = map_gene_name_str_to_genename(temp_gene_name_buf);
  find_mutation_name(temp_readid_buf, temp_mut_buf);
  KnownMutation km = map_mutation_name_to_enum(temp_mut_buf ,g);

  if (rinfo->var_id!= *prev_mut)
    {
    	// This makes S CIP calls int TB and C00007102,C00009008 fusidic as S  (should be S) and R (should be R) 
      // reset_res_var_info(rinfo);
      // *prev_mut=km; //for use in the next call to this function
      // This makes R CIP calls in TB and C00007102,C00009008 fusidic as r (should be S) and r (should be R)
      *prev_mut=rinfo->var_id; //for use in the next call to this function
      rinfo->working_current_max_sus_allele_present=0;
      rinfo->working_current_max_res_allele_present=0;
    }
  rinfo->gene=g;
  rinfo->var_id=km;
  rinfo->num_resistant_alleles=r;

  //collect min, median covg on allele and also percentage of kmers with any covg
  boolean too_short=false;

  Covg stmp_med =   median_covg_on_allele_in_specific_colour(array_nodes, 
							    num_kmers, 
							    working_ca, 
							    0, 
							    &too_short,
							    ignore_first, ignore_last);

  Covg stmp_min =   min_covg_on_allele_in_specific_colour(array_nodes, 
							  num_kmers, 
							  0, 
							  &too_short,
							  ignore_first, ignore_last);

  int stmp_perc =   percent_nonzero_on_allele_in_specific_colour(array_nodes, 
								 num_kmers, 
								 0, 
								 &too_short,
								 ignore_first, ignore_last);

  if (stmp_perc> rinfo->working_current_max_sus_allele_present)
    {
      rinfo->working_current_max_sus_allele_present=stmp_perc;
      rinfo->susceptible_allele.median_covg = stmp_med;
      rinfo->susceptible_allele.min_covg = stmp_min;
      rinfo->susceptible_allele.percent_nonzero = stmp_perc;
    }

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

      //we have not yet got all the kmers in a resistance allele
      if (rinfo->some_resistant_allele_present==false)
	{
	  too_short=false;
	  
	  int tmp_med = 
	    median_covg_on_allele_in_specific_colour(array_nodes, 
						     num_kmers, 
						     working_ca,
						     0,
						     &too_short,
						     ignore_first, ignore_last);
	  
	  int tmp_min = 
	    min_covg_on_allele_in_specific_colour(array_nodes,
						  num_kmers,
						  0,
						  &too_short,
						  ignore_first, ignore_last);
	  
	  
	  int tmp_perc = 
	    percent_nonzero_on_allele_in_specific_colour(array_nodes,
							 num_kmers,
							 0,
							 &too_short,
							 ignore_first, ignore_last);
	  //if more of the kmers of this version of this mutation
	  //i.e this version of the mutation on this background
	    //are recovered, then keep it - we want to keep the best match
	  if (tmp_perc > rinfo->working_current_max_res_allele_present)
	    {
	      rinfo->resistant_alleles[i].median_covg = tmp_med;
	      rinfo->resistant_alleles[i].min_covg = tmp_min;  
	      rinfo->resistant_alleles[i].percent_nonzero = tmp_perc;
	      //update current best
	      rinfo->working_current_max_res_allele_present = tmp_perc;
	      if (tmp_perc==100)
		{
		  // we have a complete resistance allele, no need to go further
		  rinfo->some_resistant_allele_present=true;
		}
	    }
	}
    }
}
