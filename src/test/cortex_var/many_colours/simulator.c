/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 * 
 * CORTEX project contacts:  
 *    M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 *    Z. Iqbal (zam@well.ox.ac.uk)
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
  simulator.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <CUnit.h>
#include <gsl_sf_gamma.h>
#include <gsl_rng.h>

// cortex_var headers
#include "simulator.h"
#include "db_variants.h"
#include "dB_graph.h"
#include "dB_graph_population.h"
#include "maths.h"
#include "db_complex_genotyping.h"
#include "model_selection.h"
#include "file_reader.h"


void mark_allele_with_all_1s_or_more(dBNode** allele, int len, int colour)
{
  int i;
  for (i=0; i<len; i++)
    {
      db_node_increment_coverage(allele[i], colour);
    }
}

// len is number of edges,
void update_allele(dBNode** allele, int len, int colour, int covg, int eff_read_len)
{

  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
     
  int i;
  //add covg number of reads, placing them randomly on the allele
  for (i = 0; i < covg; i++) 
    {
      int u=0;
      while (u==0)
	{
	  u = gsl_rng_uniform_int(r,len);
	}
      //now add coverage:
      int j;
      for (j=u; (j<u+eff_read_len) && (j<len); j++)
	{
	  db_node_increment_coverage(allele[j], colour);
	}
    }
     
  gsl_rng_free (r);
     

}


void zero_allele(dBNode** allele, int len, int colour_indiv, int colour_allele1, int colour_allele2, int colour_ref_minus_site)
{
  int i;
  for (i=0; i<len; i++)
  {
    db_node_set_coverage(allele[i], colour_indiv, 0);
    db_node_set_coverage(allele[i], colour_allele1, 0);
    db_node_set_coverage(allele[i], colour_allele2, 0);
    db_node_set_coverage(allele[i], colour_ref_minus_site, 0);
  }
}


void zero_path_except_two_alleles_and_ref(dBNode** allele, int len, int colour_allele1, int colour_allele2, int colour_ref_minus_site)
{
  int i,j;
  for (i=0; i<len; i++)
    {
      for (j=0; j<NUMBER_OF_COLOURS; j++)
	{
	  if (!( (j==colour_allele1) || (j==colour_allele2) || (j==colour_ref_minus_site) ))
	    {
	      db_node_set_coverage(allele[i], j, 0);
	      db_node_set_coverage(allele[i], j, 0);
	      db_node_set_coverage(allele[i], j, 0);
	      db_node_set_coverage(allele[i], j, 0);
	    }
	}
    }
}

void test(VariantBranchesAndFlanks* var, //MultiplicitiesAndOverlapsOfBiallelicVariant* var_mults, 
      GraphAndModelInfo* model_info, char* fasta,
      int colour_ref_minus_site, int colour_indiv,
      boolean use_1_and_2net, int working_colour1, int working_colour2,
      int* count_passes, int* count_fails, dBGraph* db_graph,
      char* true_ml_gt_name)
{

  int max_allele_length=100000;
  if ( (var->len_one_allele>max_allele_length)||(var->len_other_allele > max_allele_length) )
  {
    die("simulator.c exited!");
  }

  int colours_to_genotype[]={colour_indiv};

  char ml_genotype_name[100]="";
  char* ml_genotype_name_array[1];
  ml_genotype_name_array[0]=ml_genotype_name;
  char ml_but_one_genotype_name[100]="";
  char* ml_but_one_genotype_name_array[1];
  ml_but_one_genotype_name_array[0]=ml_but_one_genotype_name;

  int MIN_LLK = -999999999;
  double ml_genotype_lik=MIN_LLK;
  double ml_but_one_genotype_lik=MIN_LLK;

  calculate_max_and_max_but_one_llks_of_specified_set_of_genotypes_of_complex_site(
    colours_to_genotype, 1,
    colour_ref_minus_site, 2, 1,3,
    max_allele_length, fasta,
    AssumeUncleaned,
    &ml_genotype_lik, &ml_but_one_genotype_lik,
    ml_genotype_name_array, ml_but_one_genotype_name_array,
    false, model_info, db_graph, working_colour1, working_colour2,
    use_1_and_2net, use_1_and_2net, MIN_LLK);

  //printf("We get max lik gt %s and we expect %s\n", ml_genotype_name, true_ml_gt_name);
  if (strcmp(ml_genotype_name, true_ml_gt_name)==0) 
  {
    (*count_passes)++;
  }
  else
  {
    (*count_fails)++;
  }
}

//if you want to simulate a true hom, pass in var with both alleles the same
void simulator(int depth, int read_len, int kmer, double seq_err_per_base,
               int number_repetitions,  int colour_indiv,
               int colour_allele1, int colour_allele2, int colour_ref_minus_site,
               VariantBranchesAndFlanks* var,
               int len_genome_minus_site, zygosity true_gt,
               GraphAndModelInfo* model_info,
               char* fasta, char* true_ml_gt_name,
               int working_colour1, int working_colour2,
               boolean using_1and2_nets,
               dBGraph* db_graph)
               //dBNode** genome_minus_site
               //boolean are_the_two_alleles_identical
               //char* filelist_net1, char* filelist_net2
               //int working_colour_1net, int working_colour_2net
{

  if (NUMBER_OF_COLOURS<4)
    {
      die("Cannot run the simulator with <4 colours. Recompile.\n");
    }

  int count_passes = 0;
  int count_fails = 0;

  const gsl_rng_type * T;
  gsl_rng * r;
  
  // GLS setup:
  /* create a generator chosen by the 
     environment variable GSL_RNG_TYPE */
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);


  //put the alleles  and reference into their colours:
  mark_allele_with_all_1s_or_more(var->one_allele, var->len_one_allele,     colour_allele1);
  mark_allele_with_all_1s_or_more(var->other_allele, var->len_other_allele, colour_allele2);
  
  int i;
  for (i=0; i<number_repetitions; i++)
    {
      zero_path_except_two_alleles_and_ref(var->one_allele, var->len_one_allele, colour_allele1, colour_allele2, colour_ref_minus_site);
      zero_path_except_two_alleles_and_ref(var->other_allele, var->len_other_allele, colour_allele1, colour_allele2, colour_ref_minus_site);

      //give the each allele depth which is taken from a Poisson with mean =  (D/2) * (R-k+1)/R  * (1-k*epsilon)
      //printf("Depth %d, var->len one allele - k,er +1 = %d, and 1-kmer * seq_err = %f     \n", depth, (var->len_one_allele)-kmer+1, 1-kmer*seq_err_per_base    )    ;
      double exp_depth_on_allele1 = 0;
      if (true_gt==het)
	//het. So 1/3 of seq errors on the other allele end up on this one
      	{
	  exp_depth_on_allele1 = ((double) depth/2) *
	    ( (double)(read_len+(var->len_one_allele)-kmer+1)/read_len) * (1-kmer*seq_err_per_base);
	    // + ((double) depth/2) *
	    //( (double)(read_len+(var->len_other_allele)-kmer+1)/read_len) * (kmer*seq_err_per_base/3); //some errors from the other allele give covg here
	}
      else if (true_gt==hom_one)
	{//hom
	  exp_depth_on_allele1 = ((double) depth) * 
	    ( (double)(read_len+(var->len_one_allele)-kmer+1)/read_len) * (1-kmer*seq_err_per_base);
	}
      else if (true_gt==hom_other)
	{
	  exp_depth_on_allele1 = ((double) depth) *
	    ( (double)(read_len+(var->len_other_allele)-kmer+1)/read_len) * kmer*seq_err_per_base/3;// 1/3 of the errors on the true allele are on this one
	}
      double exp_depth_on_allele2 = 0;
      if (true_gt==het)
	//het. So 1/3 of seq errors are not a problem. So loss of covg is (1-k*(2/3)*e)
      	{
	  exp_depth_on_allele2 = exp_depth_on_allele1;
	}
      else if (true_gt==hom_one)
	{
	  exp_depth_on_allele2 = ( (double)(read_len+(var->len_one_allele)-kmer+1)/read_len) * kmer*seq_err_per_base/3;
	}
      else if (true_gt==hom_other)
	{
	  exp_depth_on_allele2 = ((double) depth) *
	    ( (double)(read_len+(var->len_other_allele)-kmer+1)/read_len) * (1-kmer*seq_err_per_base);
	}

      double exp_depth_on_ref_minus_site = (double) depth * ((double)(len_genome_minus_site-kmer+1)/read_len) * (1-kmer*seq_err_per_base);
      //printf("exp ZAMMER %f %f %f\n", exp_depth_on_allele1, exp_depth_on_allele2, exp_depth_on_ref_minus_site);
      unsigned int sampled_covg_allele1 = gsl_ran_poisson (r, exp_depth_on_allele1);
      unsigned int sampled_covg_allele2 = gsl_ran_poisson (r, exp_depth_on_allele2);
      unsigned int sampled_covg_rest_of_genome = gsl_ran_poisson (r, exp_depth_on_ref_minus_site);
      printf("Sampled covgs on alleles 1,2 and genome are  %d %d %d\n", sampled_covg_allele1, sampled_covg_allele2, sampled_covg_rest_of_genome);
      update_allele(var->one_allele, var->len_one_allele,     
		    colour_indiv, sampled_covg_allele1,read_len-kmer+1);
      update_allele(var->other_allele, var->len_other_allele, 
		    colour_indiv, sampled_covg_allele2, read_len-kmer+1);
      //update_allele(genome_minus_site,len_genome_minus_site,  colour_indiv, sampled_covg_rest_of_genome);

      test(var, model_info, fasta, colour_ref_minus_site,
	         colour_indiv, using_1and2_nets,
           working_colour1, working_colour2,
           &count_passes, &count_fails, db_graph, true_ml_gt_name);

    }

  //cleanup
  zero_allele(var->one_allele, var->len_one_allele,     colour_indiv, colour_allele1, colour_allele2, colour_ref_minus_site);
  zero_allele(var->other_allele, var->len_other_allele, colour_indiv, colour_allele1, colour_allele2, colour_ref_minus_site);

  CU_ASSERT((double)count_passes/(double)(count_passes+count_fails) > 0.9 );//actually, we could set this to ==1
  printf("Number of passes: %d, number of fails %d\n", count_passes, count_fails);


  gsl_rng_free (r);

}




