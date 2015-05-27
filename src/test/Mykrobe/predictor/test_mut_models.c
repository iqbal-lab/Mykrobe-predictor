/*
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 *
 */
/*
  test_genotype_known.c 
*/

// system headers
#include <stdlib.h>
#include <limits.h>

// third party headers
#include <CUnit.h>
#include <Basic.h>
#include <string_buffer.h>

#include "build.h"
#include "element.h"
#include "seq.h"
#include "open_hash/hash_table.h"
#include "genotyping_known.h"
#include "mut_models.h"
#include "antibiotics.h"

void test_mutation_S()
{

// typedef struct
// {
//   VarOnBackground* vob_best_sus;
//   VarOnBackground* vob_best_res;
// } Var;//corresponds to an enum - so encapsulates info across backgrounds

// typedef struct
// {
//   AlleleInfo susceptible_allele;
//   AlleleInfo resistant_alleles[60];
//   int num_resistant_alleles;
//   boolean some_resistant_allele_present;
//   int working_current_max_res_allele_present;
//   int working_current_median_covg;
//   //  int working_current_max_sus_allele_present;
//   KnownMutation var_id;
//   GeneMutationGene gene;
// }VarOnBackground;
// typedef struct
// {
//   Covg median_covg;
//   Covg median_covg_on_nonzero_nodes;
//   Covg min_covg;
//   int  percent_nonzero;

// } AlleleInfo;
	double err_rate = 0.01;
	int kmer = 15;
	int expected_covg = 45;
	
	float min_frac_to_detect_minor_pops = 0.1;
	double genome_size = 280000;
	double mean_read_length = 100;
	double bp_loaded = 28000000;
	double lambda_g =expected_covg;
	double lambda_e = expected_covg*err_rate;
	double epsilon = pow(1-err_rate, kmer);	

    AlleleInfo*	susceptible_allele = alloc_allele_info();
	susceptible_allele->median_covg= expected_covg;
	susceptible_allele->median_covg_on_nonzero_nodes= expected_covg;
	susceptible_allele->min_covg= expected_covg;
	susceptible_allele->percent_nonzero = 100;

	AlleleInfo* resistant_allele = alloc_allele_info();	
	resistant_allele->median_covg = 0;
	resistant_allele->median_covg_on_nonzero_nodes = 0;
	resistant_allele->min_covg = 0;
	resistant_allele->percent_nonzero = 0;


	VarOnBackground*  vob_best = alloc_and_init_var_on_background();
	vob_best->susceptible_allele = *susceptible_allele;
	vob_best->resistant_alleles[0] = *resistant_allele;
	vob_best->num_resistant_alleles = 1;
	vob_best->some_resistant_allele_present = false;
	vob_best->working_current_max_res_allele_present = 0;

	Var* var = alloc_var();
	var->vob_best_res = vob_best;
	var->vob_best_sus = vob_best;

	Var** vars =  malloc(sizeof(Var*)*1);
	vars[0] = var;

	ModelChoiceMethod choice = MaxAPosteriori;
    Model best_model;

     InfectionType I=
	resistotype(vars[0], err_rate, kmer, 
		    lambda_g, lambda_e, epsilon,
		    &best_model, MaxAPosteriori,
		    min_frac_to_detect_minor_pops);


      if (I==Susceptible)
	{
	  printf("S\n");
	}
      else if (I==MixedInfection)
	{
	  printf("r\n");
	}
      else if (I==Resistant)
	{
	  printf("R\n");
	}
      else
	{
	  printf("N\n");
	}	

	CU_ASSERT(I == Susceptible)
}
void test_mutation_R()
{
	double err_rate = 0.01;
	int kmer = 15;
	int expected_covg = 45;
	
	float min_frac_to_detect_minor_pops = 0.1;
	double genome_size = 280000;
	double mean_read_length = 100;
	double bp_loaded = 28000000;
	double lambda_g =expected_covg;
	double lambda_e = expected_covg*err_rate;
	double epsilon = pow(1-err_rate, kmer);	

    AlleleInfo*	resistant_allele = alloc_allele_info();
	resistant_allele->median_covg= expected_covg;
	resistant_allele->median_covg_on_nonzero_nodes= expected_covg;
	resistant_allele->min_covg= expected_covg;
	resistant_allele->percent_nonzero = 100;

	AlleleInfo* susceptible_allele = alloc_allele_info();	
	susceptible_allele->median_covg = 0;
	susceptible_allele->median_covg_on_nonzero_nodes = 0;
	susceptible_allele->min_covg = 0;
	susceptible_allele->percent_nonzero = 0;


	VarOnBackground*  vob_best = alloc_and_init_var_on_background();
	vob_best->susceptible_allele = *susceptible_allele;
	vob_best->resistant_alleles[0] = *resistant_allele;
	vob_best->num_resistant_alleles = 1;
	vob_best->some_resistant_allele_present = false;
	vob_best->working_current_max_res_allele_present = 0;

	Var* var = alloc_var();
	var->vob_best_res = vob_best;
	var->vob_best_sus = vob_best;

	Var** vars =  malloc(sizeof(Var*)*1);
	vars[0] = var;

	ModelChoiceMethod choice = MaxAPosteriori;
    Model best_model;

     InfectionType I=
	resistotype(vars[0], err_rate, kmer, 
		    lambda_g, lambda_e, epsilon,
		    &best_model, MaxAPosteriori,
		    min_frac_to_detect_minor_pops);


      if (I==Susceptible)
	{
	  printf("S\n");
	}
      else if (I==MixedInfection)
	{
	  printf("r\n");
	}
      else if (I==Resistant)
	{
	  printf("R\n");
	}
      else
	{
	  printf("N\n");
	}	

	CU_ASSERT(I == Resistant)}

void test_mutation_r()
{
	double err_rate = 0.01;
	int kmer = 15;
	int expected_covg = 45;
	
	float min_frac_to_detect_minor_pops = 0.1;
	double genome_size = 280000;
	double mean_read_length = 100;
	double bp_loaded = 28000000;
	double lambda_g =expected_covg;
	double lambda_e = expected_covg*err_rate;
	// double lambda_g = 0.572547;

	// double lambda_e = lambda_g * err_rate/3 * pow(1-err_rate, kmer-1);
	// double lambda_e = 0.001928;	
	double epsilon = pow(1-err_rate, kmer);	

    AlleleInfo*	resistant_allele = alloc_allele_info();
	resistant_allele->median_covg= 0.2 * expected_covg;
	resistant_allele->median_covg_on_nonzero_nodes= 0.2 * expected_covg;
	resistant_allele->min_covg= 0.2 * expected_covg;
	resistant_allele->percent_nonzero = 100;

	AlleleInfo* susceptible_allele = alloc_allele_info();	
	susceptible_allele->median_covg= 0.8 * expected_covg;
	susceptible_allele->median_covg_on_nonzero_nodes= 0.8 * expected_covg;
	susceptible_allele->min_covg= 0.8 * expected_covg;
	susceptible_allele->percent_nonzero = 100;


	VarOnBackground*  vob_best = alloc_and_init_var_on_background();
	vob_best->susceptible_allele = *susceptible_allele;
	vob_best->resistant_alleles[0] = *resistant_allele;
	vob_best->num_resistant_alleles = 1;
	vob_best->some_resistant_allele_present = true;
	vob_best->working_current_max_res_allele_present = expected_covg/2;

	Var* var = alloc_var();
	var->vob_best_res = vob_best;
	var->vob_best_sus = vob_best;

	Var** vars =  malloc(sizeof(Var*)*1);
	vars[0] = var;

	ModelChoiceMethod choice = MaxAPosteriori;
    Model best_model;

     InfectionType I=
	resistotype(vars[0], err_rate, kmer, 
		    lambda_g, lambda_e, epsilon,
		    &best_model, MaxAPosteriori,
		    min_frac_to_detect_minor_pops);


      if (I==Susceptible)
	{
	  printf("S\n");
	}
      else if (I==MixedInfection)
	{
	  printf("r\n");
	}
      else if (I==Resistant)
	{
	  printf("R\n");
	}
      else
	{
	  printf("N\n");
	}	

	CU_ASSERT(I == MixedInfection)
}

                // "gyrA_S84L" :{
                //         "R_per_cov": "100",
                //         "S_per_cov": "62",
                //         "R_median_cov": "3",
                //         "S_median_cov": "2",
                //         "conf": "26",
                // "induced_resistance": "Ciprofloxacin"

// {susceptible_allele = {median_covg = 2, median_covg_on_nonzero_nodes = 0, min_covg = 0, percent_nonzero = 62}, resistant_alleles = {{median_covg = 3, median\
// _covg_on_nonzero_nodes = 0, min_covg = 2, percent_nonzero = 100}, {median_covg = 0, median_covg_on_nonzero_nodes = 0, min_covg = 0, percent_nonzero = 6}, {median_\
// covg = 0, median_covg_on_nonzero_nodes = 0, min_covg = 0, percent_nonzero = 12}, {median_covg = 0, median_covg_on_nonzero_nodes = 0, min_covg = 0, percent_nonzero\
//  = 0}, {median_covg = 0, median_covg_on_nonzero_nodes = 0, min_covg = 0, percent_nonzero = 0}, {median_covg = 0, median_covg_on_nonzero_nodes = 0, min_covg = 0, p\
// ercent_nonzero = 12}, {median_covg = 0, median_covg_on_nonzero_nodes = 0, min_covg = 0, percent_nonzero = 0} <repeats 54 times>}, num_resistant_alleles = 6, some_\
// resistant_allele_present = 1 '\001', working_current_max_res_allele_present = 100, working_current_median_covg = 3, var_id = gyrA_S84L, gene = gyrA}

void test_mutation_custom_1()
{
	double err_rate = 0.01;
	int kmer = 15;
	int expected_covg = 252;
	
	float min_frac_to_detect_minor_pops = 0.1;
	double genome_size = 280000;
	double mean_read_length = 100;
	double bp_loaded = 28000000;
	// double lambda_g =expected_covg;
	// double lambda_e = expected_covg*err_rate / 3;
	double lambda_g = 3.406567;
	double lambda_e = 0.011470;
	double epsilon = pow(1-err_rate, kmer);	

    AlleleInfo*	resistant_allele = alloc_allele_info();
	resistant_allele->median_covg= 3;
	resistant_allele->median_covg_on_nonzero_nodes= 0;
	resistant_allele->min_covg= 2;
	resistant_allele->percent_nonzero = 100;

	AlleleInfo* susceptible_allele = alloc_allele_info();	
	susceptible_allele->median_covg= 2;
	susceptible_allele->median_covg_on_nonzero_nodes= 0;
	susceptible_allele->min_covg= 0;
	susceptible_allele->percent_nonzero = 62;


	VarOnBackground*  vob_best = alloc_and_init_var_on_background();
	vob_best->susceptible_allele = *susceptible_allele;
	vob_best->resistant_alleles[0] = *resistant_allele;
	vob_best->num_resistant_alleles = 1;
	vob_best->some_resistant_allele_present = true;
	vob_best->working_current_max_res_allele_present = 3;

	Var* var = alloc_var();
	var->vob_best_res = vob_best;
	var->vob_best_sus = vob_best;

	Var** vars =  malloc(sizeof(Var*)*1);
	vars[0] = var;

	ModelChoiceMethod choice = MaxAPosteriori;
    Model best_model;

     InfectionType I= resistotype(vars[0], err_rate, kmer, 
		    lambda_g, lambda_e, epsilon,
		    &best_model, MaxAPosteriori,
		    min_frac_to_detect_minor_pops);
	printf("conf %f\n", best_model.conf);
// err_rate 0.010000
// lambda_g 3.406567
// lambda_e 0.011470
//  epsilon 0.860058
// min_frac_to_detect_minor_pops 0.100000
//  LLks of S, M, R are -16.855595, -6.881425 and -3.398583
// LP of S, M, R are -33.711189, -10000005.881425 and -6.797167/////
      if (I==Susceptible)
	{
	  printf("S\n");
	}
      else if (I==MixedInfection)
	{
	  printf("r\n");
	}
      else if (I==Resistant)
	{
	  printf("R\n");
	}
      else
	{
	  printf("N\n");
	}	

	CU_ASSERT(I == Unsure)
}
