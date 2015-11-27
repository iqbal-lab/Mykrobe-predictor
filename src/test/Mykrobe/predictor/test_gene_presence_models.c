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
#include "gene_presence_models.h"
#include "gene_presence.h"
#include "mut_models.h"
void test_resistotype_gene()
{
	GeneInfo* gi = alloc_and_init_gene_info();
	gi->median_covg = 45;
	gi->median_covg_on_nonzero_nodes = 45;
	gi->percent_nonzero = 98;
	gi->num_gaps = 3;
	gi->len = 1993;
	double err_rate = 0.01;
	int kmer = 15;
	int expected_covg = 45;
	Model best_model;
	ModelChoiceMethod choice = MaxAPosteriori;
	int min_expected_kmer_recovery_for_this_gene = 80;

	double lambda_g =expected_covg;
	double lambda_e = expected_covg*err_rate;

	double epsilon = pow(1-err_rate, kmer);

	  // printf("epsilon %f \n", epsilon);  
	  // printf("lambda_g %f \n", lambda_g);   
	  // printf("lambda_e %f \n", lambda_e);   
	  // printf("err_rate %f \n", err_rate);   
	  
	  // printf("expected_covg %d \n", expected_covg);   
	  // printf("kmer_size %d \n", kmer);   
	  // printf("median %d \n", gi->median_covg);   
	  // printf("percent_nonzero %d \n", gi->percent_nonzero);   
	  // printf("len %d \n", gi->len);   
	  // printf("num_gaps %d \n", gi->num_gaps);   


	boolean genotyped_present = false;
	InfectionType I = resistotype_gene(gi, err_rate, kmer,
			       lambda_g,  lambda_e, epsilon, expected_covg, 0,
			       &best_model,
			       choice,
			       min_expected_kmer_recovery_for_this_gene,0.03,
			       &genotyped_present);
    

	CU_ASSERT(I == Resistant);
	free_gene_info(gi);
}

void test_resistotype_unsure_gene()
{


	GeneInfo* gi = alloc_and_init_gene_info();
	gi->median_covg = 2;
	gi->median_covg_on_nonzero_nodes = 2;
	gi->percent_nonzero = 82;
	gi->num_gaps = 25;
	gi->len = 1993;
	double err_rate = 0.01;
	int kmer = 15;
	int expected_covg = 100;
	Model best_model;
	ModelChoiceMethod choice = MaxAPosteriori;
	int min_expected_kmer_recovery_for_this_gene = 80;

	double lambda_g =expected_covg;
	double lambda_e = expected_covg*err_rate;

	double epsilon = pow(1-err_rate, kmer);

	boolean genotyped_present = false;
	InfectionType I = resistotype_gene(gi, err_rate, kmer,
			       lambda_g,  lambda_e, epsilon, expected_covg, 0,
			       &best_model,
			       choice,
			       min_expected_kmer_recovery_for_this_gene,0.03,
			       &genotyped_present);
    
	printf("%i\n",I );
	CU_ASSERT(I == Susceptible);
	CU_ASSERT(genotyped_present == true);
	free_gene_info(gi);
}


void test_resistotype_unsure_gene_2()
{


	GeneInfo* gi = alloc_and_init_gene_info();
	gi->median_covg = 4;
	gi->median_covg_on_nonzero_nodes = 4;
	gi->percent_nonzero = 82;
	gi->num_gaps = 25;
	gi->len = 1993;
	double err_rate = 0.01;
	int kmer = 15;
	int expected_covg = 100;
	Model best_model;
	ModelChoiceMethod choice = MaxAPosteriori;
	int min_expected_kmer_recovery_for_this_gene = 80;

	double lambda_g =expected_covg;
	double lambda_e = expected_covg*err_rate;

	double epsilon = pow(1-err_rate, kmer);

	boolean genotyped_present = false;
	InfectionType I = resistotype_gene(gi, err_rate, kmer,
			       lambda_g,  lambda_e, epsilon, expected_covg, 0,
			       &best_model,
			       choice,
			       min_expected_kmer_recovery_for_this_gene,0.03,
			       &genotyped_present);
    
	printf("%i\n",I );
	CU_ASSERT(I == MixedInfection);
	CU_ASSERT(genotyped_present == true);
	free_gene_info(gi);
}


void test_resistotype_minor_gene()
{


	GeneInfo* gi = alloc_and_init_gene_info();
	gi->median_covg = 20;
	gi->median_covg_on_nonzero_nodes = 20;
	gi->percent_nonzero = 82;
	gi->num_gaps = 25;
	gi->len = 1993;
	double err_rate = 0.01;
	int kmer = 15;
	int expected_covg = 100;
	Model best_model;
	ModelChoiceMethod choice = MaxAPosteriori;
	int min_expected_kmer_recovery_for_this_gene = 80;

	double lambda_g =expected_covg;
	double lambda_e = expected_covg*err_rate;

	double epsilon = pow(1-err_rate, kmer);

	boolean genotyped_present = false;
	InfectionType I = resistotype_gene(gi, err_rate, kmer,
			       lambda_g,  lambda_e, epsilon, expected_covg, 0,
			       &best_model,
			       choice,
			       min_expected_kmer_recovery_for_this_gene,0.03,
			       &genotyped_present);
    

	CU_ASSERT(I == MixedInfection);
	free_gene_info(gi);


}

void test_resistotype_gene_at_high_CN()
{
	GeneInfo* gi = alloc_and_init_gene_info();
	gi->median_covg = 100;
	gi->median_covg_on_nonzero_nodes = 100;
	gi->percent_nonzero = 98;
	gi->num_gaps = 3;
	gi->len = 1993;
	double err_rate = 0.01;
	int kmer = 15;
	int expected_covg = 10;
	Model best_model;
	ModelChoiceMethod choice = MaxAPosteriori;
	int min_expected_kmer_recovery_for_this_gene = 80;

	double lambda_g =expected_covg;
	double lambda_e = expected_covg*err_rate;

	double epsilon = pow(1-err_rate, kmer);

	boolean genotyped_present = false;
	InfectionType I = resistotype_gene(gi, err_rate, kmer,
			       lambda_g,  lambda_e, epsilon, expected_covg, 0,
			       &best_model,
			       choice,
			       min_expected_kmer_recovery_for_this_gene,0.03,
			       &genotyped_present);
    

	CU_ASSERT(I == Resistant);
	free_gene_info(gi);
}

void test_resistotype_gene_S()
{
	// How it might look with repeat kmers
	GeneInfo* gi = alloc_and_init_gene_info();
	gi->median_covg = 0;
	gi->median_covg_on_nonzero_nodes = 250;
	gi->percent_nonzero = 12;
	gi->num_gaps = 500;
	gi->len = 1993;
	double err_rate = 0.01;
	int kmer = 15;
	int expected_covg = 100;
	Model best_model;
	ModelChoiceMethod choice = MaxAPosteriori;
	int min_expected_kmer_recovery_for_this_gene = 80;

	double lambda_g =expected_covg;
	double lambda_e = expected_covg*err_rate;

	double epsilon = pow(1-err_rate, kmer);

	boolean genotyped_present = false;
	InfectionType I = resistotype_gene(gi, err_rate, kmer,
			       lambda_g,  lambda_e, epsilon, expected_covg, 0,
			       &best_model,
			       choice,
			       min_expected_kmer_recovery_for_this_gene,0.03,
			       &genotyped_present);
    

	CU_ASSERT(I == Susceptible);
	CU_ASSERT(genotyped_present == false);
	free_gene_info(gi);

}


void test_low_coverage_genes()
{
	// How it might look with repeat kmers
	GeneInfo* gi = alloc_and_init_gene_info();
	gi->median_covg = 16;
	gi->median_covg_on_nonzero_nodes = 16;
	gi->percent_nonzero = 2;
	gi->len = 1993;
	double err_rate = 0.01;
	int kmer = 15;
	int expected_covg = 10;
	Model best_model;
	ModelChoiceMethod choice = MaxAPosteriori;
	int min_expected_kmer_recovery_for_this_gene = 80;

	double lambda_g =expected_covg;
	double lambda_e = expected_covg*err_rate;

	double epsilon = pow(1-err_rate, kmer);

	boolean genotyped_present = false;
	InfectionType I = resistotype_gene(gi, err_rate, kmer,
			       lambda_g,  lambda_e, epsilon, expected_covg, 0,
			       &best_model,
			       choice,
			       min_expected_kmer_recovery_for_this_gene,0.03,
			       &genotyped_present);
	CU_ASSERT(I == Susceptible);
	CU_ASSERT(genotyped_present == false);

	gi->median_covg = 16;
	gi->median_covg_on_nonzero_nodes = 16;
	gi->percent_nonzero = 80;
	gi->len = 1993;

	genotyped_present = false;
	I = resistotype_gene(gi, err_rate, kmer,
			       lambda_g,  lambda_e, epsilon, expected_covg, 0,
			       &best_model,
			       choice,
			       min_expected_kmer_recovery_for_this_gene,0.03,
			       &genotyped_present);
	CU_ASSERT(I == Resistant);
	CU_ASSERT(genotyped_present == true);
	free_gene_info(gi);
}


void test_low_coverage_ont_genes()
{
	// How it might look with repeat kmers
	GeneInfo* gi = alloc_and_init_gene_info();
	gi->median_covg = 0;
	gi->median_covg_on_nonzero_nodes = 0;
	gi->percent_nonzero = 0;
	gi->len = 1993;
	double err_rate = 0.1;
	int kmer = 15;
	int expected_covg = 2;
	Model best_model;
	ModelChoiceMethod choice = MaxAPosteriori;
	int min_expected_kmer_recovery_for_this_gene = 30;

	double lambda_g =expected_covg;
	double lambda_e = expected_covg*err_rate;

	double epsilon = pow(1-err_rate, kmer);

	boolean genotyped_present = false;
	InfectionType I = resistotype_gene(gi, err_rate, kmer,
			       lambda_g,  lambda_e, epsilon, expected_covg, 0,
			       &best_model,
			       choice,
			       min_expected_kmer_recovery_for_this_gene,0.03,
			       &genotyped_present);
	CU_ASSERT(I == Susceptible);
	CU_ASSERT(genotyped_present == false);

	gi->median_covg = 2;
	gi->median_covg_on_nonzero_nodes = 2;
	gi->percent_nonzero = 50;
	gi->len = 1993;

	genotyped_present = false;
	I = resistotype_gene(gi, err_rate, kmer,
			       lambda_g,  lambda_e, epsilon, expected_covg, 0,
			       &best_model,
			       choice,
			       min_expected_kmer_recovery_for_this_gene,0.03,
			       &genotyped_present);
	CU_ASSERT(I == Resistant);
	CU_ASSERT(genotyped_present == true);

	min_expected_kmer_recovery_for_this_gene = 50;

	gi->median_covg = 1;
	gi->median_covg_on_nonzero_nodes = 1;
	gi->percent_nonzero = 44;
	gi->len = 1600;

	genotyped_present = false;
	I = resistotype_gene(gi, err_rate, kmer,
			       lambda_g,  lambda_e, epsilon, expected_covg, 0,
			       &best_model,
			       choice,
			       min_expected_kmer_recovery_for_this_gene,0.03,
			       &genotyped_present);
	CU_ASSERT(I == Resistant);
	CU_ASSERT(genotyped_present == true);
	free_gene_info(gi);

}
