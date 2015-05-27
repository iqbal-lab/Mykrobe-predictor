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

	double genome_size = 280000;
	double mean_read_length = 100;
	double bp_loaded = 28000000;
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


	InfectionType I = resistotype_gene(gi, err_rate, kmer,
			       lambda_g,  lambda_e, epsilon, expected_covg,
			       &best_model,
			       choice,
			       min_expected_kmer_recovery_for_this_gene);
    

	CU_ASSERT(I == Resistant);

}

void test_resistotype_minor_gene()
{
// typedef struct
// {
//   Covg median_covg;
//   Covg min_covg;
//   Covg median_covg_on_nonzero_nodes;
//   int num_gaps;
//   int len;
//   int  percent_nonzero;
//   StrBuf* strbuf;
//   GenePresenceGene name;
// } GeneInfo;
// """mecA"" :{
//    ""per_cov"": ""82"",
//    ""median_cov"": ""2"",
//    ""conf"": ""1171"",
//   ""induced_resistance"": ""Methicillin""
//   },"

	GeneInfo* gi = alloc_and_init_gene_info();
	gi->median_covg = 2;
	gi->median_covg_on_nonzero_nodes = 2;
	gi->percent_nonzero = 82;
	gi->num_gaps = 25;
	gi->len = 1993;
	double err_rate = 0.01;
	int kmer = 15;
	int expected_covg = 45;
	Model best_model;
	ModelChoiceMethod choice = MaxAPosteriori;
	int min_expected_kmer_recovery_for_this_gene = 80;

	double genome_size = 280000;
	double mean_read_length = 100;
	double bp_loaded = 28000000;
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
	 // epsilon 0.860058
	  // lambda_g 0.860058
	// lambda_e 0.002491
	// err_rate 0.010000
	// expected_covg 100
	// kmer_size 15
	// median 2
	// percent_nonzero 82

// LLks of S, M, R are -2298.969819, -1132.828561 and -1127.381016
                // "Methicillin": "R",

	InfectionType I = resistotype_gene(gi, err_rate, kmer,
			       lambda_g,  lambda_e, epsilon, expected_covg,
			       &best_model,
			       choice,
			       min_expected_kmer_recovery_for_this_gene);
    

	CU_ASSERT(I == MixedInfection);



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

	double genome_size = 280000;
	double mean_read_length = 100;
	double bp_loaded = 28000000;
	double lambda_g =expected_covg;
	double lambda_e = expected_covg*err_rate;

	double epsilon = pow(1-err_rate, kmer);

	InfectionType I = resistotype_gene(gi, err_rate, kmer,
			       lambda_g,  lambda_e, epsilon, expected_covg,
			       &best_model,
			       choice,
			       min_expected_kmer_recovery_for_this_gene);
    

	CU_ASSERT(I == Resistant);

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

	double genome_size = 280000;
	double mean_read_length = 100;
	double bp_loaded = 28000000;
	double lambda_g =expected_covg;
	double lambda_e = expected_covg*err_rate;

	double epsilon = pow(1-err_rate, kmer);

	InfectionType I = resistotype_gene(gi, err_rate, kmer,
			       lambda_g,  lambda_e, epsilon, expected_covg,
			       &best_model,
			       choice,
			       min_expected_kmer_recovery_for_this_gene);
    

	CU_ASSERT(I == Susceptible);

}
