/*
 * 
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 *
 */
/*
  test_db_complex_genotyping.h
*/

#ifndef TEST_DB_COMPLEX_GT_H_
#define TEST_DB_COMPLEX_GT_H_

#include "global.h"

void test_initialise_multiplicities_of_allele_nodes_wrt_both_alleles();
void test_calc_log_likelihood_of_genotype_with_complex_alleles1();
void test_calc_log_likelihood_of_genotype_with_complex_alleles2();
void test_calc_log_likelihood_of_genotype_with_complex_alleles3();
void regression_test_1_single_bubble_call_one_allele_shorter_than_k_one_very_long();
void regression_test_2_genotyping_of_PD_SNP_call();

#endif /* TEST_DB_COMPLEX_GT_H_ */
