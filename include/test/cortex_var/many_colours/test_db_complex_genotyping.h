/*
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
