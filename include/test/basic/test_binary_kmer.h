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
  test_binary_kmer.h
*/

#ifndef TEST_BIN_KMER_H_
#define TEST_BIN_KMER_H_

#include "global.h"

void test_that_bitfield_really_is_64bits();

void test_binary_kmer_assignment_operator();

void test_binary_kmer_comparison_operator();

void test_binary_kmer_less_than();

void test_binary_kmer_right_shift_one_base();

void test_binary_kmer_left_shift_one_base();

void test_seq_to_binary_kmer_and_binary_kmer_to_seq();

void test_binary_kmer_nucleotide_iterator();

void test_binary_kmer_reverse_complement();

void test_seq_reverse_complement();

void test_get_sliding_windows_from_sequence();

void test_breaking_homopolymers_in_get_sliding_windows();

#endif /* TEST_BIN_KMER_H_ */
