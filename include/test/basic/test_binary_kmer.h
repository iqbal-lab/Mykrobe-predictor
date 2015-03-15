/*
 * 
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 *
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
