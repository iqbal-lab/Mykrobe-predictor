/*
 * 
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 *
 */
/*
  test_seq.h
*/

#ifndef TEST_SEQ_H_
#define TEST_SEQ_H_

#include "global.h"

void test_read_sequence_from_fasta();

void test_read_sequence_from_fasta_when_file_has_bad_reads();

void test_read_sequence_from_fastq();

void test_read_sequence_from_long_fasta();
void test_shift_last_kmer_to_start_of_sequence();

void test_read_sequence_from_fastq_with_bad_reads_and_long_reads();

#endif /* TEST_SEQ_H_ */
