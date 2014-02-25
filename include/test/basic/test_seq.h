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
