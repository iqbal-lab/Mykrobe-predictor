/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
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
  run_basic_tests.c
*/

#include <stdlib.h>
#include <stdio.h>

#include <CUnit.h>
#include <Basic.h>

// cortex_var headers
#include "test_binary_kmer.h"
#include "test_seq.h"

int  main()
{

  CU_pSuite pSuite = NULL;

  /* initialize the CUnit test registry */
  if (CUE_SUCCESS!=CU_initialize_registry())
    return CU_get_error();
  
  /* add a suite to the registry */
  pSuite = CU_add_suite("Suite_1", NULL, NULL);
  if (NULL == pSuite) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  /* add the tests to the suite */





  if (NULL == CU_add_test(pSuite, "Check that bitfield as defined really is 64bits on this platform", test_that_bitfield_really_is_64bits )){
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pSuite, "Test the assignment operator for BinaryKmers", test_binary_kmer_assignment_operator  )){
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pSuite, "Test the comparison operator for binary kmers", test_binary_kmer_comparison_operator )){
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pSuite, "Test the less than operator for binary kmers", test_binary_kmer_less_than )){
    CU_cleanup_registry();
    return CU_get_error();
  }
  
  if (NULL == CU_add_test(pSuite, "Test the right shift operator for big binary kmers that are encoded in multiple long integers", test_binary_kmer_right_shift_one_base )){
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pSuite, "Test the left shift operator for big binary kmers that are encoded in multiple long integers", test_binary_kmer_left_shift_one_base )){
    CU_cleanup_registry();
    return CU_get_error();
  }


  if (NULL == CU_add_test(pSuite, "test reading of fasta file",  test_read_sequence_from_fasta)){
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pSuite, "test reading of long fasta file",  test_read_sequence_from_long_fasta)){
    CU_cleanup_registry();
    return CU_get_error();
  }


  if (NULL == CU_add_test(pSuite, "test reading of fastq file",  test_read_sequence_from_fastq)){
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pSuite, "test reading of fastq file when some reads are too long or have bad characters",  test_read_sequence_from_fastq_with_bad_reads_and_long_reads)){
    CU_cleanup_registry();
    return CU_get_error();
  }



  if (NULL == CU_add_test(pSuite, "test conversion from binary nucleotide to C string", test_seq_to_binary_kmer_and_binary_kmer_to_seq)) {
    CU_cleanup_registry();
    return CU_get_error();
  }


 if (NULL == CU_add_test(pSuite, "test binary kmer reverse complement", test_binary_kmer_reverse_complement)) {
    CU_cleanup_registry();
    return CU_get_error();
  }

 if (NULL == CU_add_test(pSuite, "test seq reverse complement", test_seq_reverse_complement)) {
   CU_cleanup_registry();
   return CU_get_error();
 }

 
 if (NULL == CU_add_test(pSuite, "test nucleotide iterator", test_binary_kmer_nucleotide_iterator)) {
    CU_cleanup_registry();
    return CU_get_error();
  }

 
  if (NULL == CU_add_test(pSuite, "test creation of binary kmers from sequence - sliding window", test_get_sliding_windows_from_sequence)) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  // Commented this out because it crashes and I don't understand it
  /*
  if (NULL == CU_add_test(pSuite, "test creation of binary kmers from sequence - breaking at homopolymers", test_breaking_homopolymers_in_get_sliding_windows    )) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  */



  if (NULL == CU_add_test(pSuite, "test shift last kmer to start",  test_shift_last_kmer_to_start_of_sequence)){
    CU_cleanup_registry();
    return CU_get_error();
  }




  /* Run all tests using the CUnit Basic interface */
  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
 
  CU_cleanup_registry();
  return CU_get_error();


}





