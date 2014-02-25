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
  test_binary_kmer.c
*/

// system libraries
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

// third party libraries
#include <CUnit.h>
#include <Basic.h>

// cortex_var headers
#include "binary_kmer.h"
#include "test_binary_kmer.h"


void test_that_bitfield_really_is_64bits()
{
  CU_ASSERT(sizeof(bitfield_of_64bits)==8);
}


void test_binary_kmer_assignment_operator()
{
  int i;

  bitfield_of_64bits  array[10] = {(bitfield_of_64bits) 1,(bitfield_of_64bits)2,(bitfield_of_64bits)7,(bitfield_of_64bits)25,
				   (bitfield_of_64bits)25656464646,(bitfield_of_64bits)43,(bitfield_of_64bits)4444444,(bitfield_of_64bits)101010,(bitfield_of_64bits)999,
				   (bitfield_of_64bits)32};

  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>10)
    {
      die("Had assumed you would never try going higher than 10 long longs in a BinaryKmer. Fix the test - very simple fix");
    }

  //set up your test kmer
  BinaryKmer test_kmer;
  binary_kmer_initialise_to_zero(&test_kmer);

  for (i=0 ; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
    {
      test_kmer[i]=array[i];

    }
  
  BinaryKmer assignee;
  binary_kmer_assignment_operator(assignee, test_kmer);

  for (i=0 ; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
    {
      CU_ASSERT(assignee[i]==test_kmer[i]);
    }



}

void test_binary_kmer_comparison_operator()
{

  BinaryKmer bk1;
  BinaryKmer bk2;
  BinaryKmer bk3;
  BinaryKmer bk4;
  BinaryKmer bk5;
  BinaryKmer bk6;
  BinaryKmer bk7;
  BinaryKmer bk8;

  int i;
  for (i=0; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
    {
      bk1[i]=0;
      bk2[i]=1;
      bk3[i]=0;
      bk4[i]=i+2;
      bk5[i]=i+1;
      bk6[i]=~0;
      bk7[i]=1<<i;
      bk8[i]=1<<i;
    }
  CU_ASSERT(binary_kmer_comparison_operator(bk1,bk2)==false);
  CU_ASSERT(binary_kmer_comparison_operator(bk1,bk3)==true);
  CU_ASSERT(binary_kmer_comparison_operator(bk3,bk4)==false);
  CU_ASSERT(binary_kmer_comparison_operator(bk4,bk5)==false);
  CU_ASSERT(binary_kmer_comparison_operator(bk2,bk6)==false);
  CU_ASSERT(binary_kmer_comparison_operator(bk7,bk8)==true);

}

void test_binary_kmer_less_than()
{
  
  BinaryKmer bk1;
  BinaryKmer bk2;

  int i;
  
  bitfield_of_64bits j;
  short kmer_size;
  
  for (kmer_size=31; kmer_size < NUMBER_OF_BITFIELDS_IN_BINARY_KMER*32; kmer_size=kmer_size+32)
    {
      int number_of_bitfields_fully_used = kmer_size/32;
      int number_of_bits_in_most_sig_bitfield = 2* (kmer_size-(32*number_of_bitfields_fully_used));

      for (j=0; j<10; j++)
	{
	  
	  for (i=NUMBER_OF_BITFIELDS_IN_BINARY_KMER-number_of_bitfields_fully_used-1; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
	    {
	      bk1[i]=j;
	      bk2[i]=j;
	    }
	  
	  //Now make a single difference between the two
	  for (i=NUMBER_OF_BITFIELDS_IN_BINARY_KMER-number_of_bitfields_fully_used-1; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
	    {
	      bk1[i]=(bitfield_of_64bits) 10;
	      bk2[i]=(bitfield_of_64bits) 11;
	      CU_ASSERT(binary_kmer_less_than(bk1,bk2, kmer_size)==true);
	
	if (binary_kmer_less_than(bk1,bk2, kmer_size)!=true)
	  {
      die("Breaks at kmer = %d, i=%d, j=%lu\n"
          "num of bitfields fully used is %d, and number of bits in most sig is %d\n",
          kmer_size, i, (unsigned long)j,
          number_of_bitfields_fully_used, number_of_bits_in_most_sig_bitfield);
	  }
	
	CU_ASSERT(binary_kmer_less_than(bk2,bk1, kmer_size)==false);
	
	//reset values, so now both arrays are identical
	bk1[i]=j;
	bk2[i]=j;
	
	    }
	  
	}
    }
  
}

void test_binary_kmer_right_shift_one_base()
{
  BinaryKmer test;
  BinaryKmer shifted;

  // Start with largest possible kmer size, then work down
  int kmer_size = NUMBER_OF_BITFIELDS_IN_BINARY_KMER*32-1;
  
  unsigned long bases[4] = {0,1,2,3};
  //unsigned long bases[4] = {2,2,2,2};
  
  for(; kmer_size > 0; kmer_size -= 2)
  {
    // zero
    int i;
    for(i = 0; i < NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
    {
      test[i] = 0;
      shifted[i] = 0;
    }
  
    // Fill test one base at a time
    int bases_in_most_sig_word = kmer_size % 32;

    if(bases_in_most_sig_word == 0)
    {
      bases_in_most_sig_word = 32;
    }

    int base_offset = bases_in_most_sig_word - 1;
    int word = NUMBER_OF_BITFIELDS_IN_BINARY_KMER - (kmer_size+31)/32;

    test[word] |= (bases[0] << (base_offset*2));
    base_offset--;

    for(i = 1; i < kmer_size; i++)
    {
      // Update word and offset
      if(base_offset < 0)
      {
        base_offset = 31;
        word++;
      }

      test[word] |= (bases[i % 4] << (base_offset*2));
      shifted[word] |= (bases[(i+3) % 4] << (base_offset*2));
      base_offset--;
    }

    //printf("setup\n");
    //print_kmer(test);
    //print_kmer(shifted);

    // shift and compare
    binary_kmer_right_shift_one_base(test);

    //print_kmer(test);
    //print_kmer(shifted);

    CU_ASSERT(binary_kmer_comparison_operator(test, shifted));

    // Shift off the end, compare with zero
    for(i = 0; i < kmer_size-1; i++)
    {
      binary_kmer_right_shift_one_base(test);
    }
    
    //printf("zero:\n");
    //print_kmer(test);

    for(i = 0; i < NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
    {
      CU_ASSERT(test[i] == 0);
    }
  }
}

void test_binary_kmer_left_shift_one_base()
{
  BinaryKmer test;
  BinaryKmer shifted;

  // Start with largest possible kmer size, then work down
  int kmer_size = NUMBER_OF_BITFIELDS_IN_BINARY_KMER*32-1;
  
  unsigned long bases[4] = {0,1,2,3};
  //unsigned long bases[4] = {2,2,2,2};
  
  for(; kmer_size > 0; kmer_size -= 2)
  {
    // zero
    int i;
    for(i = 0; i < NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
    {
      test[i] = 0;
      shifted[i] = 0;
    }
  
    // Fill test one base at a time
    int bases_in_most_sig_word = kmer_size % 32;

    if(bases_in_most_sig_word == 0)
    {
      bases_in_most_sig_word = 32;
    }

    int base_offset = bases_in_most_sig_word - 1;
    int word = NUMBER_OF_BITFIELDS_IN_BINARY_KMER - (kmer_size+31)/32;

    test[word] |= (bases[0] << (base_offset*2));
    base_offset--;

    for(i = 1; i < kmer_size; i++)
    {
      // Update word and offset
      if(base_offset < 0)
      {
        base_offset = 31;
        word++;
      }

      test[word] |= (bases[i % 4] << (base_offset*2));
      shifted[word] |= (bases[(i+3) % 4] << (base_offset*2));
      base_offset--;
    }

    //printf("setup\n");
    //print_kmer(test);
    //print_kmer(shifted);

    // mask last base
    test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] &= ~3;
    // shift left
    binary_kmer_left_shift_one_base(shifted, kmer_size);

    //printf("compare\n");
    //print_kmer(test);
    //print_kmer(shifted);

    CU_ASSERT(binary_kmer_comparison_operator(test, shifted));

    // Shift off the end, compare with zero
    for(i = 0; i < kmer_size-1; i++)
    {
      binary_kmer_left_shift_one_base(test, kmer_size);
    }
    
    //printf("zero:\n");
    //print_kmer(test);

    for(i = 0; i < NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
    {
      CU_ASSERT(test[i] == 0);
    }
  }
}


void test_seq_to_binary_kmer_and_binary_kmer_to_seq(){
  
  BinaryKmer kmer;
  binary_kmer_initialise_to_zero(&kmer);
  char kmer_seq[1000];
  
  seq_to_binary_kmer("ATCGCGC",7, &kmer);
 
  CU_ASSERT_STRING_EQUAL("ATCGCGC",binary_kmer_to_seq(&kmer,7,kmer_seq));

  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>1)
    {
      binary_kmer_initialise_to_zero(&kmer);
      seq_to_binary_kmer("CATCAGTGGGACATAAACCACACAGATGACACACACA",37, &kmer);
      CU_ASSERT_STRING_EQUAL("CATCAGTGGGACATAAACCACACAGATGACACACACA",binary_kmer_to_seq(&kmer,37,kmer_seq));


      binary_kmer_initialise_to_zero(&kmer);
      seq_to_binary_kmer("CATCAGTGGGACATAAACCACACAGATGACACACACACATCAGTGGGACATAAACCACACAGA",63, &kmer);
      CU_ASSERT_STRING_EQUAL("CATCAGTGGGACATAAACCACACAGATGACACACACACATCAGTGGGACATAAACCACACAGA",binary_kmer_to_seq(&kmer,63,kmer_seq));
    }

  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>2)
    {
      binary_kmer_initialise_to_zero(&kmer);
      seq_to_binary_kmer("CATCAGGGGTCAGTCAGTCACACAAAAACGTCAGTCGTGTCATCAGGGGTCAGTCAGTCACACAAAAACGTCAGTCGTGT",80, &kmer);
      CU_ASSERT_STRING_EQUAL("CATCAGGGGTCAGTCAGTCACACAAAAACGTCAGTCGTGTCATCAGGGGTCAGTCAGTCACACAAAAACGTCAGTCGTGT",binary_kmer_to_seq(&kmer,80,kmer_seq));

      binary_kmer_initialise_to_zero(&kmer);
      seq_to_binary_kmer("CATCAGGGGTCAGTCAGTCACACAAAAACGTCAGTCGTGTCATCAGGGGTCAGTCAGTCACACAAAAACGTCAGTCGTGTGGGGGGGGTTTTCAC",95, &kmer);
      CU_ASSERT_STRING_EQUAL("CATCAGGGGTCAGTCAGTCACACAAAAACGTCAGTCGTGTCATCAGGGGTCAGTCAGTCACACAAAAACGTCAGTCGTGTGGGGGGGGTTTTCAC",binary_kmer_to_seq(&kmer,95,kmer_seq));

    }



  
}

void test_binary_kmer_reverse_complement()
{
  BinaryKmer kmer, kmer_reverse;
  char seq[1000];

  // test with various different kmer sizes
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER == 1)
  {
    seq_to_binary_kmer("ATCGCGC", 7, &kmer);
    binary_kmer_reverse_complement(&kmer, 7, &kmer_reverse);
    CU_ASSERT_STRING_EQUAL("GCGCGAT",binary_kmer_to_seq(&kmer_reverse, 7, seq));

    seq_to_binary_kmer("A", 1, &kmer);
    binary_kmer_reverse_complement(&kmer, 1, &kmer_reverse);
    CU_ASSERT_STRING_EQUAL("T", binary_kmer_to_seq(&kmer_reverse, 1, seq));

    seq_to_binary_kmer("GGCCCCGCCCCGCCCCGCCCCGCCCCGCCCC",31, &kmer);
    binary_kmer_reverse_complement(&kmer, 31, &kmer_reverse);
    CU_ASSERT_STRING_EQUAL("GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCC",
                           binary_kmer_to_seq(&kmer_reverse, 31, seq));
  }
  else if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER == 2)
  {
    seq_to_binary_kmer("AAACGTAACGTAACGTAACGTTTTTTCATGGCA", 33, &kmer);
    binary_kmer_reverse_complement(&kmer, 33, &kmer_reverse);
    CU_ASSERT_STRING_EQUAL("TGCCATGAAAAAACGTTACGTTACGTTACGTTT",
                           binary_kmer_to_seq(&kmer_reverse, 33, seq));

    seq_to_binary_kmer("AAACGTAACGTAACGTAACGTTTTTTCATGGCAACGT", 37, &kmer);
    binary_kmer_reverse_complement(&kmer, 37, &kmer_reverse);
    CU_ASSERT_STRING_EQUAL("ACGTTGCCATGAAAAAACGTTACGTTACGTTACGTTT",
                           binary_kmer_to_seq(&kmer_reverse, 37, seq));
  
    seq_to_binary_kmer("ACGTTGCCATGAAAAAACGTTACGTTACGTTACGTTTAAAAAAAAAACCCCCCCCCC",
                       57, &kmer);
    binary_kmer_reverse_complement(&kmer, 57, &kmer_reverse);
    CU_ASSERT_STRING_EQUAL("GGGGGGGGGGTTTTTTTTTTAAACGTAACGTAACGTAACGTTTTTTCATGG"
                           "CAACGT",binary_kmer_to_seq(&kmer_reverse, 57, seq));

    seq_to_binary_kmer("ACGTTGCCATGAAAAAACGTTACGTTACGTTACGTTTAAAAAAAAAACCCCCCCC"
                       "CCGGGTAC", 63, &kmer);
    binary_kmer_reverse_complement(&kmer, 63, &kmer_reverse);
    CU_ASSERT_STRING_EQUAL("GTACCCGGGGGGGGGGTTTTTTTTTTAAACGTAACGTAACGTAACGTTTTT"
                           "TCATGGCAACGT",
                           binary_kmer_to_seq(&kmer_reverse, 63, seq));
  }
  else if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER == 3)
  {
    seq_to_binary_kmer("TACGTTGCCATGAAAAAACGTTACGTTACGTTACGTTTAAAAAAAAAACCCCCCC"
                       "CCCGGGTACC", 65, &kmer);
    binary_kmer_reverse_complement(&kmer, 65, &kmer_reverse);
    CU_ASSERT_STRING_EQUAL("GGTACCCGGGGGGGGGGTTTTTTTTTTAAACGTAACGTAACGTAACGTTTT"
                           "TTCATGGCAACGTA",
                           binary_kmer_to_seq(&kmer_reverse,65,seq));
      
    seq_to_binary_kmer("AACGTGTGTGCCCATACAGGAACGTGTGTGCCCATACAGGAACGTGTGTGCCCAT"
                       "ACAGGAACGTGTGTGCCCATACAGGG", 81, &kmer);
    binary_kmer_reverse_complement(&kmer, 81, &kmer_reverse);
    CU_ASSERT_STRING_EQUAL("CCCTGTATGGGCACACACGTTCCTGTATGGGCACACACGTTCCTGTATGGG"
                           "CACACACGTTCCTGTATGGGCACACACGTT",
                           binary_kmer_to_seq(&kmer_reverse, 81, seq));
  }
  else if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER == 4)
  {
    seq_to_binary_kmer("CCCCCCCCCTATGGGCACATATGGGCACATATGGGCACATATGGGCACATATG"
                       "GGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCA", 97, &kmer);
    binary_kmer_reverse_complement(&kmer, 97, &kmer_reverse);
    CU_ASSERT_STRING_EQUAL("TGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATAT"
                           "GTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATAGGGGGGGGG",
                           binary_kmer_to_seq(&kmer_reverse, 97, seq));

    seq_to_binary_kmer("CCCCCCCCTATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGG"
                       "GCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCA"
                       "CATATGGGCACATATGGGCAC", 127, &kmer);
    binary_kmer_reverse_complement(&kmer, 127, &kmer_reverse);
    CU_ASSERT_STRING_EQUAL("GTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATA"
                           "TGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCAT"
                           "ATGTGCCCATATGTGCCCATAGGGGGGGG",
                           binary_kmer_to_seq(&kmer_reverse, 127, seq));
  }
  else if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER == 5)
  {
    seq_to_binary_kmer("TGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGC"
                       "CCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATA"
                       "TGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATAGGGGGGGGG", 159,
                       &kmer);
    binary_kmer_reverse_complement(&kmer, 159, &kmer_reverse);
    CU_ASSERT_STRING_EQUAL("CCCCCCCCCTATGGGCACATATGGGCACATATGGGCACATATGGGCACATA"
                           "TGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATAT"
                           "GGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATG"
                           "GGCACA",  
                           binary_kmer_to_seq(&kmer_reverse, 159, seq));
  }
  else if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER == 6)
  {
    seq_to_binary_kmer("GTCGCCATCATCCAGGTCGCCGAACTACGGTGGTAAAGCTGGAGAGGCCAAGTTG"
                       "CTGCCGATGCTCACCTATAGGACCCTGGTGTATACACATGCATGCTTGACAGCAT"
                       "GTCCGCTCTGTGGCGCGGCTATTTCACGCTGCTCCTATGCAACCGTCCACCGATA"
                       "TCCTCCCCCAAGCATCACATTC", 187,
                       &kmer);
    binary_kmer_reverse_complement(&kmer, 187, &kmer_reverse);
    CU_ASSERT_STRING_EQUAL("GAATGTGATGCTTGGGGGAGGATATCGGTGGACGGTTGCATAGGAGCAGCG"
                           "TGAAATAGCCGCGCCACAGAGCGGACATGCTGTCAAGCATGCATGTGTATA"
                           "CACCAGGGTCCTATAGGTGAGCATCGGCAGCAACTTGGCCTCTCCAGCTTT"
                           "ACCACCGTAGTTCGGCGACCTGGATGATGGCGAC",  
                           binary_kmer_to_seq(&kmer_reverse, 187, seq));
  }
  else if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER == 7)
  {
    seq_to_binary_kmer("GTTTCTTTGTCGTGGTTGGGGTTCTTAGAATGTGCCTAGCATATTGATCCCTCGG"
                       "CTGGAACACAACAGCGTCAGGGTCCGGTATTCGTTCAGAATGGTGTACCCGATGG"
                       "GCGCATAGAAAGTTTCTGACCTAGGTGCTCCGCACGATCCGTATGTGCATTTGCA"
                       "CGTCTTGGACGGGGGAGTACCGCACCAAAAGGATCCCAACTCCTTAACTTGTAAA"
                       "TCC", 223,
                       &kmer);
    binary_kmer_reverse_complement(&kmer, 223, &kmer_reverse);
    CU_ASSERT_STRING_EQUAL("GGATTTACAAGTTAAGGAGTTGGGATCCTTTTGGTGCGGTACTCCCCCGTC"
                           "CAAGACGTGCAAATGCACATACGGATCGTGCGGAGCACCTAGGTCAGAAAC"
                           "TTTCTATGCGCCCATCGGGTACACCATTCTGAACGAATACCGGACCCTGAC"
                           "GCTGTTGTGTTCCAGCCGAGGGATCAATATGCTAGGCACATTCTAAGAACC"
                           "CCAACCACGACAAAGAAAC",  
                           binary_kmer_to_seq(&kmer_reverse, 223, seq));
  }
  else if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER == 8)
  {
    seq_to_binary_kmer("CCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATA"
                       "TGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGC"
                       "CCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATA"
                       "TGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGC"
                       "CCATATGTGCCCATATGTGCCCATATGTGCGGGGG", 255, &kmer);
    binary_kmer_reverse_complement(&kmer, 255, &kmer_reverse);
    CU_ASSERT_STRING_EQUAL("CCCCCGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACAT"
                           "ATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATA"
                           "TGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATAT"
                           "GGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATG"
                           "GGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGG",
                           binary_kmer_to_seq(&kmer_reverse, 255, seq));  
  }
}


void test_seq_reverse_complement()
{
  char out[100];
  char * seq = "AAAAAA";
  CU_ASSERT_STRING_EQUAL(seq_reverse_complement(seq, 6, out), "TTTTTT");
  
  seq = "ATAAAA";
  CU_ASSERT_STRING_EQUAL(seq_reverse_complement(seq, 6, out), "TTTTAT");
  
  seq = "CGATAAAA";
  CU_ASSERT_STRING_EQUAL(seq_reverse_complement(seq, 8, out), "TTTTATCG");
 
  seq = "CGATAAAAGG";
  CU_ASSERT_STRING_EQUAL(seq_reverse_complement(seq, 10, out), "CCTTTTATCG");
   
  seq = "";
  CU_ASSERT_STRING_EQUAL(seq_reverse_complement(seq, 0, out), "");
}



void test_binary_kmer_nucleotide_iterator(){
  
  int count = 0;
  
  void check_nucleotides(Nucleotide base){
    
    if (base == Adenine || base == Guanine || base == Cytosine || base == Thymine){
      count ++;
    }
    else{
      count --;
    }
  }
  
  nucleotide_iterator(check_nucleotides);
  
  CU_ASSERT_EQUAL(count,4);
}



void test_get_sliding_windows_from_sequence(){

    
   //----------------------------------
    KmerSlidingWindowSet * windows = malloc(sizeof(KmerSlidingWindowSet));  
    if (windows == NULL){
      die("Out of memory trying to allocate a KmerSlidingWindowSet");
    } 

    binary_kmer_alloc_kmers_set(windows,20,30);
    
    CU_ASSERT_EQUAL(windows->max_nwindows,20);

    //----------------------------------


    // ****************************************************
    // First tests do not involve any homopolymer breaking:
    // ****************************************************

  
    char * seq  = "AAAAANNTTTTGGGG";
    char qual0[15] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    char kmer_seq3[4];

    int nkmers1 = get_sliding_windows_from_sequence(seq,qual0,strlen(seq),0,3,windows,20,20, false, 0);

    CU_ASSERT_EQUAL(windows->nwindows,2);

    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 3);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),3,kmer_seq3),"AAA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[1]),3,kmer_seq3),"AAA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[2]),3,kmer_seq3),"AAA");

    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 6);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[0]),3,kmer_seq3),"TTT");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[1]),3,kmer_seq3),"TTT");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[2]),3,kmer_seq3),"TTG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[3]),3,kmer_seq3),"TGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[4]),3,kmer_seq3),"GGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[5]),3,kmer_seq3),"GGG");

    CU_ASSERT_EQUAL(nkmers1, 9);
    
    char qual1[15] = { 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 10, 20, 20, 20};

    int nkmers2 = get_sliding_windows_from_sequence(seq,qual1,strlen(seq),15,3,windows,20,20, false,0);

    CU_ASSERT_EQUAL(windows->nwindows,3);

    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 3);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),3,kmer_seq3),"AAA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[1]),3,kmer_seq3),"AAA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[2]),3,kmer_seq3),"AAA");

    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 2);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[0]),3,kmer_seq3),"TTT");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[1]),3,kmer_seq3),"TTT");

    CU_ASSERT_EQUAL((windows->window[2]).nkmers, 1);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[2].kmer[0]),3,kmer_seq3),"GGG");

    CU_ASSERT_EQUAL(nkmers2, 6);

    char qual2[15] = { 5, 10, 20, 20, 20, 20, 0, 0, 0, 20, 20, 30, 20, 20, 20};

    int nkmers3 = get_sliding_windows_from_sequence(seq,qual2,strlen(seq),15,3,windows,20,20, false,0);

    CU_ASSERT_EQUAL(windows->nwindows,2);

    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 1);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),3,kmer_seq3),"AAA");

    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 4);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[0]),3,kmer_seq3),"TTG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[1]),3,kmer_seq3),"TGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[2]),3,kmer_seq3),"GGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[3]),3,kmer_seq3),"GGG");
    
    CU_ASSERT_EQUAL(nkmers3, 5);

    char qual3[15] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int nkmers4 = get_sliding_windows_from_sequence(seq,qual3,strlen(seq),15,3,windows,20,20, false,0);

    CU_ASSERT_EQUAL(windows->nwindows,0);
    CU_ASSERT_EQUAL(nkmers4, 0);

    char kmer_seq5[6];
    int nkmers5 = get_sliding_windows_from_sequence(seq,qual0,strlen(seq),0,5,windows,20,20, false,0);
    
    CU_ASSERT_EQUAL(windows->nwindows,2);
  
    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 1);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),5,kmer_seq5),"AAAAA");

  
    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 4);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[0]),5,kmer_seq5),"TTTTG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[1]),5,kmer_seq5),"TTTGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[2]),5,kmer_seq5),"TTGGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[3]),5,kmer_seq5),"TGGGG");

    CU_ASSERT_EQUAL(nkmers5, 5);
    


    // Now try some examples with large kmers

    
    if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>2)
      {
	
	char* seq_contains_only_one_good_33mer="CCCCCCCCTANTGGGCACATATGGGCACATATGGNGCACATATGGGCACATATGGGCACATATGGNGCACATATGGNNNNNNNNNNGCACATATGGGNCANCATATGNNGCACATATGGGCACATATGGGCACATATGGGCANC";
        char qual_144_zeroes[144] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	char qual_144_30s[144] = { 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30};
	char qual_all_30s_except_for_the_one_33mer[144] = { 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 30, 30};


	char kmer_seq33[34];
	int quality_cutoff=0;
	short kmer_size=33;
	int nkmers_33 = get_sliding_windows_from_sequence(seq_contains_only_one_good_33mer,qual_144_zeroes,strlen(seq_contains_only_one_good_33mer),quality_cutoff,kmer_size,windows,40,40, false,0);

	CU_ASSERT_EQUAL(windows->nwindows,1);
	CU_ASSERT_EQUAL((windows->window[0]).nkmers, 1);
	CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),33,kmer_seq33),"GCACATATGGGCACATATGGGCACATATGGGCA");
	CU_ASSERT(nkmers_33==1);

	quality_cutoff=20;
	nkmers_33 = get_sliding_windows_from_sequence(seq_contains_only_one_good_33mer,qual_144_30s,strlen(seq_contains_only_one_good_33mer),quality_cutoff,kmer_size,windows,40,40, false,0);
	CU_ASSERT_EQUAL(windows->nwindows,1);
	CU_ASSERT_EQUAL((windows->window[0]).nkmers, 1);
	CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),33,kmer_seq33),"GCACATATGGGCACATATGGGCACATATGGGCA");
	CU_ASSERT(nkmers_33==1);

	nkmers_33 = get_sliding_windows_from_sequence(seq_contains_only_one_good_33mer,qual_all_30s_except_for_the_one_33mer,strlen(seq_contains_only_one_good_33mer),quality_cutoff,kmer_size,windows,40,40, false,0);
	CU_ASSERT_EQUAL(windows->nwindows,0);
	CU_ASSERT(nkmers_33==0);


      }



    binary_kmer_free_kmers_set(&windows);
    
    CU_ASSERT(windows == NULL);


    
}


void test_breaking_homopolymers_in_get_sliding_windows ()
{
   //----------------------------------
    KmerSlidingWindowSet * windows = malloc(sizeof(KmerSlidingWindowSet));  
    if (windows == NULL){
      die("Out of memory trying to allocate a KmerSlidingWindowSet");
    } 

    binary_kmer_alloc_kmers_set(windows,20,30);
    
    CU_ASSERT_EQUAL(windows->max_nwindows,20);

    //----------------------------------


    char * seq  = "AAAAANNTCAGAT";
    char qual0[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    char kmer_seq3[4];

    //break at homopolymers of length 4. ie you should never load any homopolymer of length >3, and when you reach one of length >=4, you don't start your next window
    //until after the end of the homopolymer
    int nkmers1 = get_sliding_windows_from_sequence(seq,qual0,strlen(seq),0,3,windows,20,20, true, 4);

    CU_ASSERT_EQUAL(windows->nwindows,2);

    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 1);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),3,kmer_seq3),"AAA");
 
    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 4);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[0]),3,kmer_seq3),"TCA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[1]),3,kmer_seq3),"CAG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[2]),3,kmer_seq3),"AGA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[3]),3,kmer_seq3),"GAT");

    CU_ASSERT_EQUAL(nkmers1, 5);

    //now try length 3, which is also the kmer_length. Should only get one window
    nkmers1 = get_sliding_windows_from_sequence(seq,qual0,strlen(seq),0,3,windows,20,20, true, 3);

    CU_ASSERT_EQUAL(windows->nwindows,1);
    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 4);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),3,kmer_seq3),"TCA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[1]),3,kmer_seq3),"CAG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[2]),3,kmer_seq3),"AGA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[3]),3,kmer_seq3),"GAT");

    CU_ASSERT_EQUAL(nkmers1, 4);

    //does it break if you enter hom cutoff of 0 or 1 or 2?

    //cutoff is 2
    nkmers1 = get_sliding_windows_from_sequence(seq,qual0,strlen(seq),0,3,windows,20,20, true, 2);

    CU_ASSERT_EQUAL(windows->nwindows,1);
    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 4);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),3,kmer_seq3),"TCA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[1]),3,kmer_seq3),"CAG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[2]),3,kmer_seq3),"AGA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[3]),3,kmer_seq3),"GAT");

    CU_ASSERT_EQUAL(nkmers1, 4);

    //cutoff is 1
    nkmers1 = get_sliding_windows_from_sequence(seq,qual0,strlen(seq),0,3,windows,20,20, true, 1);

    CU_ASSERT_EQUAL(windows->nwindows,0);
    CU_ASSERT_EQUAL(nkmers1, 0);

    //cutoff is zero
    nkmers1 = get_sliding_windows_from_sequence(seq,qual0,strlen(seq),0,3,windows,20,20, true, 0);

    CU_ASSERT_EQUAL(windows->nwindows,0);
    CU_ASSERT_EQUAL(nkmers1, 0);



    // another example:

    seq  = "AAAAANNTCTCCCCCTAGATGGGGGGGGGGGGGCCACCCNCCCCGTGATAT";
    char qual_other[51] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    char kmer_seq7[8];

    //first of all - not breaking homopolymers, only the N's cause problems:
    nkmers1 = get_sliding_windows_from_sequence(seq,qual_other,strlen(seq),0,7,windows,40,40, false, 0);

    CU_ASSERT_EQUAL(windows->nwindows,2);

    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 26);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),7,kmer_seq7),"TCTCCCC");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[1]),7,kmer_seq7),"CTCCCCC");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[2]),7,kmer_seq7),"TCCCCCT");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[3]),7,kmer_seq7),"CCCCCTA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[4]),7,kmer_seq7),"CCCCTAG");
    // ... not doing them all
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[25]),7,kmer_seq7),"GCCACCC");



    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 5);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[0]),7,kmer_seq3),"CCCCGTG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[1]),7,kmer_seq3),"CCCGTGA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[2]),7,kmer_seq3),"CCGTGAT");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[3]),7,kmer_seq3),"CGTGATA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[4]),7,kmer_seq3),"GTGATAT");

    CU_ASSERT_EQUAL(nkmers1, 31);


    //then, breaking homopolymers - no homopolymer of length 1 or greater
    nkmers1 = get_sliding_windows_from_sequence(seq,qual_other,strlen(seq),0,7,windows,40,40, true, 1);
    CU_ASSERT_EQUAL(windows->nwindows,0);
    CU_ASSERT_EQUAL(nkmers1, 0);
    //then, breaking homopolymers - no homopolymer of length 2 or greater
    nkmers1 = get_sliding_windows_from_sequence(seq,qual_other,strlen(seq),0,7,windows,40,40, true, 2);
    CU_ASSERT_EQUAL(windows->nwindows,1);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),7,kmer_seq7),"GTGATAT");
    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 1);
    CU_ASSERT_EQUAL(nkmers1, 1);
    //then, breaking homopolymers - no homopolymer of length 3 or greater
    nkmers1 = get_sliding_windows_from_sequence(seq,qual_other,strlen(seq),0,7,windows,40,40, true, 3);
    CU_ASSERT_EQUAL(windows->nwindows,2);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),7,kmer_seq7),"TAGATGG");
    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 1);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[0]),7,kmer_seq7),"GTGATAT");
    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 1);
    CU_ASSERT_EQUAL(nkmers1, 2);
    //then, breaking homopolymers - no homopolymer of length 4 or greater
    nkmers1 = get_sliding_windows_from_sequence(seq,qual_other,strlen(seq),0,7,windows,40,40, true, 4);
    CU_ASSERT_EQUAL(windows->nwindows,2);
    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 2);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),7,kmer_seq7),"TAGATGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[1]),7,kmer_seq7),"AGATGGG");

    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 1);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[0]),7,kmer_seq7),"GTGATAT");
    CU_ASSERT_EQUAL(nkmers1, 3);

    //"AAAAANNTCTCCCCCTAGATGGGGGGGGGGGGGCCACCCNCCCCGTGATAT"
    //then, breaking homopolymers - no homopolymer of length 5 or greater
    nkmers1 = get_sliding_windows_from_sequence(seq,qual_other,strlen(seq),0,7,windows,40,40, true, 5);
    CU_ASSERT_EQUAL(windows->nwindows,3);
    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 1);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),7,kmer_seq7),"TCTCCCC");

    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 3);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[0]),7,kmer_seq7),"TAGATGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[1]),7,kmer_seq7),"AGATGGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[2]),7,kmer_seq7),"GATGGGG");

    CU_ASSERT_EQUAL((windows->window[2]).nkmers, 5);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[2].kmer[0]),7,kmer_seq7),"CCCCGTG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[2].kmer[1]),7,kmer_seq7),"CCCGTGA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[2].kmer[2]),7,kmer_seq7),"CCGTGAT");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[2].kmer[3]),7,kmer_seq7),"CGTGATA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[2].kmer[4]),7,kmer_seq7),"GTGATAT");

    CU_ASSERT_EQUAL(nkmers1, 9);

    



    



    // Now try some examples with large kmers

    
    if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>2)
      {
	
	char* seq_contains_only_one_good_33mer="CCCCCCCCTANTGGGCACATATGGGCACATATGGNGCACATATGGGCACATATGGGCACATATGGNGCACATATGGNNNNNNNNNNGCACATATGGGNCANCATATGNNGCACATATGGGCACATATGGGCACATATGGGCANC";
        char qual_144_zeroes[144] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	char kmer_seq33[34];
	int quality_cutoff=0;
	short kmer_size=33;
	//let's put a homopolymer cutoff of 4 - should have no effect, as the one good 33mer contains none longer than 3
	int nkmers_33 = get_sliding_windows_from_sequence(seq_contains_only_one_good_33mer,qual_144_zeroes,strlen(seq_contains_only_one_good_33mer),quality_cutoff,kmer_size,windows,40,40, true,4);

	CU_ASSERT_EQUAL(windows->nwindows,1);
	CU_ASSERT_EQUAL((windows->window[0]).nkmers, 1);
	CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),33,kmer_seq33),"GCACATATGGGCACATATGGGCACATATGGGCA");
	CU_ASSERT(nkmers_33==1);


	//let's put a homopolymer cutoff of 3 - should find nothing
	nkmers_33 = get_sliding_windows_from_sequence(seq_contains_only_one_good_33mer,qual_144_zeroes,strlen(seq_contains_only_one_good_33mer),quality_cutoff,kmer_size,windows,40,40, true,3);

	CU_ASSERT_EQUAL(windows->nwindows,0);
	CU_ASSERT(nkmers_33==0);



      }



    binary_kmer_free_kmers_set(&windows);
    
    CU_ASSERT(windows == NULL);

  
}
