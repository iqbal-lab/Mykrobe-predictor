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
  test_pop_load_and_print.c
*/

#include <stdlib.h>

#include <CUnit.h>
#include <Basic.h>

// cortex_var headers
#include "file_reader.h"
#include "dB_graph_population.h"
#include "element.h"
#include "open_hash/hash_table.h"
#include "supernode_cmp.h"

void test_load_two_people_in_same_populations_and_print_separately_their_supernodes()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }
  if(NUMBER_OF_COLOURS < 2)
  {
    warn("Test needs >=2 colours\n");
    return;
  }

  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size = 4;
  int max_retries = 10;

  dBGraph * hash_table = hash_table_new(number_of_bits, bucket_size,
                                        max_retries, kmer_size);

  if(hash_table == NULL)
  {
    die("Unable to alloc the hash table. dead before we even started.");
  }

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/test_pop_load_and_print/two_individuals_simple.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);


  //printf("Number of bases loaded is %d",seq_loaded);
  
  CU_ASSERT(seq_loaded == 44);
  CU_ASSERT(seq_read == 44);
  CU_ASSERT(bad_reads==0);

  char** array_of_supernodes_for_person1= (char**) calloc(10,sizeof(char*));
  char** array_of_supernodes_for_person2= (char**) calloc(10,sizeof(char*));

  if ( (array_of_supernodes_for_person1==NULL) || (array_of_supernodes_for_person2==NULL))
    {
      die("Cant start - OOM\n");
    }

  //these counters are used to make sure none of the files of printed out supernodes get too big. In fact
  //they are completely unused in test code, as we don't print anything
  //long supernode_count_person1=0;
  //long supernode_count_person2=0;
  
  //this on the other hand is used in testing.
  //int number_of_supernodes_in_person_1=0;

  //print_supernode will, in debug mode, alloc memory for you in your array, and put the supernode in it


  int i;


  warn("Upgrade this test to use the NEW PRINT FUNCTION and then uncomment it.\n");
  /*

  //  db_graph_traverse_specific_person_or_pop_for_supernode_printing(&db_graph_choose_output_filename_and_print_supernode_for_specific_person_or_pop, hash_table, &supernode_count_person1, 0, 
  //								  true, array_of_supernodes_for_person1,&number_of_supernodes_in_person_1);

  //printf("PERSON 1 has %d supernodes\n", number_of_supernodes_in_person_1);
  db_graph_set_all_visited_nodes_to_status_none(hash_table);

  int number_of_supernodes_in_person_2=0;
  db_graph_traverse_specific_person_or_pop_for_supernode_printing(&db_graph_choose_output_filename_and_print_supernode_for_specific_person_or_pop, hash_table, &supernode_count_person2, 1, 
  					   true, array_of_supernodes_for_person2,&number_of_supernodes_in_person_2);
  //printf("PERSON 2 has %d supernodes\n", number_of_supernodes_in_person_2);

  CU_ASSERT(number_of_supernodes_in_person_1==2);
  CU_ASSERT(number_of_supernodes_in_person_2==2);


  //for (i=0; i<number_of_supernodes_in_person_1; i++)
  // {
  //   printf("\nPerson 1 node %d : %s\n",i,array_of_supernodes_for_person1[i]);
  // }
  //for (i=0; i<number_of_supernodes_in_person_2; i++)
  // {
  //   printf("\nPerson 2 node %d : %s\n",i,array_of_supernodes_for_person2[i]);
  // }


  //Quicksort these results:
  qsort((void*) array_of_supernodes_for_person1, number_of_supernodes_in_person_1 , sizeof(char*),  &supernode_cmp);
  qsort((void*) array_of_supernodes_for_person2, number_of_supernodes_in_person_2 , sizeof(char*),  &supernode_cmp);
  // for (i=0; i<number_of_supernodes_in_person_1; i++)
  // {
  //   printf("\nAfter sort Person 1 node %d : %s\n",i,array_of_supernodes_for_person1[i]);
  // }
  //for (i=0; i<number_of_supernodes_in_person_2; i++)
  // {
  //   printf("\nafter sort Person 2 node %d : %s\n",i,array_of_supernodes_for_person2[i]);
  // }


  //Expected results are

  char* correct_answer_person_1[] ={"AAAA", "AACGTT"};
  char* correct_answer_person_2[] ={"CCCC", "CGTCAA"};


  int j;

  CU_ASSERT(2==number_of_supernodes_in_person_1);
  for (j=0; j<number_of_supernodes_in_person_1; j++)
    {
  //   printf("\nj is %d, Person 1 should have %s and we see %s\n", j,correct_answer_person_1[j], array_of_supernodes_for_person1[j]);
      CU_ASSERT_STRING_EQUAL(correct_answer_person_1[j], array_of_supernodes_for_person1[j]);
    }


  
  CU_ASSERT(2==number_of_supernodes_in_person_2);
  for (j=0; j<number_of_supernodes_in_person_2; j++)
    {
      // printf("Person 2 should have %s and we see%s\n", correct_answer_person_2[j], array_of_supernodes_for_person2[j]);

     CU_ASSERT_STRING_EQUAL(correct_answer_person_2[j], array_of_supernodes_for_person2[j]);
    }




  */


  //cleanup

  for (i=0; i<2; i++)
    {
      free(array_of_supernodes_for_person1[i]);
    }
 for (i=0; i<2; i++)
    {
      free(array_of_supernodes_for_person2[i]);
    }

 free(array_of_supernodes_for_person1);
 free(array_of_supernodes_for_person2);


  hash_table_free(&hash_table);
}


// Three people, with slight variants at one of two loci
// This is a good test set for seeing if we find the right shared variants
// However for this test case, just check that we get the right supernodes for each person,
// and that we can find the subsets of the supernodes that they ALL share

void test_take_three_people_each_with_one_read_and_find_variants()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // FIX TO USE NEW PRINT

  /*


  //first set up the hash/graph
  int kmer_size = 5;
  int number_of_bits = 5;
  int bucket_size    = 4;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph * hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);


  if (hash_table==NULL)
    {
      die("Unable to alloc the hash table. dead before we even started.");
    }


  long long seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/test_pop_load_and_print/three_indiv_simple/three_individuals_simple.colours", &bad_reads, hash_table, NULL);

  //printf("take 3 people test Number of bases loaded is %d",seq_loaded);
  CU_ASSERT(seq_loaded == 55);
  CU_ASSERT(bad_reads==0);

  char** array_of_supernodes_for_person1= (char**) calloc(20,sizeof(char*));
  char** array_of_supernodes_for_person2= (char**) calloc(20,sizeof(char*));
  char** array_of_supernodes_for_person3= (char**) calloc(20,sizeof(char*));
  
  if ( (array_of_supernodes_for_person1==NULL) || (array_of_supernodes_for_person2==NULL) || (array_of_supernodes_for_person3==NULL) )
    {
      die("Cant start - OOM\n");
    }

  
  //these counters are used to make sure none of the files of printed out supernodes get too big. In fact
  //they are completely unused in test code, as we don't print anything
  long supernode_count_person1=0;
  long supernode_count_person2=0;
  long supernode_count_person3=0;


  int number_of_supernodes_in_person_1=0;
  //print_supernode will, in debug mode, alloc memory for you in your array, and put the supernode in it
  db_graph_traverse_specific_person_or_pop_for_supernode_printing(&db_graph_choose_output_filename_and_print_supernode_for_specific_person_or_pop, hash_table, &supernode_count_person1, 0, 
					   true, array_of_supernodes_for_person1,&number_of_supernodes_in_person_1);
  //printf("PERSON 1 has %d supernodes\n", number_of_supernodes_in_person_1);
  db_graph_set_all_visited_nodes_to_status_none(hash_table);

  int i;
  //    for (i=0; i<number_of_supernodes_in_person_1; i++)
  //{
  //  printf("SUPERNODE %s\n", array_of_supernodes_for_person1[i]);
  // }

  
  int number_of_supernodes_in_person_2=0;
  db_graph_traverse_specific_person_or_pop_for_supernode_printing(&db_graph_choose_output_filename_and_print_supernode_for_specific_person_or_pop, hash_table,&supernode_count_person2, 1, 
					   true, array_of_supernodes_for_person2,&number_of_supernodes_in_person_2);
  // printf("PERSON 2 has %d supernodes\n", number_of_supernodes_in_person_2);
  db_graph_set_all_visited_nodes_to_status_none(hash_table);

  // for (i=0; i<number_of_supernodes_in_person_2; i++)
  //{
  //  printf("SUPERNODE %s\n", array_of_supernodes_for_person2[i]);
  //  }


  int number_of_supernodes_in_person_3=0;
  db_graph_traverse_specific_person_or_pop_for_supernode_printing(&db_graph_choose_output_filename_and_print_supernode_for_specific_person_or_pop, hash_table, &supernode_count_person3, 2, 
					   true, array_of_supernodes_for_person3,&number_of_supernodes_in_person_3);
  printf("PERSON 3 has %d supernodes\n", number_of_supernodes_in_person_3);
  db_graph_set_all_visited_nodes_to_status_none(hash_table);

  //for (i=0; i<number_of_supernodes_in_person_3; i++)
  // {
  //  printf("SUPERNODE %s\n", array_of_supernodes_for_person3[i]);
  // }


  //qsort results
  qsort((void*) array_of_supernodes_for_person1, number_of_supernodes_in_person_1 , sizeof(char*),  &supernode_cmp);
  qsort((void*) array_of_supernodes_for_person2, number_of_supernodes_in_person_2 , sizeof(char*),  &supernode_cmp);
  qsort((void*) array_of_supernodes_for_person3, number_of_supernodes_in_person_3 , sizeof(char*),  &supernode_cmp);


  //Expected results are

  char* correct_answer_person_1[] ={"AAGCCTCGACAGCCATGC"};
  char* correct_answer_person_2[]={"AAGCCTCGTTCGGCCATGC"};
  char* correct_answer_person_3[]={"AAGCCTCGCTAG","GCATGGCTAG"};
  char* rev_correct_answer_person_1[]={"GCATGGCTGTCGAGGCTT"};
  char* rev_correct_answer_person_2[]={"GCATGGCCGAACGAGGCTT"};
  char* rev_correct_answer_person_3[]={"CTAGCGAGGCTT","CTAGCCATGC"};


  for (i=0; i<number_of_supernodes_in_person_1; i++)
    {
      //  printf("\ni is %d, person 1, compare %s and %s\n", i, correct_answer_person_1[i], array_of_supernodes_for_person1[i]);
      CU_ASSERT(!strcmp(correct_answer_person_1[i], array_of_supernodes_for_person1[i]) || !strcmp(rev_correct_answer_person_1[i], array_of_supernodes_for_person1[i]) );
    }
  for (i=0; i<number_of_supernodes_in_person_2; i++)
    {
      // printf("\n i is %d, person 2, compare %s and %s\n", i, correct_answer_person_2[i], array_of_supernodes_for_person2[i]);
      CU_ASSERT(!strcmp(correct_answer_person_2[i], array_of_supernodes_for_person2[i]) || !strcmp(rev_correct_answer_person_2[i], array_of_supernodes_for_person2[i]) );
    }
  for (i=0; i<number_of_supernodes_in_person_3; i++)
    {
      //printf("\n i is %d, person 3, compare %s and %s\n", i, correct_answer_person_3[i], array_of_supernodes_for_person3[i]);
      CU_ASSERT(!strcmp(correct_answer_person_3[i], array_of_supernodes_for_person3[i]) || !strcmp(rev_correct_answer_person_3[i], array_of_supernodes_for_person3[i]) );
    }



  //cleanup


  for (i=0; i<number_of_supernodes_in_person_1; i++)
    {
      free(array_of_supernodes_for_person1[i]);
    }
  for (i=0; i<number_of_supernodes_in_person_2; i++)
    {
      free(array_of_supernodes_for_person2[i]);
    }
  for (i=0; i<number_of_supernodes_in_person_3; i++)
    {
      free(array_of_supernodes_for_person3[i]);
    }

  free(array_of_supernodes_for_person1);
  free(array_of_supernodes_for_person2);
  free(array_of_supernodes_for_person3);


  hash_table_free(&hash_table);
  */

}



void test_take_two_people_sharing_an_alu_and_find_supernodes()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // FIX to use new print

  /*

  //first set up the hash/graph
  int kmer_size = 31;
  int number_of_bits = 20;
  int bucket_size    = 10;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph * hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);



  if (hash_table==NULL)
    {
      die("Unable to alloc the hash table. dead before we even started.");
    }


  long long seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/test_pop_load_and_print/two_people_sharing_alu/two_people.colours",  &bad_reads,hash_table, NULL);
  //printf("Number of bases loaded is %d",seq_loaded);
  CU_ASSERT(seq_loaded == 677);
  CU_ASSERT(bad_reads ==0);

  char** array_of_supernodes_for_person1= (char**) calloc(316,sizeof(char*));
  char** array_of_supernodes_for_person2= (char**) calloc(255,sizeof(char*));
  
  if ( (array_of_supernodes_for_person1==NULL) || (array_of_supernodes_for_person2==NULL) )
    {
      die("Cant start - OOM\n");
    }

  
  //dummy counter used for making sure don't print too many supernodes to a file
  long supernode_count_person1=0;
  long supernode_count_person2=0;
  
  int  number_of_supernodes_in_person_1=0;
  //print_supernode will, in debug mode, alloc memory for you in your array, and put the supernode in it
  db_graph_traverse_specific_person_or_pop_for_supernode_printing(&db_graph_choose_output_filename_and_print_supernode_for_specific_person_or_pop, hash_table, &supernode_count_person1, 0, 
					   true, array_of_supernodes_for_person1,&number_of_supernodes_in_person_1);
  printf("PERSON 1 has %d supernodes\n", number_of_supernodes_in_person_1);
  db_graph_set_all_visited_nodes_to_status_none(hash_table);

  int i;
  for (i=0; i<number_of_supernodes_in_person_1; i++)
    {
      printf("SUPERNODE %s\n", array_of_supernodes_for_person1[i]);
    }

  
  int number_of_supernodes_in_person_2=0;
  db_graph_traverse_specific_person_or_pop_for_supernode_printing(&db_graph_choose_output_filename_and_print_supernode_for_specific_person_or_pop, hash_table, &supernode_count_person2, 1, 
					   true, array_of_supernodes_for_person2,&number_of_supernodes_in_person_2);
 printf("PERSON 2 has %d supernodes\n", number_of_supernodes_in_person_2);
 db_graph_set_all_visited_nodes_to_status_none(hash_table);
 
 for (i=0; i<number_of_supernodes_in_person_2; i++)
   {
     printf("SUPERNODE %s\n", array_of_supernodes_for_person2[i]);
   }


   printf("\n******   TODO ********Check these supernodes are correct\n");

   

   //Now see which bits of the two supernodes the two individuals have in common



  //cleanup


  for (i=0; i<number_of_supernodes_in_person_1; i++)
    {
      free(array_of_supernodes_for_person1[i]);
    }
  for (i=0; i<number_of_supernodes_in_person_2; i++)
    {
      free(array_of_supernodes_for_person2[i]);
    }


  free(array_of_supernodes_for_person1);
  free(array_of_supernodes_for_person2);

  hash_table_free(&hash_table);

  */

}





















