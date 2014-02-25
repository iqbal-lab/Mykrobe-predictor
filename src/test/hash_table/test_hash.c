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
  test_hash.c
*/

#include <stdlib.h>

#include <CUnit.h>
#include <Basic.h>

// cortex_var headers
#include "binary_kmer.h"
#include "element.h"
#include "open_hash/hash_table.h"

void test_hash_table_find_or_insert()
{

  short kmer_size;
  long long max_key_given_kmer_size;

      

  //adds binary kmers to a hash, then tries to find them
  //will go through all possible bin kmers given the kmer_size, stepping 
  // with granulatiry step. Since for kmer_size large (eg 21) it takes to long to try ALL
  // possible 21-mers, the last argument allows you to specify the largest one to test


  BinaryKmer tmp_kmer;
  BinaryKmer tmp_kmer2;
  
  void test(short kmer_size, int num_bits, int bucket, int max_tries, int step, long long max_kmer_to_test)
    {
      
      int number_of_bits      = num_bits;
      int bucket_size         = bucket;
      long long bad_reads     = 0; 
      int max_retries         = max_tries;
      boolean found           = false;
      
      HashTable* hash_table  = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
      
      
      long long i;
          

      for (i=0; i< max_kmer_to_test; i++)
	{

	  BinaryKmer b;
	  binary_kmer_initialise_to_zero(&b);
	  b[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]=(bitfield_of_64bits) i;

	  hash_table_find_or_insert(element_get_key(&b, hash_table->kmer_size, &tmp_kmer),&found, hash_table);
	  if (found==false)
	    {
	      CU_ASSERT(binary_kmer_comparison_operator(b,*(element_get_key(&b, hash_table->kmer_size, &tmp_kmer)) ) ); 
	    }
	}
      
      Element* e=NULL;
      
      
      for (i=0; i< max_kmer_to_test; i=i+step)
	{

	  BinaryKmer b;
	  binary_kmer_initialise_to_zero(&b);
	  b[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]=(bitfield_of_64bits) i;

	  e = hash_table_find(element_get_key(&b,hash_table->kmer_size,&tmp_kmer), hash_table);
	  CU_ASSERT(e!=NULL);
	  if (e !=NULL)
	    {
	      CU_ASSERT(binary_kmer_comparison_operator(e->kmer, *element_get_key(&b, hash_table->kmer_size, &tmp_kmer2)) );
	    }
	  else
	    {
        die("Error: e is NULL for i=%lld - unable to find\n",i);
	     }
	}
      
      hash_table_free(&hash_table);
      CU_ASSERT(hash_table == NULL);
      
    }

  int num_bits;
  int bucket;
  int max_tries;
  int step;

  for (kmer_size=3; kmer_size<6; kmer_size+=2)//changed from 10 to 6
    {
      num_bits=20;//debig changed from 15 to 20
      bucket=100;
      max_tries=10;
      step=1;
      max_key_given_kmer_size = (bitfield_of_64bits) 1 <<(kmer_size*2);
      test(kmer_size, num_bits, bucket, max_tries, step, max_key_given_kmer_size); 
    }



  //kmer_size=17;
  //num_bits=16;
  //bucket=400;
  //max_tries=20;
  //step=10000000;
  //max_key_given_kmer_size = (BinaryKmer) 1 <<(kmer_size*2);
  //test(kmer_size, num_bits, bucket, max_tries, step, max_key_given_kmer_size);


  /*
  //now test with just one bucket
  kmer_size=3;
  num_bits=20;
  bucket=1;
  max_tries=10;
  step=1;
  max_key_given_kmer_size = (bitfield_of_64bits) 1 <<(kmer_size*2);
  test(kmer_size, num_bits, bucket, max_tries, step, max_key_given_kmer_size);


  //now test with many very shallow buckets
  kmer_size=3;
  num_bits=2;
  bucket=10000000;
  max_tries=10;
  step=1;
  max_key_given_kmer_size = (bitfield_of_64bits) 1 <<(kmer_size*2);
  test(kmer_size, num_bits, bucket, max_tries, step, max_key_given_kmer_size);

  */
  // this takes ages:  
  //and with large kmer but still shallow buckets
  //kmer_size=31;
  //step=200000;
  //bitfield_of_64bits max_kmer_to_test=100000;
  //test(kmer_size, num_bits, bucket, max_tries, step, max_kmer_to_test);
  




}


void test_hash_table_apply_or_insert()
{
  
  short kmer_size;
  bitfield_of_64bits max_key_given_kmer_size;
  

  //adds kmers to hash, then goes through and tries applying to them all, and checks if the apply was applied
  void test(short kmer_size, int num_bits, int bucket, int max_tries, int step, long long max_kmer_to_test)
    {
      
      int number_of_bits      = num_bits;
      int bucket_size         = bucket;
      long long bad_reads     = 0; 
      int max_retries         = max_tries;
      boolean found           = false;
      BinaryKmer tmp_kmer;

      HashTable* hash_table  = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
      
      
      long long i;
      
      //insert all possible binary kmers, prior to start of test
      for (i=0; i< max_kmer_to_test; i++)
	{
          BinaryKmer b;
          binary_kmer_initialise_to_zero(&b);
          b[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]=(bitfield_of_64bits) i;

	  hash_table_find_or_insert(element_get_key(&b, hash_table->kmer_size, &tmp_kmer),&found, hash_table);
	}
      
      
      boolean applied=false;
      
      //now test the "apply"
      for (i=0; i< max_kmer_to_test; i=i+step)
	{
	  BinaryKmer b;
          binary_kmer_initialise_to_zero(&b);
          b[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]=(bitfield_of_64bits) i;

	  applied = hash_table_apply_or_insert(element_get_key(&b,hash_table->kmer_size, &tmp_kmer),&db_node_action_set_status_pruned , hash_table);
	  CU_ASSERT(applied==true);
	  CU_ASSERT(db_node_check_status(hash_table_find(element_get_key(&b,hash_table->kmer_size, &tmp_kmer),hash_table), pruned)==true);
	}
      
      hash_table_free(&hash_table);
      CU_ASSERT(hash_table == NULL);
      
    }

  int num_bits;
  int bucket;
  int max_tries;
  int step;

  for (kmer_size=3; kmer_size<10; kmer_size+=2)
    {
      max_key_given_kmer_size = (bitfield_of_64bits) 1 <<(kmer_size*2);
      num_bits=15;
      bucket=100;
      max_tries=10;
      step=1;
      test(kmer_size, num_bits, bucket, max_tries, step, max_key_given_kmer_size);
    }

  //now test with just one bucket
  kmer_size=3;
  num_bits=20;
  bucket=1;
  max_tries=10;
  step=1;
  max_key_given_kmer_size = (bitfield_of_64bits) 1 <<(kmer_size*2);
  test(kmer_size, num_bits, bucket, max_tries, step, max_key_given_kmer_size);


  //now test with many very shallow buckets
  kmer_size=3;
  num_bits=2;
  bucket=10000000;
  max_tries=10;
  step=1;
  max_key_given_kmer_size = (bitfield_of_64bits) 1 <<(kmer_size*2);
  test(kmer_size, num_bits, bucket, max_tries, step, max_key_given_kmer_size);

  //and with large kmer but still shallow buckets
  kmer_size=31;
  step=200000;
  long long max_kmer_to_test=100000;
  test(kmer_size, num_bits, bucket, max_tries, step, max_kmer_to_test);



  
}
