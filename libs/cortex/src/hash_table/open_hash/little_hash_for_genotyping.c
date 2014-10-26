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
  open_hash/little_hash_for_genotyping.c
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "open_hash/little_hash_for_genotyping.h"
#include "hash_value.h"


LittleHashTable * little_hash_table_new(int number_bits, int bucket_size, 
				  int max_rehash_tries, short kmer_size){ 
  
  LittleHashTable *little_hash_table = malloc(sizeof(LittleHashTable));

  if (little_hash_table == NULL) {
    fprintf(stderr,"could not allocate hash table\n");
    return NULL;
    //exit(EXIT_FAILURE);
  }
  
  little_hash_table->collisions = calloc(max_rehash_tries, sizeof(long long));
  if (little_hash_table->collisions == NULL) {
    fprintf(stderr,"could not allocate memory\n");
    return NULL;
    //exit(EXIT_FAILURE);
  }
  
  little_hash_table->unique_kmers = 0;
  little_hash_table->max_rehash_tries = max_rehash_tries;
  little_hash_table->number_buckets = (long long) 1 << number_bits;
  little_hash_table->bucket_size   = bucket_size;

  //calloc is vital - we want to make sure initialised to zero
  little_hash_table->table = calloc(little_hash_table->number_buckets * little_hash_table->bucket_size, sizeof(GenotypingElement));

  if (little_hash_table->table == NULL) {
    fprintf(stderr,"could not allocate hash table of size %qd\n",little_hash_table->number_buckets * little_hash_table->bucket_size);
    return NULL;
    //exit(EXIT_FAILURE);
  }
   
  little_hash_table->next_element = calloc(little_hash_table->number_buckets, sizeof(int));
  if (little_hash_table->table == NULL) {
    fprintf(stderr,"could not allocate array of pointers for next available genotyping_element in buckets [%qd]\n",little_hash_table->number_buckets);
    return NULL;
    //exit(EXIT_FAILURE);
  }

  little_hash_table->kmer_size      = kmer_size;
  return little_hash_table;
}

void little_hash_table_free(LittleHashTable ** little_hash_table)
{ 
  free((*little_hash_table)->table);
  free((*little_hash_table)->next_element);
  free((*little_hash_table)->collisions);
  free(*little_hash_table);
  *little_hash_table = NULL;
}


// Lookup for key in bucket defined by the hash value. 
// If key is in bucket, returns true and the position of the key/element in current_pos.
// If key is not in bucket, and bucket is not full, returns the next available position in current_pos (and overflow is returned as false)
// If key is not in bucket, and bucket is full, returns overflow=true
boolean little_hash_table_find_in_bucket(Key key, long long * current_pos, boolean * overflow, LittleHashTable * little_hash_table, int rehash){


  //add the rehash to the final bitfield in the BinaryKmer
  BinaryKmer bkmer_with_rehash_added;
  binary_kmer_initialise_to_zero(&bkmer_with_rehash_added);
  binary_kmer_assignment_operator(bkmer_with_rehash_added, *key);
  bkmer_with_rehash_added[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] =   bkmer_with_rehash_added[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]+ (bitfield_of_64bits) rehash;

  int hashval = hash_value(&bkmer_with_rehash_added,little_hash_table->number_buckets);


  boolean found = false;
  int i=0;                     //position in bucket
  *overflow    = false;
  *current_pos   = (long long) hashval * little_hash_table->bucket_size;   //position in hash table

  while( (i<little_hash_table->bucket_size) &&   // still within the bucket
	 (!db_genotyping_node_check_for_flag_ALL_OFF(&little_hash_table->table[*current_pos]) )  && // not yet reached an empty space
	 (!found)
	 )
    {

    //sanity check -- to avoid out of boundary access
    if (*current_pos >= little_hash_table->number_buckets * little_hash_table->bucket_size || *current_pos<0)
      {
	die("out of bounds problem found in hash table_find_with_position\n");
      }
    
    //element found

   
    if (genotyping_element_is_key(key,little_hash_table->table[*current_pos]))
      {
	found = true;
      }
    else
      {
	(*current_pos)++;
	i++;
      }
    
    }
  

  if (i == little_hash_table->bucket_size)
    {
      *overflow = true;
    }
  

  assert(!found || !(*overflow));
  return found;
}




void little_hash_table_traverse(void (*f)(GenotypingElement *),LittleHashTable * little_hash_table){
  long long i;
 
  for(i=0;i<little_hash_table->number_buckets * little_hash_table->bucket_size;i++){
    
    if (!db_genotyping_node_check_for_flag_ALL_OFF(&little_hash_table->table[i])){
      f(&little_hash_table->table[i]);
    }
  }
}

long long little_hash_table_traverse_returning_sum(long long (*f)(GenotypingElement *),LittleHashTable * little_hash_table){
  long long i;
  long long ret=0;
  for(i=0;i<little_hash_table->number_buckets * little_hash_table->bucket_size;i++){
    if (!db_genotyping_node_check_for_flag_ALL_OFF(&little_hash_table->table[i])){
      ret += f(&little_hash_table->table[i]);
    }
  }
  return ret;
}

void little_hash_table_traverse_passing_int(void (*f)(GenotypingElement *, int*),LittleHashTable * little_hash_table, int* num){
  long long i;
  for(i=0;i<little_hash_table->number_buckets * little_hash_table->bucket_size;i++){
    if (!db_genotyping_node_check_status(&little_hash_table->table[i],unassigned)){
      f(&little_hash_table->table[i], num);
    }
  }
}
void little_hash_table_traverse_passing_ints_and_path(void (*f)(GenotypingElement *, int*, int*, dBNode**, Orientation*, Nucleotide*, char*, int),
					       LittleHashTable * little_hash_table, int* num1, int* num2, 
					       dBNode** p_n, Orientation* p_o, Nucleotide* p_lab, char* p_str, int len){
  long long i;
  for(i=0;i<little_hash_table->number_buckets * little_hash_table->bucket_size;i++){
    if (!db_genotyping_node_check_status(&little_hash_table->table[i],unassigned)){
      f(&little_hash_table->table[i], num1, num2, p_n, p_o, p_lab, p_str, len);
    }
  }
}


void little_hash_table_traverse_passing_3ints_and_big_graph_path(void (*f)(GenotypingElement *, int*, int*, int*, dBNode**, Orientation*, Nucleotide*, char*, int),
								 LittleHashTable * little_hash_table, int* num1, int* num2, int* num3, 
								 dBNode** p_n, Orientation* p_o, Nucleotide* p_lab, char* p_str, int len){
  long long i;
  for(i=0;i<little_hash_table->number_buckets * little_hash_table->bucket_size;i++){
    if (!db_genotyping_node_check_status(&little_hash_table->table[i],unassigned)){
      f(&little_hash_table->table[i], num1, num2, num3, p_n, p_o, p_lab, p_str, len);
    }
  }
}

void little_hash_table_traverse_passing_big_graph_path(void (*f)(GenotypingElement *, dBNode**, Orientation*, Nucleotide*, char*, int),
								 LittleHashTable * little_hash_table, 
								 dBNode** p_n, Orientation* p_o, Nucleotide* p_lab, char* p_str, int len){
  long long i;
  for(i=0;i<little_hash_table->number_buckets * little_hash_table->bucket_size;i++){
    if (!db_genotyping_node_check_status(&little_hash_table->table[i],unassigned)){
      f(&little_hash_table->table[i], p_n, p_o, p_lab, p_str, len);
    }
  }
}


GenotypingElement * little_hash_table_find(Key key, LittleHashTable * little_hash_table)
{
  if (little_hash_table == NULL) 
    {
      die("little_hash_table_find has been called with a NULL table! Exiting");
    }

  GenotypingElement * ret = NULL;
  long long current_pos;
  boolean overflow;
  int rehash = 0;
  boolean found; 

  do
    {
      found = little_hash_table_find_in_bucket(key,&current_pos, &overflow, little_hash_table,rehash);
      
      if (found) //then we know overflow is false - this is checked in find_in_bucket
	{
	  ret =  &little_hash_table->table[current_pos];
	}
      else if (overflow)
	{ //rehash
	  rehash++; 
	  if (rehash>little_hash_table->max_rehash_tries)
	    {
	      die("Too much rehashing!! Rehash=%d", rehash);
	    }
	}
    } while(overflow);
  
  return ret;
}


GenotypingElement * little_hash_table_find_or_insert(Key key, boolean * found,  LittleHashTable * little_hash_table){
  
  if (little_hash_table == NULL) {
    die("NULL table!");
  }
  
  GenotypingElement element;
  GenotypingElement * ret = NULL;
  int rehash = 0;
  boolean overflow; 

  long long current_pos;

  do{

    *found = little_hash_table_find_in_bucket(key,&current_pos,&overflow,little_hash_table,rehash);

    if (! *found)
      {
	if (!overflow) //it is definitely nowhere in the hashtable, so free to insert
	  {
	    //sanity check
	    if (!db_genotyping_node_check_for_flag_ALL_OFF(&little_hash_table->table[current_pos]))
	      {
		die("error trying to write on an occupied element\n");
	      }
      
	    //insert element
	    //printf("Inserting element at position %qd in bucket \n", current_pos);
	    genotyping_element_initialise(&element,key, little_hash_table->kmer_size);

	    //little_hash_table->table[current_pos] = element; //structure assignment
	    genotyping_element_assign(&(little_hash_table->table[current_pos]) , &element);
	    
	    ret = &little_hash_table->table[current_pos];
	    little_hash_table->unique_kmers++;
	    
	  }
	else
	  { //overflow -> rehashing
	    
	    rehash++;
	    if (rehash>little_hash_table->max_rehash_tries)
	      {
		die("too much rehashing!! Reserve more memory. Rehash=%d\n", rehash);
	      }
	  }
      }
    else //it is found
      {
	ret = &little_hash_table->table[current_pos];
      }
  } while (overflow);
  
  little_hash_table->collisions[rehash]++;
  return ret;
}


//this methods inserts an element in the next available bucket
//it doesn't check whether another element with the same key is present in the table
//used for fast loading when it is known that all the elements in the input have different key
GenotypingElement * little_hash_table_insert(Key key, LittleHashTable * little_hash_table){
  
  if (little_hash_table == NULL) {
    die("NULL table!");
  }
  
  GenotypingElement element;
  GenotypingElement * ret = NULL;
  int rehash = 0;
  boolean inserted = false;
  do{
    //add the rehash to the final bitfield in the BinaryKmer
    BinaryKmer bkmer_with_rehash_added;
    binary_kmer_initialise_to_zero(&bkmer_with_rehash_added);
    binary_kmer_assignment_operator(bkmer_with_rehash_added, *key);
    bkmer_with_rehash_added[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] =   bkmer_with_rehash_added[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]+ (bitfield_of_64bits) rehash;

    int hashval = hash_value(&bkmer_with_rehash_added,little_hash_table->number_buckets);
    
    if (little_hash_table->next_element[hashval] < little_hash_table->bucket_size)
      { //can insert element
	long long  current_pos   = (long long) hashval * little_hash_table->bucket_size + (long long) little_hash_table->next_element[hashval] ;   //position in hash table

	//sanity check
	if (!db_genotyping_node_check_for_flag_ALL_OFF(&little_hash_table->table[current_pos])){
	  die("Out of bounds - trying to insert new node beyond end of bucket\n");
	}
  
      
	genotyping_element_initialise(&element,key, little_hash_table->kmer_size);
	genotyping_element_assign( &(little_hash_table->table[current_pos]),  &element); 
	little_hash_table->unique_kmers++;
	little_hash_table->next_element[hashval]++;	
	ret = &little_hash_table->table[current_pos];
	inserted=true;
      }
    else
      {//rehash
	rehash++;
	if (rehash>little_hash_table->max_rehash_tries)
	  {
	    die("Too much rehashing!! Reserve more memory.  Rehash=%d\n", rehash);
	  }
      }
    
  } while (! inserted);

  return ret;
}

void little_hash_table_print_stats(LittleHashTable * little_hash_table)
{
  int k;
  printf("Collisions:\n");
  for(k=0;k<10;k++)
    {
      if (little_hash_table->collisions[k] != 0){
	printf("\t tries %i: %qd\n",k,little_hash_table->collisions[k]);
      }
    }
}



long long little_hash_table_get_unique_kmers(LittleHashTable * little_hash_table)
{
  return little_hash_table->unique_kmers;
}


long long little_hash_table_get_capacity(LittleHashTable * little_hash_table){
  return little_hash_table->number_buckets*little_hash_table->bucket_size;
}


//wipes a colour clean - for all nodes, sets covg=0, edge=0.
void little_graph_wipe_colour(int colour, LittleHashTable* little_graph)
{
  void wipe_node(GenotypingElement* node)
  {
    node->individual_edges[colour]=0;
    node->coverage[colour]=0;
  }
  little_hash_table_traverse(&wipe_node, little_graph);
}

