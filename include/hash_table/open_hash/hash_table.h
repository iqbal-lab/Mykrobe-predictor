/*
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
  hash_table.h 
  
  open hash table implementation - ie every bucket has a predifined size 
  overloads results in rehashing
  all the routines as prefixed with hash_table
*/

#ifndef OPEN_HASH_TABLE_H_
#define OPEN_HASH_TABLE_H_

#include "global.h"
#include "element.h"

typedef struct
{
  short kmer_size;
  long long number_buckets;
  int bucket_size;
  Element * table; 
  short * next_element; //keeps index of the next free element in bucket 
  long long * collisions;
  long long unique_kmers;
  int max_rehash_tries;
} HashTable;


HashTable * hash_table_new(int number_bits, int bucket_size, int max_rehash_tries, short kmer_size);

void hash_table_free(HashTable * * hash_table);

//if the key is present applies f otherwise adds a new element for kmer
boolean hash_table_apply_or_insert(Key key, void (*f)(Element*), HashTable *);

//applies f to every element of the table
void hash_table_traverse(void (*f)(Element *),HashTable *);
long long hash_table_traverse_returning_sum(long long (*f)(Element *),HashTable * hash_table);
void hash_table_traverse_passing_int(void (*f)(Element *, int*),HashTable * hash_table, int* num);
void hash_table_traverse_passing_ints_and_path(
  void (*f)(Element *, int*, int*, dBNode**, Orientation*, Nucleotide*, char*, int),
  HashTable * hash_table, int* num1, int* num2, 
  dBNode** p_n, Orientation* p_o, Nucleotide* p_lab, char* p_str, int len);

void hash_table_traverse_passing_3ints_and_path(
  void (*f)(Element *, int*, int*, int*, dBNode**, Orientation*, Nucleotide*, char*, int),
  HashTable * hash_table, int* num1, int* num2, int* num3, 
  dBNode** p_n, Orientation* p_o, Nucleotide* p_lab, char* p_str, int len);


//if the element is not in table create an element with key and adds it
Element * hash_table_find_or_insert(Key key, boolean * found, HashTable * hash_table);
Element * hash_table_insert(Key key, HashTable * hash_table);

void hash_table_print_stats(HashTable *);

long long hash_table_get_unique_kmers(HashTable *);

//return entry for kmer
Element * hash_table_find(Key key, HashTable * hash_table);

long long hash_table_get_capacity(HashTable * hash_table);

#endif /* OPEN_HASH_TABLE_H_ */
