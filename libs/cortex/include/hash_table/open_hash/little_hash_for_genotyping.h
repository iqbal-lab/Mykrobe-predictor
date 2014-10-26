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
  little_hash_for_genotyping.h
*/

#ifndef LITTLE_HASH_H_
#define LITTLE_HASH_H_

#include "global.h"
#include "genotyping_element.h"
#include "element.h"

typedef struct
{
  short kmer_size;
  long long number_buckets;
  int bucket_size;
  GenotypingElement * table; 
  short * next_element; //keeps index of the next free element in bucket 
  long long * collisions;
  long long unique_kmers;
  int max_rehash_tries;
} LittleHashTable;


LittleHashTable * little_hash_table_new(int number_bits, int bucket_size, int max_rehash_tries, short kmer_size);

void little_hash_table_free(LittleHashTable * * little_hash_table);

//if the key is present applies f otherwise adds a new element for kmer
boolean little_hash_table_apply_or_insert(Key key, void (*f)(GenotypingElement*), LittleHashTable *);

//applies f to every element of the table
void little_hash_table_traverse(void (*f)(GenotypingElement *),LittleHashTable *);
long long little_hash_table_traverse_returning_sum(long long (*f)(GenotypingElement *),LittleHashTable * little_hash_table);
void little_hash_table_traverse_passing_int(void (*f)(GenotypingElement *, int*),LittleHashTable * little_hash_table, int* num);
void little_hash_table_traverse_passing_ints_and_path(void (*f)(GenotypingElement *, int*, int*, dBNode**, Orientation*, Nucleotide*, char*, int),
					       LittleHashTable * little_hash_table, int* num1, int* num2, 
					       dBNode** p_n, Orientation* p_o, Nucleotide* p_lab, char* p_str, int len);


//note the mix of GenotypingElement and Element in the arguments - this is DELIBERATE
void little_hash_table_traverse_passing_3ints_and_big_graph_path(void (*f)(GenotypingElement *, int*, int*, int*, dBNode**, Orientation*, Nucleotide*, char*, int),
								 LittleHashTable * little_hash_table, int* num1, int* num2, int* num3, 
								 dBNode** p_n, Orientation* p_o, Nucleotide* p_lab, char* p_str, int len);

void little_hash_table_traverse_passing_big_graph_path(void (*f)(GenotypingElement *, dBNode**, Orientation*, Nucleotide*, char*, int),
						       LittleHashTable * little_hash_table, 
						       dBNode** p_n, Orientation* p_o, Nucleotide* p_lab, char* p_str, int len);
//if the element is not in table create an element with key and adds it
GenotypingElement * little_hash_table_find_or_insert(Key key, boolean * found, LittleHashTable * little_hash_table);
GenotypingElement * little_hash_table_insert(Key key, LittleHashTable * little_hash_table);

void little_hash_table_print_stats(LittleHashTable *);

long long little_hash_table_get_unique_kmers(LittleHashTable *);

//return entry for kmer
GenotypingElement * little_hash_table_find(Key key, LittleHashTable * little_hash_table);

long long little_hash_table_get_capacity(LittleHashTable * little_hash_table);

void little_graph_wipe_colour(int colour, LittleHashTable* little_graph);

#endif /* LITTLE_HASH_H_ */
