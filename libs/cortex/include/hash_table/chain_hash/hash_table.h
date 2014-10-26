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
  hash_table.h 
  
  chain hash table implementation - ie every bucket grows as more elements are stored
  bucket implemeted as priority queues
  all the routines as prefixed with hash_table
*/


#ifndef CHAIN_HASH_TABLE_H_
#define CHAIN_HASH_TABLE_H_

#include "global.h"
#include "element.h"

typedef struct
{
  short kmer_size;
  unsigned int number_buckets;
  PQueue * table; //every bucket is implemented as a priority queue
} HashTable;


HashTable * hash_table_new(unsigned int number_buckets, short kmer_size);

void hash_table_free(HashTable * * hash_table);

//if the key is present applies f otherwise adds a new element for kmer
boolean hash_table_apply_or_insert(Key key, void (*f)(Element*), HashTable *);

//applies f to every element of the table
void hash_table_traverse(void (*f)(Element *),HashTable *);

//if the element is not in table create an element with key and adds it
Element * hash_table_find_or_insert(Key key, boolean * found, HashTable * hash_table);

//return entry for kmer
Element * hash_table_find(Key key, HashTable * hash_table);

Element * hash_table_insert(Key key, HashTable * hash_table);

#endif /* CHAIN_HASH_TABLE_H_ */
