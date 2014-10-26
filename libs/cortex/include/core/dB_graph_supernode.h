/*
 * Copyright 2009-2013 Zamin Iqbal and Mario Caccamo
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
  dB_graph_supernode.h defines the interface for dealing with supernodes
*/

#ifndef DB_GRAPH_SUPERNODE_H_
#define DB_GRAPH_SUPERNODE_H_

#include "global.h"
#include "element.h"
#include "open_hash/hash_table.h"

//a Supernode object just holds an array of nodes
// and auxiliary objects.
typedef struct
{
  dBNode**     path;
  Orientation* or;
  Nucleotide*  nuc;
  char*        seq;
  boolean      is_cycle;
  int          len_sup;//number of edges. 
                       //if one node, length=0, if 2 nodes, length=1, runs from 0,1.., len_sup
  int          len_alloced;
  uint16_t     kmer_size;
} Supernode;


typedef struct
{
  Covg* covgs;
  int len;//above this, data may be garbage in general. 
  int len_alloced;
} CovgArray;

typedef struct
{
  double* llks;
  int len;//above this, data may be garbage in general. 
  int len_alloced;
} LlkArray;

Supernode*   alloc_supernode(int length);
void         init_supernode(Supernode* s, uint16_t kmer_size);
Supernode*   alloc_and_init_supernode(int length, uint16_t kmer_size);
void         free_supernode(Supernode* csup);


//covg arrays
CovgArray*   alloc_and_init_covg_array(int len);
void         free_covg_array(CovgArray* c);
void         reset_covg_array(CovgArray* c);
void         reset_used_part_of_covg_array(CovgArray* c);
void         covg_array_push(CovgArray* c, Covg val);
//likelihood arrays
LlkArray*    alloc_and_init_llk_array(int len);
void         free_llk_array(LlkArray* c);
void         reset_llk_array(LlkArray* c);
void         reset_used_part_of_llk_array(LlkArray* c);
void         llk_array_push(LlkArray* c, double val);




void supernode_assign(Supernode* sup_from, Supernode* sup_to);
void supernode_get_seq(Supernode* sup, StrBuf* seq, boolean include_first_kmer,  uint16_t kmer_size);
void supernode_copy_reverse(Supernode* sup_from, Supernode* sup_to);


#endif /* DB_GRAPH_SUPERNODE_H_ */
