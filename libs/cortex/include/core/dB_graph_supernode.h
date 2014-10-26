/*
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@tgac.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * The MIT License (MIT)
 * Copyright (c) 2009-2014 <Z. Iqbal and M. Caccamo>
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
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
