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
  priority_queue.h - routines are prefixed with priority_queue_
*/

#ifndef PQUEUE_H_
#define PQUEUE_H_

#include "global.h"
#include "element.h"

typedef struct
{
  int number_entries;
  int max_size;
  Element * elements;
} PQueue;


PQueue * pqueue_new();

boolean pqueue_apply_or_insert(Key, void (*f)(Element*), PQueue *, short kmer_size);

Element * pqueue_find_or_insert(Key key,boolean * found, PQueue *pqueue, short kmer_size, int alloc_size);

//it doesnt check for existence
Element * pqueue_insert(Key key, PQueue *pqueue, short kmer_size, int alloc_size);

void pqueue_traverse(void (*f)(Element *),PQueue *);

Element * pqueue_find(Key,PQueue *, short kmer_size);

void pqueue_free(PQueue** pqueue);
void pqueue_free_elements(PQueue* pqueue);

#endif /* PQUEUE_H_ */
