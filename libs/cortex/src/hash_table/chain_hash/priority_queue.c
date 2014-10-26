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
  chain_hash/priority_queue.c
*/

#include <stdlib.h>
#include <stdio.h>

#include "priority_queue.h"

PQueue * pqueue_new(){ 
  
  PQueue * pqueue = malloc(sizeof(PQueue));
  
  if (pqueue == NULL){
    die("Could not allocate pqueue");
  }

  pqueue->max_size = 0;
  pqueue->number_entries = 0;

  pqueue->elements = NULL;

  return pqueue;
  

}


//internal routine that reallocs space for the array represenintg the priority queue
// return true if reallocates
boolean pqueue_realloc(PQueue * pqueue, int realloc_size){
  //check if it needs to reallocate
  boolean ret = false;

  if (pqueue->number_entries+1>pqueue->max_size){
    Element * new_ptr;
    new_ptr = realloc(pqueue->elements,(pqueue->number_entries+realloc_size) * sizeof(Element));
    
    pqueue->max_size += realloc_size;
    
    if (new_ptr == NULL){
      die("priority queue: cannot allocate memory");
      //so no need to worry about orphaned pointer at this stage.
    }
    
    pqueue->elements = new_ptr;
    ret = true;
  }

  return ret;
}

//assumes the element is not already stored -- check this before calling this function
void pqueue_bubble_up(int current_index, PQueue *pqueue){
  
  while (current_index > 0){

    int parent_index = (current_index-1)/2;
    if (element_smaller(pqueue->elements[parent_index],pqueue->elements[current_index])){
      Element tmp;
      tmp = pqueue->elements[current_index];
      pqueue->elements[current_index] = pqueue->elements[parent_index];
      pqueue->elements[parent_index] = tmp;
      current_index = parent_index;
    }
    else{
      break;
    }
  }
}

boolean pqueue_apply_or_insert(Key key, void (*f)(Element*),PQueue *pqueue, short kmer_size){

  int i =0;
  boolean found = false;
  Element element;

  for(i=0;i<pqueue->number_entries;i++){
    if (element_is_key(key,pqueue->elements[i])){
      f(&pqueue->elements[i]);
      pqueue_bubble_up(i,pqueue);
      found = true;
      break;
    }
  }

  if (! found){
    pqueue_realloc(pqueue,2);
    pqueue->number_entries++;
    element_initialise(&element,key, kmer_size);
    pqueue->elements[pqueue->number_entries-1] = element; //structure assignment
  }
  return found;
}

void pqueue_traverse(void (*f)(Element *),PQueue * pqueue)
{
  int i;
  for(i=0;i<pqueue->number_entries;i++){
    f(&(pqueue->elements[i]));
  }
}


Element * pqueue_find_or_insert(Key key,boolean * found, PQueue * pqueue, short kmer_size, int alloc_size){  
  int i;
  Element element;
 
  *found = false;
  for(i=0;i<pqueue->number_entries;i++){ 
    if (element_is_key(key,pqueue->elements[i])){
      *found = true;
      return &(pqueue->elements[i]);
    }
  }
  
  if (!*found){
    pqueue_realloc(pqueue,alloc_size);
    pqueue->number_entries++;
    element_initialise(&element,key, kmer_size);
    pqueue->elements[pqueue->number_entries-1] = element; //structure assignment
  }

  return &(pqueue->elements[pqueue->number_entries-1]);  
}


Element * pqueue_insert(Key key,PQueue * pqueue, short kmer_size, int alloc_size){  
 
  Element element;
  pqueue_realloc(pqueue,alloc_size);
  pqueue->number_entries++;
  element_initialise(&element,key, kmer_size);
  pqueue->elements[pqueue->number_entries-1] = element; //structure assignment
  return &(pqueue->elements[pqueue->number_entries-1]);  
}


Element * pqueue_find(Key key,PQueue * pqueue, short kmer_size){
  
  int i;

  for(i=0;i<pqueue->number_entries;i++){ 
    if (element_is_key(key,pqueue->elements[i])){
      return &(pqueue->elements[i]);
    }
  }
  
  return NULL;
}

void pqueue_free(PQueue** pqueue){
  
  free((*pqueue)->elements);
  
  free(*pqueue);
  *pqueue = NULL;
  
}

void pqueue_free_elements(PQueue* pqueue){
  
  free(pqueue->elements);
  pqueue->elements = NULL;  
}





