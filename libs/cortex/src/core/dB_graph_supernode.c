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
  dB_graph_supernode.c  - functions for handling supernodes
*/

// system libraries
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>


// cortex_var headers
#include "binary_kmer.h"
#include "dB_graph_supernode.h"



Supernode* alloc_supernode(int length)
{
  Supernode* s = (Supernode*)calloc(sizeof(Supernode), 1 );
  if (s==NULL)
    {
      return s;
    }

  s->path= (dBNode**)     calloc(sizeof(dBNode*), length);
  if (s->path==NULL)
    {
      free(s);
      return NULL;
    }

  s->or  = (Orientation*) calloc(sizeof(Orientation),length);
  if (s->or==NULL)
    {
      free(s->path);
      free(s);
      return NULL;
    }

  s->nuc = (Nucleotide*)  calloc(sizeof(Nucleotide),length);
  if (s->nuc==NULL)
    {
      free(s->or);
      free(s->path);
      free(s);
      return NULL;
    }

  s->seq = (char*)        calloc(sizeof(char),(length+1));
  if (s->seq==NULL)
    {
      free(s->nuc);
      free(s->or);
      free(s->path);
      free(s);
      return NULL;
    }
  s->len_alloced=length;
  return s;
}

void init_supernode(Supernode* s, uint16_t kmer_size)
{
  if (s==NULL)
    {
      return;
    }
  s->is_cycle=false;
  s->len_sup=0;
  s->kmer_size = kmer_size;
}

Supernode* alloc_and_init_supernode(int length, uint16_t kmer_size)
{
  Supernode* s= alloc_supernode(length);
  if (s!=NULL)
    {
      init_supernode(s, kmer_size);
    }
  return s;
}

void free_supernode(Supernode* s)
{
  free(s->path);
  free(s->or);
  free(s->nuc);
  free(s->seq);
  free(s);
}





CovgArray*   alloc_and_init_covg_array(int len_alloced)
{
  CovgArray* ret = (CovgArray*) malloc(sizeof(CovgArray));
  if (ret==NULL)
    {
      return NULL;
    }
  ret->covgs = (Covg*) malloc(sizeof(Covg)*len_alloced);
  if (ret->covgs==NULL)
    {
      return NULL;
    }

  //  memset(ret->covgs, 0, sizeof(Covg) * len);
  ret->len_alloced = len_alloced;
  ret->len = 0;
  reset_covg_array(ret);
  return ret;
}

void free_covg_array(CovgArray* c)
{
  free(c->covgs);
  free(c);
}

void reset_covg_array(CovgArray* c)
{
  int i;
  for (i=0; i<c->len_alloced; i++)
    {
      c->covgs[i]= (Covg)0;
    }
  c->len=0;
}

void reset_used_part_of_covg_array(CovgArray* c)
{
  int i;
  for (i=0; i<c->len; i++)
    {
      c->covgs[i]= (Covg)0;
    }
  c->len=0;
}

void covg_array_push(CovgArray* c, Covg val)
{
  if (c->len+1 <= c->len_alloced)
    {
      c->covgs[c->len] = val;
      c->len = c->len+1;
    }
  else
    {
      die("Out of space in covg array - hit limit of %d\n", c->len);
    }
}


LlkArray*   alloc_and_init_llk_array(int len_alloced)
{
  LlkArray* ret = (LlkArray*) malloc(sizeof(LlkArray));
  if (ret==NULL)
    {
      die("Unable to alloc a LlkArray holder\n");
    }
  ret->llks = (double*) malloc(sizeof(double)*len_alloced);
  if (ret->llks==NULL)
    {
      die("Unable to malloc an array of llks. I don't see any value in dying gracefully here, you must be badly out of memory\n");
    }
  ret->len_alloced = len_alloced;
  ret->len = 0;
  reset_llk_array(ret);
  return ret;
}

void free_llk_array(LlkArray* la)
{
  free(la->llks);
  free(la);
}

void reset_llk_array(LlkArray* la)
{
  int i;
  for (i=0; i<la->len_alloced; i++)
    {
      la->llks[i]= (double)0;
    }
  la->len=0;
}

void reset_used_part_of_llk_array(LlkArray* la)
{
  int i;
  for (i=0; i<la->len; i++)
    {
      la->llks[i]= (double)0;
    }
  la->len=0;
}

void llk_array_push(LlkArray* la, double val)
{
  if (la->len+1 <= la->len_alloced)
    {
      la->llks[la->len] = val;
      la->len = la->len+1;
    }
  else
    {
      die("Out of space in llk array - hit limit of %d\n", la->len);
    }
}



void supernode_assign(Supernode* sup_from, Supernode* sup_to)
{
  if (sup_to->len_alloced < sup_from->len_sup)
    {
      die("Trying to copy a long supernode into a small one. Coding error.\n");
    }
  int i;
  for (i=0; i<=sup_from->len_sup; i++)
    {
      sup_to->path[i]  = sup_from->path[i];
      sup_to->or[i]    = sup_from->or[i];
      sup_to->nuc[i]   = sup_from->nuc[i];
    }
  sup_to->is_cycle     = sup_from->is_cycle;
  sup_to->len_sup      = sup_from->len_sup;
  sup_to->seq[0]       = '\0';
  strcat(sup_to->seq, sup_from->seq);
  sup_to->kmer_size = sup_from->kmer_size;
  
}



void supernode_get_seq(Supernode* sup, StrBuf* seq, boolean include_first_kmer,  uint16_t kmer_size)
{
  strbuf_reset(seq);

  char tmp[kmer_size+1];
  tmp[kmer_size]='\0';
  if (include_first_kmer==true)
    {
      if (sup->or[0]==forward)
	{
	  strbuf_append_str(seq, binary_kmer_to_seq(&(sup->path[0]->kmer), kmer_size, tmp));
	}
      else
	{
	  BinaryKmer local_copy_of_kmer;
	  binary_kmer_assignment_operator(local_copy_of_kmer, sup->path[0]->kmer);
	  BinaryKmer tmp_kmer;
	  // after the following line tmp_kmer and rev_kmer are pointing to the same B Kmer
	  BinaryKmer* rev_kmer = 
	    binary_kmer_reverse_complement(&local_copy_of_kmer, kmer_size, &tmp_kmer);
	  strbuf_append_str(seq, binary_kmer_to_seq(rev_kmer, kmer_size, tmp));
	}
    }
  strbuf_append_str(seq, sup->seq);
}


void supernode_copy_reverse(Supernode* sup_from, Supernode* sup_to)
{
  //I'm not going to check if the kmer_size is the same from from and to
  if (sup_to->len_alloced < sup_from->len_sup)
    {
      die("Trying to copy-reverse a long supernode into a small one. Coding error.\n");
    }
  

  BinaryKmer local_copy_of_kmer;
  BinaryKmer tmp_kmer;

  //local function so I don't keep buggering around on the stack
  //suppose this node goes to nodes A and B (in direction current_orientation)
  //if I was coming back from A or B in the reverse direction, in both cases the final edge
  //would be the same: let's get it
  Nucleotide db_graph_get_reverse_edge(dBNode * current_node, Orientation current_orientation)
  {
    
    binary_kmer_assignment_operator(local_copy_of_kmer, current_node->kmer);
  
    // after the following line tmp_kmer and rev_kmer are pointing to the same B Kmer
    BinaryKmer* rev_kmer = binary_kmer_reverse_complement(&local_copy_of_kmer,sup_from->kmer_size, &tmp_kmer);
    
    Nucleotide reverse_edge;
    if (current_orientation == reverse)
      {   
	reverse_edge = binary_kmer_get_last_nucleotide(&local_copy_of_kmer);
	binary_kmer_assignment_operator(local_copy_of_kmer,*rev_kmer);
      }
    else
      {
	reverse_edge = binary_kmer_get_last_nucleotide(rev_kmer);
      }

    return reverse_edge;
  }
  //end of local function


  int i;
  sup_to->is_cycle   = sup_from->is_cycle;
  sup_to->len_sup    = sup_from->len_sup;
  int len            = sup_from->len_sup;
  for (i=len; i>=0; i--)
    {
      sup_to->path[len-i]= sup_from->path[i];
      sup_to->or[len-i]  = opposite_orientation(sup_from->or[i]);//going backwards
      
      //now for reverse bases and nucleotides
      if (i==len)
	{
	  //sup_to->nuc[len-i]=Undefined;
	}
      else
	{
	  sup_to->nuc[len-i] = 
	    db_graph_get_reverse_edge(sup_from->path[i], sup_from->or[i]);
	}

    }
  sup_to->seq[0]='\0';
  nucleotides_to_string(sup_to->nuc, sup_to->len_sup-1, sup_to->seq);

}



