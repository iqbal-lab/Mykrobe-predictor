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
  seq.h
*/

#ifndef STDSEQ_H_
#define STDSEQ_H_

#include <stdio.h>

#include "global.h"

typedef struct
{
  char *name;
  int  start,end;
  char *seq;  // sequence 
  char* qual; // qualities 
} Sequence;

//returns length of sequence read

//note that fastq format dosn't support "partial" reading of a long entry -- ask Zam or Mario
//that means we can only read full fastq entries
int read_sequence_from_fastq(FILE * fp, Sequence * seq, int max_read_length, int fastq_ascii_offset);

//this routine can read long sequences (eg full chromosomes) , this is implemented by reading the sequence in chunks
int read_sequence_from_fasta(FILE * fp, Sequence * seq, int max_chunk_length, boolean new_entry, boolean * full_entry, int offset);

void alloc_sequence(Sequence * seq, int max_read_length, int max_name_length);

void free_sequence(Sequence ** );

void shift_last_kmer_to_start_of_sequence(Sequence * sequence, int length, short kmer_size);
boolean good_base(char c);

#endif /* STDSEQ_H_ */
