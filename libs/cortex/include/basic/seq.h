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
