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

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#include "seq.h"

// Read sequence from file "fp" in FASTA format.
// it returns the length of the sequence, 0 if no sequence is left in file
boolean good_base(char c);


char current_entry_name[LINE_MAX+1] = "";
int  last_end_coord;


//new_entry tells the parser to expect a new fasta entry (ie starts with >)
//full_entry tells the caller if the end of the entry has been reached
//This function can be use with big fasta entries. Successive calls on the same
//big entry remember the name of the entry in seq->name.
//seq->start, seq->end are the coordinate with respect to the full fasta entry (eg chr1 start: 1000 end: 3000)
//offset defines the length of seq->seq that is preserved, this allows to append new sequence from the entry to a prefix in seq->seq
//until max_chunk_length is reached. 
//ir returns 0 when reaches EOF (in this case full_entry is true)

int read_sequence_from_fasta(FILE *fp, Sequence * seq, int max_chunk_length,boolean new_entry, boolean * full_entry, int offset){

  char line[LINE_MAX]; //LINE_MAX is defined in limits.h
  int i;
  int j = 0; //length of sequence
  boolean good_read;

  if (fp == NULL){
    die("read_sequence_from_fasta.c: File not defined");
  }

  if (seq == NULL){
    die("read_sequence_from_fasta.c: Cannot pass a NULL pointer for seq");
  }

  if (seq->seq == NULL){
    die("Dont give me a null pointer for seq->seq - alloc memory yourself and give me that\n");
  }

  if (seq->qual == NULL){
    die("Dont give me a null pointer for seq->qual - alloc memory yourself and give me that\n");
  }

  if (seq->name == NULL){
    die("Dont give me a null pointer for seq->name - alloc memory yourself and give me that\n");
  }

 
  //get name
  if (new_entry == true){ //need a '>' followed by a name
    if (fgets(line, LINE_MAX, fp) != NULL){
      if (line[0] == '>'){
	i=1;
	while(i<LINE_MAX && ! (line[i] == '\n' || line[i] == ' ' || line[i] == '\t' || line[i] == '\r')){
	  seq->name[i-1] = line[i];
	  i++;
	}
	seq->name[i-1] = '\0';
      }
      else{
	die("Syntax error in fasta entry %s",line);
      }
    }
    

    strcpy(current_entry_name,seq->name);
    seq->start = 1; //new entry
  }
  else{ //old entry
    strcpy(seq->name,current_entry_name);
    seq->start = last_end_coord+1;
  }

    
  //get sequence

  boolean ready_to_return = false;
  long file_pointer = 0;
  long prev_file_pointer = 0;

  j=offset; //global length
  *full_entry=false;

  file_pointer = ftell(fp);
  

  while (ready_to_return==false){

    file_pointer = ftell(fp);

    if (fgets(line,LINE_MAX,fp) != NULL){

      int length = strlen(line);

      //sanity check
      
      if (line[0] == '>'){
	fseek(fp,file_pointer,SEEK_SET);
	ready_to_return = true;
	*full_entry = true;
      }
      else{
	if (j==max_chunk_length){
	  fseek(fp,prev_file_pointer,SEEK_SET); 
	  ready_to_return = true;
	  *full_entry = false;
	}
	else{
	  //check line is fine
	  i=0; //counter within line
	  while(i<length && ! (line[i] == '\n' || line[i] == ' ' || line[i] == '\t' || line[i] == '\r') && ready_to_return==false){
	    good_read = good_base(line[i]);	  
	    if (! good_read){
	      fprintf(stderr,"Invalid symbol [%c] pos:%i in entry %s\n",line[i],i,seq->name);

	    }
	    
	    seq->seq[j]  = line[i];
	    seq->qual[j] = '\0';
	    i++;
	    file_pointer++;
	    j++;
	    
	    if (j==max_chunk_length){	    
	      //check if line is not complete
	      if (i<length-1){
		ready_to_return = true;
		*full_entry = false;
		fseek(fp,file_pointer,SEEK_SET); 	      
	      }

	      //else 3 cases might happen with line complete
	      // 1. next line starts ">" -> *full_entry = true
	      // 2. next line is EOF -> *full_entry = true
	      // 3. next line is more sequence -> *full_entry = false
	    }
	  }	    
	}	    
      }
    }
    else{
      *full_entry = true;
      ready_to_return = true;
      //printf("complete line - complete entry - complete file\n");
    }
    
    prev_file_pointer = file_pointer;
  }
    
  seq->end = seq->start+j-1-offset;
  last_end_coord = seq->end;
  
  seq->seq[j]  = '\0';
  seq->qual[j] = '\0';

  return j;
}



/* Read sequence from file "fp" in FASTQ format.
   it returns the length of the sequence, 0 if no sequence is left in file
   read with bad chars are skipped -- this is a different behaviour to read_sequence_from_fasta
*/

int read_sequence_from_fastq(FILE *fp, Sequence * seq, int max_read_length, int fastq_ascii_offset){

  char line[LINE_MAX]; //LINE_MAX is defined in limits.h
  int i;
  int j = 0; //length of sequence
  int q = 0; //length of qualities
  long file_pointer;
  boolean good_read = true;
  int offset = fastq_ascii_offset; // this is usually 33, the offset to convert from the ascii in a fastq code, to quality value in standard Sanger format
  

  if (fp == NULL){
    die("File not defined");
  }

  if (seq == NULL){
    die("Cannot pass a NULL pointer for seq");
  }

  if (seq->seq == NULL){
    die("Dont give me a null pointer for seq->seq - alloc memory yourself and give me that");
  }

  if (seq->name == NULL){
    die("Dont give me a null pointer for seq->name - alloc memory yourself and give me that");
  }

  if (seq->qual == NULL){
    die("Dont give me a null pointer for seq->qual - alloc memory yourself and give me that");
  }

  boolean end_of_file=false;

  do{

    good_read=true;
    q=0;
    j=0;

    //read name of fastq entry
    if (fgets(line, LINE_MAX, fp) != NULL)
      {
	if (line[0] == '@')
	  {
	    for(i = 1;i<LINE_MAX;i++)
	      {	   
		if (line[i] == '\n' || line[i] == ' ' || line[i] == '\t' || line[i] == '\r')
		  {
		    break;
		  }
		if(i>LINE_MAX)
		  {
		    die("Name too long");
		  }
		seq->name[i-1] = line[i];
	      }
	    seq->name[i-1] = '\0';
	    
	    //read sequence 
	    
	    while (fgets(line, LINE_MAX, fp) != NULL)
	      {
		
		if ((line[0] == '+') || (line[0] == '-'))
		  { //go to get qualities
		    break;
		  }
		
		//check line is fine
		for(i=0;i<LINE_MAX;i++)  
		  {
		    if (line[i] == '\n' || line[i] == ' ' || line[i] == '\t' || line[i] == '\r')
		      {
			break; //fine but nothing to add
		      }
		    
		    if (! good_base(line[i]))
		    {
		      good_read=false;
		      fprintf(stderr,"Invalid symbol [%c] pos:%i in entry %s\n",line[i],i,seq->name);
		    }


		    if (j>=max_read_length){
		      fprintf(stdout,"read [%s] too long [%i]. Skip read\n",seq->name,j);
		      //exit(EXIT_FAILURE);
		      good_read=false;
		    }

		    seq->seq[j] = line[i];
		    j++;
		    
		    
		  }
	      }
	    
	    //read qualities -- verify position first
	    file_pointer = ftell(fp);
	    
		
	    while (fgets(line, LINE_MAX, fp) != NULL)
	      {

		if (line[0] == '@' && (j<=q) )//then we have gone on to the next read  
		  //allowing q>j in case where qualities longer than j
		  {
		    fseek(fp,file_pointer,SEEK_SET);
		    break; //goto next read
		  }
		
		
		for(i=0;i<LINE_MAX;i++)
		  {
		    if (line[i] == '\n' || line[i] == ' ' || line[i] == '\t' || line[i] == '\r'){
		      break; //fine but nothing to add
		    }
		    
		    seq->qual[q] = line[i]-offset;
		    q++;
		    
		    if (q==max_read_length)
		      {
			die("Qualities for [%s] longer than the max read length  [%i]. Exiting...",seq->name,q);
		      }
		    
		  }
		
		file_pointer = ftell(fp);
	      }
	    
	    
	    if (j!=q)
	      {
		fprintf(stdout,"Lengths of quality [%i] and sequence [%i]  strings don't coincide for this read: [%s]. Skip it\n",q,j,seq->name);
		good_read=false;
	      }
	    
	  }//if line starts with @
	else
	  {
	    die("syntax error in fastq file -- it misses @\n");
	  }
      }
    else
      {
	end_of_file=true;
      }

  } while (! good_read && !end_of_file); 

  seq->seq[j]   = '\0';
  seq->qual[q]  = '\0'; //this is not technically necessary but simplifies many checks downstream
  return j;
}


void alloc_sequence(Sequence * seq, int max_read_length, int max_name_length){
 
 
  if (seq == NULL){							
    die("Cannot pass a null seq to sequence_alloc");								
  }
  seq->name = malloc(sizeof(char) * max_name_length);
  if (seq->name == NULL){
    die("Out of memory trying to allocate string");
  }
  seq->seq  = malloc(sizeof(char) * (max_read_length+1));
  if (seq->seq == NULL){
    die("Out of memory trying to allocate string");
  }
  seq->qual  = malloc(sizeof(char) * (max_read_length+1));
  if (seq->qual == NULL){
    die("Out of memory trying to allocate string");
  }

  seq->seq[max_read_length]='\0';
  seq->qual[max_read_length]='\0';

}


boolean good_base(char c){
  boolean ret;
  if (c  != 'A' && c != 'a' && 
      c != 'C' && c != 'c' && 
      c != 'G' && c != 'g' && 
      c != 'T' && c != 't' && 
      c != 'N' && c != 'n' 
      ){
    ret = false;
  }	
  else{
    ret =  true;
  }
  
  return ret;
}


void free_sequence(Sequence ** sequence){
  free((*sequence)->name);
  free((*sequence)->seq);
  free((*sequence)->qual);
  free(*sequence);
  *sequence = NULL;
}

void shift_last_kmer_to_start_of_sequence(Sequence * sequence, int length, short kmer_size){

  int i;

  if (length-kmer_size<kmer_size){
    die("kmer_size too long");
  }

  for(i=0;i<kmer_size; i++){
    sequence->seq[i] = sequence->seq[length-kmer_size+i];   
    sequence->qual[i] = '\0';
  }

}
