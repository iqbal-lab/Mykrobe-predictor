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
   binary_kmer.c - routines to 
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "binary_kmer.h"
#include "event_encoding.h"

void binary_kmer_initialise_to_zero(BinaryKmer* bkmer)
{
  int i;
  for (i=0; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
    {
      ((*bkmer)[i])=(bitfield_of_64bits) 0;
    }
}


void binary_kmer_assignment_operator(BinaryKmer left, BinaryKmer right)
{
  int i;

  for (i=0; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
    {
      left[i]=right[i];
    }
}


void binary_kmer_set_all_bitfields(BinaryKmer assignee, bitfield_of_64bits val)
{
  int i;
  for (i=0; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
    {
      assignee[i]=val;
    }
}


//returns true if they are the same
boolean binary_kmer_comparison_operator(const BinaryKmer const left, const BinaryKmer const right)
{

  boolean they_are_the_same=true;

  int i;
  for (i=0; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
    {
      if (left[i]!=right[i])
	{
	  they_are_the_same=false;
	  break;                             //sorry Mario, I know you hate breaks
	}
      
    }

  return they_are_the_same;
}


//TODO - this wrongly says left<right when they are the same! IS really a <= operator
boolean binary_kmer_less_than(const BinaryKmer const left, const BinaryKmer const right, short kmer_size)
{
  boolean left_is_less_than_right=false;

  //need the following to work out which bits to ignore
  int number_of_bitfields_fully_used = kmer_size/32;
  //int number_of_bits_in_most_sig_bitfield = 2* (kmer_size-(32*number_of_bitfields_fully_used));

  int i;


  //start at most significant end
  // this would break if we had number_of_bitfields_fully_used==NUMBER_OF_BITFIELDS_IN_BINARY_KMER. But we can never have that as k is always odd.
  for (i=NUMBER_OF_BITFIELDS_IN_BINARY_KMER-number_of_bitfields_fully_used-1; i<NUMBER_OF_BITFIELDS_IN_BINARY_KMER ; i++)
    {
      if (left[i]<right[i])
	{
	  left_is_less_than_right=true;
	  break;                             
	}
      else if (left[i]>right[i])
	{
	  left_is_less_than_right=false;
	  break;
	}

      
    }

  return left_is_less_than_right;
  
}



// Implicit in this is the idea that you shift left,
// and mask to 0 the bits that fall off the left hand end
void binary_kmer_right_shift_one_base(BinaryKmer kmer)
{
  int i;
  for(i = NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1; i > 0; i--)
  {
    kmer[i] >>= 2;
    kmer[i] |= (kmer[i-1] << 62); // & 0x3
  }

  kmer[0] >>= 2;
}

void binary_kmer_left_shift_one_base(BinaryKmer kmer, short kmer_size)
{
  int top_word = NUMBER_OF_BITFIELDS_IN_BINARY_KMER - (kmer_size+31)/32;

  int i;
  for(i = top_word; i < NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1; i++)
  {
    kmer[i] <<= 2;
    kmer[i] |= (kmer[i+1] >> 62); // & 0x3
  }

  kmer[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] <<= 2;

  // Mask top word
  short top_bits = 2 * (kmer_size % 32); // bits in top word
  kmer[top_word] &= (~ (bitfield_of_64bits)0 >> (64 - top_bits));
}


void binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(BinaryKmer* bkmer, Nucleotide n, short kmer_size)
{

  //shift left by one base,
  binary_kmer_left_shift_one_base(*bkmer, kmer_size);

  // add new base at right hand end 
  (*bkmer)[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] |= n;

}




char reverse_char_nucleotide(char c)
{
  switch (c)
    {
    case 'A':
      return 'T';
    case 'C':
      return 'G';
    case 'G':
      return 'C';
    case 'T':
      return 'A';
    case 'a':
      return 't';
    case 'c':
      return 'g';
    case 'g':
      return 'c';
    case 't':
      return 'a';
    default:
      die("Non-existent nucleotide %c\n",c);
    }
}


//length is the length in number of bases; the char* should have one MORE base than that allocated, to hold '\0'
char * seq_reverse_complement(char * in, int length, char * out){

  int k;
  for(k=0;k<length;k++){
    out[k] = reverse_char_nucleotide(in[length-k-1]);
  }
  out[length] = '\0';
  return out;
}


/*
Nucleotide reverse_binary_nucleotide(Nucleotide n)
{
  switch (n)
    {
    case Adenine:
      return Thymine;
    case Cytosine:
      return Guanine;
    case Guanine:
      return Cytosine;
    case Thymine:
      return Adenine;
    default:
      die("Calling reverse_binary_nucleotide on non-existent nucleotide %i",n);
    }
}
*/


char * nucleotides_to_string(Nucleotide * nucleotides, int length, char * string){
  
  if (string == NULL){
    die("Seq argument cannot be NULL");
  }

  int i;
  for(i=0;i<length;i++){
    string[i] = binary_nucleotide_to_char(nucleotides[i]);
  }
  
  string[length] = '\0';
  return string;
}



//The first argument - seq - is a C string in A,C,G,T format
//The second argument - quality - is a string of qualities for the sequence, one byte per base.
//quality cutoff argument defines the threshold for quality
//return total number of kmers read
//The third argument - length - is the length in bases of the sequence.
//Final argument - homopolymer_cutoff - allows you to break a sliding window at the end of a homopolymer.
//                 If homopolymer_cutoff= n > 0, then as soon as the latest base is the n-th in a homopolymeric sequence, the window is broken, and the next
//                  window only starts when there is  a new base.
//return total number of kmers read in - ie good kmers that go into windows
int get_sliding_windows_from_sequence(char * sequence,  char * qualities, int length, char quality_cut_off, short kmer_size, KmerSlidingWindowSet * windows, 
				      int max_windows, int max_kmers, boolean break_homopolymers, int homopolymer_cutoff){  

  char first_kmer[kmer_size+1];
  first_kmer[kmer_size]='\0';

  int i=0; //current index
  int count_kmers = 0;

  if (sequence == NULL){
    die("In get_sliding_windows_from_sequence, sequence is NULL");
  }

  if (length < kmer_size || max_windows == 0 || max_kmers == 0){
    return 0;
  }

  int index_windows = 0;
  
  //loop over the bases in the sequence
  //index i is the current position in input sequence -- it nevers decreases. 
  
  int hom_ct; //count how long current homopolymer run is

  do{

    //built first kmer, ie a stretch of kmer_size good qualities bases
    int j = 0; //count how many good bases
    hom_ct=0; 
    while ((i<length) && (j<kmer_size)){

      //collects the bases in the first kmer
      first_kmer[j] = sequence[i];

      //what is current homopolymer length
      if ( (j>0) && (first_kmer[j]==first_kmer[j-1]) )
	{
	  hom_ct++;
	}
      else if (j>=0)
	{
	  hom_ct=1;
	}


      if ((char_to_binary_nucleotide(sequence[i]) == Undefined) || 
	  (quality_cut_off>0 && qualities[i]<= quality_cut_off)){
	j=0; //restart the first kmer 
      }
      else if ( (break_homopolymers==true) && (hom_ct>=homopolymer_cutoff) )
	{
	  
	  //now we may be in the middle of a very long hompoppolymer run.So we want to increment i sufficiently to hit the next base
	  int first_base_after_homopolymer=i;
	  while ( (first_base_after_homopolymer<length) && (sequence[first_base_after_homopolymer]==first_kmer[j]) )
	    {
	      first_base_after_homopolymer++;
	    }

	  i=first_base_after_homopolymer-1; //we are going to add one at the end of the loop, just below
	  j=0; //restart the first kmer

	}
      
      else{
	j++;
      }

      i++; 
    }

    if (j==kmer_size){ //ie we did not parse the entire sequence looking for a single good kmer, the first kmer
      
      count_kmers++;

      //new sliding window
      if (index_windows>=max_windows){
    die("number of windows is bigger than max_windows in "
        "get_sliding_windows_from_sequence");
	 }

      KmerSlidingWindow * current_window =&(windows->window[index_windows]);

      int index_kmers = 0;
      //do first kmer
      BinaryKmer tmp_bin_kmer;
      seq_to_binary_kmer(first_kmer,kmer_size, &tmp_bin_kmer);
      binary_kmer_assignment_operator(current_window->kmer[index_kmers] , tmp_bin_kmer);

      //do the rest --
      index_kmers++;
    
      while(i<length){
	
	if (index_kmers>=max_kmers){
	  die("Number of kmers is bigger than max_kmers in "
        "get_sliding_windows_from_sequence - second check");
	}

	Nucleotide current_base = char_to_binary_nucleotide(sequence[i]);

	if ( (i>0) && (sequence[i] == sequence[i-1]) )
	  {
	    hom_ct++;
	  }
	else 
	  {
	    hom_ct=1;
	  }

	if ((current_base == Undefined) ||
	    (quality_cut_off!=0 && qualities[i]<= quality_cut_off)){
	  i++; 
	  break;
	}
	else if ( (break_homopolymers==true) && (hom_ct>=homopolymer_cutoff) )
	  {
	    //now we may be in the middle of a very long homopolymer run.So we want to increment i sufficiently to go beyond         
	    int first_base_after_homopolymer=i;
	    while ( (first_base_after_homopolymer<length) && (sequence[first_base_after_homopolymer]==sequence[i]) )
	      {
		first_base_after_homopolymer++;
	      }

	    i=first_base_after_homopolymer; 
	    break;
	  }


	//set the kmer to previous
	binary_kmer_assignment_operator(current_window->kmer[index_kmers], current_window->kmer[index_kmers-1]);
	binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&(current_window->kmer[index_kmers]), current_base, kmer_size);

	index_kmers++;
	count_kmers++;
	i++;
      }

      current_window->nkmers = index_kmers; 
    
      index_windows++;
            
    }
  } while (i<length);
 
  windows->nwindows = index_windows;

  return count_kmers;
}






//caller passes in preallocated BinaryKmer, which is also returned in the return value
BinaryKmer* seq_to_binary_kmer(char * seq, short kmer_size, BinaryKmer* prealloced_kmer){

  //sanity checks
  if (seq==NULL)
    {
      die("Do not passs null ptr to seq_to_binary_kmer. Exiting..");
    }
  if (strlen(seq) != (unsigned)kmer_size)
    {
      die("Calling seq_to_binary_kmer with  a sequence %s of length %d, "
          "but kmer size %d, which is different. Exiting",
          seq, (int) strlen(seq),  kmer_size);
    }
  
  int j;
  binary_kmer_initialise_to_zero(prealloced_kmer);

  for(j=0;j<kmer_size;j++){

    if (char_to_binary_nucleotide(seq[j]) == Undefined){
      prealloced_kmer=NULL;
      return prealloced_kmer;
    }


    binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(prealloced_kmer, char_to_binary_nucleotide(seq[j]), kmer_size ); 
    
  }

  return prealloced_kmer;

}


  
  



//caller passes in allocated char*. This is returned and also set in 3rd argument.
//user of this method is responsible for deallocating the returned sequence
//note that the allocated space has to be kmer_size+1;
char * binary_kmer_to_seq(BinaryKmer* bkmer, short kmer_size, char * seq){
 
  BinaryKmer local_bkmer;
  binary_kmer_assignment_operator(local_bkmer, *bkmer);

  if (seq == NULL){
      die("Seq argument cannot be NULL");
    }

  int mask = 3; // 0000011 mask used to extract the two least significative bits
  int j;
  
  for(j=kmer_size-1; j>=0; j--){ //start from the back of the sequence
    
    //get translation for the two least significant bits
    seq[j] =  binary_nucleotide_to_char(local_bkmer[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] & mask);
    //shift right zam
    binary_kmer_right_shift_one_base(local_bkmer); //note this is a local copy internal to this function - not altering the original BinaryKmer

  }
  
  seq[kmer_size] = '\0';
  
  return seq;
}

// kmer and prealloc_reverse_kmer may point to the same address
// This is a highly optimised version of the above function
BinaryKmer* binary_kmer_reverse_complement(BinaryKmer* kmer, short kmer_size,
                                           BinaryKmer* prealloc_reverse_kmer)
{
  binary_kmer_initialise_to_zero(prealloc_reverse_kmer);

  BinaryKmer kmer_copy;

  // Copy
  memcpy(kmer_copy, kmer, BINARY_KMER_BYTES);

  // Complement
  int i;
  for(i = 0; i < NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
  {
    kmer_copy[i] = ~kmer_copy[i];
  }

  // Reverse
  // Loop over full bitfields
  int j, k;
  for(i = NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1; i > 0; i--)
  {
    for(j = 0; j < 32; j++)
    {
      // Shift destination left
      for(k = 0; k < NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1; k++)
      {
        (*prealloc_reverse_kmer)[k] <<= 2;
        (*prealloc_reverse_kmer)[k] |= ((*prealloc_reverse_kmer)[k+1] >> 62);
      }

      (*prealloc_reverse_kmer)[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] <<= 2;

      // Append new base
      (*prealloc_reverse_kmer)[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] |= kmer_copy[i] & 0x3;

      // Shift source right
      kmer_copy[i] >>= 2;
    }
  }

  // Do remaining bases in last bitfield [0]
  int top_bases = kmer_size % 32;

  for(i = 0; i < top_bases; i++)
  {
    // Shift destination left
    for(k = 0; k < NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1; k++)
    {
      (*prealloc_reverse_kmer)[k] <<= 2;
      (*prealloc_reverse_kmer)[k] |= ((*prealloc_reverse_kmer)[k+1] >> 62);
    }

    (*prealloc_reverse_kmer)[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] <<= 2;

    // Append new base
    (*prealloc_reverse_kmer)[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] |= kmer_copy[0] & 0x3;
    
    // Shift source right
    kmer_copy[0] >>= 2;
  }

  // Mask top word (we didn't zero at the beginning!)
  short top_bits = 2 * top_bases; // bits in top word
  (*prealloc_reverse_kmer)[0] &= (~(uint64_t)0 >> (64 - top_bits));

  return prealloc_reverse_kmer;
}

Nucleotide binary_kmer_get_first_nucleotide(BinaryKmer* kmer, short kmer_size)
{
  int number_of_bits_in_most_sig_bitfield = 2 * (kmer_size % 32);

  return ((*kmer)[0] >> (number_of_bits_in_most_sig_bitfield - 2)) & 0x3;
}

Nucleotide binary_kmer_get_last_nucleotide(BinaryKmer* kmer)
{
  return (*kmer)[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] & 0x3;
}


void binary_kmer_alloc_kmers_set(KmerSlidingWindowSet * windows, int max_windows, int max_kmers){

  if (windows == NULL){
    die("Cannot pass a NULL window to alloc");
  } 
  
  //allocate memory for the sliding windows         
  windows->max_nwindows = max_windows;

  windows->window = malloc(sizeof(KmerSlidingWindow) * max_windows);       
  if (windows->window== NULL){
    die("Out of memory trying to allocate an array of KmerSlidingWindow");
  }
  windows->nwindows = 0;
  
  //allocate memory for every every sliding window
  int w;
  for(w=0;w<max_windows;w++){
    windows->window[w].nkmers=0;
    //windows->window[w].kmer = malloc(sizeof(BinaryKmer) * max_kmers);
    windows->window[w].kmer = (BinaryKmer*) malloc(sizeof(bitfield_of_64bits)*NUMBER_OF_BITFIELDS_IN_BINARY_KMER * max_kmers);
    if (windows->window[w].kmer == NULL){
      die("binary_kmer: Out of memory trying to allocate an array of BinaryKmer");
    }      
  }      
  
}

void binary_kmer_free_kmers(KmerSlidingWindow * * kmers)
{

  free((*kmers)->kmer);
  free(*kmers); 
  *kmers = NULL;
}

void binary_kmer_free_kmers_set(KmerSlidingWindowSet * * kmers_set)
{
  int w;
  for(w=0;w<(*kmers_set)->max_nwindows;w++){
    free((*kmers_set)->window[w].kmer);
  }
  
  free((*kmers_set)->window);
  free(*kmers_set);
  *kmers_set = NULL;
}

void nucleotide_iterator(void (*f)(Nucleotide)){
  
  int i;
  for (i=0;i<4;i++){
    f(i);
  }
  
}

void nucleotide_iterator_orientation(void (*f)(Nucleotide n,Orientation o)){
  
  int i;
  for (i=0;i<4;i++){
    f(i,forward);
    f(i,reverse);
  }
  
}
