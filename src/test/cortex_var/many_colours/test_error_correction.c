/*
 * Copyright 2009-2013 Zamin Iqbal 
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
  test_error_correction.c
*/

#include <stdlib.h>

#include <CUnit.h>
#include <Basic.h>

#include "string_buffer.h"
#include "error_correction.h"
#include "file_reader.h"

void test_base_mutator()
{
  char* str="AACGT";
  StrBuf* strbuf=strbuf_new();
  strbuf_set(strbuf, str);
  CU_ASSERT(mutate_base(strbuf, 1, 0)==true);
  CU_ASSERT(strcmp(strbuf->buff, "ACCGT")==0);
  strbuf_set(strbuf, str);
  CU_ASSERT(mutate_base(strbuf, 1, 1)==true);
  CU_ASSERT(strcmp(strbuf->buff, "AGCGT")==0);
  strbuf_set(strbuf, str);
  CU_ASSERT(mutate_base(strbuf, 1, 2)==true);
  CU_ASSERT(strcmp(strbuf->buff, "ATCGT")==0);
  strbuf_set(strbuf, str);
  CU_ASSERT(mutate_base(strbuf, 1, 3)==false);

  strbuf_set(strbuf, str);
  CU_ASSERT(mutate_base(strbuf, 2, 0)==true);
  CU_ASSERT(strcmp(strbuf->buff, "AAAGT")==0);
  strbuf_set(strbuf, str);
  CU_ASSERT(mutate_base(strbuf, 2, 1)==true);
  CU_ASSERT(strcmp(strbuf->buff, "AAGGT")==0);
  strbuf_set(strbuf, str);
  CU_ASSERT(mutate_base(strbuf, 2, 2)==true);
  CU_ASSERT(strcmp(strbuf->buff, "AATGT")==0);
  strbuf_set(strbuf, str);
  CU_ASSERT(mutate_base(strbuf, 2, 3)==false);

  char* str2 = "NNAC";
  strbuf_set(strbuf, str2);
  CU_ASSERT(mutate_base(strbuf, 0, 0)==true);
  CU_ASSERT(strcmp(strbuf->buff, "ANAC")==0);
  strbuf_set(strbuf, str2);
  CU_ASSERT(mutate_base(strbuf, 0, 1)==true);
  CU_ASSERT(strcmp(strbuf->buff, "CNAC")==0);
  strbuf_set(strbuf, str2);
  CU_ASSERT(mutate_base(strbuf, 0, 2)==true);
  CU_ASSERT(strcmp(strbuf->buff, "GNAC")==0);
  strbuf_set(strbuf, str2);
  CU_ASSERT(mutate_base(strbuf, 0, 3)==true);
  CU_ASSERT(strcmp(strbuf->buff, "TNAC")==0);

  strbuf_free(strbuf);
}



void test_fix_end_if_unambiguous()
{

  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  //first set up the hash/graph
  int kmer_size = 5;
  int number_of_bits = 10;
  int bucket_size = 10;

  //*********************************************
  // Test 1 - does it fix a kmer to match the graph 
  //          when there is a unique way of doing it
  //*********************************************

  dBGraph * db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);
  
  //Load the following fasta: this is to be the trusted graph - just AAAAA
  // >
  // AAAAA


  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/error_correction/graph1.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);


  StrBuf* readbuf = strbuf_new();
  StrBuf* qualbuf = strbuf_new();
  StrBuf* kmerbuf = strbuf_init(5);
  char kmerstr[6];
  char* read1="TGACAAAAAC"; 
  char* qual1="IIIIIIIII#";
  char cutoff = 40+33;
  strbuf_set(readbuf, read1);
  strbuf_set(qualbuf, qual1);
  //at position 5 in the read we have AAAAC, which is fixable to AAAAA
  CU_ASSERT(fix_end_if_unambiguous(Right, readbuf, qualbuf, cutoff, 5, kmerbuf, kmerstr, db_graph)==true);
  CU_ASSERT(strcmp(readbuf->buff, "TGACAAAAAA")==0);
  CU_ASSERT(strcmp(qualbuf->buff, "IIIIIIIIIJ")==0);
  strbuf_set(readbuf, read1);//reset
  strbuf_set(qualbuf, qual1);

  CU_ASSERT(fix_end_if_unambiguous(Left, readbuf, qualbuf, cutoff, 5, kmerbuf, kmerstr, db_graph)==false);
  CU_ASSERT(strcmp(readbuf->buff, "TGACAAAAAC")==0);
  CU_ASSERT(strcmp(qualbuf->buff, "IIIIIIIII#")==0);

  hash_table_free(&db_graph);
  strbuf_free(readbuf);
  strbuf_free(qualbuf);
  strbuf_free(kmerbuf);




  //*********************************************
  // Test 2 - confirm it does not fix
  //          when there are two mutations which
  //          would make it in the graph
  //*********************************************

  db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);
  
  //Load the following fasta: this is to be the trusted graph. note ACGGA can be corrected in TWO ways to matc this graph
  // >
  // ACGGCTTTACGGT


  load_se_filelist_into_graph_colour(
    "../data/test/error_correction/graph2.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);


  readbuf = strbuf_new();
  qualbuf = strbuf_new();
  kmerbuf = strbuf_init(5);

  char* read2="TGACACGGAGGTACGT";//contains ACGGA 
  strbuf_set(readbuf, read2);
  char* qual2="IIIIIIII#IIIIIII";
  strbuf_set(qualbuf, qual2);

  //at position 4 in the read we have ACGGA which is fixable to ACGGC OR ACGGT
  CU_ASSERT(fix_end_if_unambiguous(Right, readbuf, qualbuf, cutoff, 4, kmerbuf, kmerstr, db_graph)==false);
  CU_ASSERT(strcmp(readbuf->buff, "TGACACGGAGGTACGT")==0);
  CU_ASSERT(strcmp(qualbuf->buff, "IIIIIIII#IIIIIII")==0);

  strbuf_set(readbuf, read2);//reset
  strbuf_set(qualbuf, qual2);//reset


  //in fact, now confirm nothing is changed anywhere in the read
  int j;
  for (j=0; j<9; j++)
    {
      CU_ASSERT(fix_end_if_unambiguous(Right, readbuf, qualbuf, cutoff, j, kmerbuf, kmerstr, db_graph)==false);
      CU_ASSERT(strcmp(readbuf->buff, "TGACACGGAGGTACGT")==0);
      CU_ASSERT(strcmp(qualbuf->buff, "IIIIIIII#IIIIIII")==0);
      strbuf_set(readbuf, read2);//reset
      strbuf_set(qualbuf, qual2);//reset
    }


  hash_table_free(&db_graph);
  strbuf_free(readbuf);
  strbuf_free(kmerbuf);
  strbuf_free(qualbuf);

  
}



void test_get_first_good_kmer_and_populate_qual_array()
{
  //pass in a read with chosen qualities, and check 
  //we correctly get int arrays describing which bases have good quals
  //and which kmers are in the graph
  //this function makes a performance choice, no abort and not bother
  //populating the kmer array as soon as we find out all the qualities are high,
  //as we immediately know we are not going to do any error-correction.


  //first set up the hash/graph
  int kmer_size = 5;
  int number_of_bits = 10;
  int bucket_size = 10;

  //*********************************************

  //*********************************************

  dBGraph * db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);
  
  //Load the following fasta: this is to be the trusted graph
  // >
  // ACGGCTTTACGGT

  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/error_correction/graph3.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  //read a small FASTQ - all these kmers ARE int he graph and have high quality
  // @zam all qual 30
  // ACGGCTTTACGGT
  // +
  // ?????????????


  SeqFile *sf = seq_file_open("../data/test/error_correction/fq_for_comparing_with_graph3.fq");
  if(sf == NULL)
    {
      // Error opening file
      fprintf(stderr, "Error: cannot read seq file '%s'\n", "../data/test/error_correction/fq_for_comparing_with_graph3.fq");
      exit(EXIT_FAILURE);
    }

  seq_next_read(sf);
  StrBuf* read_seq  = strbuf_new();
  StrBuf* read_qual = strbuf_new();
  seq_read_all_bases(sf, read_seq);
  seq_read_all_quals(sf, read_qual);
  int read_len = seq_get_length(sf);
  seq_file_close(sf);

  
  //now let's check
  int quals_good[read_len];
  set_int_array(quals_good, read_len, 1);


  char qual_cutoff=20+33;//33=ascii FASTQ sanger offset
  int first_good_kmer=-1;

  Orientation strand_first_good_kmer=reverse;//just testing
  char* readid="zam";
  ReadCorrectionDecison result = get_first_good_kmer_and_populate_qual_array(readid,read_seq, read_qual,
									     read_len-5+1, read_len,
									     quals_good,
									     qual_cutoff, &first_good_kmer,&strand_first_good_kmer,
									     db_graph, DontWorryAboutLowQualBaseUnCorrectable, false);
  CU_ASSERT(result==PrintUncorrected);
  CU_ASSERT(strand_first_good_kmer==forward);//I've not marked the ref in the graph, so everything treated as forward
  int j;
  for (j=0; j<read_len; j++)
    {
      CU_ASSERT(quals_good[j]==1);
    }
  CU_ASSERT(first_good_kmer==-1);

  //result should be independent of policy
  set_int_array(quals_good, read_len, 1);
  first_good_kmer=-1;

  strand_first_good_kmer=reverse;//just testing
  result = get_first_good_kmer_and_populate_qual_array(readid,read_seq, read_qual,
						       read_len-5+1, read_len,
						       quals_good,
						       qual_cutoff, &first_good_kmer, &strand_first_good_kmer,
						       db_graph, DiscardReadIfLowQualBaseUnCorrectable, false);
  CU_ASSERT(result==PrintUncorrected);
  CU_ASSERT(strand_first_good_kmer==forward);
  for (j=0; j<read_len; j++)
    {
      CU_ASSERT(quals_good[j]==1);
    }
  CU_ASSERT(first_good_kmer==-1);




  //now try again, but change the quality threshold so that all the qualities fail

  qual_cutoff=90+33;
  first_good_kmer=-1;
  set_int_array(quals_good, read_len, 1);
  strand_first_good_kmer=reverse;//just testing
  result = get_first_good_kmer_and_populate_qual_array(readid,read_seq, read_qual,
						       read_len-5+1, read_len,
						       quals_good,
						       qual_cutoff, &first_good_kmer, &strand_first_good_kmer,
						       db_graph, DontWorryAboutLowQualBaseUnCorrectable, false);
  CU_ASSERT(result==PrintCorrected); //this function does not know all kmers in the graph
  CU_ASSERT(strand_first_good_kmer==forward);
  for (j=0; j<read_len; j++)
    {
      CU_ASSERT(quals_good[j]==0);
    }
  CU_ASSERT(first_good_kmer==0);


  first_good_kmer=-1;
  set_int_array(quals_good, read_len, 1);
  strand_first_good_kmer=reverse;//just testing
  result = get_first_good_kmer_and_populate_qual_array(readid,read_seq, read_qual,
						       read_len-5+1, read_len,
						       quals_good,
						       qual_cutoff, &first_good_kmer,&strand_first_good_kmer,
						       db_graph, DiscardReadIfLowQualBaseUnCorrectable, false);
  CU_ASSERT(result==PrintCorrected); 
  CU_ASSERT(strand_first_good_kmer==forward);
  for (j=0; j<read_len; j++)
    {
      CU_ASSERT(quals_good[j]==0);
    }
  CU_ASSERT(first_good_kmer==0);




  //now try a more complicated read, with a kmer that is not in the graph, but HIGH quality
  //  @zam all qual 10 except30 at the base which isnot in the graph
  //  ACGGCTTGACGGT
  //  +
  //  +++++++?+++++

  //so at position 3 in the read we have GCTTG which is not in the graph, could be fixed to be in the graph,
  // BUT IS HIGH QUALITY, so should not be fixed. Upstream functions check that. This function just checks for presence.

  sf = seq_file_open("../data/test/error_correction/fq2_for_comparing_with_graph3.fq");
  if(sf == NULL)
    {
      // Error opening file
      fprintf(stderr, "Error: cannot read seq file '%s'\n", "../data/test/error_correction/fq2_for_comparing_with_graph3.fq");
      exit(EXIT_FAILURE);
    }

  seq_next_read(sf);
  seq_read_all_bases(sf, read_seq);
  seq_read_all_quals(sf, read_qual);
  read_len = seq_get_length(sf);
  seq_file_close(sf);

  first_good_kmer=-1;
  
  set_int_array(quals_good, read_len, 1);
  qual_cutoff=20+33;
  strand_first_good_kmer=reverse;//just testing
  result = get_first_good_kmer_and_populate_qual_array(readid,read_seq, read_qual,
						       read_len-5+1, read_len,
						       quals_good,
						       qual_cutoff, &first_good_kmer, &strand_first_good_kmer,
						       db_graph, DiscardReadIfLowQualBaseUnCorrectable, false);


  //NOTE - the following is correct behaviour.
  //this function just decides that it is worth trying to correct
  //as there is a kmer not in the graph. It doesn't check to see if that kmer has low quality
  CU_ASSERT(result==PrintCorrected);
  CU_ASSERT(strand_first_good_kmer==forward);

  CU_ASSERT(quals_good[0]==0);
  CU_ASSERT(quals_good[1]==0);
  CU_ASSERT(quals_good[2]==0);
  CU_ASSERT(quals_good[3]==0);
  CU_ASSERT(quals_good[4]==0);
  CU_ASSERT(quals_good[5]==0);
  CU_ASSERT(quals_good[6]==0);
  CU_ASSERT(quals_good[7]==1);
  CU_ASSERT(quals_good[8]==0);
  CU_ASSERT(quals_good[9]==0);
  CU_ASSERT(quals_good[10]==0);
  CU_ASSERT(quals_good[11]==0);
  CU_ASSERT(quals_good[12]==0);

  CU_ASSERT(first_good_kmer==0);

  //and agin changing policy - should get same answer

  set_int_array(quals_good, read_len, 1);
  qual_cutoff=20+33;
  strand_first_good_kmer=reverse;//just testing
  result = get_first_good_kmer_and_populate_qual_array(readid,read_seq, read_qual,
						       read_len-5+1, read_len,
						       quals_good,
						       qual_cutoff, &first_good_kmer,&strand_first_good_kmer,
						       db_graph, DontWorryAboutLowQualBaseUnCorrectable,
						       false);


  CU_ASSERT(result==PrintCorrected);

  CU_ASSERT(quals_good[0]==0);
  CU_ASSERT(quals_good[1]==0);
  CU_ASSERT(quals_good[2]==0);
  CU_ASSERT(quals_good[3]==0);
  CU_ASSERT(quals_good[4]==0);
  CU_ASSERT(quals_good[5]==0);
  CU_ASSERT(quals_good[6]==0);
  CU_ASSERT(quals_good[7]==1);
  CU_ASSERT(quals_good[8]==0);
  CU_ASSERT(quals_good[9]==0);
  CU_ASSERT(quals_good[10]==0);
  CU_ASSERT(quals_good[11]==0);
  CU_ASSERT(quals_good[12]==0);

  CU_ASSERT(first_good_kmer==0);







  //now shift the quality threshold to just above 30, so all the bases have low qual
  first_good_kmer=-1;
  
  set_int_array(quals_good, read_len, 1);
  qual_cutoff=31+33;
  strand_first_good_kmer=reverse;//just testing
  result = get_first_good_kmer_and_populate_qual_array(readid,read_seq, read_qual,
						       read_len-5+1, read_len,
						       quals_good,
						       qual_cutoff, &first_good_kmer,&strand_first_good_kmer,
						       db_graph, DontWorryAboutLowQualBaseUnCorrectable,
						       false);


  CU_ASSERT(result==PrintCorrected);

  CU_ASSERT(quals_good[0]==0);
  CU_ASSERT(quals_good[1]==0);
  CU_ASSERT(quals_good[2]==0);
  CU_ASSERT(quals_good[3]==0);
  CU_ASSERT(quals_good[4]==0);
  CU_ASSERT(quals_good[5]==0);
  CU_ASSERT(quals_good[6]==0);
  CU_ASSERT(quals_good[7]==0);
  CU_ASSERT(quals_good[8]==0);
  CU_ASSERT(quals_good[9]==0);
  CU_ASSERT(quals_good[10]==0);
  CU_ASSERT(quals_good[11]==0);
  CU_ASSERT(quals_good[12]==0);

  CU_ASSERT(first_good_kmer==0);



  //and again with other policy
  first_good_kmer=-1;
  
  set_int_array(quals_good, read_len, 1);
  qual_cutoff=31+33;
  strand_first_good_kmer=reverse;//just testing
  result = get_first_good_kmer_and_populate_qual_array(readid,read_seq, read_qual,
						       read_len-5+1, read_len,
						       quals_good,
						       qual_cutoff, &first_good_kmer,&strand_first_good_kmer,
						       db_graph, DiscardReadIfLowQualBaseUnCorrectable,
						       false);
  CU_ASSERT(result==PrintCorrected);

  CU_ASSERT(quals_good[0]==0);
  CU_ASSERT(quals_good[1]==0);
  CU_ASSERT(quals_good[2]==0);
  CU_ASSERT(quals_good[3]==0);
  CU_ASSERT(quals_good[4]==0);
  CU_ASSERT(quals_good[5]==0);
  CU_ASSERT(quals_good[6]==0);
  CU_ASSERT(quals_good[7]==0);
  CU_ASSERT(quals_good[8]==0);
  CU_ASSERT(quals_good[9]==0);
  CU_ASSERT(quals_good[10]==0);
  CU_ASSERT(quals_good[11]==0);
  CU_ASSERT(quals_good[12]==0);

  CU_ASSERT(first_good_kmer==0);





  //now shift the quality threshold to just below 10, so all bases have high qual
  //so now it doesn't matter about the kmers, whether they are in or out, you definitely
  //print the read uncorrected.
  
  set_int_array(quals_good, read_len, 1);
  qual_cutoff=9+33;
  first_good_kmer=-1;
  strand_first_good_kmer=reverse;//just testing
  result = get_first_good_kmer_and_populate_qual_array(readid,read_seq, read_qual,
						       read_len-5+1, read_len,
						       quals_good,
						       qual_cutoff, &first_good_kmer,&strand_first_good_kmer,
						       db_graph, DiscardReadIfLowQualBaseUnCorrectable,
						       false);


  CU_ASSERT(result==PrintUncorrected);


  CU_ASSERT(quals_good[0]==1);
  CU_ASSERT(quals_good[1]==1);
  CU_ASSERT(quals_good[2]==1);
  CU_ASSERT(quals_good[3]==1);
  CU_ASSERT(quals_good[4]==1);
  CU_ASSERT(quals_good[5]==1);
  CU_ASSERT(quals_good[6]==1);
  CU_ASSERT(quals_good[7]==1);
  CU_ASSERT(quals_good[8]==1);
  CU_ASSERT(quals_good[9]==1);
  CU_ASSERT(quals_good[10]==1);
  CU_ASSERT(quals_good[11]==1);
  CU_ASSERT(quals_good[12]==1);

  //this isnot set- does not bother to check, since all quals high
  CU_ASSERT(first_good_kmer==-1);


  //now again with other policy
  set_int_array(quals_good, read_len, 1);
  qual_cutoff=9+33;
  first_good_kmer=-1;
  strand_first_good_kmer=reverse;//just testing
  result = get_first_good_kmer_and_populate_qual_array(readid,read_seq, read_qual,
						       read_len-5+1, read_len,
						       quals_good,
						       qual_cutoff, &first_good_kmer,&strand_first_good_kmer,
						       db_graph, DontWorryAboutLowQualBaseUnCorrectable,
						       false);


  CU_ASSERT(result==PrintUncorrected);


  CU_ASSERT(quals_good[0]==1);
  CU_ASSERT(quals_good[1]==1);
  CU_ASSERT(quals_good[2]==1);
  CU_ASSERT(quals_good[3]==1);
  CU_ASSERT(quals_good[4]==1);
  CU_ASSERT(quals_good[5]==1);
  CU_ASSERT(quals_good[6]==1);
  CU_ASSERT(quals_good[7]==1);
  CU_ASSERT(quals_good[8]==1);
  CU_ASSERT(quals_good[9]==1);
  CU_ASSERT(quals_good[10]==1);
  CU_ASSERT(quals_good[11]==1);
  CU_ASSERT(quals_good[12]==1);

  //this isnot set- does not bother to check, since all quals high
  CU_ASSERT(first_good_kmer==-1);



  //now try a read where all kmers are not in the graph, but high quality
  //  @zam all qual 60, no kmers in graph
  //  ATATATATATAT
  //  +
  //  ]]]]]]]]]]]]

  //this read should be printed uncorrected, as all high quality

  sf = seq_file_open("../data/test/error_correction/fq3_for_comparing_with_graph3.fq");
  if(sf == NULL)
    {
      // Error opening file
      fprintf(stderr, "Error: cannot read seq file '%s'\n", "../data/test/error_correction/fq3_for_comparing_with_graph3.fq");
      exit(EXIT_FAILURE);
    }

  seq_next_read(sf);
  seq_read_all_bases(sf, read_seq);
  seq_read_all_quals(sf, read_qual);
  read_len = seq_get_length(sf);
  seq_file_close(sf);

  first_good_kmer=-1;  
  set_int_array(quals_good, read_len, 1);
  qual_cutoff=10+33;
  strand_first_good_kmer=reverse;//just testing
  result = get_first_good_kmer_and_populate_qual_array(readid,read_seq, read_qual,
						       read_len-5+1, read_len,
						       quals_good,
						       qual_cutoff, &first_good_kmer,&strand_first_good_kmer,
						       db_graph, DontWorryAboutLowQualBaseUnCorrectable, false);

  CU_ASSERT(result==PrintUncorrected);

  CU_ASSERT(quals_good[0]==1);
  CU_ASSERT(quals_good[1]==1);
  CU_ASSERT(quals_good[2]==1);
  CU_ASSERT(quals_good[3]==1);
  CU_ASSERT(quals_good[4]==1);
  CU_ASSERT(quals_good[5]==1);
  CU_ASSERT(quals_good[6]==1);
  CU_ASSERT(quals_good[7]==1);
  CU_ASSERT(quals_good[8]==1);
  CU_ASSERT(quals_good[9]==1);
  CU_ASSERT(quals_good[10]==1);
  CU_ASSERT(quals_good[11]==1);

  //not set
  CU_ASSERT(first_good_kmer==-1);


  //policy makes no difference:
  first_good_kmer=-1;  
  set_int_array(quals_good, read_len, 1);
  qual_cutoff=10+33;
  strand_first_good_kmer=reverse;//just testing
  result = get_first_good_kmer_and_populate_qual_array(readid,read_seq, read_qual,
						       read_len-5+1, read_len,
						       quals_good,
						       qual_cutoff, &first_good_kmer, &strand_first_good_kmer,
						       db_graph, DiscardReadIfLowQualBaseUnCorrectable, false);

  CU_ASSERT(result==PrintUncorrected);

  CU_ASSERT(quals_good[0]==1);
  CU_ASSERT(quals_good[1]==1);
  CU_ASSERT(quals_good[2]==1);
  CU_ASSERT(quals_good[3]==1);
  CU_ASSERT(quals_good[4]==1);
  CU_ASSERT(quals_good[5]==1);
  CU_ASSERT(quals_good[6]==1);
  CU_ASSERT(quals_good[7]==1);
  CU_ASSERT(quals_good[8]==1);
  CU_ASSERT(quals_good[9]==1);
  CU_ASSERT(quals_good[10]==1);
  CU_ASSERT(quals_good[11]==1);

  //not set
  CU_ASSERT(first_good_kmer==-1);





  //now try a read where all kmers are not in the graph, and NOT all high quality
  // @zam all qual 60, except one with low quality 3, but no kmers in graph
  // ATATATATATAT
  //+
  //]]]]]]]]]]]#




  sf = seq_file_open("../data/test/error_correction/fq4_for_comparing_with_graph3.fq");
  if(sf == NULL)
    {
      // Error opening file
      fprintf(stderr, "Error: cannot read seq file '%s'\n", "../data/test/error_correction/fq4_for_comparing_with_graph3.fq");
      exit(EXIT_FAILURE);
    }

  seq_next_read(sf);
  seq_read_all_bases(sf, read_seq);
  seq_read_all_quals(sf, read_qual);
  read_len = seq_get_length(sf);
  seq_file_close(sf);

  first_good_kmer=-1;  
  set_int_array(quals_good, read_len, 1);
  qual_cutoff=10+33;
  strand_first_good_kmer=reverse;//just testing
  result = get_first_good_kmer_and_populate_qual_array(readid,read_seq, read_qual,
						       read_len-5+1, read_len,
						       quals_good,
						       qual_cutoff, &first_good_kmer, &strand_first_good_kmer,
						       db_graph, DontWorryAboutLowQualBaseUnCorrectable, false);



  CU_ASSERT(result==PrintUncorrected);

  CU_ASSERT(quals_good[0]==1);
  CU_ASSERT(quals_good[1]==1);
  CU_ASSERT(quals_good[2]==1);
  CU_ASSERT(quals_good[3]==1);
  CU_ASSERT(quals_good[4]==1);
  CU_ASSERT(quals_good[5]==1);
  CU_ASSERT(quals_good[6]==1);
  CU_ASSERT(quals_good[7]==1);
  CU_ASSERT(quals_good[8]==1);
  CU_ASSERT(quals_good[9]==1);
  CU_ASSERT(quals_good[10]==1);
  CU_ASSERT(quals_good[11]==0);

  //not set - did not find a good kmer
  CU_ASSERT(first_good_kmer==-1);


  //this time the policy makes a difference

  first_good_kmer=-1;  
  set_int_array(quals_good, read_len, 1);
  qual_cutoff=10+33;
  strand_first_good_kmer=reverse;//just testing
  result = get_first_good_kmer_and_populate_qual_array(readid,read_seq, read_qual,
						       read_len-5+1, read_len,
						       quals_good,
						       qual_cutoff, &first_good_kmer, &strand_first_good_kmer,
						       db_graph, DiscardReadIfLowQualBaseUnCorrectable, false);



  CU_ASSERT(result==Discard);

  CU_ASSERT(quals_good[0]==1);
  CU_ASSERT(quals_good[1]==1);
  CU_ASSERT(quals_good[2]==1);
  CU_ASSERT(quals_good[3]==1);
  CU_ASSERT(quals_good[4]==1);
  CU_ASSERT(quals_good[5]==1);
  CU_ASSERT(quals_good[6]==1);
  CU_ASSERT(quals_good[7]==1);
  CU_ASSERT(quals_good[8]==1);
  CU_ASSERT(quals_good[9]==1);
  CU_ASSERT(quals_good[10]==1);
  CU_ASSERT(quals_good[11]==0);

  //not set also
  CU_ASSERT(first_good_kmer==-1);





  hash_table_free(&db_graph);
  strbuf_free(read_seq);
  strbuf_free(read_qual);





}


void test_error_correct_file_against_graph()
{
  //simple test, with one good kmer in middle. All low qual. Check error corrects left and right.
  // include a check of the stats collection


  //first set up the hash/graph
  int kmer_size = 5;
  int number_of_bits = 10;
  int bucket_size = 10;

  //*********************************************

  //*********************************************

  dBGraph * db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);
  
  //Load the following fasta: this is to be the trusted graph
  // >
  // ACGGCTTTACGGT



  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/error_correction/graph3.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);


  //now test this fastq
  // @only one kmer in graph, CTTTA
  // CCCCCTTTAAAAT
  // +
  // #############


  //we will start at CTTTA and move right, error correcting, and then go back to CTTTA and move left
  //
  //graph      ACGGCTTTACGGT
  // read      CCCCCTTTAAAAT
  // step1     ccccctttaCaat  << correct base at posn 9
  // step2     ccccctttaCGat  << correct base at posn 10
  // step3     ccccctttaCGGt << correct base at posn 11. 
  //                           < no need to correct posn 12
  // step4     cccGctttaCGGt  << correct base at posn 3
  // step5     ccGGctttaCGGt  << correct base at posn 2
  //                          no need to correct posn 1
  // step6     AcGGctttaCGGt  << correct base at posn 0


  char quality_cutoff= 10;
  char ascii_qual_offset = 33;
  char* outfile = "../data/test/error_correction/fq5_for_comparing_with_graph3.err_corrected.fq";
  uint64_t bases_modified_count_array[20];
  uint64_t posn_modified_count_array[20];
  int bases_modified_count_array_size=20;
  int i;
  for (i=0; i<20; i++)
    {
      bases_modified_count_array[i]=0;
      posn_modified_count_array[i]=0;
    }
  boolean add_greedy_bases_for_better_bwt_compression=false;
  int num_greedy_bases=0;
  boolean rev_comp_read_if_on_reverse_strand=false;

  error_correct_file_against_graph("../data/test/error_correction/fq5_for_comparing_with_graph3.fq", 
				   quality_cutoff, ascii_qual_offset,
				   db_graph, outfile,
				   bases_modified_count_array,
				   posn_modified_count_array,
				   bases_modified_count_array_size,
				   DontWorryAboutLowQualBaseUnCorrectable,
				   add_greedy_bases_for_better_bwt_compression, num_greedy_bases, rev_comp_read_if_on_reverse_strand);

  CU_ASSERT(bases_modified_count_array[0]==0);
  CU_ASSERT(bases_modified_count_array[1]==0);
  CU_ASSERT(bases_modified_count_array[2]==0);
  CU_ASSERT(bases_modified_count_array[3]==0);
  CU_ASSERT(bases_modified_count_array[4]==0);
  CU_ASSERT(bases_modified_count_array[5]==0);
  CU_ASSERT(bases_modified_count_array[6]==1);
  CU_ASSERT(bases_modified_count_array[7]==0);
  CU_ASSERT(bases_modified_count_array[8]==0);
  CU_ASSERT(bases_modified_count_array[9]==0);
  CU_ASSERT(bases_modified_count_array[10]==0);
  CU_ASSERT(bases_modified_count_array[11]==0);
  CU_ASSERT(bases_modified_count_array[12]==0);
  CU_ASSERT(bases_modified_count_array[13]==0);
  CU_ASSERT(bases_modified_count_array[14]==0);

  CU_ASSERT(posn_modified_count_array[0]==1);
  CU_ASSERT(posn_modified_count_array[1]==0);
  CU_ASSERT(posn_modified_count_array[2]==1);
  CU_ASSERT(posn_modified_count_array[3]==1);
  CU_ASSERT(posn_modified_count_array[4]==0);
  CU_ASSERT(posn_modified_count_array[5]==0);
  CU_ASSERT(posn_modified_count_array[6]==0);
  CU_ASSERT(posn_modified_count_array[7]==0);
  CU_ASSERT(posn_modified_count_array[8]==0);
  CU_ASSERT(posn_modified_count_array[9]==1);
  CU_ASSERT(posn_modified_count_array[10]==1);
  CU_ASSERT(posn_modified_count_array[11]==1);
  CU_ASSERT(posn_modified_count_array[12]==0);


  SeqFile* sf = seq_file_open(outfile);
  if(sf == NULL)
    {
      // Error opening file
      fprintf(stderr, "Error: cannot read seq file '%s'\n", outfile);
      exit(EXIT_FAILURE);
    }

  seq_next_read(sf);
  StrBuf* read_seq  = strbuf_new();
  StrBuf* read_qual  = strbuf_new();

  seq_read_all_bases_and_quals(sf, read_seq, read_qual);

  seq_file_close(sf);


  //graph is this          ACGGCTTTACGGT
  // input data is this    CCCCCTTTAAAAT  quals  #############
  // correction should be  ACGGCTTTACGGT  quals  ,#,,#####,,,# (at fixed bases, qual is fixed to be 1 more than cutoff)


  CU_ASSERT(strcmp(read_seq->buff, "ACGGCTTTACGGT")==0);
  CU_ASSERT(strcmp(read_qual->buff, ",#,,#####,,,#")==0);


  //Try again with an almost identical fastq, but which
  //has high quality for one of the previously corrected
  //same, but make sure does not fix high qual.
  for (i=0; i<20; i++)
    {
      bases_modified_count_array[i]=0;
      posn_modified_count_array[i]=0;
    }


  //single read is:
  //  @only one kmer in graph, CTTTA. Penultimate base has ghigh quality so should not be corrected
  //  CCCCCTTTAAAAT
  //  +
  //  ###########]#

  char* outfile2 = "../data/test/error_correction/fq6_for_comparing_with_graph3.err_corrected.fq";
  error_correct_file_against_graph("../data/test/error_correction/fq6_for_comparing_with_graph3.fq", 
				   quality_cutoff, ascii_qual_offset,
				   db_graph, outfile2,
				   bases_modified_count_array,
				   posn_modified_count_array,
				   bases_modified_count_array_size,
				   DontWorryAboutLowQualBaseUnCorrectable,
				   add_greedy_bases_for_better_bwt_compression, num_greedy_bases, rev_comp_read_if_on_reverse_strand);

  CU_ASSERT(bases_modified_count_array[0]==0);
  CU_ASSERT(bases_modified_count_array[1]==0);
  CU_ASSERT(bases_modified_count_array[2]==0);
  CU_ASSERT(bases_modified_count_array[3]==0);
  CU_ASSERT(bases_modified_count_array[4]==0);
  CU_ASSERT(bases_modified_count_array[5]==1);
  CU_ASSERT(bases_modified_count_array[6]==0);
  CU_ASSERT(bases_modified_count_array[7]==0);
  CU_ASSERT(bases_modified_count_array[8]==0);
  CU_ASSERT(bases_modified_count_array[9]==0);
  CU_ASSERT(bases_modified_count_array[10]==0);
  CU_ASSERT(bases_modified_count_array[11]==0);
  CU_ASSERT(bases_modified_count_array[12]==0);
  CU_ASSERT(bases_modified_count_array[13]==0);
  CU_ASSERT(bases_modified_count_array[14]==0);

  CU_ASSERT(posn_modified_count_array[0]==1);
  CU_ASSERT(posn_modified_count_array[1]==0);
  CU_ASSERT(posn_modified_count_array[2]==1);
  CU_ASSERT(posn_modified_count_array[3]==1);
  CU_ASSERT(posn_modified_count_array[4]==0);
  CU_ASSERT(posn_modified_count_array[5]==0);
  CU_ASSERT(posn_modified_count_array[6]==0);
  CU_ASSERT(posn_modified_count_array[7]==0);
  CU_ASSERT(posn_modified_count_array[8]==0);
  CU_ASSERT(posn_modified_count_array[9]==1);
  CU_ASSERT(posn_modified_count_array[10]==1);
  CU_ASSERT(posn_modified_count_array[11]==0);
  CU_ASSERT(posn_modified_count_array[12]==0);


  sf = seq_file_open(outfile2);
  if(sf == NULL)
    {
      // Error opening file
      fprintf(stderr, "Error: cannot read seq file '%s'\n", outfile2);
      exit(EXIT_FAILURE);
    }

  seq_next_read(sf);
  strbuf_reset(read_seq);
  strbuf_reset(read_qual);
  seq_read_all_bases_and_quals(sf, read_seq, read_qual);
  seq_file_close(sf);


  //graph is this          ACGGCTTTACGGT
  // input data is this    CCCCCTTTAAAAT    quas   ###########]#
  // correction should be  ACGGCTTTACGAT    quals  ,#,,#####,,]#


  CU_ASSERT(strcmp(read_seq->buff, "ACGGCTTTACGAT")==0);
  CU_ASSERT(strcmp(read_qual->buff, ",#,,#####,,]#")==0);


  //now, just to double check the statistics, try one fastq that contains both previous reads

  for (i=0; i<20; i++)
    {
      bases_modified_count_array[i]=0;
      posn_modified_count_array[i]=0;
    }
  char* outfile3 = "../data/test/error_correction/fq7_for_comparing_with_graph3.err_corrected.fq";
  error_correct_file_against_graph("../data/test/error_correction/fq7_for_comparing_with_graph3.fq", 
				   quality_cutoff, ascii_qual_offset,
				   db_graph, outfile3,
				   bases_modified_count_array,
				   posn_modified_count_array,
				   bases_modified_count_array_size,
				   DontWorryAboutLowQualBaseUnCorrectable,
				   add_greedy_bases_for_better_bwt_compression, num_greedy_bases, rev_comp_read_if_on_reverse_strand);

  CU_ASSERT(bases_modified_count_array[0]==0);
  CU_ASSERT(bases_modified_count_array[1]==0);
  CU_ASSERT(bases_modified_count_array[2]==0);
  CU_ASSERT(bases_modified_count_array[3]==0);
  CU_ASSERT(bases_modified_count_array[4]==0);
  CU_ASSERT(bases_modified_count_array[5]==1);
  CU_ASSERT(bases_modified_count_array[6]==1);
  CU_ASSERT(bases_modified_count_array[7]==0);
  CU_ASSERT(bases_modified_count_array[8]==0);
  CU_ASSERT(bases_modified_count_array[9]==0);
  CU_ASSERT(bases_modified_count_array[10]==0);
  CU_ASSERT(bases_modified_count_array[11]==0);
  CU_ASSERT(bases_modified_count_array[12]==0);
  CU_ASSERT(bases_modified_count_array[13]==0);
  CU_ASSERT(bases_modified_count_array[14]==0);

  CU_ASSERT(posn_modified_count_array[0]==2);
  CU_ASSERT(posn_modified_count_array[1]==0);
  CU_ASSERT(posn_modified_count_array[2]==2);
  CU_ASSERT(posn_modified_count_array[3]==2);
  CU_ASSERT(posn_modified_count_array[4]==0);
  CU_ASSERT(posn_modified_count_array[5]==0);
  CU_ASSERT(posn_modified_count_array[6]==0);
  CU_ASSERT(posn_modified_count_array[7]==0);
  CU_ASSERT(posn_modified_count_array[8]==0);
  CU_ASSERT(posn_modified_count_array[9]==2);
  CU_ASSERT(posn_modified_count_array[10]==2);
  CU_ASSERT(posn_modified_count_array[11]==1);
  CU_ASSERT(posn_modified_count_array[12]==0);


  sf = seq_file_open(outfile3);
  if(sf == NULL)
    {
      // Error opening file
      fprintf(stderr, "Error: cannot read seq file '%s'\n", outfile3);
      exit(EXIT_FAILURE);
    }

  seq_next_read(sf);
  strbuf_reset(read_seq);
  strbuf_reset(read_qual);
  seq_read_all_bases_and_quals(sf, read_seq, read_qual);
  CU_ASSERT(strcmp(read_seq->buff, "ACGGCTTTACGGT")==0);
  CU_ASSERT(strcmp(read_qual->buff, ",#,,#####,,,#")==0);

  seq_next_read(sf);
  seq_read_all_bases_and_quals(sf, read_seq, read_qual);
  CU_ASSERT(strcmp(read_seq->buff, "ACGGCTTTACGAT")==0);
  CU_ASSERT(strcmp(read_qual->buff,",#,,#####,,]#")==0);
  seq_file_close(sf);



  //last test - does the policy correctly get executed, for when to discard reads?

  // @only one kmer in graph, CTTTA.To right of CTTTA all bases high qual.Only far left is low qual, but this cannot be corrected
  // CCCCCTTTAAAAT
  // +
  // #]]]#####]]]]
  for (i=0; i<20; i++)
    {
      bases_modified_count_array[i]=0;
      posn_modified_count_array[i]=0;
    }
  char* outfile4 = "../data/test/error_correction/fq8_for_comparing_with_graph3.err_corrected.fq";
  error_correct_file_against_graph("../data/test/error_correction/fq8_for_comparing_with_graph3.fq", 
				   quality_cutoff, ascii_qual_offset,
				   db_graph, outfile4,
				   bases_modified_count_array,
				   posn_modified_count_array,
				   bases_modified_count_array_size,
				   DontWorryAboutLowQualBaseUnCorrectable,
				   add_greedy_bases_for_better_bwt_compression, num_greedy_bases, rev_comp_read_if_on_reverse_strand);

  CU_ASSERT(bases_modified_count_array[0]==1);
  CU_ASSERT(bases_modified_count_array[1]==0);
  CU_ASSERT(bases_modified_count_array[2]==0);
  CU_ASSERT(bases_modified_count_array[3]==0);
  CU_ASSERT(bases_modified_count_array[4]==0);
  CU_ASSERT(bases_modified_count_array[5]==0);
  CU_ASSERT(bases_modified_count_array[6]==0);
  CU_ASSERT(bases_modified_count_array[7]==0);
  CU_ASSERT(bases_modified_count_array[8]==0);
  CU_ASSERT(bases_modified_count_array[9]==0);
  CU_ASSERT(bases_modified_count_array[10]==0);
  CU_ASSERT(bases_modified_count_array[11]==0);
  CU_ASSERT(bases_modified_count_array[12]==0);
  CU_ASSERT(bases_modified_count_array[13]==0);
  CU_ASSERT(bases_modified_count_array[14]==0);

  CU_ASSERT(posn_modified_count_array[0]==0);
  CU_ASSERT(posn_modified_count_array[1]==0);
  CU_ASSERT(posn_modified_count_array[2]==0);
  CU_ASSERT(posn_modified_count_array[3]==0);
  CU_ASSERT(posn_modified_count_array[4]==0);
  CU_ASSERT(posn_modified_count_array[5]==0);
  CU_ASSERT(posn_modified_count_array[6]==0);
  CU_ASSERT(posn_modified_count_array[7]==0);
  CU_ASSERT(posn_modified_count_array[8]==0);
  CU_ASSERT(posn_modified_count_array[9]==0);
  CU_ASSERT(posn_modified_count_array[10]==0);
  CU_ASSERT(posn_modified_count_array[11]==0);
  CU_ASSERT(posn_modified_count_array[12]==0);


  sf = seq_file_open(outfile4);
  if(sf == NULL)
    {
      // Error opening file
      fprintf(stderr, "Error: cannot read seq file '%s'\n", outfile4);
      exit(EXIT_FAILURE);
    }

  seq_next_read(sf);
  strbuf_reset(read_seq);
  strbuf_reset(read_qual);
  seq_read_all_bases_and_quals(sf, read_seq, read_qual);
  CU_ASSERT(strcmp(read_seq->buff, "CCCCCTTTAAAAT")==0);
  CU_ASSERT(strcmp(read_qual->buff, "#]]]#####]]]]")==0);
  seq_file_close(sf);



  //and again, this time changing the policy so you want it discarded
  for (i=0; i<20; i++)
    {
      bases_modified_count_array[i]=0;
      posn_modified_count_array[i]=0;
    }
  char* outfile5 = "../data/test/error_correction/fq8_for_comparing_with_graph3.err_corrected.fa2";
  error_correct_file_against_graph("../data/test/error_correction/fq8_for_comparing_with_graph3.fq", 
				   quality_cutoff, ascii_qual_offset,
				   db_graph, outfile5,
				   bases_modified_count_array,
				   posn_modified_count_array,
				   bases_modified_count_array_size,
				   DiscardReadIfLowQualBaseUnCorrectable,
				   add_greedy_bases_for_better_bwt_compression, num_greedy_bases, rev_comp_read_if_on_reverse_strand);

  CU_ASSERT(bases_modified_count_array[0]==0);
  CU_ASSERT(bases_modified_count_array[1]==0);
  CU_ASSERT(bases_modified_count_array[2]==0);
  CU_ASSERT(bases_modified_count_array[3]==0);
  CU_ASSERT(bases_modified_count_array[4]==0);
  CU_ASSERT(bases_modified_count_array[5]==0);
  CU_ASSERT(bases_modified_count_array[6]==0);
  CU_ASSERT(bases_modified_count_array[7]==0);
  CU_ASSERT(bases_modified_count_array[8]==0);
  CU_ASSERT(bases_modified_count_array[9]==0);
  CU_ASSERT(bases_modified_count_array[10]==0);
  CU_ASSERT(bases_modified_count_array[11]==0);
  CU_ASSERT(bases_modified_count_array[12]==0);
  CU_ASSERT(bases_modified_count_array[13]==0);
  CU_ASSERT(bases_modified_count_array[14]==0);

  CU_ASSERT(posn_modified_count_array[0]==0);
  CU_ASSERT(posn_modified_count_array[1]==0);
  CU_ASSERT(posn_modified_count_array[2]==0);
  CU_ASSERT(posn_modified_count_array[3]==0);
  CU_ASSERT(posn_modified_count_array[4]==0);
  CU_ASSERT(posn_modified_count_array[5]==0);
  CU_ASSERT(posn_modified_count_array[6]==0);
  CU_ASSERT(posn_modified_count_array[7]==0);
  CU_ASSERT(posn_modified_count_array[8]==0);
  CU_ASSERT(posn_modified_count_array[9]==0);
  CU_ASSERT(posn_modified_count_array[10]==0);
  CU_ASSERT(posn_modified_count_array[11]==0);
  CU_ASSERT(posn_modified_count_array[12]==0);


  sf = seq_file_open(outfile5);
  if(sf == NULL)
    {
      // Error opening file
      fprintf(stderr, "Error: cannot read seq file '%s'\n", outfile4);
      exit(EXIT_FAILURE);
    }

  seq_next_read(sf);
  CU_ASSERT(seq_read_all_bases(sf, read_seq)==0);

  seq_file_close(sf);


  strbuf_free(read_seq);
  strbuf_free(read_qual);
  hash_table_free(&db_graph);
}




void test_take_n_greedy_random_steps()
{
    //first set up the hash/graph
  int kmer_size = 7;
  int number_of_bits = 4;
  int bucket_size = 10;

  double  avg_coverage;
  Covg min_coverage, max_coverage;

  dBGraph * db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);

 

  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  //load this fasta
  // >read 1
  // AGGAACGTCCGCCATTAGACC
  // >read 2 - differs by a SNP
  // AGGAACGTACGCCATTAGACC
  //        ^   
  //graph s just a SNP bubble, at k=7

  load_se_filelist_into_graph_colour(
				     "../data/test/error_correction/graph4.falist",
				     fq_quality_cutoff, homopolymer_cutoff,
				     remove_duplicates_se, ascii_fq_offset,
				     into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
				     &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
				     NULL, 0, &subsample_null);

  BinaryKmer tmp_kmer1;
  BinaryKmer tmp_kmer2;

  dBNode* first_node = hash_table_find(element_get_key(seq_to_binary_kmer("AGGAACG", kmer_size,&tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(first_node!=NULL);

  StrBuf* greedy_seq=strbuf_new();
  take_n_greedy_random_steps(first_node, forward, db_graph, 1, greedy_seq);  
  CU_ASSERT(strcmp(greedy_seq->buff,"T")==0);
  strbuf_reset(greedy_seq);

  take_n_greedy_random_steps(first_node, forward, db_graph, 2, greedy_seq);  
  CU_ASSERT(  (strcmp(greedy_seq->buff,"TC")==0) || (strcmp(greedy_seq->buff,"TA")==0) );
  strbuf_reset(greedy_seq);

  take_n_greedy_random_steps(first_node, forward, db_graph, 9, greedy_seq);  
  CU_ASSERT(  (strcmp(greedy_seq->buff,"TCCGCCATT")==0) || (strcmp(greedy_seq->buff,"TACGCCATT")==0) );

  strbuf_reset(greedy_seq);

  take_n_greedy_random_steps(first_node, forward, db_graph, 16, greedy_seq);  
  CU_ASSERT(  (strcmp(greedy_seq->buff,"TCCGCCATTAGACCAA")==0) || (strcmp(greedy_seq->buff,"TACGCCATTAGACCAA")==0) );//pad with A's
  strbuf_reset(greedy_seq);


  take_n_greedy_random_steps(first_node, reverse, db_graph, 5, greedy_seq);  
  CU_ASSERT(strcmp(greedy_seq->buff,"AAAAA")==0);
  strbuf_reset(greedy_seq);

  strbuf_free(greedy_seq);
  hash_table_free(&db_graph);
}

void test_reverse_comp_according_ref_pos_strand()
{
      //first set up the hash/graph
  int kmer_size = 7;
  int number_of_bits = 4;
  int bucket_size = 10;

  double  avg_coverage;
  Covg min_coverage, max_coverage;

  dBGraph * db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);

 

  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  //load this fasta
  // >read 1
  // AGGAACGTCCGCCATTAGACC
  // >read 2 - differs by a SNP
  // AGGAACGTACGCCATTAGACC
  //        ^   
  //graph s just a SNP bubble, at k=7

  load_se_filelist_into_graph_colour(
				     "../data/test/error_correction/graph4.falist",
				     fq_quality_cutoff, homopolymer_cutoff,
				     remove_duplicates_se, ascii_fq_offset,
				     into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
				     &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
				     NULL, 0, &subsample_null);

  BinaryKmer tmp_kmer1;
  BinaryKmer tmp_kmer2;

  read_ref_fasta_and_mark_strand("../data/test/error_correction/graph4.fa", db_graph);

  //arrays for stats
  uint64_t bases_modified_count_array[50];
  uint64_t posn_modified_count_array[50];
  int bases_modified_count_array_size=50;
  int i;
  for (i=0; i<50; i++)
    {
      bases_modified_count_array[i]=0;
      posn_modified_count_array[i]=0;
    }
  boolean add_greedy_bases_for_better_bwt_compression=false;
  int num_greedy_bases=0;
  boolean rev_comp_read_if_on_reverse_strand=true;

  //Now try correcting a fastq that has one read = rev complement of a section of the reference
  char* outfile = "../data/test/error_correction/fq10_revcomp_of_fq9_for_comparing_with_graph4.error_corrected.fq";

  error_correct_file_against_graph("../data/test/error_correction/fq10_revcomp_of_fq9_for_comparing_with_graph4.fq",
				   0, 33, db_graph, outfile,
				   bases_modified_count_array, posn_modified_count_array,
				   50, 
				   DontWorryAboutLowQualBaseUnCorrectable,
				   add_greedy_bases_for_better_bwt_compression, 
				   num_greedy_bases, 
				   rev_comp_read_if_on_reverse_strand);

  CU_ASSERT(bases_modified_count_array[0]==1);//1 read is uncorrected
  for (i=1; i<21; i++)
    {
      CU_ASSERT(bases_modified_count_array[i]==0);
      CU_ASSERT(posn_modified_count_array[i]==0);
    }
  SeqFile* sf = seq_file_open(outfile);
  if(sf == NULL)
    {
      // Error opening file
      fprintf(stderr, "Error: cannot read seq file '%s'\n", outfile);
      exit(EXIT_FAILURE);
    }

  seq_next_read(sf);
  StrBuf* read_seq  = strbuf_new();
  StrBuf* read_qual  = strbuf_new();
  seq_read_all_bases_and_quals(sf, read_seq, read_qual);
  //expect to get reverse complement of read
  CU_ASSERT(strcmp(read_seq->buff, "AGGAACGTCCGCCATTAGACC")==0);
  //and reverse of quals
  CU_ASSERT(strcmp(read_qual->buff, "LKJIIIIIIIIIIIIIIIIII")==0);
  seq_file_close(sf);

  hash_table_free(&db_graph);



  //second test


  //first set up the hash/graph
  kmer_size = 5;
  number_of_bits = 10;
  bucket_size = 10;

  //*********************************************

  //*********************************************

  db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);
  
  //Load the following fasta: this is to be the trusted graph
  // >
  // ACGGCTTTACGGT



  fq_quality_cutoff = 0;
  homopolymer_cutoff = 0;
  remove_duplicates_se = false;
  ascii_fq_offset = 33;
  into_colour = 0;

  files_loaded = 0;
  bad_reads = 0;
  dup_reads = 0;
  seq_loaded = 0;
  seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/error_correction/graph3.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);


  //now mark reference in reverse direction to the file i will correct
  read_ref_fasta_and_mark_strand("../data/test/error_correction/graph3_revcomp.fa", 
				 db_graph);

  //now test this fastq
  // @only one kmer in graph, CTTTA
  // CCCCCTTTAAAAT
  // +
  // #############


  //we will start at CTTTA and move right, error correcting, and then go back to CTTTA and move left
  //
  //graph      ACGGCTTTACGGT
  // read      CCCCCTTTAAAAT
  // step1     ccccctttaCaat  << correct base at posn 9
  // step2     ccccctttaCGat  << correct base at posn 10
  // step3     ccccctttaCGGt << correct base at posn 11. 
  //                           < no need to correct posn 12
  // step4     cccGctttaCGGt  << correct base at posn 3
  // step5     ccGGctttaCGGt  << correct base at posn 2
  //                          no need to correct posn 1
  // step6     AcGGctttaCGGt  << correct base at posn 0


  char quality_cutoff= 10;
  char ascii_qual_offset = 33;
  char* outfile2 = "../data/test/error_correction/fq5_for_comparing_with_graph3.err_corrected.fq";

  for (i=0; i<20; i++)
    {
      bases_modified_count_array[i]=0;
      posn_modified_count_array[i]=0;
    }
  add_greedy_bases_for_better_bwt_compression=false;
  num_greedy_bases=0;
  rev_comp_read_if_on_reverse_strand=true;
  error_correct_file_against_graph("../data/test/error_correction/fq5_for_comparing_with_graph3.fq", 
				   quality_cutoff, ascii_qual_offset,
				   db_graph, outfile2,
				   bases_modified_count_array,
				   posn_modified_count_array,
				   bases_modified_count_array_size,
				   DontWorryAboutLowQualBaseUnCorrectable,
				   add_greedy_bases_for_better_bwt_compression, 
				   num_greedy_bases, 
				   rev_comp_read_if_on_reverse_strand);

  CU_ASSERT(bases_modified_count_array[0]==0);
  CU_ASSERT(bases_modified_count_array[1]==0);
  CU_ASSERT(bases_modified_count_array[2]==0);
  CU_ASSERT(bases_modified_count_array[3]==0);
  CU_ASSERT(bases_modified_count_array[4]==0);
  CU_ASSERT(bases_modified_count_array[5]==0);
  CU_ASSERT(bases_modified_count_array[6]==1);
  CU_ASSERT(bases_modified_count_array[7]==0);
  CU_ASSERT(bases_modified_count_array[8]==0);
  CU_ASSERT(bases_modified_count_array[9]==0);
  CU_ASSERT(bases_modified_count_array[10]==0);
  CU_ASSERT(bases_modified_count_array[11]==0);
  CU_ASSERT(bases_modified_count_array[12]==0);
  CU_ASSERT(bases_modified_count_array[13]==0);
  CU_ASSERT(bases_modified_count_array[14]==0);

  CU_ASSERT(posn_modified_count_array[0]==1);
  CU_ASSERT(posn_modified_count_array[1]==0);
  CU_ASSERT(posn_modified_count_array[2]==1);
  CU_ASSERT(posn_modified_count_array[3]==1);
  CU_ASSERT(posn_modified_count_array[4]==0);
  CU_ASSERT(posn_modified_count_array[5]==0);
  CU_ASSERT(posn_modified_count_array[6]==0);
  CU_ASSERT(posn_modified_count_array[7]==0);
  CU_ASSERT(posn_modified_count_array[8]==0);
  CU_ASSERT(posn_modified_count_array[9]==1);
  CU_ASSERT(posn_modified_count_array[10]==1);
  CU_ASSERT(posn_modified_count_array[11]==1);
  CU_ASSERT(posn_modified_count_array[12]==0);


  sf = seq_file_open(outfile2);
  if(sf == NULL)
    {
      // Error opening file
      fprintf(stderr, "Error: cannot read seq file '%s'\n", outfile2);
      exit(EXIT_FAILURE);
    }

  seq_next_read(sf);
  strbuf_reset(read_seq);
  strbuf_reset(read_qual);
  seq_read_all_bases_and_quals(sf, read_seq, read_qual);
  seq_file_close(sf);


  //graph is this          ACGGCTTTACGGT
  // input data is this    CCCCCTTTAAAAT quals ###########
  // correction should be  ACGGCTTTACGGT
  //  i get rev comp of    ACGGCTTTACGGT = ACCGTAAAGCCGT and quals=reverse of   ,#,,#####,,,# - ie. #,,,#####,,#,

  CU_ASSERT(strcmp(read_seq->buff, "ACCGTAAAGCCGT")==0);
  CU_ASSERT(strcmp(read_qual->buff, "#,,,#####,,#,")==0);

  hash_table_free(&db_graph);
  strbuf_free(read_seq);
  strbuf_free(read_qual);


}
