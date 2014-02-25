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
  test_seq.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <CUnit.h>
#include <Basic.h>

// cortex_var headers
#include "seq.h"
#include "test_seq.h"


//test reading several short entries from a fasta file
void test_read_sequence_from_fasta(){

  Sequence * seq = malloc(sizeof(Sequence));
  boolean full_entry;

  if (seq == NULL){
    die("Out of memory trying to allocate Sequence\n");
  }

  //pre-allocate space where to read the sequences
  alloc_sequence(seq,200,LINE_MAX);

  int length_seq;
  FILE* fp1 = fopen("../data/test/basic/one_entry.fa", "r");
  
  if (fp1 == NULL){
    die("cannot open file:../data/test/basic/one_entry.fa\n");
  }

  // 1. Read from simple fasta:
  // >Zam
  // ACGT
  // ACGTACGTACGT
  
  length_seq = read_sequence_from_fasta(fp1,seq,1000,true,&full_entry,0);

  CU_ASSERT_EQUAL(length_seq, 16);
  CU_ASSERT_STRING_EQUAL("Zam",seq->name);
  CU_ASSERT_EQUAL(1,seq->start);
  CU_ASSERT_EQUAL(16,seq->end);
  CU_ASSERT_STRING_EQUAL("ACGTACGTACGTACGT",seq->seq);
  CU_ASSERT(full_entry);
  fclose(fp1);

  FILE* fp2 = fopen("../data/test/basic/three_entries.fa", "r");

  // 2. Read from fasta:
  //>Zam1
  //ACGT
  //ACGTACGTACGT
  //>Zam2
  //ACGT
  //ACGTACGTACGT
  //TTTTTTTT
  //>Zam3
  //ACGTNNACGTACGTACGT

  length_seq = read_sequence_from_fasta(fp2,seq,1000,true,&full_entry,0);
  
  CU_ASSERT_EQUAL(length_seq, 16);
  CU_ASSERT_STRING_EQUAL("Zam1",seq->name);
  CU_ASSERT_EQUAL(1,seq->start);
  CU_ASSERT_EQUAL(16,seq->end);
  CU_ASSERT_STRING_EQUAL("ACGTACGTACGTACGT",seq->seq);
  CU_ASSERT(full_entry);

  length_seq = read_sequence_from_fasta(fp2,seq,1000,true,&full_entry,0);
  
  CU_ASSERT_EQUAL(length_seq, 24);
  CU_ASSERT_STRING_EQUAL("Zam2",seq->name);
  CU_ASSERT_EQUAL(1,seq->start);
  CU_ASSERT_EQUAL(24,seq->end);
  CU_ASSERT_STRING_EQUAL("ACGTACGTACGTACGTTTTTTTTt",seq->seq);
  CU_ASSERT(full_entry);
   

  length_seq = read_sequence_from_fasta(fp2,seq,1000,true,&full_entry,0);

  
  CU_ASSERT_EQUAL(length_seq, 18);
  CU_ASSERT_STRING_EQUAL("Zam3",seq->name);
  CU_ASSERT_EQUAL(1,seq->start);
  CU_ASSERT_EQUAL(18,seq->end);
  CU_ASSERT_STRING_EQUAL("ACGTNNACGTACGTACGT",seq->seq);
  CU_ASSERT(full_entry);

  length_seq = read_sequence_from_fasta(fp2,seq,1000,true,&full_entry,0);

  CU_ASSERT_EQUAL(length_seq, 0);
  CU_ASSERT(full_entry);
   
  fclose(fp2);
  free_sequence(&seq);

}



//test reading from a fasta file with long entries 

void test_read_sequence_from_long_fasta(){

  Sequence * seq = malloc(sizeof(Sequence));
  boolean full_entry;

  if (seq == NULL){							
    die("Out of memory trying to allocate Sequence\n");								
  }
  //pre-allocate space where to read the sequences
  alloc_sequence(seq,200,LINE_MAX);

  int length_seq;
  FILE* fp1 = fopen("../data/test/basic/long_entries.fa", "r");
  
  if (fp1 == NULL){							
    die("cannot open file: ../data/test/basic/long_entries.fa\n");								
  }
  
  length_seq = read_sequence_from_fasta(fp1,seq,10,true,&full_entry,0);
  CU_ASSERT_EQUAL(length_seq, 10);
  CU_ASSERT_EQUAL(seq->start,1);
  CU_ASSERT_EQUAL(seq->end,10);
  CU_ASSERT_STRING_EQUAL("Mario",seq->name);
  CU_ASSERT_STRING_EQUAL("ACGTACGTAC",seq->seq);
  CU_ASSERT(full_entry==false);

  length_seq = read_sequence_from_fasta(fp1,seq,10,false,&full_entry,0);
  CU_ASSERT_EQUAL(length_seq, 10);
  CU_ASSERT_EQUAL(seq->start,11);
  CU_ASSERT_EQUAL(seq->end,20);
  CU_ASSERT_STRING_EQUAL("Mario",seq->name);
  CU_ASSERT_STRING_EQUAL("GTACGTAAAA",seq->seq);
 
  seq->seq[0]='T';
  seq->seq[1]='T';
  seq->seq[2]='T';
  length_seq = read_sequence_from_fasta(fp1,seq,10,false,&full_entry,3);
  CU_ASSERT_EQUAL(length_seq, 10);
  CU_ASSERT_STRING_EQUAL("Mario",seq->name);
  CU_ASSERT_EQUAL(seq->start, 21);
  CU_ASSERT_EQUAL(seq->end,27);
  CU_ASSERT_STRING_EQUAL("TTTAAAAAAA",seq->seq);

  //finish off the entry
  length_seq = read_sequence_from_fasta(fp1,seq,1000,false,&full_entry,0); 
  CU_ASSERT(full_entry == true);
  
  length_seq = read_sequence_from_fasta(fp1,seq,16,true,&full_entry,0); 
  CU_ASSERT(full_entry == true); 
  CU_ASSERT_EQUAL(seq->start, 1); 
  CU_ASSERT_EQUAL(seq->end,16); 
  CU_ASSERT_EQUAL(length_seq, 16); 
  CU_ASSERT_STRING_EQUAL("Pepe",seq->name);
  
  length_seq = read_sequence_from_fasta(fp1,seq,3,true,&full_entry,0); 
  CU_ASSERT(full_entry == false);
  CU_ASSERT_EQUAL(seq->start, 1);
  CU_ASSERT_EQUAL(seq->end,3);
  CU_ASSERT_EQUAL(length_seq, 3);
  CU_ASSERT_STRING_EQUAL("COCO",seq->name);
  CU_ASSERT_STRING_EQUAL("TTT",seq->seq);

  length_seq = read_sequence_from_fasta(fp1,seq,1000,false,&full_entry,0);
  CU_ASSERT(full_entry == true);
  CU_ASSERT_EQUAL(seq->start, 4);
  CU_ASSERT_EQUAL(seq->end,10);
  CU_ASSERT_EQUAL(length_seq,7);
  CU_ASSERT_STRING_EQUAL("COCO",seq->name);
  CU_ASSERT_STRING_EQUAL("TAAAATT",seq->seq);

  length_seq = read_sequence_from_fasta(fp1,seq,15,true,&full_entry,0);
  CU_ASSERT(full_entry == true);
  CU_ASSERT_EQUAL(seq->start, 1);
  CU_ASSERT_EQUAL(seq->end,15);
  CU_ASSERT_EQUAL(length_seq, 15);
  CU_ASSERT_STRING_EQUAL("CACHO",seq->name);
  CU_ASSERT_STRING_EQUAL("TTTTTTAAAGGATAT",seq->seq);


  length_seq = read_sequence_from_fasta(fp1,seq,15,true,&full_entry,0);
  CU_ASSERT(full_entry == true);
  CU_ASSERT_EQUAL(length_seq, 0);

  fclose(fp1);
  free_sequence(&seq);

}


//bad chars are reported in a message and skipt when building kmers out a sequence

void test_read_sequence_from_fasta_when_file_has_bad_reads()
{

  int length_seq;
  Sequence * seq = malloc(sizeof(Sequence));
  boolean full_entry;
  
  if (seq == NULL){							
    die("Out of memory trying to allocate Sequence\n");								
  }
  //pre-allocate space where to read the sequences
  int max_read_length=100;
  alloc_sequence(seq,max_read_length,LINE_MAX);
 
  FILE* fp2= fopen("../data/test/basic/includes_reads_that_have_bad_characters.fa", "r");

  // >read1
  // AAAAAAAAAAAA9
  // >read2
  // ¡€#9∞§¶#¶•#•#•#ª#ª#ª#ªº#º#º#º––––
  // >read3 4 c's
  // CCCC
  // >read4 10 Ts
  // TTTTTTTTTT
  // >read5
  // $
  // >read6
  // AAAAAAAAAAAAAAAAAA#A
  // >read7
  // AAA

  length_seq = read_sequence_from_fasta(fp2,seq,max_read_length,true,&full_entry,0);
  
  CU_ASSERT_EQUAL(length_seq, 13);
  CU_ASSERT_STRING_EQUAL("read1",seq->name);
  CU_ASSERT(full_entry);

  length_seq = read_sequence_from_fasta(fp2,seq,max_read_length,true,&full_entry,0);

  CU_ASSERT_EQUAL(length_seq, 63);
  CU_ASSERT_STRING_EQUAL("read2",seq->name);
  CU_ASSERT(full_entry);

  length_seq = read_sequence_from_fasta(fp2,seq,max_read_length,true,&full_entry,0);
  
  CU_ASSERT_EQUAL(length_seq, 4);
  CU_ASSERT_STRING_EQUAL("read3",seq->name);
  CU_ASSERT_STRING_EQUAL("CCCC",seq->seq);
  CU_ASSERT(full_entry);

  length_seq = read_sequence_from_fasta(fp2,seq,max_read_length,true,&full_entry,0);
  
  CU_ASSERT_EQUAL(length_seq, 10);
  CU_ASSERT_STRING_EQUAL("read4",seq->name);
  CU_ASSERT_STRING_EQUAL("TTTTTTTTTT",seq->seq);
  CU_ASSERT(full_entry);

  length_seq = read_sequence_from_fasta(fp2,seq,max_read_length,true,&full_entry,0);
  
  CU_ASSERT_EQUAL(length_seq, 1);
  CU_ASSERT_STRING_EQUAL("read5",seq->name);
  CU_ASSERT(full_entry);

  length_seq = read_sequence_from_fasta(fp2,seq,max_read_length,true,&full_entry,0);
 
  CU_ASSERT_EQUAL(length_seq, 20);
  CU_ASSERT_STRING_EQUAL("read6",seq->name);
  CU_ASSERT(full_entry);

  length_seq = read_sequence_from_fasta(fp2,seq,max_read_length,true,&full_entry,0);
 
  CU_ASSERT_EQUAL(length_seq, 3);
  CU_ASSERT_STRING_EQUAL("read7",seq->name);
  CU_ASSERT_STRING_EQUAL("AAA",seq->seq);
  CU_ASSERT(full_entry);

  length_seq = read_sequence_from_fasta(fp2,seq,max_read_length,true,&full_entry,0);
  CU_ASSERT_EQUAL(length_seq, 0);
  CU_ASSERT(full_entry);
   
  fclose(fp2);


  //now make sure we do not get trapped in an infinite loop if the last read of a file is bad

  FILE* fp3= fopen("../data/test/basic/includes_final_read_that_has_bad_characters.fa", "r");

  // >read1
  // AAAAAAAAAAAA9
  // >read2
  // ¡€#9∞§¶#¶•#•#•#ª#ª#ª#ªº#º#º#º––––
  // >read3 4 c's
  // CCCC
  // >read4 10 Ts
  // TTTTTTTTTT
  // >read5
  // $
  // >read6
  // AAAAAAAAAAAAAAAAAA#A


  length_seq = read_sequence_from_fasta(fp3,seq,max_read_length,true,&full_entry,0);

  CU_ASSERT_EQUAL(length_seq, 13);
  CU_ASSERT_STRING_EQUAL("read1",seq->name);
  CU_ASSERT(full_entry);

  length_seq = read_sequence_from_fasta(fp3,seq,max_read_length,true,&full_entry,0);
  
  CU_ASSERT_EQUAL(length_seq, 63);
  CU_ASSERT_STRING_EQUAL("read2",seq->name);
  CU_ASSERT(full_entry);

  length_seq = read_sequence_from_fasta(fp3,seq,max_read_length,true,&full_entry,0);
  
  CU_ASSERT_EQUAL(length_seq, 4);
  CU_ASSERT_STRING_EQUAL("CCCC",seq->seq);
  CU_ASSERT_STRING_EQUAL("read3",seq->name);
  CU_ASSERT(full_entry);
  
  length_seq = read_sequence_from_fasta(fp3,seq,max_read_length,true,&full_entry,0);
  
  CU_ASSERT_EQUAL(length_seq, 10);
  CU_ASSERT_STRING_EQUAL("read4",seq->name);
  CU_ASSERT_STRING_EQUAL("TTTTTTTTTT",seq->seq);
  CU_ASSERT(full_entry);

  length_seq = read_sequence_from_fasta(fp3,seq,max_read_length,true,&full_entry,0);

  CU_ASSERT_EQUAL(length_seq, 1);
  CU_ASSERT_STRING_EQUAL("read5",seq->name);
  CU_ASSERT(full_entry);

  length_seq = read_sequence_from_fasta(fp3,seq,max_read_length,true,&full_entry,0);

  
  CU_ASSERT_EQUAL(length_seq, 20);
  CU_ASSERT_STRING_EQUAL("read6",seq->name);
  CU_ASSERT(full_entry);

  length_seq = read_sequence_from_fasta(fp3,seq,max_read_length,true,&full_entry,0);

  CU_ASSERT_EQUAL(length_seq, 0);
  CU_ASSERT(full_entry);

  fclose(fp3);



  free_sequence(&seq);

}


void test_shift_last_kmer_to_start_of_sequence(){ 
  Sequence * seq = malloc(sizeof(Sequence));
  boolean full_entry;

  if (seq == NULL){
    die("Out of memory trying to allocate Sequence\n");
  }
  //pre-allocate space where to read the sequences
  alloc_sequence(seq,200,LINE_MAX);

  FILE* fp1 = fopen("../data/test/basic/long_entries.fa", "r");

  int length_seq = read_sequence_from_fasta(fp1,seq,10,true,&full_entry,0);

  CU_ASSERT(seq->seq[0]=='A');
  CU_ASSERT(seq->seq[1]=='C');
  CU_ASSERT(seq->seq[2]=='G');
  
  shift_last_kmer_to_start_of_sequence(seq, length_seq,3);
  CU_ASSERT(seq->seq[0]=='T');
  CU_ASSERT(seq->seq[1]=='A');
  CU_ASSERT(seq->seq[2]=='C');
  
  CU_ASSERT(full_entry==false);
  fclose(fp1);
  free_sequence(&seq);
}

void test_read_sequence_from_fastq(){


  int ascii_offset = 33;

  //pre-allocate space where to read the sequences
  Sequence* seq = malloc(sizeof(Sequence));
  if (seq==NULL){
    die("Out of memory trying to allocate a Sequence");
  }
  
  alloc_sequence(seq,200,LINE_MAX);
  
  int length_seq;
  FILE* fp1 = fopen("../data/test/basic/one_entry.fq", "r");

  // 1. Read from simple fasta:
  // >Zam
  // ACGT
  // +
  // &&&&

  length_seq = read_sequence_from_fastq(fp1,seq,1000, ascii_offset);
  
  CU_ASSERT_EQUAL(length_seq, 4);
  CU_ASSERT_STRING_EQUAL("Zam",seq->name);
  CU_ASSERT_STRING_EQUAL("ACGT",seq->seq);
  //  CU_ASSERT_STRING_EQUAL("&&&&",seq->qual);/Zam says - /changed this line when I changed the ual reading code to offset by 33
  CU_ASSERT((int) seq->qual[0] == 5);




  FILE* fp2 = fopen("../data/test/basic/three_entries.fq", "r");
  
  //2. Read from fastq:


  // @Zam1
  //ACGT
  //+
  //&&&&
  //@Zam2
  //AAAAAAAA
  //+
  //!((/8F+,
  //@Zam3
  //ATATATAT
  //TTTTTTTTTT
  //-
  //(((((((+AAAAAABAAA

  
  length_seq = read_sequence_from_fastq(fp2,seq,1000, ascii_offset);
  
  CU_ASSERT_EQUAL(length_seq, 4);
  CU_ASSERT_STRING_EQUAL("Zam1",seq->name);
  CU_ASSERT_STRING_EQUAL("ACGT",seq->seq);
  //  CU_ASSERT_STRING_EQUAL("&&&&",seq->qual);/Zam says - /changed this line when I changed the ual reading code to offset by 33
  CU_ASSERT((int) seq->qual[0] == 5);


  length_seq = read_sequence_from_fastq(fp2,seq,1000, ascii_offset);
  
  CU_ASSERT_EQUAL(length_seq, 8);
  CU_ASSERT_STRING_EQUAL("Zam2",seq->name);
  CU_ASSERT_STRING_EQUAL("AAAAAAAA",seq->seq);
  CU_ASSERT((int) seq->qual[0] == 0);//! is quality 0 - take 33 off its ascii code, can check on http://www.asciitable.com/
  CU_ASSERT((int) seq->qual[1] == 7);// ( is quality 7
  CU_ASSERT((int) seq->qual[2] == 7);// ( is quality 7
  CU_ASSERT((int) seq->qual[3] == 14);// / is quality 14
  CU_ASSERT((int) seq->qual[4] == 23);// 8 is quality 23
  CU_ASSERT((int) seq->qual[5] == 37);// F is quality 37
  CU_ASSERT((int) seq->qual[6] == 10);// + is quality 10
  CU_ASSERT((int) seq->qual[7] == 11);// , is quality 11



  length_seq = read_sequence_from_fastq(fp2,seq,1000, ascii_offset);
  
  CU_ASSERT_EQUAL(length_seq, 18);
  CU_ASSERT_STRING_EQUAL("Zam3",seq->name);
  CU_ASSERT_STRING_EQUAL("ATATATATTTTTTTTTTT",seq->seq);
  CU_ASSERT((int) seq->qual[0] == 7);// ( is quality 7

  length_seq = read_sequence_from_fastq(fp2,seq,1000, ascii_offset);

  CU_ASSERT_EQUAL(length_seq, 0);

  fclose(fp2);
  free_sequence(&seq);
}




void test_read_sequence_from_fastq_with_bad_reads_and_long_reads()
{

  int ascii_offset=33;

  //pre-allocate space where to read the sequences
  Sequence* seq = malloc(sizeof(Sequence));
  if (seq==NULL){
    die("Out of memory trying to allocate a Sequence");
  }

  int max_read_length=200;
  alloc_sequence(seq,max_read_length,LINE_MAX);
  

  
  int length_seq;
  
  FILE* fp1 = fopen("../data/test/basic/includes_one_read_that_is_too_long.fq", "r");
  
  // @read1
  // ACGT
  // +
  // !!!!
  // @read2
  // CCCC
  // +
  // 5555
  // @read3
  // AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
  // -
  // 0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
  // @read4
  // ACGT
  // +
  // 3333



  length_seq = read_sequence_from_fastq(fp1,seq,max_read_length, ascii_offset);
  CU_ASSERT_EQUAL(length_seq, 4);
  CU_ASSERT_STRING_EQUAL("read1",seq->name);
  CU_ASSERT_STRING_EQUAL("ACGT",seq->seq);
  CU_ASSERT((int) (seq->qual[0])==0);
  CU_ASSERT((int) (seq->qual[1])==0);
  CU_ASSERT((int) (seq->qual[2])==0);
  CU_ASSERT((int) (seq->qual[3])==0);
  
  length_seq = read_sequence_from_fastq(fp1,seq,max_read_length, ascii_offset);
  CU_ASSERT_EQUAL(length_seq, 4);
  CU_ASSERT_STRING_EQUAL("read2",seq->name);
  CU_ASSERT_STRING_EQUAL("CCCC",seq->seq);
  CU_ASSERT((int) (seq->qual[0])==20);
  CU_ASSERT((int) (seq->qual[1])==20);
  CU_ASSERT((int) (seq->qual[2])==20);
  CU_ASSERT((int) (seq->qual[3])==20);

  
  length_seq = read_sequence_from_fastq(fp1,seq,max_read_length, ascii_offset);
  CU_ASSERT_EQUAL(length_seq, 100);
  CU_ASSERT_STRING_EQUAL("read3",seq->name);
  CU_ASSERT((int) (seq->qual[0])==15); // 0 translates as ascii 48; subtract 33 and get 15
  
  

  length_seq = read_sequence_from_fastq(fp1,seq,max_read_length, ascii_offset);
  CU_ASSERT_EQUAL(length_seq, 4);
  CU_ASSERT_STRING_EQUAL("read4",seq->name);
  CU_ASSERT_STRING_EQUAL("ACGT",seq->seq);
  CU_ASSERT((int) (seq->qual[0])==18);

  
  
  fclose(fp1);
  

  FILE* fp2 = fopen("../data/test/basic/includes_reads_with_bad_characters.fq", "r");

  //@read1
  //ACGTACGTACGTACGT
  //+
  //WEW2WEW2WEW2WEWA
  //@read2
  //AAAA#5A
  //+
  //1234567
  //@read3
  //TTTT
  //+
  //3333



  length_seq = read_sequence_from_fastq(fp2,seq,max_read_length, ascii_offset);
  CU_ASSERT_EQUAL(length_seq, 16);
  CU_ASSERT_STRING_EQUAL("read1",seq->name);
  CU_ASSERT_STRING_EQUAL("ACGTACGTACGTACGT",seq->seq);


  length_seq = read_sequence_from_fastq(fp2,seq,max_read_length, ascii_offset);
  CU_ASSERT_EQUAL(length_seq, 4);
  CU_ASSERT_STRING_EQUAL("read3",seq->name);
  CU_ASSERT_STRING_EQUAL("TTTT",seq->seq);

  
  length_seq = read_sequence_from_fastq(fp2,seq,max_read_length, ascii_offset);
  CU_ASSERT_EQUAL(length_seq, 0);
  
  fclose(fp2);




  FILE* fp3 = fopen("../data/test/basic/includes_one_read_where_quality_is_longer_than_seq.fq", "r");

  //@read1
  //ACGTACGTACGTACGT
  //+
  //WEW2WEW2WEW2WEWA
  //@read2
  //AAAA#5A
  //+
  //!!!!!!!!!!!!!!!!!!!!!!
  //@read3
  //TTTT
  //+
  //3333

  length_seq = read_sequence_from_fastq(fp3,seq,max_read_length, ascii_offset);
  CU_ASSERT_EQUAL(length_seq, 16);
  CU_ASSERT_STRING_EQUAL("read1",seq->name);
  CU_ASSERT_STRING_EQUAL("ACGTACGTACGTACGT",seq->seq);

  length_seq = read_sequence_from_fastq(fp3,seq,max_read_length, ascii_offset);
  CU_ASSERT_EQUAL(length_seq, 4);
  CU_ASSERT_STRING_EQUAL("read3",seq->name);
  CU_ASSERT_STRING_EQUAL("TTTT",seq->seq);

  length_seq = read_sequence_from_fastq(fp3,seq,max_read_length, ascii_offset);
  CU_ASSERT_EQUAL(length_seq, 0);
  
  fclose(fp3);

  FILE* fp4 = fopen("../data/test/basic/includes_multiline_reads.fq", "r");

  // @read1
  // ACGT
  // +
  // @@@@
  // @read2 45 bases
  // AAAAAAAAAAAAAAA
  // CCCCCCCCCCCCCCC
  // GGGGGGGGGGGGGGG
  // +
  // 222222222222222
  // 333333333333333
  // 444444444444444
  // @read3
  // TTT
  // -
  // ggg


  length_seq = read_sequence_from_fastq(fp4,seq,max_read_length, ascii_offset);
  CU_ASSERT_EQUAL(length_seq, 4);
  CU_ASSERT_STRING_EQUAL("read1",seq->name);
  CU_ASSERT_STRING_EQUAL("ACGT",seq->seq);
  

  length_seq = read_sequence_from_fastq(fp4,seq,max_read_length, ascii_offset);
  CU_ASSERT_EQUAL(length_seq, 45);
  CU_ASSERT_STRING_EQUAL("read2",seq->name);
  CU_ASSERT_STRING_EQUAL("AAAAAAAAAAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGG",seq->seq);
  

  length_seq = read_sequence_from_fastq(fp4,seq,max_read_length, ascii_offset);
  CU_ASSERT_EQUAL(length_seq, 3);
  CU_ASSERT_STRING_EQUAL("read3",seq->name);
  CU_ASSERT_STRING_EQUAL("TTT",seq->seq);
  


  length_seq = read_sequence_from_fastq(fp4,seq,max_read_length, ascii_offset);
  CU_ASSERT_EQUAL(length_seq, 0);
  
  fclose(fp4);

  
 

  free_sequence(&seq);






}

