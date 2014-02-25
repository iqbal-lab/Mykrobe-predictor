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
  test_dB_graph_population.c
*/
 
// system libraries
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <limits.h>

// third party libraries
#include <CUnit.h>
#include <Basic.h>

// cortex_var headers
#include "dB_graph.h"
#include "element.h"
#include "binary_kmer.h"
#include "file_reader.h"
#include "test_dB_graph_population.h"
#include "dB_graph_population.h"

// there are "pure" hash table tests which know nothing of the graph. This on the other hand
// is a sanity check one level up from those
void test_hash_table_find()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size = 4;
  int max_retries = 10;

  dBGraph * db_graph = hash_table_new(number_of_bits, bucket_size,
                                      max_retries, kmer_size);

  //Load the following fasta:
  //    >read1
  //    AAAAAAAA
  //    >read2
  //    GGCT
  //    >read3
  //    TAGG

  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/graph/test_dB_graph.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0, // 0 => falist/fqlist not colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  //length of sequence read in from file
  CU_ASSERT(seq_read == 16);
  //length of sequence that gets past filtwers and into the graph
  CU_ASSERT(seq_loaded==16);

  //number of kmers
  CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph), 5);

  //bad reads
  CU_ASSERT_EQUAL(bad_reads, 0);

  //all the kmers and their reverse complements from the reads
  BinaryKmer tmp_kmer1;
  BinaryKmer tmp_kmer2;
  
  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AAA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("TTT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element3 = hash_table_find(element_get_key(seq_to_binary_kmer("GGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element4 = hash_table_find(element_get_key(seq_to_binary_kmer("GCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element5 = hash_table_find(element_get_key(seq_to_binary_kmer("GCT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element6 = hash_table_find(element_get_key(seq_to_binary_kmer("AGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element7 = hash_table_find(element_get_key(seq_to_binary_kmer("TAG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element8 = hash_table_find(element_get_key(seq_to_binary_kmer("CTA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element9 = hash_table_find(element_get_key(seq_to_binary_kmer("AGG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element10 = hash_table_find(element_get_key(seq_to_binary_kmer("CCT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);

  //some kmers that should not be in the graph

  dBNode* test_element11 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element12 = hash_table_find(element_get_key(seq_to_binary_kmer("CCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element13 = hash_table_find(element_get_key(seq_to_binary_kmer("TAT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element14 = hash_table_find(element_get_key(seq_to_binary_kmer("ATA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element15 = hash_table_find(element_get_key(seq_to_binary_kmer("TAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element16 = hash_table_find(element_get_key(seq_to_binary_kmer("ATG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element17 = hash_table_find(element_get_key(seq_to_binary_kmer("TTG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element18 = hash_table_find(element_get_key(seq_to_binary_kmer("AAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element19 = hash_table_find(element_get_key(seq_to_binary_kmer("TGA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element20 = hash_table_find(element_get_key(seq_to_binary_kmer("TCA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);

  //kmers in the graph
  CU_ASSERT(test_element1 != NULL);
  CU_ASSERT(test_element2 != NULL);
  CU_ASSERT(test_element1 == test_element2);

  CU_ASSERT(test_element3 != NULL);
  CU_ASSERT(test_element4 != NULL);
  CU_ASSERT(test_element3 == test_element4);

  CU_ASSERT(test_element5 != NULL);
  CU_ASSERT(test_element6 != NULL);  
  CU_ASSERT(test_element5 == test_element6);

  CU_ASSERT(test_element7 != NULL);
  CU_ASSERT(test_element8 != NULL);
  CU_ASSERT(test_element7 == test_element8);

  CU_ASSERT(test_element9 != NULL);
  CU_ASSERT(test_element10 != NULL);
  CU_ASSERT(test_element9 == test_element10);  

  //kmers not in the graph
  CU_ASSERT(test_element11 == NULL);
  CU_ASSERT(test_element12 == NULL);
  CU_ASSERT(test_element13 == NULL);
  CU_ASSERT(test_element14 == NULL);
  CU_ASSERT(test_element15 == NULL);
  CU_ASSERT(test_element16 == NULL);
  CU_ASSERT(test_element17 == NULL);
  CU_ASSERT(test_element18 == NULL);
  CU_ASSERT(test_element19 == NULL);
  CU_ASSERT(test_element20 == NULL);

  
  hash_table_free(&db_graph);
  CU_ASSERT(db_graph == NULL);
}

void test_tip_clipping()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size = 4;

  Covg min_coverage, max_coverage;
  double avg_coverage;
  boolean is_cycle;
  BinaryKmer tmp_kmer1;
  BinaryKmer tmp_kmer2;

  dBGraph * db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);
  
  //Load the following fasta:
  
  //>main_trunk
  //GCGTCCCAT
  //>tip
  //CGTTT

  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/graph/generates_graph_with_tip.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT_EQUAL(seq_read,14);
  CU_ASSERT_EQUAL(seq_loaded,14);

  dBNode* node1 = hash_table_find(element_get_key(seq_to_binary_kmer("TTT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);

  CU_ASSERT(node1 != NULL);
  
  int tip_length = db_graph_db_node_clip_tip_for_specific_person_or_pop(node1, 10,&db_node_action_set_status_pruned,db_graph, 0);
  
  CU_ASSERT_EQUAL(tip_length,2);

  dBNode* node2 = hash_table_find(element_get_key(seq_to_binary_kmer("GTT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,
				  db_graph);
  dBNode* node3 = hash_table_find(element_get_key(seq_to_binary_kmer("CGT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,
				  db_graph);
  dBNode* node4 = hash_table_find(element_get_key(seq_to_binary_kmer("CCA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,
				  db_graph);


  //check the change in the status
  CU_ASSERT(db_node_check_status(node1,pruned));
  CU_ASSERT(db_node_check_status(node2,pruned));
  
  //check status didn't change
  CU_ASSERT(db_node_check_status(node3,none));
  CU_ASSERT(db_node_check_status(node4,none));

  //now check that if you look to see if the cipped kmers "exist" in the person's graph, you do not find them, as there is no edge
  CU_ASSERT(db_graph_find_node_restricted_to_specific_person_or_population(element_get_key(seq_to_binary_kmer("GTT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) , 
									   db_graph, 0)==NULL);
  CU_ASSERT(db_graph_find_node_restricted_to_specific_person_or_population(element_get_key(seq_to_binary_kmer("TTT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) , 
									   db_graph, 0)==NULL);

  //check tip clip works well with db_graph_supernode (ie that the tip was really removed)
  dBNode * nodes_path[100];
  Orientation orientations_path[100];
  Nucleotide labels_path[100];
  char tmp_seq[100+1];

  int length_supernode = db_graph_supernode_for_specific_person_or_pop(node3,100,&db_node_action_set_status_visited,
								       nodes_path,orientations_path,labels_path,
								       tmp_seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								       db_graph, 0);

  CU_ASSERT_EQUAL(length_supernode,6);
  CU_ASSERT_STRING_EQUAL("GGACGC",tmp_seq);
  
  
  //check ends
  dBNode* node5 = hash_table_find(element_get_key(seq_to_binary_kmer("GCG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,
				  db_graph);
  dBNode* node6 = hash_table_find(element_get_key(seq_to_binary_kmer("CAT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,
				  db_graph);

  //check the status are correct
  
  CU_ASSERT(db_node_check_status(node1,pruned));
  CU_ASSERT(db_node_check_status(node2,pruned));

  CU_ASSERT(db_node_check_status(node3,visited));  
  CU_ASSERT(db_node_check_status(node4,visited));
  CU_ASSERT(db_node_check_status(node5,visited));
  CU_ASSERT(db_node_check_status(node6,visited));

  hash_table_free(&db_graph);
  CU_ASSERT(db_graph == NULL);
}


void test_pruning_low_coverage_nodes()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size = 4;

  BinaryKmer tmp_kmer1;
  BinaryKmer tmp_kmer2;

  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  

  //Load the following fasta:
  
  //>main_trunk
  //GCGTCCCAT
  //>tip
  //CGTTT
  
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/graph/generates_graph_with_tip.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT_EQUAL(seq_read,14);
  CU_ASSERT_EQUAL(seq_loaded,14);

  dBNode* node1 = hash_table_find(element_get_key(seq_to_binary_kmer("CCC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* node2 = hash_table_find(element_get_key(seq_to_binary_kmer("TCC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* node3 = hash_table_find(element_get_key(seq_to_binary_kmer("CCA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  

  boolean node_pruned = db_graph_db_node_prune_low_coverage_ignoring_colours(node1,1,&db_node_action_set_status_pruned,db_graph);

  CU_ASSERT(node_pruned);
  CU_ASSERT(db_node_check_status(node1,pruned));

  CU_ASSERT(!db_node_edge_exist(node3,Guanine,reverse,0));
  CU_ASSERT(db_node_edge_exist(node3,Thymine,forward,0));

  CU_ASSERT(!db_node_edge_exist(node2,Cytosine,reverse,0));
  CU_ASSERT(db_node_edge_exist(node2,Cytosine,forward,0));

  //check all arrows were removed from node1
  void nucleotide_action(Nucleotide n){
    CU_ASSERT(!db_node_edge_exist(node1,n,reverse,0));
    CU_ASSERT(!db_node_edge_exist(node1,n,forward,0));
  }
  
  nucleotide_iterator(&nucleotide_action);
  hash_table_free(&db_graph);

  CU_ASSERT(db_graph == NULL);
}

void test_get_perfect_path_in_one_colour() 
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }
 
  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size = 10;

  dBNode * path_nodes[100];
  Orientation path_orientations[100];
  Nucleotide path_labels[100];
  char tmp_seq[100];

  double  avg_coverage;
  Covg min_coverage, max_coverage;

  dBGraph * db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);

 
  //1. Sequence of tests as follows
  //         Each test loads a single specifically designed fasta file into a dB_graph.
  //         The test then picks an element in the graph, and calls get_perfect_path
  //         and checks that it gets the right sequence.
  

  // ****
  //1.1 Fasta file that generate a graph with two hairpins, and a single edge (in each rorientation) joining them.
  //  Sequence is :  ACGTAC
  // ****

  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/graph/generates_graph_with_two_self_loops.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

   CU_ASSERT_EQUAL(seq_read,6);
   CU_ASSERT_EQUAL(seq_loaded,6);
   CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),2);
   CU_ASSERT_EQUAL(bad_reads,0);
   

  // now start at GTA and get all the sequence from there to the end of the supernode, and see
  // if that is right.

  BinaryKmer tmp_kmer1;
  BinaryKmer tmp_kmer2;

  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GTA", kmer_size,&tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
   CU_ASSERT(test_element1!=NULL);
   
   boolean is_cycle=false;

   
   int test1_length = db_graph_get_perfect_path_for_specific_person_or_pop(test_element1,forward,100,
									   &db_node_action_do_nothing,
									   path_nodes,path_orientations,path_labels,
									   tmp_seq, &avg_coverage, &min_coverage, &max_coverage,
									   &is_cycle,db_graph, 0);

    CU_ASSERT(is_cycle);
    CU_ASSERT_EQUAL(test1_length,4);
    
   CU_ASSERT_STRING_EQUAL(tmp_seq,"CGTA");

   hash_table_free(&db_graph);
   CU_ASSERT(db_graph == NULL);
 
 
   /*   // **** */
   /*   //1.2 Fasta file that generate a graph with one long supernode, with a conflict at the end */
   /*   //   caused by two outward/exiting edges */
   /*   // **** */
   
   //first set up the hash/graph
   kmer_size = 3;
   number_of_bits = 4;
   bucket_size = 10;

   
   db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);

  files_loaded = 0;
  bad_reads = 0;
  dup_reads = 0;
  seq_loaded = 0;
  seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/graph/generates_graph_with_one_long_supernode_with_conflict_at_end.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

   CU_ASSERT_EQUAL(seq_read,13);
   CU_ASSERT_EQUAL(seq_loaded,13);

   CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),5);
   CU_ASSERT_EQUAL(bad_reads,0);

   test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("ACA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
   CU_ASSERT(test_element1!=NULL);

   //ACA < TGT so forward gives TT
   is_cycle=false;
   

   //go forward
   int test2_length = db_graph_get_perfect_path_for_specific_person_or_pop(test_element1,forward,100,
						&db_node_action_do_nothing,
						path_nodes,path_orientations,path_labels,
						tmp_seq, &avg_coverage, &min_coverage, &max_coverage,
						&is_cycle,db_graph, 0);

   CU_ASSERT(!is_cycle);
   CU_ASSERT_EQUAL(test2_length,2);
   CU_ASSERT_STRING_EQUAL(tmp_seq,"TT");


   //ACA < TGT so backward gives ""
   
   //go backward
   test2_length = db_graph_get_perfect_path_for_specific_person_or_pop(test_element1,reverse,100,
					    &db_node_action_do_nothing,
					    path_nodes,path_orientations,path_labels,
					    tmp_seq, &avg_coverage, &min_coverage, &max_coverage,
					    &is_cycle,db_graph, 0);

  
   CU_ASSERT(!is_cycle);
   CU_ASSERT_EQUAL(test2_length,0);

   
   CU_ASSERT_STRING_EQUAL(tmp_seq,"");
   hash_table_free(&db_graph);
   CU_ASSERT(db_graph == NULL);


   // ****
   //1.3 Fasta file that generate a graph with one long supernode, with a conflict at the end
   //   caused by two INWARD edges in the opposite direction
   // ****
   
   //first set up the hash/graph
   kmer_size = 3;
   number_of_bits = 4;
   bucket_size = 5;
   bad_reads = 0;

   db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);

  files_loaded = 0;
  bad_reads = 0;
  dup_reads = 0;
  seq_loaded = 0;
  seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/graph/generates_graph_with_one_long_supernode_with_inward_conflict_at_end.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

   CU_ASSERT_EQUAL(seq_read,13);
   CU_ASSERT_EQUAL(seq_loaded,13);
   CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),5);
   CU_ASSERT_EQUAL(bad_reads,0);
   

   test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("ACA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
   CU_ASSERT(test_element1!=NULL);

   //ACA < TGT so forward gives TT
   is_cycle=false;

   //go forward
   int test3_length = db_graph_get_perfect_path_for_specific_person_or_pop(test_element1,forward,100,
						&db_node_action_do_nothing,
						path_nodes,path_orientations,path_labels,
						tmp_seq, &avg_coverage, &min_coverage, &max_coverage,
						&is_cycle,db_graph, 0);

   CU_ASSERT(!is_cycle);

   CU_ASSERT_EQUAL(test3_length,2);
   CU_ASSERT_STRING_EQUAL(tmp_seq,"TT");

   
   //ACA < TGT so backwar gives ""
   is_cycle=false;

   //go reverse
   test3_length = db_graph_get_perfect_path_for_specific_person_or_pop(test_element1,reverse,100,
					    &db_node_action_do_nothing,
					    path_nodes,path_orientations,path_labels,
					    tmp_seq, &avg_coverage, &min_coverage, &max_coverage,
					    &is_cycle,db_graph, 0);

   CU_ASSERT(!is_cycle);

   CU_ASSERT_EQUAL(test3_length,0);
   CU_ASSERT_STRING_EQUAL(tmp_seq,"");

   hash_table_free(&db_graph);
   CU_ASSERT(db_graph == NULL);


   // ****
   //1.4 Fasta file that generate a graph with an infinite loop at a single kmer
   //
   // ****

  
   //first set up the hash/graph
   kmer_size = 3;
   number_of_bits = 8;
   bucket_size = 4;

   db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);
   
  bad_reads = 0;
  dup_reads = 0;
  seq_loaded = 0;
  seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/graph/generates_graph_with_infinite_loop.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);


   CU_ASSERT_EQUAL(seq_read,25);
   CU_ASSERT_EQUAL(seq_loaded,25);
   CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),1);
   CU_ASSERT_EQUAL(bad_reads,0);

   test_element1 = hash_table_find(seq_to_binary_kmer("AAA", kmer_size, &tmp_kmer1) ,db_graph);
   CU_ASSERT(test_element1!=NULL);

   //forward
   is_cycle=false;

   int test4_length = db_graph_get_perfect_path_for_specific_person_or_pop(test_element1,forward,100,
						&db_node_action_do_nothing,
						path_nodes,path_orientations,path_labels,
						tmp_seq, &avg_coverage, &min_coverage, &max_coverage,
						&is_cycle,db_graph, 0);
   CU_ASSERT(is_cycle);//this time it should find a cycle
   CU_ASSERT_STRING_EQUAL(tmp_seq,"A");
   CU_ASSERT_EQUAL(test4_length,1);

   //backward
   is_cycle=false;

   test4_length = db_graph_get_perfect_path_for_specific_person_or_pop(test_element1,reverse,100,
					    &db_node_action_do_nothing,
					    path_nodes,path_orientations,path_labels,
					    tmp_seq, &avg_coverage, &min_coverage, &max_coverage,
					    &is_cycle,db_graph, 0);

   CU_ASSERT(is_cycle);//this time it should find a cycle
   CU_ASSERT_EQUAL(test4_length,1);
   CU_ASSERT_STRING_EQUAL(tmp_seq,"T");

   hash_table_free(&db_graph);
   CU_ASSERT(db_graph == NULL);


   
   // ****
   // 1.5 check parameters (path nodes,labels,etc) for get_perfect_path
   //
   // ****
   
   kmer_size = 3;
   number_of_bits = 4;
   bucket_size = 4;
   db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
   
   boolean found1, found2, found3;
   dBNode * node1;
   dBNode * node2;
   dBNode * node3;
   dBNode * node4;

   
   node1 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CGT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2), &found1, db_graph);
   node2 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("GTT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found2, db_graph);
   db_node_add_edge(node1, node2, reverse,reverse, db_graph->kmer_size, 0);

   node3 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("TTA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found3, db_graph);
   db_node_add_edge(node2, node3, reverse,reverse, db_graph->kmer_size, 0);
 
   node4 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("TAG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found3, db_graph);
   db_node_add_edge(node3, node4, reverse,reverse, db_graph->kmer_size, 0);
  

   //add coverage
   db_node_update_coverage(node1, 0, 10);
   db_node_update_coverage(node2, 0,2);
   db_node_update_coverage(node3, 0,1);
   db_node_update_coverage(node4, 0,20);
  


   dBNode * nodes[10];
   Orientation orientations[10];
   Nucleotide bases[10];
   int test5_length = 0;
 

   test5_length = db_graph_get_perfect_path_for_specific_person_or_pop(node1,reverse, 10,
					    &db_node_action_set_status_visited,
					    nodes, orientations, bases,
					    tmp_seq,&avg_coverage,&min_coverage,&max_coverage,
					    &is_cycle,db_graph, 0);

  CU_ASSERT_EQUAL(test5_length,3);
  
  //check nodes
  CU_ASSERT_EQUAL(node1,nodes[0]);
  CU_ASSERT_EQUAL(node2,nodes[1]);
  CU_ASSERT_EQUAL(node3,nodes[2]);
  CU_ASSERT_EQUAL(node4,nodes[3]);

  //check labels
  CU_ASSERT_EQUAL(bases[0],Thymine);
  CU_ASSERT_EQUAL(bases[1],Adenine);
  CU_ASSERT_EQUAL(bases[2],Guanine);
  
  //check orientations
  CU_ASSERT_EQUAL(orientations[0],reverse);
  CU_ASSERT_EQUAL(orientations[1],reverse);
  CU_ASSERT_EQUAL(orientations[2],reverse);

  //check statuses
  CU_ASSERT(db_node_check_status(nodes[0], none));
  CU_ASSERT(db_node_check_status(nodes[1], visited));
  CU_ASSERT(db_node_check_status(nodes[2], visited));
  CU_ASSERT(db_node_check_status(nodes[3], none));

  //check coverage
  CU_ASSERT_EQUAL(avg_coverage,1.5);
  CU_ASSERT_EQUAL(min_coverage,1);
  CU_ASSERT_EQUAL(max_coverage,2);
  
  hash_table_free(&db_graph);
  CU_ASSERT(db_graph == NULL);
}



void test_detect_and_smooth_bubble()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  int kmer_size = 3;
  int number_of_buckets = 4;
  int bucket_size = 3;
  dBGraph * db_graph = hash_table_new(number_of_buckets,bucket_size,10,kmer_size);
  boolean found;
  BinaryKmer tmp_kmer1;
  BinaryKmer tmp_kmer2;

 
  dBNode * node1, * node2, * node3, * node4, * node5, * node6, * node7, * node8, * node9;
  Orientation orientations1[100],orientations2[100];
  Nucleotide labels1[100], labels2[100];
  dBNode * path_nodes1[100], * path_nodes2[100];

  int length1, length2;
  char seq1[101], seq2[101];

  boolean bubble;

  double avg_coverage1,avg_coverage2;
  Covg min_coverage1, max_coverage1, min_coverage2, max_coverage2;
  
  //check a perfect bubble 

  //start point
  node1 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2), &found, db_graph);

  //branch1
  node2 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);


  //modify coverage -- used below to clip one branch
  db_node_update_coverage(node2, 0, 5);
  CU_ASSERT_EQUAL(db_node_get_coverage(node2, 0),5);


  node3 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CAC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
  node4 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("ACT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);

  //end point (branch merge here)
  node5 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CTT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);

  //branch 2
  node6 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);

  db_node_update_coverage(node6, 0, 1);
  CU_ASSERT_EQUAL(db_node_get_coverage(node6, 0),1);

  node7 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CTC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
  node8 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("TCT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);

  //add 3p extension
  node9 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("TTA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
 
  //branch1
  db_node_add_edge(node1, node2, forward,forward, db_graph->kmer_size, 0);
  db_node_add_edge(node2, node3, forward,forward, db_graph->kmer_size, 0);
  db_node_add_edge(node3, node4, forward,forward, db_graph->kmer_size, 0);
  db_node_add_edge(node4, node5, forward,reverse, db_graph->kmer_size, 0);

  
  //branch2
  db_node_add_edge(node1, node6, forward,reverse, db_graph->kmer_size, 0);
  db_node_add_edge(node6, node7, reverse,forward, db_graph->kmer_size, 0);
  db_node_add_edge(node7, node8, forward,reverse, db_graph->kmer_size, 0);
  db_node_add_edge(node8, node5, reverse,reverse, db_graph->kmer_size, 0);

  CU_ASSERT(db_node_edge_exist(node1,Thymine,forward, 0));
  CU_ASSERT(db_node_edge_exist(node1,Adenine,forward, 0));

  //add 3p extension
  db_node_add_edge(node5, node9, reverse,reverse, db_graph->kmer_size, 0);


  bubble = db_graph_detect_bubble_for_specific_person_or_population(node1,forward,
								    db_graph->kmer_size+2,
								    &db_node_action_set_status_visited,
								    &length1,
								    path_nodes1,
								    orientations1,
								    labels1,
								    seq1,
								    &avg_coverage1,
								    &min_coverage1,
								    &max_coverage1,
								    &length2,
								    path_nodes2,
								    orientations2,
								    labels2,
								    seq2,
								    &avg_coverage2,
								    &min_coverage2,
								    &max_coverage2,
								    db_graph, 
								     
								    0);


  CU_ASSERT(bubble);

  CU_ASSERT(
	    ((labels1[0] == Adenine) && (labels2[0] == Thymine)) 
	    || 
	    ((labels2[0] == Adenine) && (labels1[0] == Thymine))
	    );

  CU_ASSERT_EQUAL(labels1[1],Cytosine);
  CU_ASSERT_EQUAL(labels1[2],Thymine);
  CU_ASSERT_EQUAL(labels1[3],Thymine);

  CU_ASSERT_EQUAL(labels2[1],Cytosine);
  CU_ASSERT_EQUAL(labels2[2],Thymine);
  CU_ASSERT_EQUAL(labels2[3],Thymine);



  //check extensions
  
  Orientation orientations5p[100];
  Orientation orientations3p[100];
  Nucleotide labels_flank5p[100];
  Nucleotide labels_flank3p[100];
  boolean is_cycle5p, is_cycle3p;
  int length_flank5p,length_flank3p;
  dBNode * nodes5p[100];
  dBNode * nodes3p[100];
  char tmp_seq[101];

  double avg_coverage;
  Covg min_coverage, max_coverage;
  
  length_flank5p = db_graph_get_perfect_path_for_specific_person_or_pop(node1, reverse,100,
									&db_node_action_set_status_visited,
									nodes5p,orientations5p, labels_flank5p,
									tmp_seq,&avg_coverage,&min_coverage,&max_coverage,
									&is_cycle5p,db_graph, 0);

  CU_ASSERT_EQUAL(length_flank5p,0);
  CU_ASSERT_EQUAL(nodes5p[0],node1);


  length_flank3p = db_graph_get_perfect_path_for_specific_person_or_pop(node5, reverse,100,
									&db_node_action_set_status_visited,
									nodes3p,orientations3p, labels_flank3p,
									tmp_seq,&avg_coverage,&min_coverage,&max_coverage,
									&is_cycle3p,db_graph, 0);
 
  
  CU_ASSERT_EQUAL(length_flank3p,1);
  CU_ASSERT_EQUAL(nodes3p[0],node5);
  CU_ASSERT_EQUAL(nodes3p[1],node9);
  CU_ASSERT_EQUAL(labels_flank3p[0],Adenine);

  //check statuses
  CU_ASSERT(db_node_check_status(node1,none));
  CU_ASSERT(db_node_check_status(node2,visited));
  CU_ASSERT(db_node_check_status(node3,visited));
  CU_ASSERT(db_node_check_status(node4,visited));
  CU_ASSERT(db_node_check_status(node6,visited));
  CU_ASSERT(db_node_check_status(node7,visited));
  CU_ASSERT(db_node_check_status(node8,visited));
  //this last node shouldn't be marked as visited
  CU_ASSERT(db_node_check_status(node5,none));


  //check the bubble from node5 perspective
  bubble = db_graph_detect_bubble_for_specific_person_or_population(node5,forward,db_graph->kmer_size+1,
								    &db_node_action_set_status_visited,
								    &length1,path_nodes1,orientations1,labels1,
								    seq1,&avg_coverage1,&min_coverage1,&max_coverage1,
								    &length2,path_nodes2,orientations2,labels2,
								    seq2,&avg_coverage2,&min_coverage2,&max_coverage2,
								    db_graph,   0);
  

  CU_ASSERT(bubble);
  CU_ASSERT(((labels1[0] == Adenine) && (labels2[0] == Thymine)) || (((labels2[0] == Adenine) && (labels1[0] == Thymine))));
  

  //test smooth bubble

  //remove a branch
  //this one should fail -- because of limit 0 on min coverage
  boolean removed = db_graph_db_node_smooth_bubble_for_specific_person_or_pop(node1,forward,db_graph->kmer_size+1,1,
									      &db_node_action_set_status_pruned,db_graph,   0);
  CU_ASSERT(!removed);


  //this one should work
  removed = db_graph_db_node_smooth_bubble_for_specific_person_or_pop(node1,forward,db_graph->kmer_size+1,10,
								      &db_node_action_set_status_pruned,db_graph,   0);
  
  CU_ASSERT(removed);
  //check that correct path was removed

  //the first node shouldn't be marked
  CU_ASSERT(db_node_check_status(node1,none));
  CU_ASSERT(db_node_check_status(node2,visited));
  CU_ASSERT(db_node_check_status(node3,visited));
  CU_ASSERT(db_node_check_status(node4,visited));
  CU_ASSERT(db_node_check_status(node6,pruned));
  CU_ASSERT(db_node_check_status(node7,pruned));
  CU_ASSERT(db_node_check_status(node8,pruned));
  //the last node shouldn't be marked 
  CU_ASSERT(db_node_check_status(node5,none));

  CU_ASSERT(!db_node_edge_exist(node1,Thymine,forward, 0));
  CU_ASSERT(db_node_edge_exist(node1,Adenine,forward, 0));

  //check arrows by getting a supernode

  dBNode * nodes_path[100];
  Orientation orientations_path[100];
  Nucleotide labels_path[100];
  boolean is_cycle;

  db_graph_supernode_for_specific_person_or_pop(node2,100,&db_node_action_set_status_visited,
						nodes_path,orientations_path,labels_path,
						tmp_seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
						db_graph, 0);
  
  CU_ASSERT_STRING_EQUAL(tmp_seq,"ACTTA");
  
  
  db_graph_supernode_for_specific_person_or_pop(node7,100,&db_node_action_set_status_visited,
						nodes_path,orientations_path,labels_path,
						tmp_seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
						db_graph, 0);
  
  CU_ASSERT_STRING_EQUAL(tmp_seq,"");
  

  hash_table_free(&db_graph);
  CU_ASSERT(db_graph == NULL);

  //test bubbles with unequal branch sizes
   
  kmer_size = 21;
  number_of_buckets = 4; //2**4 buckets
  bucket_size = 40;
  db_graph = hash_table_new(number_of_buckets,bucket_size,10,kmer_size);

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/graph/generate_bubble_with_unequal_branch_sizes.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT_EQUAL(seq_read,343);

  //fetch the kmer where the path splits (bubble appears) GCCAACCATGCCTGTTAAGGG
  node1 = hash_table_find(element_get_key(seq_to_binary_kmer("GCCAACCATGCCTGTTAAGGG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
 
  CU_ASSERT(node1!=NULL);

  bubble = db_graph_detect_bubble_for_specific_person_or_population(node1,reverse,100,
								    &db_node_action_set_status_visited,
								    &length1,path_nodes1,orientations1,labels1,
								    seq1,&avg_coverage1,&min_coverage1,&max_coverage1,
								    &length2,path_nodes2,orientations2,labels2,
								    seq2,&avg_coverage2,&min_coverage2,&max_coverage2,
								    db_graph, 0);
  

  CU_ASSERT(bubble);
  
  CU_ASSERT((length1 == 21 && length2==22) || (length2 == 21 && length1==22));

  if (length1==21)
    {
      CU_ASSERT_EQUAL(labels1[0],Guanine);
      CU_ASSERT_STRING_EQUAL(seq1,"GGGTTATTTTTCTAGAGAGTT");
      CU_ASSERT_EQUAL(labels2[0],Thymine);
      CU_ASSERT_STRING_EQUAL(seq2,"TGGGTTATTTTTCTAGAGAGTT");
    }
  else{
    CU_ASSERT_EQUAL(labels2[0],Guanine);
    CU_ASSERT_STRING_EQUAL(seq2,"GGGTTATTTTTCTAGAGAGTT");
    CU_ASSERT_EQUAL(labels1[0],Thymine);
    CU_ASSERT_STRING_EQUAL(seq1,"TGGGTTATTTTTCTAGAGAGTT");
  }
    
  hash_table_free(&db_graph);
  CU_ASSERT(db_graph == NULL);

}



void test_db_graph_db_node_has_precisely_n_edges_with_status_in_one_colour()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size = 4;
  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);

  boolean found;
  dBNode * node1, * node2, * node3, * node4;
  dBNode  * next_node[4];
  Orientation next_orientation[4];
  Nucleotide next_base[4];

  BinaryKmer tmp_kmer1;
  BinaryKmer tmp_kmer2;

  node1 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2), &found, db_graph);
  node2 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2), &found, db_graph);
  node3 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2), &found, db_graph);

 

  db_node_add_edge(node1, node2, forward,reverse, db_graph->kmer_size, 0);
  db_node_add_edge(node1, node3, forward,forward, db_graph->kmer_size, 0);

  boolean one_edge1 = db_graph_db_node_has_precisely_n_edges_with_status_for_specific_person_or_pop(node1,forward,none,1,
												    next_node,next_orientation,next_base,db_graph, 0);


  CU_ASSERT_EQUAL(one_edge1,false);

  db_node_set_status(node2,visited);

  boolean one_edge2 = db_graph_db_node_has_precisely_n_edges_with_status_for_specific_person_or_pop(node1,forward,none,1,
												    next_node,next_orientation,next_base,db_graph, 0);

 
  CU_ASSERT_EQUAL(one_edge2,true);
  CU_ASSERT_EQUAL(node3,next_node[0]);
  CU_ASSERT_EQUAL(next_orientation[0],forward);
  CU_ASSERT_EQUAL(next_base[0],Guanine);
  

  node4 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2), &found, db_graph);
  db_node_add_edge(node1, node4, forward, forward, db_graph->kmer_size, 0);
   
  boolean one_edge3 = db_graph_db_node_has_precisely_n_edges_with_status_for_specific_person_or_pop(node1,forward,none,2,
												    next_node,next_orientation,next_base,db_graph, 0);
  

  
  CU_ASSERT_EQUAL(one_edge3,true);

  //observe that the tests below require to know in wich order the edges are visited

  CU_ASSERT_EQUAL(node4,next_node[0]);
  CU_ASSERT_EQUAL(node3,next_node[1]);

  CU_ASSERT_EQUAL(next_base[0],Adenine);
  CU_ASSERT_EQUAL(next_base[1],Guanine);
  
  CU_ASSERT_EQUAL(next_orientation[0],forward);
  CU_ASSERT_EQUAL(next_orientation[1],forward);

  hash_table_free(&db_graph);
  CU_ASSERT(db_graph == NULL);
}


// DEV: is this code needed? Zam: It needs fixing - it's functionality I want to work - both th main function and this test
/* - ready to be uncommented when N50 cde is upgraded to use db_graph_supernode
void test_get_N50()
{

  //first set up the hash/graph
  int kmer_size = 5;
  int number_of_bits = 5;
  int bucket_size    = 4;
  int max_retries = 10;

  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  int seq_length=0;
  //long long count_kmers = 0;
  long long bad_reads=0; long long dup_reads=0;
  boolean remove_duplicates_single_endedly=false; 
  boolean break_homopolymers=false;
  int homopolymer_cutoff=0;


  //1. One fasta file containing two reads:
  //  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA which becomes just a one node supernode
  //  CACGGTATTATTTACCT which is a 13 node supernode.
  // so the N50 is 13
  
   int max_read_length=70;

   seq_length = load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/graph/n50_example1.fa",&bad_reads, &dup_reads, max_read_length,  
										   remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff, 
										   db_graph, 0);
   

   CU_ASSERT_EQUAL(seq_length,58);
   //CU_ASSERT_EQUAL(count_kmers,14);
   CU_ASSERT_EQUAL(bad_reads,0);

   CU_ASSERT(db_graph_get_N50_of_supernodes(db_graph, 0)==13);

   //clean-up
   hash_table_free(&db_graph);



   // ==========================================================================
   // example 2:

   // >read1
   // AAAAAA
   // >read2
   // CCCCCC
   // >read3
   // ACGTA
   // >read4
   // CTATG
   // >read5
   // TTTAT
   // >read6
   // GCTTA
   // >read7
   // AGGCT
   // >read8
   // CACGGTATT
   // 7 singleton super´nodes, and one 5-node supernode. So N50 is 1
   
   //first set up the hash/graph
   db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);


   //count_kmers=0;
   bad_reads=0;
   seq_length=0;
   seq_length = load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/graph/n50_example2.fa",&bad_reads, &dup_reads, max_read_length,  
										   remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff, 
										   db_graph, 0);

   CU_ASSERT_EQUAL(seq_length,46);
   //  CU_ASSERT_EQUAL(count_kmers,12);
   CU_ASSERT_EQUAL(bad_reads,0);


   CU_ASSERT(db_graph_get_N50_of_supernodes(db_graph, 0)==1); //remember this leaves the nodes all visited
 
   //clean-up
   hash_table_free(&db_graph);

}
*/



void test_is_condition_true_for_all_nodes_in_supernode()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size = 10;

  BinaryKmer tmp_kmer1;
  BinaryKmer tmp_kmer2;

  dBGraph * db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);

 
  //1. Sequence of tests as follows
  //         Each test loads a single specifically designed fasta file into a dB_graph.
  

  // ****
  //1.1 Fasta file that generate a graph with two hairpins, and a single edge (in each rorientation) joining them.
  //  Sequence is :  ACGTAC
  // ****

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/graph/generates_graph_with_two_self_loops.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT_EQUAL(seq_read,6);
  CU_ASSERT_EQUAL(seq_loaded,6);
  CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),2);
  CU_ASSERT_EQUAL(bad_reads,0);


  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GTA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(test_element1!=NULL);
  dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("ACG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(test_element2!=NULL);

  int limit = 50;
  dBNode * nodes_path[limit];
  Orientation orientations_path[limit];
  Nucleotide labels_path[limit];
  char seq[limit+1];
  int length_path;
  double avg_coverage;
  boolean is_cycle;
  Covg min_coverage, max_coverage;

  //  printf("Check if all nodes have status none - this should be true\n");
  CU_ASSERT(db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,						   
								  &db_node_check_status_none,
								  &db_node_action_do_nothing, 
								  nodes_path, orientations_path, labels_path, &length_path,
								  seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								  db_graph, 0)==true);
  CU_ASSERT(db_node_check_status(test_element1, none));
  //  printf("set status of one node to visited\n");
  db_node_set_status(test_element1, visited);
  CU_ASSERT(db_node_check_status(test_element1, visited));
  CU_ASSERT(db_node_check_status(test_element2, none));
  // printf("Checkif all nodes have status none - this should not be true - confirm this\n");
  CU_ASSERT(!db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,
								   &db_node_check_status_none,
								   &db_node_action_do_nothing, 
								   nodes_path, orientations_path, labels_path, &length_path, 
								   seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								   db_graph, 0));


  //  printf("double check statuses of nodes unchanged\n");
  CU_ASSERT(db_node_check_status(test_element1, visited));
  CU_ASSERT(db_node_check_status(test_element2, none));

  // printf("Check if all nodes have status visited - this should not be true - confirm this\n");
  CU_ASSERT(!db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,
								   &db_node_check_status_visited,
								   &db_node_action_do_nothing, 
								   nodes_path, orientations_path, labels_path, &length_path,
								   seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								   db_graph, 0));


  //check again, but this time, set all nodes to visited in the process
  //printf("Set all nodes to visited while checking them all\n");
  CU_ASSERT(!db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50, 
								   &db_node_check_status_visited,
								   &db_node_action_set_status_visited_or_visited_and_exists_in_reference, 
								   nodes_path, orientations_path, labels_path, &length_path, 
								   seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								   db_graph, 0));
  //and now this time it SHOULD BE TRUE
  CU_ASSERT(db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50, 
								  &db_node_check_status_visited,
								  &db_node_action_do_nothing, 
								  nodes_path, orientations_path, labels_path, &length_path, 
								  seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								  db_graph, 0));

  // and should still be true
  CU_ASSERT(db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,
								  &db_node_check_status_visited,
								  &db_node_action_do_nothing, 
								  nodes_path, orientations_path, labels_path, &length_path, 
								  seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								  db_graph, 0));
  

  //printf("Nodes currently all visited. Set all nodes to none while checking them all to confirm that they are currently all visited (before doing the set-to-none)\n");
  CU_ASSERT(db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,
								  &db_node_check_status_visited,
								  &db_node_action_set_status_none, 
								  nodes_path, orientations_path, labels_path, &length_path,
								  seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								  db_graph, 0));

  CU_ASSERT(db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,
								  &db_node_check_status_none,
								  &db_node_action_do_nothing, 
								  nodes_path, orientations_path, labels_path, &length_path, 
								  seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								  db_graph, 0));
  
  db_node_set_status(test_element2, exists_in_reference);
  CU_ASSERT(!db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,
								   &db_node_check_status_none,
								   &db_node_action_do_nothing, 
								   nodes_path, orientations_path, labels_path, &length_path, 
								   seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								   db_graph, 0));

  CU_ASSERT(!db_graph_is_condition_true_for_all_nodes_in_supernode(test_element2, 50,
								   &db_node_check_status_none,
								   &db_node_action_do_nothing, 
								   nodes_path, orientations_path, labels_path, &length_path,
								   seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								   db_graph, 0));

   
  hash_table_free(&db_graph);
  
  
}


void test_db_graph_supernode_for_specific_person_or_pop()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }
  if(NUMBER_OF_COLOURS < 2)
  {
    warn("Test needs >= 2 colours\n");
    return;
  }

  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 8;
  int bucket_size = 8;
  int max_retries = 10;
  
  dBGraph * hash_table = hash_table_new(number_of_bits, bucket_size,
                                        max_retries, kmer_size);

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/supernode/one_person_one_long_supernode_with_conflict_at_end.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT(seq_loaded==13);
  CU_ASSERT(bad_reads==0);

  //we have loaded a fasta containing these supernodes ACATT, TTC, TTG
  // reads are: ACATT, ATTC, ATTG
  //note we have loaded one person only
 

  int max_expected_supernode_length=100;
  dBNode * nodes_path[max_expected_supernode_length];
  Orientation orientations_path[max_expected_supernode_length];
  Nucleotide labels_path[max_expected_supernode_length];
  char seq[max_expected_supernode_length+hash_table->kmer_size+1];
  BinaryKmer tmp_kmer1, tmp_kmer2;

  //get the supernodes one at a time. Remember length is the number of edges between nodes.

  dBNode* test_elem1 = hash_table_find(element_get_key(seq_to_binary_kmer("ACA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table);
  dBNode* test_elem2 = hash_table_find(element_get_key(seq_to_binary_kmer("CAT", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table);
  dBNode* test_elem3 = hash_table_find(element_get_key(seq_to_binary_kmer("ATT", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table);
  CU_ASSERT(!(test_elem1 == NULL));
  CU_ASSERT(!(test_elem2 == NULL));
  CU_ASSERT(!(test_elem3 == NULL));


  CU_ASSERT(db_node_check_status(test_elem1,none));
  CU_ASSERT(db_node_check_status(test_elem2,none));
  CU_ASSERT(db_node_check_status(test_elem3,none));

  //coverage variables
  double avg_coverage=0;
  Covg max_coverage = 0;
  Covg min_coverage = 0;
  boolean is_cycle=false;

  int length = db_graph_supernode_for_specific_person_or_pop(test_elem1,max_expected_supernode_length, &db_node_action_set_status_visited_or_visited_and_exists_in_reference,
							     nodes_path,orientations_path, labels_path, seq,
							     &avg_coverage, &min_coverage, &max_coverage, &is_cycle,
							     hash_table, 0); //0 because we are looking at person 0 - the only person in the graph

  CU_ASSERT(length==2);
  CU_ASSERT(avg_coverage==1);
  CU_ASSERT(min_coverage==1);
  CU_ASSERT(max_coverage==1); //remember coverage only measured on internal nodes of a supernode.
  CU_ASSERT(is_cycle==false);

  CU_ASSERT_STRING_EQUAL(seq, "TT"); //does not put initial kmer "ACA" into the string, just the edges
  CU_ASSERT(nodes_path[0]==test_elem1);
  CU_ASSERT(nodes_path[1]==test_elem2);
  CU_ASSERT(nodes_path[2]==test_elem3);
  CU_ASSERT(orientations_path[0]==forward);
  CU_ASSERT(orientations_path[1]==reverse);
  CU_ASSERT(orientations_path[2]==reverse);
  CU_ASSERT(labels_path[0]==Thymine);
  CU_ASSERT(labels_path[1]==Thymine);
	     
  CU_ASSERT(db_node_check_status(test_elem1, visited));
  CU_ASSERT(db_node_check_status(test_elem2, visited));
  CU_ASSERT(db_node_check_status(test_elem3, visited));


  hash_table_free(&hash_table);


  // ******** try another example ******************
  

  hash_table = hash_table_new(number_of_bits, bucket_size,
                              max_retries, kmer_size);

  bad_reads = 0;
  dup_reads = 0;
  seq_loaded = 0;
  seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/test_pop_load_and_print/two_individuals_simple.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  //person 1:
  //>person1_read1
  //AAAAAAAAAAAAAAA
  //>person1_read2
  //ACGTT     <<<<<<<<<<<<<<<<Note this generates supernode (AAC)GTT due to a hairpin
  //person 2:
  //>person2_read1
  //GGGGGGGGGGGGGGGGGG
  //>person2_read2
  //TTGACG   <<<<<<<<<<<<<<<< Also has a hairpin

  test_elem1 = hash_table_find(element_get_key(seq_to_binary_kmer("AAA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table);
  test_elem2 = hash_table_find(element_get_key(seq_to_binary_kmer("ACG", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table);
  test_elem3 = hash_table_find(element_get_key(seq_to_binary_kmer("CGT", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table);
  dBNode* test_elem4 = hash_table_find(element_get_key(seq_to_binary_kmer("GTT", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table);
  dBNode* test_elem5 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table);

  dBNode* test_elem6 = hash_table_find(element_get_key(seq_to_binary_kmer("TTG", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table);

  dBNode* test_elem7 = hash_table_find(element_get_key(seq_to_binary_kmer("TGA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table);
  dBNode* test_elem8 = hash_table_find(element_get_key(seq_to_binary_kmer("GAC", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table);
  dBNode* test_elem9 = hash_table_find(element_get_key(seq_to_binary_kmer("ACG", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table);


  CU_ASSERT(!(test_elem1 == NULL));
  CU_ASSERT(!(test_elem2 == NULL));
  CU_ASSERT(!(test_elem3 == NULL));
  CU_ASSERT(!(test_elem4 == NULL));
  CU_ASSERT(!(test_elem5 == NULL));
  CU_ASSERT(!(test_elem6 == NULL));
  CU_ASSERT(!(test_elem7 == NULL));
  CU_ASSERT(!(test_elem8 == NULL));
  CU_ASSERT(!(test_elem9 == NULL));

  avg_coverage=0;
  min_coverage=0;
  max_coverage=0;
  is_cycle=false;

  length = db_graph_supernode_for_specific_person_or_pop(test_elem1,max_expected_supernode_length, &db_node_action_set_status_visited_or_visited_and_exists_in_reference,
							 nodes_path,orientations_path, labels_path, seq,
							 &avg_coverage, &min_coverage, &max_coverage, &is_cycle,
							 hash_table, 0); //person 0 int his array is person in person1.fa - sorry, offset by 1


  CU_ASSERT(length==1); //just one edge - self loop
  CU_ASSERT_STRING_EQUAL(seq, "A");
  CU_ASSERT(nodes_path[0]==test_elem1);
  CU_ASSERT(nodes_path[1]==test_elem1);
  CU_ASSERT(orientations_path[0]==forward);
  CU_ASSERT(orientations_path[1]==forward);
  CU_ASSERT(labels_path[0]==Adenine);

  CU_ASSERT(avg_coverage==0);
  CU_ASSERT(min_coverage==0);
  CU_ASSERT(max_coverage==0);//there are not interior nodes of supernode
  CU_ASSERT(is_cycle==true);

  CU_ASSERT(db_node_check_status(test_elem1, visited)==true); 
  CU_ASSERT(db_node_check_status(test_elem2, none)==true);
  CU_ASSERT(db_node_check_status(test_elem3, none)==true);
  CU_ASSERT(db_node_check_status(test_elem4, none)==true);
  CU_ASSERT(db_node_check_status(test_elem5, none)==true);
  CU_ASSERT(db_node_check_status(test_elem6, none)==true);
  CU_ASSERT(db_node_check_status(test_elem7, none)==true);
  CU_ASSERT(db_node_check_status(test_elem8, none)==true);
  CU_ASSERT(db_node_check_status(test_elem9, none)==true);

  db_graph_set_all_visited_nodes_to_status_none(hash_table);
  CU_ASSERT(db_node_check_status(test_elem1, none)==true);


  avg_coverage=0;
  min_coverage=0;
  max_coverage=0;
  is_cycle=false;

  length = db_graph_supernode_for_specific_person_or_pop(test_elem2,max_expected_supernode_length,&db_node_action_set_status_visited_or_visited_and_exists_in_reference,
							 nodes_path,orientations_path, labels_path, seq,
							 &avg_coverage, &min_coverage, &max_coverage, &is_cycle,
							 hash_table, 0); //person 0 int his array is person in person1.fa - sorry, offset by 1


  CU_ASSERT(length==3);
  CU_ASSERT_STRING_EQUAL(seq, "GTT");
  CU_ASSERT(nodes_path[0]==test_elem4);
  CU_ASSERT(nodes_path[1]==test_elem2);
  CU_ASSERT(nodes_path[2]==test_elem2);//hairpin!
  CU_ASSERT(nodes_path[3]==test_elem4);
  CU_ASSERT(orientations_path[0]==forward);
  CU_ASSERT(orientations_path[1]==forward);
  CU_ASSERT(orientations_path[2]==reverse);
  CU_ASSERT(orientations_path[3]==reverse);
  CU_ASSERT(labels_path[0]==Guanine);
  CU_ASSERT(labels_path[1]==Thymine);
  CU_ASSERT(labels_path[2]==Thymine);
  CU_ASSERT(avg_coverage==2);
  CU_ASSERT(min_coverage==2);
  CU_ASSERT(max_coverage==2);//note that the single interior node ALSO exists in the graph of person2, but that should not incrememtn this coverage to 3.
  CU_ASSERT(is_cycle==false);
  
  CU_ASSERT(db_node_check_status(test_elem1, none)==true);
  CU_ASSERT(db_node_check_status(test_elem2, visited)==true);
  CU_ASSERT(db_node_check_status(test_elem3, visited)==true);
  CU_ASSERT(db_node_check_status(test_elem4, visited)==true); 
  CU_ASSERT(db_node_check_status(test_elem5, none)==true);
  CU_ASSERT(db_node_check_status(test_elem6, none)==true);
  CU_ASSERT(db_node_check_status(test_elem7, none)==true);
  CU_ASSERT(db_node_check_status(test_elem8, none)==true);
  CU_ASSERT(db_node_check_status(test_elem9, visited)==true);

  
  //now try again with a different initial node within the same supernode. Should return the same as before.
  //In general it doesn't guarantee to trvaerse the supernode in the same direction; we know which was it will go however

  avg_coverage=0;
  min_coverage=0;
  max_coverage=0;
  is_cycle=false;

  length = db_graph_supernode_for_specific_person_or_pop(test_elem3,max_expected_supernode_length,&db_node_action_set_status_visited_or_visited_and_exists_in_reference,
							 nodes_path,orientations_path, labels_path, seq,
							 &avg_coverage, &min_coverage, &max_coverage, &is_cycle,
							 hash_table, 0); //person 0 int his array is person in person1.fa - sorry, offset by 1
  

  CU_ASSERT(length==3);
  CU_ASSERT_STRING_EQUAL(seq, "GTT");
  CU_ASSERT(nodes_path[0]==test_elem4);
  CU_ASSERT(nodes_path[1]==test_elem2);
  CU_ASSERT(nodes_path[2]==test_elem2);//hairpin!
  CU_ASSERT(nodes_path[3]==test_elem4);
  CU_ASSERT(orientations_path[0]==forward);
  CU_ASSERT(orientations_path[1]==forward);
  CU_ASSERT(orientations_path[2]==reverse);
  CU_ASSERT(orientations_path[3]==reverse);
  CU_ASSERT(labels_path[0]==Guanine);
  CU_ASSERT(labels_path[1]==Thymine);
  CU_ASSERT(labels_path[2]==Thymine);
  CU_ASSERT(avg_coverage==2);
  CU_ASSERT(min_coverage==2);
  CU_ASSERT(max_coverage==2);//note that the single interior node ALSO exists in the graph of person2, but that should not incrememtn this coverage to 3.
  CU_ASSERT(is_cycle==false);
  
  CU_ASSERT(db_node_check_status(test_elem1, none)==true);
  CU_ASSERT(db_node_check_status(test_elem2, visited)==true);
  CU_ASSERT(db_node_check_status(test_elem3, visited)==true);
  CU_ASSERT(db_node_check_status(test_elem4, visited)==true); 
  CU_ASSERT(db_node_check_status(test_elem5, none)==true);
  CU_ASSERT(db_node_check_status(test_elem6, none)==true);
  CU_ASSERT(db_node_check_status(test_elem7, none)==true);
  CU_ASSERT(db_node_check_status(test_elem8, none)==true);
  CU_ASSERT(db_node_check_status(test_elem9, visited)==true);


  length = db_graph_supernode_for_specific_person_or_pop(test_elem4,max_expected_supernode_length,&db_node_action_set_status_visited_or_visited_and_exists_in_reference,
							 nodes_path,orientations_path, labels_path, seq,
							 &avg_coverage, &min_coverage, &max_coverage, &is_cycle,
							 hash_table, 0); //person 0 int his array is person in person1.fa - sorry, offset by 1
  
  
  CU_ASSERT(length==3);
  CU_ASSERT_STRING_EQUAL(seq, "GTT");
  CU_ASSERT(nodes_path[0]==test_elem4);
  CU_ASSERT(nodes_path[1]==test_elem2);
  CU_ASSERT(nodes_path[2]==test_elem2);//hairpin!
  CU_ASSERT(nodes_path[3]==test_elem4);
  CU_ASSERT(orientations_path[0]==forward);
  CU_ASSERT(orientations_path[1]==forward);
  CU_ASSERT(orientations_path[2]==reverse);
  CU_ASSERT(orientations_path[3]==reverse);
  CU_ASSERT(labels_path[0]==Guanine);
  CU_ASSERT(labels_path[1]==Thymine);
  CU_ASSERT(labels_path[2]==Thymine);
  CU_ASSERT(avg_coverage==2);
  CU_ASSERT(min_coverage==2);
  CU_ASSERT(max_coverage==2);//note that the single interior node ALSO exists in the graph of person2, but that should not incrememtn this coverage to 3.
  CU_ASSERT(is_cycle==false);
  
  CU_ASSERT(db_node_check_status(test_elem1, none)==true);
  CU_ASSERT(db_node_check_status(test_elem2, visited)==true);
  CU_ASSERT(db_node_check_status(test_elem3, visited)==true);
  CU_ASSERT(db_node_check_status(test_elem4, visited)==true); 
  CU_ASSERT(db_node_check_status(test_elem5, none)==true);
  CU_ASSERT(db_node_check_status(test_elem6, none)==true);
  CU_ASSERT(db_node_check_status(test_elem7, none)==true);
  CU_ASSERT(db_node_check_status(test_elem8, none)==true);
  CU_ASSERT(db_node_check_status(test_elem9, visited)==true);

  //and now a trap! This node is not there in person 1, but is in person 2
  length = db_graph_supernode_for_specific_person_or_pop(test_elem6,max_expected_supernode_length,&db_node_action_set_status_visited_or_visited_and_exists_in_reference,
							 nodes_path,orientations_path, labels_path, seq,
							 &avg_coverage, &min_coverage, &max_coverage, &is_cycle,
							 hash_table, 0);   
  CU_ASSERT(length==0);



  length = db_graph_supernode_for_specific_person_or_pop(test_elem6,max_expected_supernode_length,&db_node_action_set_status_visited_or_visited_and_exists_in_reference,
							 nodes_path,orientations_path, labels_path, seq,
							 &avg_coverage, &min_coverage, &max_coverage, &is_cycle,
							 hash_table,  1); 
  CU_ASSERT(length==3);
  CU_ASSERT_STRING_EQUAL(seq, "CAA");  

  hash_table_free(&hash_table);
}


void test_is_supernode_end()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 8;
  int bucket_size = 8;
  int max_retries = 10;
  BinaryKmer tmp_kmer1, tmp_kmer2;

  dBGraph * hash_table = hash_table_new(number_of_bits, bucket_size,
                                        max_retries, kmer_size);

 
  //1. Sequence of tests as follows
  //         Each test loads a single specifically designed fasta file into a graph/hash table.
  //         The test then picks an element in the graph, and calls get_seq_from_elem_to_end_of_supernode
  //         and checks that it gets the right sequence.
  

  // ****
  //1.1 Fasta file that generate a graph with two hairpins, and a single edge (in each rorientation) joining them.
  //  Sequence is :  ACGTAC
  // ****

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/supernode/one_person_two_self_loops.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT(seq_loaded==6);
  CU_ASSERT(seq_read==6);
  CU_ASSERT(seq_read==6);
  CU_ASSERT(bad_reads==0);

  //GTA is not the end of a supernode in either direction
  dBNode* query_node = hash_table_find(element_get_key(seq_to_binary_kmer("GTA",hash_table->kmer_size, &tmp_kmer1), hash_table->kmer_size, &tmp_kmer2), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node, forward, 0, hash_table));
  CU_ASSERT(!db_node_is_supernode_end(query_node, reverse, 0, hash_table));

  //CGT is not the end of a supernode in either direction
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("CGT",hash_table->kmer_size, &tmp_kmer1), hash_table->kmer_size, &tmp_kmer2), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node, forward, 0, hash_table));
  CU_ASSERT(!db_node_is_supernode_end(query_node, reverse, 0, hash_table));

  //ACG is not  the end of a supernode in the reverse direction
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("ACG",hash_table->kmer_size, &tmp_kmer1), hash_table->kmer_size, &tmp_kmer2), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node, forward, 0, hash_table));
  CU_ASSERT(!db_node_is_supernode_end(query_node, reverse, 0, hash_table));

  hash_table_free(&hash_table);


  // ****
  //1.2 Fasta file that generate a graph with one long supernode, with a conflict at the end
  //   caused by two outward/exiting edges 
  // ****


  //first set up the hash/graph
  kmer_size = 3;
  number_of_bits = 4;
  bucket_size    = 4;
  bad_reads = 0;
  max_retries=10;

  hash_table = hash_table_new(number_of_bits, bucket_size,
                              max_retries, kmer_size);

  bad_reads = 0;
  dup_reads = 0;
  seq_loaded = 0;
  seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/supernode/one_person_one_long_supernode_with_conflict_at_end.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT(seq_loaded==13);
  CU_ASSERT(seq_read==13);
  CU_ASSERT(bad_reads==0);

  //ACA IS  the end of a supernode in the reverse direction
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("ACA",hash_table->kmer_size, &tmp_kmer1), hash_table->kmer_size, &tmp_kmer2), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node, forward, 0, hash_table));
  CU_ASSERT(db_node_is_supernode_end(query_node, reverse, 0, hash_table));
  

  //ATG is not a supernode end
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("ATG",hash_table->kmer_size, &tmp_kmer1), hash_table->kmer_size, &tmp_kmer2), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node, forward, 0, hash_table));
  CU_ASSERT(!db_node_is_supernode_end(query_node, reverse, 0, hash_table));


  //ATT IS  the end of a supernode in the reverse direction
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("ATT",hash_table->kmer_size, &tmp_kmer1), hash_table->kmer_size, &tmp_kmer2), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node, forward, 0, hash_table));
  CU_ASSERT(db_node_is_supernode_end(query_node, reverse, 0, hash_table));
  
  //TTC is a supernode end in both directions
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("TTC",hash_table->kmer_size, &tmp_kmer1), hash_table->kmer_size, &tmp_kmer2), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(db_node_is_supernode_end(query_node, forward, 0, hash_table));
  CU_ASSERT(db_node_is_supernode_end(query_node, reverse, 0, hash_table));

  //TTG is a supernode end in both directions
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("TTG",hash_table->kmer_size, &tmp_kmer1), hash_table->kmer_size, &tmp_kmer2), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(db_node_is_supernode_end(query_node, forward, 0, hash_table));
  CU_ASSERT(db_node_is_supernode_end(query_node, reverse, 0, hash_table));

  hash_table_free(&hash_table);

  // ****
  //1.3 Fasta file that generate a graph with one long supernode, with a conflict at the end
  //   caused by two INWARD edges in the opposite direction
  // ****

  //first set up the hash/graph
  kmer_size = 3;
  number_of_bits = 4;
  bucket_size = 4;
  bad_reads = 0;
  max_retries = 10;

  hash_table = hash_table_new(number_of_bits, bucket_size,
                              max_retries, kmer_size);

  bad_reads = 0;
  dup_reads = 0;
  seq_loaded = 0;
  seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/supernode/one_person_one_long_supernode_with_inward_conflict_at_end.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT(seq_loaded==13);
  CU_ASSERT(seq_read==13);
  CU_ASSERT(bad_reads==0);

  //ACA IS  the end of a supernode in the reverse direction
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("ACA",hash_table->kmer_size, &tmp_kmer1), hash_table->kmer_size, &tmp_kmer2), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node, forward, 0, hash_table));
  CU_ASSERT(db_node_is_supernode_end(query_node, reverse, 0, hash_table));
  
  //ATG is not a supernode end
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("ATG",hash_table->kmer_size, &tmp_kmer1), hash_table->kmer_size, &tmp_kmer2), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node, forward, 0, hash_table));
  CU_ASSERT(!db_node_is_supernode_end(query_node, reverse, 0, hash_table));


  //ATT IS  the end of a supernode in the reverse direction
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("ATT",hash_table->kmer_size, &tmp_kmer1), hash_table->kmer_size, &tmp_kmer2), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node, forward, 0, hash_table));
  CU_ASSERT(db_node_is_supernode_end(query_node, reverse, 0, hash_table));
  
  //TTC is a supernode end in both directions
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("TTC",hash_table->kmer_size, &tmp_kmer1), hash_table->kmer_size, &tmp_kmer2), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(db_node_is_supernode_end(query_node, forward, 0, hash_table));
  CU_ASSERT(db_node_is_supernode_end(query_node, reverse, 0, hash_table));

  //AAT is a supernode end in reverse directions
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("AAT",hash_table->kmer_size, &tmp_kmer1), hash_table->kmer_size, &tmp_kmer2), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node, forward, 0, hash_table));
  CU_ASSERT(db_node_is_supernode_end(query_node, reverse, 0, hash_table));

  //TAA is a supernode end in both directions
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("TAA",hash_table->kmer_size, &tmp_kmer1), hash_table->kmer_size, &tmp_kmer2), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(db_node_is_supernode_end(query_node, forward, 0, hash_table));
  CU_ASSERT(db_node_is_supernode_end(query_node, reverse, 0, hash_table));

  hash_table_free(&hash_table);

  // ****
  //1.4 Fasta file that generate a graph with an infinite loop at a single kmer
  //   
  // ****

  
  //first set up the hash/graph
  kmer_size = 3;
  number_of_bits = 4;
  bucket_size    = 4;
  bad_reads = 0;
  max_retries=10;

  hash_table = hash_table_new(number_of_bits, bucket_size,
                              max_retries, kmer_size);

  bad_reads = 0;
  dup_reads = 0;
  seq_loaded = 0;
  seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/supernode/one_person_infiniteloop.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT(bad_reads==0);
  CU_ASSERT(seq_loaded==25);
  CU_ASSERT(seq_read==25);

  //AAA is a supernode end in both directions
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("AAA",hash_table->kmer_size, &tmp_kmer1), hash_table->kmer_size, &tmp_kmer2), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(db_node_is_supernode_end(query_node, forward, 0, hash_table));
  CU_ASSERT(db_node_is_supernode_end(query_node, reverse, 0, hash_table));

  hash_table_free(&hash_table);

}


void test_getting_stats_of_how_many_indivduals_share_a_node()
{ 
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  //       >person1_read1
  //        AAAAAAAAAAAAAAA
  //        >person1_read2
  //        ACGTT
  //
  //
  //       >person2_read1
  //       GGGGGGGGGGGGGGGGGG
  //       >person2_read2
  //       TTGACG


  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size = 4;
  int max_retries = 10;

  dBGraph * hash_table = hash_table_new(number_of_bits, bucket_size,
                                        max_retries, kmer_size);
  
  if (hash_table==NULL)
    {
      die("unable to alloc the hash table. dead before we even started. OOM");
    }

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/test_pop_load_and_print/two_individuals_simple.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  //printf("Number of bases loaded is %d",seq_loaded);
  
  CU_ASSERT(seq_loaded == 44);
  CU_ASSERT(bad_reads==0);

  int stats[6]; //array of ints. stats[0] - number of kmers shared by noone; stats[1], number shared by one person, etc
  stats[0]=0; stats[1]=0; stats[2]=0; stats[3]=0; stats[4]=0; stats[5]=0; 

  int* array[3];
  array[0]=&stats[0];
  array[1]=&stats[1];
  array[2]=&stats[2];

  int num_people=2;

  db_graph_traverse_to_gather_statistics_about_people(&find_out_how_many_individuals_share_this_node_and_add_to_statistics , hash_table, array, num_people);

  CU_ASSERT(stats[0]==0);
  CU_ASSERT(stats[1]==6);
  CU_ASSERT(stats[2]==1);

  hash_table_free(&hash_table);
 
}


void test_get_min_and_max_covg_of_nodes_in_supernode()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size = 4;
  int max_retries = 10;
  BinaryKmer tmp_kmer1, tmp_kmer2;

  dBGraph * hash_table = hash_table_new(number_of_bits, bucket_size,
                                        max_retries, kmer_size);

  
  if (hash_table==NULL)
    {
      die("unable to alloc the hash table. dead before we even started. OOM");
    }


  //start with an example with just oner person - note this is an example with the same kmer occuring twice in one read :) - CCG
  // so need to make sure that coverga eof that kmer is not double-incremented for the two times in that read

  // >r1
  //AACCGG
  //>r2
  //AAC
  //>r3
  //AAC
  //>r4
  //AAC
  //>r5
  //ACC
  //>r6
  //CCG
  //>r7
  //CCG

  //supernode will be
  // AAC --> ACC --> CCG ___|
  // TTG <-- TGG <-- GGC <--|
  // hairpin :-)

  //double coverness of hairpin must not confuse read coverage

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/coverage/one_person.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  //printf("Number of bases loaded is %d",seq_loaded);
  CU_ASSERT(seq_loaded == 24);
  CU_ASSERT(seq_read == 24);
  CU_ASSERT(bad_reads==0);

  dBNode* query_node = hash_table_find(element_get_key(seq_to_binary_kmer("AAC",hash_table->kmer_size, &tmp_kmer1), hash_table->kmer_size, &tmp_kmer2), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(db_node_get_coverage(query_node, 0)==4);
  //printf("Expect covg aac is 4, but is %d", db_node_get_coverage(query_node, 0));
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("ACC",hash_table->kmer_size, &tmp_kmer1), hash_table->kmer_size, &tmp_kmer2), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(db_node_get_coverage(query_node, 0)==2);
  //printf("Expect covg acc is 2, but is %d", db_node_get_coverage(query_node, 0));
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("CCG",hash_table->kmer_size, &tmp_kmer1), hash_table->kmer_size, &tmp_kmer2), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(db_node_get_coverage(query_node, 0)==4); //NOTE this is only 4 because we cound kmer coverage. CCG occurs twice in read r1, once on each strand
  //printf("Expect covg ccg is 4, but is %d", db_node_get_coverage(query_node, 0));

  hash_table_free(&hash_table);
}


void test_db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  int length_of_arrays=14;
  int number_of_nodes_to_load=7;


  //first set up the hash/graph
  int kmer_size = 7;
  int number_of_bits = 8;
  int bucket_size = 10;
  int max_retries = 10;

  dBGraph * db_graph = hash_table_new(number_of_bits, bucket_size,
                                      max_retries, kmer_size);

  if (db_graph==NULL)
    {
      die("unable to alloc the hash table. dead before we even started. OOM");
    }

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/one_person_for_testing_array_loading.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  //>one read
  //AATAGACGCCCACACCTGATAGACCCCACACTCTAA

  
  CU_ASSERT(seq_loaded==36);
  
  FILE* chrom_fptr = fopen("../data/test/pop_graph/one_person.fa", "r");
  if (chrom_fptr==NULL)
    {
      die("Cannot open ./data/test/pop_graph/one_person.fa\n");
    }

  //*************************************
  // malloc and initialising
  //*************************************
  dBNode**     chrom_path_array        = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays); //everything these dBNode*'s point to will be in the hash table - ie owned by the hash table
  Orientation* chrom_orientation_array = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays); 
  Nucleotide*  chrom_labels         = (Nucleotide*) malloc(sizeof(Nucleotide)*length_of_arrays);
  char*        chrom_string            = (char*) malloc(sizeof(char)*length_of_arrays+1); //+1 for \0

  int n;
  for (n=0; n<length_of_arrays; n++)
    {
      chrom_path_array[n]=NULL;
      chrom_orientation_array[n]=forward;
      chrom_labels[n]=Undefined;
      chrom_string[n]='N';
    }
  chrom_string[0]='\0';


  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence\n");
  }
  alloc_sequence(seq,number_of_nodes_to_load+db_graph->kmer_size+1,LINE_MAX);
  seq->seq[0]='\0';


  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      die("Failed to malloc kmer sliding window . Exit.\n");
    }
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*1000);    //*(number_of_nodes_to_load + db_graph->kmer_size));
  if (kmer_window->kmer==NULL)
    {
      die("Failed to malloc kmer_window->kmer");
    }
  
  kmer_window->nkmers=0;

  // ***********************************************************
  //********** end of malloc and initialise ********************



  int ret = db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta(chrom_fptr, number_of_nodes_to_load, 0, 
													     length_of_arrays,
													     chrom_path_array, chrom_orientation_array, chrom_labels, chrom_string,
													     seq, kmer_window, 
													     true, false,
													     db_graph);
  CU_ASSERT(ret==number_of_nodes_to_load);

  char tmp_seqzam[db_graph->kmer_size+1];
  tmp_seqzam[db_graph->kmer_size]='\0';


  int i;
  //length of arrays is 14, number of nodes to load is 7
  for(i=0; i<length_of_arrays-number_of_nodes_to_load; i++)
    {
      CU_ASSERT(chrom_path_array[i]==NULL);
    }
  CU_ASSERT_STRING_EQUAL( "AATAGAC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[7]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ATAGACG", binary_kmer_to_seq(element_get_kmer(chrom_path_array[8]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "GCGTCTA", binary_kmer_to_seq(element_get_kmer(chrom_path_array[9]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "AGACGCC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[10]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "GACGCCC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[11]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACGCCCA", binary_kmer_to_seq(element_get_kmer(chrom_path_array[12]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CGCCCAC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[13]), db_graph->kmer_size, tmp_seqzam) );


  //one more batch, then array is full,
  ret = db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta(chrom_fptr, number_of_nodes_to_load, number_of_nodes_to_load, 
													 length_of_arrays,
													 chrom_path_array, chrom_orientation_array, chrom_labels, chrom_string,
													 seq, kmer_window, 
													 false, false,
													 db_graph);
  CU_ASSERT(ret==number_of_nodes_to_load);



  CU_ASSERT_STRING_EQUAL( "AATAGAC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[0]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ATAGACG", binary_kmer_to_seq(element_get_kmer(chrom_path_array[1]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "GCGTCTA", binary_kmer_to_seq(element_get_kmer(chrom_path_array[2]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "AGACGCC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[3]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "GACGCCC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[4]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACGCCCA", binary_kmer_to_seq(element_get_kmer(chrom_path_array[5]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CGCCCAC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[6]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "GCCCACA", binary_kmer_to_seq(element_get_kmer(chrom_path_array[7]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CCCACAC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[8]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CCACACC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[9]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "AGGTGTG", binary_kmer_to_seq(element_get_kmer(chrom_path_array[10]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACACCTG", binary_kmer_to_seq(element_get_kmer(chrom_path_array[11]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CACCTGA", binary_kmer_to_seq(element_get_kmer(chrom_path_array[12]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACCTGAT", binary_kmer_to_seq(element_get_kmer(chrom_path_array[13]), db_graph->kmer_size, tmp_seqzam) );
  

  //from now on, it is always true that LAST time was not a new fasta entry, so penultimate argument is TRUE
  ret = db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta(chrom_fptr, number_of_nodes_to_load, number_of_nodes_to_load, 
													 length_of_arrays,
													 chrom_path_array, chrom_orientation_array, chrom_labels, chrom_string,
													 seq, kmer_window, 
													 false, true,
                                                                                                         db_graph);

  CU_ASSERT(ret==number_of_nodes_to_load);


  CU_ASSERT_STRING_EQUAL( "GCCCACA", binary_kmer_to_seq(element_get_kmer(chrom_path_array[0]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CCCACAC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[1]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CCACACC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[2]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "AGGTGTG", binary_kmer_to_seq(element_get_kmer(chrom_path_array[3]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACACCTG", binary_kmer_to_seq(element_get_kmer(chrom_path_array[4]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CACCTGA", binary_kmer_to_seq(element_get_kmer(chrom_path_array[5]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACCTGAT", binary_kmer_to_seq(element_get_kmer(chrom_path_array[6]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CCTGATA", binary_kmer_to_seq(element_get_kmer(chrom_path_array[7]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CTATCAG", binary_kmer_to_seq(element_get_kmer(chrom_path_array[8]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "TCTATCA", binary_kmer_to_seq(element_get_kmer(chrom_path_array[9]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "GATAGAC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[10]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ATAGACC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[11]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "GGGTCTA", binary_kmer_to_seq(element_get_kmer(chrom_path_array[12]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "AGACCCC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[13]), db_graph->kmer_size, tmp_seqzam) );
  

  //and again

  ret = db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta(chrom_fptr, number_of_nodes_to_load, number_of_nodes_to_load, 
													 length_of_arrays,
													 chrom_path_array, chrom_orientation_array, chrom_labels, chrom_string,
													 seq, kmer_window, 
													 false, true,
                                                                                                         db_graph);

  CU_ASSERT(ret==number_of_nodes_to_load);


  CU_ASSERT_STRING_EQUAL( "CCTGATA", binary_kmer_to_seq(element_get_kmer(chrom_path_array[0]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CTATCAG", binary_kmer_to_seq(element_get_kmer(chrom_path_array[1]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "TCTATCA", binary_kmer_to_seq(element_get_kmer(chrom_path_array[2]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "GATAGAC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[3]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ATAGACC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[4]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "GGGTCTA", binary_kmer_to_seq(element_get_kmer(chrom_path_array[5]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "AGACCCC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[6]), db_graph->kmer_size, tmp_seqzam) );

  CU_ASSERT_STRING_EQUAL( "GACCCCA", binary_kmer_to_seq(element_get_kmer(chrom_path_array[7]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACCCCAC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[8]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CCCCACA", binary_kmer_to_seq(element_get_kmer(chrom_path_array[9]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CCCACAC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[10]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "AGTGTGG", binary_kmer_to_seq(element_get_kmer(chrom_path_array[11]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CACACTC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[12]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACACTCT", binary_kmer_to_seq(element_get_kmer(chrom_path_array[13]), db_graph->kmer_size, tmp_seqzam) );


  //now - does it cope with hitting the end of the entry before getting the required number of nodes

  ret = db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta(chrom_fptr, number_of_nodes_to_load, number_of_nodes_to_load, 
													 length_of_arrays,
													 chrom_path_array, chrom_orientation_array, chrom_labels, chrom_string,
													 seq, kmer_window, 
													 false, true,
                                                                                                         db_graph);

  CU_ASSERT(ret==2);

  CU_ASSERT_STRING_EQUAL( "GACCCCA", binary_kmer_to_seq(element_get_kmer(chrom_path_array[0]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACCCCAC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[1]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CCCCACA", binary_kmer_to_seq(element_get_kmer(chrom_path_array[2]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CCCACAC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[3]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "AGTGTGG", binary_kmer_to_seq(element_get_kmer(chrom_path_array[4]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CACACTC", binary_kmer_to_seq(element_get_kmer(chrom_path_array[5]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACACTCT", binary_kmer_to_seq(element_get_kmer(chrom_path_array[6]), db_graph->kmer_size, tmp_seqzam) );

  CU_ASSERT_STRING_EQUAL( "CACTCTA", binary_kmer_to_seq(element_get_kmer(chrom_path_array[7]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACTCTAA", binary_kmer_to_seq(element_get_kmer(chrom_path_array[8]), db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT(chrom_path_array[9]==NULL);
  CU_ASSERT(chrom_path_array[10]==NULL);
  CU_ASSERT(chrom_path_array[11]==NULL);
  CU_ASSERT(chrom_path_array[12]==NULL);
  CU_ASSERT(chrom_path_array[13]==NULL);

  //cleanup

  free(kmer_window->kmer);
  free(kmer_window);
  free(chrom_path_array);
  free(chrom_orientation_array);
  free(chrom_string);
  free(chrom_labels);
  free_sequence(&seq);
  hash_table_free(&db_graph);
  fclose(chrom_fptr);

}


// A fundamental assumption in our implementation is that the fasta whose path
// you are following is long, much longer than the longest supernode you are
// going to find. This is fine when running using reference chromosomes. However
// in tests, this is a pain, and so I am using fasta's for these tests that are
// padded at the end with N's.

void test_db_graph_make_reference_path_based_sv_calls_null_test_1()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // ===========================================================================
  // 1. NULL test. Load short fake reference, and load a person whose sequence
  //    is identical. This should find nothing.
  // ===========================================================================


  //first set up the hash/graph
  int kmer_size = 7;
  int number_of_bits = 8;
  int bucket_size = 10;
  int max_retries = 10;

  dBGraph * hash_table = hash_table_new(number_of_bits, bucket_size,
                                        max_retries, kmer_size);

  if (hash_table==NULL)
    {
      die("unable to alloc the hash table. dead before we even started. OOM");
    }

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/one_person_with_Ns_on_end.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  //>one read
  //AATAGACGCCCACACCTGATAGACCCCACAC 

  // The de Bruijn graph of this in 7mers is as follows:
  // a line of 8 nodes, with a loop of length 16 off the 9th node. ie all nodes
  // have 1-in and 1-out, except the 9th node, CCCACAC, which has 2 ins and
  // 2 outs. So it's a supernode of length 8, a loop-supernode of length 16,
  // our algorithm should look at the supernode starting with AATAGAC, then the
  // supernode starting at CCCACAC, then stop

  
  CU_ASSERT(seq_read==1462);
  CU_ASSERT(seq_loaded==31);

  FILE* chrom_fptr = fopen("../data/test/pop_graph/first_person_with_one_read_and_Ns_on_end.fa", "r");
  if (chrom_fptr==NULL)
    {
      die("Cannot open ./data/test/pop_graph/first_person_with_one_read_and_Ns_on_end.fa\n");
    }

  int min_fiveprime_flank_anchor = 2;
  int min_threeprime_flank_anchor= 3;// I want to try to attach anchors, and the supernodes are quite short, so for this test only ask for (ridiculously) small anchors
  int max_anchor_span = 20;
  int length_of_arrays=40;
  int min_covg =1;
  int max_covg = 10;
  int max_expected_size_of_supernode=20;

  
  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, 0, 
							1,
							min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
							max_expected_size_of_supernode, length_of_arrays, hash_table, NULL,
							0, NULL, NULL, NULL, NULL, NULL, &make_reference_path_based_sv_calls_condition_always_true, &action_set_flanks_and_branches_to_be_ignored,
							&print_no_extra_info, NULL, NoIdeaWhatCleaning);

  CU_ASSERT(ret==0);

  hash_table_free(&hash_table);
  fclose(chrom_fptr);

}



void test_db_graph_make_reference_path_based_sv_calls_null_test_2()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // ===========================================================================
  // 2. Harder NULL test. Reference=an ALU and load a person whose sequence is
  //    also that ALU. This should find nothing. Our implementation of
  //    db_graph_make_reference_path_based_sv_calls assumes we have a very large
  //    ref fasta, so we load big arrays of bases which are much longer than we
  //    expect any of the supernodes to be. In this case, we know the supernode
  //    is as long as the Alu.  
  // ===========================================================================

  int kmer_size = 31;
  int number_of_bits = 8;
  int bucket_size = 10;
  int max_retries = 10;

  dBGraph* hash_table = hash_table_new(number_of_bits, bucket_size,
                                       max_retries, kmer_size);

  if (hash_table==NULL)
    {
      die("unable to alloc the hash table. dead before we even started. OOM");
    }

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/test_pop_load_and_print/two_people_sharing_alu/just_one_of_the_two_people.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  //   >7SLRNA#SINE/Alu  plus GTTCAGAG at start and GTCAGCGTAG at end
  //   GTTCAGAGGCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG
  //   GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCC
  //   GATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG
  //   CGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAA
  //   ACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGT
  //   GAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTT
  //   AAAAAAAAAAAAAAAAAAAAGTCAGCCGTAG
  // plus loads of N's on end
   
  
  CU_ASSERT(seq_read==5158);//length of input sequence
  CU_ASSERT(seq_loaded==339);//amount loaded (ie removing Ns in this case)
  
  FILE* chrom_fptr = fopen("../data/test/pop_graph/test_pop_load_and_print/two_people_sharing_alu/person1.fa", "r");
  if (chrom_fptr==NULL)
    {
      die("Cannot open ../data/test/pop_graph/test_pop_load_and_print/two_people_sharing_alu/person1.fa");
    }

  int min_fiveprime_flank_anchor = 10;
  int min_threeprime_flank_anchor= 10;
  int max_anchor_span = 370;
  int length_of_arrays=740;
  int min_covg =1;
  int max_covg = 10;
  int max_expected_size_of_supernode=370;
  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, 0, 
							1,
							min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
							max_expected_size_of_supernode, length_of_arrays, hash_table, NULL,
							0, NULL, NULL, NULL, NULL, NULL, &make_reference_path_based_sv_calls_condition_always_true, &action_set_flanks_and_branches_to_be_ignored,
							&print_no_extra_info, NULL, NoIdeaWhatCleaning);
  
  CU_ASSERT(ret==0);

  hash_table_free(&hash_table);
  fclose(chrom_fptr);
}


void test_db_graph_make_reference_path_based_sv_calls_null_test_3()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // ===========================================================================
  // 3. Reference = Alu-NNNNN- same Alu, and person is identical to reference.
  //    This should find nothing. Note this really is testing something
  //    different. Since the reference is basically the same sequence repeated
  //    twice (with N's in between), there is a risk that we would match the
  //    5' anchor of the (one and only) supernode at the start of the reference,
  //    and the 3' anchor at the end of the 2nd copy of the repeat.
  //    To be clear, look at this:
  //        Ref:  ZamZamZamNNNNNNNNNNNNNNNNZamZamZam
  //        Indiv ZamZamZamNNNNNNNNNNNNNNNNZamZamZam
  //      So we look at supernode ZamZamzZam. We could, if we implemented things
  //      wrong, do this:
  //
  //        Ref: ZamZamZamNNNNNNNNNNNNNNNNZamZamZam
  //             Zam............................Zam
  //              ^                              ^
  //             5' anchor                    3' anchor    match start of
  //     supernode ZamZamZam with start of reference, and end of supernode with
  //     end of reference this would be a  bug - failing to realise that the
  //     entire supernode matches exactly the reference on the first instance of
  //     the repeat ZamZamZam
  //
  //    In fact, I found exactly this bug through this test.
  // ===========================================================================


  int kmer_size = 31;
  int number_of_bits = 8;
  int bucket_size = 10;
  int max_retries = 10;

  dBGraph* hash_table = hash_table_new(number_of_bits, bucket_size,
                                       max_retries, kmer_size);
  
  if (hash_table==NULL)
    {
      die("unable to alloc the hash table. dead before we even started. OOM");
    }

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/variations/one_person_is_alu_Ns_then_same_alu.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 0, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);
  
  //   >7SLRNA#SINE/Alu  
  // GTTCAGAGGCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG
  // GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCC
  // GATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG
  // CGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAA
  // ACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGT
  // GAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTT
  // AAAAAAAAAAAAAAAAAAAAGTCAGCCGTAGNNNNNNNNNNNNNNNNNNN
  // GTTCAGAGGCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG
  // GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCC
  // GATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG
  // CGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAA
  // ACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGT
  // GAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTT
  // AAAAAAAAAAAAAAAAAAAAGTCAGCCGTAG


  CU_ASSERT(seq_read==697);
  CU_ASSERT(seq_loaded==678);//not a typo - removing Ns

  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/one_person_aluNsalu.fa", "r");
  if (chrom_fptr==NULL)
    {
      die("Cannot open ../data/test/pop_graph/variations/one_person_aluNsalu.fa");
    }

  int min_fiveprime_flank_anchor = 10;
  int min_threeprime_flank_anchor= 10;
  int max_anchor_span = 370;
  int length_of_arrays= 740;
  int min_covg =1;
  int max_covg = 10;
  int max_expected_size_of_supernode=370;
  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, 0, 
							1,
							min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
							max_expected_size_of_supernode, length_of_arrays, hash_table, NULL,
							0, NULL, NULL, NULL, NULL, NULL, &make_reference_path_based_sv_calls_condition_always_true, &action_set_flanks_and_branches_to_be_ignored,
							&print_no_extra_info, NULL, NoIdeaWhatCleaning);

  
  CU_ASSERT(ret==0);

  hash_table_free(&hash_table);
  fclose(chrom_fptr);
}


void test_db_graph_make_reference_path_based_sv_calls_null_test_4()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // ===========================================================================
  // 4. Reference = 10kb of chromosome 1. Individual is that same sequence. This
  //    should find nothing.
  // ===========================================================================

  int kmer_size = 31;
  int number_of_bits = 13;
  int bucket_size = 10;
  int max_retries = 10;

  dBGraph* hash_table = hash_table_new(number_of_bits, bucket_size,
                                       max_retries, kmer_size);
  
  if (hash_table==NULL)
    {
      die("unable to alloc the hash table. dead before we even started. OOM");
    }

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/variations/one_person_is_10kb_of_chrom1.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 0, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);
  
  CU_ASSERT(seq_read==16320);
  
  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/person_1kb_chrom1.fa", "r");
  if (chrom_fptr==NULL)
    {
      die("Cannot open ../data/test/pop_graph/variations/person_1kb_chrom1.fa");
    }

  int min_fiveprime_flank_anchor = 20;
  int min_threeprime_flank_anchor= 10;
  int max_anchor_span = 10000;
  int length_of_arrays=20000;
  int min_covg =1;
  int max_covg = 10;
  int max_expected_size_of_supernode=10000;
  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, 0, 
							1,
							min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
							max_expected_size_of_supernode, length_of_arrays, hash_table, NULL,
							0, NULL, NULL, NULL, NULL, NULL, &make_reference_path_based_sv_calls_condition_always_true, &action_set_flanks_and_branches_to_be_ignored,
							&print_no_extra_info, NULL, NoIdeaWhatCleaning);

  CU_ASSERT(ret==0);

  hash_table_free(&hash_table);
  fclose(chrom_fptr);

}


void test_db_graph_make_reference_path_based_sv_calls_null_test_5()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // ===========================================================================
  // Reference is two copies of a single sequence, tandem repeat of about 36
  // bases. Individual is the same, with an Alu inserted between. Since the
  // supernode in the individual has one copy of th repeat, and then the Alu,
  // and then stops, you should be unable to find an anchor at the Alu-end of
  // the supernode. So this should find nothing.
  // ===========================================================================


  int kmer_size = 31;
  int number_of_bits = 15;
  int bucket_size = 10;
  int max_retries = 10;

  dBGraph* hash_table = hash_table_new(number_of_bits, bucket_size,
                                       max_retries, kmer_size);

  if (hash_table==NULL)
    {
      die("unable to alloc the hash table. dead before we even started. OOM");
    }

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/variations/two_people_one_with_alu_insertion.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/person_without_alu.fa", "r");
  if (chrom_fptr==NULL)
    {
      die("Cannot open ../data/test/pop_graph/variations/person_without_alu.fa");
    }


  int min_fiveprime_flank_anchor = 3;
  int min_threeprime_flank_anchor= 3;
  int max_anchor_span =400;
  int length_of_arrays=800;
  int min_covg =1;
  int max_covg = 100000;
  int max_expected_size_of_supernode=400;



  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr,  1, 
							0,
							min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
							max_expected_size_of_supernode, length_of_arrays, hash_table, NULL,
							0, NULL, NULL, NULL, NULL, NULL, &make_reference_path_based_sv_calls_condition_always_true, &action_set_flanks_and_branches_to_be_ignored,
							&print_no_extra_info, NULL, NoIdeaWhatCleaning);
  
  CU_ASSERT(ret==0);

  hash_table_free(&hash_table);
  fclose(chrom_fptr);

}

//numbering of tests goes back to 1 - these unti tests are for cases where we DO expect to find a variant
void test_db_graph_make_reference_path_based_sv_calls_test_1()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // ===========================================================================
  // 1. Reference = short sequence
  //    Individual is idential except for one base change
  // ===========================================================================



  int kmer_size = 7;
  int number_of_bits = 8;
  int bucket_size = 10;
  int max_retries = 10;

  dBGraph* hash_table = hash_table_new(number_of_bits, bucket_size,
                                       max_retries, kmer_size);
  
  if (hash_table==NULL)
    {
      die("unable to alloc the hash table. dead before we even started. OOM");
    }

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/variations/two_people_short_seq_with_one_base_difference.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT(seq_read==192);
  
  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/second_person_same_short_seq_one_base_diff.fa", "r");
  if (chrom_fptr==NULL)
    {
      die("Cannot open ../data/test/pop_graph/variations/second_person_same_short_seq_one_base_diff.fa");
    }

  int min_fiveprime_flank_anchor = 5;
  int min_threeprime_flank_anchor= 5;
  int max_anchor_span =40;
  int length_of_arrays=80;
  int min_covg =1;
  int max_covg = 10;
  int max_expected_size_of_supernode=40;


  //this time, we expect results. So  malloc some arrays to hold the results. Only use this option for testing, as in general use, we expect to find millions
  // of variants

  char** return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  char** return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  char** return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  char** return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      die("Failed to alloc return_something_array - cannot start test\n");
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 30);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 30);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 30);
  return_branch2_array[0] = (char*) malloc( sizeof(char) * 30);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      die("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';
  
  int return_variant_start_coords_array[2];
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  int* return_variant_start_coords_array_ptr[2];
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);


  //now do the test!!!
  FILE* fp = fopen("../data/tempfiles_can_be_deleted/temp_outputfile_trustedpath_sv_caller_test1", "w");
  //indiv in colour 0, ref in colour 1
  int ret = 
    db_graph_make_reference_path_based_sv_calls_in_subgraph_defined_by_func_of_colours(chrom_fptr, &element_get_colour0, &element_get_covg_colour0, 1,
										       min_fiveprime_flank_anchor, min_threeprime_flank_anchor, 
										       max_anchor_span, min_covg, max_covg, 
										       max_expected_size_of_supernode, length_of_arrays, hash_table, fp,
										       1, return_flank5p_array, return_trusted_branch_array, 
										       return_branch2_array, return_flank3p_array, 
										       return_variant_start_coords_array_ptr, 
										       &make_reference_path_based_sv_calls_condition_always_true_for_func_of_colours, 
										       &action_set_flanks_and_branches_to_be_ignored,
										       &print_no_extra_info, NULL, NoIdeaWhatCleaning, 0);//last 0 - start numbering variant at 0
  
  fclose(fp);

  CU_ASSERT(ret==1);
  CU_ASSERT_STRING_EQUAL("AATAGACGCCCACACCTGATAG", return_flank5p_array[0]);
  CU_ASSERT_STRING_EQUAL("TAGCCACA", return_trusted_branch_array[0]);
  CU_ASSERT_STRING_EQUAL("AAGCCACA", return_branch2_array[0]);
  CU_ASSERT_STRING_EQUAL("CTGTACTTGTA", return_flank3p_array[0]);

  CU_ASSERT(return_variant_start_coords_array[0]==23);

  //cleanup
  hash_table_free(&hash_table);
  fclose(chrom_fptr);
  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);

}



void test_db_graph_make_reference_path_based_sv_calls_test_2()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // ===========================================================================
  // 2. Reference = Alu-NNNNN- same Alu, and person is identical to reference
  //    except for a single base difference. This should find the single variant!
  //    It is also an example where we traverse the supernode in the reverse
  //    direction to that in which is appears in our array. This is therefore a
  //    good unit test to ensure that our code works in this case.
  //    Note that the coverage is double on a lot of these nodes because we have
  //    two copies of the Alu. But that doesn't in any way add credence to the
  //    variant call - it's paralogous sequence that doubles the coverage.
  // ===========================================================================

  int kmer_size = 31;
  int number_of_bits = 8;
  int bucket_size = 10;
  int max_retries = 10;

  dBGraph* hash_table = hash_table_new(number_of_bits, bucket_size,
                                       max_retries, kmer_size);

  if (hash_table==NULL)
    {
      die("unable to alloc the hash table. dead before we even started. OOM");
    }

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/variations/two_people_both_alu_Ns_alu_with_one_base_difference.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/one_person_aluNsalu_PLUS_SINGLE_BASE_CHANGE.fa", "r");
  if (chrom_fptr==NULL)
    {
      die("Cannot open ../data/test/pop_graph/variations/one_person_aluNsalu_PLUS_SINGLE_BASE_CHANGE.fa");
    }

  int min_fiveprime_flank_anchor = 10;
  int min_threeprime_flank_anchor= 10;
  int max_anchor_span =350;
  int length_of_arrays=700;
  int min_covg =1;
  int max_covg = 10;
  int max_expected_size_of_supernode=350;


  char** return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  char** return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  char** return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  char** return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      die("Failed to alloc return_something_array - cannot start test\n");
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 400);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 400);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 400);
  return_branch2_array[0] = (char*) malloc( sizeof(char) * 400);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      die("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';


  int return_variant_start_coords_array[2];  
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  int* return_variant_start_coords_array_ptr[2];
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);


  FILE* fp = fopen("../data/tempfiles_can_be_deleted/temp_outputfile_trustedpath_sv_caller_test2", "w");

  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, 0, 
							1,
							min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
							max_expected_size_of_supernode, length_of_arrays, hash_table, fp,
							1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array,return_variant_start_coords_array_ptr,
							&make_reference_path_based_sv_calls_condition_always_true, &action_set_flanks_and_branches_to_be_ignored,
							&print_no_extra_info, NULL, NoIdeaWhatCleaning);
  fclose(fp);

  CU_ASSERT(ret==1);
  CU_ASSERT_STRING_EQUAL("GTTCAGAGGCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG", return_flank5p_array[0]);
  CU_ASSERT_STRING_EQUAL("GCTGTAGTGCGCTATGCCGATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAGCGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAAACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGTGAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTTAAAAAAAAAAAAAAAAAAAAGTCAGCCGTAG", return_flank3p_array[0]);
  CU_ASSERT_STRING_EQUAL("CTGGGAGGATCGCTTGAGTCCAGGAGTTCTGG", return_trusted_branch_array[0]);
  CU_ASSERT_STRING_EQUAL("GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGG", return_branch2_array[0]);
  CU_ASSERT(return_variant_start_coords_array[0]==59);

  hash_table_free(&hash_table);
  fclose(chrom_fptr);

  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);

}


void test_db_graph_make_reference_path_based_sv_calls_test_3()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // ===========================================================================
  // 3. Reference is just a short sequence. Individual is missing 2 bases in the
  //    middle--> we should find this. Note this test is also interesting
  //    because the reference and individual both contain a closed loop in their
  //    de Bruijn graph
  // ===========================================================================


  int kmer_size = 7;
  int number_of_bits = 8;
  int bucket_size = 10;
  int max_retries = 10;

  dBGraph* hash_table = hash_table_new(number_of_bits, bucket_size,
                                       max_retries, kmer_size);

  if (hash_table==NULL)
    {
      die("unable to alloc the hash table. dead before we even started. OOM");
    }

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/variations/two_people_one_with_2_bases_missing.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/person_without_2_bases_missing.fa", "r");
  if (chrom_fptr==NULL)
    {
      die("Cannot open ../data/test/pop_graph/variations/person_without_2_bases_missing.fa");
    }


  int min_fiveprime_flank_anchor = 2;
  int min_threeprime_flank_anchor= 2;
  int max_anchor_span =40;
  int length_of_arrays=80;
  int min_covg =1;
  int max_covg = 1000;
  int max_expected_size_of_supernode=40;


  char** return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  char** return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  char** return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  char** return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      die("Failed to alloc return_something_array - cannot start test\n");
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 80);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 80);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 80);
  return_branch2_array[0] = (char*) malloc( sizeof(char) * 80);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      die("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';


  int return_variant_start_coords_array[2];    
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  int* return_variant_start_coords_array_ptr[2];
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);


  FILE* fp = fopen("../data/tempfiles_can_be_deleted/temp_outputfile_trustedpath_sv_caller_test3", "w");


  //individual has 2 missing bases, reference does not
  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, 0, 
							1,
							min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
							max_expected_size_of_supernode, length_of_arrays, hash_table, fp,
							1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, return_variant_start_coords_array_ptr,
							&make_reference_path_based_sv_calls_condition_always_true, &action_set_flanks_and_branches_to_be_ignored,
							&print_no_extra_info, NULL, NoIdeaWhatCleaning);
  
  fclose(fp);

  CU_ASSERT(ret==1);
  CU_ASSERT_STRING_EQUAL("CCACACCTGA", return_flank5p_array[0]);
  CU_ASSERT_STRING_EQUAL("ACAC", return_flank3p_array[0]);
  CU_ASSERT_STRING_EQUAL("TAGACCCC", return_trusted_branch_array[0]);
  CU_ASSERT_STRING_EQUAL("GACCCC", return_branch2_array[0]);
  CU_ASSERT(return_variant_start_coords_array[0]==20);


  hash_table_free(&hash_table);
  fclose(chrom_fptr);

  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);

}


void test_db_graph_make_reference_path_based_sv_calls_test_4()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // ===========================================================================
  // Reference is just a short sequence. Individual has an extra 2 bases in the
  // middle--> we should find this
  // ===========================================================================

  int kmer_size = 7;
  int number_of_bits = 8;
  int bucket_size    = 10;
  int max_retries=10;

  dBGraph* hash_table = hash_table_new(number_of_bits, bucket_size,
                                       max_retries, kmer_size);

  if (hash_table==NULL)
    {
      die("unable to alloc the hash table. dead before we even started. OOM");
    }

  // We use the same two people as last time, but swap their roles

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/variations/two_people_one_with_2_bases_missing.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/person_with_2_bases_missing.fa", "r");
  if (chrom_fptr==NULL)
    {
      die("Cannot open ../data/test/pop_graph/variations/person_with_2_bases_missing.fa");
    }


  int min_fiveprime_flank_anchor = 2;
  int min_threeprime_flank_anchor= 2;
  int max_anchor_span =40;
  int length_of_arrays=80;
  int min_covg =1;
  int max_covg = 1000;
  int max_expected_size_of_supernode=40;


  char** return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  char** return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  char** return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  char** return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      die("Failed to alloc return_something_array - cannot start test\n");
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 80);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 80);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 80);
  return_branch2_array[0] = (char*) malloc( sizeof(char) * 80);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      die("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';


  int return_variant_start_coords_array[2];      
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  int* return_variant_start_coords_array_ptr[2];
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);

  FILE* fp = fopen("../data/tempfiles_can_be_deleted/temp_outputfile_trustedpath_sv_caller_test4", "w");



  //cf previous test - we are using individua with index 1, not 0 as last time. ie we have swapped which is individual and which is reference
  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr,  1, 
							1,
							min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
							max_expected_size_of_supernode, length_of_arrays, hash_table, fp,
							1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, return_variant_start_coords_array_ptr,
							&make_reference_path_based_sv_calls_condition_always_true, &action_set_flanks_and_branches_to_be_ignored,
							&print_no_extra_info, NULL, NoIdeaWhatCleaning);

  fclose(fp);

  CU_ASSERT(ret==1);
  CU_ASSERT_STRING_EQUAL("CCACACCTGA", return_flank5p_array[0]);
  CU_ASSERT_STRING_EQUAL("ACAC", return_flank3p_array[0]);
  CU_ASSERT_STRING_EQUAL("TAGACCCC", return_branch2_array[0]);
  CU_ASSERT_STRING_EQUAL("GACCCC", return_trusted_branch_array[0]);
  CU_ASSERT(return_variant_start_coords_array[0]==20);


  hash_table_free(&hash_table);
  fclose(chrom_fptr);

  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);
}


void test_db_graph_make_reference_path_based_sv_calls_test_5()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // ===========================================================================
  // 5. Reference is a single sequence which is a single supernode. Individual
  //    is the same, with an Alu inserted in the middle. Should be able to find
  //    anchors at both ends this time, and find the insertion.
  //    Note that by missing the alu in the reference, it has a load of nodes
  //    that won't be seen in the individual, so for example here we will ask
  //    for a min 3' flank of 8. This actually means 31+8 bases of agreement.
  // ===========================================================================


  int kmer_size = 31;
  int number_of_bits = 15;
  int bucket_size = 10;
  int max_retries = 10;

  dBGraph* hash_table = hash_table_new(number_of_bits, bucket_size,
                                       max_retries, kmer_size);

  if (hash_table==NULL)
    {
      die("unable to alloc the hash table. dead before we even started. OOM");
    }

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/variations/two_people_one_with_alu_inserted_mid_supernode.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/person_with_one_supernode_and_without_alu.fa", "r");
  if (chrom_fptr==NULL)
    {
      die("Cannot open ../data/test/pop_graph/variations/person_with_one_supernode_and_without_alu.fa");
    }


  int   min_fiveprime_flank_anchor = 3;
  int   min_threeprime_flank_anchor=  8;
  int   max_anchor_span =500;
  int   length_of_arrays=1000;
  int   min_covg =1;
  int   max_covg = 100000;
  int   max_expected_size_of_supernode=500;


  char** return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  char** return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  char** return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  char** return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      die("Failed to alloc return_something_array - cannot start test\n");
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 500);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 500);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 500);
  return_branch2_array[0] = (char*) malloc( sizeof(char) *500);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      die("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';


  int return_variant_start_coords_array[2];        
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  int* return_variant_start_coords_array_ptr[2];  
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);


  
  FILE* fp = fopen("../data/tempfiles_can_be_deleted/temp_outputfile_trustedpath_sv_caller_test5", "w");
  //indiv colour 0, ref col1
  int ret = db_graph_make_reference_path_based_sv_calls_in_subgraph_defined_by_func_of_colours(chrom_fptr, &element_get_colour0, &element_get_covg_colour0, 1,
											       min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
											       max_expected_size_of_supernode, length_of_arrays, hash_table, fp,
											       1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, 
											       return_variant_start_coords_array_ptr, &make_reference_path_based_sv_calls_condition_always_true_for_func_of_colours, 
											       &action_set_flanks_and_branches_to_be_ignored,
											       &print_no_extra_info, NULL, NoIdeaWhatCleaning,0); //last 0 - start numbering variant at 0
  fclose(fp);

  CU_ASSERT(ret==1);

  //first base of Alu insert = first base of displaced sequence. ie we had xxxxGyyyy and then we insert Alu: xxxx[G-Alu]Gyyyy. So first base of Alu look slike is in flank5p
  CU_ASSERT_STRING_EQUAL("TGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATGGTG", return_flank5p_array[0]);
  CU_ASSERT_STRING_EQUAL("TAGCCAGG", return_flank3p_array[0]);
  CU_ASSERT_STRING_EQUAL("AAACTCCGTCTGTACTAATAGTACAAAAAT", return_trusted_branch_array[0]);
  CU_ASSERT_STRING_EQUAL("CCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCCGATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAGCGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAAACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGTGAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTTAAAAAAAAAAAAAAAAAAAAGAAACTCCGTCTGTACTAATAGTACAAAAAT", return_branch2_array[0]);
  CU_ASSERT(return_variant_start_coords_array[0]==40); //note the insertion happens at coordinate 39, but the first inserted base is the same as what would be there anyway - a G. 
  
  hash_table_free(&hash_table);
  fclose(chrom_fptr);
  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);

}



void test_db_graph_make_reference_path_based_sv_calls_test_6()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // ===========================================================================
  // 6. Reverse the roles of previous test. Reference has an Alu, individual
  //    does not. Should be able to find anchors at both ends this time, and
  //    find the deletion
  // ===========================================================================


  int kmer_size = 31;
  int number_of_bits = 15;
  int bucket_size = 10;
  int max_retries = 10;

  dBGraph* hash_table = hash_table_new(number_of_bits, bucket_size,
                                       max_retries, kmer_size);

  if (hash_table==NULL)
    {
      die("unable to alloc the hash table. dead before we even started. OOM");
    }

  // We use the same two people as last time, but swap their roles

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/variations/two_people_one_with_alu_inserted_mid_supernode.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);


  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/person_with_alu_in_middle_of_supernode.fa", "r");
  if (chrom_fptr==NULL)
    {
      die("Cannot open ../data/test/pop_graph/variations/person_with_alu_in_middle_of_supernode.fa");
    }


  int  min_fiveprime_flank_anchor = 3;
  int  min_threeprime_flank_anchor= 8;
  int  max_anchor_span =500;
  int  length_of_arrays=1000;
  int  min_covg =1;
  int  max_covg = 100000;
  int  max_expected_size_of_supernode=500;


  char**  return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  char**  return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  char**  return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  char**  return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      die("Failed to alloc return_something_array - cannot start test\n");
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 500);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 500);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 500);
  return_branch2_array[0] = (char*) malloc( sizeof(char) *500);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      die("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';

  int return_variant_start_coords_array[2];          
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  int* return_variant_start_coords_array_ptr[2];    
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);
  
  FILE* fp = fopen("../data/tempfiles_can_be_deleted/temp_outputfile_trustedpath_sv_caller_test6", "w");

  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr,  1, 
							0,
							min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
							max_expected_size_of_supernode, length_of_arrays, hash_table, fp,
							1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, 
							return_variant_start_coords_array_ptr, &make_reference_path_based_sv_calls_condition_always_true, 
							&action_set_flanks_and_branches_to_be_ignored,
							&print_no_extra_info, NULL, NoIdeaWhatCleaning);
  fclose(fp);


  CU_ASSERT(ret==1);

  //first base of Alu insert = first base of displaced sequence. ie we had xxxxGyyyy and then we insert Alu: xxxx[G-Alu]Gyyyy. So first base of Alu look slike is in flank5p
  CU_ASSERT_STRING_EQUAL("TGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATGGTG", return_flank5p_array[0]);
  CU_ASSERT_STRING_EQUAL("TAGCCAGG", return_flank3p_array[0]);
  CU_ASSERT_STRING_EQUAL("AAACTCCGTCTGTACTAATAGTACAAAAAT", return_branch2_array[0]);
  CU_ASSERT_STRING_EQUAL("CCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCCGATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAGCGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAAACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGTGAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTTAAAAAAAAAAAAAAAAAAAAGAAACTCCGTCTGTACTAATAGTACAAAAAT", return_trusted_branch_array[0]);
  CU_ASSERT(return_variant_start_coords_array[0]==40); //note the insertion happens at coordinate 39, but the first inserted base is the same as what would be there anyway - a G. 


  hash_table_free(&hash_table);
  fclose(chrom_fptr);
  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);
}


void test_db_graph_make_reference_path_based_sv_calls_test_7()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // ===========================================================================
  // 7. Reference is an Alu with a different Alu inserted in the middle.
  //    Individual is just the first Alu - ie they have a deletion of an Alu
  //    from within an Alu.
  // ===========================================================================

  int  kmer_size = 31;
  int  number_of_bits = 15;
  int  bucket_size = 10;
  int  max_retries = 10;

  dBGraph* hash_table = hash_table_new(number_of_bits, bucket_size,
                                       max_retries, kmer_size);
  
  if (hash_table==NULL)
    {
      die("unable to alloc the hash table. dead before we even started. OOM");
    }

  // We use the same two people as last time, but swap their roles

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/variations/two_people_one_is_alu_other_has_2nd_alu_inserted.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/person_with_alu_in_middle_of_alu.fa", "r");
  if (chrom_fptr==NULL)
    {
      die("Cannot open ../data/test/pop_graph/variations/person_with_alu_in_middle_of_alu.fa");
    }


  
 //    reference is
 //   > AluJo#SINE/Alu inserted in middle of 7SLRNA#SINE/Alu
 //   GCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG
 //   GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCC
 //   GATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG
 //   GCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGA
 //   GGCGGGAGGATCGCTTGAGCCCAGGAGTTCGAGACCAGCCTGGGCAACAT
 //   AGCGAGACCCCGTCTCTACAAAAAATACAAAAATTAGCCGGGCGTGGTGG
 //   CGCGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGGATCGCT
 //   TGAGCCCAGGAGTTCGAGGCTGCAGTGAGCTATGATCGCGCCACTGCACT
 //   CCAGCCTGGGCGACAGAGCGAGACCCTGTCTCAAAAAAAAAAAAAAAAAA
 //   AAAAAAAAAAAA
 //   CGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAA
 //   ACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGT
 //   GAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTT
 //   AAAAAAAAAAAAAAAAAAAA
 //   ..plus N's on the end
    
    
 //   and individual is
    
 //   >7SLRNA#SINE/Alu
 //   GCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG
 //   GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCC
 //   GATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG
 //   CGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAA
 //   ACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGT
 //   GAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTT
 //   AAAAAAAAAAAAAAAAAAAA
    
    


  int  min_fiveprime_flank_anchor = 3;
  int  min_threeprime_flank_anchor= 3;
  int  max_anchor_span =500;
  int  length_of_arrays=1000;
  int  min_covg =1;
  int  max_covg = 100000;
  int  max_expected_size_of_supernode=500;


  char** return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  char** return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  char** return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  char** return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      die("Failed to alloc return_something_array - cannot start test\n");
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 500);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 500);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 500);
  return_branch2_array[0] = (char*) malloc( sizeof(char) *500);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      die("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';


  int return_variant_start_coords_array[2];          
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  int* return_variant_start_coords_array_ptr[2];    
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);

  FILE* fp = fopen("../data/tempfiles_can_be_deleted/temp_outputfile_trustedpath_sv_caller_test7", "w");



  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr,  1, 
							0,
							min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
							max_expected_size_of_supernode, length_of_arrays, hash_table, fp,
							1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, 
							return_variant_start_coords_array_ptr, &make_reference_path_based_sv_calls_condition_always_true, 
							&action_set_flanks_and_branches_to_be_ignored,
							&print_no_extra_info, NULL, NoIdeaWhatCleaning);


  fclose(fp);


  CU_ASSERT(ret==1);

  CU_ASSERT_STRING_EQUAL("GCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCCGATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG", return_flank5p_array[0]);
  CU_ASSERT_STRING_EQUAL("GAACCGGCCCAGGTCGGAAACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGTGAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTTAAAAAAAAAAAAAAAAAAAA", return_flank3p_array[0]);
  CU_ASSERT_STRING_EQUAL("CGGGGGACCACCAGGTTGCCTAAGGAGGGGT", return_branch2_array[0]);
  CU_ASSERT_STRING_EQUAL("GCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGAGGATCGCTTGAGCCCAGGAGTTCGAGACCAGCCTGGGCAACATAGCGAGACCCCGTCTCTACAAAAAATACAAAAATTAGCCGGGCGTGGTGGCGCGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGGATCGCTTGAGCCCAGGAGTTCGAGGCTGCAGTGAGCTATGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACCCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGGGGGACCACCAGGTTGCCTAAGGAGGGGT", return_trusted_branch_array[0]);
  CU_ASSERT(return_variant_start_coords_array[0]==151); 

  hash_table_free(&hash_table);
  fclose(chrom_fptr);
  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);
}


void test_db_graph_make_reference_path_based_sv_calls_test_8()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // ===========================================================================
  // 8. Reference is 10 kb of chromosome 1 plus 1kb of sequence inserted
  //    mid-supernode, and person is identical, except for that 1kb is missing 
  // ===========================================================================

  int  kmer_size = 31;
  int  number_of_bits = 15;
  int  bucket_size = 10;
  int  max_retries = 10;

  dBGraph* hash_table = hash_table_new(number_of_bits, bucket_size,
                                       max_retries, kmer_size);

  if (hash_table==NULL)
    {
      die("unable to alloc the hash table. dead before we even started. OOM");
    }

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/variations/two_people_one_with_1kb_deletion.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);


  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/first_person_10kb_chrom1_plus_1kb_inserted_mid_supernode.fa", "r");
  if (chrom_fptr==NULL)
    {
      die("Cannot open ../data/test/pop_graph/variations/first_person_10kb_chrom1_plus_1kb_inserted_mid_supernode.fa");
    }

  // Let's be clear about what this test looks like. 

  //       The individual is 10kb of chromosome 1, which has the following supernodes

//>node_0  cylcle:false length:4546 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:2 fst:2 lst_coverage:2 lstr:1 lstf:2
//TGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATGGTGAAACTCCGTCTGTACTAATAGTACAAAAATTAGCCAGGCGTGATGGTGTGCACCTGTAGTCCTTGCTACTCAGAAAGCTGAGGCAGGAGAATCGCTTGTACCCAGGAGGCAGAGGTTGCAGTGAGCAGAGATTGTGCCACTGCACTCCATCCTGGGTGACAGAGTGCTATGAGTCACCACACCTGGTATGAGCCACCGTGCCTGGCCCACAATGACTTTTACACATGTTGTTAAATCATCTTACAGATTTTATAATTTGGGGGAAGAAAAGTTTTACTAAATTGTCTTTTAATGGAAACTCTACAAGAACCAGAATCTTTGCTTTGTTCACTTATGTATCCATTCCTAGGCCTAGAAAAATGTCTGACACATAGCGGCAATTATTCATTGAATAAATGGACCCAGGGATAGTACATTAGCTATGCTATATGCATACATTAAAGATGTAGATTATCGACTTTCAAAAGATAATTAATGTAACTTCTTACTGCTTCTGAACATGTTTGTGAGTTATATTGCTGAGGGACCTTTATCTTCTCATTCTTTCATCTTAACCCAGTGTTATAAAATTGAAATCACCAATATTATTCCATATCTAAAATTAATATCTACCTTGTAAAAAATATCACTCTGCTGCATTTGAGAATAGACTTTTTAGGTAATAATGATGCAATCCATAGGGTTTTTTGGGGGCACAGAGGGATTCATGCTAACAGAACATTTTATTTTCTATTTTCCCAGAGCTGTAAAACATGAAATTACGGTAGTATAAGGCATATTTTTACTCTTTTTATAATTTTTTCTAAAAAAAATTAGTGTTTGTTCCCTATATAACTTTTAACTTTATAGGTAAATATTTGTTTCTTTCAGCTCCAGTTTTATGTGAAATAGAGTTTTCAGATTTATGTAGCATGGAAAGTTTTAATACGTCAGAGTTACTGATTTTTGCCAATCATTTTCTCAATTATTTCTTTTTTATCTTTAGTTGATTTTTTTGTAGTGACACATTTTGTTTCTAGTCTCATTTCCTTTTGTTTATATTCTATGTATATTTCATTTTTGGTTACTATGAGAATTACATATAACATCCTAGAGTTATAACATTTTAATTTGAATTTATTTCAACTTAAGTTCAATCACATACCAAAATTCTACTGCTATATATATAGCTCTACTCTTTTTATGTTATTGATGTGACAAATTATATCTTTATTCATTGTATACCAGCTAACAGATTTACAATTACATTTTATGCATTTGCCTTTTAAATTATGTAGAAAATAAAAAGCAGAGTTACAAACCAAAATTACAATAGGACTGTTTTTATGTTTGTTTATGTATTTACCTTTACCAGAGAGCTTTGTATATTCATACAGCTTGCTTATTTACTTATATAGTTATTGCCTAGAGTTCATTTATTTCAACCTGAAGGACTTAACACTTCCTGAATGTCAAATTCAGGGATAAATGGATTTTTTTCAGTTTTAAAAAAAAATCCGGAAATGTCTTAATTTCTCCTTCATTTTTGAAGGATAAGTTTTCCAGCTATATATTTCTCAATTGACAGGTTTCTTCATTATTTTAAATATATAATCCACTGCCTACTGGCCTTCAAGGTTTCTGCTGAGAAATCAGCTGCTAATGTTATCTGGATCCCTATCTGTGAGAGTTGCTCTTCTCTCTGAGTTTTCAACATTCTCCCATTATCTTTTTTTTGTTTGTTTTTGAGACAAATAATTGTACATATTCATGGGATACAGAGTGATATTTTGATACATGTATACAATGCCCAGTGATCAAATAAGGATAATTAGCGTATCCATCACCTCAAATAGTTGTCATTTATTTGTATTGTGAACAGTCAACATTCTTTCTTCTAGTTTTTTAAATTTATAAACATTTAAATTTTATTACAGAAATTTAAATTTTTTGATTCTGAAAAAGTCATATATGTATGCAACATTTTTTATCGTTTATTTATATATTTATGCATCTTTCCTTTTAGTTTTGACAGAGATTTTCTATTTTATCATTATTTCAAAAGAACTCTTACCTGTATTTATTTATCAAGTATATTTCCCTTGTTTTTTCCTAGTATATTAATTTATTTACTTATCTTCTAAAAATCCTCCATATAATCTGTTTATTTTGTTTCCTTTCTATAATTTCTTCAATAATTAGTTCTGTTCTATTTTCCATTAAAATATTTAAATCTTGTATGAATTTTTGTCAGATTAGAAATTTAGGGCGTTTCTTAATTTCTCTATACTCTAGCTTTTGACTTTTTTTTTCTGACCTAAGAGGTATTTAGAGCACATTTTAGATTTTTTATTTTGACTAATCATTTAAAATGTATACTAATCTTCAATTTAAATAAAAAACTGGTCTATAGTGACAAAAATTACAAATGAGCCTAACTAATAAATTATCAGCTGTGTTTATATGTATAAGCATGCACAGATTTTGGTAAATATGTACATAGTATATTGGTGAGCTTATTTTTATCATTCTTAACTCATTGTGTAGTCTAAACGTTGGGGAAAAAATAACATACAATAATCAGATGGTGTGAATAAGAAAATTGTTCTAATGTTTGTAAACCAAGCAACTGTTTTAACTGCTCCCCTCTTCCTGATTGACTTCTAAAAGGGATTGATCCATATTGGGTCCTATCATGTACGTCACGGTATAACATCTCCAGCTATAAAATGGAAATTTGAGAATAACTTTGCTGCTACTCAGATATATTTTATTTCAAAAACATACACTAAGGTGTTGCTGTTGGATCTTTCCAAAAACATATTCACACAGAACTTTCAATCACACTGAGCCATATTTGAACAATCTTTCAAGGTCAGCTCTGGCATAAGCTAACATTATACCATTTAACTCAGAAATTTCTTTAGTATTTGATTAATGGGTTTATGTTTGATATGTAATGTAATTTTCTAATACTAAATCAAGTGGTAATTTTGTTAGTCAAGTTGATTTAGTGGCTTGGGAAGAAAGCTTTTAATGTTCCCCTAATTTTTCTTTCCTTTGACATGATCCTTCACATGTCTTATTTTGCTTAGTGATTTTTCTTTTTTTTTTTTTTTTTTTGAGACAGGGTCTTACTCTACCACCCAGGCTTGAGTGCAGTGGTGCAATCACAGCTCATTGCAGCCTTGACCTCCCAGACTCAAGCTATTCTTCCACCTCAGCCTCCCAAGTAGCTGGTACTACAGGCACATGCCACCAAACTTGGCTAATTTTTGTATTTTTTGTAGAGACAGAGTTTTGCCAAATTCTCAGGCTGGTCTGGAATTTCTGGGCTCAAGTAATCCTGCCTTGGCCTCCCAACATGCTGATATTACAGACATAAGCCACAGTACCTGGCCAGTTTTCTTTTTAAAAAAATCTATTGGTTATTAATTTGAAGCCTTCCTTTTCATAGCTGTGCTCCTTAATTGGGAGCAAACATGAATGGACCACAACTTAGCCAATTTTCTATATACGATCTTTGCCATCCTAATTTAAAGGAATATTAATTCTTTCTTTTCCTCTTTCATTCCACAAACCTGTATTGACTACATCTAAGTTCTAAATGGTGCACTGGATGTTGAAAAAGTTGATGATGAGCAAGAACAAAATTCCTGCTTTCAGGAGACTTACAGTTCAATATGGGAAATATAATTTGTTAAAATATAAAAGTGCAATTGTGTTACATGCTGTACGAAGTACATGTTGACATGTGAGCATATAATAAATGGGCTGGAGGCCAGAGGATTGCCAAAGAGAATGGGCCTCCTGCTGAGATGAAAAGTTGAGCAGGGATTAGTTGGCAAAAGTGGAGGGACGATCCTTTCTAGGCAGGAGGAAGAACATGTACAGAATCTCTGAGGTGTGATGCGACAAAGTCTATATAAAAAACTGAAGAAAGGTCTAATATGGCTTAAATACAGAAGCTAGTAGGAGAGGAGTCGAAAAGAGGCTGGAGAAGTAGAAAGTGTCTGCATTCTGCAGGAACTTATATTGTATAAAATATATTGTATATTGTATAAAAAGAATTTCTCTTTATTCTAAGTGCAATGTGAAGCCAATGAAGTGCTTTAAACAGGTGATGTGATTTGATTGAATTTATTACTTCACTTAACAAATATTCATTACATGCCCACTGTTTGTCAGATATTGCTCTAGCCCCTGGTGATACAGTAGGGAATAAAACAGGCAAAAATCCCTGTCCTCTTGCAGCTTATAATGGACTGCAATGTTTAATATGTCAGAGGAGGTCCACGGAGGAGTGACTTCTAAGCAAGAATCTGAAAAAAATGAGGATATCTAAGGAGGGAACAAATGGTTCAAAAGCCCTATAATTGCAAGCAGGCATGATGAAGCAATTGCAGTTGTCCTGACTCTCAACACCGTGGAACTCAAAGGAGATGGAAAGATTCCTTCTCTCCCTCATATATTTTCTCTCTTTCTGTCTATATATATAGAATATGAGACATTTCCCTAATCATTATGT
//>node_1  cylcle:false length:1696 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:2 fst:1 lst_coverage:2 lstr:1 lstf:2
//ACATAATGATTAGGGAAATGTCTCATATTCTTCTACTCAGAAATAAGCAATATAGCAATTACTGTTTTTTACATTTTACAGTTACAGTTTCAGAGAAAGTTTGATATTTATCTAAAATTTTTCAATGTATGAACTTTTTCATTTGACAAACCATAATTGTACATATTCTTGGGATACAGAGTGATATTTTCTTTACATGTATAGAATGTGTAGTGATCAAATCAGGGTAATTTCCACTAATTTAAAATGCCACCTTTATGTTATTGTAATTTATATATATACTATATATATACACACACACATATATATATATACATGTCCACATACAGTGTGTGTGTGCACATGTACACACATGCATATGTGTATATAATGCCCAGTATAAGCAATGTGCACAAATAAAATTAGCTAACAGAGATAGTATAGAGTGAGAGGAGAGGCAGATTAATCTTTGAGGAAAAGCACAATTTTATAGCTGAATGGAGAAAGCTGAGGTGGTTTCTAAGATGGAGAATAAGACGAAAAATGTAAGTACGTTGTCTGACTGAATTCAAGAAAGAAGGGTAAAAGAGAAGAAAGTAGTGGTCTTATCATTAAATGCCACAGAGAGGTAAAGATAAAGACAACATATTGTTTTGGGTTTAGTAATTTAAGGGTTACCAAATTCCGTTTTGGAGGAGGAACAGATTCCATGTCCACTAGAATGGAATGAACAAGAAATGGAGGAGGAAAATAGGTAGTTTTTCAAAAGTTTTCAAAAATATGAAAAGAAGAAATGAAGTGGTACTTGGAAGAGATTGTTGAAATGGGAGAGACTATGGTGGCTTGTTTAGAAGCAGTTGAGATAGATCCAATTGAGATAGAGATATTGACTATATAAACAAAAGAATGACAAATTAATAGTGTAATGGATAACTTGACTTTGGCAAATATTGTGAATTTTTGTGAAAGTACAACTAAAAGGCAATGTCACTCCAATAATCACCAGAGTAATCAATTTGCTTATTGCTGTCCCTTTAAATATAGTTCTCTGGTATCAACTAACATGTTTTTAACTAATGATGCTTCTTAAAGAAAAGGGAAAAGACCTTTTTCTTTCTTTCAGTCTTCAATGATTCACTGCTTCATCTCGCTCCACCAAAGATAAATGAAATCTACATCTCTTATACATTAACAATGCATGACAATTTACAAATAGCTAAATTTTTGGAGCTAACTTTAAGTACCTGAATGGAATTTAATCAACCCACTAATCTCCTTCTCACTTCTCAGTTATTTATCAAGTTTATGTCAAGGGACAAGGAAAAATTATCCAAACATTGTTTAAAACAATCATCATTAATTAGTAACACTTATCCAGGGGGGTTTTTAACCTTTCCCCCACTCAAGGATTATTCTAATGTCAGAGTAGAATAAAAAATAAGTGCAGTGATGCTGACTCTTCCAAGCTTAACATTTCTCACAAGTCAATTAGCTTTGTACTGGGAGGAGGGCGTGAAGGGCTGCTTGCGGTAGTTGTGTAGCAGCAGCACAATGGCCGCAGACAAGGAAAACAGTTTCTAGGAATTCCTCGTATATAATTTTATATTTTTGACAAGATTAATGACCCATGCTCCCTTCCTCTCCATTTCTTTTTTTGGAATTCTGTTGGTATGTAGTTACTATATTTTATTAAAGGAAATTAGCCTTATCTCTTATT
//>node_2  cylcle:false length:566 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:2 fst:1 lst_coverage:2 lstr:1 lstf:2
//AGGAAATTAGCCTTATCTCTTATTATATTTTTTATGACCTTCAAAGTAGTGTCTCTGCTTAAAAGTGTACCCTGGCCGGGCGTGGTGGCTCACACCTGTAATTCCAGCACTTTGGGAGGCCGAGGCGGGTGGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTGTACTAAAAATACAAAAAATTAGCAGGGCATAGTGGCGGGCGCCTGTAGTCCCAGCTACTCAGGAGGCTCAGGCAGGAGAATGGCGTGAACCCGGGAGACGGAGCTTGCGGTGAGCTGAGATCGCACCGCTGCACTCCAGCCTGGGCGACAGAGCAAGACTCCGTCTCAAAAAAAAAAAAAAAGTATACCCTGAGGCACACATCAAGCGACATGTAGAGTTCATAAATTCTGGCCAAATGGTCATACCTCAAACCTCATCAGCAGTAAGGCTCTTTACTTGCACTGACAAATATGAACGCTGGGGAATTTGGAAATGATATATAATATATAATATTATATATATAATAGATATATAATATATAATATATATAATACATAT
//>node_3  cylcle:false length:2989 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:1 fst:0 lst_coverage:2 lstr:1 lstf:2
//AATTTGTTACAGCAGTAATAGAACAAGTGGTTATCCATATGAGGCAAATTAGATTGGATACCTATCTCCAATAGAAATCAATTCAAGGTGAATTCCAGGAAAATACTTAAAACATTTAGATTAAAAATAAATGAGAATTTTTGTTACTTTTGGTAGGTCATAGAACCAAGAAAAACAAACATTAAGGAGGAAAAATGAACATATGACTACATCAAAATATAAAGCTTCTCTATTTGGAAGATATCATAAGGTGACAAATCATAAACTGTAATATTTACAACATATATATAAGTGAATAAATATACATTTAGAATATATATGAACTCCCAAAAATCAACAGGAAAAATAAGACATAGAACAAGCAAAATGCATAAACAAAAGAAGGCAAAACAAAAATAATGACTCATAATTATATGAAAAGAAGCTCATCTTCATAGATGAGCAGATAAATGCAAATTAAAACCACCCTGAGATGCTTTTTACATCCATGAGCCTGATAAAAGTTAGAGTCTAAAAGTAATAACAAAGATGGGAAGTAATAGAAAATCTTGTCCATTACTGGTTAAAGTATAAACTGATACAGCTACTTTATAGAATATTACATTATAGAATAAAGTTGTGAGTATGTATATGCAGTGACTCAGCATATTCATTGCTAGTATGTACTCAAGAGAAACTTACAGGAGTGGACTAGGAAGTAAATACAAAATGATTACAACATTGTTTGTTATATCAAAAAATAAAAAAGACACCCAATTTTCCAGCAAAAAAAATAAGTAAAAATAAATCCTGGTGTATTCTAACAATGGAATAATATATAGCCATTAAAATAAATCAACTATTACTGTACATATGAATGTAAATATCAGCAAAACATATTGTTTAGTGAAAAACTAAGAAGCTGAAGAAGAATATATACAATATGGTTACATTTATATGAAGTCCAAAAACTTGCAAAATAAAGAAATGTATTTAGAAATAGATTCACATGTGAGAAAACTAGAAGAAAATTAATGAAAGGATAAGAGGGATAGCAGTAATTCTGAGTAGTTGAGGGAATTTCAATTGGAAAAAAATAATATCATATTCTTTAAGTCAGGTAGTGGGTATTAGCATTTGTTTTACCATCGTTCTTTATTCTTATAGCTACACTATATATTTTCAATGTATTTAATGTATTTTTTGCATAATTAAATATTATGCAATAAAAATGAGAAAACAAAAAAGTAGAAAATGATAAATTACAATAAAGAAATGGAGAAAAAATTATAATCTAGTTGAGTAATGGTATATTACATAGCTATTTTCTTAAGTAGATGTATGTACATGATGTATGCACGATTGTACATACATGTTCTTAATTATATATAAATATATATGTACATATTTTTAATATAAAATACTAAACAAAGTACACCAAAATATTAGCTCCTATGTTAGTGAGATAATGTTTTGTTTTTTTGTATTTTAAGTTTTACATAGTAGGTGTATTTTTCTGTTTTCATACTGCTATAAAGAACTGCCCAAGACTGGGTAATTTATAAAGGAAAGAAGTTTAATTGGCTCACCGTTCAGCACAGCTTGGGAGGCCTCAGGAAATCTACAATCATGGCGGAAGACAAAGAGGAAGCAAGCCAGCTTCTTCGCAAGGCAGCATGAAGAAGTGCCGAGCAAAGGGGAAAGAATCCCTTATAAAACCATCAAATCTCGTGAGAACTCACTATCACAAGAACAGCACAGGGGAAACTGCCCCCATGATTCAATTACCTCCACCTGGTCTCTCCCTTGACCTGTGGGGATTATGGGGACTATGGGGATTACAATTCAAGATGAGATTCAGGTGGGGATACAAAGCCTAACCATATCAGTAGGCATGTATTGAATTTTAAACTCAGAGAAAAATACTAGTGTTTTTATAGGATTCTTACTAAAGAAAAACCAGAAAGTAATAAACCATCTACGCTAAGACATAAAATTCAGTTGTTTAGTTACAAGATAGAATGTGGCCTTGTAAGAAAGCAAATTAACTTCTAACATACAAAGCCTTAGAGAAGATTCAAGTGACTGACGGATCTTAAACAGAGCTATTATTACAACTCGAACTGCAGTAAAATATCCTCAGCAACATAGATGTGTGTGTTTCACTAGTCAGAGCAATACAAATTTAATGAAACTCCACTGGTGGTGTTTTTAATCAGACAATTTCTGAAGATGTCCTGGCTTATTCACAGATGCAAGCCAAATCTCTAGAAGAGTACCATAATAAGAAAAAAAAGAATACAGGCAATTGAGAGCTGTTCCAAAGTTTAGGGAGTTTTTGTAAGGAATTAATAAATAAAAATGTTCTTGAAAGACAGAAATTAATATGCAGTTCATACTGTCAGAATTGCAGGCAATTTATCAAAGTCCCCTAATCCTCCAAAATCGCTATTTTTTTTTTGACACACACTTTACAGTACAGAAGAAAATGTCTCCGGCAATAAATCACAAAGTTAAAATTACCTAGTCTACAATTAACTACACAGTGATGGTAAATCATTTTCTACCAAAAGAAAGAAATGTCTTGTCTATTCAGGTTCTGCTCTACTTAAAAGTTTTCCTTGTTGGCGAGCAAGTGGTTAGAAAATTATATTTTATACGTACATTCAGCTTAACTATCATTCAGCTCAGGAAGATGACTCAGGGCCTTATCCATACCTTCAAGTTTGCTCTTAGCAAGTAATTGTTTCAGTATCTATATCAAAAATGGCTTAAGCCTGCAACATGTTTCTGAATGATTAACAAGGTGATAGTCAGTTCTTCATTGAATCCTGGATGCTTTATTTTTCTTAATAAGAGGAATTCATATGGATCAGCTAGAAAAAAATTAAGAGGAAAATCACATGGAAAGTTATATATTATATATCTATTATATATAATATATATCTATTACATATTATATATTGTATATCTATTACATATATATTATATATGTATTATATATATTATA
//>node_4  cylcle:true length:1691 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:2 fst:2 lst_coverage:2 lstr:2 lstf:2
//TGGCCAGGCTGGTCTCGAACTCCTGACCTCAAATGATCCATCTACCTCGGCCTCCCAAAGTGCTGGAATTACAGATGTGAGCCACAATGCCCGGCCTTATTTTCTACAACTTTGGTAACTTTAGCATATACCCCAAATCTGTAAGACATAATATTATAATTCAAATGCAACTCATGGCTTCTCTTTGTACTCTTTCTCTAGCTTTTGAATTATTTATTCTAATACCAGTTTTAATTCTGACACAAAATCATGGGAGTTCTAATCAAAATCCAACCTTTTATCATAAAAACTATGAAGAAATTATGAGTAGAATTTAAAAAGGAAAATAGGCCTATTAATTAGATTTGTCTTTGTAGCATTTAACTCTATAATAAATAATATTTTATGCCTATGAGTCCCCAACAAAGCCTCCAGCTTCTATTTAGATATAAACTGTAAAAGTCACTACTGGATCCACAAGCAAGACTATGGTAAATAAATTTCTCCACCTAACCAGCTTCTTTTACATGATGTTACATGTTTCTTTTGTTTTTTCATTTTGGCAAATATTGATTGTCATCTTCGTGTTTGTCTATGTCCTAAGTGCTGGGATACAGAATCTGAAAAGATGGACACAGGACCTGCCTTCAAGTTCACCCCCTTTTTTTTTTTTTTTTGAGATGCAGTTTTGCTCTTGTCACCCAGGCTGGAGTGTACTGGTGAGATCTCTGCTCACTGCAACCTCCACCTTCAGGGTTCAAGTGATTCTCCTGCCTCAGCCTCCCAAGTAGCTGTGATTACAGGTCCCAGCCACCACGCCTAGCTAATTTTTGTATTTTTAGTAGAGACAGCGTTTCATCATGCTGGTCAGGCTGGTCTCGAACTCCTAACCTCAGGTAGTCGACCCACCTCGGCCTCCCACAGTGCTGAGATTACAGGCATGAGCCACCACGCCCTGCTAGGAGTTCACGCTTTAGTTGGGGAAAATATACAATAAGCAAGCCAGTTTTTAAAATGAGAACTGCAATTAGAGTTAAATGCTACAAAGACAAACTCACAGGAAGATGGGATGTAGAATGATAAGGCTCTCAGAATAGTAAGAGAAACTATTGCTTCTTACGATGTTTGTCTTTCTTTGTATCGGTGCTCAGCTGAGTCTGCAGTGCTTCAGAGGCAGCTTTCATTTTATAAAAATCTATGATTTCTCCTTCCAGTTTTTTTTTCTCTTCCTCGAGCTTCCTTATCTCCTCCTGTTGAATCATTTTAAGATGCTCGAACTTGTCCTGCAGCTGTGAAACCAATGTGCAGTTGTGACACCAAAGCAGTGTGGCTGAACACCTAAAAGAATACGCTTTTTTTCTGATTATCAAACAAACCCAAATCATCACAGTAGACCACGATCTTAATAACAATCTCAAAAACTCAGGAGTAAACACTCAGATATGGAATTTTTCTTTTCTTTCTTTTTTCCTTTTATAAGATGGAGTCTCACTCTGTTGCCCAGGCTGGAGTGCACTGGTGCGATCTCAGCTCACTGCAACCTCCATCTCCCAGTTCAAGTGATTCTCCTGCCTCAGCCTCTTGAGTAGCTGGGACTATAGGCATGCACCACCACTACAGGCGTGTGCCACCACACCTGGCTAATTTTTGTATTTTTAGTAGAGATGGGGTTTTGCCATGATGGCCAGGCTGGTCTCGAACTCCTGACCTCA
//>node_5  cylcle:false length:462 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:2 fst:2 lst_coverage:1 lstr:0 lstf:1
//TGGCCAGGCTGGTCTCGAACTCCTGACCTCAGGTGATCCTCCCGCTTTGGCCTCCCAAAGACTTTTTTTTTTTTTTTAATATAGAGACAAGTTCTCAGTACGTTGCCCAGGCTGGTCTCAAACTCCTGAGCTCAAGTGATCCTCCCACCTCAGCTTCCCAAAGTGCTGGGACTGACTGGATGCAGTGGCTCATGCTTGTAAACTCAGCACTTTGGGAGGCCAAGGTGGGAGGATCGCTTGAGCCCAGGAGTTCAAGACCAGACTGGGTGATATAACACAATAGTCAACTTCAACAGGAGAGAGAATCTGTAAACTTGAATATAGATCTTCCGAAATTATCCAGTCAGTGGACAGAGAAAAAAAGAATAAAAGAGAGAAAAGAAGGCTGGGTGTGGTGGCTCAAGCCTGTAATCCCAACACTTTGGGAGGCCGAGGCAGGCAGATTAAGAGGTCAGGAGTTCA
//>node_6  cylcle:false length:116 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:2 fst:1 lst_coverage:2 lstr:1 lstf:2
//AATAAGAGATAAGGCTAATTTCCTTTAATAATAATAAAATCCTTTAATAAAAATATAAAGGAATAATATAATAATTTTCTTTAATAAAATATAATAAGAGATAAGGCTAATTTCCT
//>node_7  cylcle:false length:41 min_coverage:2 max_coverage:2 average_coverage: 2.00 fst_coverage:2 fst:2 lst_coverage:2 lstr:2 lstf:1
//ATATATAATATATAATATATATAATACATATATAATATATA
//>node_8  cylcle:false length:38 min_coverage:2 max_coverage:2 average_coverage: 2.00 fst_coverage:2 fst:2 lst_coverage:2 lstr:2 lstf:1
//AAAATATAATAAGAGATAAGGCTAATTTCCTTTAATAA
//>node_9  cylcle:false length:48 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:2 fst:1 lst_coverage:2 lstr:1 lstf:2
//TATAATATATATAATACATATATAATATATAATATATATAATACATAT
//>node_10  cylcle:false length:100 min_coverage:2 max_coverage:2 average_coverage: 2.00 fst_coverage:2 fst:2 lst_coverage:2 lstr:2 lstf:1
//AGAATATGAGACATTTCCCTAATCATTATGTGTAATTACAATTACATATATATATGTAATTGTAATTACACATAATGATTAGGGAAATGTCTCATATTCT


//The reference is the same as the individual, except it has about 1kb inserted in node 0. To be clearer, the reverse complement of node 0 is what is in the 10kb without insertion:

//GTGTGGTGACTCATAGCACTCTGTCACCCAGGATGGAGTGCAGTGGCACAATCTCTGCTCACTGCAACCTCTGCCTCCTGGGTACAAGCGATTCTCCTGCCTCAGCTTTCTGAGTAGCAAGGACTACAGGTGCACACCATCACGCCTGGCTAATTTTTGTACTATTAGTACAGACGGAGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCAGGGTCCATTTATTCAATGAATAATTGCCGCTATGTGTCAGACATTTTTCTAGGCCTAGGAATGGATACATAAGTGAACAAAGCAAAGATTCTGGTTCTTGTAGAGTTTCCATTAAAAGACAATTTAGTAAAACTTTTCTTCCCCCAAATTATAAAATCTGTAAGATGATTTAACAACATGTGTAAAAGTCATTGTGGGCCAGGCACGGTGGCTCATACCAGATTTTTTACAAGGTAGATATTAATTTTAGATATGGAATAATATTGGTGATTTCAATTTTATAACACTGGGTTAAGATGAAAGAATGAGAAGATAAAGGTCCCTCAGCAATATAACTCACAAACATGTTCAGAAGCAGTAAGAAGTTACATTAATTATCTTTTGAAAGTCGATAATCTACATCTTTAATGTATGCATATAGCATAGCTAATGTACTATCCCTAAAGTTAAAAGTTATATAGGGAACAAACACTAATTTTTTTTAGAAAAAATTATAAAAAGAGTAAAAATATGCCTTATACTACCGTAATTTCATGTTTTACAGCTCTGGGAAAATAGAAAATAAAATGTTCTGTTAGCATGAATCCCTCTGTGCCCCCAAAAAACCCTATGGATTGCATCATTATTACCTAAAAAGTCTATTCTCAAATGCAGCAGAGTGATTAACCAAAAATGAAATATACATAGAATATAAACAAAAGGAAATGAGACTAGAAACAAAATGTGTCACTACAAAAAAATCAACTAAAGATAAAAAAGAAATAATTGAGAAAATGATTGGCAAAAATCAGTAACTCTGACGTATTAAAACTTTCCATGCTACATAAATCTGAAAACTCTATTTCACATAAAACTGGAGCTGAAAGAAACAAATATTTACCTATTTTTTATTTTCTACATAATTTAAAAGGCAAATGCATAAAATGTAATTGTAAATCTGTTAGCTGGTATACAATGAATAAAGATATAATTTGTCACATCAATAACATAAAAAGAGTAGAGCTATATATATAGCAGTAGAATTTTGGTATGTGATTGAACTTAAGTTGAAATAAATTCAAATTAAAATGTTATAACTCTAGGATGTTATATGTAATTCTCATAGAGACATTTCCGGATTTTTTTTTAAAACTGAAAAAAATCCATTTATCCCTGAATTTGACATTCAGGAAGTGTTAAGTCCTTCAGGTTGAAATAAATGAACTCTAGGCAATAACTATATAAGTAAATAAGCAAGCTGTATGAATATACAAAGCTCTCTGGTAAAGGTAAATACATAAACAAACATAAAAACAGTCCTATTGTAATTTTGGTTTGTAACTCTGCCTCAAAAACAAACAAAAAAAAGATAATGGGAGAATGTTGAAAACTCAGAGAGAAGAGCAACTCTCACAGATAGGGATCCAGATAACATTAGCAGCTGATTTCTCAGCAGAAACCTTGAAGGCCAGTAGGCAGTGGATTATATATTTAAAATAATGAAGAAACCTGTCAATTGAGAAATATATAGCTGGAAAACTTATCCTTCAAAAATGAAGGAGAAATTAATATGACTTTTTCAGAATCAAAAAATTTAAATTTCTGTAATAAAATTTAAATGTTTATAAATTTAAAAAACTAGAAGAAAGAATGTTGACTGTTCACAATACAAATAAATGACAACTATTTGAGGTGATGGATACGCTAATTATCCTTATTTGATCACTGGGCATTGTATACATGTATCAAAATATCACTCTGTATCCCATGAATATGTACAATTATTTGTAAGAAATTATAGAAAGGAAACAAAATAAACAGATTATATGGAGGATTTTTAGAAGATAAGTAAATAAATTAATATACTAGGAAAAAACAAGGGAAATATACTTGATAAATAAATACAGGTAAGAGTTCTTTTGAAATAATGATAAAATAGAAAATCTCTGTCAAAACTAAAAGGAAAGATGCATAAATATATAAATAAACGATAAAAAATGTTGCATACATCTATAGACCAGTTTTTTATTTAAATTGAAGATTAGTATACATTTTAAATGATTAGTCAAAATAAAAAATCTAAAATGTGCTCTAAATACCTCTTAGGTCAGAAAAAAAAAGTCAAAAGCTAGAGTATAGAGAAATTAAGAAACGCCCTAAATTTCTAATCTGACAAAAATTCATACAAGATTTAAATATTTTAATGGAAAATAGAACAGAACTAATTATTGAACAGTTGCTTGGTTTACAAACATTAGAACAATTTTCTTATTCACACCATCTGATTATTGTATGTTATTTTTTCCCCAACGTTTAGACTACACAATGAGTTAAGAATGATAAAAATAAGCTCACCAATATACTATGTACATATTTACCAAAATCTGTGCATGCTTATACATATAAACACAGCTGATAATTTATTAGTTAGGCTCATTTGTAATTTTTGTCAATGGCTCAGTGTGATTGAAAGTTCTGTGTGAATATGTTTTTGGAAAGATCCAACAGCAACACCTTAGTGTATGTTTTTGAAATAAAATATATCTGAGTAGCAGCAAAGTTATTCTCAAATTTCCATTTTATAGCTGGAGATGTTATACCGTGACGTACATGATAGGACCCAATATGGATCAATCCCTTTTAGAAGTCAATCAGGAAGAGGGGAGCAGTTAATGAAGGATCATGTCAAAGGAAAGAAAAATTAGGGGAACATTAAAAGCTTTCTTCCCAAGCCACTAAATCAACTTGACTAACAAAATTACCACTTGATTTAGTATTAGAAAATTACATTACATATCAAACATAAACCCATTAATCAAATACTAAAGAAATTTCTGAGTTAAATGGTATAATGTTAGCTTATGCCAGAGCTGACCTTGAAAGATTGTTCAAATTCTGTCTCTACAAAAAATACAAAAATTAGCCAAGTTTGGTGGCATGTGCCTGTAGTACCAGCTACTTGGGAGGCTGAGGTGGAAGAATAGCTTGAGTCTGGGAGGTCAAGGCTGCAATGAGCTGTGA  TTGCACCACTGCACTCAAGCCTGGGTGGTAGAGTAAGACCCTGTCTCAAAAAAAAAAAAAAAAAAAGAAAAATCACTAAGCAAAATAAGACATGGTATATAGAAAATTGGCTAAGTTGTGGTCCATTCATGTTTGCTCCCAATTAAGGAGCACAGCTATGAAAAGGAAGGCTTCAAATTAATAACCAATAGATTTTTTTAAAAAGAAAACTGGCCAGGTACTGTGGCTTATGTCTGTAATATCAGCATGTTGGGAGGCCAAGGCAGGATTACTTGAGCCCAGAAATTCCAGACCAGCCTGAGAATTTGGCAAAACTCGTACAGCATGTAACACAATTGCACTTTTATATTTTAACAAATTATATTTCCCATATTGAACTGTAAGTCTCCTGAAAGCAGGAATTTTGTTCTTGCTCATCATCAACTTTTTCAACATCCAGTGCACCATTTAGAACTTAGATGTAGTCAATACAGGTTTGTGGAATGAAAGAGGAAAAGAAAGAATTAATATTCCTTTAAATTAGGATGGCAAAGATCTTTAAGCCATATTAGACCTTTCTTCAGTTTTTTATATAGACTTTGTCGCATCACACCTCAGAGATTCTGTACATGTTCTTCCTCCTGCCTAGAAAGGATCGTCCCTCCACTTTTGCCAACTAATCCCTGCTCAACTTTTCATCTCAGCAGGAGGCCCATTCTCTTTGGCAATCCTCTGGCCTCCAGCCCATTTATTATATGCTCACATGTCAACATGTACTACAGTGGGCATGTAATGAATATTTGTTAAGTGAAGTAATAAATTCAATCAAATCACATCACCTGTTTAAAGCACTTCATTGGCTTCACATTGCACTTAGAATAAAGAGAAATTCTTTTTATACAATATACAATATATTTTATACAATATAAGTTCCTGCAGAATGCAGACACTTTCTACTTCTCCAGCCTCTTTTCGACTCCTCTCCTACTAGCTTCTGTAAATTGCTTCATCATGCCTGCTTGCAATTATAGGGCTTTTGAACCATTTGTTCCCTCCTTAGATATCCTCATTTTTTTCAGATTCTTGCTTAGAAGTCACTCCTCCGTGGACCTCCTCTGACATATTAAACATTGCAGTCCATTATAAGCTGCAAGAGGACAGGGATTTTTGCCTGTTTTATTCCCTACTGTATCACCAGGGGCTAGAGCAATATCTGACAAACATAATGATTAGGGAAATGTCTCATATTCTATATATATAGACAGAAAGAGAGAAAATATATGAGGGAGAGAAGGAATCTTTCCATCTCCTTTGAGTTCCACGGTGTTGAGAGTCAGGACAACTGC

//  and the insert happens just after TTGGGAGGCTGAGGTGGAAGAATAGCTTGAGTCTGGGAGGTCAAGGCTGCAATGAGCTGTGA, where I have put a space above


//The inserted sequence is

// ACCTCAATCTATAACCACTACCTTCTGGGTCCTTTCTAAAAATTGACAAATAATAATCATATATAATTAATGTACAATGTTATGTTTCAATACATGTTTGCATTGTGGAATAATTAAATCAAGCTACTTGGCATGTCAGTAACATCACATGCTTAGCATTTTTGTAGTGAAAATATTTAAAATCTACTCTTTTAGCAATTTTGAAATATACAATACAGTACTTACTCACTTAACATCATTGATTGGTCCTCAGAAACTGCAACTTGGAGTGAAAAGATGTATAAAGAAACCAATTTTCCCATAGGCTAATAGCTATAAATAAGAGTTAGATTCTTACAGCATATTTTTGGTCACAAAATATCACCAAACTTCTAAATAAAGACCAAAACACTTCAAATATTAAACATTGAAATAAATATGAGCTTTGCATACATTTAAGAAAGATTAATAAAAACAAGTAAGATAATTATTTGCCCAATTATTTCATTCAGGGTTGGGGAGACTGGAGTCTGTGCTGGAAGCTCAGGGCTCAAGCTGGGCAACAGCCCTGGACAGGATGCCATCCCACTGCAGGATGGCTCACACATGCCCACAGCCACTCAGCCTGGGACCATTTGGACACAGCAATTAACCTTACCTGCATGTCTTTGTGGGGAGGAAACCAGAGTGCTTAGAAAAACCCATGCAGACAGACACAGAGCAAACATGCAAACCTCACAAAGATATTGTTTCTTCTGTCACCTGTGCTTTTGGGTCATATTCAAGAAATCATTAACCAAATAAAAGTCGTGGAGCTTTTCCCTATGTTTTCTTTTAGTAGTTTTATAGTTTCAGGTCTTACATTTAACTCCTTAATCCATTTTGATTTTTGCATATGGTGTGAGATAAGCTTCTGGTTTCATTCTCCCACATGTGGATATCCAGTTCTCTGAACACCATATATGGAAGAGACTGTCATTTCCTCATGATATGTTCCTGGCACTTTTGTTGAAATCAATTGACCATAGATCTGTGGGTTTATTTCTGGCTTTTTATTCTGTTCCATTGGCCAATGTACCTGTGTTTATGCTTGTGCCTTGCTGTTTTGATTATTATAGCTTTATAATATGTTTTGAAATCAGGTAGTGTGATGCCTCCATCTTTGCTTTTTATGCTCAAGATAGTTTGGATATTCAGAGTGTTTTATGGTTCCATATACAT
// which I took from chromosome 2.






//NOTE - the inserted sequence starts with an A. The sequence that follows the inserted sequence also starts with an A. So our algorithm will add that A to the flanking region.
// We allow for this is checking results match expectations



  int  min_fiveprime_flank_anchor = 21;
  int  min_threeprime_flank_anchor= 21;
  int  max_anchor_span =8000;
  int  length_of_arrays=16000;
  int  min_covg =1;
  int  max_covg = 100000000;
  int  max_expected_size_of_supernode=8000;


  char** return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  char** return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  char** return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  char** return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      die("Failed to alloc return_something_array - cannot start test\n");
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 5600);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 5600);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 5600);
  return_branch2_array[0] = (char*) malloc( sizeof(char) *5600);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      die("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';


  int return_variant_start_coords_array[2];            
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  int* return_variant_start_coords_array_ptr[2];      
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);

  FILE* fp = fopen("../data/tempfiles_can_be_deleted/temp_outputfile_trustedpath_sv_caller_test8", "w");





  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, 0, 
							1,
							min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
							max_expected_size_of_supernode, length_of_arrays, hash_table,  fp,
							1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, 
							return_variant_start_coords_array_ptr, &make_reference_path_based_sv_calls_condition_always_true, 
							&action_set_flanks_and_branches_to_be_ignored,
							&print_no_extra_info, NULL, NoIdeaWhatCleaning);


  fclose(fp);


  CU_ASSERT(ret==1);

  CU_ASSERT_STRING_EQUAL("CATAATGATTAGGGAAATGTCTCATATTCTATATATATAGACAGAAAGAGAGAAAATATATGAGGGAGAGAAGGAATCTTTCCATCTCCTTTGAGTTCCACGGTGTTGAGAGTCAGGACAACTGCAATTGCTTCATCATGCCTGCTTGCAATTATAGGGCTTTTGAACCATTTGTTCCCTCCTTAGATATCCTCATTTTTTTCAGATTCTTGCTTAGAAGTCACTCCTCCGTGGACCTCCTCTGACATATTAAACATTGCAGTCCATTATAAGCTGCAAGAGGACAGGGATTTTTGCCTGTTTTATTCCCTACTGTATCACCAGGGGCTAGAGCAATATCTGACAAACAGTGGGCATGTAATGAATATTTGTTAAGTGAAGTAATAAATTCAATCAAATCACATCACCTGTTTAAAGCACTTCATTGGCTTCACATTGCACTTAGAATAAAGAGAAATTCTTTTTATACAATATACAATATATTTTATACAATATAAGTTCCTGCAGAATGCAGACACTTTCTACTTCTCCAGCCTCTTTTCGACTCCTCTCCTACTAGCTTCTGTATTTAAGCCATATTAGACCTTTCTTCAGTTTTTTATATAGACTTTGTCGCATCACACCTCAGAGATTCTGTACATGTTCTTCCTCCTGCCTAGAAAGGATCGTCCCTCCACTTTTGCCAACTAATCCCTGCTCAACTTTTCATCTCAGCAGGAGGCCCATTCTCTTTGGCAATCCTCTGGCCTCCAGCCCATTTATTATATGCTCACATGTCAACATGTACTTCGTACAGCATGTAACACAATTGCACTTTTATATTTTAACAAATTATATTTCCCATATTGAACTGTAAGTCTCCTGAAAGCAGGAATTTTGTTCTTGCTCATCATCAACTTTTTCAACATCCAGTGCACCATTTAGAACTTAGATGTAGTCAATACAGGTTTGTGGAATGAAAGAGGAAAAGAAAGAATTAATATTCCTTTAAATTAGGATGGCAAAGATCGTATATAGAAAATTGGCTAAGTTGTGGTCCATTCATGTTTGCTCCCAATTAAGGAGCACAGCTATGAAAAGGAAGGCTTCAAATTAATAACCAATAGATTTTTTTAAAAAGAAAACTGGCCAGGTACTGTGGCTTATGTCTGTAATATCAGCATGTTGGGAGGCCAAGGCAGGATTACTTGAGCCCAGAAATTCCAGACCAGCCTGAGAATTTGGCAAAACTCTGTCTCTACAAAAAATACAAAAATTAGCCAAGTTTGGTGGCATGTGCCTGTAGTACCAGCTACTTGGGAGGCTGAGGTGGAAGAATAGCTTGAGTCTGGGAGGTCAAGGCTGCAATGAGCTGTGA", return_flank5p_array[0]);

  CU_ASSERT_STRING_EQUAL("GAGTAAGACCCTGTCTCAAAAAAAAAAAAAAAAAAAGAAAAATCACTAAGCAAAATAAGACATGTGAAGGATCATGTCAAAGGAAAGAAAAATTAGGGGAACATTAAAAGCTTTCTTCCCAAGCCACTAAATCAACTTGACTAACAAAATTACCACTTGATTTAGTATTAGAAAATTACATTACATATCAAACATAAACCCATTAATCAAATACTAAAGAAATTTCTGAGTTAAATGGTATAATGTTAGCTTATGCCAGAGCTGACCTTGAAAGATTGTTCAAATATGGCTCAGTGTGATTGAAAGTTCTGTGTGAATATGTTTTTGGAAAGATCCAACAGCAACACCTTAGTGTATGTTTTTGAAATAAAATATATCTGAGTAGCAGCAAAGTTATTCTCAAATTTCCATTTTATAGCTGGAGATGTTATACCGTGACGTACATGATAGGACCCAATATGGATCAATCCCTTTTAGAAGTCAATCAGGAAGAGGGGAGCAGTTAAAACAGTTGCTTGGTTTACAAACATTAGAACAATTTTCTTATTCACACCATCTGATTATTGTATGTTATTTTTTCCCCAACGTTTAGACTACACAATGAGTTAAGAATGATAAAAATAAGCTCACCAATATACTATGTACATATTTACCAAAATCTGTGCATGCTTATACATATAAACACAGCTGATAATTTATTAGTTAGGCTCATTTGTAATTTTTGTCACTATAGACCAGTTTTTTATTTAAATTGAAGATTAGTATACATTTTAAATGATTAGTCAAAATAAAAAATCTAAAATGTGCTCTAAATACCTCTTAGGTCAGAAAAAAAAAGTCAAAAGCTAGAGTATAGAGAAATTAAGAAACGCCCTAAATTTCTAATCTGACAAAAATTCATACAAGATTTAAATATTTTAATGGAAAATAGAACAGAACTAATTATTGAAGAAATTATAGAAAGGAAACAAAATAAACAGATTATATGGAGGATTTTTAGAAGATAAGTAAATAAATTAATATACTAGGAAAAAACAAGGGAAATATACTTGATAAATAAATACAGGTAAGAGTTCTTTTGAAATAATGATAAAATAGAAAATCTCTGTCAAAACTAAAAGGAAAGATGCATAAATATATAAATAAACGATAAAAAATGTTGCATACATATATGACTTTTTCAGAATCAAAAAATTTAAATTTCTGTAATAAAATTTAAATGTTTATAAATTTAAAAAACTAGAAGAAAGAATGTTGACTGTTCACAATACAAATAAATGACAACTATTTGAGGTGATGGATACGCTAATTATCCTTATTTGATCACTGGGCATTGTATACATGTATCAAAATATCACTCTGTATCCCATGAATATGTACAATTATTTGTCTCAAAAACAAACAAAAAAAAGATAATGGGAGAATGTTGAAAACTCAGAGAGAAGAGCAACTCTCACAGATAGGGATCCAGATAACATTAGCAGCTGATTTCTCAGCAGAAACCTTGAAGGCCAGTAGGCAGTGGATTATATATTTAAAATAATGAAGAAACCTGTCAATTGAGAAATATATAGCTGGAAAACTTATCCTTCAAAAATGAAGGAGAAATTAAGACATTTCCGGATTTTTTTTTAAAACTGAAAAAAATCCATTTATCCCTGAATTTGACATTCAGGAAGTGTTAAGTCCTTCAGGTTGAAATAAATGAACTCTAGGCAATAACTATATAAGTAAATAAGCAAGCTGTATGAATATACAAAGCTCTCTGGTAAAGGTAAATACATAAACAAACATAAAAACAGTCCTATTGTAATTTTGGTTTGTAACTCTGCTTTTTATTTTCTACATAATTTAAAAGGCAAATGCATAAAATGTAATTGTAAATCTGTTAGCTGGTATACAATGAATAAAGATATAATTTGTCACATCAATAACATAAAAAGAGTAGAGCTATATATATAGCAGTAGAATTTTGGTATGTGATTGAACTTAAGTTGAAATAAATTCAAATTAAAATGTTATAACTCTAGGATGTTATATGTAATTCTCATAGTAACCAAAAATGAAATATACATAGAATATAAACAAAAGGAAATGAGACTAGAAACAAAATGTGTCACTACAAAAAAATCAACTAAAGATAAAAAAGAAATAATTGAGAAAATGATTGGCAAAAATCAGTAACTCTGACGTATTAAAACTTTCCATGCTACATAAATCTGAAAACTCTATTTCACATAAAACTGGAGCTGAAAGAAACAAATATTTACCTATAAAGTTAAAAGTTATATAGGGAACAAACACTAATTTTTTTTAGAAAAAATTATAAAAAGAGTAAAAATATGCCTTATACTACCGTAATTTCATGTTTTACAGCTCTGGGAAAATAGAAAATAAAATGTTCTGTTAGCATGAATCCCTCTGTGCCCCCAAAAAACCCTATGGATTGCATCATTATTACCTAAAAAGTCTATTCTCAAATGCAGCAGAGTGATATTTTTTACAAGGTAGATATTAATTTTAGATATGGAATAATATTGGTGATTTCAATTTTATAACACTGGGTTAAGATGAAAGAATGAGAAGATAAAGGTCCCTCAGCAATATAACTCACAAACATGTTCAGAAGCAGTAAGAAGTTACATTAATTATCTTTTGAAAGTCGATAATCTACATCTTTAATGTATGCATATAGCATAGCTAATGTACTATCCCTGGGTCCATTTATTCAATGAATAATTGCCGCTATGTGTCAGACATTTTTCTAGGCCTAGGAATGGATACATAAGTGAACAAAGCAAAGATTCTGGTTCTTGTAGAGTTTCCATTAAAAGACAATTTAGTAAAACTTTTCTTCCCCCAAATTATAAAATCTGTAAGATGATTTAACAACATGTGTAAAAGTCATTGTGGGCCAGGCACGGTGGCTCATACCAGGTGTGGTGACTCATAGCACTCTGTCACCCAGGATGGAGTGCAGTGGCACAATCTCTGCTCACTGCAACCTCTGCCTCCTGGGTACAAGCGATTCTCCTGCCTCAGCTTTCTGAGTAGCAAGGACTACAGGTGCACACCATCACGCCTGGCTAATTTTTGTACTATTAGTACAGACGGAGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCA", 
			 return_flank3p_array[0]);

  CU_ASSERT_STRING_EQUAL("TTGCACCACTGCACTCAAGCCTGGGTGGTA", return_branch2_array[0]);

  CU_ASSERT_STRING_EQUAL("CCTCAATCTATAACCACTACCTTCTGGGTCCTTTCTAAAAATTGACAAATAATAATCATATATAATTAATGTACAATGTTATGTTTCAATACATGTTTGCATTGTGGAATAATTAAATCAAGCTACTTGGCATGTCAGTAACATCACATGCTTAGCATTTTTGTAGTGAAAATATTTAAAATCTACTCTTTTAGCAATTTTGAAATATACAATACAGTACTTACTCACTTAACATCATTGATTGGTCCTCAGAAACTGCAACTTGGAGTGAAAAGATGTATAAAGAAACCAATTTTCCCATAGGCTAATAGCTATAAATAAGAGTTAGATTCTTACAGCATATTTTTGGTCACAAAATATCACCAAACTTCTAAATAAAGACCAAAACACTTCAAATATTAAACATTGAAATAAATATGAGCTTTGCATACATTTAAGAAAGATTAATAAAAACAAGTAAGATAATTATTTGCCCAATTATTTCATTCAGGGTTGGGGAGACTGGAGTCTGTGCTGGAAGCTCAGGGCTCAAGCTGGGCAACAGCCCTGGACAGGATGCCATCCCACTGCAGGATGGCTCACACATGCCCACAGCCACTCAGCCTGGGACCATTTGGACACAGCAATTAACCTTACCTGCATGTCTTTGTGGGGAGGAAACCAGAGTGCTTAGAAAAACCCATGCAGACAGACACAGAGCAAACATGCAAACCTCACAAAGATATTGTTTCTTCTGTCACCTGTGCTTTTGGGTCATATTCAAGAAATCATTAACCAAATAAAAGTCGTGGAGCTTTTCCCTATGTTTTCTTTTAGTAGTTTTATAGTTTCAGGTCTTACATTTAACTCCTTAATCCATTTTGATTTTTGCATATGGTGTGAGATAAGCTTCTGGTTTCATTCTCCCACATGTGGATATCCAGTTCTCTGAACACCATATATGGAAGAGACTGTCATTTCCTCATGATATGTTCCTGGCACTTTTGTTGAAATCAATTGACCATAGATCTGTGGGTTTATTTCTGGCTTTTTATTCTGTTCCATTGGCCAATGTACCTGTGTTTATGCTTGTGCCTTGCTGTTTTGATTATTATAGCTTTATAATATGTTTTGAAATCAGGTAGTGTGATGCCTCCATCTTTGCTTTTTATGCTCAAGATAGTTTGGATATTCAGAGTGTTTTATGGTTCCATATACATATTGCACCACTGCACTCAAGCCTGGGTGGTA", return_trusted_branch_array[0]);

  CU_ASSERT(return_variant_start_coords_array[0]==6722); 



  hash_table_free(&hash_table);
  fclose(chrom_fptr);
  
  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);


}

void test_db_graph_make_reference_path_based_sv_calls_test_9()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // ===========================================================================
  // 9. Identical to previous test, but each person has identical 600 lines of
  //    chromosome 12 added before the sequence that was there in previous test
  //    This is purely to test that the start coordinate of the variant is found
  //    correctly, even when you have to load more and more sequence into the
  //    array.
  // ===========================================================================
  

  int  kmer_size = 31;
  int  number_of_bits = 15;
  int  bucket_size = 10;
  int  max_retries = 10;

  dBGraph* hash_table = hash_table_new(number_of_bits, bucket_size,
                                       max_retries, kmer_size);

  if (hash_table==NULL)
    {
      die("unable to alloc the hash table. dead before we even started. OOM");
    }

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/variations/two_people_one_with_1kb_deletion_both_with_600lineschrom12_beforehand.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/first_person_600lineschrom12_then_10kb_chrom1_plus_1kb_inserted_mid_supernode.fa", "r");
  if (chrom_fptr==NULL)
    {
      die("Cannot open ../data/test/pop_graph/variations/first_person_600lineschrom12_then_10kb_chrom1_plus_1kb_inserted_mid_supernode.fa");
    }


  int min_fiveprime_flank_anchor = 21;
  int  min_threeprime_flank_anchor= 21;
  int  max_anchor_span =8000;
  int  length_of_arrays=16000;
  int  min_covg =1;
  int  max_covg = 100000000;
  int  max_expected_size_of_supernode=8000;


  char** return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  char** return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  char** return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  char** return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      die("Failed to alloc return_something_array - cannot start test\n");
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 5600);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 5600);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 5600);
  return_branch2_array[0] = (char*) malloc( sizeof(char) *5600);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      die("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';


  int return_variant_start_coords_array[2];              
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  int* return_variant_start_coords_array_ptr[2];        
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);

  FILE* fp = fopen("../data/tempfiles_can_be_deleted/temp_outputfile_trustedpath_sv_caller_test9", "w");

  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, 0, 
							1,
							min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
							max_expected_size_of_supernode, length_of_arrays, hash_table,  fp,
							1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, 
							return_variant_start_coords_array_ptr, &make_reference_path_based_sv_calls_condition_always_true, 
							&action_set_flanks_and_branches_to_be_ignored,
							&print_no_extra_info, NULL, NoIdeaWhatCleaning);


  fclose(fp);


  CU_ASSERT(ret==1);

  CU_ASSERT_STRING_EQUAL("CATAATGATTAGGGAAATGTCTCATATTCTATATATATAGACAGAAAGAGAGAAAATATATGAGGGAGAGAAGGAATCTTTCCATCTCCTTTGAGTTCCACGGTGTTGAGAGTCAGGACAACTGCAATTGCTTCATCATGCCTGCTTGCAATTATAGGGCTTTTGAACCATTTGTTCCCTCCTTAGATATCCTCATTTTTTTCAGATTCTTGCTTAGAAGTCACTCCTCCGTGGACCTCCTCTGACATATTAAACATTGCAGTCCATTATAAGCTGCAAGAGGACAGGGATTTTTGCCTGTTTTATTCCCTACTGTATCACCAGGGGCTAGAGCAATATCTGACAAACAGTGGGCATGTAATGAATATTTGTTAAGTGAAGTAATAAATTCAATCAAATCACATCACCTGTTTAAAGCACTTCATTGGCTTCACATTGCACTTAGAATAAAGAGAAATTCTTTTTATACAATATACAATATATTTTATACAATATAAGTTCCTGCAGAATGCAGACACTTTCTACTTCTCCAGCCTCTTTTCGACTCCTCTCCTACTAGCTTCTGTATTTAAGCCATATTAGACCTTTCTTCAGTTTTTTATATAGACTTTGTCGCATCACACCTCAGAGATTCTGTACATGTTCTTCCTCCTGCCTAGAAAGGATCGTCCCTCCACTTTTGCCAACTAATCCCTGCTCAACTTTTCATCTCAGCAGGAGGCCCATTCTCTTTGGCAATCCTCTGGCCTCCAGCCCATTTATTATATGCTCACATGTCAACATGTACTTCGTACAGCATGTAACACAATTGCACTTTTATATTTTAACAAATTATATTTCCCATATTGAACTGTAAGTCTCCTGAAAGCAGGAATTTTGTTCTTGCTCATCATCAACTTTTTCAACATCCAGTGCACCATTTAGAACTTAGATGTAGTCAATACAGGTTTGTGGAATGAAAGAGGAAAAGAAAGAATTAATATTCCTTTAAATTAGGATGGCAAAGATCGTATATAGAAAATTGGCTAAGTTGTGGTCCATTCATGTTTGCTCCCAATTAAGGAGCACAGCTATGAAAAGGAAGGCTTCAAATTAATAACCAATAGATTTTTTTAAAAAGAAAACTGGCCAGGTACTGTGGCTTATGTCTGTAATATCAGCATGTTGGGAGGCCAAGGCAGGATTACTTGAGCCCAGAAATTCCAGACCAGCCTGAGAATTTGGCAAAACTCTGTCTCTACAAAAAATACAAAAATTAGCCAAGTTTGGTGGCATGTGCCTGTAGTACCAGCTACTTGGGAGGCTGAGGTGGAAGAATAGCTTGAGTCTGGGAGGTCAAGGCTGCAATGAGCTGTGA", return_flank5p_array[0]);

  CU_ASSERT_STRING_EQUAL("GAGTAAGACCCTGTCTCAAAAAAAAAAAAAAAAAAAGAAAAATCACTAAGCAAAATAAGACATGTGAAGGATCATGTCAAAGGAAAGAAAAATTAGGGGAACATTAAAAGCTTTCTTCCCAAGCCACTAAATCAACTTGACTAACAAAATTACCACTTGATTTAGTATTAGAAAATTACATTACATATCAAACATAAACCCATTAATCAAATACTAAAGAAATTTCTGAGTTAAATGGTATAATGTTAGCTTATGCCAGAGCTGACCTTGAAAGATTGTTCAAATATGGCTCAGTGTGATTGAAAGTTCTGTGTGAATATGTTTTTGGAAAGATCCAACAGCAACACCTTAGTGTATGTTTTTGAAATAAAATATATCTGAGTAGCAGCAAAGTTATTCTCAAATTTCCATTTTATAGCTGGAGATGTTATACCGTGACGTACATGATAGGACCCAATATGGATCAATCCCTTTTAGAAGTCAATCAGGAAGAGGGGAGCAGTTAAAACAGTTGCTTGGTTTACAAACATTAGAACAATTTTCTTATTCACACCATCTGATTATTGTATGTTATTTTTTCCCCAACGTTTAGACTACACAATGAGTTAAGAATGATAAAAATAAGCTCACCAATATACTATGTACATATTTACCAAAATCTGTGCATGCTTATACATATAAACACAGCTGATAATTTATTAGTTAGGCTCATTTGTAATTTTTGTCACTATAGACCAGTTTTTTATTTAAATTGAAGATTAGTATACATTTTAAATGATTAGTCAAAATAAAAAATCTAAAATGTGCTCTAAATACCTCTTAGGTCAGAAAAAAAAAGTCAAAAGCTAGAGTATAGAGAAATTAAGAAACGCCCTAAATTTCTAATCTGACAAAAATTCATACAAGATTTAAATATTTTAATGGAAAATAGAACAGAACTAATTATTGAAGAAATTATAGAAAGGAAACAAAATAAACAGATTATATGGAGGATTTTTAGAAGATAAGTAAATAAATTAATATACTAGGAAAAAACAAGGGAAATATACTTGATAAATAAATACAGGTAAGAGTTCTTTTGAAATAATGATAAAATAGAAAATCTCTGTCAAAACTAAAAGGAAAGATGCATAAATATATAAATAAACGATAAAAAATGTTGCATACATATATGACTTTTTCAGAATCAAAAAATTTAAATTTCTGTAATAAAATTTAAATGTTTATAAATTTAAAAAACTAGAAGAAAGAATGTTGACTGTTCACAATACAAATAAATGACAACTATTTGAGGTGATGGATACGCTAATTATCCTTATTTGATCACTGGGCATTGTATACATGTATCAAAATATCACTCTGTATCCCATGAATATGTACAATTATTTGTCTCAAAAACAAACAAAAAAAAGATAATGGGAGAATGTTGAAAACTCAGAGAGAAGAGCAACTCTCACAGATAGGGATCCAGATAACATTAGCAGCTGATTTCTCAGCAGAAACCTTGAAGGCCAGTAGGCAGTGGATTATATATTTAAAATAATGAAGAAACCTGTCAATTGAGAAATATATAGCTGGAAAACTTATCCTTCAAAAATGAAGGAGAAATTAAGACATTTCCGGATTTTTTTTTAAAACTGAAAAAAATCCATTTATCCCTGAATTTGACATTCAGGAAGTGTTAAGTCCTTCAGGTTGAAATAAATGAACTCTAGGCAATAACTATATAAGTAAATAAGCAAGCTGTATGAATATACAAAGCTCTCTGGTAAAGGTAAATACATAAACAAACATAAAAACAGTCCTATTGTAATTTTGGTTTGTAACTCTGCTTTTTATTTTCTACATAATTTAAAAGGCAAATGCATAAAATGTAATTGTAAATCTGTTAGCTGGTATACAATGAATAAAGATATAATTTGTCACATCAATAACATAAAAAGAGTAGAGCTATATATATAGCAGTAGAATTTTGGTATGTGATTGAACTTAAGTTGAAATAAATTCAAATTAAAATGTTATAACTCTAGGATGTTATATGTAATTCTCATAGTAACCAAAAATGAAATATACATAGAATATAAACAAAAGGAAATGAGACTAGAAACAAAATGTGTCACTACAAAAAAATCAACTAAAGATAAAAAAGAAATAATTGAGAAAATGATTGGCAAAAATCAGTAACTCTGACGTATTAAAACTTTCCATGCTACATAAATCTGAAAACTCTATTTCACATAAAACTGGAGCTGAAAGAAACAAATATTTACCTATAAAGTTAAAAGTTATATAGGGAACAAACACTAATTTTTTTTAGAAAAAATTATAAAAAGAGTAAAAATATGCCTTATACTACCGTAATTTCATGTTTTACAGCTCTGGGAAAATAGAAAATAAAATGTTCTGTTAGCATGAATCCCTCTGTGCCCCCAAAAAACCCTATGGATTGCATCATTATTACCTAAAAAGTCTATTCTCAAATGCAGCAGAGTGATATTTTTTACAAGGTAGATATTAATTTTAGATATGGAATAATATTGGTGATTTCAATTTTATAACACTGGGTTAAGATGAAAGAATGAGAAGATAAAGGTCCCTCAGCAATATAACTCACAAACATGTTCAGAAGCAGTAAGAAGTTACATTAATTATCTTTTGAAAGTCGATAATCTACATCTTTAATGTATGCATATAGCATAGCTAATGTACTATCCCTGGGTCCATTTATTCAATGAATAATTGCCGCTATGTGTCAGACATTTTTCTAGGCCTAGGAATGGATACATAAGTGAACAAAGCAAAGATTCTGGTTCTTGTAGAGTTTCCATTAAAAGACAATTTAGTAAAACTTTTCTTCCCCCAAATTATAAAATCTGTAAGATGATTTAACAACATGTGTAAAAGTCATTGTGGGCCAGGCACGGTGGCTCATACCAGGTGTGGTGACTCATAGCACTCTGTCACCCAGGATGGAGTGCAGTGGCACAATCTCTGCTCACTGCAACCTCTGCCTCCTGGGTACAAGCGATTCTCCTGCCTCAGCTTTCTGAGTAGCAAGGACTACAGGTGCACACCATCACGCCTGGCTAATTTTTGTACTATTAGTACAGACGGAGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCA", 
			 return_flank3p_array[0]);

  CU_ASSERT_STRING_EQUAL("TTGCACCACTGCACTCAAGCCTGGGTGGTA", return_branch2_array[0]);

  CU_ASSERT_STRING_EQUAL("CCTCAATCTATAACCACTACCTTCTGGGTCCTTTCTAAAAATTGACAAATAATAATCATATATAATTAATGTACAATGTTATGTTTCAATACATGTTTGCATTGTGGAATAATTAAATCAAGCTACTTGGCATGTCAGTAACATCACATGCTTAGCATTTTTGTAGTGAAAATATTTAAAATCTACTCTTTTAGCAATTTTGAAATATACAATACAGTACTTACTCACTTAACATCATTGATTGGTCCTCAGAAACTGCAACTTGGAGTGAAAAGATGTATAAAGAAACCAATTTTCCCATAGGCTAATAGCTATAAATAAGAGTTAGATTCTTACAGCATATTTTTGGTCACAAAATATCACCAAACTTCTAAATAAAGACCAAAACACTTCAAATATTAAACATTGAAATAAATATGAGCTTTGCATACATTTAAGAAAGATTAATAAAAACAAGTAAGATAATTATTTGCCCAATTATTTCATTCAGGGTTGGGGAGACTGGAGTCTGTGCTGGAAGCTCAGGGCTCAAGCTGGGCAACAGCCCTGGACAGGATGCCATCCCACTGCAGGATGGCTCACACATGCCCACAGCCACTCAGCCTGGGACCATTTGGACACAGCAATTAACCTTACCTGCATGTCTTTGTGGGGAGGAAACCAGAGTGCTTAGAAAAACCCATGCAGACAGACACAGAGCAAACATGCAAACCTCACAAAGATATTGTTTCTTCTGTCACCTGTGCTTTTGGGTCATATTCAAGAAATCATTAACCAAATAAAAGTCGTGGAGCTTTTCCCTATGTTTTCTTTTAGTAGTTTTATAGTTTCAGGTCTTACATTTAACTCCTTAATCCATTTTGATTTTTGCATATGGTGTGAGATAAGCTTCTGGTTTCATTCTCCCACATGTGGATATCCAGTTCTCTGAACACCATATATGGAAGAGACTGTCATTTCCTCATGATATGTTCCTGGCACTTTTGTTGAAATCAATTGACCATAGATCTGTGGGTTTATTTCTGGCTTTTTATTCTGTTCCATTGGCCAATGTACCTGTGTTTATGCTTGTGCCTTGCTGTTTTGATTATTATAGCTTTATAATATGTTTTGAAATCAGGTAGTGTGATGCCTCCATCTTTGCTTTTTATGCTCAAGATAGTTTGGATATTCAGAGTGTTTTATGGTTCCATATACATATTGCACCACTGCACTCAAGCCTGGGTGGTA", return_trusted_branch_array[0]);


  CU_ASSERT(return_variant_start_coords_array[0]==42662); 
  
  hash_table_free(&hash_table);
  fclose(chrom_fptr);
  
  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);

}




void test_get_covg_of_nodes_in_one_but_not_other_of_two_arrays()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // First set up the hash/graph
  int kmer_size = 5;
  int number_of_bits = 8;
  int bucket_size = 8;
  int max_retries = 10;

  BinaryKmer tmp_kmer1, tmp_kmer2;

  dBGraph * hash_table = hash_table_new(number_of_bits, bucket_size,
                                        max_retries, kmer_size);

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/one_person_two_reads.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  //>read1
  //AAAACGAAAAAATTCGAG
  //>read2
  //TCGAGAAAACGAGA

  //We will make two arrays of nodes, one corresponding to each read, and compare these two arrays. There is no reason to be associated with a read, this is just to make the test easier.
  // Where the nodes come from is irrelevant

  dBNode* n1 = hash_table_find(element_get_key(seq_to_binary_kmer("AAAAC", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table); //seen in both reads
  dBNode* n2 = hash_table_find(element_get_key(seq_to_binary_kmer("AAACG", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table); //covg 2
  dBNode* n3 = hash_table_find(element_get_key(seq_to_binary_kmer("AACGA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table); //covg 2
  dBNode* n4 = hash_table_find(element_get_key(seq_to_binary_kmer("ACGAA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table); //covg 1
  dBNode* n5 = hash_table_find(element_get_key(seq_to_binary_kmer("CGAAA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table); //covg 1
  dBNode* n6 = hash_table_find(element_get_key(seq_to_binary_kmer("GAAAA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table); //seen in both reads
  dBNode* n7 = hash_table_find(element_get_key(seq_to_binary_kmer("AAAAA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table); //covg 2
  dBNode* n8 = hash_table_find(element_get_key(seq_to_binary_kmer("AAAAA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table); //covg 2
  dBNode* n9 = hash_table_find(element_get_key(seq_to_binary_kmer("AAAAT", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table); //covg1 
  dBNode* n10 = hash_table_find(element_get_key(seq_to_binary_kmer("AAATT", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table); //covg 1
  dBNode* n11 = hash_table_find(element_get_key(seq_to_binary_kmer("AATTC", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table); //covg 1
  dBNode* n12 = hash_table_find(element_get_key(seq_to_binary_kmer("ATTCG", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table); //covg 1
  dBNode* n13 = hash_table_find(element_get_key(seq_to_binary_kmer("TTCGA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table); //covg1
  dBNode* n14 = hash_table_find(element_get_key(seq_to_binary_kmer("TCGAG", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table); // seen in both reads

  dBNode* e1 = hash_table_find(element_get_key(seq_to_binary_kmer("TCGAG", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table); //seen in both reads
  dBNode* e2 = hash_table_find(element_get_key(seq_to_binary_kmer("CGAGA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table); //covg 2
  dBNode* e3 = hash_table_find(element_get_key(seq_to_binary_kmer("GAGAA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table); //covg 1
  dBNode* e4 = hash_table_find(element_get_key(seq_to_binary_kmer("AGAAA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table); //covg 1
  dBNode* e5 = hash_table_find(element_get_key(seq_to_binary_kmer("GAAAA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table); //seen in both reads
  dBNode* e6 = hash_table_find(element_get_key(seq_to_binary_kmer("AAAAC", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2), hash_table); //seen in both reads
  //note we haven't gone all the way to the end of the second read



  CU_ASSERT(n1 != NULL);
  CU_ASSERT(n2 != NULL);
  CU_ASSERT(n3 != NULL);
  CU_ASSERT(n4 != NULL);
  CU_ASSERT(n5 != NULL);
  CU_ASSERT(n6 != NULL);
  CU_ASSERT(n7 != NULL);
  CU_ASSERT(n8 != NULL);
  CU_ASSERT(n9 != NULL);
  CU_ASSERT(n10 != NULL);
  CU_ASSERT(n11 != NULL);
  CU_ASSERT(n12 != NULL);
  CU_ASSERT(n13 != NULL);
  CU_ASSERT(n14 != NULL);

  CU_ASSERT(e1 != NULL);
  CU_ASSERT(e2 != NULL);
  CU_ASSERT(e3 != NULL);
  CU_ASSERT(e4 != NULL);
  CU_ASSERT(e5 != NULL);
  CU_ASSERT(e6 != NULL);


  

  dBNode* array1[] = {n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14};
  dBNode* array2[] = {e1,e2,e3,e4,e5,e6};


  int num_nodes_in_array1_not_array2 = 0;
  int num_nodes_in_array2_not_array1 = 0;

  Covg coverages_of_nodes_in_array1_not_array2[14];
  Covg coverages_of_nodes_in_array2_not_array1[6];
  
  Covg* ptrs_to_coverages_of_nodes_in_array1_not_array2[14];
  Covg* ptrs_to_coverages_of_nodes_in_array2_not_array1[6];

  dBNode* working_array1[20];
  dBNode* working_array2[20];

  int i;
  for (i=0; i<14; i++)
    {
      coverages_of_nodes_in_array1_not_array2[i]=0;
      ptrs_to_coverages_of_nodes_in_array1_not_array2[i]=&coverages_of_nodes_in_array1_not_array2[i];
      
    }
  for (i=0; i<6; i++)
    {
      coverages_of_nodes_in_array2_not_array1[i]=0;
      ptrs_to_coverages_of_nodes_in_array2_not_array1[i]=&coverages_of_nodes_in_array2_not_array1[i];

    }

  get_covg_of_nodes_in_one_but_not_other_of_two_arrays(array1, array2, 14,6,&num_nodes_in_array1_not_array2 , &num_nodes_in_array2_not_array1, 
						       ptrs_to_coverages_of_nodes_in_array1_not_array2, ptrs_to_coverages_of_nodes_in_array2_not_array1, 
						       working_array1, working_array2,0);

  CU_ASSERT(num_nodes_in_array1_not_array2==10 );
  CU_ASSERT(num_nodes_in_array2_not_array1==3 );

  //we also get the coverages of these nodes, in some random order (depends on the mem addresses of the nodes)
  //Note that in advance, we don't know how many such nodes there will be, if any
  qsort(coverages_of_nodes_in_array1_not_array2, num_nodes_in_array1_not_array2, sizeof(int), int_cmp); 
  qsort(coverages_of_nodes_in_array2_not_array1, num_nodes_in_array2_not_array1, sizeof(int), int_cmp); 

  CU_ASSERT(coverages_of_nodes_in_array1_not_array2[0]==1);
  CU_ASSERT(coverages_of_nodes_in_array1_not_array2[1]==1);
  CU_ASSERT(coverages_of_nodes_in_array1_not_array2[2]==1);
  CU_ASSERT(coverages_of_nodes_in_array1_not_array2[3]==1);
  CU_ASSERT(coverages_of_nodes_in_array1_not_array2[4]==1);
  CU_ASSERT(coverages_of_nodes_in_array1_not_array2[5]==1);
  CU_ASSERT(coverages_of_nodes_in_array1_not_array2[6]==1);
  CU_ASSERT(coverages_of_nodes_in_array1_not_array2[7]==2);
  CU_ASSERT(coverages_of_nodes_in_array1_not_array2[8]==2);
  CU_ASSERT(coverages_of_nodes_in_array1_not_array2[9]==2);

  CU_ASSERT(coverages_of_nodes_in_array2_not_array1[0]==1);
  CU_ASSERT(coverages_of_nodes_in_array2_not_array1[1]==1);
  CU_ASSERT(coverages_of_nodes_in_array2_not_array1[2]==2);

  hash_table_free(&hash_table);
}




void test_apply_to_all_nodes_in_path_defined_by_fasta()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // First set up the hash/graph
  int kmer_size = 5;
  int number_of_bits = 10;
  int bucket_size = 8;
  int max_retries = 10;

  dBGraph *db_graph = hash_table_new(number_of_bits, bucket_size,
                                     max_retries, kmer_size);

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/variations/two_people_short_seq_with_one_base_difference.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  // Load two people, one of whom contains this:

  // >read 1
  // AATAGACGCCCACACCTGATAGAAGCCACACTGTACTTGTANNNNNNNNNNNN
  // NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN



  char** results_array= (char**) malloc(sizeof(char*)*1000);
  if (results_array==NULL)
    {
      die("Cannot alloc results_array in test of apply fn - exit.");
    }
  int i;
  for(i=0; i< 1000; i++)
    {
      results_array[i]=(char*) malloc(sizeof(char)*6);
      if (results_array[i]==NULL)
	{
	  die("Cannot alloc results array - quit");
	}
      results_array[i][0]='\0';
    }

  
  char tmpseq[db_graph->kmer_size+1];
  int count=0;





  void test_func(dBNode* node)
    {
      if ( (node!=NULL) ) 
	{
	  BinaryKmer k;
	  binary_kmer_assignment_operator(k, *element_get_kmer(node));
	  binary_kmer_to_seq(&k, db_graph->kmer_size, tmpseq);
	  //printf("%s\n", tmpseq);
	  strcat(results_array[count],tmpseq);
	  count++;
	}
      else
	{
	  //printf("Hit a null kmer at count %d\n", count);
	}

    }

  FILE* fasta_fptr = fopen("../data/test/pop_graph/variations/first_person_short_seq.fa", "r");
  apply_to_all_nodes_in_path_defined_by_fasta(&test_func, fasta_fptr, 10, db_graph);
  fclose(fasta_fptr);

  CU_ASSERT_STRING_EQUAL(results_array[0], "AATAG");
  CU_ASSERT_STRING_EQUAL(results_array[1], "ATAGA");
  CU_ASSERT_STRING_EQUAL(results_array[2], "GTCTA");
  CU_ASSERT_STRING_EQUAL(results_array[3], "AGACG");
  CU_ASSERT_STRING_EQUAL(results_array[4], "GACGC");
  CU_ASSERT_STRING_EQUAL(results_array[5], "ACGCC");
  CU_ASSERT_STRING_EQUAL(results_array[6], "CGCCC");
  CU_ASSERT_STRING_EQUAL(results_array[7], "GCCCA");
  CU_ASSERT_STRING_EQUAL(results_array[8], "CCCAC");
  CU_ASSERT_STRING_EQUAL(results_array[9], "CCACA");
  CU_ASSERT_STRING_EQUAL(results_array[10],"CACAC");
  CU_ASSERT_STRING_EQUAL(results_array[11],"ACACC");
  CU_ASSERT_STRING_EQUAL(results_array[12], "AGGTG");
  CU_ASSERT_STRING_EQUAL(results_array[13], "ACCTG");
  CU_ASSERT_STRING_EQUAL(results_array[36], "TACAA");


  
  // Now load some more fastas, including one which has some N's in the middle
  // Make sure we can get the right kmers when we follow that path

  files_loaded = 0;
  bad_reads = 0;
  dup_reads = 0;
  seq_loaded = 0;
  seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/variations/two_people_both_alu_Ns_alu_with_one_base_difference.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  //cleanup results

  for(i=0; i< 1000; i++)
    {
      results_array[i][0]='\0';
    }
  count =0;

  fasta_fptr = fopen("../data/test/pop_graph/variations/one_person_aluNsalu.fa", "r");

  /*
    >7SLRNA#SINE/Alu 
    GTTCAGAGGCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG
    GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCC
    GATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG
    CGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAA
    ACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGT
    GAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTT
    AAAAAAAAAAAAAAAAAAAAGTCAGCCGTAGNNNNNNNNNNNNNNNNNNN
    GTTCAGAGGCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG
    GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCC
    GATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG
    CGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAA
    ACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGT
    GAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTT
    AAAAAAAAAAAAAAAAAAAAGTCAGCCGTAG
  */

  
  apply_to_all_nodes_in_path_defined_by_fasta(&test_func, fasta_fptr, 10, db_graph);
  fclose(fasta_fptr);

  CU_ASSERT_STRING_EQUAL(results_array[0], "GTTCA");

  CU_ASSERT_STRING_EQUAL(results_array[29], "CTGTA");

  CU_ASSERT_STRING_EQUAL(results_array[334], "CGTAG");
  CU_ASSERT_STRING_EQUAL(results_array[335], "GTTCA");


  for(i=0; i< 1000; i++)
    {
      free(results_array[i]);
    }
  free(results_array);
  hash_table_free(&db_graph);
}

void test_does_this_path_exist_in_this_colour()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  if (NUMBER_OF_COLOURS < 3)
  {
    warn("Test is not configured for NUMBER_OF_COLOURS < 3");
    return;
  }

  //first set up the hash/graph                                                                                                                                                                            
  int kmer_size = 7;
  int number_of_bits = 10;
  int bucket_size = 8;
  int max_retries = 10;
  
  dBGraph * db_graph = hash_table_new(number_of_bits, bucket_size,
                                      max_retries, kmer_size);

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/three_colours.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 1, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  // Annoyingly, I have called these file colour1 and colour2, but in the
  // hash table they are colours 0 and 1. Sorry for this. From here on, I
  // use colour0 colour1 - the original filenames irrelevant
  /* in colour 0 we have

>read1
AAGCATAGACACAGTGAGAGAC
>read2
CCCCCTTTTCCCCCTTTTCCCC

and in colour1 we have

>read1
AAGCATAGACACAGTGAGA
>read3
TTTGTGTTTTGTGTG

and in colour2  - which is NOT loaded into th hgraph, we have
>read4 bridges read2 and read3 - so if you mix colours 0 and 1, you'll get
one long contig joining read2 and read3

TTCCCCTTTGTGTTTTGTGTG

  */
  
  //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------
  int max_read_length = 100;
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence");
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);
  
  //We are going to load all the bases into a single sliding window 
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
  {
    die("Failed to malloc kmer sliding window in test_does_this_path_exist_in_this_colour. Exit.");
  }
  
  
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-db_graph->kmer_size-1));
  if (kmer_window->kmer==NULL)
  {
    die("Failed to malloc kmer_window->kmer in test_does_this_path_exist_in_this_colour. Exit.");
  }
  kmer_window->nkmers=0;
  
  //create file reader
  int file_reader(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
long long ret;
int offset = 0;
if (new_entry == false){
die("new_entry must be true in hsi test function");
}
ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);

return ret;
  }
  
  
  dBNode* array_nodes[50];//in fact there are 43 17-mers in the first line of the fasta
  Orientation array_or[50];
  
  //end of intialisation 
  
  
  //so let's check read1 - should see it in colour1 and colour2, and the union of colours 1 and 2
  
  FILE* fp = fopen("../data/test/pop_graph/colour0.fa", "r");
  if (fp==NULL)
  {
    die("Cannot open ../data/test/pop_graph/colour0.fa");
  }

  boolean f_entry=true;
  int len_array = align_next_read_to_graph_and_return_node_array(fp, 50, array_nodes, array_or, false, &f_entry, file_reader, seq, kmer_window, db_graph, 0);
  CU_ASSERT(f_entry==true);
  //now the test
  
  CU_ASSERT(does_this_path_exist_in_this_colour(array_nodes, array_or, len_array, &element_get_colour0, db_graph)==true);
  CU_ASSERT(does_this_path_exist_in_this_colour(array_nodes, array_or, len_array, &element_get_colour1, db_graph)==true);
  
  
  //now read2:
  f_entry=true;
  len_array = align_next_read_to_graph_and_return_node_array(fp, 50, array_nodes, array_or, false, &f_entry, file_reader, seq, kmer_window, db_graph, 0);
  CU_ASSERT(f_entry=true);
  CU_ASSERT(does_this_path_exist_in_this_colour(array_nodes, array_or, len_array, &element_get_colour0, db_graph)==true);
  CU_ASSERT(does_this_path_exist_in_this_colour(array_nodes, array_or, len_array, &element_get_colour1, db_graph)==false);
  
  //now try the other fasta file and check read1 and 3
  fclose(fp);
  fp = fopen("../data/test/pop_graph/colour1.fa", "r");
  if (fp==NULL)
  {
    die("Cannot open ../data/test/pop_graph/colour1.fa");
  }
  
  
  //this read is in both colours
  f_entry=true;
  len_array = align_next_read_to_graph_and_return_node_array(fp, 50, array_nodes, array_or, false, &f_entry, file_reader, seq, kmer_window, db_graph, 0);
  CU_ASSERT(f_entry=true);
  CU_ASSERT(does_this_path_exist_in_this_colour(array_nodes, array_or, len_array, &element_get_colour0, db_graph)==true);
  CU_ASSERT(does_this_path_exist_in_this_colour(array_nodes, array_or, len_array, &element_get_colour1, db_graph)==true);
  //this read is in colour1 only
  f_entry=true;
  len_array = align_next_read_to_graph_and_return_node_array(fp, 50, array_nodes, array_or, false, &f_entry, file_reader, seq, kmer_window, db_graph, 0);
  CU_ASSERT(f_entry=true);
  CU_ASSERT(does_this_path_exist_in_this_colour(array_nodes, array_or, len_array, &element_get_colour0, db_graph)==false);
  CU_ASSERT(does_this_path_exist_in_this_colour(array_nodes, array_or, len_array, &element_get_colour1, db_graph)==true);
  
  //now for the final test, take a read which is there in the union of two colours, but not in either
  fclose(fp);
  fp = fopen("../data/test/pop_graph/colour2.fa", "r");
  if (fp==NULL)
  {
    die("Cannot open ../data/test/pop_graph/colour2.fa");
  }

  f_entry=true;
  len_array = align_next_read_to_graph_and_return_node_array(fp, 50, array_nodes, array_or, false, &f_entry,file_reader, seq, kmer_window, db_graph, 0);
  CU_ASSERT(f_entry=true);
  CU_ASSERT(does_this_path_exist_in_this_colour(array_nodes, array_or, len_array, &element_get_colour0, db_graph)==false);
  CU_ASSERT(does_this_path_exist_in_this_colour(array_nodes, array_or, len_array, &element_get_colour1, db_graph)==false);
  CU_ASSERT(does_this_path_exist_in_this_colour(array_nodes, array_or, len_array, &element_get_colour_union_of_all_colours, db_graph)==true);
  
  
  
  free(kmer_window->kmer);
  free(kmer_window);
  free_sequence(&seq);
  hash_table_free(&db_graph);
}



//this is a regression test to check that there is no problem with either overflowing ints,
//or dumping causing segfaults. I'm not checking content of the covg distribution.
void test_dump_covg_distribution()
{

  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  //first set up the hash/graph
  int kmer_size = 5;
  int number_of_bits = 3;
  int bucket_size = 2;
  int max_retries = 10;

  dBGraph * hash_table = hash_table_new(number_of_bits, bucket_size,
                                        max_retries, kmer_size);

  if(hash_table == NULL)
  {
    die("Unable to alloc the hash table in test_dump_covg_distribution. dead before we even started.");
  }

  // Read FASTA sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/file_for_testing_dump_covg.fa.filelist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, hash_table, 0, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  //now let's up the covg on one of the nodes
  BinaryKmer tmp_kmer1, tmp_kmer2;
  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AAAAA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,hash_table);
  CU_ASSERT(test_element1!=NULL);

  uint64_t i;
  for (i=0; i<=8000000000; i++)
    {
      db_node_update_coverage(test_element1, 0, 1);
    }
  
  //dump binary, dump covg
  db_graph_dump_single_colour_binary_of_colour0("../data/tempfiles_can_be_deleted/binary_for_cvg_distrib.ctx", 
						&db_node_condition_always_true,
						hash_table, NULL, BINVERSION);

  
  db_graph_get_covg_distribution("../data/tempfiles_can_be_deleted/covg_distribution",
				 hash_table, 
				 0, 
				 &db_node_check_status_not_pruned);
  
  hash_table_free(&hash_table);


  

}
