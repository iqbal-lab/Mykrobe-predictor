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
  test_file_reader.c
*/

// system headers
#include <stdlib.h>
#include <limits.h>

// third party headers
#include <CUnit.h>
#include <Basic.h>
#include <string_buffer.h>

// cortex_var headers
#include "file_reader.h"
#include "dB_graph_population.h"
#include "element.h"
#include "seq.h"
#include "open_hash/hash_table.h"
#include "test_file_reader.h"
#include "graph_info.h"

void test_dump_load_sv_trio_binary()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER == 1)
  {
    int kmer_size = 3;
    int number_of_bits_pre = 4;
    int number_of_bits_post = 8;
    int bucket_size = 5;
    int max_retries = 10;

    dBGraph *db_graph_pre, *db_graph_post;

    db_graph_pre = hash_table_new(number_of_bits_pre, bucket_size,
                                  max_retries, kmer_size);

    //we need the following arguments for the API but we will not use them -
    // for duplicate removal and homopolymer breaking
    int fq_quality_cutoff = 20;
    int homopolymer_cutoff = 0;
    boolean remove_duplicates_se = false;
    char ascii_fq_offset = 33;
    int into_colour = 0;

    unsigned int files_loaded = 0;
    unsigned long long bad_reads = 0, dup_reads = 0;
    unsigned long long seq_read = 0, seq_loaded = 0;

    load_se_filelist_into_graph_colour(
      "../data/test/graph/test_dB_graph.falist",
      fq_quality_cutoff, homopolymer_cutoff,
      remove_duplicates_se, ascii_fq_offset,
      into_colour, db_graph_pre, 0,
      &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
      NULL, 0, &subsample_null);

    CU_ASSERT(seq_read == 16);
    CU_ASSERT(seq_loaded == 16);

    //hash_table_traverse(&print_node_binary,db_graph_pre);
    GraphInfo* ginfo = graph_info_alloc_and_init();

    graph_info_set_seq(ginfo, 0, seq_read);
    graph_info_set_mean_readlen(ginfo, 0, 5);
    graph_info_set_specific_colour_to_cleaned_against_pool(ginfo, 0, "zammo.ctx", 401);
    graph_info_set_seq_err(ginfo, 0, 0.456);
    graph_info_set_remv_low_cov_sups(ginfo, 0, 178);
     
    db_graph_dump_binary("../data/tempfiles_can_be_deleted/dump_cortex_var_graph.ctx", 
  		&db_node_condition_always_true, db_graph_pre, ginfo, BINVERSION);


    hash_table_free(&db_graph_pre);
    CU_ASSERT(db_graph_pre==NULL);

    db_graph_post = hash_table_new(number_of_bits_post, bucket_size, 10, kmer_size);
    int num_cols_in_binary = -1;

    graph_info_initialise(ginfo);

    long long seq_length_post
      = load_multicolour_binary_from_filename_into_graph(
          "../data/tempfiles_can_be_deleted/dump_cortex_var_graph.ctx",
          db_graph_post, ginfo, &num_cols_in_binary);

    CU_ASSERT(num_cols_in_binary==NUMBER_OF_COLOURS);
    CU_ASSERT(ginfo->mean_read_length[0]==5);
    CU_ASSERT(ginfo->total_sequence[0]==seq_read);

    CU_ASSERT(ginfo->seq_err[0]-0.456<0.0001);
    CU_ASSERT(ginfo->cleaning[0]->tip_clipping==false);
    CU_ASSERT(ginfo->cleaning[0]->remv_low_cov_sups==true);
    CU_ASSERT(ginfo->cleaning[0]->remv_low_cov_sups_thresh==178);
    CU_ASSERT(ginfo->cleaning[0]->remv_low_cov_nodes==false);
    CU_ASSERT(ginfo->cleaning[0]->remv_low_cov_nodes_thresh==-1);
    CU_ASSERT(strcmp(ginfo->cleaning[0]->name_of_graph_against_which_was_cleaned, "zammo.ctx (colour 401)")==0);
    CU_ASSERT(ginfo->cleaning[0]->len_name_of_graph_against_which_was_cleaned==strlen("zammo.ctx (colour 401)"));
    // load_multicolour_binary_data_from_filename_into_graph returns total number
    // of unique kmers loaded, times kmer_length
    CU_ASSERT_EQUAL(seq_length_post,15);
    CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph_post),5);

    BinaryKmer tmp_kmer1, tmp_kmer2;

    // all the kmers and their reverse complements from the reads
    dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AAA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("TTT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element3 = hash_table_find(element_get_key(seq_to_binary_kmer("GGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element4 = hash_table_find(element_get_key(seq_to_binary_kmer("GCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element5 = hash_table_find(element_get_key(seq_to_binary_kmer("GCT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element6 = hash_table_find(element_get_key(seq_to_binary_kmer("AGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element7 = hash_table_find(element_get_key(seq_to_binary_kmer("TAG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element8 = hash_table_find(element_get_key(seq_to_binary_kmer("CTA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element9 = hash_table_find(element_get_key(seq_to_binary_kmer("AGG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element10 = hash_table_find(element_get_key(seq_to_binary_kmer("CCT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);


    // kmers that should not be in the graph
    dBNode* test_element11 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element12 = hash_table_find(element_get_key(seq_to_binary_kmer("CCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element13 = hash_table_find(element_get_key(seq_to_binary_kmer("TAT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element14 = hash_table_find(element_get_key(seq_to_binary_kmer("ATA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element15 = hash_table_find(element_get_key(seq_to_binary_kmer("TAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element16 = hash_table_find(element_get_key(seq_to_binary_kmer("ATG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element17 = hash_table_find(element_get_key(seq_to_binary_kmer("TTG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element18 = hash_table_find(element_get_key(seq_to_binary_kmer("AAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element19 = hash_table_find(element_get_key(seq_to_binary_kmer("TGA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element20 = hash_table_find(element_get_key(seq_to_binary_kmer("TCA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);


    CU_ASSERT(test_element1 != NULL);
    CU_ASSERT(test_element2 != NULL);
    CU_ASSERT(test_element1 == test_element2);
    // Checking that we count KMER coverage. ie this kmer occurs many times in the same read
    CU_ASSERT(db_node_get_coverage(test_element1,0)==6);

    CU_ASSERT(test_element3 != NULL);
    CU_ASSERT(test_element4 != NULL);
    CU_ASSERT(test_element3 == test_element4);
    CU_ASSERT(db_node_get_coverage(test_element3,0)==1);

    CU_ASSERT(test_element5 != NULL);
    CU_ASSERT(test_element6 != NULL);
    CU_ASSERT(test_element5 == test_element6);
    CU_ASSERT(db_node_get_coverage(test_element5,0)==1);

    CU_ASSERT(test_element7 != NULL);
    CU_ASSERT(test_element8 != NULL);
    CU_ASSERT(test_element7 == test_element8);
    CU_ASSERT(db_node_get_coverage(test_element7,0)==1);

    CU_ASSERT(test_element9 != NULL);
    CU_ASSERT(test_element10 != NULL);
    CU_ASSERT(test_element9 == test_element10);
    CU_ASSERT(db_node_get_coverage(test_element9,0)==1);

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


    // Check arrows

    Nucleotide base;

    // AAA -A-> AAA
    CU_ASSERT(db_node_has_precisely_one_edge(test_element1, forward,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Adenine);

    CU_ASSERT(db_node_has_precisely_one_edge(test_element1, reverse,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Thymine);

    //GGC -T-> GCT
    CU_ASSERT(db_node_has_precisely_one_edge(test_element3, reverse,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Thymine);


    //TAG -G-> AGG
    CU_ASSERT(db_node_has_precisely_one_edge(test_element7, reverse,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Guanine);

    //AGC -C-> GCC
    CU_ASSERT(db_node_has_precisely_one_edge(test_element6, forward,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Cytosine);

    //CCT -A-> CTA
    CU_ASSERT(db_node_has_precisely_one_edge(test_element10, reverse,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Adenine);

    //add one extra arrows by had -- this breaks the graph - it is only to check that the arrows get set correctly

    add_edges(test_element1,0,0x02);
    CU_ASSERT(db_node_has_precisely_one_edge(test_element1, forward,&base, 0)==false);


    CU_ASSERT(db_node_has_precisely_one_edge(test_element1, reverse,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Thymine);

    CU_ASSERT(db_node_edge_exist(test_element1, Adenine, forward, 0)==true);
    CU_ASSERT(db_node_edge_exist(test_element1, Cytosine, forward, 0)==true);
    CU_ASSERT(db_node_edge_exist(test_element1, Guanine, forward, 0)==false);
    CU_ASSERT(db_node_edge_exist(test_element1, Thymine, forward, 0)==false);



    add_edges(test_element1,0,0x20);
    CU_ASSERT(db_node_has_precisely_one_edge(test_element1, reverse,&base, 0)==false);

    CU_ASSERT(db_node_edge_exist(test_element1, Adenine, forward, 0)==true);
    CU_ASSERT(db_node_edge_exist(test_element1, Cytosine, forward, 0)==true);
    CU_ASSERT(db_node_edge_exist(test_element1, Guanine, forward, 0)==false);
    CU_ASSERT(db_node_edge_exist(test_element1, Thymine, forward, 0)==false);
    CU_ASSERT(db_node_edge_exist(test_element1, Adenine, reverse, 0)==false);
    CU_ASSERT(db_node_edge_exist(test_element1, Cytosine, reverse, 0)==true);
    CU_ASSERT(db_node_edge_exist(test_element1, Guanine, reverse, 0)==false);
    CU_ASSERT(db_node_edge_exist(test_element1, Thymine, reverse, 0)==true);


    hash_table_free(&db_graph_post);
    graph_info_free(ginfo);
    CU_ASSERT(db_graph_post == NULL);
  }

  // Now try the same thing with big kmers

  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER == 2)
  {
    int kmer_size = 33;
    int number_of_bits_pre = 10;
    int number_of_bits_post = 12;
    int bucket_size = 10;
    int max_retries = 10;

    dBGraph *db_graph_pre, *db_graph_post;

    db_graph_pre = hash_table_new(number_of_bits_pre, bucket_size,
                                  max_retries, kmer_size);

    int fq_quality_cutoff = 20;
    int homopolymer_cutoff = 0;
    boolean remove_duplicates_se = false;
    char ascii_fq_offset = 33;
    int into_colour = 0;

    unsigned int files_loaded = 0;
    unsigned long long bad_reads = 0, dup_reads = 0;
    unsigned long long seq_read = 0, seq_loaded = 0;

    load_se_filelist_into_graph_colour(
      "../data/test/graph/person2.falist",
      fq_quality_cutoff, homopolymer_cutoff,
      remove_duplicates_se, ascii_fq_offset,
      into_colour, db_graph_pre, 0,
      &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
      NULL, 0, &subsample_null);

    // > 6 unique 33-mers
    // TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACC
    // > 13 unique 33-mers
    // ACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAAC
    // > 12 unique 33-mers
    // GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT

    GraphInfo* ginfo = graph_info_alloc_and_init();

    graph_info_initialise(ginfo);
    // Not bothering to set mean read length this time
    graph_info_set_seq(ginfo, 0, seq_read);
    db_graph_dump_binary("../data/tempfiles_can_be_deleted/dump_cortex_var_graph_2.ctx",
      &db_node_condition_always_true, db_graph_pre, ginfo, BINVERSION);
    //hash_table_traverse(&print_node_binary,db_graph_pre);

    hash_table_free(&db_graph_pre);
    CU_ASSERT(db_graph_pre==NULL);

    db_graph_post = hash_table_new(number_of_bits_post,bucket_size,max_retries,kmer_size);

    int num_cols_in_binary=-1;

    graph_info_initialise(ginfo);
    
    load_multicolour_binary_from_filename_into_graph(
          "../data/tempfiles_can_be_deleted/dump_cortex_var_graph_2.ctx",
          db_graph_post, ginfo, &num_cols_in_binary);

    CU_ASSERT(ginfo->mean_read_length[0]==0);//because we did not set it
    CU_ASSERT(num_cols_in_binary==NUMBER_OF_COLOURS);
    CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph_post),31);

    BinaryKmer tmp_kmer1, tmp_kmer2;

    //some the kmers and their reverse complements from the reads
    dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("TAACCCTAACCCTAACCCTAACCCTAACCCTAA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element3 = hash_table_find(element_get_key(seq_to_binary_kmer("AACCCTAACCCTAACCCTAACCCTAACCCTAAC", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element4 = hash_table_find(element_get_key(seq_to_binary_kmer("GTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTT", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element5 = hash_table_find(element_get_key(seq_to_binary_kmer("ACCCTAACCCTAACCCTAACCCTAACCCTAACC", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element6 = hash_table_find(element_get_key(seq_to_binary_kmer("GGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGT", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element7 = hash_table_find(element_get_key(seq_to_binary_kmer("CCCTAACCCTAACCCTAACCCTAACCCTAACCC", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element8 = hash_table_find(element_get_key(seq_to_binary_kmer("GGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element9 = hash_table_find(element_get_key(seq_to_binary_kmer("CCTAACCCTAACCCTAACCCTAACCCTAACCCT", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element10 = hash_table_find(element_get_key(seq_to_binary_kmer("AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGG", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element11 = hash_table_find(element_get_key(seq_to_binary_kmer("CTAACCCTAACCCTAACCCTAACCCTAACCCTA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element12 = hash_table_find(element_get_key(seq_to_binary_kmer("TAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAG", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);

      //from read 2:
    dBNode* test_element13 = hash_table_find(element_get_key(seq_to_binary_kmer("ACCCTAACCCTAACCCTAACCCCTAACCCTAAC", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element14 = hash_table_find(element_get_key(seq_to_binary_kmer("GTTAGGGTTAGGGGTTAGGGTTAGGGTTAGGGT", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element15 = hash_table_find(element_get_key(seq_to_binary_kmer("CCCTAACCCTAACCCTAACCCCTAACCCTAACC", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element16 = hash_table_find(element_get_key(seq_to_binary_kmer("GGTTAGGGTTAGGGGTTAGGGTTAGGGTTAGGG", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);

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

    CU_ASSERT(test_element11 != NULL);
    CU_ASSERT(test_element12 != NULL);
    CU_ASSERT(test_element11 == test_element12);


    CU_ASSERT(test_element13 != NULL);
    CU_ASSERT(test_element14 != NULL);
    CU_ASSERT(test_element13 == test_element14);

    CU_ASSERT(test_element15 != NULL);
    CU_ASSERT(test_element16 != NULL);
    CU_ASSERT(test_element15 == test_element16);

    //now check arrows
    Nucleotide base;

    //Note we know that read1 forms a supernode that is a loop.

    // TAACCCTAACCCTAACCCTAACCCTAACCCTAA ----- C ----> AACCCTAACCCTAACCCTAACCCTAACCCTAAC
    CU_ASSERT(db_node_has_precisely_one_edge(test_element1, forward,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Cytosine);

    CU_ASSERT(db_node_has_precisely_one_edge(test_element1, reverse,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Guanine);

    // AACCCTAACCCTAACCCTAACCCTAACCCTAAC ----C  ----> ACCCTAACCCTAACCCTAACCCTAACCCTAACC
    CU_ASSERT(db_node_has_precisely_one_edge(test_element3, forward,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Cytosine);

    CU_ASSERT(db_node_has_precisely_one_edge(test_element3, reverse,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Adenine);

    // ACCCTAACCCTAACCCTAACCCTAACCCTAACC ---C -----> CCCTAACCCTAACCCTAACCCTAACCCTAACCC
    CU_ASSERT(db_node_has_precisely_one_edge(test_element5, forward,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Cytosine);

    CU_ASSERT(db_node_has_precisely_one_edge(test_element5, reverse,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Thymine);

    //OK. Looks good.

    hash_table_free(&db_graph_post);
    graph_info_free(ginfo);
  }


  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER == 1)
  {
    // Finally a test case which found a bug in binary read/write which no
    // other test case found
    int kmer_size = 17;
    int number_of_bits_pre = 10;
    int number_of_bits_post = 10;
    int bucket_size = 30;
    int max_retries = 10;

    dBGraph *db_graph_pre, *db_graph_post;

    db_graph_pre = hash_table_new(number_of_bits_pre, bucket_size,
                                  max_retries, kmer_size);

    int fq_quality_cutoff = 20;
    int homopolymer_cutoff = 0;
    boolean remove_duplicates_se = false;
    char ascii_fq_offset = 33;
    int into_colour = 0;

    unsigned int files_loaded = 0;
    unsigned long long bad_reads = 0, dup_reads = 0;
    unsigned long long seq_read = 0, seq_loaded = 0;

    load_se_filelist_into_graph_colour(
      "../data/test/graph/person3.falist",
      fq_quality_cutoff, homopolymer_cutoff, remove_duplicates_se, ascii_fq_offset,
      into_colour, db_graph_pre, 0,
      &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
      NULL, 0, &subsample_null);

    // >read1 overlaps human chrom 1
    // TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACC
    // > read 2 overlaps human chrom 1
    // ACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAAC
    // > read 3 does not
    // GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
    // > read 3 does not
    // GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
    // > read 3 does not
    // GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
    // > read 4 does not, but has too low coverage
    // TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

    GraphInfo* ginfo = graph_info_alloc_and_init();

    graph_info_initialise(ginfo);
    graph_info_set_seq(ginfo, 0, seq_read);
    db_graph_dump_binary("../data/tempfiles_can_be_deleted/dump_cortex_var_graph_3.ctx", &db_node_condition_always_true, db_graph_pre, ginfo, BINVERSION);
    hash_table_free(&db_graph_pre);
    CU_ASSERT(db_graph_pre==NULL);

    db_graph_post = hash_table_new(number_of_bits_post,bucket_size,max_retries,kmer_size);

    int num_cols_in_binary = -1;

    graph_info_initialise(ginfo);
    load_multicolour_binary_from_filename_into_graph(
      "../data/tempfiles_can_be_deleted/dump_cortex_var_graph_3.ctx",
      db_graph_post, ginfo,&num_cols_in_binary);

    CU_ASSERT(num_cols_in_binary==NUMBER_OF_COLOURS);
    CU_ASSERT(ginfo->total_sequence[0]==seq_read);

    BinaryKmer tmp_kmer1, tmp_kmer2;

    //now try to traverse a supernode. This is effectively a regressiontest for a bug in graph/element.c: print_binary/read_binary
    dBNode *test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("TAACCCTAACCCTAACC", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);


    // this node is in the middle of this supernode: CCCTAACCCTAACCCTAACCC. You can't extend forward because there is one copy in the reads with a T afterwards
    // and one with a C afterwards. And you can' extend the other way because the first kmer CCCTAACCCTAACCCTA has two arrows in. ie is preceded by two different bases (A and T) in
    // different reads

    CU_ASSERT(test_element1 != NULL);

    dBNode * path_nodes[100];
    Orientation path_orientations[100];
    Nucleotide path_labels[100];
    char path_string[100];
    int limit=100;
    double avg_covg;
    Covg min_covg, max_covg;
    boolean is_cycle;

    int len = db_graph_supernode_for_specific_person_or_pop(
                test_element1, limit, &db_node_action_do_nothing,
                path_nodes, path_orientations, path_labels, path_string,
                &avg_covg, &min_covg, &max_covg, &is_cycle, db_graph_post, 0);

    CU_ASSERT(len==4);
    CU_ASSERT(is_cycle==false);

    CU_ASSERT( (!strcmp(path_string, "AGGG") ) || (!strcmp(path_string, "ACCC"))) ;

    hash_table_free(&db_graph_post);
    graph_info_free(ginfo);
  }
}



void test_load_singlecolour_binary()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not written for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // Build graph
  int kmer_size = 5;
  int number_of_bits = 7;
  int bucket_size = 5;
  int max_retries = 10;

  dBGraph* db_graph_pre = hash_table_new(number_of_bits, bucket_size,
                                         max_retries, kmer_size);

  /*
  >read
  AACGTTCC
  >read
  AACGTTCC
  >zam
  GTTCA
  >read1 - will repeat 300 times
  AAAAAAA
  >read1
  AAAAAAA

  note this contains 5 unique kmers.read 1 seems to contain 4 unique kmers
  but only contains 3, as AACGTT is just one kmer looped back on itself

  We need the following arguments for the API but we will not use them -
  for duplicate removal and homopolymer breaking
  */

  int fq_quality_cutoff = 20;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_read = 0, seq_loaded = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/graph/test_dumping_by_graph_and_reload_by_sv_trio.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph_pre, 0,
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  GraphInfo* ginfo = graph_info_alloc_and_init();

  graph_info_set_seq(ginfo, 0, seq_read);
  //not setting mean read len, so will not test it
  db_graph_dump_single_colour_binary_of_colour0(
						"../data/tempfiles_can_be_deleted/dump_single_colour_cortex_var_graph.ctx",
						&db_node_condition_always_true, db_graph_pre, ginfo, BINVERSION);

  hash_table_free(&db_graph_pre);
  CU_ASSERT(db_graph_pre == NULL);

  dBGraph* db_graph_post = hash_table_new(number_of_bits, bucket_size,
                                          10, kmer_size);

  // zero it, and see if you get your data back
  graph_info_initialise(ginfo);
  load_single_colour_binary_data_from_filename_into_graph(
							  "../data/tempfiles_can_be_deleted/dump_single_colour_cortex_var_graph.ctx", 
							  db_graph_post, ginfo, 
							  true,0, false,0, false);
  
  CU_ASSERT(ginfo->total_sequence[0]==2809);
  CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph_post), 5);


  BinaryKmer tmp_kmer1, tmp_kmer2;

  // All the nodes and their rev complements from the graph
  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AACGT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("ACGTT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element3 = hash_table_find(element_get_key(seq_to_binary_kmer("CGTTC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element4 = hash_table_find(element_get_key(seq_to_binary_kmer("GAACG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element5 = hash_table_find(element_get_key(seq_to_binary_kmer("GTTCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element6 = hash_table_find(element_get_key(seq_to_binary_kmer("GGAAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element7 = hash_table_find(element_get_key(seq_to_binary_kmer("GTTCA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element8 = hash_table_find(element_get_key(seq_to_binary_kmer("TGAAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element9 = hash_table_find(element_get_key(seq_to_binary_kmer("AAAAA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element10 = hash_table_find(element_get_key(seq_to_binary_kmer("TTTTT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);

  // nodes that should not be in the graph
  dBNode* test_element11 = hash_table_find(element_get_key(seq_to_binary_kmer("ATATA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element12 = hash_table_find(element_get_key(seq_to_binary_kmer("TGGGG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element13 = hash_table_find(element_get_key(seq_to_binary_kmer("AATAG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element14 = hash_table_find(element_get_key(seq_to_binary_kmer("CTCTC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element15 = hash_table_find(element_get_key(seq_to_binary_kmer("GGCGG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element16 = hash_table_find(element_get_key(seq_to_binary_kmer("GGGGA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element17 = hash_table_find(element_get_key(seq_to_binary_kmer("TACTA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);


  CU_ASSERT(test_element1 != NULL);
  CU_ASSERT(test_element2 != NULL);
  CU_ASSERT(test_element1 == test_element2);
  CU_ASSERT(db_node_get_coverage(test_element1,0)==4);

  CU_ASSERT(test_element3 != NULL);
  CU_ASSERT(test_element4 != NULL);
  CU_ASSERT(test_element3 == test_element4);
  CU_ASSERT(db_node_get_coverage(test_element3,0)==2);

  CU_ASSERT(test_element5 != NULL);
  CU_ASSERT(test_element6 != NULL);
  CU_ASSERT(test_element5 == test_element6);
  CU_ASSERT(db_node_get_coverage(test_element5,0)==2);

  CU_ASSERT(test_element7 != NULL);
  CU_ASSERT(test_element8 != NULL);
  CU_ASSERT(test_element7 == test_element8);
  CU_ASSERT(db_node_get_coverage(test_element7,0)==1);


  CU_ASSERT(test_element9 != NULL);
  CU_ASSERT(test_element10 != NULL);
  CU_ASSERT(test_element9 == test_element10);

  CU_ASSERT(test_element11 == NULL);
  CU_ASSERT(test_element12 == NULL);
  CU_ASSERT(test_element13 == NULL);
  CU_ASSERT(test_element14 == NULL);
  CU_ASSERT(test_element15 == NULL);
  CU_ASSERT(test_element16 == NULL);
  CU_ASSERT(test_element17 == NULL);

  hash_table_free(&db_graph_post);
  graph_info_free(ginfo);
}



void test_load_individual_binaries_into_sv_trio()
{
  if(NUMBER_OF_COLOURS < 3)
  {
    warn("Test not configured for NUMBER_OF_COLOURS < 3\n");
    return;
  }
  else if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }
  else
  {
    //prepare a hash table

    int kmer_size = 31;
    int number_of_bits = 10;
    int bucket_size = 41;
    int max_retries = 82;

    dBGraph* db_graph = hash_table_new(number_of_bits, bucket_size,
                                       max_retries, kmer_size);

    //first dump 3 single-colour graphs

    //Person 0: created binary from this fasta
    //  >read1 taken from an Alu
    //  GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCC

    //Person 1:
    //  > read1 different line from same Alu as person1, so will be supernode of its own
    //  GATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG

    //Person2:
    //  > read1 matches first kmer of person 0, followed by a final A not G
    //  GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGA

    int fq_quality_cutoff = 0;
    int homopolymer_cutoff = 0;
    boolean remove_duplicates_se = false;
    char ascii_fq_offset = 33;
    int into_colour = 0;

    unsigned int files_loaded = 0;
    unsigned long long bad_reads = 0, dup_reads = 0;
    unsigned long long seq_read = 0;

    unsigned long long seq_loaded0 = 0, seq_loaded1 = 0, seq_loaded2 = 0;

    // Load Person0
    load_se_filelist_into_graph_colour(
      "../data/test/graph/test_dumping_by_graph_and_reload_by_sv_trio_person0.falist",
      fq_quality_cutoff, homopolymer_cutoff,
      remove_duplicates_se, ascii_fq_offset,
      into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
      &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded0,
      NULL, 0, &subsample_null);

    GraphInfo* ginfo = graph_info_alloc_and_init();

    graph_info_set_seq(ginfo, 0, seq_loaded0);
    db_graph_dump_single_colour_binary_of_colour0(
      "../data/tempfiles_can_be_deleted/test_dumping_by_graph_and_reload_by_sv_trio_person0_kmer31.ctx",
      &db_node_condition_always_true, db_graph, ginfo, BINVERSION);
    hash_table_free(&db_graph);


    // Load Person1
    db_graph = hash_table_new(number_of_bits, bucket_size,
                              max_retries, kmer_size);

    load_se_filelist_into_graph_colour(
      "../data/test/graph/test_dumping_by_graph_and_reload_by_sv_trio_person1.falist",
      fq_quality_cutoff, homopolymer_cutoff,
      remove_duplicates_se, ascii_fq_offset,
      into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
      &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded1,
      NULL, 0, &subsample_null);

    graph_info_initialise(ginfo);
    graph_info_set_seq(ginfo, 0, seq_loaded1);
    db_graph_dump_single_colour_binary_of_colour0(
      "../data/tempfiles_can_be_deleted/test_dumping_by_graph_and_reload_by_sv_trio_person1_kmer31.ctx",
      &db_node_condition_always_true, db_graph, ginfo, BINVERSION);
    hash_table_free(&db_graph);

    // Load Person2
    db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

    load_se_filelist_into_graph_colour(
      "../data/test/graph/test_dumping_by_graph_and_reload_by_sv_trio_person2.falist",
      fq_quality_cutoff, homopolymer_cutoff,
      remove_duplicates_se, ascii_fq_offset,
      into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
      &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded2,
      NULL, 0, &subsample_null);

    graph_info_initialise(ginfo);
    graph_info_set_seq(ginfo, 0, seq_loaded2);
    db_graph_dump_single_colour_binary_of_colour0(
      "../data/tempfiles_can_be_deleted/test_dumping_by_graph_and_reload_by_sv_trio_person2_kmer31.ctx",
      &db_node_condition_always_true, db_graph, ginfo, BINVERSION);
    hash_table_free(&db_graph);

    db_graph = hash_table_new(number_of_bits, bucket_size, max_retries, kmer_size);
    graph_info_initialise(ginfo);
    int first_colour = 0;
    load_population_as_binaries_from_graph(
      "../data/test/graph/trio_filelist_for_testing_loading_singlecolour_bins_into_multicol_bin.colours",
      first_colour, true, db_graph, ginfo, false, 0, false);

    CU_ASSERT(ginfo->total_sequence[0] == seq_loaded0);
    CU_ASSERT(ginfo->total_sequence[1] == seq_loaded1);
    CU_ASSERT(ginfo->total_sequence[2] == seq_loaded2);

    CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph), 41);


    //start with a kmer that should be in person0 and person 2 only
    BinaryKmer tmp_kmer1, tmp_kmer2;

    dBNode* test_element1_person0 = db_graph_find_node_restricted_to_specific_person_or_population(
      element_get_key(seq_to_binary_kmer("GTGGGAGGATCGCTTGAGTCCAGGAGTTCTG", kmer_size, &tmp_kmer1),
                      kmer_size, &tmp_kmer2), db_graph, 0);

      //look for reverse complement of GTGGGAGGATCGCTTGAGTCCAGGAGTTCTG
    dBNode* test_element2_person0 = db_graph_find_node_restricted_to_specific_person_or_population(
      element_get_key(seq_to_binary_kmer("CAGAACTCCTGGACTCAAGCGATCCTCCCAC", kmer_size, &tmp_kmer1),
                      kmer_size, &tmp_kmer2), db_graph, 0);

    dBNode* test_element1_person1 = db_graph_find_node_restricted_to_specific_person_or_population(
      element_get_key(seq_to_binary_kmer("GTGGGAGGATCGCTTGAGTCCAGGAGTTCTG", kmer_size, &tmp_kmer1),
                      kmer_size, &tmp_kmer2), db_graph, 1);

      //look for reverse complement of GTGGGAGGATCGCTTGAGTCCAGGAGTTCTG
    dBNode* test_element2_person1 = db_graph_find_node_restricted_to_specific_person_or_population(
      element_get_key(seq_to_binary_kmer("CAGAACTCCTGGACTCAAGCGATCCTCCCAC", kmer_size, &tmp_kmer1),
                      kmer_size, &tmp_kmer2), db_graph, 1);

    dBNode* test_element1_person2 = db_graph_find_node_restricted_to_specific_person_or_population(
      element_get_key(seq_to_binary_kmer("GTGGGAGGATCGCTTGAGTCCAGGAGTTCTG", kmer_size, &tmp_kmer1),
                      kmer_size, &tmp_kmer2), db_graph, 2);

      //look for reverse complement of GTGGGAGGATCGCTTGAGTCCAGGAGTTCTG
    dBNode* test_element2_person2 = db_graph_find_node_restricted_to_specific_person_or_population(
      element_get_key(seq_to_binary_kmer("CAGAACTCCTGGACTCAAGCGATCCTCCCAC", kmer_size, &tmp_kmer1),
                      kmer_size, &tmp_kmer2), db_graph, 2);

    CU_ASSERT(test_element1_person0 != NULL);
    CU_ASSERT(test_element2_person0 != NULL);
    CU_ASSERT(test_element1_person0 == test_element2_person0);
    CU_ASSERT(test_element1_person1 == NULL);
    CU_ASSERT(test_element2_person1 == NULL);
    CU_ASSERT(test_element1_person2 != NULL);
    CU_ASSERT(test_element2_person2 != NULL);
    CU_ASSERT(test_element1_person2 == test_element2_person0);
    CU_ASSERT(test_element1_person0 == test_element1_person2);

    // Now some kmers that exist in person 2 only

    dBNode* test_element3_person0
      = db_graph_find_node_restricted_to_specific_person_or_population(
          element_get_key(seq_to_binary_kmer("GATCGGGTGTCCGCACTAAGTTCGGCATCAA",
                          kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),
          db_graph, 0);

    dBNode* test_element3_person1
      = db_graph_find_node_restricted_to_specific_person_or_population(
          element_get_key(seq_to_binary_kmer("GATCGGGTGTCCGCACTAAGTTCGGCATCAA",
                          kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),
        db_graph, 1);

    dBNode* test_element3_person2
      = db_graph_find_node_restricted_to_specific_person_or_population(
          element_get_key(seq_to_binary_kmer("GATCGGGTGTCCGCACTAAGTTCGGCATCAA",
                          kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),
          db_graph, 2);


    CU_ASSERT(test_element3_person0==NULL);
    CU_ASSERT(test_element3_person1!=NULL);
    CU_ASSERT(test_element3_person2==NULL);

    if(test_element3_person1 == NULL)
    {
      warn("%s:%i: test_element3_person1 is NULL", __FILE__, __LINE__);
    }

    if(test_element1_person2 == NULL)
    {
      // Prepare a hash table
      int kmer_size = 31;
      int number_of_bits = 0;
      int bucket_size = 41;
      int max_retries = 82;
      
      dBGraph* db_graph = hash_table_new(number_of_bits, bucket_size,
                                         max_retries, kmer_size);
      

      //first dump 3 single-colour graphs
      
      //Person 0: created binary from this fasta
      //  >read1 taken from an Alu
      //  GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCC
      
      //Person 1:
      //  > read1 different line from same Alu as person1, so will be supernode of its own
      //  GATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG
      
      //Person2:
      //  > read1 matches first kmer of person 0, followed by a final A not G
      //  GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGA

      int fq_quality_cutoff = 0;
      int homopolymer_cutoff = 0;
      boolean remove_duplicates_se = false;
      char ascii_fq_offset = 33;
      int into_colour = 0;

      unsigned long long bad_reads = 0, dup_reads = 0;
      unsigned long long seq_read = 0;
      unsigned long long seq_loaded0 = 0, seq_loaded1 = 0, seq_loaded2 = 0;

      // Load Person0
      load_se_filelist_into_graph_colour(
        "../data/test/graph/test_dumping_by_graph_and_reload_by_sv_trio_person0.falist",
        fq_quality_cutoff, homopolymer_cutoff,
        remove_duplicates_se, ascii_fq_offset,
        into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
        &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded0,
        NULL, 0, &subsample_null);

      //GraphInfo* ginfo = graph_info_alloc_and_init();

      graph_info_set_seq(ginfo, 0, seq_loaded0);
      db_graph_dump_single_colour_binary_of_colour0(
        "../data/tempfiles_can_be_deleted/test_dumping_by_graph_and_reload_by_sv_trio_person0_kmer31.ctx", 
				&db_node_condition_always_true, db_graph, ginfo, BINVERSION);
      hash_table_free(&db_graph);


      db_graph = hash_table_new(number_of_bits, bucket_size,
                                max_retries, kmer_size);
     
      load_se_filelist_into_graph_colour(
      "../data/test/graph/test_dumping_by_graph_and_reload_by_sv_trio_person1.falist",
      fq_quality_cutoff, homopolymer_cutoff,
      remove_duplicates_se, ascii_fq_offset,
      into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
      &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded1,
      NULL, 0, &subsample_null);

      graph_info_initialise(ginfo);
      graph_info_set_seq(ginfo, 0, seq_loaded1);

      db_graph_dump_single_colour_binary_of_colour0(
        "../data/tempfiles_can_be_deleted/test_dumping_by_graph_and_reload_by_sv_trio_person1_kmer31.ctx", 
        &db_node_condition_always_true, db_graph, ginfo, BINVERSION);
     hash_table_free(&db_graph);


     db_graph = hash_table_new(number_of_bits, bucket_size,
                               max_retries, kmer_size);
     
     load_se_filelist_into_graph_colour(
      "../data/test/graph/test_dumping_by_graph_and_reload_by_sv_trio_person2.falist",
      fq_quality_cutoff, homopolymer_cutoff,
      remove_duplicates_se, ascii_fq_offset,
      into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
      &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded2,
      NULL, 0, &subsample_null);

      graph_info_initialise(ginfo);
      graph_info_set_seq(ginfo, 0, seq_loaded2);

      db_graph_dump_single_colour_binary_of_colour0(
        "../data/tempfiles_can_be_deleted/test_dumping_by_graph_and_reload_by_sv_trio_person2_kmer31.ctx", 
        &db_node_condition_always_true, db_graph, ginfo, BINVERSION);
      hash_table_free(&db_graph);


      db_graph = hash_table_new(number_of_bits, bucket_size,
                                max_retries, kmer_size);

      graph_info_initialise(ginfo);

      int first_colour = 0;
      load_population_as_binaries_from_graph(
        "../data/test/graph/trio_filelist_for_testing_loading_singlecolour_bins_into_multicol_bin.colours",
        first_colour, true, db_graph, ginfo,
        false, 0, false);

      CU_ASSERT(ginfo->total_sequence[0]==seq_loaded0);
      CU_ASSERT(ginfo->total_sequence[1]==seq_loaded1);
      CU_ASSERT(ginfo->total_sequence[2]==seq_loaded2);
      
     CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph), 41);
     
     
      // Start with a kmer that should be in person0 and person 2 only
      BinaryKmer tmp_kmer1, tmp_kmer2;
     
      dBNode* test_element1_person0
        = db_graph_find_node_restricted_to_specific_person_or_population(
					  element_get_key(seq_to_binary_kmer("GTGGGAGGATCGCTTGAGTCCAGGAGTTCTG", 
                                               kmer_size, &tmp_kmer1),
                            kmer_size, &tmp_kmer2),
            db_graph, 0);
      
      //look for reverse complement of GTGGGAGGATCGCTTGAGTCCAGGAGTTCTG
      dBNode* test_element2_person0
        = db_graph_find_node_restricted_to_specific_person_or_population(
            element_get_key(seq_to_binary_kmer("CAGAACTCCTGGACTCAAGCGATCCTCCCAC",
                                               kmer_size, &tmp_kmer1),
                            kmer_size, &tmp_kmer2),
            db_graph, 0);
      
      dBNode* test_element1_person1
        = db_graph_find_node_restricted_to_specific_person_or_population(
            element_get_key(seq_to_binary_kmer("GTGGGAGGATCGCTTGAGTCCAGGAGTTCTG",
                                               kmer_size, &tmp_kmer1),
                            kmer_size, &tmp_kmer2),
            db_graph, 1);
      
      //look for reverse complement of GTGGGAGGATCGCTTGAGTCCAGGAGTTCTG
      dBNode* test_element2_person1
        = db_graph_find_node_restricted_to_specific_person_or_population(
            element_get_key(seq_to_binary_kmer("CAGAACTCCTGGACTCAAGCGATCCTCCCAC",
                                               kmer_size, &tmp_kmer1),
                            kmer_size, &tmp_kmer2),
            db_graph, 1);
      
      
      dBNode* test_element1_person2
        = db_graph_find_node_restricted_to_specific_person_or_population(
            element_get_key(seq_to_binary_kmer("GTGGGAGGATCGCTTGAGTCCAGGAGTTCTG",
                                               kmer_size, &tmp_kmer1),
                            kmer_size, &tmp_kmer2),
            db_graph, 2);
      
      //look for reverse complement of GTGGGAGGATCGCTTGAGTCCAGGAGTTCTG
      dBNode* test_element2_person2
        = db_graph_find_node_restricted_to_specific_person_or_population(
            element_get_key(seq_to_binary_kmer("CAGAACTCCTGGACTCAAGCGATCCTCCCAC",
                                               kmer_size, &tmp_kmer1),
                            kmer_size, &tmp_kmer2),
            db_graph, 2);

      CU_ASSERT(test_element1_person0 !=NULL);
      CU_ASSERT(test_element2_person0 !=NULL);
      CU_ASSERT(test_element1_person0==test_element2_person0);
      CU_ASSERT(test_element1_person1 ==NULL);
      CU_ASSERT(test_element2_person1 ==NULL);
      CU_ASSERT(test_element1_person2 !=NULL);
      CU_ASSERT(test_element2_person2 !=NULL);
      CU_ASSERT(test_element1_person2==test_element2_person0);
      CU_ASSERT(test_element1_person0==test_element1_person2);
      
      
      //Now some kmers that exist in person 2 only
      
      dBNode* test_element3_person0
        = db_graph_find_node_restricted_to_specific_person_or_population(
            element_get_key(seq_to_binary_kmer("GATCGGGTGTCCGCACTAAGTTCGGCATCAA",
                            kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),
            db_graph, 0);

      dBNode* test_element3_person1
      = db_graph_find_node_restricted_to_specific_person_or_population(
          element_get_key(seq_to_binary_kmer("GATCGGGTGTCCGCACTAAGTTCGGCATCAA",
                          kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),
          db_graph, 1);

      dBNode* test_element3_person2
        = db_graph_find_node_restricted_to_specific_person_or_population(
            element_get_key(seq_to_binary_kmer("GATCGGGTGTCCGCACTAAGTTCGGCATCAA",
                            kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),
            db_graph, 2);
      
      
      CU_ASSERT(test_element3_person0==NULL);
      CU_ASSERT(test_element3_person1!=NULL);
      CU_ASSERT(test_element3_person2==NULL);
      
      if (test_element3_person1==NULL)
      {
        warn("%s:%i: test_element3_person1 is NULL", __FILE__, __LINE__);
      }
      if (test_element1_person2==NULL)
      {
        warn("%s:%i: test_element1_person2 is NULL", __FILE__, __LINE__);
      }
      
      //now check the edges
      
      //in person1 we know this kmer GATCGGGTGTCCGCACTAAGTTCGGCATCAA has an Thymine edge. This is in the forward direction with respect to the node, as the reverse complement
      //  of GATCGGGTGTCCGCACTAAGTTCGGCATCAA starts with an T, and so is bigger
      Nucleotide base;
      CU_ASSERT(db_node_has_precisely_one_edge(test_element3_person1, forward,&base, 1)==true);
      CU_ASSERT(db_node_edge_exist(test_element3_person1, Thymine, forward, 1)==true);
      
      
      //and the coup de grace. Does this kmer GTGGGAGGATCGCTTGAGTCCAGGAGTTCTG node have two different edges, for person 0 and person 2. ie does person 0 have precisely a G,
      // and person 2 an A?
      
      CU_ASSERT(db_node_has_precisely_one_edge(test_element1_person0, reverse,&base, 0)==true);
      CU_ASSERT(db_node_edge_exist(test_element1_person0, Guanine, reverse, 0)==true);
      CU_ASSERT(db_node_has_precisely_one_edge(test_element1_person2, reverse,&base, 2)==true);
      CU_ASSERT(db_node_edge_exist(test_element1_person2, Adenine, reverse, 2)==true);
      
      
      //from now on we are pretty happy. It's loading the right nodes, and putting the edges in the right place.
      //the following is not ideal from a modular code point of view - I'm going to use code in db_graph_population to get supernodes, and see if they are what I expect.
      
      
      int max_expected_supernode_length=30;
      dBNode * nodes_path[max_expected_supernode_length];
      Orientation orientations_path[max_expected_supernode_length];
      Nucleotide labels_path[max_expected_supernode_length];
      char seq[max_expected_supernode_length+db_graph->kmer_size+1];
      
      double avg_coverage = 0;
      Covg max_coverage = 0;
      Covg min_coverage = 0;
      boolean is_cycle = false;
      
      
      db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop(
        test_element1_person0, reverse, max_expected_supernode_length, Guanine,
        &db_node_action_do_nothing, nodes_path, orientations_path, labels_path,
        seq, &avg_coverage, &min_coverage, &max_coverage, &is_cycle, db_graph,
        0);

      CU_ASSERT_STRING_EQUAL("GGCTGTAGTGCGCTATGCC", seq);

      db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop(
        test_element1_person0, reverse, max_expected_supernode_length, Adenine,
        &db_node_action_do_nothing, nodes_path,orientations_path, labels_path,
        seq, &avg_coverage, &min_coverage, &max_coverage, &is_cycle, db_graph,
        2);

      CU_ASSERT_STRING_EQUAL("A", seq);

      graph_info_free(ginfo);
      hash_table_free(&db_graph);
    }

    // Now check the edges

    // In person1 we know this kmer GATCGGGTGTCCGCACTAAGTTCGGCATCAA has an
    // Thymine edge. This is in the forward direction with respect to the node,
    // as the reverse complement of GATCGGGTGTCCGCACTAAGTTCGGCATCAA starts
    // with an T, and so is bigger
    Nucleotide base;
    CU_ASSERT(db_node_has_precisely_one_edge(test_element3_person1, forward,&base, 1)==true);
    CU_ASSERT(db_node_edge_exist(test_element3_person1, Thymine, forward, 1)==true);


    //and the coup de grace. Does this kmer GTGGGAGGATCGCTTGAGTCCAGGAGTTCTG node
    // have two different edges, for person 0 and person 2. ie does person 0
    // have precisely a G, and person 2 an A?

    CU_ASSERT(db_node_has_precisely_one_edge(test_element1_person0, reverse,&base, 0)==true);
    CU_ASSERT(db_node_edge_exist(test_element1_person0, Guanine, reverse, 0)==true);
    CU_ASSERT(db_node_has_precisely_one_edge(test_element1_person2, reverse,&base, 2)==true);
    CU_ASSERT(db_node_edge_exist(test_element1_person2, Adenine, reverse, 2)==true);


    // From now on we are pretty happy. It's loading the right nodes, and
    // putting the edges in the right place. The following is not ideal from a
    // modular code point of view - I'm going to use code in db_graph_population
    // to get supernodes, and see if they are what I expect.

    int max_expected_supernode_length = 30;
    dBNode * nodes_path[max_expected_supernode_length];
    Orientation orientations_path[max_expected_supernode_length];
    Nucleotide labels_path[max_expected_supernode_length];
    char seq[max_expected_supernode_length+db_graph->kmer_size+1];

    double avg_coverage = 0;
    Covg max_coverage = 0;
    Covg min_coverage = 0;
    boolean is_cycle = false;


    db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop(test_element1_person0, reverse, max_expected_supernode_length, Guanine, &db_node_action_do_nothing,
      nodes_path,orientations_path, labels_path, seq,
      &avg_coverage, &min_coverage, &max_coverage, &is_cycle,
      db_graph, 0);

    CU_ASSERT_STRING_EQUAL("GGCTGTAGTGCGCTATGCC", seq);


    db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop(test_element1_person0, reverse, max_expected_supernode_length, Adenine, &db_node_action_do_nothing,
      nodes_path,orientations_path, labels_path, seq,
      &avg_coverage, &min_coverage, &max_coverage, &is_cycle,
      db_graph, 2);

    CU_ASSERT_STRING_EQUAL("A", seq);

    hash_table_free(&db_graph);
    graph_info_free(ginfo);
  }
}

void test_coverage_is_correctly_counted_on_loading_from_file()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // First set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size = 10;

  dBGraph * db_graph = hash_table_new(number_of_bits, bucket_size,
                                      10, kmer_size);

  int fq_quality_cutoff = 20;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_read = 0, seq_loaded = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/graph/file_to_test_covg_of_reads.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0,
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  /*
  >occurs once
  AATAATAATAATAATAATAATAATAAT
  >occurs twice
  GGGCAGTCTCT
  >occurs twice
  GGGCAGTCTCT
  >also occurs once
  TTTTTTTTTT
  */

  CU_ASSERT_EQUAL(seq_read, 59);
  CU_ASSERT_EQUAL(seq_loaded, 59);

  CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),12);
  //CU_ASSERT_EQUAL(bad_reads,0);

  BinaryKmer tmp_kmer;
  BinaryKmer tmp_kmer2;

  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("TTT", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(db_node_get_coverage(test_element1,0)==8);
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AAT", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(db_node_get_coverage(test_element1,0)==9);

  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(db_node_get_coverage(test_element1,0)==2);
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GGC", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(db_node_get_coverage(test_element1,0)==2);
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GCA", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(db_node_get_coverage(test_element1,0)==2);
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("CAG", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(db_node_get_coverage(test_element1,0)==2)
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AGT", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(db_node_get_coverage(test_element1,0)==2);
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GTC", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(db_node_get_coverage(test_element1,0)==2);
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("TCT", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(db_node_get_coverage(test_element1,0)==4);

  hash_table_free(&db_graph);
}


void test_getting_sliding_windows_where_you_break_at_kmers_not_in_db_graph()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // Create graph
  int kmer_size = 17;
  int number_of_bits = 10;
  int bucket_size = 30;

  dBGraph* db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);

  // Load file
  int fq_quality_cutoff = 20;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_read = 0, seq_loaded = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/graph/person3.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0,
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  //OK - we have graph. Now test getting sliding windows from
  // 1. a sequence that is all in the graph
  // 2. a sequence that has one bad base in the middle
  // 3. total garbage sequence

  FILE* fp = fopen("../data/test/graph/person3_with_errors_extended_file.fq", "r");

  if(fp == NULL)
  {
    die("Cannot open ../data/test/graph/person3_with_errors_extended_file.fq\n"
"in test_getting_sliding_windows_where_you_break_at_kmers_not_in_db_graph");
  }


  //allocations
  int max_read_length=100;
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence");
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);


  //max_read_length/(kmer_size+1) is the worst case for the number of sliding windows, ie a kmer follow by a low-quality/bad base
  int max_windows = max_read_length/(kmer_size+1);

  //number of possible kmers in a 'perfect' read
  int max_kmers   = max_read_length-kmer_size+1;


  //----------------------------------
  //preallocate the space of memory used to keep the sliding_windows. NB: this space of memory is reused for every call -- with the view
  //to avoid memory fragmentation
  //NB: this space needs to preallocate memory for orthogonal situations:
  //    * a good read -> few windows, many kmers per window
  //    * a bad read  -> many windows, few kmers per window
  //----------------------------------
  KmerSlidingWindowSet * windows = malloc(sizeof(KmerSlidingWindowSet));
  if (windows == NULL){
    die("Out of memory trying to allocate a KmerArraySet");
  }
  //allocate memory for the sliding windows
  binary_kmer_alloc_kmers_set(windows, max_windows, max_kmers);




  // GET READ 1 - this si full of sequencing errors, and it should be impossible to find any kmers that are in the graph
  int len = read_sequence_from_fastq(fp, seq, max_read_length, ascii_fq_offset);

  get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq, seq->qual, len, fq_quality_cutoff,
    windows, max_windows, max_kmers, db_graph);
  CU_ASSERT(windows->nwindows==0);



  // GET READ 2 - this one lies entirely in the graph
  len = read_sequence_from_fastq(fp, seq, max_read_length, ascii_fq_offset);
  get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq, seq->qual, len, fq_quality_cutoff,
    windows, max_windows, max_kmers, db_graph);
  CU_ASSERT(windows->nwindows==1);
  CU_ASSERT((windows->window[0]).nkmers=28);

  BinaryKmer test_kmer;
  seq_to_binary_kmer("ACCCTAACCCTAACCCT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true );
  seq_to_binary_kmer("CCCTAACCCTAACCCTA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[1],test_kmer)==true );
  seq_to_binary_kmer("CCTAACCCTAACCCTAA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[2],test_kmer)==true );

  seq_to_binary_kmer("CTAACCCTAACCCTAAC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[3],test_kmer)==true );
  seq_to_binary_kmer("TAACCCTAACCCTAACC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[4],test_kmer)==true );
  seq_to_binary_kmer("AACCCTAACCCTAACCC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[5],test_kmer)==true );
  seq_to_binary_kmer("ACCCTAACCCTAACCCC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[6],test_kmer)==true );

  seq_to_binary_kmer("CCCTAACCCTAACCCCT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[7],test_kmer)==true );
  seq_to_binary_kmer("CCTAACCCTAACCCCTA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[8],test_kmer)==true );
  seq_to_binary_kmer("CTAACCCTAACCCCTAA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[9],test_kmer)==true );


  seq_to_binary_kmer("TAACCCTAACCCCTAAC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[10],test_kmer)==true );
  seq_to_binary_kmer("AACCCTAACCCCTAACC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[11],test_kmer)==true );
  seq_to_binary_kmer("ACCCTAACCCCTAACCC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[12],test_kmer)==true );

  seq_to_binary_kmer("CCCTAACCCCTAACCCT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[13],test_kmer)==true );
  seq_to_binary_kmer("CCTAACCCCTAACCCTA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[14],test_kmer)==true );
  seq_to_binary_kmer("CTAACCCCTAACCCTAA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[15],test_kmer)==true );

  // ..not going all the way to the end.
  // move on to next read



  // GET READ 3 - lots of errors again
  len = read_sequence_from_fastq(fp, seq, max_read_length,ascii_fq_offset);
  get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq, seq->qual, len, fq_quality_cutoff,
    windows, max_windows, max_kmers, db_graph);
  CU_ASSERT(windows->nwindows==0);


  // GET fourth read - this one has a single base in the middle which means kmers containing it won't be in thr graph. Otherwise the same as read3 in person3.fa

  len = read_sequence_from_fastq(fp, seq, max_read_length,ascii_fq_offset);
  get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq, seq->qual, len, fq_quality_cutoff,
    windows, max_windows, max_kmers, db_graph);
  CU_ASSERT(windows->nwindows==2);
  CU_ASSERT((windows->window[0]).nkmers==3);

  seq_to_binary_kmer("GGGGCGGGGCGGGGCGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true );
  seq_to_binary_kmer("GGGCGGGGCGGGGCGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[1],test_kmer)==true );
  seq_to_binary_kmer("GGCGGGGCGGGGCGGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[2],test_kmer)==true );


  //and in the second window, after the interrupting T:

  CU_ASSERT((windows->window[1]).nkmers==8);
  seq_to_binary_kmer("GGGGCGGGGCCCCCTCA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[0],test_kmer)==true );
  seq_to_binary_kmer("GGGCGGGGCCCCCTCAC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[1],test_kmer)==true );
  seq_to_binary_kmer("GGCGGGGCCCCCTCACA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[2],test_kmer)==true );
  seq_to_binary_kmer("GCGGGGCCCCCTCACAC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[3],test_kmer)==true );
  seq_to_binary_kmer("CGGGGCCCCCTCACACA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[4],test_kmer)==true );
  seq_to_binary_kmer("GGGGCCCCCTCACACAC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[5],test_kmer)==true );
  seq_to_binary_kmer("GGGCCCCCTCACACACA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[6],test_kmer)==true );
  seq_to_binary_kmer("GGCCCCCTCACACACAT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[7],test_kmer)==true );


  // GET FIFTH READ - lies entirely in graph, so should just get one window

  len = read_sequence_from_fastq(fp, seq, max_read_length,ascii_fq_offset);
  get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq, seq->qual, len, fq_quality_cutoff,
    windows, max_windows, max_kmers, db_graph);
  CU_ASSERT(windows->nwindows==1);
  CU_ASSERT((windows->window[0]).nkmers==28);

  seq_to_binary_kmer("GGGGCGGGGCGGGGCGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true );
  seq_to_binary_kmer("GGGCGGGGCGGGGCGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[1],test_kmer)==true );
  seq_to_binary_kmer("GGCGGGGCGGGGCGGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[2],test_kmer)==true );
  seq_to_binary_kmer("GCGGGGCGGGGCGGGGC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[3],test_kmer)==true );
  seq_to_binary_kmer("CGGGGCGGGGCGGGGCG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[4],test_kmer)==true );
  seq_to_binary_kmer("GGGGCGGGGCGGGGCGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[5],test_kmer)==true );

  //... etc ...

  seq_to_binary_kmer("GGCCCCCTCACACACAT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[27],test_kmer)==true );


  // READ 6 - entirely in graph
  len = read_sequence_from_fastq(fp, seq, max_read_length,ascii_fq_offset);
  get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq, seq->qual, len, fq_quality_cutoff,
    windows, max_windows, max_kmers, db_graph);
  CU_ASSERT(windows->nwindows==1);
  CU_ASSERT((windows->window[0]).nkmers==15);

  seq_to_binary_kmer("TTTTTTTTTTTTTTTTT", kmer_size, &test_kmer);
  int i;
  for (i=0; i<15; i++)
  {
    CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[i],test_kmer)==true );
  }

  // READ 7 - first character wrong. Vital test this one. Slightly different code path if the first kmer of all is bad.
  len = read_sequence_from_fastq(fp, seq, max_read_length,ascii_fq_offset);
  get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq, seq->qual, len, fq_quality_cutoff,
    windows, max_windows, max_kmers, db_graph);
  CU_ASSERT(windows->nwindows==1);
  CU_ASSERT((windows->window[0]).nkmers==27);

  seq_to_binary_kmer("GGGCGGGGCGGGGCGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true );
  seq_to_binary_kmer("GGCGGGGCGGGGCGGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[1],test_kmer)==true );
  seq_to_binary_kmer("GCGGGGCGGGGCGGGGC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[2],test_kmer)==true );
  seq_to_binary_kmer("CGGGGCGGGGCGGGGCG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[3],test_kmer)==true );

  // etc

  seq_to_binary_kmer("GGCCCCCTCACACACAT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[26],test_kmer)==true );



  // READ 7 - errors spaces such that precisely two kmers can be pulled out, each in their own window

  len = read_sequence_from_fastq(fp, seq, max_read_length,ascii_fq_offset);
  get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq, seq->qual, len, fq_quality_cutoff,
    windows, max_windows, max_kmers, db_graph);
  CU_ASSERT(windows->nwindows==2);
  CU_ASSERT((windows->window[0]).nkmers==1);

  seq_to_binary_kmer("GGGCGGGGCGGGGCGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true );


  CU_ASSERT((windows->window[1]).nkmers==1);
  seq_to_binary_kmer("GGGGCGGGGCCCCCTCA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[0],test_kmer)==true );



  // last read contains two kmers that are in the graph, but which have no edge between them. This function is not sensitive to that, and should just get one window

  len = read_sequence_from_fastq(fp, seq, max_read_length,ascii_fq_offset);
  get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq, seq->qual, len, fq_quality_cutoff,
    windows, max_windows, max_kmers, db_graph);

  CU_ASSERT(windows->nwindows==1);
  CU_ASSERT((windows->window[0]).nkmers==5);

  seq_to_binary_kmer("TTTTTTTTTTTTTTTTT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true );
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[1],test_kmer)==true );
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[2],test_kmer)==true );
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[3],test_kmer)==true );
  seq_to_binary_kmer("TTTTTTTTTTTTTTTTA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[4],test_kmer)==true );


  free_sequence(&seq);
  binary_kmer_free_kmers_set(&windows);
  hash_table_free(&db_graph);
}




void test_get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  int kmer_size = 17;
  int number_of_bits = 10;
  int bucket_size = 30;

  dBGraph* db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);

  // Load file
  int fq_quality_cutoff = 20;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_read = 0, seq_loaded = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/graph/person3.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0,
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  //OK - we have graph. Now test getting sliding windows from
  // 1. a sequence that is all in the graph
  // 2. a sequence that has one bad base in the middle
  // 3. total garbage sequence

  FILE* fp = fopen("../data/test/graph/person3_with_errors_extended_file.fq", "r");
  if(fp == NULL)
  {
    die("Cannot open ../data/test/graph/person3_with_errors_extended_file.fq\n"
        "in test_getting_sliding_windows_where_you_break_at_kmers_not_in_db_graph\n");
  }


  //allocations
  int max_read_length=100;
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence");
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);


  //max_read_length/(kmer_size+1) is the worst case for the number of sliding windows, ie a kmer follow by a low-quality/bad base
  int max_windows = max_read_length/(kmer_size+1);

  //number of possible kmers in a 'perfect' read
  int max_kmers   = max_read_length-kmer_size+1;



  //----------------------------------
  //preallocate the space of memory used to keep the sliding_windows. NB: this space of memory is reused for every call -- with the view
  //to avoid memory fragmentation
  //NB: this space needs to preallocate memory for orthogonal situations:
  //    * a good read -> few windows, many kmers per window
  //    * a bad read  -> many windows, few kmers per window
  //----------------------------------
  KmerSlidingWindowSet * windows = malloc(sizeof(KmerSlidingWindowSet));
  if (windows == NULL){
    die("Out of memory trying to allocate a KmerArraySet");
  }
  //allocate memory for the sliding windows
  binary_kmer_alloc_kmers_set(windows, max_windows, max_kmers);


  // GET READ 1 - this si full of sequencing errors, and it should be impossible to find any kmers that are in the graph
  int len = read_sequence_from_fastq(fp, seq, max_read_length,ascii_fq_offset);

  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, fq_quality_cutoff,
    windows, max_windows, max_kmers, db_graph, 0);
  CU_ASSERT(windows->nwindows==0);



  // GET READ 2 - this one lies entirely in the graph
  len = read_sequence_from_fastq(fp, seq, max_read_length,ascii_fq_offset);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, fq_quality_cutoff,
    windows, max_windows, max_kmers, db_graph, 0);
  CU_ASSERT(windows->nwindows==1);
  CU_ASSERT((windows->window[0]).nkmers=28);

  BinaryKmer test_kmer;
  seq_to_binary_kmer("ACCCTAACCCTAACCCT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true );
  seq_to_binary_kmer("CCCTAACCCTAACCCTA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[1],test_kmer)==true );
  seq_to_binary_kmer("CCTAACCCTAACCCTAA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[2],test_kmer)==true );

  seq_to_binary_kmer("CTAACCCTAACCCTAAC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[3],test_kmer)==true );
  seq_to_binary_kmer("TAACCCTAACCCTAACC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[4],test_kmer)==true );
  seq_to_binary_kmer("AACCCTAACCCTAACCC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[5],test_kmer)==true );
  seq_to_binary_kmer("ACCCTAACCCTAACCCC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[6],test_kmer)==true );

  seq_to_binary_kmer("CCCTAACCCTAACCCCT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[7],test_kmer)==true );
  seq_to_binary_kmer("CCTAACCCTAACCCCTA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[8],test_kmer)==true );
  seq_to_binary_kmer("CTAACCCTAACCCCTAA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[9],test_kmer)==true );


  seq_to_binary_kmer("TAACCCTAACCCCTAAC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[10],test_kmer)==true );
  seq_to_binary_kmer("AACCCTAACCCCTAACC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[11],test_kmer)==true );
  seq_to_binary_kmer("ACCCTAACCCCTAACCC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[12],test_kmer)==true );

  seq_to_binary_kmer("CCCTAACCCCTAACCCT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[13],test_kmer)==true );
  seq_to_binary_kmer("CCTAACCCCTAACCCTA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[14],test_kmer)==true );
  seq_to_binary_kmer("CTAACCCCTAACCCTAA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[15],test_kmer)==true );

  // ..not going all the way to the end.
  // move on to next read



  // GET READ 3 - lots of errors again
  len = read_sequence_from_fastq(fp, seq, max_read_length,ascii_fq_offset);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, fq_quality_cutoff,
    windows, max_windows, max_kmers, db_graph, 0);
  CU_ASSERT(windows->nwindows==0);


  // GET fourth read - this one has a single base in the middle which means kmers containing it won't be in thr graph. Otherwise the same as read3 in person3.fa

  len = read_sequence_from_fastq(fp, seq, max_read_length, ascii_fq_offset);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, fq_quality_cutoff,
   windows, max_windows, max_kmers, db_graph, 0);
  CU_ASSERT(windows->nwindows==2);
  CU_ASSERT((windows->window[0]).nkmers==3);

  seq_to_binary_kmer("GGGGCGGGGCGGGGCGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true );
  seq_to_binary_kmer("GGGCGGGGCGGGGCGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[1],test_kmer)==true );
  seq_to_binary_kmer("GGCGGGGCGGGGCGGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[2],test_kmer)==true );

  //and in the second window, after the interrupting T:

  CU_ASSERT((windows->window[1]).nkmers==8);
  seq_to_binary_kmer("GGGGCGGGGCCCCCTCA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[0],test_kmer)==true );
  seq_to_binary_kmer("GGGCGGGGCCCCCTCAC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[1],test_kmer)==true );
  seq_to_binary_kmer("GGCGGGGCCCCCTCACA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[2],test_kmer)==true );
  seq_to_binary_kmer("GCGGGGCCCCCTCACAC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[3],test_kmer)==true );
  seq_to_binary_kmer("CGGGGCCCCCTCACACA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[4],test_kmer)==true );
  seq_to_binary_kmer("GGGGCCCCCTCACACAC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[5],test_kmer)==true );
  seq_to_binary_kmer("GGGCCCCCTCACACACA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[6],test_kmer)==true );
  seq_to_binary_kmer("GGCCCCCTCACACACAT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[7],test_kmer)==true );


  // GET FIFTH READ - lies entirely in graph, so should just get one window

  len = read_sequence_from_fastq(fp, seq, max_read_length, ascii_fq_offset);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, fq_quality_cutoff,
    windows, max_windows, max_kmers, db_graph, 0);
  CU_ASSERT(windows->nwindows==1);
  CU_ASSERT((windows->window[0]).nkmers==28);

  seq_to_binary_kmer("GGGGCGGGGCGGGGCGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true );
  seq_to_binary_kmer("GGGCGGGGCGGGGCGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[1],test_kmer)==true );
  seq_to_binary_kmer("GGCGGGGCGGGGCGGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[2],test_kmer)==true );
  seq_to_binary_kmer("GCGGGGCGGGGCGGGGC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[3],test_kmer)==true );
  seq_to_binary_kmer("CGGGGCGGGGCGGGGCG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[4],test_kmer)==true );
  seq_to_binary_kmer("GGGGCGGGGCGGGGCGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[5],test_kmer)==true );

  //... etc ...

  seq_to_binary_kmer("GGCCCCCTCACACACAT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[27],test_kmer)==true );


  // READ 6 - entirely in graph
  len = read_sequence_from_fastq(fp, seq, max_read_length, ascii_fq_offset);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, fq_quality_cutoff,
    windows, max_windows, max_kmers, db_graph, 0);
  CU_ASSERT(windows->nwindows==1);
  CU_ASSERT((windows->window[0]).nkmers==15);

  seq_to_binary_kmer("TTTTTTTTTTTTTTTTT", kmer_size, &test_kmer);
  int i;
  for (i=0; i<15; i++)
  {
    CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[i],test_kmer)==true );
  }

  // READ 7 - first character wrong. Vital test this one. Slightly different code path if the first kmer of all is bad.
  len = read_sequence_from_fastq(fp, seq, max_read_length, ascii_fq_offset);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, fq_quality_cutoff,
    windows, max_windows, max_kmers, db_graph, 0);
  CU_ASSERT(windows->nwindows==1);
  CU_ASSERT((windows->window[0]).nkmers==27);

  seq_to_binary_kmer("GGGCGGGGCGGGGCGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true );
  seq_to_binary_kmer("GGCGGGGCGGGGCGGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[1],test_kmer)==true );
  seq_to_binary_kmer("GCGGGGCGGGGCGGGGC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[2],test_kmer)==true );
  seq_to_binary_kmer("CGGGGCGGGGCGGGGCG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[3],test_kmer)==true );

  // etc

  seq_to_binary_kmer("GGCCCCCTCACACACAT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[26],test_kmer)==true );


  // READ 7 - errors spaced apart

  len = read_sequence_from_fastq(fp, seq, max_read_length, ascii_fq_offset);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, fq_quality_cutoff,
    windows, max_windows, max_kmers, db_graph, 0);
  CU_ASSERT(windows->nwindows==2);

  //last read - we shoud get two separate windows, as the kmers have no edge joining them
  len = read_sequence_from_fastq(fp, seq, max_read_length, ascii_fq_offset);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, fq_quality_cutoff,
    windows, max_windows, max_kmers, db_graph, 0);
  CU_ASSERT(windows->nwindows==2);

  CU_ASSERT((windows->window[0]).nkmers==4);
  seq_to_binary_kmer("TTTTTTTTTTTTTTTTT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true );
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[1],test_kmer)==true );
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[2],test_kmer)==true );
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[3],test_kmer)==true );


  CU_ASSERT((windows->window[1]).nkmers==1);
  seq_to_binary_kmer("TTTTTTTTTTTTTTTTA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[0],test_kmer)==true );

  free_sequence(&seq);
  binary_kmer_free_kmers_set(&windows);
  hash_table_free(&db_graph);
}





// DEV: propose this is removed or updated
//      - read_fastq_and_print_reads_that_lie_in_graph() isn't used anywhere else
//      - Doesn't use new seq_file sequence reader
//
//Assumption is that you use a bunch of fastq to build a graph, then clean it.
//You then want access to a set of fasta files that correspond to the good reads only.
void test_dumping_of_clean_fasta()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // Construct hash table
  int kmer_size = 17;
  int number_of_bits = 10;
  int bucket_size = 30;

  dBGraph* db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);

  // Read sequence
  int fq_quality_cutoff = 20;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_read = 0, seq_loaded = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/graph/person3.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0,
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  // OK, we now have a graph. Let's see if dumping a clean fasta gives the
  // right answers.

  FILE* fptr = fopen("../data/test/graph/person3_with_errors.fq", "r");
  if(fptr == NULL)
  {
    die("Cannot open ../data/test/graph/person3_with_errors.fq in\n"
        "test_dumping_of_clean_fasta. Exiting.");
  }

  int file_reader(FILE * fp, Sequence * seq, int max_read_length,
                  boolean new_entry, boolean * full_entry)
  {
    * full_entry = true;

    if (new_entry!= true){
      die("new_entry has to be true for fastq");
    }

    return read_sequence_from_fastq(fp,seq,max_read_length, ascii_fq_offset);
  }


  int max_read_length = 100;

  //alloc array to hold results
  char** array_of_reads= (char**) calloc(30,sizeof(char*));
  int i;
  for (i=0; i<30; i++)
  {
    array_of_reads[i]= (char*)calloc(100,sizeof(char));
  }

  int number_of_reads = 0;
  long long fastq_bad_reads = 0;

  read_fastq_and_print_reads_that_lie_in_graph(fptr, stdout, &file_reader,
    &fastq_bad_reads, max_read_length, db_graph,
					       //false, NULL, NULL);
    true, array_of_reads, &number_of_reads);

  CU_ASSERT_STRING_EQUAL(array_of_reads[0], "ACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAAC");
  CU_ASSERT_STRING_EQUAL(array_of_reads[1], "GGGGCGGGGCGGGGCGGGG");
  CU_ASSERT_STRING_EQUAL(array_of_reads[2], "GGGGCGGGGCCCCCTCACACACAT");

  CU_ASSERT_STRING_EQUAL(array_of_reads[3], "GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT");
  CU_ASSERT_STRING_EQUAL(array_of_reads[4], "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");

  for (i=0; i<30; i++)
  {
    free(array_of_reads[i]);
  }
  free(array_of_reads);
  hash_table_free(&db_graph);

  /*> read
ACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAAC
Examing this:  AGGGCGGGGCGGGGCAGGGCAGGGCGGAGCCCACTCACACACAT
Examing this:  GGGGCGGGGCGGGGCGGGGTGGGGCGGGGCCCCCTCACACACAT
Examing this:  GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
> read
GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
Examing this:  TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
> read
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  */


}




void test_loading_of_paired_end_reads_removing_duplicates()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // First set up the hash/graph
  int kmer_size = 21;
  int number_of_bits = 10;
  int bucket_size = 10;

  dBGraph * db_graph = hash_table_new(number_of_bits, bucket_size,
                                      10, kmer_size);

  // first a test where the file contains no duplicates and you do not try to
  // remove duplicates - does all the data get loaded?

  // Read in sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_pe = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int file_pairs_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_read = 0, seq_loaded = 0;

  load_pe_filelists_into_graph_colour(
    "../data/test/graph/paired_end_file1_1.fqlist",
    "../data/test/graph/paired_end_file1_2.fqlist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_pe, ascii_fq_offset,
    into_colour, db_graph, 0,
    &file_pairs_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  // Because there are various N's in the sequence, we lose some bases
  CU_ASSERT(seq_read == 720);
  CU_ASSERT(seq_loaded == 665);
  CU_ASSERT(db_graph->unique_kmers == 243);

  // Then a test where you try to remove duplicates from files that don't
  // contain any
  hash_table_free(&db_graph);
  db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);

  remove_duplicates_pe = true;

  // Reset counters
  file_pairs_loaded = 0;
  bad_reads = 0;
  dup_reads = 0;
  seq_read = 0;
  seq_loaded = 0;

  load_pe_filelists_into_graph_colour(
    "../data/test/graph/paired_end_file1_1.fqlist",
    "../data/test/graph/paired_end_file1_2.fqlist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_pe, ascii_fq_offset,
    into_colour, db_graph, 0,
    &file_pairs_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT(seq_read == 720);
  CU_ASSERT(seq_loaded == 665);
  CU_ASSERT(db_graph->unique_kmers == 243);
  CU_ASSERT(dup_reads == 0);

  hash_table_free(&db_graph);

  // Now try a pair of files containing one duplicate pair of reads
  // (corresponding mates duplicated in both files) and also one read that is
  // duplicated only in the _1 file only the former pair of reads should be
  // discarded

  db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);

  // Reset counters
  file_pairs_loaded = 0;
  bad_reads = 0;
  dup_reads = 0;
  seq_read = 0;
  seq_loaded = 0;

  load_pe_filelists_into_graph_colour(
    "../data/test/graph/paired_end_file2_with_dup_1.fqlist",
    "../data/test/graph/paired_end_file2_with_dup_2.fqlist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_pe, ascii_fq_offset,
    into_colour, db_graph, 0,
    &file_pairs_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  // As before, but with four more 36bp reads
  // All the sequence we have read is loaded except for 2 reads, plus the effect
  // of the Ns
  CU_ASSERT(seq_read == 864);
  CU_ASSERT(seq_loaded == 737);
  CU_ASSERT(db_graph->unique_kmers == 245);
  CU_ASSERT(dup_reads == 2);

  hash_table_free(&db_graph);

  // now take a pair of files where there is an original mate pair, one proper
  // duplicate, and then 3 other mate-pairs that might be confused for
  // duplicates - one read is a dup of the original, and the other is reverse
  // complement of the other.
  // Of these, only one pair is what we want to call a duplicate. However our
  // filter cannot distinguish between the following cases
  //   1. a second mate pair that is identical to a previous mate pair
  //   2. Pair1_read1 is identical to Pair2_read1, and Pair2_read2 is identical
  //      to Pair3_read1.
  // the final mate pair in our file has one mate which is the same as a
  // left_mate in one read pair, and a right mate which is the same as the right
  // mate in another read pair. So it isn't really a dup, but we discard it
  // anyway as we cannot tell the difference.

  db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);

  // Reset counters
  file_pairs_loaded = 0;
  bad_reads = 0;
  dup_reads = 0;
  seq_read = 0;
  seq_loaded = 0;

  load_pe_filelists_into_graph_colour(
    "../data/test/graph/paired_end_file3_with_dups_1.fqlist",
    "../data/test/graph/paired_end_file3_with_dups_2.fqlist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_pe, ascii_fq_offset,
    into_colour, db_graph, 0,
    &file_pairs_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  // Five 36bp reads, left and right
  CU_ASSERT(seq_read == 360);
  CU_ASSERT(seq_loaded == 216);
  CU_ASSERT(dup_reads == 4);
  CU_ASSERT(file_pairs_loaded == 1);

  hash_table_free(&db_graph);
}



void test_loading_of_single_ended_reads_removing_duplicates()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // First set up the hash/graph
  int kmer_size = 21;
  int number_of_bits = 10;
  int bucket_size = 10;

  dBGraph * db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);

  // first a test where the file contains no duplicates and you do not try to
  // remove duplicates - does all the data get loaded?
  int fq_quality_cutoff = 1;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_read = 0, seq_loaded = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/graph/paired_end_file1_1.fqlist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0,
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  // Just lose 1 read plus part of another due to Ns
  CU_ASSERT(dup_reads == 0);
  CU_ASSERT(seq_read == 360);
  CU_ASSERT(seq_loaded == 311);
  CU_ASSERT(db_graph->unique_kmers == 105);

  hash_table_free(&db_graph);

  // then a test where you try to remove duplicates from files that don't contain any

  db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);

  remove_duplicates_se = true;
  seq_read = 0;
  seq_loaded = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/graph/paired_end_file1_1.fqlist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0,
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT(dup_reads == 0);
  CU_ASSERT(seq_read == 360);
  CU_ASSERT(seq_loaded == 311);
  CU_ASSERT(db_graph->unique_kmers == 105);

  hash_table_free(&db_graph);

  // then a test where you try to remove duplicates from files that do contain dupes

  db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);

  seq_read = 0;
  seq_loaded = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/graph/paired_end_file2_with_dup_1.fqlist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0,
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT(dup_reads == 2);
  // As before, but with four more 36bp reads
  CU_ASSERT(seq_read == 432);
  // Basically same as previous fastq - the extra reads are ignored
  CU_ASSERT(seq_loaded == 311);

  hash_table_free(&db_graph);
}

void test_load_seq_into_array()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  int length_of_arrays = 120;
  boolean expecting_new_fasta_entry = true;
  int max_kmer_size_used_in_this_test = 31;

  //============================================================================
  // Example 1. Load fasta containing one read, GGGG into the graph. Then use
  // same fasta and try to pull out the array of nodes corresponding to that
  // path through the graph

  // First set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size = 4;
  int max_retries = 10;

  dBGraph* db_graph = hash_table_new(number_of_bits, bucket_size,
                                     max_retries, kmer_size);

  if(db_graph == NULL)
  {
    die("Unable to alloc the hash table.");
  }

  // Read sequence
  int fq_quality_cutoff = 20;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_loaded = 0, seq_read = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/simplepeople_onlyperson1.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 1, // 1 => colourlist not falist/fqlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  // yes, only 4 bases in this fasta!
  CU_ASSERT(seq_read == 4);
  CU_ASSERT(seq_loaded == 4);


  dBNode**     path_nodes        = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays);
  Orientation* path_orientations = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays);
  Nucleotide*  path_labels       = (Nucleotide*) malloc(sizeof(Nucleotide)*length_of_arrays);
  char*        path_string       = (char*) malloc(sizeof(char)*(length_of_arrays+ max_kmer_size_used_in_this_test +1) ); //+1 for \0

  if ( (path_nodes==NULL) || (path_orientations==NULL) || (path_labels==NULL) || (path_string==NULL))
  {
    die("Failed to allocate arrays for test");
  }

  int i=0;

  //initialise
  for (i=0; i<length_of_arrays; i++)
  {
    path_nodes[i]=NULL;
    path_orientations[i]=forward;
    path_labels[i]=Undefined;
  }
  for (i=0; i< length_of_arrays+max_kmer_size_used_in_this_test; i++)
  {
    path_string[i]='Z';
  }
  path_string[length_of_arrays+max_kmer_size_used_in_this_test]='\0';

  char tmp_seq[5000];


  int num_of_nodes_to_read=2;

  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence\n");
  }
  alloc_sequence(seq, num_of_nodes_to_read+db_graph->kmer_size, 10);

  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*1000);
  kmer_window->nkmers=0;

  FILE* fptr = fopen("../data/test/pop_graph/simple1.fa", "r");

  if (fptr==NULL)
  {
    die("Error: Cannot open ../data/test/pop_graph/simple1.fa");
  }

  int retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string,
    seq, kmer_window, expecting_new_fasta_entry, db_graph);


  int offset=length_of_arrays-num_of_nodes_to_read;


  CU_ASSERT(retvalue==2);
  CU_ASSERT_STRING_EQUAL("CCC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset]==reverse);
  CU_ASSERT_STRING_EQUAL("CCC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+1]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+1]==reverse);
  CU_ASSERT_STRING_EQUAL("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZG", path_string);
  CU_ASSERT(path_labels[offset]==Guanine);
  CU_ASSERT(path_labels[offset+1]==Undefined);


  //we should now hit the end of the file, but it should not affect anything in our arrays
  CU_ASSERT(load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph)==0);
  CU_ASSERT_STRING_EQUAL("CCC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset]==reverse);
  CU_ASSERT_STRING_EQUAL("CCC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+1]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+1]==reverse);
  CU_ASSERT_STRING_EQUAL("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZG", path_string);
  CU_ASSERT(path_labels[offset]==Guanine);
  CU_ASSERT(path_labels[offset+1]==Undefined);


  fclose(fptr);


  //============================================================================
  // Example 2.
  // Fasta file containing GGNG. We expect to get back two nodes, both
  // containing the dummy kmer TTT, signifying that they contain an N


  // deliberately do not clean up db_graph. Should make no difference what other
  // paths/nodes/edges there are in th graph, so long as the fasta whose path we
  // are following has itself been loaded - we are only ever following one path.

  seq_read = 0;
  seq_loaded = 0;
  
  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/simplepeople_onlyperson2.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 1, // 1 => colourlist not falist/fqlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT(seq_read==4);
  CU_ASSERT(seq_loaded==0);

  fptr = fopen("../data/test/pop_graph/simple2.fa", "r");

  if (fptr==NULL)
  {
    die("Cannot open ../data/test/pop_graph/simple2.fa");
  }

  num_of_nodes_to_read = 2;
  expecting_new_fasta_entry = true;
  offset=length_of_arrays-num_of_nodes_to_read; // 120 - 2

  // Initialise
  for(i = 0; i < length_of_arrays; i++)
  {
    path_nodes[i]=NULL;
    path_orientations[i]=forward;
    path_labels[i]=Undefined;
  }
  for (i=0; i< length_of_arrays+max_kmer_size_used_in_this_test; i++)
  {
    path_string[i]='Z';
  }
  path_string[length_of_arrays+max_kmer_size_used_in_this_test]='\0';



  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);


  CU_ASSERT(retvalue==2);
  CU_ASSERT(path_nodes[offset]==NULL);
  CU_ASSERT(path_orientations[offset]==forward);
  CU_ASSERT(path_labels[offset]==Undefined);
  CU_ASSERT(path_string[offset]=='G');
  CU_ASSERT_STRING_EQUAL("G", path_string+offset);


  CU_ASSERT(path_nodes[offset+1]==NULL);
  CU_ASSERT(path_orientations[offset+1]==forward);
  CU_ASSERT(path_labels[offset+1]==Undefined);


  //we should now hit the end of the file, but it should not affect anything in our arrays
  CU_ASSERT(load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph)
   ==0);
  CU_ASSERT(path_nodes[offset]==NULL);
  CU_ASSERT(path_orientations[offset]==forward);
  CU_ASSERT(path_labels[offset]==Undefined);
  CU_ASSERT(path_string[offset]=='G');

  CU_ASSERT(path_nodes[offset+1]==NULL);
  CU_ASSERT(path_orientations[offset+1]==forward);
  CU_ASSERT(path_labels[offset+1]==Undefined);

  fclose(fptr);


  //============================================================================
  // Example 3: read contains only GNGNGNGNGNG

  seq_read = 0;
  seq_loaded = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/simplepeople_onlyperson3.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 1, // 1 => colourlist not falist/fqlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT(seq_read == 11);
  CU_ASSERT(seq_loaded == 0);

  fptr = fopen("../data/test/pop_graph/simple3.fa", "r");

  if (fptr==NULL)
  {
    die("Cannot open ../data/test/pop_graph/simple3.fa");
  }

  num_of_nodes_to_read = 9;

  //need more space in seq
  free_sequence(&seq);
  seq = malloc(sizeof(Sequence));

  if(seq == NULL)
  {
    die("Out of memory trying to allocate Sequence\n");
  }

  alloc_sequence(seq, num_of_nodes_to_read+db_graph->kmer_size, 10);

  expecting_new_fasta_entry=true;
  offset = length_of_arrays-num_of_nodes_to_read; //- db_graph->kmer_size+1


  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays,
    path_nodes, path_orientations, path_labels,
    path_string, seq, kmer_window,
    expecting_new_fasta_entry, db_graph);

  CU_ASSERT(retvalue==num_of_nodes_to_read);

  for (i=0; i< num_of_nodes_to_read; i++)
  {
    CU_ASSERT(path_nodes[offset+i]==NULL);
    CU_ASSERT(path_orientations[offset+i]==forward);
    CU_ASSERT(path_labels[offset+i]==Undefined);
  }


  CU_ASSERT(path_string[offset]   == 'N');
  CU_ASSERT(path_string[offset+1] == 'G');
  CU_ASSERT(path_string[offset+2] == 'N');
  CU_ASSERT(path_string[offset+3] == 'G');
  CU_ASSERT(path_string[offset+4] == 'N');
  CU_ASSERT(path_string[offset+5] == 'G');
  CU_ASSERT(path_string[offset+6] == 'N');
  CU_ASSERT(path_string[offset+7] == 'G');
  CU_ASSERT(path_string[offset+8]=='\0');

  CU_ASSERT(load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry,
    db_graph)==0);

  for (i=0; i< num_of_nodes_to_read; i++)
  {

    CU_ASSERT(path_nodes[offset+i]==NULL);
    CU_ASSERT(path_orientations[offset+i]==forward);
    CU_ASSERT(path_labels[offset+i]==Undefined);
  }

  CU_ASSERT(path_string[offset]   == 'N');
  CU_ASSERT(path_string[offset+1] == 'G');
  CU_ASSERT(path_string[offset+2] == 'N');
  CU_ASSERT(path_string[offset+3] == 'G');
  CU_ASSERT(path_string[offset+4] == 'N');
  CU_ASSERT(path_string[offset+5] == 'G');
  CU_ASSERT(path_string[offset+6] == 'N');
  CU_ASSERT(path_string[offset+7] == 'G');
  CU_ASSERT(path_string[offset+8]=='\0');

  fclose(fptr);


  //============================================================================
  // Example 4: fasta file contains: ACGCGCGTTTACG

  seq_read = 0;
  seq_loaded = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/simplepeople_onlyperson4.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 1, // 1 => colourlist not falist/fqlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT(seq_loaded==13);
  CU_ASSERT(seq_read==13);

  fptr = fopen("../data/test/pop_graph/simple4.fa", "r");
  if (fptr==NULL)
  {
    die("Cannot open ../data/test/pop_graph/simple4.fa\n");
  }


  num_of_nodes_to_read=11;

  //need more space in seq
  free_sequence(&seq);
  seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence\n");
  }
  alloc_sequence(seq, num_of_nodes_to_read+db_graph->kmer_size, 10);

  free(path_nodes);
  free(path_orientations);
  free(path_labels);
  free(path_string);

  path_nodes        = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays);
  path_orientations = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays);
  path_labels       = (Nucleotide*) malloc(sizeof(Nucleotide)*length_of_arrays);
  path_string       = (char*) malloc(sizeof(char)*(length_of_arrays+ max_kmer_size_used_in_this_test +1) ); //+1 for \0

  if ( (path_nodes==NULL) || (path_orientations==NULL) || (path_labels==NULL) || (path_string==NULL))
  {
    die("Failed to allocate arrays for test\n");
  }

  //initialise
  for (i=0; i<length_of_arrays; i++)
  {
    path_nodes[i]=NULL;
    path_orientations[i]=forward;
    path_labels[i]=Undefined;
  }
  for (i=0; i< length_of_arrays+max_kmer_size_used_in_this_test; i++)
  {
    path_string[i]='Z';
  }
  path_string[length_of_arrays+max_kmer_size_used_in_this_test]='\0';


  expecting_new_fasta_entry=true;
  offset=length_of_arrays-num_of_nodes_to_read;


  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);

  CU_ASSERT(retvalue==num_of_nodes_to_read);


  // here is where we test the path_string
  CU_ASSERT_STRING_EQUAL("CGCGTTTACG", path_string+offset);



  CU_ASSERT(path_nodes[offset]!=NULL);
  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[offset]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset]==forward);


  CU_ASSERT_STRING_EQUAL("CGC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+1]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+1]==forward);
  CU_ASSERT(path_labels[offset]==Cytosine);

  CU_ASSERT_STRING_EQUAL("CGC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+2]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+2]==reverse);
  CU_ASSERT(path_labels[offset+1]==Guanine);

  CU_ASSERT_STRING_EQUAL("CGC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+3]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+3]==forward);
  CU_ASSERT(path_labels[offset+2]==Cytosine);

  CU_ASSERT_STRING_EQUAL("CGC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+4]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+4]==reverse);
  CU_ASSERT(path_labels[offset+3]==Guanine);

  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+5]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+5]==reverse);
  CU_ASSERT(path_labels[offset+4]==Thymine);

  CU_ASSERT_STRING_EQUAL("AAC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+6]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+6]==reverse);
  CU_ASSERT(path_labels[offset+5]==Thymine);

  //note that path_nodes[offset+7] is NOT NULL. even though it encoded a TTT. This is because we check for 1's in ALL NUMBER_OF_BITFIELDS_IN_BINARY_KMER*64 bits of the BinaryKmer,
  //whereas when you read a TTT you only get 1's in the first 6 bits.
  CU_ASSERT_STRING_EQUAL("AAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+7]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+7]==reverse);
  CU_ASSERT(path_labels[offset+6]==Thymine);

  CU_ASSERT_STRING_EQUAL("TAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+8]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+8]==reverse);
  CU_ASSERT(path_labels[offset+7]==Adenine );

  CU_ASSERT_STRING_EQUAL("GTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+9]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+9]==reverse);
  CU_ASSERT(path_labels[offset+8]==Cytosine);

  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+10]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+10]==forward);
  CU_ASSERT(path_labels[offset+9]==Guanine);


  //now reach the end of file
  CU_ASSERT(load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry,
    db_graph)==0);
  //arrays should be unaffected


  CU_ASSERT_STRING_EQUAL("CGCGTTTACG", path_string+offset);

  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[offset]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset]==forward);

  CU_ASSERT_STRING_EQUAL("CGC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+1]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+1]==forward);
  CU_ASSERT(path_labels[offset]==Cytosine);

  CU_ASSERT_STRING_EQUAL("CGC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+2]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+2]==reverse);
  CU_ASSERT(path_labels[offset+1]==Guanine);

  CU_ASSERT_STRING_EQUAL("CGC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+3]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+3]==forward);
  CU_ASSERT(path_labels[offset+2]==Cytosine);

  CU_ASSERT_STRING_EQUAL("CGC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+4]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+4]==reverse);
  CU_ASSERT(path_labels[offset+3]==Guanine);

  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+5]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+5]==reverse);
  CU_ASSERT(path_labels[offset+4]==Thymine);

  CU_ASSERT_STRING_EQUAL("AAC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+6]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+6]==reverse);
  CU_ASSERT(path_labels[offset+5]==Thymine);

  //note that path_nodes[offset+7] is NOT NULL. even though it encoded a TTT. This is because we check for 1's in ALL 64 bits of the kmer,
  //whereas when you read a TTT you only get 1's in the first 6 bits.
  CU_ASSERT_STRING_EQUAL("AAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+7]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+7]==reverse);
  CU_ASSERT(path_labels[offset+6]==Thymine);

  CU_ASSERT_STRING_EQUAL("TAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+8]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+8]==reverse);
  CU_ASSERT(path_labels[offset+7]==Adenine );

  CU_ASSERT_STRING_EQUAL("GTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+9]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+9]==reverse);
  CU_ASSERT(path_labels[offset+8]==Cytosine);

  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+10]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+10]==forward);
  CU_ASSERT(path_labels[offset+9]==Guanine);

  fclose(fptr);
  free_sequence(&seq);


  //============================================================================
  // Example 5: Fasta contains a chunk of chromosome 1

  // >from chromosome 1
  // AACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAAC
  // CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAAACCCTAAACCCTAACCCTAAC
  // CCTAACCCTAACCCTAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCTAACCCCTAACCCTAACCC
  // TAACCCTACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCCTAACCCTAACCCTAACCCTA
  // ACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTCGCGGTACCCTC

  // no N's in this

  // Cleanup, as we now want a new graph
  hash_table_free(&db_graph);


  kmer_size = 31;
  number_of_bits = 7;
  bucket_size    = 5;
  bad_reads = 0;
  max_retries=10;

  db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  if (db_graph==NULL)
  {
    die("Unable to alloc the hash table.\n");
  }

  seq_read = 0;
  seq_loaded = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/simplepeople_onlyperson5.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 1, // 1 => colourlist not falist/fqlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT(seq_loaded==360);
  CU_ASSERT(seq_read==360);

  fptr = fopen("../data/test/pop_graph/simple5.fa", "r");

  if (fptr==NULL)
  {
    die("Cannot open ../data/test/pop_graph/simple5.fa\n");
  }


  num_of_nodes_to_read=330;
  length_of_arrays=400;


  seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence\n");
  }
  alloc_sequence(seq, num_of_nodes_to_read+db_graph->kmer_size, 10);

  free(path_nodes);
  free(path_orientations);
  free(path_labels);
  free(path_string);
  path_nodes       = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays);
  path_orientations = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays);
  path_labels         = (Nucleotide*) malloc(sizeof(Nucleotide)*length_of_arrays);
  path_string            = (char*) malloc(sizeof(char)*(length_of_arrays+ max_kmer_size_used_in_this_test +1) ); //+1 for \0

  if ( (path_nodes==NULL) || (path_orientations==NULL) || (path_labels==NULL) || (path_string==NULL))
  {
    die("Failed to allocate arrays for test\n");
  }

  //initialise
  for (i=0; i<length_of_arrays; i++)
  {
    path_nodes[i]=NULL;
    path_orientations[i]=forward;
    path_labels[i]=Undefined;
  }
  for (i=0; i< length_of_arrays+max_kmer_size_used_in_this_test; i++)
  {
    path_string[i]='Z';
  }
  path_string[length_of_arrays+max_kmer_size_used_in_this_test]='\0';


  expecting_new_fasta_entry=true;
  offset=length_of_arrays-num_of_nodes_to_read;


  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);

  CU_ASSERT(retvalue==num_of_nodes_to_read);

  //path_string should be the sequence loaded, but not the first kmer - ie matches the edges
  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCTAACCCCTAACCCTAACCCTAACCCTACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTCGCGGTACCCTC", path_string+offset);


  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCCCT", binary_kmer_to_seq(element_get_kmer(path_nodes[offset]), db_graph->kmer_size, tmp_seq));

  CU_ASSERT_STRING_EQUAL("ACCCTAACCCTAACCCTAACCCTAACCCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+1]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+1]==forward);
  CU_ASSERT(path_labels[offset]==Adenine);

  CU_ASSERT_STRING_EQUAL("CCCTAACCCTAACCCTAACCCTAACCCCTAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+2]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+2]==forward);
  CU_ASSERT(path_labels[offset+1]==Adenine);


  CU_ASSERT_STRING_EQUAL("CCTAACCCTAACCCTAACCCTAACCCCTAAC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+3]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+3]==forward);
  CU_ASSERT(path_labels[offset+2]==Cytosine);


  CU_ASSERT_STRING_EQUAL("CTAACCCTAACCCTAACCCTAACCCCTAACC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+4]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+4]==forward);
  CU_ASSERT(path_labels[offset+3]==Cytosine);

  CU_ASSERT_STRING_EQUAL("GGGTTAGGGGTTAGGGTTAGGGTTAGGGTTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+5]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+5]==reverse);
  CU_ASSERT(path_labels[offset+4]==Cytosine);

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCCTAACCCT", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+6]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+6]==forward);
  CU_ASSERT(path_labels[offset+5]==Thymine);

  CU_ASSERT_STRING_EQUAL("ACCCTAACCCTAACCCTAACCCCTAACCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+7]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+7]==forward);
  CU_ASSERT(path_labels[offset+6]==Adenine);

  CU_ASSERT_STRING_EQUAL("CCCTAACCCTAACCCTAACCCCTAACCCTAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+8]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+8]==forward);
  CU_ASSERT(path_labels[offset+7]==Adenine);

  CU_ASSERT_STRING_EQUAL("CCTAACCCTAACCCTAACCCCTAACCCTAAC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+9]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+9]==forward);
  CU_ASSERT(path_labels[offset+8]==Cytosine);

  //and one node near the end
  CU_ASSERT_STRING_EQUAL("CCTAACCCTAACCCCTAACCCTAACCCTAAC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+308]),db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+308]==forward);
  CU_ASSERT(path_labels[offset+307]==Cytosine);

  fclose(fptr);


  //============================================================================
  // Example 6. Fasta contains a small chunk of chromosome 1 with an N inserted

  // >r
  // AACCCTAACCCTAACCCTAACCCTAACCCCTAACCNTAACCCTAACCCTAACCCTAACCTAACCCTAA

  // Deliberately do not cleanup, as we do not want a new graph - we want to be
  // sure we can follow our path without other stuff in graqph distracting

  seq_read = 0;
  seq_loaded = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/simplepeople_onlyperson6.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 1, // 1 => colourlist not falist/fqlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT(seq_read==68);
  CU_ASSERT(seq_loaded==67);//one N

  fptr = fopen("../data/test/pop_graph/simple6.fa", "r");

  if (fptr==NULL)
  {
    die("Cannot open ../data/test/pop_graph/simple6.fa\n");
  }


  num_of_nodes_to_read=38;//num of bases is 68
  expecting_new_fasta_entry=true;
  offset=length_of_arrays-num_of_nodes_to_read;

  //initialise
  for (i=0; i<length_of_arrays; i++)
  {
    path_nodes[i]=NULL;
    path_orientations[i]=forward;
    path_labels[i]=Undefined;
  }
  for (i=0; i< length_of_arrays+max_kmer_size_used_in_this_test; i++)
  {
    path_string[i]='Z';
  }
  path_string[length_of_arrays+max_kmer_size_used_in_this_test]='\0';

  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);

  CU_ASSERT(retvalue==num_of_nodes_to_read);

  CU_ASSERT_STRING_EQUAL("AACCNTAACCCTAACCCTAACCCTAACCTAACCCTAA", path_string+offset);

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCCCT", binary_kmer_to_seq(element_get_kmer(path_nodes[offset]), db_graph->kmer_size, tmp_seq));

  CU_ASSERT_STRING_EQUAL("ACCCTAACCCTAACCCTAACCCTAACCCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+1]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+1]==forward);
  CU_ASSERT(path_labels[offset]==Adenine);

  CU_ASSERT_STRING_EQUAL("CCCTAACCCTAACCCTAACCCTAACCCCTAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+2]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+2]==forward);
  CU_ASSERT(path_labels[offset+1]==Adenine);

  CU_ASSERT_STRING_EQUAL("CCTAACCCTAACCCTAACCCTAACCCCTAAC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+3]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+3]==forward);
  CU_ASSERT(path_labels[offset+2]==Cytosine);

  CU_ASSERT_STRING_EQUAL("CTAACCCTAACCCTAACCCTAACCCCTAACC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+4]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+4]==forward);
  CU_ASSERT(path_labels[offset+3]==Cytosine);


  for (i=0; i<31; i++)
  {
    CU_ASSERT(path_nodes[offset+5+i]==NULL);
  }

  CU_ASSERT_STRING_EQUAL("TAACCCTAACCCTAACCCTAACCTAACCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+36]),  db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+36]==forward);
  CU_ASSERT(path_labels[offset+35]==Undefined);  // This is important!!!

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCTAACCCTAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+37]),  db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+37]==forward);
  CU_ASSERT(path_labels[offset+36]==Adenine);



  //============================================================================
  // Example 7. Fasta contains a small chunk of chromosome 1 with several
  // consecutive N's inserted

  // >zammo
  // AACCCTAACCCTAACCCTAACCCTAACCCCTAACCNNNNNTAACCCTAACCCTAACCCTAACCTAACCCTAA


  //Deliberately do not cleanup, as we do not want a new graph - we want to be sure we can follow our path without other stuff in graqph distracting

  seq_read = 0;
  seq_loaded = 0;
  
  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/simplepeople_onlyperson7.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 1, // 1 => colourlist not falist/fqlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);
  
  CU_ASSERT(seq_read==72);
  CU_ASSERT(seq_loaded==67);//5 N's

  fptr = fopen("../data/test/pop_graph/simple7.fa", "r");
  if (fptr==NULL)
  {
    die("Cannot open ../data/test/pop_graph/simple7.fa\n");
  }

  num_of_nodes_to_read=42;//number of bases is 72
  expecting_new_fasta_entry=true;
  offset=length_of_arrays-num_of_nodes_to_read;

  //initialise
  for (i=0; i<length_of_arrays; i++)
  {
    path_nodes[i]=NULL;
    path_orientations[i]=forward;
    path_labels[i]=Undefined;
  }
  for (i=0; i< length_of_arrays+max_kmer_size_used_in_this_test; i++)
  {
    path_string[i]='Z';
  }
  path_string[length_of_arrays+max_kmer_size_used_in_this_test]='\0';


  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays,
   path_nodes, path_orientations, path_labels,
   path_string, seq, kmer_window,
   expecting_new_fasta_entry, db_graph);


  CU_ASSERT(retvalue==num_of_nodes_to_read);

  CU_ASSERT_STRING_EQUAL("AACCNNNNNTAACCCTAACCCTAACCCTAACCTAACCCTAA", path_string+offset);

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCCCT", binary_kmer_to_seq(element_get_kmer(path_nodes[offset]), db_graph->kmer_size, tmp_seq));

  CU_ASSERT_STRING_EQUAL("ACCCTAACCCTAACCCTAACCCTAACCCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+1]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+1]==forward);
  CU_ASSERT(path_labels[offset]==Adenine);

  CU_ASSERT_STRING_EQUAL("CCCTAACCCTAACCCTAACCCTAACCCCTAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+2]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+2]==forward);
  CU_ASSERT(path_labels[offset+1]==Adenine);

  CU_ASSERT_STRING_EQUAL("CCTAACCCTAACCCTAACCCTAACCCCTAAC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+3]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+3]==forward);
  CU_ASSERT(path_labels[offset+2]==Cytosine);

  CU_ASSERT_STRING_EQUAL("CTAACCCTAACCCTAACCCTAACCCCTAACC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+4]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+4]==forward);
  CU_ASSERT(path_labels[offset+3]==Cytosine);

  for (i=0; i<35; i++)
  {
    CU_ASSERT(path_nodes[offset+5+i]==NULL);
  }

  CU_ASSERT_STRING_EQUAL("TAACCCTAACCCTAACCCTAACCTAACCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+40]),  db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+40]==forward);
  CU_ASSERT(path_labels[offset+39]==Undefined);  // this is important!!!

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCTAACCCTAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+41]),  db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+41]==forward);
  CU_ASSERT(path_labels[offset+40]==Adenine);


  fclose(fptr);





  //============================================================================
  // Example 8: Fasta entry starts with Ns

  // >z
  // NNNNAAACGT



  //cleanup, as we now want a new graph - changing kmer size
  hash_table_free(&db_graph);


  kmer_size = 5;
  number_of_bits = 7;
  bucket_size    = 5;
  bad_reads = 0;
  max_retries=10;

  db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  if (db_graph==NULL)
  {
    die("Unable to alloc the hash table.\n");
  }


  num_of_nodes_to_read=6;
  length_of_arrays=40;
  offset=length_of_arrays-num_of_nodes_to_read;
  expecting_new_fasta_entry=true;

  free_sequence(&seq);
  seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence\n");
  }
  alloc_sequence(seq, num_of_nodes_to_read+db_graph->kmer_size, 10);

  free(path_nodes);
  free(path_orientations);
  free(path_labels);
  free(path_string);
  path_nodes       = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays);
  path_orientations = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays);
  path_labels         = (Nucleotide*) malloc(sizeof(Nucleotide)*length_of_arrays);
  path_string            = (char*) malloc(sizeof(char)*(length_of_arrays+ max_kmer_size_used_in_this_test +1) ); //+1 for \0

  if ( (path_nodes==NULL) || (path_orientations==NULL) || (path_labels==NULL) || (path_string==NULL))
  {
    die("%s:%i: Failed to allocate arrays for tests", __FILE__, __LINE__);
  }


  // initialise
  for (i=0; i<length_of_arrays; i++)
  {
    path_nodes[i]=NULL;
    path_orientations[i]=forward;
    path_labels[i]=Undefined;
  }
  for (i=0; i< length_of_arrays+max_kmer_size_used_in_this_test; i++)
  {
    path_string[i]='Z';
  }
  path_string[length_of_arrays+max_kmer_size_used_in_this_test]='\0';

  seq_read = 0;
  seq_loaded = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/simplepeople_onlyperson8.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 1, // 1 => colourlist not falist/fqlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT(seq_read==10);
  CU_ASSERT(seq_loaded==6);

  fptr = fopen("../data/test/pop_graph/simple8.fa", "r");

  if (fptr==NULL)
  {
    die("Cannot open ../data/test/pop_graph/simple8.fa\n");
  }

  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);

  CU_ASSERT(retvalue==num_of_nodes_to_read);

  CU_ASSERT_STRING_EQUAL("AACGT", path_string+offset);

  CU_ASSERT(path_nodes[offset]==NULL);
  CU_ASSERT(path_nodes[offset+1]==NULL);
  CU_ASSERT(path_nodes[offset+2]==NULL);
  CU_ASSERT(path_nodes[offset+3]==NULL);


  CU_ASSERT_STRING_EQUAL("AAACG", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+4]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+4]==forward);
  CU_ASSERT(path_labels[offset+3]==Undefined);

  CU_ASSERT_STRING_EQUAL("AACGT", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+5]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+5]==forward);
  CU_ASSERT(path_labels[offset+4]==Thymine);

  fclose(fptr);


  //============================================================================
  // Example 9:

  // >zam
  // AAAANAANAAANAAAANCGTCT

  num_of_nodes_to_read=18;//22 bases
  length_of_arrays=40;
  offset=length_of_arrays-num_of_nodes_to_read;
  expecting_new_fasta_entry=true;

  free_sequence(&seq);
  seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence\n");
  }
  alloc_sequence(seq, num_of_nodes_to_read+db_graph->kmer_size, 10);

  seq_read = 0;
  seq_loaded = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/simplepeople_onlyperson9.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 1, // 1 => colourlist not falist/fqlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT(seq_read == 22);
  CU_ASSERT(seq_loaded == 5); // kmer is 5 and there are many Ns

  fptr = fopen("../data/test/pop_graph/simple9.fa", "r");
  if (fptr==NULL)
  {
    die("Cannot open ../data/test/pop_graph/simple9.fa\n");
  }

  //initialise
  for (i=0; i<length_of_arrays; i++)
  {
    path_nodes[i]=NULL;
    path_orientations[i]=forward;
    path_labels[i]=Undefined;
  }
  for (i=0; i< length_of_arrays+max_kmer_size_used_in_this_test; i++)
  {
    path_string[i]='Z';
  }
  path_string[length_of_arrays+max_kmer_size_used_in_this_test]='\0';


  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);


  CU_ASSERT(retvalue==num_of_nodes_to_read);
  CU_ASSERT_STRING_EQUAL("AANAAANAAAANCGTCT", path_string+offset);

  CU_ASSERT(path_nodes[offset]==NULL);
  CU_ASSERT(path_nodes[offset+1]==NULL);
  CU_ASSERT(path_nodes[offset+2]==NULL);
  CU_ASSERT(path_nodes[offset+3]==NULL);
  CU_ASSERT(path_nodes[offset+4]==NULL);
  CU_ASSERT(path_nodes[offset+5]==NULL);
  CU_ASSERT(path_nodes[offset+6]==NULL);
  CU_ASSERT(path_nodes[offset+7]==NULL);
  CU_ASSERT(path_nodes[offset+8]==NULL);
  CU_ASSERT(path_nodes[offset+9]==NULL);
  CU_ASSERT(path_nodes[offset+10]==NULL);
  CU_ASSERT(path_nodes[offset+11]==NULL);
  CU_ASSERT(path_nodes[offset+12]==NULL);
  CU_ASSERT(path_nodes[offset+13]==NULL);
  CU_ASSERT(path_nodes[offset+14]==NULL);
  CU_ASSERT(path_nodes[offset+15]==NULL);
  CU_ASSERT(path_nodes[offset+16]==NULL);


  CU_ASSERT_STRING_EQUAL("AGACG", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+17]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+4]==forward);
  CU_ASSERT(path_labels[offset+3]==Undefined);

  free_sequence(&seq);
  free(path_nodes);
  free(path_orientations);
  free(path_labels);
  free(path_string);
  hash_table_free(&db_graph);

  fclose(fptr);


  // So far we have tested that when we call this function once, we get the
  // right nodes/Nucleotides etc. What remains is to check that we can call it
  // in a loop and it will successfully keep reading the next chunk of the file,
  // and putting it in the right place in the array. Remember we are doing this:

  // array (........, X,X,X,X,...X)        The X's are the result of a previous
  // call. Say there are 100 of them then push the array 50 places to the left
  // and call load_seq_into_array, loading another 50 into the right-hand-end
  // Does everything go into the right index of the array? Do we get the right
  // nodes at the start of the second batch? At that point we will read an extra
  // base, but need the previous k-1 to construct the binary kmer. What happens
  // if we had N's at the end of the previous batch?

  // ===========================================================================
  // Example 10 - simplest example of loading batches. We will have an array of
  // length 8, and load 4 bases at a time.

  // >zam
  // ACGTACGTACGTACGT


  kmer_size = 3;
  number_of_bits = 7;
  bucket_size    = 5;
  bad_reads = 0;
  max_retries=10;

  db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  if (db_graph==NULL)
  {
    die("Unable to alloc the hash table.\n");
  }


  num_of_nodes_to_read=2;//let's get the first 4 bases first
  length_of_arrays=8;

  seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence\n");
  }
  alloc_sequence(seq, length_of_arrays+db_graph->kmer_size, 10);


  path_nodes       = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays);
  path_orientations = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays);
  path_labels         = (Nucleotide*) malloc(sizeof(Nucleotide)*length_of_arrays);
  path_string            = (char*) malloc(sizeof(char)*(length_of_arrays+ max_kmer_size_used_in_this_test +1) ); //+1 for \0

  if ( (path_nodes==NULL) || (path_orientations==NULL) || (path_labels==NULL) || (path_string==NULL))
  {
    die("Failed to allocate arrays for test\n");
  }

  //initialise
  for (i=0; i<length_of_arrays; i++)
  {
    path_nodes[i]=NULL;
    path_orientations[i]=forward;
    path_labels[i]=Undefined;
  }
  for (i=0; i< length_of_arrays+max_kmer_size_used_in_this_test; i++)
  {
    path_string[i]='Z';
  }
  path_string[length_of_arrays+max_kmer_size_used_in_this_test]='\0';


  seq_read = 0;
  seq_loaded = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/simplepeople_onlyperson10.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 1, // 1 => colourlist not falist/fqlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT(seq_read == 16);
  CU_ASSERT(seq_loaded == 16);

  fptr = fopen("../data/test/pop_graph/simple10.fa", "r");

  if (fptr==NULL)
  {
    die("Cannot open ../data/test/pop_graph/simple10.fa\n");
  }

  //initialise the arrays
  for (i=0; i<length_of_arrays; i++)
  {
    path_nodes[i]=NULL;
    path_orientations[i]=forward;
    path_labels[i]=Undefined;
  }

  //when we load the first chunk, we expect the start of an entry..
  expecting_new_fasta_entry=true;
  //remember offset is not passed into the function call, it just tells us where we expect the answers to be
  offset=length_of_arrays-num_of_nodes_to_read;

  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);

  CU_ASSERT(retvalue==num_of_nodes_to_read);

  for (i=0; i<offset; i++)
  {
    CU_ASSERT(path_nodes[i]==NULL);
  }

  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[offset]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset]==forward);
  CU_ASSERT_STRING_EQUAL("T", path_string+offset);

  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+1]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+1]==reverse);
  CU_ASSERT(path_labels[offset]==Thymine);


  //now we want to read another 4 nodes
  num_of_nodes_to_read=4;

  //Now move the data we have currently in our arrays along by Number_of_nodes_to_read
  for (i=0; i<length_of_arrays-num_of_nodes_to_read; i++)
  {
    path_nodes[i]       = path_nodes[i+num_of_nodes_to_read];
    path_orientations[i]= path_orientations[i+num_of_nodes_to_read];
    path_labels[i]      = path_labels[i+num_of_nodes_to_read];
    path_string[i]      = path_string[i+num_of_nodes_to_read];
  }

  //check we have moved it correctly! Arrays are length 8, and we have loaded 4 bases=2 kmers so far
  for (i=0; i<2; i++)
  {
    CU_ASSERT(path_nodes[i]==NULL);
    CU_ASSERT(path_orientations[i]==forward);
    CU_ASSERT(path_labels[i]==Undefined);
  }

  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[2]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT_STRING_EQUAL("T", path_string+offset-num_of_nodes_to_read);//offset

  CU_ASSERT(path_orientations[2]==forward);
  CU_ASSERT(path_labels[1]==Undefined);
  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[3]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[3]==reverse);
  CU_ASSERT(path_labels[2]==Thymine);

  //move final kmer to the start of seq.
  for (i=0; i<db_graph->kmer_size; i++)
  {
    seq->seq[i]=seq->seq[num_of_nodes_to_read-db_graph->kmer_size+i];
  }

  //and of course offset alters:
  offset=offset-num_of_nodes_to_read;

  //OK - now ready to load next batch, this time of 4 nodes, remembering that now we no longer expect a new fasta entry.
  expecting_new_fasta_entry=false;
  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);
  CU_ASSERT(retvalue==num_of_nodes_to_read);


  CU_ASSERT_STRING_EQUAL("TACGT", path_string+offset);


  //first four entries in array unchanged
  //check first 4 entries
  for (i=0; i<2; i++)
  {
    CU_ASSERT(path_nodes[i]==NULL);
    CU_ASSERT(path_orientations[i]==forward);
    CU_ASSERT(path_labels[i]==Undefined);
  }


  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[2]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[2]==forward);
  CU_ASSERT(path_labels[1]==Undefined);
  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[3]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[3]==reverse);
  CU_ASSERT(path_labels[2]==Thymine);

  //check last 4
  CU_ASSERT_STRING_EQUAL("GTA", binary_kmer_to_seq(element_get_kmer(path_nodes[4]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[4]==forward);
  CU_ASSERT(path_labels[3]==Adenine);

  CU_ASSERT_STRING_EQUAL("GTA", binary_kmer_to_seq(element_get_kmer(path_nodes[5]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[5]==reverse);
  CU_ASSERT(path_labels[4]==Cytosine);

  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[6]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[6]==forward);
  CU_ASSERT(path_labels[5]==Guanine);

  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[7]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[7]==reverse);
  CU_ASSERT(path_labels[6]==Thymine);


  free_sequence(&seq);
  free(path_nodes);
  free(path_orientations);
  free(path_labels);
  free(path_string);
  hash_table_free(&db_graph);


// Example 11: Harder case of loading batches, this time where batch size much
// bigger than kmer_size

// First 20 lines of chromosome 1, contains 1140 bases, but I want to not have
// to worry about that and just keep loading chunks until I run out.
// Kmer_size 31. Take an array of length 800, and load 400 nodes at a time.

//  >1 dna:chromosome chromosome:NCBI36:1:1:247249719:1
//  TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC
//  TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCT
//  AACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCT
//  AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCC
//  TAACCCTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACCCCAAC
//  CCCAACCCCAACCCCAACCCCAACCCTAACCCCTAACCCTAACCCTAACCCTACCCTAAC
//  CCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCCTAACCCTAACCCTAACCCTA
//  ACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTCGCGGTACCCTC
//  AGCCGGCCCGCCCGCCCGGGTCTGACCTGAGGAGAACTGTGCTCCGCCTTCAGAGTACCA
//  CCGAAATCTGTGCAGAGGACAACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCTGTGCT
//  GAGGAGAACGCAACTCCGCCGGCGCAGGCGCAGAGAGGCGCGCCGCGCCGGCGCAGGCGC
//  AGACACATGCTAGCGCGTCGGGGTGGAGGCGTGGCGCAGGCGCAGAGAGGCGCGCCGCGC
//  CGGCGCAGGCGCAGAGACACATGCTACCGCGTCCAGGGGTGGAGGCGTGGCGCAGGCGCA
//  GAGAGGCGCACCGCGCCGGCGCAGGCGCAGAGACACATGCTAGCGCGTCCAGGGGTGGAG
//  GCGTGGCGCAGGCGCAGAGACGCAAGCCTACGGGCGGGGGTTGGGGGGGCGTGTGTTGCA
//  GGAGCAAAGTCGCACGGCGCCGGGCTGGGGCGGGGGGAGGGTGGCGCCGTGCACGCGCAG
//  AAACTCACGTCACGGTGGCGCGGCGCAGAGACGGGTAGAACCTCAGTAATCCGAAAAGCC
//  GGGATCGACCGCCCCTTGCTTGCAGCCGGGCACTACAGGACCCGCTTGCTCACGGTGCTG
//  TGCCAGGGCGCCCCCTGCTGGCGACTAGGGCAACTGCAGGGCTCTCTTGCTTAGAGTGGT
//


  kmer_size = 31;
  number_of_bits = 12;
  bucket_size    = 20;
  bad_reads = 0;
  max_retries=10;

  db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  if (db_graph==NULL)
  {
    die("Unable to alloc the hash table. \n");
  }


  num_of_nodes_to_read=400;
  length_of_arrays=800;;

  seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence\n");
  }
  // 40 is max length of read name
  alloc_sequence(seq, length_of_arrays+db_graph->kmer_size, 40);


  path_nodes             = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays);
  path_orientations      = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays);
  path_labels            = (Nucleotide*) malloc(sizeof(Nucleotide)*length_of_arrays);
  path_string            = (char*) malloc(sizeof(char)*(length_of_arrays+ max_kmer_size_used_in_this_test +1) ); //+1 for \0

  if ( (path_nodes==NULL) || (path_orientations==NULL) || (path_labels==NULL) || (path_string==NULL))
  {
    die("%s:%i: Failed to allocate arrays for test", __FILE__, __LINE__);
  }

  //initialise
  for (i=0; i<length_of_arrays; i++)
  {
    path_nodes[i]=NULL;
    path_orientations[i]=forward;
    path_labels[i]=Undefined;
  }
  for (i=0; i< length_of_arrays+max_kmer_size_used_in_this_test; i++)
  {
    path_string[i]='Z';
  }
  path_string[length_of_arrays+max_kmer_size_used_in_this_test]='\0';

  seq_read = 0;
  seq_loaded = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/simplepeople_onlyperson11.colours",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 1, // 1 => colourlist not falist/fqlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT(seq_read == 1140);
  CU_ASSERT(seq_loaded == 1140);

  fptr = fopen("../data/test/pop_graph/simple11.fa", "r");
  if (fptr==NULL)
  {
    die("Cannot open ../data/test/pop_graph/simple10.fa\n");
  }

  // Initialise the arrays
  for (i=0; i<length_of_arrays; i++)
  {
    path_nodes[i]=NULL;
    path_orientations[i]=forward;
    path_labels[i]=Undefined;
    path_string[i]='N';
  }

  // When we load the first chunk, we expect the start of an entry..
  expecting_new_fasta_entry=true;

  // Remember offset is not passed into the function call, it just tells us
  // where we expect the answers to be
  offset=length_of_arrays-num_of_nodes_to_read;

  // == first batch load ==
  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);

  CU_ASSERT(retvalue==num_of_nodes_to_read);

  for (i=0; i<offset; i++)
  {
    CU_ASSERT(path_nodes[i]       ==NULL);
    CU_ASSERT(path_orientations[i]==forward);
    CU_ASSERT(path_labels[i]      ==Undefined);
    CU_ASSERT(path_string[i]      =='N');
  }


  //take first 430 bases of file = first 400 nodes, MINUS the first 31-mer
  CU_ASSERT_STRING_EQUAL( "AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCTAACCCCTAACCCTAACCCTAACCCTACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCC" ,path_string+offset);


  CU_ASSERT_STRING_EQUAL("AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset]), db_graph->kmer_size, tmp_seq)); //reverse complement of first kmer in fasta - TAACCCTAACCCTAACCCTAACCCTAACCCT
  CU_ASSERT(path_orientations[offset]==reverse);

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+1]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+1]==forward);
  CU_ASSERT(path_labels[offset]==Adenine);

  CU_ASSERT_STRING_EQUAL("ACCCTAACCCTAACCCTAACCCTAACCCTAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+2]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+2]==forward);
  CU_ASSERT(path_labels[offset+1]==Adenine);

  CU_ASSERT_STRING_EQUAL("CCCTAACCCTAACCCTAACCCTAACCCTAAC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+3]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+3]==forward);
  CU_ASSERT(path_labels[offset+2]==Cytosine);

  CU_ASSERT_STRING_EQUAL("AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+60]), db_graph->kmer_size, tmp_seq)); //first kmer on second line of fasta
  CU_ASSERT(path_orientations[offset+60]==reverse);
  CU_ASSERT(path_labels[offset+59]==Thymine);

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCCCT", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+120]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+120]==forward);
  CU_ASSERT(path_labels[offset+119]==Thymine);

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+180]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+180]==forward);
  CU_ASSERT(path_labels[offset+179]==Adenine);

  CU_ASSERT_STRING_EQUAL("GGTTAGGGTTAGGGTTTAGGGTTTAGGGTTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+240]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+240]==reverse);
  CU_ASSERT(path_labels[offset+239]==Cytosine);

  CU_ASSERT_STRING_EQUAL("CCTAACCCTAACCCTAACCCTAACCCTAACC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+360]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+360]==forward);
  CU_ASSERT(path_labels[offset+359]==Cytosine);


  CU_ASSERT_STRING_EQUAL("CTAACCCTAACCCTAACCCTAACCCTAACCC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+399]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+399]==forward);
  CU_ASSERT(path_labels[offset+398]==Cytosine);


  // First batch correctly loaded. Now load another 400, and ensure the
  // transition is OK. num_of_nodes_to_read is already = 400
  // From now on, each new base will be a new node.


  //push everything along in the arrays:
  for (i=0; i<length_of_arrays-num_of_nodes_to_read; i++)
  {
    path_nodes[i]       = path_nodes[i+num_of_nodes_to_read];
    path_orientations[i]= path_orientations[i+num_of_nodes_to_read];
    path_labels[i]      = path_labels[i+num_of_nodes_to_read];
    path_string[i]      = path_string[i+num_of_nodes_to_read];
  }


  //move final kmer to the start of seq.

  // Note there are two separate complications relating to how far we push seq
  // to the left
  // 1. When you start, to get n nodes, you need k + n-1 bases. After that, you
  //    only need n
  // 2. After the first batch, you always want to leave a kmer's worth of bases
  //    in seq before loading the next batch, to allow calculation of the edge
  //    between the last node of first batch, and first of next batch.

  // This time, you don't have an extra kmer of bases in seq from the last time
  // around.
  for (i=0; i<db_graph->kmer_size; i++)
  {
    // because num_nodes = num_bases-k+1; so num_bases=num_nodes+k-1;
    // so last k bases are the num_nodes,...,num_nodes+k-1 -th bases.
    // But we count from 0...
    seq->seq[i]=seq->seq[num_of_nodes_to_read-1+i];
  }

  seq->seq[db_graph->kmer_size]='\0';

  // Double check, compare with what used to be in path_nodes[offset+399] -
  // see above
  CU_ASSERT_STRING_EQUAL("CTAACCCTAACCCTAACCCTAACCCTAACCC", seq->seq);

  // Load in another 400
  expecting_new_fasta_entry=false;

  // == second batch load ==
  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);
  CU_ASSERT(retvalue==num_of_nodes_to_read);


  // First the obvious - these nodes that were between offset and offset+399
  // should now be between 0 and 399 this is just checking my test essentially

  CU_ASSERT_STRING_EQUAL("AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset]), db_graph->kmer_size, tmp_seq)); //reverse complement of first kmer in fasta - TAACCCTAACCCTAACCCTAACCCTAACCCT
  CU_ASSERT(path_orientations[0]==reverse);

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[1]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[1]==forward);
  CU_ASSERT(path_labels[1]==Adenine);

  CU_ASSERT_STRING_EQUAL("ACCCTAACCCTAACCCTAACCCTAACCCTAA", binary_kmer_to_seq(element_get_kmer(path_nodes[2]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[2]==forward);
  CU_ASSERT(path_labels[1]==Adenine);

  CU_ASSERT_STRING_EQUAL("CCCTAACCCTAACCCTAACCCTAACCCTAAC", binary_kmer_to_seq(element_get_kmer(path_nodes[3]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[3]==forward);
  CU_ASSERT(path_labels[2]==Cytosine);

  CU_ASSERT_STRING_EQUAL("AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA", binary_kmer_to_seq(element_get_kmer(path_nodes[60]), db_graph->kmer_size, tmp_seq)); //first kmer on second line of fasta
  CU_ASSERT(path_orientations[60]==reverse);
  CU_ASSERT(path_labels[59]==Thymine);

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCCCT", binary_kmer_to_seq(element_get_kmer(path_nodes[120]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[120]==forward);
  CU_ASSERT(path_labels[119]==Thymine);

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[180]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[180]==forward);
  CU_ASSERT(path_labels[179]==Adenine);

  CU_ASSERT_STRING_EQUAL("GGTTAGGGTTAGGGTTTAGGGTTTAGGGTTA", binary_kmer_to_seq(element_get_kmer(path_nodes[240]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[240]==reverse);
  CU_ASSERT(path_labels[239]==Cytosine);

  CU_ASSERT_STRING_EQUAL("CCTAACCCTAACCCTAACCCTAACCCTAACC", binary_kmer_to_seq(element_get_kmer(path_nodes[360]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[360]==forward);
  CU_ASSERT(path_labels[359]==Cytosine);

  CU_ASSERT_STRING_EQUAL("CTAACCCTAACCCTAACCCTAACCCTAACCC", binary_kmer_to_seq(element_get_kmer(path_nodes[399]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[399]==forward);
  CU_ASSERT(path_labels[398]==Cytosine);

  // now the new nodes :-)
  // with path_string test the whole lot in one go

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC"
                         "TAACCCTAACCCTAACCCTAACCCAACCCTAACCCTAACCCTAACCCTAACCC"
                         "TAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACC"
                         "CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAA"
                         "CCCTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACCC"
                         "CAACCCCAACCCCAACCCCAACCCCAACCCTAACCCCTAACCCTAACCCTAAC"
                         "CCTACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCCTA"
                         "ACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCC"
                         "TAACCCTAACCCTCGCGGTACCCTCAGCCGGCCCGCCCGCCCGGGTCTGACCT"
                         "GAGGAGAACTGTGCTCCGCCTTCAGAGTACCACCGAAATCTGTGCAGAGGACA"
                         "ACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCTGTGCTGAGGAGAACGCAAC"
                         "TCCGCCGGCGCAGGCGCAGAGAGGCGCGCCGCGCCGGCGCAGGCGCAGACACA"
                         "TGCTAGCGCGTCGGGGTGGAGGCGTGGCGCAGGCGCAGAGAGGCGCGCCGCGC"
                         "CGGCGCAGGCGCAGAGACACATGCTACCGCGTCCAGGGGTGGAGGCGTGGCGC"
                         "AGGCGCAGAGAGGCGCACCGCGCCGGCGCAGGCGCAGAGACACATGCTAGCGC"
                         "GTCC", path_string);

  CU_ASSERT_STRING_EQUAL("AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA", binary_kmer_to_seq(element_get_kmer(path_nodes[400]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[400]==reverse);

  // == these next two asserts are VERY IMPORTANT ==
  // they test whether we get the correct edges across points where we load more
  // nodes into the array

  CU_ASSERT(path_labels[399]==Thymine);
  CU_ASSERT(path_string[399]=='T'); //iqbal

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[401]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[401]==forward);
  CU_ASSERT(path_labels[400]==Adenine);
  CU_ASSERT(path_string[400]=='A');

  CU_ASSERT_STRING_EQUAL("ACCCTAACCCTAACCCTAACCCTAACCCTAA", binary_kmer_to_seq(element_get_kmer(path_nodes[402]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[402]==forward);
  CU_ASSERT(path_labels[401]==Adenine);
  CU_ASSERT(path_string[401]=='A');

  CU_ASSERT_STRING_EQUAL("CCCTAACCCTAACCCTAACCCTAACCCTAAC", binary_kmer_to_seq(element_get_kmer(path_nodes[403]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[403]==forward);
  CU_ASSERT(path_labels[402]==Cytosine);
  CU_ASSERT(path_string[402]=='C');


  CU_ASSERT_STRING_EQUAL("ACCCTAACCCTAACCCTAACCCCTAACCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[420]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[420]==forward);
  CU_ASSERT(path_labels[419]==Adenine);
  CU_ASSERT(path_string[419]=='A');

  CU_ASSERT_STRING_EQUAL("GAGAGGCGCACCGCGCCGGCGCAGGCGCAGA", binary_kmer_to_seq(element_get_kmer(path_nodes[780]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[780]==forward);
  CU_ASSERT(path_labels[779]==Adenine);
  CU_ASSERT(path_string[779]=='A');


  CU_ASSERT_STRING_EQUAL("CGCAGGCGCAGAGACACATGCTAGCGCGTCC", binary_kmer_to_seq(element_get_kmer(path_nodes[799]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[799]==forward);
  CU_ASSERT(path_labels[798]==Cytosine);
  CU_ASSERT(path_string[798]=='C');


  //prepare to load next 400

  //push everything along in the arrays:
  for (i=0; i<length_of_arrays-num_of_nodes_to_read; i++)
  {
    path_nodes[i]       = path_nodes[i+num_of_nodes_to_read];
    path_orientations[i]= path_orientations[i+num_of_nodes_to_read];
    path_labels[i]      = path_labels[i+num_of_nodes_to_read];
    path_string[i]      = path_string[i+num_of_nodes_to_read];
  }


  // Move final kmer to the start of seq. Note final means at the end of seq,
  // which has got num_of_nodes_to_read PLUS ONE in it
  // (you preloaded with the previous kmer)
  for (i=0; i<db_graph->kmer_size; i++)
  {
    //num_bases=num_nodes+k-1, but we have an extra base,  so last k bases are
    // the num_nodes+1,....num_nodes+k -th bases. But we count from 0...
    seq->seq[i]=seq->seq[num_of_nodes_to_read+i];
  }
  seq->seq[db_graph->kmer_size]='\0';


  //load the next 400 - or at least try to. This time there are only
  // 1140-(800+31-1)=310 nodes before you run out of file.
  expecting_new_fasta_entry=false;


  // === third batch load ===
  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);
  CU_ASSERT(retvalue==310);

  //check we have the transition correct, between previously loaded nodes and new ones

  //this node used to be at index 799:
   CU_ASSERT_STRING_EQUAL("CGCAGGCGCAGAGACACATGCTAGCGCGTCC", binary_kmer_to_seq(element_get_kmer(path_nodes[399]), db_graph->kmer_size, tmp_seq));
   CU_ASSERT(path_orientations[399]==forward);
   CU_ASSERT(path_labels[398]==Cytosine);
   CU_ASSERT(path_string[398]=='C');


  //finally, just check the last one. What is it's index? Well, we asked for 400, but only got 310. So we expect to find last one at 799-90=709.

   CU_ASSERT_STRING_EQUAL("ACCACTCTAAGCAAGAGAGCCCTGCAGTTGC", binary_kmer_to_seq(element_get_kmer(path_nodes[709]), db_graph->kmer_size, tmp_seq));
   CU_ASSERT(path_orientations[709]==reverse);
   CU_ASSERT(path_labels[708]==Thymine);
   CU_ASSERT(path_string[708]=='T');


   free_sequence(&seq);
   free(path_nodes);
   free(path_orientations);
   free(path_labels);
   free(path_string);
   hash_table_free(&db_graph);


   free(kmer_window->kmer);
   free(kmer_window);

 }



 void test_align_next_read_to_graph_and_return_node_array()
 {
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // Set up db_graph
  int kmer_size = 17;
  int number_of_bits = 10;
  int bucket_size = 30;

  dBGraph* db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);

  // Read sequence
  int fq_quality_cutoff = 20;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_read = 0, seq_loaded = 0;

  load_se_filelist_into_graph_colour(
    "../data/test/graph/person3.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0,
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  int max_read_length = 200;

  //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);

  //We are going to load all the bases into a single sliding window
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
  {
    die("%s:%i: Failed to malloc kmer sliding window", __FILE__, __LINE__);
  }

  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-db_graph->kmer_size-1));
  if (kmer_window->kmer==NULL)
  {
    die("%s:%i: Failed to malloc kmer_window->kmer", __FILE__, __LINE__);
  }

  kmer_window->nkmers=0;


  //end of intialisation



  //create file reader
  int file_reader(FILE * fp, Sequence * seq, int max_read_length,
                  boolean new_entry, boolean * full_entry)
  {
    long long ret;
    int offset = 0;
    if (new_entry == false){
      warn("new_entry must be true in hsi test function");
    }
    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);

    return ret;
  }


  /*
  Now we know person3.fa is all in the graph in colour 0. Now let's just try to get our array of nodes:
  person3 looks like this

  >read1 overlaps human chrom 1
  TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACC
  > read 2 overlaps human chrom 1
  ACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAAC
  > read 3 does not
  GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
  > read 3 does not
  GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
  > read 3 does not
  GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
  > read 4 does not, but has too low coverage
  TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  >extra read
  TTTTTTTTTTTTTTTTAAA
  */

  FILE* fp = fopen("../data/test/graph/person3.fa", "r");
  if (fp==NULL)
  {
    die("%s:%i: Cannot open ../data/test/graph/person3.fa", __FILE__, __LINE__);
  }


  dBNode* array_nodes[200];//in fact there are 43 17-mers in the first line of the fasta
  Orientation array_or[200];
  int colour = 0;
  boolean f_entry = true;

  int num_kmers
    = align_next_read_to_graph_and_return_node_array(
        fp, max_read_length, array_nodes, array_or, true, &f_entry, file_reader,
        seq, kmer_window, db_graph, colour);

  CU_ASSERT(f_entry==true);
  CU_ASSERT(num_kmers==43);
  BinaryKmer test_kmer, test_kmer_rev;
  seq_to_binary_kmer("TAACCCTAACCCTAACC", 17, &test_kmer);
  binary_kmer_reverse_complement(&test_kmer, 17, &test_kmer_rev);
  CU_ASSERT(array_or[0]==reverse);
  CU_ASSERT(binary_kmer_comparison_operator(array_nodes[0]->kmer,test_kmer_rev)==true );

  seq_to_binary_kmer("AACCCTAACCCTAACCC", 17, &test_kmer);
  binary_kmer_reverse_complement(&test_kmer, 17, &test_kmer_rev);
  CU_ASSERT(array_or[1]==forward);
  CU_ASSERT(binary_kmer_comparison_operator(array_nodes[1]->kmer,test_kmer)==true );

  seq_to_binary_kmer("TAACCCTAACCCTAACC", 17, &test_kmer);
  binary_kmer_reverse_complement(&test_kmer, 17, &test_kmer_rev);
  CU_ASSERT(array_or[42]==reverse);
  CU_ASSERT(binary_kmer_comparison_operator(array_nodes[42]->kmer,test_kmer_rev)==true );


  //get the next read
  f_entry=true;
  num_kmers = align_next_read_to_graph_and_return_node_array(fp, max_read_length, array_nodes, array_or, true, &f_entry, file_reader,
   seq, kmer_window, db_graph, colour);

  CU_ASSERT(f_entry==true);
  CU_ASSERT(num_kmers==45-17+1);//next read is 45 bases long

  seq_to_binary_kmer("ACCCTAACCCTAACCCT", 17, &test_kmer);
  binary_kmer_reverse_complement(&test_kmer, 17, &test_kmer_rev);
  CU_ASSERT(array_or[0]==forward);
  CU_ASSERT(binary_kmer_comparison_operator(array_nodes[0]->kmer,test_kmer)==true );

  seq_to_binary_kmer("CCCTAACCCTAACCCTA", 17, &test_kmer);
  binary_kmer_reverse_complement(&test_kmer, 17, &test_kmer_rev);
  CU_ASSERT(array_or[1]==forward);
  CU_ASSERT(binary_kmer_comparison_operator(array_nodes[1]->kmer,test_kmer)==true );

  seq_to_binary_kmer("CCTAACCCTAACCCTAA", 17, &test_kmer);
  binary_kmer_reverse_complement(&test_kmer, 17, &test_kmer_rev);
  CU_ASSERT(array_or[2]==forward);
  CU_ASSERT(binary_kmer_comparison_operator(array_nodes[2]->kmer,test_kmer)==true );

  seq_to_binary_kmer("CTAACCCTAACCCTAAC", 17, &test_kmer);
  binary_kmer_reverse_complement(&test_kmer, 17, &test_kmer_rev);
  CU_ASSERT(array_or[28]==forward);
  CU_ASSERT(binary_kmer_comparison_operator(array_nodes[28]->kmer,test_kmer)==true );



  free(kmer_window->kmer);
  free(kmer_window);
  free_sequence(&seq);
  hash_table_free(&db_graph);
}



void test_read_next_variant_from_full_flank_file()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // Load a single file containing a SNP between two reads, and use detect_vars
  // and trusted SV caller to dump a pair of files.  Then load these dumped
  // files and check we get what we expect.

  // First set up the hash/graph
  int kmer_size = 5;
  int number_of_bits = 6;
  int bucket_size = 10;

  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);

  // Read sequence
  int fq_quality_cutoff = 20;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_read = 0, seq_loaded = 0;


  //local function
  int this_fasta_file_reader(FILE * fp, Sequence * seq, int max_read_length,
                             boolean new_entry, boolean * full_entry)
  {
    long long ret;
    int offset = 0;
    if (new_entry == false){
      offset = db_graph->kmer_size;
    }
    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);

    return ret;
  }

  //mallocing
  int max_read_length = 50;
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);

  Sequence * seq_inc_prev_kmer = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
  }
  alloc_sequence(seq_inc_prev_kmer,max_read_length+kmer_size,LINE_MAX);


  //We are going to load all the bases into a single sliding window
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
  {
    die("%s:%i: Failed to malloc kmer sliding window", __FILE__, __LINE__);
  }


  //  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-db_graph->kmer_size-1));
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length+1));
  if (kmer_window->kmer==NULL)
  {
    die("%s:%i: Failed to malloc kmer_window->kmer", __FILE__, __LINE__);
  }
  kmer_window->nkmers=0;

  //end mallocing


  // this fasta
  // >read1
  // ACGCAGGCGCTGATGAGCTCCCGTTAGATAGT
  // >read2
  // ACGCAGGCGCTTATGAGCTCCCGTTAGATAGT

  // if I break this into flank, branches etc
  //  ACGCAGGCGCT G ATGAGCTCCCGTTAGATAGT
  //  ACGCAGGCGCT T ATGAGCTCCCGTTAGATAGT

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/variations/one_person_with_SNP.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0,
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  FILE* fout_bubble = fopen("../data/tempfiles_can_be_deleted/tmp_test.bubbles", "w");
  if (fout_bubble==NULL)
  {
    die("Unable to open: ../data/tempfiles_can_be_deleted/tmp_test.bubbles");
  }

  int max_branch_len=10;

  /*  GraphInfo* ginfo = graph_info_alloc_and_init();
  graph_info_set_seq(ginfo, 0, seq_read);
  graph_info_set_mean_readlen(ginfo,0, 32);
  GraphAndModelInfo model_info;
  float repeat_geometric_param_mu = 0.8;//not used in this
  int num_chroms_in_expt=2;
  initialise_model_info(&model_info, ginfo, 50, repeat_geometric_param_mu, 
  -1, 2, EachColourADiploidSample, AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice);*/

  db_graph_detect_vars(fout_bubble, max_branch_len,db_graph, &detect_vars_condition_always_true,
		       &db_node_action_set_status_visited, &db_node_action_set_status_visited,
		       &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs, &print_no_extra_info,
		       false, NULL, NULL, &db_node_condition_always_true);
  fclose(fout_bubble);

  //output will look like this

  // >var_1_5p_flank
  // AGCTCAT
  // >branch_1_1
  // AAGCGC
  // >branch_1_2
  // CAGCGC
  // >var_1_3p_flank
  // CTGCGT

  // These are the kmers in the VariantFlanksAndBranches object as generated by detect_vars.
  //  5p flank
  //  AGCTC GCTCA CTCAT            << completely matches what is printed out
  //  allele 1
  //  CTCAT TCATA CATAA ATAAG TAAGC AAGCG AGCGC   <<<< note first and last nodes of two branches are SAME
  //  allele 2
  //  CTCAT TCATC CATCA ATCAG TCAGC CAGCG AGCGC
  //  3p flank
  //  AGCGC GCGCC AGGCG CAGGC CCTGC CGCAG ACGCA   <<< note starts on node same last node of the two branches



  //Now check that we can read back in what we have printed out.
  /*
  dBNode* flank5p[20];
  dBNode* one_allele[20];
  dBNode* other_allele[20];
  dBNode* flank3p[20];
  Orientation flank5p_or[20];
  Orientation one_allele_or[20];
  Orientation other_allele_or[20];
  Orientation flank3p_or[20];
  int len_flank5p=0;
  int len_one_allele=0;
  int len_other_allele=0;
  int len_flank3p=0;
  */


  VariantBranchesAndFlanks* var = alloc_VariantBranchesAndFlanks_object(40,40,40,40,kmer_size);
  if (var==NULL)
  {
    die("Failed to malloc var");
  }

  //FILE* var_fptr =  fopen("tmp_test_read_next_variant_from_full_flank_file_bubble.fff", "r");
  FILE* var_fptr = fopen("../data/tempfiles_can_be_deleted/tmp_test.bubbles", "r");



  read_next_variant_from_full_flank_file(var_fptr, max_read_length,
    var, db_graph, &this_fasta_file_reader, seq, seq_inc_prev_kmer, kmer_window);

  CU_ASSERT(var->len_flank5p==2 );
  CU_ASSERT(var->len_one_allele==6);
  CU_ASSERT(var->len_other_allele==6);
  CU_ASSERT(var->len_flank3p==6);
  CU_ASSERT(strcmp(var->var_name,"var_1")==0);
  // printf("ZAMZAM\n\n len flank 5p is %d\nvariant is %s\n%s\n%s\n%s\n\n", len_flank5p);
  BinaryKmer tmp_kmer, tmp_kmer_rev;
  BinaryKmer tmp_kmer2, tmp_kmer_rev2;

  seq_to_binary_kmer("AGCTC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank5p[0]->kmer, tmp_kmer));
  CU_ASSERT(var->flank5p_or[0]==forward);

  seq_to_binary_kmer("GCTCA", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank5p[1]->kmer, tmp_kmer));
  CU_ASSERT(var->flank5p_or[1]==forward);

  seq_to_binary_kmer("CTCAT", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank5p[2]->kmer, tmp_kmer_rev));
  CU_ASSERT(var->flank5p_or[2]==reverse);


  //now the alleles
  // reminder....
  // >var_1_5p_flank
  // AGCTCAT
  // >branch_1_1
  // AAGCGC
  // >branch_1_2
  // CAGCGC
  // >var_1_3p_flank
  // CTGCGT



  CU_ASSERT(var->len_one_allele   == 6);
  CU_ASSERT(var->len_other_allele == 6);
  //detect_vars will print alleles in alphabetical order of first base
  seq_to_binary_kmer("CTCAT", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  seq_to_binary_kmer("CTCAT", 5, &tmp_kmer2);
  binary_kmer_reverse_complement(&tmp_kmer2, 5, &tmp_kmer_rev2);

  //char zamzam[db_graph->kmer_size+1];
  //printf("\n\n\n IQIQIQ   start of one allele is %s\n", binary_kmer_to_seq(&(var->one_allele[0]->kmer), db_graph->kmer_size, zamzam));
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[0]  ->kmer, tmp_kmer_rev) );
  CU_ASSERT(var->one_allele_or[0]  ==reverse);
  CU_ASSERT(binary_kmer_comparison_operator(var->other_allele[0]->kmer, tmp_kmer_rev2) );
  CU_ASSERT(var->other_allele_or[0]==reverse);

  seq_to_binary_kmer("TCATA", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  seq_to_binary_kmer("TCATC", 5, &tmp_kmer2);
  binary_kmer_reverse_complement(&tmp_kmer2, 5, &tmp_kmer_rev2);
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[1]->kmer, tmp_kmer_rev) );
  CU_ASSERT(var->one_allele_or[1]==reverse);
  CU_ASSERT(binary_kmer_comparison_operator(var->other_allele[1]->kmer, tmp_kmer_rev2) );
  CU_ASSERT(var->other_allele_or[1]==reverse);


  seq_to_binary_kmer("CATAA", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  seq_to_binary_kmer("CATCA", 5, &tmp_kmer2);
  binary_kmer_reverse_complement(&tmp_kmer2, 5, &tmp_kmer_rev2);
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[2]->kmer, tmp_kmer) );
  CU_ASSERT(var->one_allele_or[2]==forward);
  CU_ASSERT(binary_kmer_comparison_operator(var->other_allele[2]->kmer, tmp_kmer2) );
  CU_ASSERT(var->other_allele_or[2]==forward);


  seq_to_binary_kmer("ATAAG", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  seq_to_binary_kmer("ATCAG", 5, &tmp_kmer2);
  binary_kmer_reverse_complement(&tmp_kmer2, 5, &tmp_kmer_rev2);
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[3]->kmer, tmp_kmer) );
  CU_ASSERT(var->one_allele_or[3]==forward);
  CU_ASSERT(binary_kmer_comparison_operator(var->other_allele[3]->kmer, tmp_kmer2) );
  CU_ASSERT(var->other_allele_or[3]==forward);


  seq_to_binary_kmer("TAAGC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  seq_to_binary_kmer("TCAGC", 5, &tmp_kmer2);
  binary_kmer_reverse_complement(&tmp_kmer2, 5, &tmp_kmer_rev2);
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[4]->kmer, tmp_kmer_rev) );
  CU_ASSERT(var->one_allele_or[4]==reverse);
  CU_ASSERT(binary_kmer_comparison_operator(var->other_allele[4]->kmer, tmp_kmer_rev2) );
  CU_ASSERT(var->other_allele_or[4]==reverse);


  seq_to_binary_kmer("AAGCG", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  seq_to_binary_kmer("CAGCG", 5, &tmp_kmer2);
  binary_kmer_reverse_complement(&tmp_kmer2, 5, &tmp_kmer_rev2);
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[5]->kmer, tmp_kmer) );
  CU_ASSERT(var->one_allele_or[5]==forward);
  CU_ASSERT(binary_kmer_comparison_operator(var->other_allele[5]->kmer, tmp_kmer2) );
  CU_ASSERT(var->other_allele_or[5]==forward);


  seq_to_binary_kmer("AGCGC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  seq_to_binary_kmer("AGCGC", 5, &tmp_kmer2);
  binary_kmer_reverse_complement(&tmp_kmer2, 5, &tmp_kmer_rev2);
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[6]->kmer, tmp_kmer) );
  CU_ASSERT(var->one_allele_or[6]==forward);
  CU_ASSERT(binary_kmer_comparison_operator(var->other_allele[6]->kmer, tmp_kmer2) );
  CU_ASSERT(var->other_allele_or[6]==forward);


  //finally the 3p flank

  //char zamzam[db_graph->kmer_size+1];
  //printf("\n\n\n IQIQIQ   start of 3p flank is %s\n", binary_kmer_to_seq(&(var->flank3p[0]->kmer), db_graph->kmer_size, zamzam));


  seq_to_binary_kmer("AGCGC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank3p[0]->kmer, tmp_kmer));
  CU_ASSERT(var->flank3p_or[0]==forward);


  seq_to_binary_kmer("GCGCC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank3p[1]->kmer, tmp_kmer));
  CU_ASSERT(var->flank3p_or[1]==forward);

  seq_to_binary_kmer("CGCCT", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank3p[2]->kmer, tmp_kmer_rev));
  CU_ASSERT(var->flank3p_or[2]==reverse);

  seq_to_binary_kmer("GCCTG", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank3p[3]->kmer, tmp_kmer_rev));
  CU_ASSERT(var->flank3p_or[3]==reverse);


  seq_to_binary_kmer("CCTGC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank3p[4]->kmer, tmp_kmer));
  CU_ASSERT(var->flank3p_or[4]==forward);

  seq_to_binary_kmer("CTGCG", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank3p[5]->kmer, tmp_kmer_rev));
  CU_ASSERT(var->flank3p_or[5]==reverse);

  seq_to_binary_kmer("TGCGT", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank3p[6]->kmer, tmp_kmer_rev));
  CU_ASSERT(var->flank3p_or[6]==reverse);


  hash_table_free(&db_graph);
  free(kmer_window->kmer);
  free(kmer_window);
  free_sequence(&seq);
  free_sequence(&seq_inc_prev_kmer);
  free_VariantBranchesAndFlanks_object(var);
  //  graph_info_free(ginfo);
}


void test_read_next_variant_from_full_flank_file_2()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  /*

  Test an example where both branches are shorter than the kmer

  >read1
  ACGCAGGCATGAGCTCCCGTTAGATAGT
  >read2
  ACGCAGGCATGGATGAGCTCCCGTTAGATAGT
  with k=5 we get these branches

  >var_1_5p_flank length:7
  AGCTCAT
  >branch_1_1 length:6
  CCATGC
  >branch_1_2 length:2
  GC
  >var_1_3p_flank length:6
  CTGCGT
   */


   //first set up the hash/graph
  int kmer_size = 5;
  int number_of_bits = 6;
  int bucket_size = 10;

  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size, 10, kmer_size);

  // Read sequence
  int fq_quality_cutoff = 20;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_read = 0, seq_loaded = 0;


  //local function
  int this_fasta_file_reader(FILE * fp, Sequence * seq, int max_read_length,
   boolean new_entry, boolean * full_entry)
  {
    long long ret;
    int offset = 0;
    if (new_entry == false){
      offset = db_graph->kmer_size;
    }
    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);

    return ret;
  }

  //mallocing
  int max_read_length = 50;
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence");
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);

  Sequence * seq_inc_prev_kmer = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence");
  }
  alloc_sequence(seq_inc_prev_kmer,max_read_length+kmer_size,LINE_MAX);


  //We are going to load all the bases into a single sliding window
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
  {
    die("Failed to malloc kmer sliding window in a test. Exit.");
  }


  //  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-db_graph->kmer_size-1));
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length+1));
  if (kmer_window->kmer==NULL)
  {
    die("Failed to malloc kmer_window->kmer in test . Exit.");
  }
  kmer_window->nkmers=0;

  //end mallocing

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/variations/one_person_with_SNP_one_branch_shorter_than_k_if_k_is_5.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0,
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  FILE* fout_bubble = fopen("../data/tempfiles_can_be_deleted/tmp_test.bubbles", "w");
  
  if(fout_bubble==NULL)
  {
    die("Unable to open: ../data/tempfiles_can_be_deleted/tmp_test.bubbles");
  }

  int max_branch_len=10;

  db_graph_detect_vars(fout_bubble, max_branch_len,db_graph, &detect_vars_condition_always_true,
		       &db_node_action_set_status_visited, &db_node_action_set_status_visited,
		       &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs, &print_no_extra_info,
		       false, NULL, NULL, db_node_condition_always_true);
  fclose(fout_bubble);

  //output will look like this

  // >var_1_5p_flank length:7
  // AGCTCAT
  // >branch_1_1 length:6
  // CCATGC
  // >branch_1_2 length:2
  // GC
  // >var_1_3p_flank length:6
  // CTGCGT


  // These are the kmers in the VariantFlanksAndBranches object as generated by detect_vars.
  //  5p flank
  // AGCTC  GCTCA  ATGAG   (fw fw rev)
  // branch 1
  // ATGAG GATGA CATCC ATCCA ATGGA CATGG CATGC  (rev rev fw fw rev rev fw)
  // branch2
  // ATGAG CATGA CATGC  (rev rev fw)
  // 3p flank
  // CATGC ATGCC AGGCA CAGGC CCTGC CGCAG ACGCA (fw fw rev rev fw rev rev)


  //Now check that we can read back in what we have printed out.
  /*
  dBNode* flank5p[20];
  dBNode* one_allele[20];
  dBNode* other_allele[20];
  dBNode* flank3p[20];
  Orientation flank5p_or[20];
  Orientation one_allele_or[20];
  Orientation other_allele_or[20];
  Orientation flank3p_or[20];
  int len_flank5p=0;
  int len_one_allele=0;
  int len_other_allele=0;
  int len_flank3p=0;
  */

  VariantBranchesAndFlanks* var = alloc_VariantBranchesAndFlanks_object(40,40,40,40,kmer_size);
  if (var==NULL)
  {
    die("Failed to malloc var 2");
  }


  //FILE* var_fptr =  fopen("tmp_test_read_next_variant_from_full_flank_file_bubble.fff", "r");
  FILE* var_fptr =  fopen("../data/tempfiles_can_be_deleted/tmp_test.bubbles", "r");


  read_next_variant_from_full_flank_file(var_fptr, max_read_length,
    var, db_graph, &this_fasta_file_reader, seq, seq_inc_prev_kmer, kmer_window);

  CU_ASSERT(var->len_flank5p==2 );
  CU_ASSERT(var->len_one_allele==6);
  CU_ASSERT(var->len_other_allele==2);
  CU_ASSERT(var->len_flank3p==6);
  // printf("ZAMZAM\n\n len flank 5p is %d\nvariant is %s\n%s\n%s\n%s\n\n", len_flank5p);
  BinaryKmer tmp_kmer, tmp_kmer_rev;
  BinaryKmer tmp_kmer2, tmp_kmer_rev2;


  seq_to_binary_kmer("AGCTC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank5p[0]->kmer, tmp_kmer));
  CU_ASSERT(var->flank5p_or[0]==forward);

  seq_to_binary_kmer("GCTCA", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank5p[1]->kmer, tmp_kmer));
  CU_ASSERT(var->flank5p_or[1]==forward);

  seq_to_binary_kmer("CTCAT", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank5p[2]->kmer, tmp_kmer_rev));
  CU_ASSERT(var->flank5p_or[2]==reverse);


  //now the alleles
  // branch 1
  // ATGAG GATGA CATCC ATCCA ATGGA CATGG CATGC  (rev rev fw fw rev rev fw)
  // branch2
  // ATGAG CATGA CATGC  (rev rev fw)


  //detect_vars will print alleles in alphabetical order of first base
  seq_to_binary_kmer("ATGAG", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  seq_to_binary_kmer("ATGAG", 5, &tmp_kmer2);
  binary_kmer_reverse_complement(&tmp_kmer2, 5, &tmp_kmer_rev2);

  //char zamzam[db_graph->kmer_size+1];
  //printf("\n\n\n IQIQIQ   start of one allele is %s\n", binary_kmer_to_seq(&(var->one_allele[0]->kmer), db_graph->kmer_size, zamzam));
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[0]  ->kmer, tmp_kmer) );
  CU_ASSERT(var->one_allele_or[0]  ==reverse);
  CU_ASSERT(binary_kmer_comparison_operator(var->other_allele[0]->kmer, tmp_kmer2) );
  CU_ASSERT(var->other_allele_or[0]==reverse);

  seq_to_binary_kmer("GATGA", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  seq_to_binary_kmer("CATGA", 5, &tmp_kmer2);
  binary_kmer_reverse_complement(&tmp_kmer2, 5, &tmp_kmer_rev2);
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[1]->kmer, tmp_kmer) );
  CU_ASSERT(var->one_allele_or[1]==reverse);
  CU_ASSERT(binary_kmer_comparison_operator(var->other_allele[1]->kmer, tmp_kmer2) );
  CU_ASSERT(var->other_allele_or[1]==reverse);


  seq_to_binary_kmer("CATCC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  seq_to_binary_kmer("CATGC", 5, &tmp_kmer2);
  binary_kmer_reverse_complement(&tmp_kmer2, 5, &tmp_kmer_rev2);
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[2]->kmer, tmp_kmer) );
  CU_ASSERT(var->one_allele_or[2]==forward);
  CU_ASSERT(binary_kmer_comparison_operator(var->other_allele[2]->kmer, tmp_kmer2) );
  CU_ASSERT(var->other_allele_or[2]==forward);


  seq_to_binary_kmer("ATCCA", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[3]->kmer, tmp_kmer) );
  CU_ASSERT(var->one_allele_or[3]==forward);

  seq_to_binary_kmer("ATGGA", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[4]->kmer, tmp_kmer) ); ///ZAM THIS IS WHERE T GOES
  CU_ASSERT(var->one_allele_or[4]==reverse);

  seq_to_binary_kmer("CATGG", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[5]->kmer, tmp_kmer) );
  CU_ASSERT(var->one_allele_or[5]==reverse);

  seq_to_binary_kmer("CATGC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[6]->kmer, tmp_kmer) );
  CU_ASSERT(var->one_allele_or[6]==forward);


  //finally the 3p flank
  // CATGC ATGCC AGGCA CAGGC CCTGC CGCAG ACGCA (fw fw rev rev fw rev rev)

  seq_to_binary_kmer("CATGC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank3p[0]->kmer, tmp_kmer));
  CU_ASSERT(var->flank3p_or[0]==forward);


  seq_to_binary_kmer("ATGCC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank3p[1]->kmer, tmp_kmer));
  CU_ASSERT(var->flank3p_or[1]==forward);

  seq_to_binary_kmer("AGGCA", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank3p[2]->kmer, tmp_kmer));
  CU_ASSERT(var->flank3p_or[2]==reverse);

  seq_to_binary_kmer("CAGGC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank3p[3]->kmer, tmp_kmer));
  CU_ASSERT(var->flank3p_or[3]==reverse);


  seq_to_binary_kmer("CCTGC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank3p[4]->kmer, tmp_kmer));
  CU_ASSERT(var->flank3p_or[4]==forward);

  seq_to_binary_kmer("CGCAG", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank3p[5]->kmer, tmp_kmer));
  CU_ASSERT(var->flank3p_or[5]==reverse);

  seq_to_binary_kmer("ACGCA", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank3p[6]->kmer, tmp_kmer));
  CU_ASSERT(var->flank3p_or[6]==reverse);


  hash_table_free(&db_graph);
  free(kmer_window->kmer);
  free(kmer_window);
  free_sequence(&seq);
  free_sequence(&seq_inc_prev_kmer);
  free_VariantBranchesAndFlanks_object(var);
}


void test_read_next_variant_from_full_flank_file_3()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  /*
  Test an example where both branches are shorter than the kmer,
  this time the BOTH will be <k

  >read1
  ACTATCTAACGGGAGCTCATCCAGCGT
  >read2
  ACTATCTAACGGGAGCTCATCCCAGCGT

  with k=5 we get these branches

  >var_1_5p_flank length:7
  ACGCTGG
  >branch_1_1 length:3
  ATG
  >branch_1_2 length:4
  GATG
  >var_1_3p_flank length:4
  AGCT

  */

  //first set up the hash/graph
  int kmer_size = 5;
  int number_of_bits = 10;
  int bucket_size = 100;

  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);

  // Read in sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_read = 0, seq_loaded = 0;

  //local function
  int this_fasta_file_reader(FILE * fp, Sequence * seq, int max_read_length,
                             boolean new_entry, boolean * full_entry)
  {
    long long ret;
    int offset = 0;
    if (new_entry == false){
      offset = db_graph->kmer_size;
    }

    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);

    return ret;
  }

  //mallocing
  int max_read_length = 50;
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence");
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);

  Sequence * seq_inc_prev_kmer = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence");
  }
  alloc_sequence(seq_inc_prev_kmer,max_read_length+kmer_size,LINE_MAX);


  //We are going to load all the bases into a single sliding window
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
  {
    die("Failed to malloc kmer sliding window in a test. Exit.");
  }

  //  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-db_graph->kmer_size-1));
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length+1));

  if (kmer_window->kmer==NULL)
  {
    die("Failed to malloc kmer_window->kmer in test . Exit.");
  }

  kmer_window->nkmers = 0;


  //end mallocing

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/variations/one_person_with_SNP_both_branches_shorter_than_k_if_k_is_5.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0,
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);


  FILE* fout_bubble = fopen("../data/tempfiles_can_be_deleted/tmp_test.bubbles", "w");
  if (fout_bubble==NULL)
  {
    die("Unable to open: ../data/tempfiles_can_be_deleted/tmp_test.bubbles");
  }

  int max_branch_len=10;

  db_graph_detect_vars(fout_bubble, max_branch_len,db_graph,
                       &detect_vars_condition_always_true,
                       &db_node_action_set_status_visited,
                       &db_node_action_set_status_visited,
                       &element_get_colour_union_of_all_colours,
                       &element_get_covg_union_of_all_covgs,
                       &print_no_extra_info, false, NULL, NULL, &db_node_condition_always_true);
  fclose(fout_bubble);

  // output will look like this
  // >var_1_5p_flank length:7
  // ACGCTGG
  // >branch_1_1 length:3
  // ATG
  // >branch_1_2 length:4
  // GATG
  // >var_1_3p_flank length:4
  // AGCT


  // These are the kmers in the VariantFlanksAndBranches object as generated by detect_vars.
  // 5p flank
  // ACGCT CAGCG CCAGC (fw rev rev)
  // branch1
  // CCAGC CTGGA ATCCA CATCC (rev fw rev rev)
  // branch2
  // CCAGC CCCAG TCCCA ATCCC CATCC (rev rev rev rev rev)
  // 3p flank
  // CATCC GATGA ATGAG GCTCA AGCTC (rev fw fw rev rev)


  //Now check that we can read back in what we have printed out.
  /*
  dBNode* flank5p[20];
  dBNode* one_allele[20];
  dBNode* other_allele[20];
  dBNode* flank3p[20];
  Orientation flank5p_or[20];
  Orientation one_allele_or[20];
  Orientation other_allele_or[20];
  Orientation flank3p_or[20];
  int len_flank5p=0;
  int len_one_allele=0;
  int len_other_allele=0;
  int len_flank3p=0;
  */

  VariantBranchesAndFlanks* var = alloc_VariantBranchesAndFlanks_object(40,40,40,40,kmer_size);
  if (var==NULL)
  {
    die("Failed to malloc var 3");
  }


  //FILE* var_fptr =  fopen("tmp_test_read_next_variant_from_full_flank_file_bubble.fff", "r");
  FILE* var_fptr = fopen("../data/tempfiles_can_be_deleted/tmp_test.bubbles", "r");



  read_next_variant_from_full_flank_file(var_fptr, max_read_length,
    var, db_graph, &this_fasta_file_reader, seq, seq_inc_prev_kmer, kmer_window);


  CU_ASSERT(var->len_flank5p==2 );
  CU_ASSERT(var->len_one_allele==3);
  CU_ASSERT(var->len_other_allele==4);
  CU_ASSERT(var->len_flank3p==4);
  // printf("ZAMZAM\n\n len flank 5p is %d\nvariant is %s\n%s\n%s\n%s\n\n", len_flank5p);
  BinaryKmer tmp_kmer, tmp_kmer_rev;
  BinaryKmer tmp_kmer2, tmp_kmer_rev2;


  // These are the kmers in the VariantFlanksAndBranches object as generated by detect_vars.
  // 5p flank
  // ACGCT CAGCG CCAGC (fw rev rev)

  seq_to_binary_kmer("ACGCT", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank5p[0]->kmer, tmp_kmer));
  CU_ASSERT(var->flank5p_or[0]==forward);

  seq_to_binary_kmer("CAGCG", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank5p[1]->kmer, tmp_kmer));
  CU_ASSERT(var->flank5p_or[1]==reverse);

  seq_to_binary_kmer("CCAGC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank5p[2]->kmer, tmp_kmer));
  CU_ASSERT(var->flank5p_or[2]==reverse);


  //now the alleles
  // branch1
  // CCAGC CTGGA ATCCA CATCC (rev fw rev rev)
  // branch2
  // CCAGC CCCAG TCCCA ATCCC CATCC (rev rev rev rev rev)

  //detect_vars will print alleles in alphabetical order of first base
  seq_to_binary_kmer("CCAGC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  seq_to_binary_kmer("CCAGC", 5, &tmp_kmer2);
  binary_kmer_reverse_complement(&tmp_kmer2, 5, &tmp_kmer_rev2);

  //char zamzam[db_graph->kmer_size+1];
  //printf("\n\n\n IQIQIQ   start of one allele is %s\n", binary_kmer_to_seq(&(var->one_allele[0]->kmer), db_graph->kmer_size, zamzam));
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[0]  ->kmer, tmp_kmer) );
  CU_ASSERT(var->one_allele_or[0]  ==reverse);
  CU_ASSERT(binary_kmer_comparison_operator(var->other_allele[0]->kmer, tmp_kmer2) );
  CU_ASSERT(var->other_allele_or[0]==reverse);

  seq_to_binary_kmer("CTGGA", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  seq_to_binary_kmer("CCCAG", 5, &tmp_kmer2);
  binary_kmer_reverse_complement(&tmp_kmer2, 5, &tmp_kmer_rev2);
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[1]->kmer, tmp_kmer) );
  CU_ASSERT(var->one_allele_or[1]==forward);
  CU_ASSERT(binary_kmer_comparison_operator(var->other_allele[1]->kmer, tmp_kmer2) );
  CU_ASSERT(var->other_allele_or[1]==reverse);


  seq_to_binary_kmer("ATCCA", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  seq_to_binary_kmer("TCCCA", 5, &tmp_kmer2);
  binary_kmer_reverse_complement(&tmp_kmer2, 5, &tmp_kmer_rev2);
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[2]->kmer, tmp_kmer) );
  CU_ASSERT(var->one_allele_or[2]==reverse);
  CU_ASSERT(binary_kmer_comparison_operator(var->other_allele[2]->kmer, tmp_kmer2) );
  CU_ASSERT(var->other_allele_or[2]==reverse);

  seq_to_binary_kmer("CATCC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  seq_to_binary_kmer("ATCCC", 5, &tmp_kmer2);
  binary_kmer_reverse_complement(&tmp_kmer2, 5, &tmp_kmer_rev2);
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[3]->kmer, tmp_kmer) );
  CU_ASSERT(var->one_allele_or[3]==reverse);
  CU_ASSERT(binary_kmer_comparison_operator(var->other_allele[3]->kmer, tmp_kmer2) );
  CU_ASSERT(var->other_allele_or[3]==reverse);


  seq_to_binary_kmer("CATCC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->other_allele[4]->kmer, tmp_kmer) );
  CU_ASSERT(var->one_allele_or[3]==reverse);


  //finally the 3p flank
  // 3p flank
  // CATCC GATGA ATGAG GCTCA AGCTC (rev fw fw rev rev)


  seq_to_binary_kmer("CATCC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank3p[0]->kmer, tmp_kmer));
  CU_ASSERT(var->flank3p_or[0]==reverse);


  seq_to_binary_kmer("GATGA", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank3p[1]->kmer, tmp_kmer));
  CU_ASSERT(var->flank3p_or[1]==forward);

  seq_to_binary_kmer("ATGAG", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank3p[2]->kmer, tmp_kmer));
  CU_ASSERT(var->flank3p_or[2]==forward);

  seq_to_binary_kmer("GCTCA", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank3p[3]->kmer, tmp_kmer));
  CU_ASSERT(var->flank3p_or[3]==reverse);


  seq_to_binary_kmer("AGCTC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank3p[4]->kmer, tmp_kmer));
  CU_ASSERT(var->flank3p_or[4]==reverse);


  hash_table_free(&db_graph);
  free(kmer_window->kmer);
  free(kmer_window);
  free_sequence(&seq);
  free_sequence(&seq_inc_prev_kmer);
  free_VariantBranchesAndFlanks_object(var);
}


void test_read_next_variant_from_full_flank_file_4()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  /*
  Test an example where both branches are shorter than the kmer,
  this time the BOTH will be <k, AND there is no 3p flank printed

  >read1
  ACTATCTAACGGGAGCTCATCCAGC
  >read2
  ACTATCTAACGGGAGCTCATCCCAGC

  with k=5 we get these branches

  >var_1_5p_flank length:9
  AGCTCATCC
  >branch_1_1 length:3
  AGC
  >branch_1_2 length:4
  CAGC
  >var_1_3p_flank length:0
  */

  // first set up the hash/graph
  //for christs sake dont bugger around with the mem height and width -
  // affects whether it calls it in fw or reverse dir

  int kmer_size = 5;
  int number_of_bits = 6;
  int bucket_size = 10;

  dBGraph * db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);

  // Read in sequence
  int fq_quality_cutoff = 20;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_read = 0, seq_loaded = 0;

  //local function
  int this_fasta_file_reader(FILE * fp, Sequence * seq, int max_read_length,
                             boolean new_entry, boolean * full_entry)
  {
    long long ret;
    int offset = 0;
    if (new_entry == false){
      offset = db_graph->kmer_size;
    }
    ret = read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);

    return ret;
  }

  //mallocing
  int max_read_length = 50;
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence");
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);

  Sequence * seq_inc_prev_kmer = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence");
  }
  alloc_sequence(seq_inc_prev_kmer,max_read_length+kmer_size,LINE_MAX);


  //We are going to load all the bases into a single sliding window
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
  {
    die("Failed to malloc kmer sliding window in a test. Exit.");
  }


  //  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-db_graph->kmer_size-1));
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length+1));
  if (kmer_window->kmer==NULL)
  {
    die("Failed to malloc kmer_window->kmer in test. Exit.");
  }
  kmer_window->nkmers=0;



  //end mallocing

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/variations/one_person_with_SNP_both_branches_"
      "shorter_than_k_if_k_is_5_and_no_3p_flank.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0,
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  FILE* fout_bubble = fopen("../data/tempfiles_can_be_deleted/tmp_test.bubbles", "w");
  if (fout_bubble==NULL)
  {
    die("Unable to open: ../data/tempfiles_can_be_deleted/tmp_test.bubbles");
  }

  int max_branch_len=10;

  db_graph_detect_vars(fout_bubble, max_branch_len,db_graph,
		       &detect_vars_condition_always_true, &db_node_action_set_status_visited,
		       &db_node_action_set_status_visited, &element_get_colour_union_of_all_colours,
		       &element_get_covg_union_of_all_covgs, &print_no_extra_info,
		       false, NULL, NULL, &db_node_condition_always_true);
  fclose(fout_bubble);

  // output will look like this
  // >var_1_5p_flank length:9
  // AGCTCATCC
  // >branch_1_1 length:3
  // AGC
  // >branch_1_2 length:4
  // CAGC
  // >var_1_3p_flank length:0




  // These are the kmers in the VariantFlanksAndBranches object as generated by detect_vars.
  // 5p flank
  // AGCTC GCTCA ATGAG GATGA CATCC (fw fw rev rev fw)
  // branch 1
  // CATCC ATCCA CTGGA CCAGC (fw fw rev fw)
  // branch2
  // CATCC ATCCC TCCCA CCCAG CCAGC (fw fw fw fw fw )
  // 3p flank
  // CCAGC  (fw )



  //Now check that we can read back in what we have printed out.
  /*
  dBNode* flank5p[20];
  dBNode* one_allele[20];
  dBNode* other_allele[20];
  dBNode* flank3p[20];
  Orientation flank5p_or[20];
  Orientation one_allele_or[20];
  Orientation other_allele_or[20];
  Orientation flank3p_or[20];
  int len_flank5p=0;
  int len_one_allele=0;
  int len_other_allele=0;
  int len_flank3p=0;
  */

  VariantBranchesAndFlanks* var = alloc_VariantBranchesAndFlanks_object(40,40,40,40,kmer_size);
  if (var==NULL)
  {
    die("Failed to malloc var 4");
  }

  //FILE* var_fptr =  fopen("tmp_test_read_next_variant_from_full_flank_file_bubble.fff", "r");
  FILE* var_fptr = fopen("../data/tempfiles_can_be_deleted/tmp_test.bubbles", "r");

  read_next_variant_from_full_flank_file(var_fptr, max_read_length,
    var, db_graph, &this_fasta_file_reader,
    seq, seq_inc_prev_kmer, kmer_window);

  //printf("lengths are %d %d %d %d\n", var->len_flank5p, var->len_one_allele,
  //  var->len_other_allele, var->len_flank3p);
  CU_ASSERT(var->len_flank5p==4 );
  CU_ASSERT(var->len_one_allele==3);
  CU_ASSERT(var->len_other_allele==4);
  CU_ASSERT(var->len_flank3p==0);
  // printf("ZAMZAM\n\n len flank 5p is %d\nvariant is %s\n%s\n%s\n%s\n\n", len_flank5p);
  BinaryKmer tmp_kmer, tmp_kmer_rev;
  BinaryKmer tmp_kmer2, tmp_kmer_rev2;


  // These are the kmers in the VariantFlanksAndBranches object as generated by detect_vars.

  // 5p flank
  // AGCTC GCTCA ATGAG GATGA CATCC (fw fw rev rev fw)

  seq_to_binary_kmer("AGCTC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank5p[0]->kmer, tmp_kmer));
  CU_ASSERT(var->flank5p_or[0]==forward);

  seq_to_binary_kmer("GCTCA", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank5p[1]->kmer, tmp_kmer));
  CU_ASSERT(var->flank5p_or[1]==forward);

  seq_to_binary_kmer("ATGAG", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank5p[2]->kmer, tmp_kmer));
  CU_ASSERT(var->flank5p_or[2]==reverse);

  seq_to_binary_kmer("GATGA", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank5p[3]->kmer, tmp_kmer));
  CU_ASSERT(var->flank5p_or[3]==reverse);

  seq_to_binary_kmer("CATCC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank5p[4]->kmer, tmp_kmer));
  CU_ASSERT(var->flank5p_or[4]==forward);

  //now the alleles
  // branch 1
  // CATCC ATCCA CTGGA CCAGC (fw fw rev fw)
  // branch2
  // CATCC ATCCC TCCCA CCCAG CCAGC (fw fw fw fw fw )
  // 3p flank
  // CCAGC  (fw )

  //detect_vars will print alleles in alphabetical order of first base
  seq_to_binary_kmer("CATCC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  seq_to_binary_kmer("CATCC", 5, &tmp_kmer2);
  binary_kmer_reverse_complement(&tmp_kmer2, 5, &tmp_kmer_rev2);

  //char zamzam[db_graph->kmer_size+1];
  //printf("\n\n\n IQIQIQ   start of one allele is %s\n",
  //  binary_kmer_to_seq(&(var->one_allele[0]->kmer), db_graph->kmer_size, zamzam));
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[0]  ->kmer, tmp_kmer) );
  CU_ASSERT(var->one_allele_or[0]  ==forward);
  CU_ASSERT(binary_kmer_comparison_operator(var->other_allele[0]->kmer, tmp_kmer2) );
  CU_ASSERT(var->other_allele_or[0]==forward);

  seq_to_binary_kmer("ATCCA", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  seq_to_binary_kmer("ATCCC", 5, &tmp_kmer2);
  binary_kmer_reverse_complement(&tmp_kmer2, 5, &tmp_kmer_rev2);
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[1]->kmer, tmp_kmer) );
  CU_ASSERT(var->one_allele_or[1]==forward);
  CU_ASSERT(binary_kmer_comparison_operator(var->other_allele[1]->kmer, tmp_kmer2) );
  CU_ASSERT(var->other_allele_or[1]==forward);

  seq_to_binary_kmer("CTGGA", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  seq_to_binary_kmer("TCCCA", 5, &tmp_kmer2);
  binary_kmer_reverse_complement(&tmp_kmer2, 5, &tmp_kmer_rev2);
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[2]->kmer, tmp_kmer) );
  CU_ASSERT(var->one_allele_or[2]==reverse);
  CU_ASSERT(binary_kmer_comparison_operator(var->other_allele[2]->kmer, tmp_kmer2) );
  CU_ASSERT(var->other_allele_or[2]==forward);

  seq_to_binary_kmer("CCAGC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  seq_to_binary_kmer("CCCAG", 5, &tmp_kmer2);
  binary_kmer_reverse_complement(&tmp_kmer2, 5, &tmp_kmer_rev2);
  CU_ASSERT(binary_kmer_comparison_operator(var->one_allele[3]->kmer, tmp_kmer) );
  CU_ASSERT(var->one_allele_or[3]==forward);
  CU_ASSERT(binary_kmer_comparison_operator(var->other_allele[3]->kmer, tmp_kmer2) );
  CU_ASSERT(var->other_allele_or[3]==forward);


  seq_to_binary_kmer("CCAGC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->other_allele[4]->kmer, tmp_kmer) );
  CU_ASSERT(var->one_allele_or[3]==forward);


  //finally the 3p flank
  // 3p flank
  // CCAGC  (fw )
  seq_to_binary_kmer("CCAGC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(var->flank3p[0]->kmer, tmp_kmer));
  CU_ASSERT(var->flank3p_or[0]==forward);

  hash_table_free(&db_graph);
  free(kmer_window->kmer);
  free(kmer_window);
  free_sequence(&seq);
  free_sequence(&seq_inc_prev_kmer);
  free_VariantBranchesAndFlanks_object(var);
}



void _test_getting_readlength_dist_with_file(
  char fq_quality_cutoffs[4], char homopolymer_cutoffs[4], // num_of_tests = 4
  int kmerlist[3], // num_of_kmers = 3
  boolean remove_duplicates,
  char *file1, char* file2,
  unsigned long correct_answers_readlens[4][3][90])
{
  int num_of_tests = 4;
  int num_of_kmers = 3;
  int readlen_distrib_arrlen = 90;
  unsigned long readlen_counts[90];

  // Read in settings
  char ascii_fq_offset = 33;
  int into_colour = 0;

  int number_of_bits = 10;
  int bucket_size = 100;
  int max_retries = 82;

  boolean pe = (file2 != NULL);

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_read = 0, seq_loaded;

  int i,j,k;

  for(i = 0; i < num_of_tests; i++)
  {
    for(k = 0; k < num_of_kmers; k++)
    {
      int kmer_size = kmerlist[k];

      if(kmer_size < 32 * NUMBER_OF_BITFIELDS_IN_BINARY_KMER &&
         kmer_size > 32 * (NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1))
      {
        int fq_quality_cutoff = fq_quality_cutoffs[i];
        int homopolymer_cutoff = homopolymer_cutoffs[i];

        for(j = 0; j < readlen_distrib_arrlen; j++)
          readlen_counts[j] = 0;

	/*
        printf("Running test: kmer: %i; fq_quality_cutoff: %i; "
               "homopolymer_cutoff: %i; remove_dupes: %s\n",
               kmer_size, fq_quality_cutoff, homopolymer_cutoff,
               (remove_duplicates == true ? "on" : "off"));
	*/
        dBGraph *db_graph = hash_table_new(number_of_bits, bucket_size,
                                           max_retries, kmer_size);

        if(pe)
        {
          load_pe_filelists_into_graph_colour(
          file1, file2,
          fq_quality_cutoff, homopolymer_cutoff,
          remove_duplicates, ascii_fq_offset,
          into_colour, db_graph, 0,
          &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
          readlen_counts, readlen_distrib_arrlen, &subsample_null);
        }
        else
        {
          load_se_filelist_into_graph_colour(
            file1,
            fq_quality_cutoff, homopolymer_cutoff,
            remove_duplicates, ascii_fq_offset,
            into_colour, db_graph, 0,
            &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
            readlen_counts, readlen_distrib_arrlen, &subsample_null);
        }

        for(j = 0; j < readlen_distrib_arrlen; j++)
        {
          CU_ASSERT(readlen_counts[j] == correct_answers_readlens[i][k][j]);

          if(readlen_counts[j] != correct_answers_readlens[i][k][j])
          {
            warn("%s:%i: %s test: %i; k[%i]: %i; j: %i; got: %lu; correct: %lu\n",
                 __FILE__, __LINE__, (pe ? "pe" : "se"),
                 i, k, kmer_size, j,
                 readlen_counts[j], correct_answers_readlens[i][k][j]);
          }
        }

        hash_table_free(&db_graph);
      }
    }
  }
}

void test_getting_readlength_distribution()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 2)
  {
    warn("Tests only written for MAXK=31 and MAXK=63\n");
    return;
  }

  char *file1, *file2;

  int kmerlist[] = {15,19,39};

  // For each kmer try the following filter settings
  char fq_quality_cutoffs[] = {0,10,10,0};
  char homopolymer_cutoffs[] = {0,0,5,5};

  int num_of_tests = 4;
  int num_of_kmers = 3;
  int readlen_distrib_arrlen = 90;
  unsigned long correct_answers_readlens[4][3][90];

  // ===========================================================================
  // 1. Check that with single-ended parsing of fastq, when we get the read
  //   length distribution it correctly does so when there are low quality bases
  //   below threshold, homopolymers, Ns
  // ===========================================================================

  int i, j, k;

  for(i = 0; i < num_of_tests; i++)
    for(k = 0; k < num_of_kmers; k++)
      for(j = 0; j < readlen_distrib_arrlen; j++)
        correct_answers_readlens[i][k][j] = 0;

  //
  // TEST 0: fq_quality_cutoff=0, homopolymer_cutoff=0
  //

  // kmer=15, fq_quality_cutoff=0, homopolymer_cutoff=0
  // read1: 41,26,17
  // read2: 95
  // read2: 95
  // read3: 41,39,17
  correct_answers_readlens[0][0][17] = 2;
  correct_answers_readlens[0][0][26] = 1;
  correct_answers_readlens[0][0][39] = 1;
  correct_answers_readlens[0][0][41] = 2;
  correct_answers_readlens[0][0][89] = 2; // 95 dumped into last bin

  // kmer=19, fq_quality_cutoff=0, homopolymer_cutoff=0
  // read1: 41,26
  // read2: 95
  // read2: 95
  // read3: 41,39
  correct_answers_readlens[0][1][26] = 1;
  correct_answers_readlens[0][1][39] = 1;
  correct_answers_readlens[0][1][41] = 2;
  correct_answers_readlens[0][1][89] = 2; // 95 dumped into last bin

  // kmer=39, fq_quality_cutoff=0, homopolymer_cutoff=0
  // read1: 41
  // read2: 95
  // read2: 95
  // read3: 41,39
  correct_answers_readlens[0][2][39] = 1;
  correct_answers_readlens[0][2][41] = 2;
  correct_answers_readlens[0][2][89] = 2; // 95 dumped into last bin

  //
  // TEST 1: fq_quality_cutoff=10, homopolymer_cutoff=0
  //

  // kmer=15, fq_quality_cutoff=10, homopolymer_cutoff=0
  // read1: 39,18,17
  // read2: 41,36
  // read2: 41,36
  // read3: 16,39,16
  correct_answers_readlens[1][0][16] = 2;
  correct_answers_readlens[1][0][17] = 1;
  correct_answers_readlens[1][0][18] = 1;
  correct_answers_readlens[1][0][36] = 2;
  correct_answers_readlens[1][0][39] = 2;
  correct_answers_readlens[1][0][41] = 2;

  // kmer=19, fq_quality_cutoff=10, homopolymer_cutoff=0
  // read1: 39
  // read2: 41,36
  // read2: 41,36
  // read3: 39
  correct_answers_readlens[1][1][36] = 2;
  correct_answers_readlens[1][1][39] = 2;
  correct_answers_readlens[1][1][41] = 2;

  // kmer=39, fq_quality_cutoff=10, homopolymer_cutoff=0
  // read1: 39
  // read2: 41
  // read2: 41
  // read3: 39
  correct_answers_readlens[1][2][39] = 2;
  correct_answers_readlens[1][2][41] = 2;

  //
  // TEST 2: fq_quality_cutoff=10, homopolymer_cutoff=5
  //

  // kmer=15, fq_quality_cutoff=10, homopolymer_cutoff=5
  // read1: 23,18,17
  // read2: 29,36
  // read2: 29,36
  // read3: 16,39,16
  correct_answers_readlens[2][0][16] = 2;
  correct_answers_readlens[2][0][17] = 1;
  correct_answers_readlens[2][0][18] = 1;
  correct_answers_readlens[2][0][23] = 1;
  correct_answers_readlens[2][0][29] = 2;
  correct_answers_readlens[2][0][36] = 2;
  correct_answers_readlens[2][0][39] = 1;

  // kmer=19, fq_quality_cutoff=10, homopolymer_cutoff=5
  // read1: 23
  // read2: 29,36
  // read2: 29,36
  // read3: 39
  correct_answers_readlens[2][1][23] = 1;
  correct_answers_readlens[2][1][29] = 2;
  correct_answers_readlens[2][1][36] = 2;
  correct_answers_readlens[2][1][39] = 1;

  // kmer=39, fq_quality_cutoff=10, homopolymer_cutoff=5
  // read1: 
  // read2: 
  // read2: 
  // read3: 39
  correct_answers_readlens[2][2][39] = 1;

  //
  // TEST 3: fq_quality_cutoff=0, homopolymer_cutoff=5
  //

  // kmer=15, fq_quality_cutoff=0, homopolymer_cutoff=5
  // read1: 23,26,17
  // read2: 83
  // read2: 83
  // read3: 41,39,17
  correct_answers_readlens[3][0][17] = 2;
  correct_answers_readlens[3][0][23] = 1;
  correct_answers_readlens[3][0][26] = 1;
  correct_answers_readlens[3][0][39] = 1;
  correct_answers_readlens[3][0][41] = 1;
  correct_answers_readlens[3][0][83] = 2;

  // kmer=19, fq_quality_cutoff=0, homopolymer_cutoff=5
  // read1: 23,26
  // read2: 83
  // read2: 83
  // read3: 41,39
  correct_answers_readlens[3][1][23] = 1;
  correct_answers_readlens[3][1][26] = 1;
  correct_answers_readlens[3][1][39] = 1;
  correct_answers_readlens[3][1][41] = 1;
  correct_answers_readlens[3][1][83] = 2;

  // kmer=39, fq_quality_cutoff=0, homopolymer_cutoff=5
  // read1: 
  // read2: 83
  // read2: 83
  // read3: 41,39
  correct_answers_readlens[3][2][39] = 1;
  correct_answers_readlens[3][2][41] = 1;
  correct_answers_readlens[3][2][83] = 2;

  boolean remove_duplicates = false;

  // Test single-ended FASTQ without remove_dups
  file1 = "../data/test/pop_graph/readlen_distrib/readlen_distrib.fqlist";

  _test_getting_readlength_dist_with_file(fq_quality_cutoffs, homopolymer_cutoffs,
                                          kmerlist, remove_duplicates,
                                          file1, NULL,
                                          correct_answers_readlens);

  // Test single-ended SAMs without remove_dups
  file1 = "../data/test/pop_graph/readlen_distrib/readlen_distrib.samlist";

  _test_getting_readlength_dist_with_file(fq_quality_cutoffs, homopolymer_cutoffs,
                                          kmerlist, remove_duplicates,
                                          file1, NULL,
                                          correct_answers_readlens);

  // Test single-ended BAMs without remove_dups
  file1 = "../data/test/pop_graph/readlen_distrib/readlen_distrib.bamlist";

  _test_getting_readlength_dist_with_file(fq_quality_cutoffs, homopolymer_cutoffs,
                                          kmerlist, remove_duplicates,
                                          file1, NULL,
                                          correct_answers_readlens);

  // Now test with remove_dupes
  remove_duplicates = true;

  // Correct for duplicates in single end reads
  // test 0 kmer=15,19,39
  correct_answers_readlens[0][0][89] = 1;
  correct_answers_readlens[0][1][89] = 1;
  correct_answers_readlens[0][2][89] = 1;
  // test 1 kmer=15
  correct_answers_readlens[1][0][36] = 1;
  correct_answers_readlens[1][0][41] = 1;
  // test 1 kmer=19
  correct_answers_readlens[1][1][36] = 1;
  correct_answers_readlens[1][1][41] = 1;
  // test 1 kmer=39
  correct_answers_readlens[1][2][41] = 1;
  // test 2 kmer=15
  correct_answers_readlens[2][0][29] = 1;
  correct_answers_readlens[2][0][36] = 1;
  // test 2 kmer=19
  correct_answers_readlens[2][1][29] = 1;
  correct_answers_readlens[2][1][36] = 1;
  // test 3 kmer=15,19,39
  correct_answers_readlens[3][0][83] = 1;
  correct_answers_readlens[3][1][83] = 1;
  correct_answers_readlens[3][2][83] = 1;

  // Test single-ended FASTQ with remove_dups
  file1 = "../data/test/pop_graph/readlen_distrib/readlen_distrib.fqlist";

  _test_getting_readlength_dist_with_file(fq_quality_cutoffs, homopolymer_cutoffs,
                                          kmerlist, remove_duplicates,
                                          file1, NULL,
                                          correct_answers_readlens);

  // Test single-ended SAMs with remove_dups
  file1 = "../data/test/pop_graph/readlen_distrib/readlen_distrib.samlist";

  _test_getting_readlength_dist_with_file(fq_quality_cutoffs, homopolymer_cutoffs,
                                          kmerlist, remove_duplicates,
                                          file1, NULL,
                                          correct_answers_readlens);

  // Test single-ended BAMs with remove_dups
  file1 = "../data/test/pop_graph/readlen_distrib/readlen_distrib.bamlist";

  _test_getting_readlength_dist_with_file(fq_quality_cutoffs, homopolymer_cutoffs,
                                          kmerlist, remove_duplicates,
                                          file1, NULL,
                                          correct_answers_readlens);

  // ===========================================================================
  // 2. Check that with PAIRED-END parsing of fastq, when we get the read length
  //    distribution it correctly does so when there are low quality bases below
  //    threshold, homopolymers, Ns
  // ===========================================================================

  // We are going to load 2 files. One is the same as the previous fastq we
  // loaded, with an extra copy of the third read, and the other is its mate.
  // file1: read1, read2, read3, copy of read3
  // file2: read4, read5, copy of read5, copy of read5
  // so actually only the final pair will be removed by dup removal

  for(i = 0; i < num_of_tests; i++)
    for(k = 0; k < num_of_kmers; k++)
      for(j = 0; j < readlen_distrib_arrlen; j++)
        correct_answers_readlens[i][k][j] = 0;

  //
  // TEST 0: fq_quality_cutoff=0, homopolymer_cutoff=0
  //

  // kmer=15, fq_quality_cutoff=0, homopolymer_cutoff=0
  // read1: 41,26,17 | 99
  // read2: 95       | 99
  // read3: 41,39,17 | 99
  // read3: 41,39,17 | 99
  correct_answers_readlens[0][0][17] = 3;
  correct_answers_readlens[0][0][26] = 1;
  correct_answers_readlens[0][0][39] = 2;
  correct_answers_readlens[0][0][41] = 3;
  correct_answers_readlens[0][0][89] = 5;

  // kmer=19, fq_quality_cutoff=0, homopolymer_cutoff=0
  // read1: 41,26 | 99
  // read2: 95    | 99
  // read3: 41,39 | 99
  // read3: 41,39 | 99
  correct_answers_readlens[0][1][26] = 1;
  correct_answers_readlens[0][1][39] = 2;
  correct_answers_readlens[0][1][41] = 3;
  correct_answers_readlens[0][1][89] = 5;

  // kmer=39, fq_quality_cutoff=0, homopolymer_cutoff=0
  // read1: 41    | 99
  // read2: 95    | 99
  // read3: 41,39 | 99
  // read3: 41,39 | 99
  correct_answers_readlens[0][2][39] = 2;
  correct_answers_readlens[0][2][41] = 3;
  correct_answers_readlens[0][2][89] = 5;

  //
  // TEST 1: fq_quality_cutoff=10, homopolymer_cutoff=0
  //

  // kmer=15, fq_quality_cutoff=10, homopolymer_cutoff=0
  // read1: 39,18,17 | 45,36
  // read2: 41,36    | 99
  // read3: 16,39,16 | 99
  // read3: 16,39,16 | 99
  correct_answers_readlens[1][0][16] = 4;
  correct_answers_readlens[1][0][17] = 1;
  correct_answers_readlens[1][0][18] = 1;
  correct_answers_readlens[1][0][36] = 2;
  correct_answers_readlens[1][0][39] = 3;
  correct_answers_readlens[1][0][41] = 1;
  correct_answers_readlens[1][0][45] = 1;
  correct_answers_readlens[1][0][89] = 3; // 3x99s get dumped into final bin

  // kmer=19, fq_quality_cutoff=10, homopolymer_cutoff=0
  // read1: 39    | 45,36
  // read2: 41,36 | 99
  // read3: 39    | 99
  // read3: 39    | 99
  correct_answers_readlens[1][1][36] = 2;
  correct_answers_readlens[1][1][39] = 3;
  correct_answers_readlens[1][1][41] = 1;
  correct_answers_readlens[1][1][45] = 1;
  correct_answers_readlens[1][1][89] = 3; // 3x99s get dumped into final bin

  // kmer=39, fq_quality_cutoff=10, homopolymer_cutoff=0
  // read1: 39 | 45
  // read2: 41 | 99
  // read3: 39 | 99
  // read3: 39 | 99
  correct_answers_readlens[1][2][39] = 3;
  correct_answers_readlens[1][2][41] = 1;
  correct_answers_readlens[1][2][45] = 1;
  correct_answers_readlens[1][2][89] = 3; // 3x99s get dumped into final bin

  //
  // TEST 2: fq_quality_cutoff=10, homopolymer_cutoff=5
  //

  // kmer=15, fq_quality_cutoff=10, homopolymer_cutoff=5
  // read1: 23,18,17 | 45,36
  // read2: 29,36    | 99
  // read3: 16,39,16 | 99
  // read3: 16,39,16 | 99
  correct_answers_readlens[2][0][16] = 4;
  correct_answers_readlens[2][0][17] = 1;
  correct_answers_readlens[2][0][18] = 1;
  correct_answers_readlens[2][0][23] = 1;
  correct_answers_readlens[2][0][29] = 1;
  correct_answers_readlens[2][0][36] = 2;
  correct_answers_readlens[2][0][39] = 2;
  correct_answers_readlens[2][0][45] = 1;
  correct_answers_readlens[2][0][89] = 3; // 3x99s get dumped into final bin

  // kmer=19, fq_quality_cutoff=10, homopolymer_cutoff=5
  // read1: 23    | 45,36
  // read2: 29,36 | 99
  // read3: 39    | 99
  // read3: 39    | 99
  correct_answers_readlens[2][1][23] = 1;
  correct_answers_readlens[2][1][29] = 1;
  correct_answers_readlens[2][1][36] = 2;
  correct_answers_readlens[2][1][39] = 2;
  correct_answers_readlens[2][1][45] = 1;
  correct_answers_readlens[2][1][89] = 3; // 3x99s get dumped into final bin

  // kmer=39, fq_quality_cutoff=10, homopolymer_cutoff=5
  // read1:    | 45
  // read2:    | 99
  // read3: 39 | 99
  // read3: 39 | 99
  correct_answers_readlens[2][2][39] = 2;
  correct_answers_readlens[2][2][45] = 1;
  correct_answers_readlens[2][2][89] = 3; // 3x99s get dumped into final bin

  //
  // TEST 3: fq_quality_cutoff=0, homopolymer_cutoff=5
  //

  // kmer=15, fq_quality_cutoff=0, homopolymer_cutoff=5
  // read1: 23,26,17 | 99
  // read2: 83       | 99
  // read3: 41,39,17 | 99
  // read3: 41,39,17 | 99
  correct_answers_readlens[3][0][17] = 3;
  correct_answers_readlens[3][0][23] = 1;
  correct_answers_readlens[3][0][26] = 1;
  correct_answers_readlens[3][0][39] = 2;
  correct_answers_readlens[3][0][41] = 2;
  correct_answers_readlens[3][0][83] = 1;
  correct_answers_readlens[3][0][89] = 4; // 4x99s get dumped into final bin

  // kmer=19, fq_quality_cutoff=0, homopolymer_cutoff=5
  // read1: 23,26 | 99
  // read2: 83    | 99
  // read3: 41,39 | 99
  // read3: 41,39 | 99
  correct_answers_readlens[3][1][23] = 1;
  correct_answers_readlens[3][1][26] = 1;
  correct_answers_readlens[3][1][39] = 2;
  correct_answers_readlens[3][1][41] = 2;
  correct_answers_readlens[3][1][83] = 1;
  correct_answers_readlens[3][1][89] = 4; // 4x99s get dumped into final bin

  // kmer=39, fq_quality_cutoff=0, homopolymer_cutoff=5
  // read1:       | 99
  // read2: 83    | 99
  // read3: 41,39 | 99
  // read3: 41,39 | 99
  correct_answers_readlens[3][2][39] = 2;
  correct_answers_readlens[3][2][41] = 2;
  correct_answers_readlens[3][2][83] = 1;
  correct_answers_readlens[3][2][89] = 4; // 4x99s get dumped into final bin

  remove_duplicates = false;

  // Test paired-ended FASTQ without remove_dups
  file1 = "../data/test/pop_graph/readlen_distrib/readlen_distrib.pair1.fqlist";
  file2 = "../data/test/pop_graph/readlen_distrib/readlen_distrib.pair2.fqlist";

  _test_getting_readlength_dist_with_file(fq_quality_cutoffs, homopolymer_cutoffs,
                                          kmerlist, remove_duplicates,
                                          file1, file2,
                                          correct_answers_readlens);

  // Test paired-ended SAMs without remove_dups
  file1 = "../data/test/pop_graph/readlen_distrib/readlen_distrib.pair1.samlist";
  file2 = "../data/test/pop_graph/readlen_distrib/readlen_distrib.pair2.samlist";

  _test_getting_readlength_dist_with_file(fq_quality_cutoffs, homopolymer_cutoffs,
                                          kmerlist, remove_duplicates,
                                          file1, file2,
                                          correct_answers_readlens);

  // Test paired-ended BAMs without remove_dups
  file1 = "../data/test/pop_graph/readlen_distrib/readlen_distrib.pair1.bamlist";
  file2 = "../data/test/pop_graph/readlen_distrib/readlen_distrib.pair2.bamlist";

  _test_getting_readlength_dist_with_file(fq_quality_cutoffs, homopolymer_cutoffs,
                                          kmerlist, remove_duplicates,
                                          file1, file2,
                                          correct_answers_readlens);

  // Now test with remove_dupes
  remove_duplicates = true;

  // Correct for duplicates in pair end reads
  // test 0 kmer=15
  correct_answers_readlens[0][0][17] = 2;
  correct_answers_readlens[0][0][39] = 1;
  correct_answers_readlens[0][0][41] = 2;
  correct_answers_readlens[0][0][89] = 4;
  // test 0 kmer=19
  correct_answers_readlens[0][1][39] = 1;
  correct_answers_readlens[0][1][41] = 2;
  correct_answers_readlens[0][1][89] = 4;
  // test 0 kmer=39
  correct_answers_readlens[0][2][39] = 1;
  correct_answers_readlens[0][2][41] = 2;
  correct_answers_readlens[0][2][89] = 4;
  // test 1 kmer=15
  correct_answers_readlens[1][0][16] = 2;
  correct_answers_readlens[1][0][39] = 2;
  correct_answers_readlens[1][0][89] = 2;
  // test 1 kmer=19
  correct_answers_readlens[1][1][39] = 2;
  correct_answers_readlens[1][1][89] = 2;
  // test 1 kmer=39
  correct_answers_readlens[1][2][39] = 2;
  correct_answers_readlens[1][2][89] = 2;
  // test 2 kmer=15
  correct_answers_readlens[2][0][16] = 2;
  correct_answers_readlens[2][0][39] = 1;
  correct_answers_readlens[2][0][89] = 2;
  // test 2 kmer=19
  correct_answers_readlens[2][1][39] = 1;
  correct_answers_readlens[2][1][89] = 2;
  // test 2 kmer=39
  correct_answers_readlens[2][2][39] = 1;
  correct_answers_readlens[2][2][89] = 2;
  // test 3 kmer=15
  correct_answers_readlens[3][0][17] = 2;
  correct_answers_readlens[3][0][39] = 1;
  correct_answers_readlens[3][0][41] = 1;
  correct_answers_readlens[3][0][89] = 3;
  // test 3 kmer=19
  correct_answers_readlens[3][1][39] = 1;
  correct_answers_readlens[3][1][41] = 1;
  correct_answers_readlens[3][1][89] = 3;
  // test 3 kmer=39
  correct_answers_readlens[3][2][39] = 1;
  correct_answers_readlens[3][2][41] = 1;
  correct_answers_readlens[3][2][89] = 3;

  // Test paired-ended FASTQ with remove_dups
  file1 = "../data/test/pop_graph/readlen_distrib/readlen_distrib.pair1.fqlist";
  file2 = "../data/test/pop_graph/readlen_distrib/readlen_distrib.pair2.fqlist";

  _test_getting_readlength_dist_with_file(fq_quality_cutoffs, homopolymer_cutoffs,
                                          kmerlist, remove_duplicates,
                                          file1, file2,
                                          correct_answers_readlens);

  // Test paired-ended SAMs with remove_dups
  file1 = "../data/test/pop_graph/readlen_distrib/readlen_distrib.pair1.samlist";
  file2 = "../data/test/pop_graph/readlen_distrib/readlen_distrib.pair2.samlist";

  _test_getting_readlength_dist_with_file(fq_quality_cutoffs, homopolymer_cutoffs,
                                          kmerlist, remove_duplicates,
                                          file1, file2,
                                          correct_answers_readlens);

  // Test paired-ended BAMs with remove_dups
  file1 = "../data/test/pop_graph/readlen_distrib/readlen_distrib.pair1.bamlist";
  file2 = "../data/test/pop_graph/readlen_distrib/readlen_distrib.pair2.bamlist";

  _test_getting_readlength_dist_with_file(fq_quality_cutoffs, homopolymer_cutoffs,
                                          kmerlist, remove_duplicates,
                                          file1, file2,
                                          correct_answers_readlens);

}


void test_loading_binary_data_iff_it_overlaps_a_fixed_colour()
{
  if(NUMBER_OF_COLOURS <= 1)
  {
    warn("This test is redundant with only one colour\n");
    return;
  }

  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  int kmer_size = 3;
  int number_of_bits_pre = 4;
  int number_of_bits_post = 8;
  int bucket_size = 5;
  int max_retries = 10;
  //int max_chunk_len_reading_from_fasta = 200;

  dBGraph* db_graph_pre = hash_table_new(number_of_bits_pre, bucket_size,
                                         max_retries, kmer_size);

  dBGraph * db_graph_post;

  // We need the following arguments for the API but we will not use them -
  // for duplicate removal and homopolymer breaking
  int fq_quality_cutoff = 20;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_read = 0, seq_loaded = 0;

  /*
  >read1
  GATCGGGTGT
  >read1 copy
  GATCGGGTGT
  >read2
  GGCT
  >read2 copy
  GGCT
  >read3
  TAGG
  >read3 copy
  TAGG
  >read4=read 1 with a single base error
  GATCGGGAGT
  */

  load_se_filelist_into_graph_colour(
    "../data/test/graph/test_loading_binarynodes_if_overlap_a_colour.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph_pre, 0,
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  CU_ASSERT(seq_read == 46);
  CU_ASSERT(seq_loaded == 46);

  GraphInfo* ginfo = graph_info_alloc_and_init();
  graph_info_set_seq(ginfo, 0, seq_read);

  // Dump a single colour binary just of this graph
  db_graph_dump_single_colour_binary_of_colour0(
    "../data/tempfiles_can_be_deleted/dump_cortex_var_graph.singlecol.ctx",
    &db_node_condition_always_true, db_graph_pre, ginfo, BINVERSION);


  // OK, we now have dumped a binary corresponding to colour 0.
  // Now let's clean up, removing the bubble created by the single base error on
  // the 2nd copy of read 4
  db_graph_remove_low_coverage_nodes_ignoring_colours(1, db_graph_pre);

  // and dump a clean graph,
  db_graph_dump_binary("../data/tempfiles_can_be_deleted/dump_cortex_var_graph.clean.ctx",
                       &db_node_check_status_not_pruned, db_graph_pre, ginfo,
                       BINVERSION);

  //cleanup, before starting all over
  hash_table_free(&db_graph_pre);
  CU_ASSERT(db_graph_pre == NULL);

  //Now load the clean graph, so the "dirty" nodes are not even there
  db_graph_post = hash_table_new(number_of_bits_post, bucket_size,
                                 10, kmer_size);

  int num_cols_in_binary = -1;
  graph_info_initialise(ginfo);

  load_multicolour_binary_from_filename_into_graph(
    "../data/tempfiles_can_be_deleted/dump_cortex_var_graph.clean.ctx",
    db_graph_post, ginfo, &num_cols_in_binary);

  CU_ASSERT(num_cols_in_binary == NUMBER_OF_COLOURS);
  CU_ASSERT(ginfo->mean_read_length[0] == 0);
  CU_ASSERT(ginfo->total_sequence[0] == seq_read);

  // OK, finally we are ready for our real test. Load the singlecolour uncleaned
  // binary into colour1, but only load the bits that overlap the cleaned graph
  // (ie colour 0)
  load_single_colour_binary_data_from_filename_into_graph(
    "../data/tempfiles_can_be_deleted/dump_cortex_var_graph.singlecol.ctx",
    db_graph_post, ginfo, 
    false, 1, true, 0, false);



  BinaryKmer tmp_kmer1, tmp_kmer2;

  //all the kmers and their reverse complements from the cleaned graph - these should now be in colours 0 and 1
  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GAT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("ATC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element3 = hash_table_find(element_get_key(seq_to_binary_kmer("TCG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element4 = hash_table_find(element_get_key(seq_to_binary_kmer("CGA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element5 = hash_table_find(element_get_key(seq_to_binary_kmer("CGG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element6 = hash_table_find(element_get_key(seq_to_binary_kmer("CCG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element7 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element8 = hash_table_find(element_get_key(seq_to_binary_kmer("CCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element9 = hash_table_find(element_get_key(seq_to_binary_kmer("GGT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element10 = hash_table_find(element_get_key(seq_to_binary_kmer("ACC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element11 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element12 = hash_table_find(element_get_key(seq_to_binary_kmer("CCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element13 = hash_table_find(element_get_key(seq_to_binary_kmer("GGT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element14 = hash_table_find(element_get_key(seq_to_binary_kmer("ACC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element15 = hash_table_find(element_get_key(seq_to_binary_kmer("GTG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element16 = hash_table_find(element_get_key(seq_to_binary_kmer("CAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element17 = hash_table_find(element_get_key(seq_to_binary_kmer("TGT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element18 = hash_table_find(element_get_key(seq_to_binary_kmer("ACA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element19 = hash_table_find(element_get_key(seq_to_binary_kmer("GGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element20 = hash_table_find(element_get_key(seq_to_binary_kmer("GCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element21 = hash_table_find(element_get_key(seq_to_binary_kmer("GCT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element22 = hash_table_find(element_get_key(seq_to_binary_kmer("AGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element23 = hash_table_find(element_get_key(seq_to_binary_kmer("TAG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element24 = hash_table_find(element_get_key(seq_to_binary_kmer("CTA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element25 = hash_table_find(element_get_key(seq_to_binary_kmer("AGG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element26 = hash_table_find(element_get_key(seq_to_binary_kmer("CCT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);


  // kmers that should not be in the graph, in either colour 0 or 1
  dBNode* test_element27 = hash_table_find(element_get_key(seq_to_binary_kmer("GGA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element28 = hash_table_find(element_get_key(seq_to_binary_kmer("TCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element29 = hash_table_find(element_get_key(seq_to_binary_kmer("GAG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  // dBNode* test_element30 = hash_table_find(element_get_key(seq_to_binary_kmer("CTC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element31 = hash_table_find(element_get_key(seq_to_binary_kmer("AGT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element32 = hash_table_find(element_get_key(seq_to_binary_kmer("ACT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);


  CU_ASSERT(test_element1 != NULL);
  CU_ASSERT(test_element2 != NULL);
  CU_ASSERT(test_element1 == test_element2);
  CU_ASSERT(db_node_get_coverage(test_element1,0)==6);
  CU_ASSERT(db_node_get_coverage(test_element1,1)==6);
  CU_ASSERT(get_edge_copy(*test_element1,0)==get_edge_copy(*test_element1,1));

  CU_ASSERT(test_element3 != NULL);
  CU_ASSERT(test_element4 != NULL);
  CU_ASSERT(test_element3 == test_element4);
  CU_ASSERT(db_node_get_coverage(test_element3,0)==3);
  CU_ASSERT(db_node_get_coverage(test_element3,1)==3);
  CU_ASSERT(get_edge_copy(*test_element3,0)==get_edge_copy(*test_element3,1));

  CU_ASSERT(test_element5 != NULL);
  CU_ASSERT(test_element6 != NULL);
  CU_ASSERT(test_element5 == test_element6);
  CU_ASSERT(db_node_get_coverage(test_element5,0)==3);
  CU_ASSERT(db_node_get_coverage(test_element5,1)==3);
  CU_ASSERT(get_edge_copy(*test_element5,0)==get_edge_copy(*test_element5,1));

  CU_ASSERT(test_element7 != NULL);
  CU_ASSERT(test_element8 != NULL);
  CU_ASSERT(test_element7 == test_element8);
  CU_ASSERT(db_node_get_coverage(test_element7,0)==3);
  CU_ASSERT(db_node_get_coverage(test_element7,1)==3);
  CU_ASSERT(get_edge_copy(*test_element7,0)==get_edge_copy(*test_element7,1));

  CU_ASSERT(test_element9 != NULL);
  CU_ASSERT(test_element10 != NULL);
  CU_ASSERT(test_element9 == test_element10);
  CU_ASSERT(db_node_get_coverage(test_element9,0)==2);
  CU_ASSERT(db_node_get_coverage(test_element9,1)==2);
  CU_ASSERT(get_edge_copy(*test_element9,0)==get_edge_copy(*test_element9,1));

  CU_ASSERT(test_element11 != NULL);
  CU_ASSERT(test_element12 != NULL);
  CU_ASSERT(test_element11 == test_element12);
  CU_ASSERT(db_node_get_coverage(test_element11,0)==3);
  CU_ASSERT(db_node_get_coverage(test_element11,1)==3);
  CU_ASSERT(get_edge_copy(*test_element11,0)==get_edge_copy(*test_element11,1));

  CU_ASSERT(test_element13 != NULL);
  CU_ASSERT(test_element14 != NULL);
  CU_ASSERT(test_element13 == test_element14);
  CU_ASSERT(db_node_get_coverage(test_element13,0)==2);
  CU_ASSERT(db_node_get_coverage(test_element13,1)==2);
  CU_ASSERT(get_edge_copy(*test_element13,0)==get_edge_copy(*test_element13,1));

  CU_ASSERT(test_element15 != NULL);
  CU_ASSERT(test_element16 != NULL);
  CU_ASSERT(test_element15 == test_element16);
  CU_ASSERT(db_node_get_coverage(test_element15,0)==2);
  CU_ASSERT(db_node_get_coverage(test_element15,1)==2);
  CU_ASSERT(get_edge_copy(*test_element15,0)==get_edge_copy(*test_element15,1));

  CU_ASSERT(test_element17 != NULL);
  CU_ASSERT(test_element18 != NULL);
  CU_ASSERT(test_element17 == test_element18);
  CU_ASSERT(db_node_get_coverage(test_element17,0)==2);
  CU_ASSERT(db_node_get_coverage(test_element17,1)==2);
  CU_ASSERT(get_edge_copy(*test_element17,0)==get_edge_copy(*test_element17,1));

  CU_ASSERT(test_element19 != NULL);
  CU_ASSERT(test_element20 != NULL);
  CU_ASSERT(test_element19 == test_element20);
  CU_ASSERT(db_node_get_coverage(test_element19,0)==2);
  CU_ASSERT(db_node_get_coverage(test_element19,1)==2);
  CU_ASSERT(get_edge_copy(*test_element19,0)==get_edge_copy(*test_element19,1));

  CU_ASSERT(test_element21 != NULL);
  CU_ASSERT(test_element22 != NULL);
  CU_ASSERT(test_element21 == test_element22);
  CU_ASSERT(db_node_get_coverage(test_element21,0)==2);
  CU_ASSERT(db_node_get_coverage(test_element21,1)==2);
  CU_ASSERT(get_edge_copy(*test_element21,0)==get_edge_copy(*test_element21,1));

  CU_ASSERT(test_element23 != NULL);
  CU_ASSERT(test_element24 != NULL);
  CU_ASSERT(test_element23 == test_element24);
  CU_ASSERT(db_node_get_coverage(test_element23,0)==2);
  CU_ASSERT(db_node_get_coverage(test_element23,1)==2);
  CU_ASSERT(get_edge_copy(*test_element23,0)==get_edge_copy(*test_element23,1));

  CU_ASSERT(test_element25 != NULL);
  CU_ASSERT(test_element26 != NULL);
  CU_ASSERT(test_element25 == test_element26);
  CU_ASSERT(db_node_get_coverage(test_element25,0)==2);
  CU_ASSERT(db_node_get_coverage(test_element25,1)==2);
  CU_ASSERT(get_edge_copy(*test_element25,0)==get_edge_copy(*test_element25,1));

  //these nodes should just not be there
  CU_ASSERT(test_element27 == NULL);
  CU_ASSERT(test_element28 == NULL);
  CU_ASSERT(test_element29 == NULL);
  CU_ASSERT(test_element31 == NULL);
  CU_ASSERT(test_element32 == NULL);

  hash_table_free(&db_graph_post);
  graph_info_free(ginfo);

  CU_ASSERT(db_graph_post == NULL);
}


// DEV: remove test CU_ASSERT(num_cols_in_binary == NUMBER_OF_COLOURS); ?
void test_load_binversion5_binary()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER == 1)
  {
    int kmer_size = 3;
    int number_of_bits_pre = 4;
    int number_of_bits_post = 8;
    int bucket_size = 5;
    int max_retries = 10;

    dBGraph* db_graph_pre = hash_table_new(number_of_bits_pre, bucket_size,
                                           max_retries, kmer_size);
    dBGraph* db_graph_post;

    // Read sequence
    int fq_quality_cutoff = 20;
    int homopolymer_cutoff = 0;
    boolean remove_duplicates_se = false;
    char ascii_fq_offset = 33;
    int into_colour = 0;

    unsigned int files_loaded = 0;
    unsigned long long bad_reads = 0, dup_reads = 0;
    unsigned long long seq_read = 0, seq_loaded = 0;

    load_se_filelist_into_graph_colour(
      "../data/test/graph/test_dB_graph.falist",
      fq_quality_cutoff, homopolymer_cutoff,
      remove_duplicates_se, ascii_fq_offset,
      into_colour, db_graph_pre, 0,
      &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
      NULL, 0, &subsample_null);

    CU_ASSERT(seq_read == 16);
    CU_ASSERT(seq_loaded == 16);

    GraphInfo* ginfo = graph_info_alloc_and_init();

    graph_info_set_seq(ginfo, 0, seq_read);
    graph_info_set_mean_readlen(ginfo, 0, 5);

    db_graph_dump_binary(
      "../data/tempfiles_can_be_deleted/dump_cortex_var_graph.ctx.binversion5",
      &db_node_condition_always_true, db_graph_pre, ginfo, 5);

    hash_table_free(&db_graph_pre);

    CU_ASSERT(db_graph_pre == NULL);

    db_graph_post = hash_table_new(number_of_bits_post, bucket_size,
                                   10, kmer_size);

    int num_cols_in_binary = -1;

    graph_info_initialise(ginfo);
    long long seq_length_post = 
      load_multicolour_binary_from_filename_into_graph(
						       "../data/tempfiles_can_be_deleted/dump_cortex_var_graph.ctx.binversion5",
						       db_graph_post, ginfo, &num_cols_in_binary);
    
    CU_ASSERT(num_cols_in_binary == NUMBER_OF_COLOURS);
    CU_ASSERT(ginfo->mean_read_length[0] == 5);
    CU_ASSERT(ginfo->total_sequence[0] == seq_read);

    // Load_multicolour_binary_data_from_filename_into_graph returns total number
    // of unique kmers loaded, times kmer_length
    CU_ASSERT_EQUAL(seq_length_post,15);
    CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph_post),5);

    BinaryKmer tmp_kmer1, tmp_kmer2;

    // All the kmers and their reverse complements from the reads
    dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AAA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("TTT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element3 = hash_table_find(element_get_key(seq_to_binary_kmer("GGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element4 = hash_table_find(element_get_key(seq_to_binary_kmer("GCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element5 = hash_table_find(element_get_key(seq_to_binary_kmer("GCT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element6 = hash_table_find(element_get_key(seq_to_binary_kmer("AGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element7 = hash_table_find(element_get_key(seq_to_binary_kmer("TAG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element8 = hash_table_find(element_get_key(seq_to_binary_kmer("CTA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element9 = hash_table_find(element_get_key(seq_to_binary_kmer("AGG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element10 = hash_table_find(element_get_key(seq_to_binary_kmer("CCT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);


    // Kmers that should not be in the graph
    dBNode* test_element11 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element12 = hash_table_find(element_get_key(seq_to_binary_kmer("CCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element13 = hash_table_find(element_get_key(seq_to_binary_kmer("TAT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element14 = hash_table_find(element_get_key(seq_to_binary_kmer("ATA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element15 = hash_table_find(element_get_key(seq_to_binary_kmer("TAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element16 = hash_table_find(element_get_key(seq_to_binary_kmer("ATG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element17 = hash_table_find(element_get_key(seq_to_binary_kmer("TTG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element18 = hash_table_find(element_get_key(seq_to_binary_kmer("AAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element19 = hash_table_find(element_get_key(seq_to_binary_kmer("TGA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode* test_element20 = hash_table_find(element_get_key(seq_to_binary_kmer("TCA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);


    CU_ASSERT(test_element1 != NULL);
    CU_ASSERT(test_element2 != NULL);
    CU_ASSERT(test_element1 == test_element2);
    // Checking that we count KMER coverage.
    // ie this kmer occurs many times in the same read
    CU_ASSERT(db_node_get_coverage(test_element1,0)==6);

    CU_ASSERT(test_element3 != NULL);
    CU_ASSERT(test_element4 != NULL);
    CU_ASSERT(test_element3 == test_element4);
    CU_ASSERT(db_node_get_coverage(test_element3,0)==1);

    CU_ASSERT(test_element5 != NULL);
    CU_ASSERT(test_element6 != NULL);
    CU_ASSERT(test_element5 == test_element6);
    CU_ASSERT(db_node_get_coverage(test_element5,0)==1);

    CU_ASSERT(test_element7 != NULL);
    CU_ASSERT(test_element8 != NULL);
    CU_ASSERT(test_element7 == test_element8);
    CU_ASSERT(db_node_get_coverage(test_element7,0)==1);

    CU_ASSERT(test_element9 != NULL);
    CU_ASSERT(test_element10 != NULL);
    CU_ASSERT(test_element9 == test_element10);
    CU_ASSERT(db_node_get_coverage(test_element9,0)==1);

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


    // Check arrows
    Nucleotide base;

    // AAA -A-> AAA
    CU_ASSERT(db_node_has_precisely_one_edge(test_element1, forward,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Adenine);

    CU_ASSERT(db_node_has_precisely_one_edge(test_element1, reverse,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Thymine);

    //GGC -T-> GCT
    CU_ASSERT(db_node_has_precisely_one_edge(test_element3, reverse,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Thymine);


    //TAG -G-> AGG
    CU_ASSERT(db_node_has_precisely_one_edge(test_element7, reverse,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Guanine);

    //AGC -C-> GCC
    CU_ASSERT(db_node_has_precisely_one_edge(test_element6, forward,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Cytosine);

    //CCT -A-> CTA
    CU_ASSERT(db_node_has_precisely_one_edge(test_element10, reverse,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Adenine);

    // Add one extra arrows by had -- this breaks the graph -
    // it is only to check that the arrows get set correctly

    add_edges(test_element1,0,0x02);
    CU_ASSERT(db_node_has_precisely_one_edge(test_element1, forward,&base, 0)==false);

    CU_ASSERT(db_node_has_precisely_one_edge(test_element1, reverse,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Thymine);

    CU_ASSERT(db_node_edge_exist(test_element1, Adenine, forward, 0)==true);
    CU_ASSERT(db_node_edge_exist(test_element1, Cytosine, forward, 0)==true);
    CU_ASSERT(db_node_edge_exist(test_element1, Guanine, forward, 0)==false);
    CU_ASSERT(db_node_edge_exist(test_element1, Thymine, forward, 0)==false);

    add_edges(test_element1,0,0x20);
    CU_ASSERT(db_node_has_precisely_one_edge(test_element1, reverse,&base, 0)==false);

    CU_ASSERT(db_node_edge_exist(test_element1, Adenine, forward, 0)==true);
    CU_ASSERT(db_node_edge_exist(test_element1, Cytosine, forward, 0)==true);
    CU_ASSERT(db_node_edge_exist(test_element1, Guanine, forward, 0)==false);
    CU_ASSERT(db_node_edge_exist(test_element1, Thymine, forward, 0)==false);
    CU_ASSERT(db_node_edge_exist(test_element1, Adenine, reverse, 0)==false);
    CU_ASSERT(db_node_edge_exist(test_element1, Cytosine, reverse, 0)==true);
    CU_ASSERT(db_node_edge_exist(test_element1, Guanine, reverse, 0)==false);
    CU_ASSERT(db_node_edge_exist(test_element1, Thymine, reverse, 0)==true);

    hash_table_free(&db_graph_post);
    graph_info_free(ginfo);
    CU_ASSERT(db_graph_post == NULL);
  }

  // Now try the same thing with big kmers

  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER == 2)
  {
    int kmer_size = 33;
    int number_of_bits_pre = 10;
    int number_of_bits_post = 12;
    int bucket_size = 10;
    int max_retries = 10;

    dBGraph *db_graph_pre = hash_table_new(number_of_bits_pre, bucket_size,
                                           max_retries, kmer_size);

    // Read sequence
    int fq_quality_cutoff = 20;
    int homopolymer_cutoff = 0;
    boolean remove_duplicates_se = false;
    char ascii_fq_offset = 33;
    int into_colour = 0;

    unsigned int files_loaded = 0;
    unsigned long long bad_reads = 0, dup_reads = 0;
    unsigned long long seq_read = 0, seq_loaded = 0;

    load_se_filelist_into_graph_colour(
      "../data/test/graph/person2.falist",
      fq_quality_cutoff, homopolymer_cutoff,
      remove_duplicates_se, ascii_fq_offset,
      into_colour, db_graph_pre, 0,
      &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
      NULL, 0, &subsample_null);

    /*
    > 6 unique 33-mers
    TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACC
    > 13 unique 33-mers
    ACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAAC
    > 12 unique 33-mers
    GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
    */

    // Not bothering to set mean read length this time
    GraphInfo* ginfo = graph_info_alloc_and_init();

    graph_info_initialise(ginfo);
    graph_info_set_seq(ginfo, 0, seq_read);

    db_graph_dump_binary(
      "../data/tempfiles_can_be_deleted/dump_cortex_var_graph_2.ctx.binversion5",
      &db_node_condition_always_true, db_graph_pre, ginfo, 5);

    hash_table_free(&db_graph_pre);
    CU_ASSERT(db_graph_pre==NULL);

    dBGraph *db_graph_post = hash_table_new(number_of_bits_post, bucket_size,
                                            max_retries, kmer_size);

    graph_info_initialise(ginfo);

    int num_cols_in_binary = -1;

    load_multicolour_binary_from_filename_into_graph(
      "../data/tempfiles_can_be_deleted/dump_cortex_var_graph_2.ctx.binversion5",
      db_graph_post, ginfo, &num_cols_in_binary);


    CU_ASSERT(ginfo->mean_read_length[0]==0);//because we did not set it
    CU_ASSERT(num_cols_in_binary==NUMBER_OF_COLOURS);
    CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph_post),31);

    BinaryKmer tmp_kmer1, tmp_kmer2;

    //some the kmers and their reverse complements from the reads
    dBNode *test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("TAACCCTAACCCTAACCCTAACCCTAACCCTAA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode *test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode *test_element3 = hash_table_find(element_get_key(seq_to_binary_kmer("AACCCTAACCCTAACCCTAACCCTAACCCTAAC", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode *test_element4 = hash_table_find(element_get_key(seq_to_binary_kmer("GTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTT", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode *test_element5 = hash_table_find(element_get_key(seq_to_binary_kmer("ACCCTAACCCTAACCCTAACCCTAACCCTAACC", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode *test_element6 = hash_table_find(element_get_key(seq_to_binary_kmer("GGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGT", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode *test_element7 = hash_table_find(element_get_key(seq_to_binary_kmer("CCCTAACCCTAACCCTAACCCTAACCCTAACCC", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode *test_element8 = hash_table_find(element_get_key(seq_to_binary_kmer("GGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode *test_element9 = hash_table_find(element_get_key(seq_to_binary_kmer("CCTAACCCTAACCCTAACCCTAACCCTAACCCT", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode *test_element10 = hash_table_find(element_get_key(seq_to_binary_kmer("AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGG", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode *test_element11 = hash_table_find(element_get_key(seq_to_binary_kmer("CTAACCCTAACCCTAACCCTAACCCTAACCCTA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode *test_element12 = hash_table_find(element_get_key(seq_to_binary_kmer("TAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAG", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);

    //from read 2:
    dBNode *test_element13 = hash_table_find(element_get_key(seq_to_binary_kmer("ACCCTAACCCTAACCCTAACCCCTAACCCTAAC", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode *test_element14 = hash_table_find(element_get_key(seq_to_binary_kmer("GTTAGGGTTAGGGGTTAGGGTTAGGGTTAGGGT", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode *test_element15 = hash_table_find(element_get_key(seq_to_binary_kmer("CCCTAACCCTAACCCTAACCCCTAACCCTAACC", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
    dBNode *test_element16 = hash_table_find(element_get_key(seq_to_binary_kmer("GGTTAGGGTTAGGGGTTAGGGTTAGGGTTAGGG", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);

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

    CU_ASSERT(test_element11 != NULL);
    CU_ASSERT(test_element12 != NULL);
    CU_ASSERT(test_element11 == test_element12);

    CU_ASSERT(test_element13 != NULL);
    CU_ASSERT(test_element14 != NULL);
    CU_ASSERT(test_element13 == test_element14);

    CU_ASSERT(test_element15 != NULL);
    CU_ASSERT(test_element16 != NULL);
    CU_ASSERT(test_element15 == test_element16);

    // Now check arrows
    // Note we know that read1 forms a supernode that is a loop.
    Nucleotide base;

    // TAACCCTAACCCTAACCCTAACCCTAACCCTAA ----- C ----> AACCCTAACCCTAACCCTAACCCTAACCCTAAC
    CU_ASSERT(db_node_has_precisely_one_edge(test_element1, forward,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Cytosine);

    CU_ASSERT(db_node_has_precisely_one_edge(test_element1, reverse,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Guanine);

    // AACCCTAACCCTAACCCTAACCCTAACCCTAAC ----C  ----> ACCCTAACCCTAACCCTAACCCTAACCCTAACC
    CU_ASSERT(db_node_has_precisely_one_edge(test_element3, forward,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Cytosine);

    CU_ASSERT(db_node_has_precisely_one_edge(test_element3, reverse,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Adenine);

    // ACCCTAACCCTAACCCTAACCCTAACCCTAACC ---C -----> CCCTAACCCTAACCCTAACCCTAACCCTAACCC
    CU_ASSERT(db_node_has_precisely_one_edge(test_element5, forward,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Cytosine);

    CU_ASSERT(db_node_has_precisely_one_edge(test_element5, reverse,&base, 0)==true);
    CU_ASSERT_EQUAL(base,Thymine);

    // OK. Looks good.
    hash_table_free(&db_graph_post);
    graph_info_free(ginfo);
  }

  // Finally a test case which found a bug in binary read/write which no other
  // test case found
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER == 1)
  {
    int kmer_size = 17;
    int number_of_bits_pre = 10;
    int number_of_bits_post = 10;
    int bucket_size = 30;
    int max_retries = 10;

    dBGraph *db_graph_pre = hash_table_new(number_of_bits_pre, bucket_size,
                                           max_retries, kmer_size);

    // Read sequence
    int fq_quality_cutoff = 20;
    int homopolymer_cutoff = 0;
    boolean remove_duplicates_se = false;
    char ascii_fq_offset = 33;
    int into_colour = 0;

    unsigned int files_loaded = 0;
    unsigned long long bad_reads = 0, dup_reads = 0;
    unsigned long long seq_read = 0, seq_loaded = 0;

    load_se_filelist_into_graph_colour(
      "../data/test/graph/person3.falist",
      fq_quality_cutoff, homopolymer_cutoff,
      remove_duplicates_se, ascii_fq_offset,
      into_colour, db_graph_pre, 0,
      &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
      NULL, 0, &subsample_null);

    /*
    >read1 overlaps human chrom 1
    TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACC
    > read 2 overlaps human chrom 1
    ACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAAC
    > read 3 does not
    GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
    > read 3 does not
    GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
    > read 3 does not
    GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
    > read 4 does not, but has too low coverage
    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    */

    GraphInfo* ginfo = graph_info_alloc_and_init();

    graph_info_initialise(ginfo);
    graph_info_set_seq(ginfo, 0, seq_read);

    db_graph_dump_binary("../data/tempfiles_can_be_deleted/dump_cortex_var_graph_3.ctx.binversion5",
                         &db_node_condition_always_true, db_graph_pre, ginfo, 5);

    hash_table_free(&db_graph_pre);
    CU_ASSERT(db_graph_pre == NULL);

    dBGraph *db_graph_post = hash_table_new(number_of_bits_post, bucket_size,
                                            max_retries, kmer_size);

    graph_info_initialise(ginfo);

    int num_cols_in_binary = -1;

    load_multicolour_binary_from_filename_into_graph(
      "../data/tempfiles_can_be_deleted/dump_cortex_var_graph_3.ctx",
      db_graph_post, ginfo, &num_cols_in_binary);

    CU_ASSERT(num_cols_in_binary == NUMBER_OF_COLOURS);
    CU_ASSERT(ginfo->total_sequence[0] == seq_read);

    BinaryKmer tmp_kmer1, tmp_kmer2;

    // Now try to traverse a supernode. This is effectively a regression test for
    // a bug in graph/element.c: print_binary/read_binary
    dBNode *test_element1
      = hash_table_find(element_get_key(seq_to_binary_kmer("TAACCCTAACCCTAACC",
                                                           kmer_size, &tmp_kmer1),
                                        kmer_size, &tmp_kmer2), db_graph_post);

    // this node is in the middle of this supernode: CCCTAACCCTAACCCTAACCC. You
    // can't extend forward because there is one copy in the reads with a T
    // afterwards and one with a C afterwards. And you can' extend the other way
    // because the first kmer CCCTAACCCTAACCCTA has two arrows in. ie is preceded
    // by two different bases (A and T) in different reads

    CU_ASSERT(test_element1 != NULL);

    dBNode * path_nodes[100];
    Orientation path_orientations[100];
    Nucleotide path_labels[100];
    char path_string[100];
    int limit = 100;
    double avg_covg;
    Covg min_covg, max_covg;
    boolean is_cycle;

    int len
      = db_graph_supernode_for_specific_person_or_pop(
          test_element1, limit, &db_node_action_do_nothing,
          path_nodes, path_orientations, path_labels, path_string,
          &avg_covg, &min_covg, &max_covg, &is_cycle, db_graph_post,
          0);

    CU_ASSERT(len == 4);
    CU_ASSERT(is_cycle == false);

    CU_ASSERT(!strcmp(path_string, "AGGG") || !strcmp(path_string, "ACCC"));

    hash_table_free(&db_graph_post);
    graph_info_free(ginfo);

  }
}
