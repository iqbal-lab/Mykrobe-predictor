/*
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
  test_seq_error_estimation.c
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <CUnit.h>
#include <Basic.h>

// cortex_var headers
#include "maths.h"
#include "element.h"
#include "file_reader.h"
#include "seq_error_rate_estimation.h"
#include "test_seq_error_estimation.h"

void test_estimate_seq_error_rate_for_one_colour_from_snp_allele_fasta()
{
  if (NUMBER_OF_COLOURS < 11)
    {
      warn("This test requires 11 colours. Please recompile\n");
      return;
    }

  //first set up the hash/graph
  int kmer_size = 7;
  int number_of_bits = 10;
  int bucket_size = 4;
  int max_retries = 10;

  dBGraph *db_graph = hash_table_new(number_of_bits, bucket_size,
                                     max_retries, kmer_size);

  // Genome length 89bp, basically load 20 copies of the same genome

  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_read = 0, seq_loaded = 0;

  unsigned long readlen_distrib_arrlen = 200;
  unsigned long *readlen_distrib
    = (unsigned long*) malloc(sizeof(unsigned long) * readlen_distrib_arrlen);

  if(readlen_distrib == NULL)
  {
    die("Unable to malloc array to hold readlen distirbution! Exiting.");
  }

  unsigned long i;
  for(i = 0; i < readlen_distrib_arrlen; i++)
  {
    readlen_distrib[i] = 0;
  }

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/seq_error_estimation/list_sample1.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    readlen_distrib, readlen_distrib_arrlen, &subsample_null);


  // Initialise a graph_info object
  GraphInfo* ginfo = graph_info_alloc_and_init();

  unsigned long mean_readlen = calculate_mean_ulong(readlen_distrib,
                                                    readlen_distrib_arrlen);

  graph_info_update_mean_readlen_and_total_seq(ginfo, 0, mean_readlen, seq_loaded);  
  
  long double default_err_rate = 0.01;
  estimate_seq_error_rate_from_snps_for_each_colour(
    "../data/test/pop_graph/seq_error_estimation/test_sample1.falist", ginfo,
    db_graph, -1, //89, //genome_size not used
    default_err_rate, NULL);

  CU_ASSERT_DOUBLE_EQUAL(ginfo->seq_err[0],0.0, 0.0001);



  hash_table_free(&db_graph);
  graph_info_free(ginfo);
  free(readlen_distrib);
}



void test_estimate_seq_error_rate_for_one_colour_from_snp_allele_fasta_test2()
{
  if (NUMBER_OF_COLOURS < 11)
    {
      warn("This test assumes 11 colours. Please recompile\n");
      return;
    }

  // First set up the hash/graph
  int kmer_size = 31;
  int number_of_bits = 13;
  int bucket_size = 100;
  int max_retries = 10;

  dBGraph *db_graph = hash_table_new(number_of_bits, bucket_size,
                                     max_retries, kmer_size);

  // Genome length 770bp, basically load 100 copies of the same genome,
  // plus two read errors
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned int files_loaded = 0;
  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_read = 0, seq_loaded = 0;

  unsigned long readlen_distrib_arrlen = 5000;
  unsigned long *readlen_distrib
    = (unsigned long*) malloc(sizeof(unsigned long) * readlen_distrib_arrlen);

  if(readlen_distrib == NULL)
  {
    die("Unable to malloc array to hold readlen distirbution! Exiting.");
  }

  unsigned long i;
  for(i = 0; i < readlen_distrib_arrlen; i++)
  {
    readlen_distrib[i] = 0;
  }

  load_se_filelist_into_graph_colour(
    "../data/test/pop_graph/seq_error_estimation/list_sample2.falist",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, 0, // 0 => falist/fqlist; 1 => colourlist
    &files_loaded, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    readlen_distrib, readlen_distrib_arrlen, &subsample_null);

  // Initialise a graph_info object
  GraphInfo* ginfo = graph_info_alloc_and_init();

  unsigned long mean_readlen = calculate_mean_ulong(readlen_distrib,
                                                    readlen_distrib_arrlen);

  graph_info_update_mean_readlen_and_total_seq(ginfo, 0, mean_readlen, seq_loaded);  
  
  long double default_err_rate = 0.01;

  estimate_seq_error_rate_from_snps_for_each_colour(
    "../data/test/pop_graph/seq_error_estimation/test_sample2.falist",
    ginfo, db_graph, -1, //770, //genome_size not used
    default_err_rate, NULL);


  CU_ASSERT_DOUBLE_EQUAL(ginfo->seq_err[0],0.0909, 0.01);


  hash_table_free(&db_graph);
  graph_info_free(ginfo);
  free(readlen_distrib);
}
