/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 * **********************************************************************
 *
 * This file is part of myKrobe.
 *
 * myKrobe is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * myKrobe is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with myKrobe.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */
/*
  test_genotype_known.c 
*/

// system headers
#include <stdlib.h>
#include <limits.h>

// third party headers
#include <CUnit.h>
#include <Basic.h>
#include <string_buffer.h>

#include "build.h"
#include "element.h"
#include "seq.h"
#include "open_hash/hash_table.h"
#include "genotyping_known.h"

void test_get_next_mutation_allele_info()
{
  
  uint16_t kmer_size = 31;
  int number_of_bits = 10;
  int bucket_size = 100;
  int max_retries = 10;


  dBGraph *db_graph= hash_table_new(number_of_bits, bucket_size,
				    max_retries, kmer_size);

  uint64_t* kmer_covg_array = calloc(100, sizeof(uint64_t));
  int max_read_length = 120;
  uint64_t* readlen_array = calloc(max_read_length, sizeof(uint64_t));

  StrBuf* list = strbuf_create("../data/test/myKrobe/predictor/mutations/sample1.fa.list");
  unsigned long long  num_bases = build_unclean_graph(db_graph, 
						      list, true,
						      kmer_size,
						      readlen_array, max_read_length,
						      kmer_covg_array, 100,
						      false, 0);

  FILE* fp = fopen("../data/test/myKrobe/predictor/mutations/some_snps1.fa", "r");
  if (fp==NULL)
    {
      die("Cannot open this file: ../data/test/myKrobe/predictor/mutations/some_snps1.fa");
    }
  
  ResVarInfo* rvi = alloc_and_init_res_var_info();


  //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence");
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);
  
  //We are going to load all the bases into a single sliding window 
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      die("Failed to malloc kmer sliding window");
    }
  

  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-db_graph->kmer_size+1));
  if (kmer_window->kmer==NULL)
    {
      die("Failed to malloc kmer_window->kmer");
    }
  kmer_window->nkmers=0;
  

  CovgArray* working_ca = alloc_and_init_covg_array(500);
  //end of intialisation 
	  
	  
  //create file readers
  int file_reader_fasta(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;
    if (new_entry == false){
      offset = db_graph->kmer_size;
    }
    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);
    
    return ret;
  }

  dBNode** array_nodes = (dBNode**) malloc(sizeof(dBNode*)*500);
  Orientation* array_or =(Orientation*)  malloc(sizeof(Orientation)*500);
  if ( (array_nodes==NULL) || (array_or==NULL))
    {
      die("Cannot alloc array of nodes or of orientations");
    }
  StrBuf* temp_rid = strbuf_new();
  StrBuf* temp_mut = strbuf_new();
  StrBuf* temp_gene = strbuf_new();
  int ignore = 1;
  get_next_mutation_allele_info(fp, db_graph, rvi,
				seq, kmer_window,
				&file_reader_fasta,
				array_nodes, array_or, working_ca, max_read_length,
				temp_rid, temp_mut, temp_gene,
				ignore, ignore);

  CU_ASSERT(rvi->susceptible_allele.median_covg==9);
  CU_ASSERT(rvi->susceptible_allele.min_covg==9);
  CU_ASSERT(rvi->susceptible_allele.percent_nonzero==100);

  CU_ASSERT(rvi->resistant_alleles[0].median_covg==0);
  CU_ASSERT(rvi->resistant_alleles[0].min_covg==0);
  CU_ASSERT(rvi->resistant_alleles[0].percent_nonzero==0);



  get_next_mutation_allele_info(fp, db_graph, rvi,
				seq, kmer_window,
				&file_reader_fasta,
				array_nodes, array_or, working_ca, max_read_length,
				temp_rid, temp_mut, temp_gene,
				ignore, ignore);


  CU_ASSERT(rvi->susceptible_allele.median_covg==0);
  CU_ASSERT(rvi->susceptible_allele.min_covg==0);
  CU_ASSERT(rvi->susceptible_allele.percent_nonzero==0);

  CU_ASSERT(rvi->resistant_alleles[0].median_covg==2);
  CU_ASSERT(rvi->resistant_alleles[0].min_covg==2);
  CU_ASSERT(rvi->resistant_alleles[0].percent_nonzero==100);


  get_next_mutation_allele_info(fp, db_graph, rvi,
				seq, kmer_window,
				&file_reader_fasta,
				array_nodes, array_or, working_ca, max_read_length,
				temp_rid, temp_mut, temp_gene,
				ignore, ignore);


  CU_ASSERT(rvi->susceptible_allele.median_covg==0);
  CU_ASSERT(rvi->susceptible_allele.min_covg==0);
  CU_ASSERT(rvi->susceptible_allele.percent_nonzero==0);

  CU_ASSERT(rvi->resistant_alleles[0].median_covg==0);
  CU_ASSERT(rvi->resistant_alleles[0].min_covg==0);
  CU_ASSERT(rvi->resistant_alleles[0].percent_nonzero==0);



  strbuf_free(temp_rid);
  strbuf_free(temp_mut);
  strbuf_free(temp_gene);
  free_res_var_info(rvi);
  free(array_nodes);
  free(array_or);
  free_covg_array(working_ca);
  free(kmer_window->kmer);
  free(kmer_window);
  free_sequence(&seq);
  free(kmer_covg_array);
  free(readlen_array);
  strbuf_free(list);
  hash_table_free(&db_graph);
}



