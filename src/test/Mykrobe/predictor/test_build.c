/*
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 *
 */
/*
  test_build.c 
*/

// system headers
#include <stdlib.h>
#include <limits.h>

// third party headers
#include <CUnit.h>
#include <Basic.h>
#include <string_buffer.h>

// cortex_var headers
#include "build.h"
#include "element.h"
#include "seq.h"
#include "open_hash/hash_table.h"

void test_build_unclean_graph()
{
  
  uint16_t kmer_size = 31;
  int number_of_bits = 19;
  int bucket_size = 100;
  int max_retries = 10;

  dBGraph *db_graph= hash_table_new(number_of_bits, bucket_size,
				    max_retries, kmer_size);

  uint64_t* kmer_covg_array = calloc(100, sizeof(uint64_t));
  uint64_t* readlen_array = calloc(120, sizeof(uint64_t));

  StrBuf* list = strbuf_create("../data/test/myKrobe/predictor/test1.bam.list");
  unsigned long long  num_bases = build_unclean_graph(db_graph, 
						      list, 
						      true,
						      kmer_size,
						      readlen_array, 120,
						      kmer_covg_array, 100,
						      false, 0);

  
  CU_ASSERT(num_bases == 346933015);
  CU_ASSERT(kmer_covg_array[1]==4100361);
  CU_ASSERT(kmer_covg_array[2]==67309);
  CU_ASSERT(kmer_covg_array[3]==3482);

  free(kmer_covg_array);
  free(readlen_array);
  strbuf_free(list);
  hash_table_free(&db_graph);
}



