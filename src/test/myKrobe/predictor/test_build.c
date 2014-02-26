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
						      kmer_size,
						      readlen_array, 120,
						      kmer_covg_array, 100);

  
  CU_ASSERT(num_bases == 346933015);
  CU_ASSERT(kmer_covg_array[1]==4100361);
  CU_ASSERT(kmer_covg_array[2]==67309);
  CU_ASSERT(kmer_covg_array[3]==3482);

  free(kmer_covg_array);
  free(readlen_array);
  strbuf_free(list);
  hash_table_free(&db_graph);
}



