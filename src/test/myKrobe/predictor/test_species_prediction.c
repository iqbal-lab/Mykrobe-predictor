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
  test_species_prediction.c 
*/

#include "test_species_prediction.h"
#include "species.h"
// system headers
#include <stdlib.h>
#include <limits.h>

// third party headers
#include <CUnit.h>
#include <Basic.h>
#include <string_buffer.h>

void test_get_species_info()
{
  // Creating the graph of the sample
  uint16_t kmer_size = 31;
  int number_of_bits = 10;
  int bucket_size = 100;
  int max_retries = 10;
  dBGraph *db_graph= hash_table_new(number_of_bits, bucket_size,
            max_retries, kmer_size);

  int max_gene_len = 1500;
  uint64_t* kmer_covg_array = calloc(150, sizeof(uint64_t));
  uint64_t* readlen_array = calloc(max_gene_len, sizeof(uint64_t));

  StrBuf* list = strbuf_create("../data/test/myKrobe/predictor/gene_presence/sample1.fa.list");
  unsigned long long  num_bases = build_unclean_graph(db_graph, 
                  list, 
                  kmer_size,
                  readlen_array, max_gene_len,
                  kmer_covg_array, 150);

  Staph_species species_known = Aureus;
  Staph_species species_assigned = get_species(db_graph,max_gene_len);



  CU_ASSERT(species_assigned == species_known);

  free(readlen_array);
  strbuf_free(list);
  free(kmer_covg_array);
  hash_table_free(&db_graph);
}