/*
 * Copyrightf 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 *
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
#include "build.h"

void test_get_species_info()
{
  // Creating the graph of the sample
  uint16_t kmer_size = 31;
  int number_of_bits = 10;
  int bucket_size = 100;
  int max_retries = 10;
  dBGraph *db_graph= hash_table_new(number_of_bits, bucket_size,
            max_retries, kmer_size);

  int max_gene_len = 10000;
  uint64_t* kmer_covg_array = calloc(150, sizeof(uint64_t));
  uint64_t* readlen_array = calloc(max_gene_len, sizeof(uint64_t));

  StrBuf* list = strbuf_create("../data/test/myKrobe/predictor/species_assignment/species_ref_list");
  unsigned long long  num_bases = build_unclean_graph(db_graph, 
						      list, true,
						      kmer_size,
						      readlen_array, max_gene_len,
						      kmer_covg_array, 150,
						      false,0);

  StrBuf* install_dir = strbuf_create("../");

  Staph_species species_assigned = get_species(db_graph,max_gene_len, install_dir, 0,0);
  CU_ASSERT(species_assigned == Aureus);

  free(readlen_array);
  strbuf_free(list);
  strbuf_free(install_dir);
  free(kmer_covg_array);
  hash_table_free(&db_graph);
}
