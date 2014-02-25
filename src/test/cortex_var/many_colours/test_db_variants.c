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
  test_db_variants.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include <CUnit.h>

// cortex_var headers
#include "element.h"
#include "file_reader.h"
#include "dB_graph.h"
#include "dB_graph_population.h"
#include "cmd_line.h"
#include "graph_info.h"
#include "db_differentiation.h"

void test_count_reads_on_allele_in_specific_colour()
{
  dBNode* node_array[10];
  boolean too_short=false;

  int i;
  for (i=0; i<10; i++)
    {
      node_array[i]=new_element();
    }

  // start simple - all nodes same covg
  for (i=0; i<10; i++)
    {
      db_node_set_coverage(node_array[i], 0,5);
    }
  
  //number of edges in allele of 10 nodes is 9
  int colour0 = 0;
  CU_ASSERT(count_reads_on_allele_in_specific_colour(node_array, 9, colour0,&too_short)==5);
  CU_ASSERT(too_short==false);

  //add one jump
  db_node_set_coverage(node_array[1], colour0, 10);
  db_node_set_coverage(node_array[2], colour0, 10);
  CU_ASSERT(count_reads_on_allele_in_specific_colour(node_array, 9, colour0,&too_short)==10);
  CU_ASSERT(too_short==false);
  
  //add a drop back to 5
  db_node_set_coverage(node_array[3], colour0, 5);
  CU_ASSERT(count_reads_on_allele_in_specific_colour(node_array, 9, colour0,&too_short)==10);
  CU_ASSERT(too_short==false);

  //add a second jump, which is an isolated node with higher covg - ignore it!
  db_node_set_coverage(node_array[5], colour0, 8);
  CU_ASSERT(count_reads_on_allele_in_specific_colour(node_array, 9, colour0,&too_short)==10);
  CU_ASSERT(too_short==false);


  //add a third jump, not isolated
  db_node_set_coverage(node_array[6], colour0, 8);
  db_node_set_coverage(node_array[7], colour0, 8);
  CU_ASSERT(count_reads_on_allele_in_specific_colour(node_array, 9, colour0,&too_short)==13);
  CU_ASSERT(too_short==false);


  
  //what if allele is short?
  CU_ASSERT(count_reads_on_allele_in_specific_colour(node_array, 1, colour0,&too_short)==0);
  CU_ASSERT(too_short==true);


  for (i=0; i<10; i++)
    {
      free_element(&node_array[i]);
    }

}


void test_get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller()
{
  /*
  // Looks like this test doesn't rely on a particular kmer size
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }
  */

  // recall in get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller,
  // theta/2 is the expected depth of covg on the allele - Dl_i/2R in the language of our paper

  double ret_het;
  double ret_hom1;
  double ret_hom2;

  //1. Simple het
  ret_het = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(het, 0.01, 10,10, 20,20);
  ret_hom1 = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(hom_one, 0.01, 10,10, 20,20);
  ret_hom2 = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(hom_other, 0.01, 10,10, 20,20);
  CU_ASSERT( (ret_het>ret_hom1) && (ret_het>ret_hom2) );
  CU_ASSERT(ret_hom1 == ret_hom2);

  //2. simple hom
  ret_het = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(het, 0.01, 10,0, 20,20);
  ret_hom1 = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(hom_one, 0.01, 10,0, 20,20);
  ret_hom2 = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(hom_other, 0.01, 10,0, 20,20);
  CU_ASSERT( (ret_het<ret_hom1) && (ret_hom1>ret_hom2) );

  ret_het = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(het, 0.01, 0,10, 20,20);
  ret_hom1 = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(hom_one, 0.01, 0,10, 20,20);
  ret_hom2 = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(hom_other, 0.01, 0,10, 20,20);
  CU_ASSERT( (ret_het<ret_hom2) && (ret_hom2>ret_hom1) );

  //3. if we have coverage one one allele only, but high covg, we call confident hom. but as covg decreases, confidence decreases
  double t1,t2,t3,t4,t5;
  t1 = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(hom_one, 0.01, 10,0, 20,20);
  t2 = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(hom_one, 0.01, 8,0, 20,20);
  t3 = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(hom_one, 0.01, 6,0, 20,20);
  t4 = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(hom_one, 0.01, 4,0, 20,20);
  t5 = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(hom_one, 0.01, 1,0, 20,20);
  CU_ASSERT(t5<t4);
  CU_ASSERT(t4<t3);
  CU_ASSERT(t3<t2);
  CU_ASSERT(t2<t1);



  //4. Confidence in a het increases as we increase covg on one of the alleles up from 1
  t1 = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(het, 0.01, 10,1, 20,20);
  t2 = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(het, 0.01, 10,3, 20,20);
  t3 = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(het, 0.01, 10,5, 20,20);
  t4 = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(het, 0.01, 10,10, 20,20);
  CU_ASSERT(t4>t3);
  CU_ASSERT(t3>t2);
  CU_ASSERT(t2>t1);
  
  //5. as we increase the sequencing error rate, confidence in hets with low covg on one allele go down
  t1 = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(het,     0.04, 10,1, 20,20);
  t2 = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(hom_one, 0.04, 10,1, 20,20);
  t3 = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(het,     0.01, 10,1, 20,20);
  t4 = get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(hom_one, 0.01, 10,1, 20,20);

  //higher sequencing error gives lower CONFIDENCE (gap between likelihood of het and hom models)
  CU_ASSERT(t1-t2 < t3-t4);


}





