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

//#include "test_dB_graph.h"
#include "test_dB_graph_node.h"
#include "test_dB_graph_population.h"
#include "test_pop_element.h"
#include "test_pop_load_and_print.h"
#include "test_pop_supernode_consensus.h"
#include "test_db_genotyping.h"
#include "test_db_variants.h"
#include "test_model_selection.h"
#include <test_file_reader.h>
#include <test_db_complex_genotyping.h>
#include <test_genome_complexity.h>
#include <test_seq_error_estimation.h>
#include "test_error_correction.h"
#include <CUnit.h>
#include <Basic.h>

int  main()
{

  CU_pSuite pPopGraphSuite = NULL;

  /* initialize the CUnit test registry */
  if (CUE_SUCCESS!=CU_initialize_registry())
    return CU_get_error();

  /* add a suite to the registry */
  pPopGraphSuite = CU_add_suite("Test main functionality of Cortex", NULL, NULL);
  if (NULL == pPopGraphSuite) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  /* add the tests to the suites */



  if (NULL == CU_add_test(pPopGraphSuite, "Test element - get edge copy", test_get_edge_copy)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Test assignment operator for element", test_element_assign)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Test setting and checking of element status",test_element_status_set_and_checks )) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Test element - get coverage of node for specific person", test_get_coverage)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Test element -increment coverage for different people", test_increment_coverage)) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pPopGraphSuite, "Test function checking number of edges to nodes with specified status",test_db_graph_db_node_has_precisely_n_edges_with_status_in_one_colour    )) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Sanity test for hash_table_find - adding kmers to graph and finding them", test_hash_table_find)) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pPopGraphSuite, "Test dumping and reloading of multicolour binary", test_dump_load_sv_trio_binary)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Test dump Binary format version5, and reload", test_load_binversion5_binary)) {
    CU_cleanup_registry();
    return CU_get_error();
  }




  if (NULL == CU_add_test(pPopGraphSuite, "Test loading of a singlecolour binary",test_load_singlecolour_binary )) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Test loading of a single-colour binary where nodes are only loaded if they already exist in a predefined \"clean\" colour, and if so, only taking the edges that overlap both colours", test_loading_binary_data_iff_it_overlaps_a_fixed_colour  )) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Test that coverage of kmers are correcty counted as reads are loaded from a file",test_coverage_is_correctly_counted_on_loading_from_file )) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  
  // comment out while debugging if (NULL == CU_add_test(pPopGraphSuite, "Regression test: integer overflow and dumping of covergae distribution does not segfault",test_dump_covg_distribution )) //{
  //    CU_cleanup_registry();
  //  return CU_get_error();
  // }
  


  if (NULL == CU_add_test(pPopGraphSuite, "Test getting sliding windows, breaking when kmer is not in graph - internal function",test_getting_sliding_windows_where_you_break_at_kmers_not_in_db_graph )) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Test getting sliding windows breaking when kner or edge is not in graph - internal function",test_get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph )) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pPopGraphSuite, "Test function for getting length distribution of filtered/effective reads after N's, low quality, PCR duplicates and homopolymers have been cut or filtered",test_getting_readlength_distribution )) {
    CU_cleanup_registry();
    return CU_get_error();
  }



  if (NULL == CU_add_test(pPopGraphSuite, "Test checking if a path lies in a (function of) colours", test_does_this_path_exist_in_this_colour  )) 
    {
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pPopGraphSuite, "Test dumping of cleaned fasta by aligning against a cleaned graph",test_dumping_of_clean_fasta )) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Test removal of PCR duplicate reads from paired-end fastq",test_loading_of_paired_end_reads_removing_duplicates )) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Test removal of PCR duplicate reads from single-ended fastq", test_loading_of_single_ended_reads_removing_duplicates)) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pPopGraphSuite, "Test loading of three binaries as dumped by graph into sv_trio as separate people", test_load_individual_binaries_into_sv_trio )) {
    CU_cleanup_registry();
    return CU_get_error();
  }


  if (NULL == CU_add_test(pPopGraphSuite, "Test utility function for getting supernode containing a given node",   test_db_graph_supernode_for_specific_person_or_pop)) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pPopGraphSuite, "Test that can identify supernode ends",   test_is_supernode_end)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Test that can get the min and max coverage of a supernode",   test_get_min_and_max_covg_of_nodes_in_supernode)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Test function for checking if condition is true for all nodes in supernode",   test_is_condition_true_for_all_nodes_in_supernode)) {
    CU_cleanup_registry();
    return CU_get_error();
  }


  if (NULL == CU_add_test(pPopGraphSuite, "Test error correction via topology of graph - tip clipping", test_tip_clipping )) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Test removal of low-coverage nodes",  test_pruning_low_coverage_nodes)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  //  if (NULL == CU_add_test(pPopGraphSuite, "Test (currently unused) function for smoothing bubbles",  test_detect_and_smooth_bubble   )) {
  //CU_cleanup_registry();
  //return CU_get_error();
  //}
  if (NULL == CU_add_test(pPopGraphSuite, "Test utility function for applying some other function to all nodes in a path defined by a fasta",   test_apply_to_all_nodes_in_path_defined_by_fasta)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
    if (NULL == CU_add_test(pPopGraphSuite, "Test get_perfect_path (maximal path without branches) in single colour of graph", test_get_perfect_path_in_one_colour     )) {
    CU_cleanup_registry();
    return CU_get_error();
  }


  if (NULL == CU_add_test(pPopGraphSuite, "Loading two people, one fasta each, in same population, and print supernodes. \n\tEach person has 2 reads, and we have manually determined what the supernodes should be", test_load_two_people_in_same_populations_and_print_separately_their_supernodes)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Load three people in one population, each consists of one toy chromosome.\n\t Find their individual supernodes correctly", test_take_three_people_each_with_one_read_and_find_variants))
    {
      CU_cleanup_registry();
      return CU_get_error();
    }
  if (NULL == CU_add_test(pPopGraphSuite, "Load two people sharing an Alu, differing by sequene before and after the Alu. Find supernodes.", test_take_two_people_sharing_an_alu_and_find_supernodes))
    {
      CU_cleanup_registry();
      return CU_get_error();
    }
  if (NULL == CU_add_test(pPopGraphSuite, "Load two people in one population, and test that given a kmer, we can find the first node in the supernode in which that kmer lies.", test_find_first_node_in_supernode))
    {
      CU_cleanup_registry();
      return CU_get_error();
    }

  if (NULL == CU_add_test(pPopGraphSuite, "Load two people in one population, and test that given a node in a person's graoh, you can find the next node in the supernode.", test_find_next_node_in_supernode))
    {
      CU_cleanup_registry();
      return CU_get_error();
    }

  if (NULL == CU_add_test(pPopGraphSuite, "Load two people in one population, and test that can pull out required subsection of a supernode.", test_correctly_find_subsection_of_supernode))
    {
      CU_cleanup_registry();
      return CU_get_error();
    }
    if (NULL == CU_add_test(pPopGraphSuite, "Load two people in one population, and test that for each separately, can find the best sub_supernode, requiring certain people_coverage", test_find_best_subsection_of_supernode_with_just_two_people))
    {
     CU_cleanup_registry();
     return CU_get_error();
    }

    if (NULL == CU_add_test(pPopGraphSuite, "Check that correctly get stats on how many kmers are shared by 1,2,3,... people in a population", test_getting_stats_of_how_many_indivduals_share_a_node))
    {
     CU_cleanup_registry();
     return CU_get_error();
    }
    
    if (NULL == CU_add_test(pPopGraphSuite, "Unit test of utilty function to produce array of nodes corresponding to a given path through graph, as specified by a fasta file",  test_load_seq_into_array))
      {
    	CU_cleanup_registry();
    	return CU_get_error();
     }
    

    if (NULL == CU_add_test(pPopGraphSuite, "Unit test of wrapper function for above utility function",  test_db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }
    if (NULL == CU_add_test(pPopGraphSuite, "Test utility function to get next read from file and return it as an array of nodes",  test_align_next_read_to_graph_and_return_node_array))
      {
    	CU_cleanup_registry();
    	return CU_get_error();
     }



    if (NULL == CU_add_test(pPopGraphSuite, "Simple null test of new algorithm to find variants using trusted paths AND supernodes. Take a person whose genome is a single read, and a reference genome which is identical.\n\tAlgorithm should find no variants!",  test_db_graph_make_reference_path_based_sv_calls_null_test_1))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }
    if (NULL == CU_add_test(pPopGraphSuite, "Simple null test of new algorithm to find variants using trusted paths AND supernodes. Take a person whose genome is a single Alu, and a reference genome which is identical.\n\t Algorithm should find no variants.",  test_db_graph_make_reference_path_based_sv_calls_null_test_2))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }

    if (NULL == CU_add_test(pPopGraphSuite, "Another null test of the Path Divergence/supernode SV caller. Reference and individual are identical again.\n\t In this case they are both equal to an Alu, some N's, and then the same Alu again. Algorithm should find no variants.\n\t Note this is testing something new - since the reference is essentially a repeat of an Alu, there is the risk that our algorithm would map the 3prime anchor onto the wrong copy of the repeat,\n\t and falsely call an insertion between the start of the first copy, and the end of the second.",  test_db_graph_make_reference_path_based_sv_calls_null_test_3))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }
   if (NULL == CU_add_test(pPopGraphSuite, "Another null test of the Path Divergence/supernode SV caller. Reference and individual are identical - this time they are both equal to around 10kb taken from the start of Human chromosome 1.\n\t Algorithm should find nothing.",  test_db_graph_make_reference_path_based_sv_calls_null_test_4))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }
   if (NULL == CU_add_test(pPopGraphSuite, "A targeted null test of the Path Divergence/supernode SV caller. The reference is a tandem repeat of about 36 bases; the individual has an Alu inserted between the two copies of the repeat. \n\tSince the supernode in the individual has one copy of the repeat, then the Alu, and then stops, the algorithm should fail to attach a 3prime anchor, and fail to call the insertion.",  
			   test_db_graph_make_reference_path_based_sv_calls_null_test_5))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }

   if (NULL == CU_add_test(pPopGraphSuite, "First simple (positive) test of the Path Divergence SV caller. Reference and individual are identical except for a single SNP. \n\tTest that caller finds the SNP, correctly gives the coordinates in the reference, and gives the sequence of the 5prime, and 3prime anchors, and the two alternative branches (reference and individual).",  test_db_graph_make_reference_path_based_sv_calls_test_1))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }

   if (NULL == CU_add_test(pPopGraphSuite, "Test of the Path Divergence SV caller. Reference has two copies of the same Alu (genome is just those two Alu's plus some N's - nothing else) and the individual differs by a single base in one copy of the Alu. \n\tAlgorithm should call the SNP despite being in the Alu. \n\t(Further, this is an example where we walk the supernode in the opposite direction to that in which it lies in our array - this therefore covers the other main code path in our implementation.)",  test_db_graph_make_reference_path_based_sv_calls_test_2))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }

   if (NULL == CU_add_test(pPopGraphSuite, "Test of the Path Divergence SV caller. Reference is a short sequence, and individual has a 2-base deletion.\n\tCaller should spot this and correctly identify the coordinates in the reference, the flanking regions, and the sequence for the two branches",  test_db_graph_make_reference_path_based_sv_calls_test_3))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }
   if (NULL == CU_add_test(pPopGraphSuite, "Test of the Path Divergence SV caller. Reverse the roles of reference and individual in the previous test. Caller should spot the 2 base insertion\n\tIdentify the coordinates of variant start, and the relevant flank and branch sequences",  test_db_graph_make_reference_path_based_sv_calls_test_4))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }
   if (NULL == CU_add_test(pPopGraphSuite, "Test of the Path Divergence SV caller. Reference is a single supernode, and individual has a complete Alu inserted within the supernode. Find the insertion.\n\t",  test_db_graph_make_reference_path_based_sv_calls_test_5))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }
   if (NULL == CU_add_test(pPopGraphSuite, "Test of the Path Divergence SV caller. Reverse roles of reference/individual in previous test, and call deletion of Alu.",  test_db_graph_make_reference_path_based_sv_calls_test_6))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }
   if (NULL == CU_add_test(pPopGraphSuite, "Test of the Path Divergence SV caller. Reference has an Alu inserted within another Alu. Individual lacks the insertion.\n\tCall the Alu insertion and correctly identify the coordinates of the variant, flanking sequences, and sequences of the two branches",  test_db_graph_make_reference_path_based_sv_calls_test_7))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }
   if (NULL == CU_add_test(pPopGraphSuite, "Test of the Path Divergence SV caller. Reference is about 10kb ofHuman chromosome 1, with about 1 kb of Human chromosome 2 inserted in the middle. Individual lacks the 1kb insertion.\n\tCall the deletion, and correctly find maximal flanking regions, and identfy the sequence of the two branches",  test_db_graph_make_reference_path_based_sv_calls_test_8))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }

   if (NULL == CU_add_test(pPopGraphSuite, "Test of the Path Divergence SV caller. Identical to above test, but with 30kb of Human chromosome 12 added before sequence of above test\n\tThis test is purely to check that the code which loads chunks of data from the Path Divergence fasta bit by bit, works properly, and allows accurate determining of the correct variant location in the trusted fasta file.",  test_db_graph_make_reference_path_based_sv_calls_test_9))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }

   if (NULL == CU_add_test(pPopGraphSuite, "Test utility function for finding coverages of nodes that lie  on one allele but not the other, in a variant.",  test_get_covg_of_nodes_in_one_but_not_other_of_two_arrays))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }
   if (NULL == CU_add_test(pPopGraphSuite, "Test utility function for calculating read-arrivals on an allele",  test_count_reads_on_allele_in_specific_colour))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }
   if (NULL == CU_add_test(pPopGraphSuite, "Test calculation of genotype likelihoods for sites called by Bubble Caller", test_get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller ))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }

//   if (NULL == CU_add_test(pPopGraphSuite, "Test classification of bubble as variant or repeat, comparing Hardy-Weinberg model with a repeat model", test_get_log_bayesfactor_varmodel_over_repeatmodel ))
  //    {
//	CU_cleanup_registry();
//	return CU_get_error();
  //    }



   if (NULL == CU_add_test(pPopGraphSuite, "Test utility function for comparing overlap of alleles with each other and selves", test_initialise_multiplicities_of_allele_nodes_wrt_both_alleles ))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }



   if (NULL == CU_add_test(pPopGraphSuite, "Test reading of variant call output (full flank file)", test_read_next_variant_from_full_flank_file))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }


   if (NULL == CU_add_test(pPopGraphSuite, "Test reading of variant call output (full flank file), when one branch < k long", test_read_next_variant_from_full_flank_file_2))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }


  
   if (NULL == CU_add_test(pPopGraphSuite, "Test reading of variant call output (full flank file), when both branches < k long", test_read_next_variant_from_full_flank_file_3))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }




   if (NULL == CU_add_test(pPopGraphSuite, "Test reading of variant call output (full flank file), when both branches < k long and zero-length 3p flank", test_read_next_variant_from_full_flank_file_4))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }




   //if (NULL == CU_add_test(pPopGraphSuite, "Test function for estimating genome complexity - first test", test_count_reads_where_snp_makes_clean_bubble1 ))
   //   {
//	CU_cleanup_registry();
//	return CU_get_error();
  //    }






   /*

   if (NULL == CU_add_test(pPopGraphSuite, "Regression test case 1 - genotyping of bubble with one branch <k and one long branch", regression_test_1_single_bubble_call_one_allele_shorter_than_k_one_very_long ))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }
 

     if (NULL == CU_add_test(pPopGraphSuite, "Regression test case 2 - genotyping of a PD call", regression_test_2_genotyping_of_PD_SNP_call))
    {
  	CU_cleanup_registry();
  	return CU_get_error();
    }

  

   if (NULL == CU_add_test(pPopGraphSuite, "Test algorithm for genotyping of complex site at a simple site (repeat for different coverages and sequencing error rates, simulating real coverage coording to our model, 100 iterations each time)", test_calc_log_likelihood_of_genotype_with_complex_alleles1 ))
      {
	CU_cleanup_registry();
	return CU_get_error();
	}


  

   
   if (NULL == CU_add_test(pPopGraphSuite, "Test genotyping of complex sites - test with two HLA-B alleles", test_calc_log_likelihood_of_genotype_with_complex_alleles2 ))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }


   if (NULL == CU_add_test(pPopGraphSuite, "Test genotyping of complex sites - test with two HLA-B alleles, using 1-net and 2-net error model", test_calc_log_likelihood_of_genotype_with_complex_alleles3 ))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }



   if (NULL == CU_add_test(pPopGraphSuite, "Test estimation of sequencing error rate from fasta of pairs of SNP alleles (first allele known not to be present, second known to be homozygous present) - very basic test", test_estimate_seq_error_rate_for_one_colour_from_snp_allele_fasta ))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }

   if (NULL == CU_add_test(pPopGraphSuite, "Test estimation of sequencing error rate from fasta of pairs of SNP alleles (first allele known not to be present, second known to be homozygous present) - another test", test_estimate_seq_error_rate_for_one_colour_from_snp_allele_fasta_test2 ))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }

   */



   if (NULL == CU_add_test(pPopGraphSuite, "Test utility function for mutating specific bases in a string", test_base_mutator ))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }

   if (NULL == CU_add_test(pPopGraphSuite, "Test utility function for error-correcting reads to match a trusted graph", test_fix_end_if_unambiguous ))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }

   if (NULL == CU_add_test(pPopGraphSuite, "Test utility function for checking kmers and qualities in a read that might need correction", test_get_first_good_kmer_and_populate_qual_array ))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }

   if (NULL == CU_add_test(pPopGraphSuite, "Test error correction of fastq files against a graph", test_error_correct_file_against_graph))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }

   if (NULL == CU_add_test(pPopGraphSuite, "Test error correction utility function for greedy random walk", test_take_n_greedy_random_steps))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }

   if (NULL == CU_add_test(pPopGraphSuite, "Test error correction can output corrected reads in positive strand direction if desired", 
			   test_reverse_comp_according_ref_pos_strand))
      {
	CU_cleanup_registry();
	return CU_get_error();
      }



 

  /* Run all tests using the CUnit Basic interface */
  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
 
  CU_cleanup_registry();
  return CU_get_error();


}

