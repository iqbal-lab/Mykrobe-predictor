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
  test_dB_graph_population.h
*/

#ifndef TEST_DB_GRAPH_POP_H_
#define TEST_DB_GRAPH_POP_H_

#include "global.h"

void test_hash_table_find();
void test_tip_clipping();
void test_pruning_low_coverage_nodes();
void test_get_perfect_path_in_one_colour() ;
void test_detect_and_smooth_bubble();
void test_db_graph_db_node_has_precisely_n_edges_with_status_in_one_colour();
void test_is_supernode_end();
//void test_get_N50();
void test_is_condition_true_for_all_nodes_in_supernode();
void test_getting_stats_of_how_many_indivduals_share_a_node();
void test_get_min_and_max_covg_of_nodes_in_supernode();
void test_db_graph_supernode_for_specific_person_or_pop();
void test_db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta();
void test_db_graph_make_reference_path_based_sv_calls();
void test_db_graph_make_reference_path_based_sv_calls_null_test_1();
void test_db_graph_make_reference_path_based_sv_calls_null_test_2();
void test_db_graph_make_reference_path_based_sv_calls_null_test_3();
void test_db_graph_make_reference_path_based_sv_calls_null_test_4();
void test_db_graph_make_reference_path_based_sv_calls_null_test_5();
void test_db_graph_make_reference_path_based_sv_calls_test_1();
void test_db_graph_make_reference_path_based_sv_calls_test_2();
void test_db_graph_make_reference_path_based_sv_calls_test_3();
void test_db_graph_make_reference_path_based_sv_calls_test_4();
void test_db_graph_make_reference_path_based_sv_calls_test_5();
void test_db_graph_make_reference_path_based_sv_calls_test_6();
void test_db_graph_make_reference_path_based_sv_calls_test_7();
void test_db_graph_make_reference_path_based_sv_calls_test_8();
void test_db_graph_make_reference_path_based_sv_calls_test_9();
void test_get_covg_of_nodes_in_one_but_not_other_of_two_arrays();
void test_apply_to_all_nodes_in_path_defined_by_fasta();
void test_does_this_path_exist_in_this_colour();
void test_dump_covg_distribution();

#endif /* TEST_DB_GRAPH_POP_H_ */
