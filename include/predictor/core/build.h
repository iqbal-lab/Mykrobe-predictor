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
  build.h - build the de Bruijn graph
*/

#ifndef BUILD_H_
#define BUILD_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <inttypes.h>

#include "element.h"
#include "open_hash/hash_table.h"
#include "seq.h"
#include "file_reader.h"

void set_all_coverages_to_zero(dBGraph* dbg, int colour);

unsigned long long build_unclean_graph(dBGraph* db_graph, StrBuf* path, boolean is_list, uint16_t kmer,
				       uint64_t* readlen_distrib, int readlen_distrib_len,
				       uint64_t* kmer_covg_array, int len_kmer_covg_array,
				       boolean only_load_pre_existing_kmers, int into_colour);


void db_graph_get_covg_distribution_array(dBGraph* db_graph, int colour,
					  boolean (*condition)(dBNode* elem),
					  uint64_t* kmer_covg_array, 
					  uint32_t len_kmer_covg_array);

//if we return -1 we have failed due to OOM
int choose_cleaning_threshold(uint64_t* kmer_covg_array, 
			      uint32_t len_kmer_covg_array,
			      int expected_depth);   //passing D_eff as integer - good enough


boolean clean_graph(dBGraph* db_graph, 
		    uint64_t* kmer_covg_array, uint32_t len, 
		    int d_eff, int max_expected_sup);

                                                                          
dBNode * db_graph_get_next_node_for_specific_person_or_pop(dBNode * current_node, 
							   Orientation current_orientation, 
							   Orientation * next_orientation,
							   Nucleotide edge, 
							   Nucleotide * reverse_edge,
							   dBGraph * db_graph, 
							   int index);

//This function does not  check that it there is such an edge in the specified person/colour - but it does check if the target node is in the specific person.
//if you want to be sure the dge exists in that colour, then check it before calling this function
//The last argument allows you to apply the operation to some subgraph - eg you might take the unuiion of colours 2 and 3, or of all colours.
//If for example you wanted to get the next node in the graph irrespective of colour, get_colour would return the union (bitwise AND) of all edges in a node.
dBNode * db_graph_get_next_node_in_subgraph_defined_by_func_of_colours(dBNode * current_node, 
								       Orientation current_orientation, 
								       Orientation * next_orientation,
								       Nucleotide edge, 
								       Nucleotide * reverse_edge,
								       dBGraph * db_graph, 
								       Edges (*get_colour)(const dBNode*)
								       );


int db_graph_db_node_clip_tip_with_orientation_for_specific_person_or_pop(dBNode * node, 
									  Orientation orientation, 
									  int limit,
									  void (*node_action)(dBNode * node),dBGraph * db_graph, int index);
int db_graph_db_node_clip_tip_for_specific_person_or_pop(dBNode * node, int limit,
							 void (*node_action)(dBNode * node),
							 dBGraph * db_graph, int index);

int db_graph_db_node_clip_tip_with_orientation_in_subgraph_defined_by_func_of_colours(dBNode * node, 
										      Orientation orientation, 
										      int limit,
										      void (*node_action)(dBNode * node),
										      dBGraph * db_graph, 
										      Edges (*get_colour)(const dBNode*),
										      void (*apply_reset_to_specific_edge_in_colour)(dBNode*, Orientation, Nucleotide),
										      void (*apply_reset_to_colour)(dBNode*)
										      );


// clip a tip in the graph (the tip starts in node)
// limit is max length for tip
// node_action is applied to all the elements in the tip
// returns the length of the tip (0 means no length)
int db_graph_db_node_clip_tip_in_subgraph_defined_by_func_of_colours(dBNode * node, int limit,
								     void (*node_action)(dBNode * node),
								     dBGraph * db_graph, 
								     Edges (*get_colour)(const dBNode*),
								     void (*apply_reset_to_specific_edge_in_colour)(dBNode*, Orientation, Nucleotide),
								     void (*apply_reset_to_colour)(dBNode*)
								     );


boolean db_graph_db_node_prune_low_coverage(dBNode * node, Covg coverage,
					    void (*node_action)(dBNode * node),
					    dBGraph * db_graph, 
					    Covg (*sum_of_covgs_in_desired_colours)(const Element *), 
					    Edges (*get_edge_of_interest)(const Element*),
					    void (*apply_reset_to_specified_edges)(dBNode*, Orientation, Nucleotide),
					    void (*apply_reset_to_specified_edges_2)(dBNode*) );

boolean db_graph_db_node_prune_without_condition(dBNode * node, 
						 void (*node_action)(dBNode * node),
						 dBGraph * db_graph, 
						 Edges (*get_edge_of_interest)(const Element*),
						 void (*apply_reset_to_specified_edges)(dBNode*, Orientation, Nucleotide),
						 void (*apply_reset_to_specified_edges_2)(dBNode*) );


boolean db_graph_db_node_prune_low_coverage_ignoring_colours(dBNode * node, Covg coverage,
							     void (*node_action)(dBNode * node),
							     dBGraph * db_graph);

int db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop(
									 dBNode *node, 
									 Orientation orientation, 
									 int limit, 
									 Nucleotide fst_nucleotide,
									 void (*node_action)(dBNode *node),
									 dBNode **path_nodes, 
									 Orientation *path_orientations, 
									 Nucleotide *path_labels,
									 char *seq, 
									 double *avg_coverage, 
									 Covg *min_coverage, 
									 Covg *max_coverage,
									 boolean *is_cycle, 
									 dBGraph *db_graph, 
									 int index);

int db_graph_get_perfect_path_with_first_edge_in_subgraph_defined_by_func_of_colours(
										     dBNode *node, 
										     Orientation orientation, 
										     int limit, 
										     Nucleotide fst_nucleotide,
										     void (*node_action)(dBNode *node), 
										     dBNode **path_nodes,
										     Orientation *path_orientations, 
										     Nucleotide *path_labels, 
										     char *seq,
										     double *avg_coverage, 
										     Covg *min_coverage, 
										     Covg *max_coverage,
										     boolean *is_cycle, 
										     dBGraph *db_graph, 
										     Edges (*get_colour)(const dBNode*),
										     Covg (*get_covg)(const dBNode*));





// perfect path -- no conflict no cycle -- returns length
// path: node_0 edge_0 node_1 edge_1 ... node_n-1 edge_n-1 node_n
// path_nodes is a n+1 array from 0..n with all the nodes in the path
// path_orientations is n+1 array from 0..n with the orientations of the node in the path
// path labels is n array from 0..n-1 with all the labels for the edges (see above)
// node_action only applied to internal nodes (not node_0 and node_n)
// seq is a string with all the labels concatenated (NB: it doesn't contain the kmer in the first node)
// avg_coverage, min_coverge, max_coverge -> refers to the internal nodes only

int db_graph_get_perfect_path_for_specific_person_or_pop(
							 dBNode *node, 
							 Orientation orientation, 
							 int limit,
							 void (*node_action)(dBNode *node), dBNode **path_nodes,
							 Orientation *path_orientations, Nucleotide *path_labels, char *seq,
							 double *avg_coverage, Covg *min_coverage, Covg *max_coverage,
							 boolean *is_cycle, dBGraph *db_graph, int index);


// perfect path -- no conflict no cycle -- returns length
// path: node_0 edge_0 node_1 edge_1 ... node_n-1 edge_n-1 node_n
// path_nodes is a n+1 array from 0..n with all the nodes in the path
// path_orientations is n+1 array from 0..n with the orientations of the node in the path
// path labels is n array from 0..n-1 with all the labels for the edges (see above)
// node_action only applied to internal nodes (not node_0 and node_n)
// seq is a string with all the labels concatenated (NB: it doesn't contain the kmer in the first node)
// avg_coverage, min_coverge, max_coverge -> refers to the internal nodes only

int db_graph_get_perfect_path_in_subgraph_defined_by_func_of_colours(
  dBNode *node, Orientation orientation, int limit,
  void (*node_action)(dBNode *node), dBNode **path_nodes,
  Orientation *path_orientations, Nucleotide *path_labels, char *seq,
  double *avg_coverage, Covg *min_coverage, Covg *max_coverage,
  boolean *is_cycle, dBGraph *db_graph, 
  Edges (*get_colour)(const dBNode*), Covg (*get_covg)(const dBNode*));




// it returns the supernode containing 'node'  
// string has to support limit+1 (+1 as you need a space for the \0 at the end)
// node_action has to be idempotent as it can be applied to the same node twice!!
// supernode_str returns the string made of the labels of the path (doesn't include first kmer). 
// returns the length of supernode

int db_graph_supernode_for_specific_person_or_pop(
  dBNode *node, int limit, void (*node_action)(dBNode *node), 
  dBNode **path_nodes, Orientation *path_orientations,
  Nucleotide *path_labels, char *supernode_str, 
  double *avg_coverage, Covg *min, Covg *max, boolean *is_cycle, 
  dBGraph *db_graph, int index);


// it returns the supernode containing 'node'  
// string has to support limit+1 (+1 as you need a space for the \0 at the end)
// node_action has to be idempotent as it can be applied to the same node twice!!
// supernode_str returns the string made of the labels of the path (doesn't include first kmer). 
// returns the length of supernode 
// last two arguments mean you find supernode in a subgraph of the main graph defined by the edges as returned by get_colour.
int db_graph_supernode_in_subgraph_defined_by_func_of_colours(
  dBNode * node,int limit,void (*node_action)(dBNode * node), 
  dBNode * * path_nodes, Orientation * path_orientations,
  Nucleotide * path_labels, char * supernode_str, 
  double * avg_coverage, Covg * min, Covg * max,
  boolean * is_cycle, dBGraph * db_graph, 
  Edges (*get_colour)(const dBNode*),
  Covg (*get_covg)(const dBNode*));




// identical code to db_graph_supernode_for_specific_person_or_pop but this returns the index of the query node (first argument) in th supernode array
// will make a significant performance gain compared with getting the path_nodes array back and then searching it

int db_graph_supernode_returning_query_node_posn_for_specific_person_or_pop(
  dBNode *node,int limit,void (*node_action)(dBNode *node), 
  dBNode **path_nodes, Orientation *path_orientations, Nucleotide *path_labels,
  char *supernode_str, double *avg_coverage, Covg *min, Covg *max,
  boolean *is_cycle, int *query_node_posn, dBGraph *db_graph, int index);


void db_graph_clip_tips_for_specific_person_or_pop(dBGraph * db_graph, int index);


void db_graph_clip_tips_in_subgraph_defined_by_func_of_colours(dBGraph * db_graph,
							       Edges (*get_colour)(const dBNode*),
							       void (*apply_reset_to_specific_edge_in_colour)(dBNode*, Orientation, Nucleotide),
							       void (*apply_reset_to_colour)(dBNode*));

void apply_reset_to_specific_edge_in_union_of_all_colours(dBNode* node, Orientation or, Nucleotide nuc);

void apply_reset_to_all_edges_in_union_of_all_colours(dBNode* node );


void db_graph_clip_tips_in_union_of_all_colours(dBGraph* db_graph);


void traverse_hash_collecting_sums(
  void (*f)(Element *, Covg*, Covg*, Covg*, Covg*, Covg*, Covg*, Covg*),
  HashTable * hash_table, 
  Covg* pgf, Covg* cox, Covg* data1, Covg* data2, Covg* data3, Covg* data4, Covg* data5);



 //NOTE THIS LOOKS AT MAX COVG ON SUPERNODE
//if the node has covg <= coverage (arg2) and its supernode has length <=kmer+1 AND all the interiro nodes of the supernode have this low covg, then
//prune the whole of the interior of the supernode
//note the argument supernode_len is set to -1 if the passed in node has status != none
// returns true if the supernode is pruned, otherwise false
// if returns false, cannot assume path_nodes etc contain anything useful
boolean db_graph_remove_supernode_containing_this_node_if_looks_like_induced_by_error(
  dBNode* node, Covg coverage, dBGraph * db_graph, int max_expected_sup_len,
  Covg (*sum_of_covgs_in_desired_colours)(const Element *), 
  Edges (*get_edge_of_interest)(const Element*), 
  void (*apply_reset_to_specified_edges)(dBNode*, Orientation, Nucleotide), 
  void (*apply_reset_to_specified_edges_2)(dBNode*),
  dBNode** path_nodes, Orientation* path_orientations, Nucleotide* path_labels, 
  char* supernode_str, int* supernode_len);






//ESTIMATES NUMBER OF READS
//depth of covg = bp of covg/genome length
boolean db_graph_remove_supernode_containing_this_node_if_more_likely_error_than_sampling(dBNode* node, int num_haploid_chroms, 
											  double total_depth_of_covg, int read_len, 
											  double error_rate_per_base,
											  int max_length_allowed_to_prune,
											  dBGraph* db_graph, int max_expected_sup,
											  Covg (*sum_of_covgs_in_desired_colours)(const Element *), 
											  Edges (*get_edge_of_interest)(const Element*), 
											  void (*apply_reset_to_specified_edges)(dBNode*, Orientation, Nucleotide), 
											  void (*apply_reset_to_specified_edges_2)(dBNode*),
											  dBNode** path_nodes, 
											  Orientation* path_orientations, 
											  Nucleotide* path_labels,
											  char* supernode_string, int* supernode_len,
											  Covg covg_thresh);



// traverse graph. At each node, if covg <= arg1, get its supernode. If that supernode length is <= kmer-length, and ALL interior nodes have covg <= arg1 
// then prune the node, and the interior nodes of the supernode.
// returns the number of pruned supernodes
long long db_graph_remove_errors_considering_covg_and_topology(
  Covg coverage, dBGraph * db_graph,
  Covg (*sum_of_covgs_in_desired_colours)(const Element *), 
  Edges (*get_edge_of_interest)(const Element*), 
  void (*apply_reset_to_specified_edges)(dBNode*, Orientation, Nucleotide), 
  void (*apply_reset_to_specified_edges_2)(dBNode*),
  int max_expected_sup);




void db_graph_traverse_with_array(void (*f)(HashTable*, Element *, int**, int, int),
                                  HashTable * hash_table, int** array,
                                  int length_of_array, int index);

void db_graph_traverse_with_array_of_uint64(void (*f)(HashTable*, Element *, uint64_t*, uint32_t, int),
                                            HashTable * hash_table,
                                            uint64_t* array,
                                            uint32_t length_of_array, int colour);



void db_graph_get_covg_distribution(char* filename, dBGraph* db_graph, 
                                    int index, boolean (*condition)(dBNode* elem));



//TODO - get rid of this

dBNode* db_graph_get_next_node_in_supernode_for_specific_person_or_pop(dBNode* node, Orientation orientation, Orientation* next_orientation, int index, dBGraph* db_graph);


// This assumes the given fasta has already been loaded into the graph. All we do here is get the sequence of nodes
// that correspond to the sequence of bases in the file.

// Note that we have to handle the N's we expect to see in a reference fasta. The most reliable thing to do, for these purposes, is to put a NULL node
// in the array for every kmer in the fasta that contains an N

//Note we need to know whether we expect a new entry this time, AND whether there was a new fasta entry last time:
// When about to load a new set of nodes, decide whether to preload seq with the last kmer of the previous set. Do this iff you are not expecting the start of a new entry.
// Only tricky bit is determining the index of the last kmer within seq (from the last time around), when you want to move it to the start of seq just before loading a new set of nodes.
// If last time, we preloaded seq with a kmer at the start (indices 0 to k-1), and we loaded N nodes, then the last kmer in seq is from indices N to N+k-1
// If last time we did not preload seq, and we loaded N nodes, then the last kmer is from N-1 to N+k-2.
// This boils down to knowing whether the last time was a new fasta entry


// Note what happens as you repeatedly call this, with number_of_nodes_to_load = length_of_arrays/2. The first time, you load nodes into the back (right hand) 
// end of the array, and the front end is empty.
//    The next time, these are pushed to the front, and new ones are loaded into the back.
//  the first time you call this, expecting_new_fasta_entry=true,  last_time_was_not_start_of_entry=false
//  the next time,                expecting_new_fasta_entry=false, last_time_was_not_start_of_entry=false,
//  after that                    expecting_new_fasta_entry=false, last_time_was_not_start_of_entry=true

// returns the number of nodes loaded into the array
int db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta
                                                  (FILE* chrom_fptr, 
						   int number_of_nodes_to_load, int number_of_nodes_loaded_last_time,
						   int length_of_arrays, 
						   dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, char* path_string,
						   Sequence* seq, KmerSlidingWindow* kmer_window, 
						   boolean expecting_new_fasta_entry, boolean last_time_was_not_start_of_entry, 
						   dBGraph* db_graph ) ;


int db_node_addr_cmp(const void *a, const void *b);


void get_coverage_from_array_of_nodes(dBNode** array, int length,
                                      Covg *min_coverage, Covg *max_coverage,
                                      double* avg_coverage, Covg *mode_coverage,
                                      double* percent_nodes_having_modal_value,
                                      int index);

void get_coverage_from_array_of_nodes_in_subgraph_defined_by_func_of_colours(
  dBNode** array, int length,  Covg *min_coverage, Covg *max_coverage,
  double* avg_coverage, Covg *mode_coverage,
  double *percent_nodes_having_modal_value,
  Edges (*get_colour)(const dBNode*),
  Covg (*get_covg)(const dBNode*));





boolean does_this_path_exist_in_this_colour(dBNode** array_nodes, Orientation* array_orientations,  int len, Edges (*get_colour)(const dBNode*), dBGraph* db_graph );


//check all edges in graph
long long db_graph_health_check(boolean fix, dBGraph * db_graph);

//clean off edges that point nowhere - caused when you dump a subgraph
long long db_graph_clean_orphan_edges(dBGraph * db_graph);





#endif
