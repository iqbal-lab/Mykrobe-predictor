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
  element.h  - defines the interface for the de Bruijn graph node. The
  implementation is complemented by a hash table that stores every node indexed
  by kmers (BinaryKmers). 

  The element routines, ie the one required by hash_table/priority queue, are
  prefixed with element_

  The de Bruijn based routines are prefixed with db_node_
*/

#ifndef ELEMENT_H_
#define ELEMENT_H_

#include <stdio.h>
#include <inttypes.h>

#include "global.h"
#include "binary_kmer.h"

// type definitions

typedef char Edges;


typedef enum
  {
    unassigned   = 0,
    none         = 1,
    visited      = 2,
    pruned       = 3,
    exists_in_reference = 4,
    visited_and_exists_in_reference = 5,
    to_be_dumped = 6, //to be dumped as binary 
    read_start_forward = 7,//used when removing duplicate reads
    read_start_reverse = 8,//used when removing duplicate reads
    read_start_forward_and_reverse = 9,//used when removing duplicate reads  
    ignore_this_node = 10,
    in_desired_genotype = 11,
    special_visited = 12,
    special_none = 13,
    special_pruned = 14,
    fw_strand = 15,
    rv_strand = 16,
  } NodeStatus;


typedef enum
  {
    neither =0,
    one     =1,
    two     =2,
    both    =3,
  } AlleleStatus;


typedef struct{
  BinaryKmer kmer;
  Covg       coverage[NUMBER_OF_COLOURS];
  Edges      individual_edges[NUMBER_OF_COLOURS];
  char       status; // will cast a NodeStatus to char
  char       allele_status;
} Element;


typedef Element dBNode;
typedef BinaryKmer* Key;


Element* new_element();
void free_element(Element** element);
void element_assign(Element* e1, Element* e2);


// utility function for getting the desired edge char,
// by specifying the appropriate colour

// gets copy of edge
Edges get_edge_copy(const Element e, int colour);
Edges get_union_of_edges(Element e);
Edges element_get_colour_union_of_all_colours(const Element*);

Edges element_get_colour0(const Element* e);
Edges element_get_colour1(const Element* e);
Edges element_get_last_colour(const Element* e);

Covg element_get_covg_union_of_all_covgs(const dBNode*);
Covg element_get_covg_for_colourlist(const dBNode* e, int* colour_list,
                                     int list_len);

Covg element_get_covg_colour0(const dBNode* e);

#if NUMBER_OF_COLOURS > 1
Covg element_get_covg_colour1(const dBNode* e);
#endif

Covg element_get_covg_last_colour(const dBNode* e);

void add_edges(Element* e, int colour, Edges edge_char);
void set_edges(Element* e, int colour, Edges edge_char);
void reset_one_edge(Element* e, Orientation orientation, Nucleotide nucleotide, int colour);

int element_get_number_of_people_or_pops_containing_this_element(Element* e);


boolean element_smaller(Element e1, Element e2);
BinaryKmer* element_get_kmer(Element *e);
boolean element_is_key(Key key, Element e);
Key element_get_key(BinaryKmer*,short kmer_size, Key preallocated_key);
void element_initialise(Element *,Key, short kmer_size);
void element_initialise_kmer_covgs_edges_and_status_to_zero(Element * e);

void element_set_kmer(Element *e, Key kmer, short kmer_size);


// reverse orientation
Orientation opposite_orientation(Orientation);
Orientation db_node_get_orientation(BinaryKmer*, dBNode *, short kmer_size);

// add an edge between nodes -- NB: it adds both edges: forward and reverse
boolean db_node_add_edge(dBNode *, dBNode *, Orientation, Orientation, short kmer_size, int colour); 


// returns yes if the label defined by the nucleotide coresponds to an 
// outgoing edge in the side defined by the orientation.   
boolean db_node_edge_exist(dBNode *e, Nucleotide n, Orientation o, int colour);

// final argument f is a function that returns an Edge that is a function of the
// different colured edges in a node. e.g. we might be interested in the union
// of all the coloured edges in the graph, or just the colour/edge for the first
// person, or the union of all edges except that corresponding to the reference.
boolean db_node_edge_exist_within_specified_function_of_coloured_edges(
  dBNode *element, Nucleotide base,
  Orientation orientation, Edges (*f)(const Element* ));


//returns the label of the first outgoing edge -- leaving from the side 
//defined by orientation. 
boolean db_node_has_precisely_one_edge(dBNode *e, Orientation o, Nucleotide *n, int colour);


boolean db_node_has_precisely_one_edge_in_subgraph_defined_by_func_of_colours(
  dBNode * node, Orientation orientation, Nucleotide * nucleotide, 
	Edges (*get_colour)(const dBNode*) );

boolean db_node_has_precisely_one_edge_in_union_graph_over_all_people(
  dBNode * node, Orientation orientation, Nucleotide * nucleotide);

boolean db_node_has_precisely_two_edges(dBNode * node, Orientation orientation,
                                        Nucleotide *n1, Nucleotide *n2, int colour);

void db_node_reset_all_edges_for_all_people_and_pops_to_zero(Element* e);

//forgets about the edges
void db_node_reset_edges(dBNode *e, int colour);

void db_node_reset_edge(dBNode *e, Orientation o, Nucleotide n, int colour);



//TODO - maybe do not need to export this:
void db_node_reset_specified_edges(dBNode * node, Orientation orientation,
                                   Nucleotide nucleotide,
                                   void (*f)(dBNode*, Orientation, Nucleotide));


//check that the edges are 0's
boolean db_node_edges_reset(dBNode *node, int colour);

boolean db_node_check_status(dBNode *node, NodeStatus status);
boolean db_node_check_allele_status(dBNode * node, AlleleStatus status);
boolean db_node_check_status_not_pruned(dBNode *node);
boolean db_node_check_status_not_pruned_or_visited(dBNode *node);
boolean db_node_check_status_to_be_dumped(dBNode *node);
boolean db_node_check_status_is_not_visited(dBNode *node);


void db_node_set_status(dBNode *node,NodeStatus status);
void db_node_set_allele_status(dBNode *node, AlleleStatus status);
void db_node_set_status_to_none(dBNode *node);



//actions and conditions 

void db_node_action_set_status_none(dBNode *node);
void db_node_action_set_status_of_unpruned_to_none(dBNode *node);

void db_node_action_set_status_pruned(dBNode *node);
void db_node_action_set_status_visited(dBNode *node);
void db_node_action_set_status_visited_unless_marked_to_be_ignored(dBNode * node);
void db_node_action_set_status_special_visited(dBNode *node);
boolean db_node_check_status_special(dBNode*node);
void db_node_action_specialise_status(dBNode *node);
void db_node_action_unspecialise_status(dBNode *node);

void db_node_action_set_status_ignore_this_node(dBNode *node);

void db_node_action_set_status_visited_or_visited_and_exists_in_reference(dBNode *node);

void db_node_action_unset_status_visited_or_visited_and_exists_in_reference(dBNode *node);

void db_node_action_unset_status_visited_or_visited_and_exists_in_reference_or_ignore_this_node(dBNode * node);


void db_node_action_do_nothing(dBNode * node);

boolean db_node_check_status_none(dBNode * node);
boolean db_node_check_for_flag_ALL_OFF(dBNode * node);


boolean db_node_check_status_visited(dBNode * node);

boolean db_node_check_status_exists_in_reference(dBNode * node);

boolean db_node_check_status_visited_and_exists_in_reference(dBNode * node);

boolean db_node_check_status_is_not_exists_in_reference(dBNode * node);

boolean db_node_check_status_is_not_visited_or_visited_and_exists_in_reference(dBNode * node);

boolean db_node_condition_is_not_marked_to_be_ignored(dBNode* node);

boolean db_node_condition_always_true(dBNode* node);





void db_node_increment_coverage(dBNode* e, int colour);
void db_node_update_coverage(dBNode* e, int colour, long update);
Covg db_node_get_coverage_tolerate_null(const dBNode* e, int colour);
Covg db_node_get_coverage(const dBNode* e, int colour);
void db_node_set_coverage(dBNode* e, int colour, Covg covg);

Covg db_node_get_coverage_in_subgraph_defined_by_func_of_colours(
  const dBNode* e, Covg (*get_covg)(const dBNode*));



//check if node doesn't have any edges in a given orientation
boolean db_node_is_blunt_end(dBNode * node, Orientation orientation, int colour);

boolean db_node_is_blunt_end_in_subgraph_given_by_func_of_colours(
  dBNode * node, Orientation orientation,  Edges (*get_colour)(const dBNode*));


boolean db_node_is_this_node_in_this_person_or_populations_graph(dBNode* node, int colour);


boolean db_node_is_this_node_in_subgraph_defined_by_func_of_colours(
  dBNode* node, Edges (*get_colour)(const dBNode*));


//functions for binary format
void db_node_print_multicolour_binary(FILE *fp, dBNode *node);

void db_node_print_single_colour_binary_of_colour0(FILE *fp, dBNode *node);
void db_node_print_single_colour_binary_of_specified_colour(FILE *fp, dBNode *node, int colour);

//reading multicolour binaries
boolean db_node_read_multicolour_binary(FILE *fp, short kmer_size, dBNode *node,
                                        int num_colours_in_binary,
                                        int binversion_in_binheader);


// read a binary for an individual person, as dumped by the target "graph"
// load this data into given colour
boolean db_node_read_single_colour_binary(FILE *fp, short kmer_size, dBNode *node,
                                          int colour, int binversion_in_binheader);


boolean db_node_check_read_start(dBNode*node, Orientation ori);

void db_node_set_read_start_status(dBNode* node, Orientation ori);


#endif /* ELEMENT_H_ */
