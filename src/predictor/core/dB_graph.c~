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
  dB_graph.c - implementation
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// cortex_var headers
#include "binary_kmer.h"
#include "element.h"
#include "dB_graph.h"
#include "dB_graph_population.h"


//This gets the next node in the graph, and does not care about whether it
//is there for any specific person or population
//it just looks to see if it is there at all, and if so, gets it
// Also does not care if there is an edge joining this node to the one it returns.
// You simply tell it to add a certain nucleotide, see if there is such a node, and if so, return it.

dBNode * db_graph_get_next_node(dBNode * current_node, Orientation current_orientation, 
			       Orientation * next_orientation,
			       Nucleotide edge, Nucleotide * reverse_edge,dBGraph * db_graph){

  BinaryKmer local_copy_of_kmer;
  binary_kmer_assignment_operator(local_copy_of_kmer, current_node->kmer);
  
  BinaryKmer tmp_kmer;
  dBNode * next_node=NULL;

  // after the following line tmp_kmer and rev_kmer are pointing to the same B Kmer
  BinaryKmer* rev_kmer = binary_kmer_reverse_complement(&local_copy_of_kmer,db_graph->kmer_size, &tmp_kmer);

  
  if (current_orientation == reverse){   
    *reverse_edge = binary_kmer_get_last_nucleotide(&local_copy_of_kmer);
    binary_kmer_assignment_operator(local_copy_of_kmer,*rev_kmer);
  }
  else{
    *reverse_edge = binary_kmer_get_last_nucleotide(rev_kmer);
  }


  binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&local_copy_of_kmer, edge, db_graph->kmer_size);

   //get node from table
  next_node = hash_table_find(element_get_key(&local_copy_of_kmer,db_graph->kmer_size, &tmp_kmer),db_graph);

  if (next_node != NULL){
    *next_orientation = db_node_get_orientation(&local_copy_of_kmer,next_node,db_graph->kmer_size);
  }
  else
    {
      // debug
      char tmpzamseq[db_graph->kmer_size+1];
      warn("Cannot find %s so get a NULL node\n", binary_kmer_to_seq(&tmp_kmer, db_graph->kmer_size, tmpzamseq));
    }

  return next_node;

}



void db_graph_set_all_visited_nodes_to_status_none(dBGraph* hash_table)
{
  hash_table_traverse(&db_node_set_status_to_none, hash_table);
}














