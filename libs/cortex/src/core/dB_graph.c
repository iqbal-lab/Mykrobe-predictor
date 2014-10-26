/*
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@tgac.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * The MIT License (MIT)
 * Copyright (c) 2009-2014 <Z. Iqbal and M. Caccamo>
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
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
#include "build.h"


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














