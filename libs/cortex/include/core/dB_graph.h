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
  dB_graph.h 

  all the routines as prefixed with db_graph
*/

#ifndef DB_GRAPH_H_
#define DB_GRAPH_H_

#include <stdio.h>
  
#include "global.h"
#include "open_hash/hash_table.h"

typedef HashTable dBGraph;

//pays no attention to whether there is an edge joining current_node to the node you would get by adding this nucleotide.
//just checksto see if such a node is in the graph
dBNode * db_graph_get_next_node(dBNode * current_node, Orientation current_orientation, 
				Orientation * next_orientation,
				Nucleotide edge, Nucleotide * reverse_edge,dBGraph * db_graph);




//Functions applying to whole graph

int db_graph_clip_tip(dBNode * node, int limit,dBGraph * db_graph);

void db_graph_set_all_visited_nodes_to_status_none(dBGraph* hash_table);


#endif /* DB_GRAPH_H_ */
