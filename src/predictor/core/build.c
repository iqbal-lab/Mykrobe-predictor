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
  build.c - build the de Bruijn graph
*/

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
#include "build.h"
#include "graph_info.h"


void set_all_coverages_to_zero(dBGraph* dbg, int colour)
{
    void reset_covg(dBNode * node)
    {
      if (node!=NULL)
	{
	  node->coverage[colour]=0;
	}
    }
    hash_table_traverse(&reset_covg, dbg);
}

double estimate_err_rate(StrBuf* path, boolean is_list)
{
  return 0.005;
}

//if boolean is_list==true, then path=list of fastq (or bams)
//if boolean is_list==false, then path is a fastq file  (or bam)
unsigned long long build_unclean_graph(dBGraph* db_graph, 
				       StrBuf* path, 
				       boolean is_list, 
				       uint16_t kmer,
				       uint64_t* readlen_distrib, 
				       int readlen_distrib_len,
				       uint64_t* kmer_covg_array, 
				       int len_kmer_covg_array,
				       boolean only_load_pre_existing_kmers,
				       int into_colour,
				       boolean (*subsample_function)() )
{
  int ascii_fq_offset = 33;
  int qual_thresh = 10;
  int homopolymer_cutoff=0;

  unsigned int num_files_loaded=0;
  unsigned long long total_bad_reads=0; 
  unsigned long long total_dup_reads=0;
  unsigned long long total_bases_read=0; 
  unsigned long long total_bases_loaded=0;


  if (is_list==true)
    {
      //TODO either put in error-handling, or modify function to always succeed - easy to do - do it on next iteration.
      load_se_filelist_into_graph_colour(path->buff,
					 qual_thresh, 
					 homopolymer_cutoff, 
					 false,
					 ascii_fq_offset,
					 into_colour, 
					 db_graph, 
					 0, // 0 => filelist not colourlist
					 &num_files_loaded, 
					 &total_bad_reads, 
					 &total_dup_reads,
					 &total_bases_read, 
					 &total_bases_loaded,
					 readlen_distrib, 
					 readlen_distrib_len, 
					 subsample_function,
					 only_load_pre_existing_kmers);

    }
  else
    {
      num_files_loaded=1;
      load_se_seq_data_into_graph_colour(path->buff,
					 qual_thresh, 
					 homopolymer_cutoff, 
					 false,
					 ascii_fq_offset, 
					 into_colour, 
					 db_graph,
					 &total_bad_reads, 
					 &total_dup_reads,
					 &total_bases_read, 
					 &total_bases_loaded,
					 readlen_distrib, 
					 readlen_distrib_len, 
					 subsample_function,
					 only_load_pre_existing_kmers);
					 
    }
  db_graph_get_covg_distribution_array(db_graph, 0, &db_node_condition_always_true,
				       kmer_covg_array, len_kmer_covg_array);

  return total_bases_loaded;
}

void db_graph_get_covg_distribution_array(dBGraph* db_graph, int colour,
					  boolean (*condition)(dBNode* elem),
					  uint64_t * kmer_covg_array, 
					  uint32_t len_kmer_covg_array)
{
  if (kmer_covg_array==NULL)
    {
      return;
    }
  uint32_t i;


  for(i = 0; i < len_kmer_covg_array; i++)
    kmer_covg_array[i] = (Covg) 0;

  void bin_covg_and_add_to_array(HashTable* htable, Element *e,
                                 uint64_t* arr, uint32_t len, int colour)
  {
    if(condition(e)==true)
    {
      uint64_t bin = MIN(e->coverage[colour], len-1);
      arr[bin]++;
    }
  }
  
  db_graph_traverse_with_array_of_uint64(&bin_covg_and_add_to_array, db_graph,
					 kmer_covg_array, len_kmer_covg_array, colour);


}


//if we return -1 we have failed due to OOM
int choose_cleaning_threshold(uint64_t* kmer_covg_array, 
			      uint32_t len_kmer_covg_array,
			      int expected_depth)//passing D_eff as integer - good enough
{
  
  /* debug
  int p;
  for (p=0; p<=30; p++)
    {
      printf("%d\t%" PRIu64 "\n", p, kmer_covg_array[p]);
    }
  */

  double* d1 = calloc(len_kmer_covg_array+1, sizeof(double));
  double* d2 = calloc(len_kmer_covg_array+1, sizeof(double));

  if ( (d1==NULL) || (d2==NULL) )
    {
      return -1;
    }
  
  int i;
  int firstPosD1=0;
  boolean found_first_posd1=false;

  //kmer_covg_array has values for covg=0 (there will be none) and covg=1, etc
  //I want d1 to start with index1 (ie leave index=0 at value 0.
  for (i=2; i<len_kmer_covg_array; i++)
    {
      if (kmer_covg_array[i]>0)
	{
	  d1[i-1]=(double) kmer_covg_array[i-1]/ (double)kmer_covg_array[i];
	}
      else
	{
	  d1[i-1]=(double) COVG_MAX;
	}
      if ( (found_first_posd1==false) && (d1[i-1]<1) )
	{
	  found_first_posd1=true;
	  firstPosD1=i-1;
	}
      //printf("d1 %d\t%f\n", i-1, d1[i-1]);
    }

  int firstPosD2=0;
  boolean found_first_posd2=false;

  for (i=2; i<len_kmer_covg_array; i++)
    {
      if (d1[i]>0)
	{
	  d2[i-1]=d1[i-1]/d1[i];
	}
      else
	{
	  d2[i-1]=(double) COVG_MAX;
	}
      if ( (found_first_posd2==false) && (d2[i-1]<=1) )
	{
	  found_first_posd2=true;
	  firstPosD2=i-1;
	}
      //printf("d2 %d\t%f\n", i-1, d2[i-1]);
    }

  if ( (firstPosD1!=0) && (firstPosD1< 0.75*expected_depth) )
    {
      //printf("Use D1 threshold %d\n", firstPosD1);
      return firstPosD1;
    }
  else if (firstPosD2!=0)
    {
      //printf("Use D2 threshold %d\n", firstPosD2);
      return firstPosD2;
    }
  else
    {
      int ed = (expected_depth+0.5)/2; 
      if (ed>1)
	{
	  //printf("Use threshold ed %d\n", ed);
	  return ed;
	} 
      else
	{
	  //printf("Fall back on threshold 1\n");
	  return 1;
	}
    }
}

boolean clean_graph(dBGraph* db_graph, 
		    uint64_t* kmer_covg_array, uint32_t len, 
		    int d_eff, int max_expected_sup)
{
  int thresh = choose_cleaning_threshold(kmer_covg_array, len, d_eff); 
  if (thresh==-1)
    {
      return false;
    }
  //printf("Chose cleaning threshold %d\n", thresh);
  db_graph_clip_tips_in_union_of_all_colours(db_graph);

  //printf("Remove low coverage supernodes covg (<= %d) \n", thresh);
  db_graph_remove_errors_considering_covg_and_topology(thresh,
						       db_graph, 
						       &element_get_covg_union_of_all_covgs, 
						       &element_get_colour_union_of_all_colours,
						       &apply_reset_to_specific_edge_in_union_of_all_colours, 
						       &apply_reset_to_all_edges_in_union_of_all_colours,
						       max_expected_sup);

  return true;
}






                                                                          
//This function does not  check that it there is such an edge in the specified person/colour - but it does check if the target node is in the specific person.
//if you want to be sure the dge exists in that colour, then check it before calling this function
dBNode * db_graph_get_next_node_for_specific_person_or_pop(dBNode * current_node, Orientation current_orientation, 
							   Orientation * next_orientation,
							   Nucleotide edge, Nucleotide * reverse_edge,dBGraph * db_graph, int index)
{

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

  if (next_node != NULL)
    {
      *next_orientation = db_node_get_orientation(&local_copy_of_kmer,next_node,db_graph->kmer_size);
    }
  else
    {
      //no else
    }


  //need to check the node is in this person's graph
  if (! (db_node_is_this_node_in_this_person_or_populations_graph(next_node, index)))
    {
      return NULL;
    }
 
  return next_node;
  

}



//This function does not  check that it there is such an edge in the specified person/colour - but it does check if the target node is in the specific person.
//if you want to be sure the dge exists in that colour, then check it before calling this function
//The last argument allows you to apply the operation to some subgraph - eg you might take the unuiion of colours 2 and 3, or of all colours.
//If for example you wanted to get the next node in the graph irrespective of colour, get_colour would return the union (bitwise AND) of all edges in a node.
dBNode * db_graph_get_next_node_in_subgraph_defined_by_func_of_colours(dBNode * current_node, Orientation current_orientation, 
								       Orientation * next_orientation,
								       Nucleotide edge, Nucleotide * reverse_edge,dBGraph * db_graph, 
								       Edges (*get_colour)(const dBNode*)
								       )
{

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

  if (next_node != NULL)
    {
      *next_orientation = db_node_get_orientation(&local_copy_of_kmer,next_node,db_graph->kmer_size);
    }
  else
    {
      //no else
    }


  //need to check the node is in the specified subgraph graph
  if (! (db_node_is_this_node_in_subgraph_defined_by_func_of_colours(next_node, get_colour)) )
    {
      return NULL;
    }
 
  return next_node;
  

}


//MARIO Mario - noode seems to call or use this?
boolean db_graph_db_node_has_precisely_n_edges_with_status_for_specific_person_or_pop(
  dBNode * node,Orientation orientation,NodeStatus status,int n,
  dBNode * * next_node, Orientation * next_orientation, Nucleotide * next_base,
  dBGraph * db_graph, int index)
{

  if ( (n>4)||(n<0))
    {
      warn("calling db_graph_db_node_has_precisely_n_edges_with_status with n=%d", n);
      return false;
    }

  int count = 0;
  boolean ret = false;

  void check_edge(Nucleotide base){

    dBNode * tmp_next_node;
    Orientation tmp_next_orientation;
    Nucleotide rev_base;

   

    if (db_node_edge_exist(node,base,orientation, index)){

      tmp_next_node = db_graph_get_next_node_for_specific_person_or_pop(node,orientation,&tmp_next_orientation,
									base,&rev_base,db_graph, index);

      if (db_node_check_status(tmp_next_node,status)){
	next_base[count] = base;
	next_node[count] = tmp_next_node;
	
	next_orientation[count] = tmp_next_orientation;
	count++;
      }
      
    }
    
  }
  
  
  nucleotide_iterator(&check_edge);

  if (count == n){
    ret = true;
  }

  return ret;

}


int db_graph_db_node_clip_tip_with_orientation_for_specific_person_or_pop(
  dBNode * node, Orientation orientation, int limit,
  void (*node_action)(dBNode * node),dBGraph * db_graph, int index)
{ 

  Nucleotide nucleotide, reverse_nucleotide;
  int length = 0;
  int i;

  dBNode** nodes = (dBNode**) malloc(sizeof(dBNode*)*limit);
  if (nodes==NULL)
    {
      die("Unable to malloc array of nodes for tip clipping\n");
    }
  Orientation next_orientation;
  dBNode * next_node;
  char seq[db_graph->kmer_size+1];
  
  //starting in a blunt end also prevents full loops 
  if (db_node_is_blunt_end(node, opposite_orientation(orientation), index)){
    
    boolean join_main_trunk = false;

    while(db_node_has_precisely_one_edge(node, orientation, &nucleotide, index)) {
  
       nodes[length] = node;

       next_node = db_graph_get_next_node_for_specific_person_or_pop(node,orientation,&next_orientation,nucleotide,&reverse_nucleotide,db_graph, index);
       
       if(next_node == NULL){
	 die("dB_graph_db_node_clip_tip_with_orientation_for_specific_person_or_pop: "
       "didnt find node in hash table: %s",
       binary_kmer_to_seq(element_get_kmer(node), db_graph->kmer_size, seq));
       }
       
       length ++;

       
       if (length>limit){
	 break;
       }

       //want to stop when we join another trunk
       if (!db_node_has_precisely_one_edge(next_node, opposite_orientation(next_orientation), &nucleotide, index)){
	 join_main_trunk = true;
	 break;
       }

       //keep track of node
       
       node = next_node;
       orientation = next_orientation;        
     }
  
     if (! join_main_trunk){
       length = 0;
     }
     else{//clear edges and mark nodes as pruned
       for(i=0;i<length;i++){

	 if (DEBUG){
	   printf("CLIPPING node: %s\n",binary_kmer_to_seq(element_get_kmer(nodes[i]),db_graph->kmer_size,seq));
	 }
	 
	 node_action(nodes[i]);
	 //perhaps we want to move this inside the node action?
	 db_node_reset_edges(nodes[i], index);
       }

       if (DEBUG){
	printf("RESET %c BACK\n",binary_nucleotide_to_char(reverse_nucleotide));
      }
       db_node_reset_edge(next_node,opposite_orientation(next_orientation),reverse_nucleotide, index);
     }

  }

  free(nodes);
  return length;
}




// clip a tip in the graph (the tip starts in node)
// limit is max length for tip
// node_action is applied to all the elements in the tip
// returns the length of the tip (0 means no length)
int db_graph_db_node_clip_tip_for_specific_person_or_pop(dBNode * node, int limit,
							 void (*node_action)(dBNode * node),
							 dBGraph * db_graph, int index)
{

  int length_tip = 0;

  
  length_tip = db_graph_db_node_clip_tip_with_orientation_for_specific_person_or_pop(node,forward,limit,node_action,db_graph, index);
  
  if (length_tip==0){
    length_tip = db_graph_db_node_clip_tip_with_orientation_for_specific_person_or_pop(node,reverse,limit,node_action,db_graph, index);
    
  }
  
  return length_tip;
}




int db_graph_db_node_clip_tip_with_orientation_in_subgraph_defined_by_func_of_colours(dBNode * node, Orientation orientation, int limit,
										       void (*node_action)(dBNode * node),dBGraph * db_graph, 
										       Edges (*get_colour)(const dBNode*),
										      void (*apply_reset_to_specific_edge_in_colour)(dBNode*, Orientation, Nucleotide),
										       void (*apply_reset_to_colour)(dBNode*)
										       )
{ 

  Nucleotide nucleotide, reverse_nucleotide;
  int length = 0;
  int i;
  dBNode** nodes=(dBNode**) malloc(sizeof(dBNode*)*limit);

  Orientation next_orientation;
  dBNode * next_node;
  char seq[db_graph->kmer_size+1];
  
  //starting in a blunt end also prevents full loops 
  if (db_node_is_blunt_end_in_subgraph_given_by_func_of_colours(node, opposite_orientation(orientation), get_colour))
    {
      boolean join_main_trunk = false;

      while(db_node_has_precisely_one_edge_in_subgraph_defined_by_func_of_colours(node, orientation, &nucleotide, get_colour)) 
	{
	  nodes[length] = node;
	  
	  next_node = db_graph_get_next_node_in_subgraph_defined_by_func_of_colours(node,orientation,&next_orientation,nucleotide,&reverse_nucleotide,db_graph, get_colour);
	  
	  if(next_node == NULL){
	    die("dB_graph_db_node_clip_tip_with_orientation_for_specific_person_or_pop: "
          "didnt find node in hash table: %s",
          binary_kmer_to_seq(element_get_kmer(node), db_graph->kmer_size, seq));
	  }
	  
	  length ++;
	  
	  
	  if (length>limit){
	    break;
	  }
	  
	  //want to stop when we join another trunk
	  if (!db_node_has_precisely_one_edge_in_subgraph_defined_by_func_of_colours(next_node, opposite_orientation(next_orientation), &nucleotide, get_colour)){
	    join_main_trunk = true;
	    break;
	  }

	  //keep track of node
	  
	  node = next_node;
	  orientation = next_orientation;        
	}
      
      if (! join_main_trunk)
	{
	  length = 0;
	}
      else
	{//clear edges and mark nodes as pruned
	  for(i=0;i<length;i++)
	    {
	      if (DEBUG){
		printf("CLIPPING node: %s\n",binary_kmer_to_seq(element_get_kmer(nodes[i]),db_graph->kmer_size,seq));
	      }
	      
	      node_action(nodes[i]);
	      //perhaps we want to move this inside the node action?
	      apply_reset_to_colour(nodes[i]);
	    }
	  
	  if (DEBUG){
	    printf("RESET %c BACK\n",binary_nucleotide_to_char(reverse_nucleotide));
	  }
	  apply_reset_to_specific_edge_in_colour(next_node,opposite_orientation(next_orientation),reverse_nucleotide);
     }
      
    }
  free(nodes);
  return length;
}




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
								     )
{

  int length_tip = 0;

  
  length_tip = db_graph_db_node_clip_tip_with_orientation_in_subgraph_defined_by_func_of_colours(node,forward,limit,node_action,db_graph, get_colour, apply_reset_to_specific_edge_in_colour, apply_reset_to_colour );
  
  if (length_tip==0){
    length_tip = db_graph_db_node_clip_tip_with_orientation_in_subgraph_defined_by_func_of_colours(node,reverse,limit,node_action,db_graph, get_colour, apply_reset_to_specific_edge_in_colour, apply_reset_to_colour);
    
  }
  
  return length_tip;
}






// 1. the argument sum_of_covgs_in_desired_colours allows you to choose which colour (or union of colours) you want to apply this to
// eg you might want  to "remove" (as defined by the action) any node that has coverage <= your threshold in the UNION of all colours, or in colour red or whatever
// SO - this func returns the sum of coverages in the colours you care about
// 2. the argument get_edge_of_interest is a function that gets the "edge" you are interested in - may be a single edge/colour from the graph, or might be a union of some edges
// 3 Pass apply_reset_to_specified_edges which applies reset_one_edge to whichever set of edges you care about, 
// 4 Pass apply_reset_to_specified_edges_2 which applies db_node_reset_edges to whichever set of edges you care about, 
boolean db_graph_db_node_prune_low_coverage(dBNode * node, Covg coverage,
					    void (*node_action)(dBNode * node),
					    dBGraph * db_graph, 
					    Covg (*sum_of_covgs_in_desired_colours)(const Element *), 
					    Edges (*get_edge_of_interest)(const Element*),
					    void (*apply_reset_to_specified_edges)(dBNode*, Orientation, Nucleotide),
					    void (*apply_reset_to_specified_edges_2)(dBNode*) )
{

  boolean ret = false;

  if (sum_of_covgs_in_desired_colours(node)<=coverage){
    ret = true;
  
    void nucleotide_action(Nucleotide nucleotide){
      Orientation next_orientation;
      Nucleotide reverse_nucleotide;
      dBNode * next_node;


      //remove evidence in adjacent nodes where they point to this one
      if (db_node_edge_exist_within_specified_function_of_coloured_edges(node,nucleotide,forward, get_edge_of_interest)){
	next_node = db_graph_get_next_node(node,forward,&next_orientation,nucleotide,&reverse_nucleotide,db_graph);
	db_node_reset_specified_edges(next_node,opposite_orientation(next_orientation),reverse_nucleotide, apply_reset_to_specified_edges);
	  
      }
	
      if (db_node_edge_exist_within_specified_function_of_coloured_edges(node,nucleotide,reverse, get_edge_of_interest)){
	next_node = db_graph_get_next_node(node,reverse,&next_orientation,nucleotide,&reverse_nucleotide,db_graph);
	db_node_reset_specified_edges(next_node,opposite_orientation(next_orientation),reverse_nucleotide, apply_reset_to_specified_edges);
	
      }
    }

    nucleotide_iterator(&nucleotide_action);

    node_action(node);
    //remove all edges from this node, in colours we care about
    apply_reset_to_specified_edges_2(node);   

  }

  return ret;
}



boolean db_graph_db_node_prune_without_condition(dBNode * node, 
						 void (*node_action)(dBNode * node),
						 dBGraph * db_graph, 
						 Edges (*get_edge_of_interest)(const Element*),
						 void (*apply_reset_to_specified_edges)(dBNode*, Orientation, Nucleotide),
						 void (*apply_reset_to_specified_edges_2)(dBNode*) )
{

  boolean ret = true;
  
  void nucleotide_action(Nucleotide nucleotide){
    Orientation next_orientation;
    Nucleotide reverse_nucleotide;
    dBNode * next_node;
    

    //remove evidence in adjacent nodes where they point to this one
    if (db_node_edge_exist_within_specified_function_of_coloured_edges(node,nucleotide,forward, get_edge_of_interest)){
      next_node = db_graph_get_next_node(node,forward,&next_orientation,nucleotide,&reverse_nucleotide,db_graph);
      db_node_reset_specified_edges(next_node,opposite_orientation(next_orientation),reverse_nucleotide, apply_reset_to_specified_edges);
      
    }
    
    if (db_node_edge_exist_within_specified_function_of_coloured_edges(node,nucleotide,reverse, get_edge_of_interest)){
      next_node = db_graph_get_next_node(node,reverse,&next_orientation,nucleotide,&reverse_nucleotide,db_graph);
      db_node_reset_specified_edges(next_node,opposite_orientation(next_orientation),reverse_nucleotide, apply_reset_to_specified_edges);
      
    }
  }
  
  nucleotide_iterator(&nucleotide_action);
  
  node_action(node);
  //remove all edges from this node, in colours we care about
  apply_reset_to_specified_edges_2(node);   
  
  return ret;
}


boolean db_graph_db_node_prune_low_coverage_ignoring_colours(dBNode * node, Covg coverage,
							     void (*node_action)(dBNode * node),
							     dBGraph * db_graph)
{

  Edges get_edge_of_interest(const Element* node)
  {
    return get_union_of_edges(*node);
  }
  
  Edges set_to_zero(Edges edge)
  {
    return 0;
  }

  void apply_reset_to_specified_edges(dBNode* node , Orientation or , Nucleotide nuc)
  {
      int i;
      for (i=0; i< NUMBER_OF_COLOURS; i++)
	{
	  reset_one_edge(node, or, nuc, i);
	}
  }
  
  void apply_reset_to_specified_edges_2(dBNode* node)
  {
      int i;
      for (i=0; i< NUMBER_OF_COLOURS; i++)
	{
	  db_node_reset_edges(node, i);
	}
  }

  return db_graph_db_node_prune_low_coverage(node, coverage, node_action, db_graph,
					     &element_get_covg_union_of_all_covgs,
					     &get_edge_of_interest,
					     &apply_reset_to_specified_edges,
					     &apply_reset_to_specified_edges_2);
}


// computes a perfect path starting from a node and an edge
// ie the starting node can have  multiple exits
// perfect path -- no conflict no cycle -- returns length
// path: node_0 edge_0 node_1 edge_1 ... node_n-1 edge_n-1 node_n
// path_nodes is a n+1 array from 0..n with all the nodes in the path
// path_orientations is n+1 array from 0..n with the orientations of the node in the path
// path labels is n array from 0..n-1 with all the labels for the edges (see above)
// node_action only applied to internal nodes (not node_0 and node_n)
// seq is a string with all the labels concatenated (NB: it doesn't contain the kmer in the first node)
// avg_coverage, min_coverge, max_coverge -> refers to the internal nodes only


int db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop(
  dBNode *node, Orientation orientation, int limit, 
  Nucleotide fst_nucleotide,
  void (*node_action)(dBNode *node),
  dBNode **path_nodes, Orientation *path_orientations, Nucleotide *path_labels,
  char *seq, double *avg_coverage, Covg *min_coverage, Covg *max_coverage,
  boolean *is_cycle, dBGraph *db_graph, int index)
{

  Orientation  current_orientation,next_orientation;
  dBNode *current_node = NULL;
  dBNode *next_node = NULL;
  Nucleotide nucleotide,rev_nucleotide,nucleotide2;
  int length =0;
  char tmp_seq[db_graph->kmer_size+1];
  tmp_seq[db_graph->kmer_size] = '\0';
  Covg sum_coverage = 0;
  Covg coverage  = 0;

  //sanity checks
  if (node == NULL)
    {
      die("db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop: "
          "can't pass a null node");
    }
  else if (! (db_node_is_this_node_in_this_person_or_populations_graph(node, index)))
    {
      //printf("\nThis node is not in the graph of this person - in db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop\n");
      return false;
    }

  current_node        = node;
  current_orientation = orientation;

  *is_cycle = false;
  
  path_nodes[0]         = node;
  path_orientations[0]  = orientation;  
  *max_coverage         = 0;
  *min_coverage         = COVG_MAX;
 
  if (DEBUG){
    printf("\n ZAM Node %i in path: %s\n", length, binary_kmer_to_seq(element_get_kmer(current_node),db_graph->kmer_size,tmp_seq));
  }
    
  
  //first edge defined
  nucleotide = fst_nucleotide;

  do{ 
    if (length>0){
      node_action(current_node);
      sum_coverage += coverage;
      *max_coverage = MAX(coverage, *max_coverage);
      *min_coverage = MIN(coverage, *min_coverage);
    }

    //this will return NULL if the next node is not in the person's graph. It does NOT check if the edge is in the graph...
    next_node =  db_graph_get_next_node_for_specific_person_or_pop(current_node,current_orientation,&next_orientation,nucleotide,&rev_nucleotide,db_graph, index);

      


    //sanity check
    if(next_node == NULL)
    {
    	die("db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop: "
          "didnt find node in hash table: %s %c %s",
      binary_kmer_to_seq(element_get_kmer(current_node),db_graph->kmer_size,tmp_seq),
      binary_nucleotide_to_char(nucleotide), current_orientation == forward ? "forward" : "reverse");
    }


    path_labels[length]        = nucleotide;
    seq[length]                = binary_nucleotide_to_char(nucleotide);
    coverage                   = db_node_get_coverage(next_node, index);
    
    length++;

    path_nodes[length]         = next_node;
    path_orientations[length]  = next_orientation;
    
    if (DEBUG)
    {
      binary_kmer_to_seq(element_get_kmer(next_node), db_graph->kmer_size, tmp_seq);
      printf("\n Length of path so far is  %i : %s\n", length, seq);
    }
    
    current_node        = next_node;
    current_orientation = next_orientation;
    
  } while (length<(limit-1) && 
	   !((next_node == node) && (next_orientation == orientation)) && //loop
	   db_node_has_precisely_one_edge(next_node, opposite_orientation(next_orientation), &nucleotide2, index) && //multiple entries
	   db_node_has_precisely_one_edge(current_node, current_orientation, &nucleotide, index)); //has one next edge only
  
  
  if ((next_node == node) && (next_orientation == orientation)){
    *is_cycle = true;
  }

  //debug
  /*
  if ((next_node == node) && (next_orientation == orientation))
    {
      printf("stopped becaise of loop\n");
    }
  else if (!db_node_has_precisely_one_edge(next_node,opposite_orientation(next_orientation),&nucleotide2, index))
    {
      printf("stopped because next node has >1 edge in - multple entries\n");
    }
  else if (!db_node_has_precisely_one_edge(current_node, current_orientation, &nucleotide, index))
    {
      printf("stopped because current node has >1 edge\n");
    }
  else if (length>=limit)
    {
      printf("stopped becase limit exceeded: length %d and limit %d\n", length, limit);
    }
  else
    {
      printf("WARNING IMPOSSIBLE");
    }
  */
  

  if (length>=limit)
    {
      die("Stopped becase supernode length limit exceeded: length %d and limit %d", length, limit);
    }
  
   seq[length] = '\0';
  *avg_coverage = (length-1<=0) ? 0 : (double) sum_coverage/(double) (length-1);

  if (*min_coverage == COVG_MAX)
    {
      *min_coverage = 0;
    }

  return length;
  
}



// computes a perfect path starting from a node and an edge
// ie the starting node can have  multiple exits
// perfect path -- no conflict no cycle -- returns length
// path: node_0 edge_0 node_1 edge_1 ... node_n-1 edge_n-1 node_n
// path_nodes is a n+1 array from 0..n with all the nodes in the path
// path_orientations is n+1 array from 0..n with the orientations of the node in the path
// path labels is n array from 0..n-1 with all the labels for the edges (see above)
// node_action only applied to internal nodes (not node_0 and node_n)
// seq is a string with all the labels concatenated (NB: it doesn't contain the kmer in the first node)
// avg_coverage, min_coverge, max_coverge -> refers to the internal nodes only

//get_colour returns an edge which is a function of the edges in a node. eg union of all edges, or of colours 1 and 2
// get covg returns coverge desired - most ikely is sum of covgs in each of colours which are summed/unioned in get_colour
int db_graph_get_perfect_path_with_first_edge_in_subgraph_defined_by_func_of_colours(
										     dBNode *node, Orientation orientation, int limit, Nucleotide fst_nucleotide,
										     void (*node_action)(dBNode *node), dBNode **path_nodes,
										     Orientation *path_orientations, Nucleotide *path_labels, char *seq,
										     double *avg_coverage, Covg *min_coverage, Covg *max_coverage,
										     boolean *is_cycle, dBGraph *db_graph, 
										     Edges (*get_colour)(const dBNode*),
										     Covg (*get_covg)(const dBNode*))
{

  Orientation  current_orientation,next_orientation;
  dBNode * current_node = NULL;
  dBNode * next_node = NULL;
  Nucleotide nucleotide,rev_nucleotide,nucleotide2;
  int length =0;
  char tmp_seq[db_graph->kmer_size+1];
  tmp_seq[db_graph->kmer_size]='\0';
  Covg sum_coverage = 0;
  Covg coverage  = 0;

  //sanity checks
  if (node == NULL)
    {
      die("db_graph_get_perfect_path_with_first_edge_in_subgraph_defined_by_func_of_colours:"
          "can't pass a null node");
    }
  else if (! (db_node_is_this_node_in_subgraph_defined_by_func_of_colours(node, get_colour)))
    {
      warn("This node is not in the graph of this person - in db_graph_get_perfect_path_\n");
      return false;
    }

  current_node        = node;
  current_orientation = orientation;

  *is_cycle = false;
  
  path_nodes[0]         = node;
  path_orientations[0]  = orientation;  
  *max_coverage         = 0;
  *min_coverage         = COVG_MAX;
 
  if (DEBUG){
    printf("\n ZAM Node %i in path: %s\n", length, binary_kmer_to_seq(element_get_kmer(current_node),db_graph->kmer_size,tmp_seq));
  }
    
  
  //first edge defined
  nucleotide = fst_nucleotide;

  do{ 
    if (length>0){
      node_action(current_node);
      sum_coverage += coverage;
      *max_coverage = MAX(coverage, *max_coverage);
      *min_coverage = MIN(coverage, *min_coverage);
    }

    //this will return NULL if the next node is not in the subgraph specified by get_colour. 
    next_node =  db_graph_get_next_node_in_subgraph_defined_by_func_of_colours(current_node,current_orientation,&next_orientation,nucleotide,&rev_nucleotide,db_graph, get_colour);

      


    //sanity check
    if(next_node == NULL)
      {
	die("db_graph_get_perfect_path_with_first_edge_in_subgraph_defined_by_func_of_colours: "
      "didnt find node in hash table: %s %c %s",
      binary_kmer_to_seq(element_get_kmer(current_node),db_graph->kmer_size,tmp_seq),
      binary_nucleotide_to_char(nucleotide),
      current_orientation == forward ? "forward" : "reverse");
      }


    path_labels[length]        = nucleotide;
    seq[length]                = binary_nucleotide_to_char(nucleotide);
    coverage                   = db_node_get_coverage_in_subgraph_defined_by_func_of_colours(next_node, get_covg);
    
    length++;

    path_nodes[length]         = next_node;
    path_orientations[length]  = next_orientation;
    
    if (DEBUG){
      printf("\n Length of path so far is  %i : %s\n", length, binary_kmer_to_seq(element_get_kmer(next_node),db_graph->kmer_size,tmp_seq));
    }
    
    current_node        = next_node;
    current_orientation = next_orientation;
    
  } while (length<(limit-1) && 
	   !((next_node == node) && (next_orientation == orientation)) && //loop
	   db_node_has_precisely_one_edge_in_subgraph_defined_by_func_of_colours(next_node, opposite_orientation(next_orientation), &nucleotide2, get_colour) && //multiple entries
	   db_node_has_precisely_one_edge_in_subgraph_defined_by_func_of_colours(current_node, current_orientation, &nucleotide, get_colour)); //has one next edge only
  
  
  if ((next_node == node) && (next_orientation == orientation)){
    *is_cycle = true;
  }

  
  
   seq[length] = '\0';
  *avg_coverage = (length-1<=0) ? 0 : (double) sum_coverage/(double) (length-1);

  if (*min_coverage == COVG_MAX)
    {
      *min_coverage = 0;
    }

  if (length>=limit)
    {
      die("This genome has a supernode longer than the maximum %d you have specified with --max_var_len. Re-run with higher --max_var_len\n", limit);
    }
  return length;
  
}






// perfect path -- no conflict no cycle -- returns length
// path: node_0 edge_0 node_1 edge_1 ... node_n-1 edge_n-1 node_n
// path_nodes is a n+1 array from 0..n with all the nodes in the path
// path_orientations is n+1 array from 0..n with the orientations of the node in the path
// path labels is n array from 0..n-1 with all the labels for the edges (see above)
// node_action only applied to internal nodes (not node_0 and node_n)
// seq is a string with all the labels concatenated (NB: it doesn't contain the kmer in the first node)
// avg_coverage, min_coverge, max_coverge -> refers to the internal nodes only

int db_graph_get_perfect_path_for_specific_person_or_pop(
  dBNode *node, Orientation orientation, int limit,
  void (*node_action)(dBNode *node), dBNode **path_nodes,
  Orientation *path_orientations, Nucleotide *path_labels, char *seq,
  double *avg_coverage, Covg *min_coverage, Covg *max_coverage,
  boolean *is_cycle, dBGraph *db_graph, int index)
{

  int length =0;
  Nucleotide nucleotide;

  //sanity check
  if (node == NULL){
    die("db_graph_get_perfect_path_for_specific_person_or_pop: can't pass a null node");
  }

  path_nodes[0]         = node;
  path_orientations[0]  = orientation;  
 

  if (db_node_has_precisely_one_edge(node, orientation, &nucleotide, index))
    {
    
      length= db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop(node,orientation, limit, nucleotide,
										   node_action,
										   path_nodes,path_orientations,path_labels,
										   seq,avg_coverage,min_coverage,max_coverage,
										   is_cycle,db_graph, index);
    }
  else{
    *max_coverage         = 0;
    *min_coverage         = 0;
    seq[0] = '\0';
  }

  return length;

}


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
  Edges (*get_colour)(const dBNode*), Covg (*get_covg)(const dBNode*))
{

  int length =0;
  Nucleotide nucleotide;

  //sanity check
  if (node == NULL){
    die("db_graph_get_perfect_path_in_subgraph_defined_by_func_of_colours: can't pass a null node");
  }

  path_nodes[0]         = node;
  path_orientations[0]  = orientation;  
 

  if (db_node_has_precisely_one_edge_in_subgraph_defined_by_func_of_colours(node, orientation, &nucleotide, get_colour))
    {
    
      length= db_graph_get_perfect_path_with_first_edge_in_subgraph_defined_by_func_of_colours(node,orientation, limit, nucleotide,
											       node_action,
											       path_nodes,path_orientations,path_labels,
											       seq,avg_coverage,min_coverage,max_coverage,
											       is_cycle,db_graph, 
											       get_colour, get_covg);
    }
  else{
    *max_coverage         = 0;
    *min_coverage         = 0;
    seq[0] = '\0';
  }

  return length;

}







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
  dBGraph *db_graph, int index)
{

  //use the allocated space as a temporary space
  dBNode **nodes_reverse = path_nodes;
  Orientation *orientations_reverse = path_orientations;
  Nucleotide *labels_reverse = path_labels;
     
  boolean is_cycler;
  int length_reverse;
  int length = 0;
  
  Covg minr,maxr;
  double avg_coverager;


  //compute the reverse path until the end of the supernode
  //return is_cycle_reverse == true if the path closes a loop    
  
  length_reverse
    = db_graph_get_perfect_path_for_specific_person_or_pop(
        node, reverse, limit, &db_node_action_do_nothing, 
        nodes_reverse, orientations_reverse, labels_reverse, 
        supernode_str, &avg_coverager, &minr, &maxr, 
        &is_cycler, db_graph, index);
  
  if (length_reverse>0){
    //let's re do the last step, we need to do that because the last node could have had multiple entries
    
    Nucleotide label;
    Orientation next_orientation;
    
    dBNode * lst_node = db_graph_get_next_node(nodes_reverse[length_reverse-1],orientations_reverse[length_reverse-1],
					       &next_orientation, labels_reverse[length_reverse-1],&label,db_graph);
    
    //sanity check
    if (lst_node != nodes_reverse[length_reverse]){
      die("db_graph_supernode broken!");
    }
    
    
    length = db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop(nodes_reverse[length_reverse],
										  opposite_orientation(orientations_reverse[length_reverse]),
										  limit,label,
										  node_action,
										  path_nodes,path_orientations,path_labels,
										  supernode_str,avg_coverage,min,max,
										  is_cycle,db_graph, index);
    
  }
  else{
    length = db_graph_get_perfect_path_for_specific_person_or_pop(node,forward,
								   limit,
								   node_action,
								   path_nodes,path_orientations,path_labels,
								   supernode_str,avg_coverage,min,max,
								   is_cycle,db_graph, index);
  }
  
  
  

  //apply action to the fst and last node
  node_action(path_nodes[0]);
  node_action(path_nodes[length]);
  
  return length;
}


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
  Covg (*get_covg)(const dBNode*))
{

  
  //use the allocated space as a temporary space
  dBNode * * nodes_reverse = path_nodes;
  Orientation * orientations_reverse = path_orientations;
  Nucleotide * labels_reverse = path_labels;
     
  boolean is_cycler;
  int length_reverse;
  int length = 0;
  
  Covg minr, maxr;
  double avg_coverager;


  //compute the reverse path until the end of the supernode
  //return is_cycle_reverse == true if the path closes a loop    
  
  length_reverse = db_graph_get_perfect_path_in_subgraph_defined_by_func_of_colours(node,reverse,limit,&db_node_action_do_nothing,
										    nodes_reverse,orientations_reverse,labels_reverse,
										    supernode_str,&avg_coverager,&minr,&maxr,
										    &is_cycler,db_graph, get_colour, get_covg);
  
  if (length_reverse>0){
    //let's re do the last step, we need to do that because the last node could have had multiple entries
    
    Nucleotide label;
    Orientation next_orientation;
    
    dBNode * lst_node = db_graph_get_next_node_in_subgraph_defined_by_func_of_colours(nodes_reverse[length_reverse-1],orientations_reverse[length_reverse-1],
										      &next_orientation, labels_reverse[length_reverse-1],&label,db_graph, get_colour);
    
    //sanity check
    if (lst_node != nodes_reverse[length_reverse]){
      die("db_graph_supernode broken!");
    }
    
    
    length = db_graph_get_perfect_path_with_first_edge_in_subgraph_defined_by_func_of_colours(nodes_reverse[length_reverse],
											      opposite_orientation(orientations_reverse[length_reverse]),
											      limit,label,
											      node_action,
											      path_nodes,path_orientations,path_labels,
											      supernode_str,avg_coverage,min,max,
											      is_cycle,db_graph, get_colour, get_covg);
    
  }
  else{
    length = db_graph_get_perfect_path_in_subgraph_defined_by_func_of_colours(node,forward,
									      limit,
									      node_action,
									      path_nodes,path_orientations,path_labels,
									      supernode_str,avg_coverage,min,max,
									      is_cycle,db_graph, get_colour, get_covg);
  }
  
  
  

  //apply action to the fst and last node
  node_action(path_nodes[0]);
  node_action(path_nodes[length]);
  
  return length;
}




// identical code to db_graph_supernode_for_specific_person_or_pop but this returns the index of the query node (first argument) in th supernode array
// will make a significant performance gain compared with getting the path_nodes array back and then searching it

int db_graph_supernode_returning_query_node_posn_for_specific_person_or_pop(
  dBNode *node,int limit,void (*node_action)(dBNode *node), 
  dBNode **path_nodes, Orientation *path_orientations, Nucleotide *path_labels,
  char *supernode_str, double *avg_coverage, Covg *min, Covg *max,
  boolean *is_cycle, int *query_node_posn, dBGraph *db_graph, int index)
{

  //use the allocated space as a temporary space
  dBNode * * nodes_reverse = path_nodes;
  Orientation * orientations_reverse = path_orientations;
  Nucleotide * labels_reverse = path_labels;
     
  boolean is_cycler;
  int length_reverse;
  int length = 0;
  
  Covg minr,maxr;
  double avg_coverager;


  //compute the reverse path until the end of the supernode
  //return is_cycle_reverse == true if the path closes a loop    
  
  length_reverse = db_graph_get_perfect_path_for_specific_person_or_pop(node,reverse,limit,&db_node_action_do_nothing,
									nodes_reverse,orientations_reverse,labels_reverse,
									supernode_str,&avg_coverager,&minr,&maxr,
									&is_cycler,db_graph, index);

  *query_node_posn = length_reverse;
  
  
  if (length_reverse>0){
    //let's re do the last step, we need to do that because the last node could have had multiple entries
    
    Nucleotide label;
    Orientation next_orientation;
    
    dBNode * lst_node = db_graph_get_next_node(nodes_reverse[length_reverse-1],orientations_reverse[length_reverse-1],
					       &next_orientation, labels_reverse[length_reverse-1],&label,db_graph);
    
    //sanity check
    if (lst_node != nodes_reverse[length_reverse]){
      die("db_graph_supernode broken!");
    }
    
    
    length = db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop(nodes_reverse[length_reverse],
										  opposite_orientation(orientations_reverse[length_reverse]),
										  limit,label,
										  node_action,
										  path_nodes,path_orientations,path_labels,
										  supernode_str,avg_coverage,min,max,
										  is_cycle,db_graph, index);
    
  }
  else{
    length = db_graph_get_perfect_path_for_specific_person_or_pop(node,forward,
								   limit,
								   node_action,
								   path_nodes,path_orientations,path_labels,
								   supernode_str,avg_coverage,min,max,
								   is_cycle,db_graph, index);
  }
  
  
  

  //apply action to the fst and last node
  node_action(path_nodes[0]);
  node_action(path_nodes[length]);
  
  return length;
}



void db_graph_clip_tips_for_specific_person_or_pop(dBGraph * db_graph, int index)
{
  
  void clip_tips(dBNode * node){
    
    //use max length k+1, which is what you would get with a single base error - a bubble of that length
    if (db_node_check_status_none(node)){
      db_graph_db_node_clip_tip_for_specific_person_or_pop(node, 1+db_graph->kmer_size,&db_node_action_set_status_pruned,db_graph, index);
    }
  }

  hash_table_traverse(&clip_tips,db_graph);
  
}





void db_graph_clip_tips_in_subgraph_defined_by_func_of_colours(dBGraph * db_graph,
							       Edges (*get_colour)(const dBNode*),
							       void (*apply_reset_to_specific_edge_in_colour)(dBNode*, Orientation, Nucleotide),
							       void (*apply_reset_to_colour)(dBNode*))
{
  
  void clip_tips(dBNode * node){
    
    //use max length k+1, which is what you would get with a single base error - a bubble of that length
    if (db_node_check_status_none(node)){
      db_graph_db_node_clip_tip_in_subgraph_defined_by_func_of_colours(node, 1+db_graph->kmer_size,&db_node_action_set_status_pruned,db_graph, 
								       get_colour, apply_reset_to_specific_edge_in_colour, apply_reset_to_colour);
    }
  }

  hash_table_traverse(&clip_tips,db_graph);
  
}


void apply_reset_to_specific_edge_in_union_of_all_colours(dBNode* node, Orientation or, Nucleotide nuc)
{
  int j;
  for (j=0; j<NUMBER_OF_COLOURS; j++)
    {
      reset_one_edge(node, or, nuc, j);
    }
}

void apply_reset_to_all_edges_in_union_of_all_colours(dBNode* node )
{
  int j;
  for (j=0; j<NUMBER_OF_COLOURS; j++)
    {
      db_node_reset_edges(node, j);
    }
}


void db_graph_clip_tips_in_union_of_all_colours(dBGraph* db_graph)
{
  db_graph_clip_tips_in_subgraph_defined_by_func_of_colours(db_graph,
							    &element_get_colour_union_of_all_colours,
							    &apply_reset_to_specific_edge_in_union_of_all_colours,
							    &apply_reset_to_all_edges_in_union_of_all_colours);
}







void traverse_hash_collecting_sums(
  void (*f)(Element *, Covg*, Covg*, Covg*, Covg*, Covg*, Covg*, Covg*),
  HashTable * hash_table, 
  Covg* pgf, Covg* cox, Covg* data1, Covg* data2, Covg* data3, Covg* data4, Covg* data5)
{
  long long i;
  for(i=0;i<hash_table->number_buckets * hash_table->bucket_size;i++){
    if (!db_node_check_status(&hash_table->table[i],unassigned)){
      f(&hash_table->table[i], pgf, cox, data1, data2, data3, data4, data5);
    }
  }
}



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
  char* supernode_str, int* supernode_len)
{

  boolean is_supernode_pruned=true;

  if (node==NULL)
    {
      die("Called db_graph_remove_supernode_containing_this_node_if_looks_like_induced_by_error with a NULL node. Programming error\n");
    }

  if ( 
      (db_node_check_status(node, none)==false)//don't touch stuff that is visited or pruned, or whatever
      ||
      (sum_of_covgs_in_desired_colours(node)==0)
       )
    {
      *supernode_len=-1;//caller can check this
      is_supernode_pruned=false;
      return is_supernode_pruned;
    }


  if ( (sum_of_covgs_in_desired_colours(node)>0) && (sum_of_covgs_in_desired_colours(node)<=coverage) )
      {
	//get the supernode, setting nodes to visited
	double avg_cov;
	Covg min_cov;
	Covg max_cov;
	boolean is_cycle;
	
	//length_sup is the number of edges in the supernode
	int length_sup =  db_graph_supernode_in_subgraph_defined_by_func_of_colours(node,max_expected_sup_len,
										    &db_node_action_set_status_visited,
										    path_nodes, path_orientations, path_labels, supernode_str,
										    &avg_cov,&min_cov, &max_cov, &is_cycle,
										    db_graph, 
										    get_edge_of_interest,
										    sum_of_covgs_in_desired_colours);

	*supernode_len=length_sup;

	if (length_sup <=1)
	  {
	    is_supernode_pruned=false;

	    return is_supernode_pruned;//do nothing. This supernode has no interior, is just 1 or 2 nodes, so cannot prune it
	  }
	else// if (length_sup <= 2*db_graph->kmer_size +2)
	  {
	    int i;
	    //to look like an error, must all have actual coverage, caused by an actual errored read, BUT must have low covg, <=threshold
	    boolean interior_nodes_look_like_error=true;
	    for (i=1; (i<=length_sup-1) && (interior_nodes_look_like_error==true); i++)
	      {
		if (sum_of_covgs_in_desired_colours(path_nodes[i])>coverage)
		  {
		    interior_nodes_look_like_error=false;
		  }
	      }

	    if (interior_nodes_look_like_error==true)
	      {

		for (i=1; (i<=length_sup-1); i++)
		  {

		    db_graph_db_node_prune_low_coverage(
							path_nodes[i], coverage,
							&db_node_action_set_status_pruned, db_graph,
							sum_of_covgs_in_desired_colours,
							get_edge_of_interest, apply_reset_to_specified_edges,
							apply_reset_to_specified_edges_2);
		  }
	      }
	    else//interior nodes do not look like error
	      {
		//don't prune - some interior ode has high coverage
	    	is_supernode_pruned=false;
	      }
	  
	  }
      }
  else//debug only
    {
      is_supernode_pruned=false;
    }
    return is_supernode_pruned;
}






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
											  Covg covg_thresh)
{

  boolean condition_error_more_likely(dBNode** p_nodes, int len_sp, int num_hap_chroms, 
				      double  total_dep_of_covg, int rd_len, double err_rate_per_base)
  {
    boolean too_short = false;
    Covg cov
      = count_reads_on_allele_in_specific_func_of_colours(
          p_nodes, len_sp, sum_of_covgs_in_desired_colours, &too_short);
    
    if (too_short==true)
      {
	return false;
      }

    /*
    // log of      dpois(c, D_over_R*e*k/3)  * exp(-D_over_R*e*len/3)
    double llk_cov_under_error_model = -total_dep_of_covg*err_rate_per_base*k/(3*rd_len) 
                                       + cov*log(total_dep_of_covg*err_rate_per_base*k/(3*rd_len)) 
                                       - total_dep_of_covg*err_rate_per_base*sp_len/(3*rd_len)
                                       - gsl_sf_lnfact(cov);

    // dpois(c, total_covg*f*(R-k+1)/R)*exp(-total_covg*f*len/R)
    double llk_cov_under_pop_var_model = -total_dep_of_covg*f
    */

    if ((len_sp <= db_graph->kmer_size +1 ) && (cov<=covg_thresh))
      {
	return true;
      }
    else if ((len_sp > db_graph->kmer_size +1 ) && (len_sp <= 2* (db_graph->kmer_size +1) ) && (cov<= covg_thresh+1))
      {
	return true;
      }
    else
      {
	return false;
      }


  }

  boolean is_supernode_pruned=true;

  if (node==NULL)
    {
       die("Called db_graph_remove_supernode_containing_this_node_if_more_likely_error_than_sampling  with a NULL node. Programming error\n");
    }

  if (db_node_check_status(node, none)==false)//don't touch stuff that is visited or pruned, or whatever
    {
      //printf("Ignore as status us not none\n\n");
      *supernode_len=-1;//caller can check this
      is_supernode_pruned=false;
      return is_supernode_pruned;
    }


  if (sum_of_covgs_in_desired_colours(node)>0)
      {
	//get the supernode, setting nodes to visited
	double avg_cov;
	Covg min_cov;
	Covg max_cov;
	boolean is_cycle;
	
	//length_sup is the number of edges in the supernode
	int length_sup =  db_graph_supernode_in_subgraph_defined_by_func_of_colours(node,max_expected_sup,
										    &db_node_action_set_status_visited,
										    path_nodes, path_orientations, path_labels, supernode_string,
										    &avg_cov,&min_cov, &max_cov, &is_cycle,
										    db_graph, 
										    get_edge_of_interest,
										    sum_of_covgs_in_desired_colours);

	*supernode_len=length_sup;

        if (length_sup > max_length_allowed_to_prune)
	  {
	    is_supernode_pruned=false;
	    return is_supernode_pruned;//do nothing. This is too long 
	  }
	else if (length_sup <=1)
	  {
	    is_supernode_pruned=false;
	    return is_supernode_pruned;//do nothing. This supernode has no interior, is just 1 or 2 nodes, so cannot prune it
	  }
	else if ( condition_error_more_likely(path_nodes, length_sup, num_haploid_chroms, 
					      total_depth_of_covg, read_len, error_rate_per_base)==true)
	  {
	    int i;
	    for (i=1; (i<=length_sup-1); i++)
	      {
		
		db_graph_db_node_prune_without_condition(path_nodes[i], 
							 &db_node_action_set_status_pruned,
							 db_graph,
							 get_edge_of_interest, apply_reset_to_specified_edges, apply_reset_to_specified_edges_2
							 );
	      }
	  }
      }
  else//debug only
    {
      //printf("OK - query node has too much covg to consider removing\n");
      is_supernode_pruned=false;
    }
    return is_supernode_pruned;

}




// traverse graph. At each node, if covg <= arg1, get its supernode. If that supernode length is <= kmer-length, and ALL interior nodes have covg <= arg1 
// then prune the node, and the interior nodes of the supernode.
// returns the number of pruned supernodes
long long db_graph_remove_errors_considering_covg_and_topology(
  Covg coverage, dBGraph * db_graph,
  Covg (*sum_of_covgs_in_desired_colours)(const Element *), 
  Edges (*get_edge_of_interest)(const Element*), 
  void (*apply_reset_to_specified_edges)(dBNode*, Orientation, Nucleotide), 
  void (*apply_reset_to_specified_edges_2)(dBNode*),
  int max_expected_sup)
{

  dBNode**     path_nodes        = (dBNode**) malloc(sizeof(dBNode*)* (max_expected_sup+1)); 
  Orientation* path_orientations = (Orientation*) malloc(sizeof(Orientation)*(max_expected_sup+1)); 
  Nucleotide*  path_labels       = (Nucleotide*) malloc(sizeof(Nucleotide)*(max_expected_sup+1));
  char*        supernode_string  = (char*) malloc(sizeof(char)*(max_expected_sup+2)); //2nd +1 for \0

  if ( (path_nodes==NULL) || (path_orientations==NULL) || (path_labels==NULL) || (supernode_string==NULL) )
    {
      die("Cannot malloc arrays for db_graph_remove_errors_considering_covg_and_topology");
    }


  long long prune_supernode_if_it_looks_like_is_induced_by_errors(dBNode* node)
  {

    int len;
    boolean is_sup_pruned;
    
    is_sup_pruned = db_graph_remove_supernode_containing_this_node_if_looks_like_induced_by_error(node, coverage, db_graph, max_expected_sup,
												  sum_of_covgs_in_desired_colours,
												  get_edge_of_interest,
												  apply_reset_to_specified_edges, 
												  apply_reset_to_specified_edges_2,
												  path_nodes, path_orientations, path_labels,supernode_string,&len);
    if (is_sup_pruned==true)
      {
	return 1;
      }
    else
      {
	return 0;
      }

  }
  

  long long number_of_pruned_supernodes  = hash_table_traverse_returning_sum(&prune_supernode_if_it_looks_like_is_induced_by_errors, db_graph);
  hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);


  free(path_nodes);
  free(path_orientations);
  free(path_labels);
  free(supernode_string);

  return number_of_pruned_supernodes;
}






boolean db_node_is_supernode_end(dBNode * element,Orientation orientation, int edge_index, dBGraph* db_graph)
{
  if (element==NULL)
    {
      //printf("Sending null pointer to db_node_is_supernode_end\n");
      return false;
    }
  char edges = get_edge_copy(*element, edge_index);
  
  if (orientation == reverse)
    {
      //shift along so the 4 most significant bits become the 4 least - we've passed our argument by copy so not altering the original
      edges >>= 4;
    }

  edges &= 15; // AND with 00001111 so that we only look at the 4 least significant bits 
   
  //is supernode end EITHER if it has 0 edges out
  if (edges == 0)
    {
      return true;
    }
  else if ((edges != 1) && (edges != 2) && (edges != 4) && (edges != 8))  // or if it has too many edges, ie >1
    {
      return true;
    }
  else  //or next node has more than one arrow in
    {
      Orientation next_orientation=forward;
      
      //we know this element has only one edge out. What is it? The function db_node_has_precisely_one_edge fills the answer into argument 3
      Nucleotide nuc, reverse_nuc;
      if (db_node_has_precisely_one_edge(element, orientation, &nuc, edge_index))
	{

	  dBNode* next_node = db_graph_get_next_node_for_specific_person_or_pop(element, orientation, &next_orientation, nuc, &reverse_nuc,db_graph, edge_index);

	  if ( (next_node==element) && (next_orientation==orientation) )
	    {
	      return true; //infinite self-loop
	    }
	  else if (db_node_has_precisely_one_edge(next_node, opposite_orientation(next_orientation), &nuc, edge_index))
	    {
	      return false;//successor node also has only one edge coming in
	    }
	  else //successor has multiple edges
	    {
	      return true;      
	    }
	}
    }

  //to keep compiler happy
  die("We have got the end of of the is_supernode_end function - should not reach here");
  return true;
}


// wrapper for hash_table_find, which allows you to look in the hash table
// specifically for nodes related to a specific colour

dBNode *db_graph_find_node_restricted_to_specific_person_or_population(Key key, dBGraph * hash_table, int index)
{

  dBNode *e = hash_table_find(key, hash_table);

  //ASSUMING read length is strictly greater than kmer-length
  //then you should never see a kmer (node) which is unconnected to another (even if the other is itself)
  //ie to check if this kmer is seen in a given person/pop, it is enough to check if it has an edge

  if (db_node_is_this_node_in_this_person_or_populations_graph(e, index)==false)
    {
      return NULL;
    }
  else
    {
      return e;
    }
  
}



void db_graph_traverse_with_array(void (*f)(HashTable*, Element *, int**, int, int),
                                  HashTable * hash_table, int** array,
                                  int length_of_array, int index)
{
  uint64_t i;
  uint64_t num_elements = (uint64_t)hash_table->number_buckets * hash_table->bucket_size;

  for(i = 0; i < num_elements; i++)
  {
    if(!db_node_check_status(&hash_table->table[i], unassigned))
    {
      f(hash_table, &hash_table->table[i], array, length_of_array, index);
    }
  }
}


void db_graph_traverse_with_array_of_uint64(void (*f)(HashTable*, Element *, uint64_t*, uint32_t, int),
                                            HashTable * hash_table,
                                            uint64_t* array,
                                            uint32_t length_of_array, int colour)
{
  uint64_t i;
  uint64_t num_elements = (uint64_t)hash_table->number_buckets * hash_table->bucket_size;

  for(i = 0; i < num_elements; i++)
  {
    if(!db_node_check_status(&hash_table->table[i], unassigned))
    {
      f(hash_table, &(hash_table->table[i]), array, length_of_array, colour);
    }
  }
}



void db_graph_get_covg_distribution(char* filename, dBGraph* db_graph, 
                                    int index, boolean (*condition)(dBNode* elem))
{
  uint32_t i;

  FILE* fout = fopen(filename, "w");

  if(fout == NULL)
    die("Cannot open %s\n", filename);

  uint32_t covgs_len = 10001;
  // uint64_t instead of Covg since these could get large
  uint64_t* covgs = (uint64_t*) malloc(sizeof(uint64_t) * covgs_len);

  if(covgs == NULL)
    die("Could not alloc array to hold covg distrib\n");

  for(i = 0; i < covgs_len; i++)
    covgs[i] = 0;

  void bin_covg_and_add_to_array(HashTable* htable, Element *e,
                                 uint64_t* arr, uint32_t len, int colour)
  {
    if(condition(e)==true)
    {
      uint64_t bin = MIN(e->coverage[colour], len-1);
      arr[bin]++;
    }
  }
  
  db_graph_traverse_with_array_of_uint64(&bin_covg_and_add_to_array, db_graph,
					 covgs, covgs_len, index);

  fprintf(fout, "KMER_COVG\tFREQUENCY\n");
    for(i = 0; i < covgs_len; i++)
    fprintf(fout, "%u\t%" PRIu64 "\n", i, covgs[i]);

  fclose(fout);
  free(covgs);
}






//TODO - get rid of this

dBNode* db_graph_get_next_node_in_supernode_for_specific_person_or_pop(dBNode* node, Orientation orientation, Orientation* next_orientation, int index, dBGraph* db_graph)
{
  char tmp_seq[db_graph->kmer_size+1];
  tmp_seq[db_graph->kmer_size]='\0';

  if (! (db_node_is_this_node_in_this_person_or_populations_graph(node, index)))
    {
      //printf("\nThis node is not in the graph of this person\n");
      return NULL;
    }
  else if (!db_node_check_status_not_pruned(node))
    {
      printf("ignore pruned node");
      //don't waste time with pruned nodes.
      return NULL;
    }
  else if (db_node_is_supernode_end(node, orientation, index, db_graph))
    {
      if (DEBUG)
	{
	  printf("this node is at the end of the supernode, in this orientation, so cant return the next one\n");
	}
      return NULL;
    }

  Nucleotide nucleotide_for_only_edge, reverse_nucleotide_for_only_edge;

  
  db_node_has_precisely_one_edge(node, orientation, &nucleotide_for_only_edge, index);//gives us nucleotide
  
  dBNode* next_node =  db_graph_get_next_node_for_specific_person_or_pop(node, orientation, next_orientation, nucleotide_for_only_edge,&reverse_nucleotide_for_only_edge, db_graph, index);
  
  if(next_node == NULL){
    die("dB_graph: didnt find node in hash table: %s", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size, tmp_seq));
  }

  if (DEBUG)
    {
      BinaryKmer tmp_kmer;
      printf("TRY TO ADD %c - next node %s\n",binary_nucleotide_to_char(nucleotide_for_only_edge),
	     next_orientation == forward ? binary_kmer_to_seq(element_get_kmer(next_node),db_graph->kmer_size, tmp_seq) :  
	     binary_kmer_to_seq(binary_kmer_reverse_complement(element_get_kmer(next_node),db_graph->kmer_size, &tmp_kmer),db_graph->kmer_size, tmp_seq));
      
    }
    

  //check for multiple entry edges 
  Nucleotide nuc;
  if (db_node_has_precisely_one_edge(next_node,opposite_orientation(*next_orientation), &nuc, index))
    {
    }
  else
    {
      //double check
      if (node ==NULL)
	{
	  die("programming error. returning null node when my model in my head says impossible");
	}
      return node; //we have gone as far as we can go - the next node has multiple entries. So we are now at the first node of the supernode
    }
    
    
  //loop
  if ((next_node == node) && (*next_orientation == orientation))
    {      
      //double check
      if (node ==NULL)
	{
          die("programming error. returning null node when my model in my head says impossible");
        }

	return node; //we have a kmer that loops back on itself
    }
  
  return next_node;

}




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
						   dBGraph* db_graph ) 
{

 if (length_of_arrays%2 !=0)
    {
      die(
"Must only call db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta "
"with even length_of_arrays");
    }

  if (number_of_nodes_to_load>length_of_arrays)
    {
      die("Insufficient space in arrays, length %d, to load %d nodes, in "
          "db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta", 
	        length_of_arrays, number_of_nodes_to_load);
    }


  //push everything in the arrays left by number_of_nodes_to_load. 
  int i;
  for(i=0; i<length_of_arrays-number_of_nodes_to_load; i++)
    {
      path_nodes[i]        = path_nodes[i+number_of_nodes_to_load];
      path_orientations[i] = path_orientations[i+number_of_nodes_to_load];
      path_labels[i]       = path_labels[i+number_of_nodes_to_load];
      path_string[i]       = path_string[i+number_of_nodes_to_load];
    } 
  for(i=length_of_arrays-number_of_nodes_to_load; i<length_of_arrays; i++)
    {
      path_nodes[i]        = NULL;
      path_orientations[i] = forward;
      path_labels[i]     =  Undefined;
      path_string[i]     = 'X';//so we spot a problem if we start seeing X's.
    }
  path_string[length_of_arrays-number_of_nodes_to_load]='\0';

  //move a kmer-worth of bases to front of seq, if appropriate
  if(!expecting_new_fasta_entry)
    {

      //sanity
      if (number_of_nodes_loaded_last_time==0)
	{
	  die("Not expecting a new fasta entry, but did not load any nodes last time - programming error");
	}
	
      // When about to load a new set of nodes, decide whether to preload seq with the last kmer of the previous set. Do this iff you are not expecting the start of a new entry.
      // Only tricky bit is determining the index of the last kmer within seq (from the last time around), when you want to move it to the start of seq just before loading a new set of nodes. 
      // If last time, we preloaded seq with a kmer at the start (indices 0 to k-1), and we loaded N nodes, then the last kmer in seq is from indices N to N+k-1
      // If last time we did not preload seq, and we loaded N nodes, then the last kmer is from N-1 to N+k-2.
      // This boils down to knowing whether the last time was a new fasta entry
      
      
      if (last_time_was_not_start_of_entry==true) //which will be true most of the time
	{
	  for (i=0; i<db_graph->kmer_size; i++)
	    {
	      seq->seq[i]=seq->seq[number_of_nodes_loaded_last_time+i]; 
	    }
	}
      else
	{
	  for (i=0; i<db_graph->kmer_size; i++)
	    {
	      seq->seq[i]=seq->seq[number_of_nodes_loaded_last_time-1+i]; 
	    }
	}
      seq->seq[db_graph->kmer_size]='\0';
    }

  return load_seq_into_array(chrom_fptr, number_of_nodes_to_load, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);

  
}


int db_node_addr_cmp(const void *a, const void *b)
{
  const dBNode** ia = (const dBNode **)a; // casting pointer types
  const dBNode** ib = (const dBNode **)b;
  return *ia  - *ib;
  //address comparison: returns negative if b > a
  // and positive if a > b 
  
}


void get_coverage_from_array_of_nodes(dBNode** array, int length,
                                      Covg *min_coverage, Covg *max_coverage,
                                      double* avg_coverage, Covg *mode_coverage,
                                      double* percent_nodes_having_modal_value,
                                      int index)
{
  Covg sum_coverage = 0;

  *max_coverage = 0;
  *min_coverage = COVG_MAX;

  int i;
  Covg coverages[length];

  for (i=0; i< length; i++)
    {
      Covg this_covg = 0;
      
      if (array[i]!= NULL)
	{
	  this_covg = db_node_get_coverage(array[i], index);
	}
      
      coverages[i]=this_covg; //will use this later, for the mode

      sum_coverage += this_covg;
      *max_coverage = MAX(this_covg, *max_coverage);
      *min_coverage = MIN(this_covg, *min_coverage);
      //printf("i is %d, this node has coverage %d, min is %d, max is %d\n", i, this_covg, *min_coverage, *max_coverage);

    }  

  if (*min_coverage==COVG_MAX)
    {
      *min_coverage=0;
    }
  *avg_coverage = (double)sum_coverage/length;


  qsort( coverages, length, sizeof(int), int_cmp);
  Covg covg_seen_most_often = coverages[0];
  int number_of_nodes_with_covg_seen_most_often=1;
  int current_run_of_identical_adjacent_covgs=1;

  for (i=1; i< length; i++)
    {
      if (coverages[i]==coverages[i-1])
	{
	  current_run_of_identical_adjacent_covgs++;
	}
      else
	{
	  current_run_of_identical_adjacent_covgs=1;
	}

      if (current_run_of_identical_adjacent_covgs > number_of_nodes_with_covg_seen_most_often)
	{
	  number_of_nodes_with_covg_seen_most_often = current_run_of_identical_adjacent_covgs;
	  covg_seen_most_often=coverages[i];
	}
    }

  *mode_coverage = covg_seen_most_often;
  *percent_nodes_having_modal_value = 100* number_of_nodes_with_covg_seen_most_often/length;

}


void get_coverage_from_array_of_nodes_in_subgraph_defined_by_func_of_colours(
  dBNode** array, int length,  Covg *min_coverage, Covg *max_coverage,
  double* avg_coverage, Covg *mode_coverage,
  double *percent_nodes_having_modal_value,
  Edges (*get_colour)(const dBNode*),
  Covg (*get_covg)(const dBNode*))
{
  Covg sum_coverage = 0;

  *max_coverage = 0;
  *min_coverage = COVG_MAX;

  int i;
  Covg coverages[length];

  for (i=0; i< length; i++)
    {
      Covg this_covg = 0;
      
      if (array[i]!= NULL)
	{
	  this_covg = db_node_get_coverage_in_subgraph_defined_by_func_of_colours(array[i], get_covg);
	}
      
      coverages[i]=this_covg; //will use this later, for the mode

      sum_coverage += this_covg;
      *max_coverage = MAX(this_covg, *max_coverage);
      *min_coverage = MIN(this_covg, *min_coverage);
      //printf("i is %d, this node has coverage %d, min is %d, max is %d\n", i, this_covg, *min_coverage, *max_coverage);

    }  

  if (*min_coverage==COVG_MAX)
    {
      *min_coverage=0;
    }
  *avg_coverage = sum_coverage/length;


  qsort( coverages, length, sizeof(int), int_cmp);
  Covg covg_seen_most_often = coverages[0];
  int number_of_nodes_with_covg_seen_most_often=1;
  int current_run_of_identical_adjacent_covgs=1;

  for (i=1; i< length; i++)
    {
      if (coverages[i]==coverages[i-1])
	{
	  current_run_of_identical_adjacent_covgs++;
	}
      else
	{
	  current_run_of_identical_adjacent_covgs=1;
	}

      if (current_run_of_identical_adjacent_covgs > number_of_nodes_with_covg_seen_most_often)
	{
	  number_of_nodes_with_covg_seen_most_often = current_run_of_identical_adjacent_covgs;
	  covg_seen_most_often=coverages[i];
	}
    }

  *mode_coverage = covg_seen_most_often;
  *percent_nodes_having_modal_value = 100* number_of_nodes_with_covg_seen_most_often/length;

}







boolean does_this_path_exist_in_this_colour(dBNode** array_nodes, Orientation* array_orientations,  int len, Edges (*get_colour)(const dBNode*), dBGraph* db_graph )
{
  boolean ret = true;
  int i;

  if (array_nodes[0]==NULL)
    {
      return false;
    }

  for (i=1; i<len; i++)
    {
      //get last base in kmer
      Nucleotide n;

      if (array_nodes[i]==NULL)
	{
	  return false;
	}

      if (array_orientations[i]==forward)
	{
	  n = binary_kmer_get_last_nucleotide(&(array_nodes[i]->kmer));
	}
      else
	{
	  BinaryKmer tmp_kmer;
	  n = binary_kmer_get_last_nucleotide(binary_kmer_reverse_complement(&(array_nodes[i]->kmer), db_graph->kmer_size, &tmp_kmer));
	}
      
      if (!(db_node_edge_exist_within_specified_function_of_coloured_edges(array_nodes[i-1], n, array_orientations[i-1], get_colour)) )
	{
	  ret=false;
	  break;
	}
      
    }

  return ret;
  
  
}





//check all edges in graph
long long db_graph_health_check(boolean fix, dBGraph * db_graph){
  dBNode * next_node;
  Nucleotide reverse_nucleotide;
  Orientation next_orientation;
  long long count_nodes=0;
  char tmp_seq[db_graph->kmer_size+1]; 
  tmp_seq[db_graph->kmer_size]='\0';

  void check_node_with_orientation(dBNode * node, Orientation orientation){

    int j;
    for (j=0; j< NUMBER_OF_COLOURS; j++)
      {

	void check_base(Nucleotide n)
	{     
	  if (db_node_edge_exist(node,n,orientation, j)==true){

	    next_node = db_graph_get_next_node_for_specific_person_or_pop(node,orientation,&next_orientation,n,&reverse_nucleotide,db_graph, j);
	    
	    if(next_node == NULL)
	      {
		printf("Health check problem -  didnt find node in hash table: %s %c %s\n",
		       binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,tmp_seq),
		       binary_nucleotide_to_char(n), 
		       orientation == forward ? "forward" : "reverse");
		if (fix)
		  {
		    db_node_reset_edge(node,orientation,n, j);
		  }
	      }
	    else
	      {
		if (db_node_edge_exist(next_node,reverse_nucleotide,opposite_orientation(next_orientation), j)==false)
		  {
		    printf("Health check problem - inconsitency return arrow missing: %s %c %s\n",
			   binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,tmp_seq),
			   binary_nucleotide_to_char(n), 
			   orientation == forward ? "forward" : "reverse");
		    if (fix)
		      {
			db_node_reset_edge(node,orientation, n, j);
		      }
		  }
	      }
	    
	  }
	}

	nucleotide_iterator(&check_base);
      }
  }


  BinaryKmer zerokmer;
  binary_kmer_initialise_to_zero(&zerokmer);
  uint64_t  count_num_zero_kmers= (uint64_t) 0;

  void check_node(dBNode * node){
    check_node_with_orientation(node,forward);
    check_node_with_orientation(node,reverse);
    if (binary_kmer_comparison_operator(zerokmer, node->kmer)==true)
      {
	count_num_zero_kmers++;
      }
    count_nodes++;
  }


  hash_table_traverse(&check_node,db_graph); 
  if (count_num_zero_kmers>1)
    {
      printf("Error: found %" PRIu64 " zero kmers in the graph\n", count_num_zero_kmers);
    }
  printf("%qd nodes checked\n",count_nodes);
  return count_nodes;
}


//clean off edges that point nowhere - caused when you dump a subgraph
long long db_graph_clean_orphan_edges(dBGraph * db_graph){
  dBNode * next_node;
  Nucleotide reverse_nucleotide;
  Orientation next_orientation;
  long long count_nodes=0;

  void check_node_with_orientation(dBNode * node, Orientation orientation){

    int j;
    for (j=0; j< NUMBER_OF_COLOURS; j++)
      {

	void check_base(Nucleotide n)
	{     
	  if (db_node_edge_exist(node,n,orientation, j)==true){

	    next_node = db_graph_get_next_node_for_specific_person_or_pop(node,orientation,&next_orientation,n,&reverse_nucleotide,db_graph, j);
	    
	    if(next_node == NULL)
	      {
		db_node_reset_edge(node,orientation,n, j);
	      }
	    else
	      {
		if (db_node_edge_exist(next_node,reverse_nucleotide,opposite_orientation(next_orientation), j)==false)
		  {
		    db_node_reset_edge(node,orientation, n, j);
		  }
	      }
	    
	  }
	}

	nucleotide_iterator(&check_base);
      }
  }


  void check_node(dBNode * node){
    check_node_with_orientation(node,forward);
    check_node_with_orientation(node,reverse);
    count_nodes++;
  }


  hash_table_traverse(&check_node,db_graph); 
  return count_nodes;
}





//wipes a colour clean - for all nodes, sets covg=0, edge=0.
void db_graph_wipe_colour(int colour, dBGraph* db_graph)
{
  void wipe_node(dBNode* node)
  {
    node->individual_edges[colour]=0;
    node->coverage[colour]=0;
  }
  hash_table_traverse(&wipe_node, db_graph);
}




//wipes two colours clean - for all nodes, sets covg=0, edge=0 - but only traverses the hash once
void db_graph_wipe_two_colours_in_one_traversal(int colour1, int colour2, dBGraph* db_graph)
{
  void wipe_node(dBNode* node)
  {
    node->individual_edges[colour1]=0;
    node->coverage[colour1]=0;
    node->individual_edges[colour2]=0;
    node->coverage[colour2]=0;
  }
  hash_table_traverse(&wipe_node, db_graph);
}

void db_graph_print_supernodes_defined_by_func_of_colours(char * filename_sups, char* filename_sings, int max_length, 
							  dBGraph * db_graph, Edges (*get_colour)(const dBNode*), Covg (*get_covg)(const dBNode*),
							  void (*print_extra_info)(dBNode**, Orientation*, int, FILE*)){

  boolean do_we_print_singletons=true;//singletons are supernodes consisting of ONE node.

  FILE * fout1; //file to which we will write all supernodes which are longer than 1 node in fasta format
  fout1= fopen(filename_sups, "w"); 
  if (fout1==NULL)
    {
      die("Cannot open file %s in db_graph_print_supernodes_defined_by_func_of_colours",
          filename_sups);
    }

  FILE * fout2=NULL; //file to which we will write all "singleton" supernodes, that are just  1 node, in fasta format
  if ( strcmp(filename_sings, "")==0 )
    {
      //      printf("Only printing supernodes consisting of >1 node "
      //     "(ie contigs longer than %d bases)", db_graph->kmer_size);
    
      do_we_print_singletons = false;
    }
  else
    {
      fout2= fopen(filename_sings, "w"); 
      if (fout2==NULL)
	{
	  die("Cannot open file %s in db_graph_print_supernodes_defined_by_func_of_colours",
        filename_sings);
	}
    }


  int count_nodes=0;
  
  dBNode * *    path_nodes;
  Orientation * path_orientations;
  Nucleotide *  path_labels;
  char * seq;
  boolean is_cycle;
  double avg_coverage=0;
  Covg min_covg = 0, max_covg = 0;

  path_nodes        = calloc(max_length,sizeof(dBNode*));
  path_orientations = calloc(max_length,sizeof(Orientation));
  path_labels       = calloc(max_length,sizeof(Nucleotide));
  seq               = calloc(max_length+1+db_graph->kmer_size,sizeof(char));
  
  

  long long count_kmers = 0;
  long long count_sing  = 0;

  void print_supernode(dBNode * node){
    count_kmers++;
    char name[100];

    if (db_node_check_status_none(node) == true){
      int length = db_graph_supernode_in_subgraph_defined_by_func_of_colours(node,max_length,&db_node_action_set_status_visited,
									     path_nodes,path_orientations,path_labels,
									     seq,&avg_coverage,&min_covg,&max_covg,&is_cycle,
									     db_graph, get_colour, get_covg);

      if (length>0){	
	sprintf(name,"node_%i",count_nodes);

	print_ultra_minimal_fasta_from_path(fout1,name,length,
					    path_nodes[0],path_orientations[0], seq,
					    db_graph->kmer_size,true);
	if (length==max_length){
	  printf("contig length equals max length [%i] for node_%i\n",max_length,count_nodes);
	}
	//fprintf(fout1, "extra information:\n");
	print_extra_info(path_nodes, path_orientations, length, fout1);
	count_nodes++;
      }
      else{
	count_sing++;
	if (do_we_print_singletons==true)
	  {
	    sprintf(name,"node_%qd",count_sing);
	    print_ultra_minimal_fasta_from_path(fout2,name,length,
						path_nodes[0],path_orientations[0],seq,
						db_graph->kmer_size,true);
	    //fprintf(fout2, "extra information:\n");
	    print_extra_info(path_nodes, path_orientations, length, fout2);
	  }
	

      }
    
    }
  }
  
  hash_table_traverse(&print_supernode,db_graph); 
  printf("%qd nodes visted [%qd singletons]\n",count_kmers,count_sing);

  free(path_nodes);
  free(path_orientations);
  free(path_labels);
  free(seq);
  fclose(fout1);
  if (fout2 != NULL)
    {
      fclose(fout2);
    }
}



void print_ultra_minimal_fasta_from_path(FILE *fout,
					 char * name,
					 int length,
					 dBNode * fst_node,
					 Orientation fst_orientation,
					 // dBNode * lst_node,
					 //Orientation lst_orientation,
					 char * string, //labels of paths
					 int kmer_size,
					 boolean include_first_kmer)
{

  if (fout==NULL)
    {
      die("Exiting - have passed a null file pointer to "
          "print_ultra_minimal_fasta_from_path\n");
    }
  
  
  if (include_first_kmer==false)
    {
      fprintf(fout,">%s length:%i kmer:%d\n%s\n", name, length, kmer_size, string);
    }
  else
    {
      if ( fst_node==NULL )
	{
	  die("WARNING - print_ultra_minimal_fasta_from_path command has been "
        "given a NULL node as first node, and needs to dereference it.");
	}
      else
	{
	  char fst_seq[kmer_size+1];
	  fst_seq[kmer_size]='\0';
	  BinaryKmer fst_kmer; BinaryKmer tmp_kmer;
	  
	  if (fst_orientation==reverse)
	    {
	      binary_kmer_assignment_operator(fst_kmer, *(  binary_kmer_reverse_complement(element_get_kmer(fst_node),kmer_size, &tmp_kmer) ) );
	    } 
	  else
	    {
	      binary_kmer_assignment_operator(fst_kmer, *(element_get_kmer(fst_node)) );
	    }
	  
	  binary_kmer_to_seq(&fst_kmer,kmer_size,fst_seq);
	  
	  fprintf(fout,">%s length:%i INFO:KMER=%d\n", name,length+kmer_size, kmer_size);
	  fprintf(fout,"%s", fst_seq);
	  fprintf(fout,"%s\n",string);
      
	}
      
    }

}


void print_standard_extra_supernode_info(dBNode** node_array,
                                         Orientation* or_array,
                                         int len, FILE* fout)
{
  // Let the compiler know that we're deliberately ignoring a parameter
  (void)or_array;
  
  int col;
  for (col=0; col<NUMBER_OF_COLOURS; col++)
    {
      
      fprintf(fout, "Covg in Colour %d:\n", col);
      int i;
      for (i=0; i<len; i++)
	{
	  if (node_array[i]!=NULL)
	    {
	      fprintf(fout, "%d ", node_array[i]->coverage[col]);
	    }
	  else
	    {
	      fprintf(fout, "0 ");
	    }
	}
      fprintf(fout, "\n");
      
    }
  
}


void print_no_extra_supernode_info(dBNode** node_array, Orientation* or_array,
                                   int len, FILE* fout)
{
  // Let the compiler know that we are deliberately doing nothing with our
  // paramters
  (void)node_array;
  (void)or_array;
  (void)len;
  (void)fout;
}

