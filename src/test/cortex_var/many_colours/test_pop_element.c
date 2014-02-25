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
  test_pop_element.c
*/

#include <stdlib.h>
#include <stdio.h>

#include <CUnit.h>
#include <Basic.h>

// cortex_var headers
#include "element.h"
#include "open_hash/hash_table.h"
#include "test_pop_element.h"

void test_get_edge_copy()
{

  Element* my_element=new_element();

  my_element->individual_edges[0]=1;
  if (NUMBER_OF_COLOURS>1)
    {
      my_element->individual_edges[1]=3;
    }

  Edges edges=get_edge_copy( *my_element, 0);
  CU_ASSERT(edges==1);
  if (NUMBER_OF_COLOURS>1)
    {
      edges=get_edge_copy( *my_element, 1);
      CU_ASSERT(edges==3);
    }
  free_element(&my_element);




  my_element=new_element();
 
  my_element->individual_edges[0]=2;

  if (NUMBER_OF_COLOURS>1)
    {
      my_element->individual_edges[1]=4;
    }

  edges=get_edge_copy( *my_element, 0);
  CU_ASSERT(edges==2);
  if (NUMBER_OF_COLOURS>1)
    {
      edges=get_edge_copy( *my_element, 1);
      CU_ASSERT(edges==4);
    }
  free_element(&my_element);


}


void test_increment_coverage()
{
  dBNode* e = new_element();
  db_node_increment_coverage(e,0);
  CU_ASSERT(e->coverage[0]==1);
  db_node_increment_coverage(e,0);
  CU_ASSERT(e->coverage[0]==2);
  db_node_increment_coverage(e,0);
  CU_ASSERT(e->coverage[0]==3);
  db_node_increment_coverage(e,0);
  CU_ASSERT(e->coverage[0]==4);
  db_node_increment_coverage(e,0);
  CU_ASSERT(e->coverage[0]==5);
  db_node_increment_coverage(e,0);
  CU_ASSERT(e->coverage[0]==6);

  if (NUMBER_OF_COLOURS>1)
    {
      CU_ASSERT(e->coverage[1]==0);
      db_node_increment_coverage(e,1);
      CU_ASSERT(e->coverage[1]==1);
    }
  if (NUMBER_OF_COLOURS>2)
    {
      CU_ASSERT(e->coverage[2]==0);
    }

 
  
  free_element(&e);
}

void test_get_coverage()
{
  dBNode* e = new_element();
  db_node_increment_coverage(e,0);
  CU_ASSERT(db_node_get_coverage(e,0)==1);
  db_node_increment_coverage(e,0);
  CU_ASSERT(db_node_get_coverage(e,0)==2);
  db_node_increment_coverage(e,0);
  CU_ASSERT(db_node_get_coverage(e,0)==3);
  db_node_increment_coverage(e,0);
  CU_ASSERT(db_node_get_coverage(e,0)==4);
  db_node_increment_coverage(e,0);
  CU_ASSERT(db_node_get_coverage(e,0)==5);
  
  e->coverage[0]=120;
  CU_ASSERT(db_node_get_coverage(e,0)==120);
  
  free_element(&e);
}


void test_element_status_set_and_checks()
{
  dBNode* e = new_element();
  
  db_node_set_status(e, visited);
  CU_ASSERT(db_node_check_status(e, visited)==true);
  CU_ASSERT(db_node_check_status(e, none)==false);
  CU_ASSERT(db_node_check_status(e, pruned) ==false);

  db_node_set_status(e, none);
  CU_ASSERT(db_node_check_status(e, none)==true);
  CU_ASSERT(db_node_check_status(e, visited)==false);
  CU_ASSERT(db_node_check_status(e, pruned)==false);



  db_node_set_status(e, pruned);
  CU_ASSERT(db_node_check_status(e, pruned)==true);
  CU_ASSERT(db_node_check_status(e, none)==false);
  CU_ASSERT(db_node_check_status(e, visited)==false);


  db_node_set_status_to_none(e);
  CU_ASSERT(db_node_check_status(e, none)==true);


  db_node_action_set_status_none(e);
  CU_ASSERT(db_node_check_status(e, none)==true);

  db_node_action_set_status_visited(e);
  CU_ASSERT(db_node_check_status(e, visited)==true);
  db_node_action_set_status_visited_or_visited_and_exists_in_reference(e);
  CU_ASSERT(db_node_check_status(e, visited)==true);
  db_node_action_unset_status_visited_or_visited_and_exists_in_reference(e);
  CU_ASSERT(db_node_check_status(e, none)==true);

  
  db_node_set_status(e,exists_in_reference);
  CU_ASSERT(db_node_check_status(e, exists_in_reference)==true);
  db_node_action_set_status_visited_or_visited_and_exists_in_reference(e);
  CU_ASSERT(db_node_check_status(e, visited_and_exists_in_reference)==true);
  db_node_action_unset_status_visited_or_visited_and_exists_in_reference(e);
  CU_ASSERT(db_node_check_status(e, exists_in_reference)==true);


  db_node_set_status(e,exists_in_reference);
  CU_ASSERT(!db_node_check_status_is_not_exists_in_reference(e));
  db_node_set_status(e,visited_and_exists_in_reference);
  CU_ASSERT(db_node_check_status_is_not_exists_in_reference(e)); //note we are being sematically strict here - the STATUS is not exists_in_reference, even though the kmer IS in the reference
  db_node_set_status(e,none);
  CU_ASSERT(db_node_check_status_is_not_exists_in_reference(e));
  db_node_set_status(e,visited);
  CU_ASSERT(db_node_check_status_is_not_exists_in_reference(e));
  db_node_set_status(e,pruned);
  CU_ASSERT(db_node_check_status_is_not_exists_in_reference(e));

  db_node_set_status_to_none(e);
  CU_ASSERT(db_node_check_status_is_not_visited_or_visited_and_exists_in_reference(e));
  db_node_set_status(e,visited);
  CU_ASSERT(!db_node_check_status_is_not_visited_or_visited_and_exists_in_reference(e));
  db_node_set_status(e, visited_and_exists_in_reference);
  CU_ASSERT(!db_node_check_status_is_not_visited_or_visited_and_exists_in_reference(e));
  

  free_element(&e);


  
}


void test_element_assign()
{
  Element e1, e2;

  BinaryKmer b1;
  binary_kmer_initialise_to_zero(&b1);
  b1[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] = (bitfield_of_64bits) 1;

  BinaryKmer b2;
  binary_kmer_initialise_to_zero(&b2);
  b2[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] = (bitfield_of_64bits) 3;

  
  element_initialise(&e1, &b1, 31);
  element_initialise(&e2, &b2, 31);

  set_edges(&e1, 0, 1);
  set_edges(&e2, 0, 2);
  
  db_node_set_status(&e1, visited);
  db_node_set_status(&e2, pruned);

  db_node_update_coverage(&e1, 0, 101);
  db_node_update_coverage(&e2, 0, 202);
  
  element_assign(&e1, &e2);

  CU_ASSERT(db_node_get_coverage(&e1,0)==202);
  CU_ASSERT(get_edge_copy(e1,0)==2 );
  CU_ASSERT(binary_kmer_comparison_operator(e1.kmer,b2) );
  CU_ASSERT(e1.status==pruned );
}
