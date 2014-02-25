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
  test_pop_supernode_consensus.h
*/

#ifndef TEST_POP_SUPERNODE_CONSENSUS_H_
#define TEST_POP_SUPERNODE_CONSENSUS_H_

#include "global.h"

void test_find_first_node_in_supernode();
void test_find_next_node_in_supernode();
void test_correctly_find_subsection_of_supernode();
void test_find_best_subsection_of_supernode_with_just_two_people();
void test_get_population_consensus_supernode();

#endif /* TEST_POP_SUPERNODE_CONSENSUS_H_ */
