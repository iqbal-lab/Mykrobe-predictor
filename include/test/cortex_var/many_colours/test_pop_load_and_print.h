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
  test_pop_load_and_print.h
*/

#ifndef TEST_POP_LOAD_AND_PRINT_H_
#define TEST_POP_LOAD_AND_PRINT_H_

#include "global.h"

void test_load_two_people_in_same_populations_and_print_separately_their_supernodes();
void test_take_three_people_each_with_one_read_and_find_variants();
void test_take_two_people_sharing_an_alu_and_find_supernodes();
void test_printing_supernode_with_chromosome_intersections_simple();
void  test_printing_supernode_with_chromosome_intersections_simple_alu_example();
void  test_printing_supernode_with_chromosome_intersections_simple_alu_example_2();
void test_printing_of_supernode_that_might_be_an_inversion_simple();

#endif /* TEST_POP_LOAD_AND_PRINT_H_ */
