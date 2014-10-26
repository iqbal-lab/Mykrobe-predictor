/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 * 
 * CORTEX project contacts:  
 *              M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 *              Z. Iqbal (zam@well.ox.ac.uk)
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
  maths.h
*/

#ifndef MATHS_H_
#define MATHS_H_

#include "global.h"
#include "element.h"
#include "dB_graph.h"

float log_factorial(int number);
float log_factorial_ll(long long number);
int min_of_ints(int a, int b);
int max_of_ints(int a, int b);
uint64_t calculate_mean_uint64_t(uint64_t* array, uint64_t len);
long long calculate_mean(long long* array, long long len);

//void set_int_array_to_zero(int* array, int len);

#endif /* MATHS_H_ */
