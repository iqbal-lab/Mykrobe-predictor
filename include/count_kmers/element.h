/*
 * 
 * CORTEX project contacts:  
 *    M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 *    Z. Iqbal (zam@well.ox.ac.uk)
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
  element.h defines the interface to count kmers in a set of reads
*/

#ifndef ELEMENT_H_
#define ELEMENT_H_

#include <stdio.h>

#include "global.h"
#include "binary_kmer.h"

typedef enum{
  unassigned = 0,
  none       = 1
} Status;


typedef struct{
	BinaryKmer kmer;
	Status status;
	long long count;
} Element;

typedef BinaryKmer Key;

boolean element_smaller(Element,Element);

boolean db_node_check_status(Element *, Status status);

boolean element_is_key(Key,Element,short kmer_size);

Key element_get_key(BinaryKmer,short kmer_size);

void element_initialise(Element *,Key, short kmer_size);

void element_increment_count(Element*);

void element_print(FILE * file, Element* e,short kmer_size, char * prefix);

#endif /* ELEMENT_H_ */
