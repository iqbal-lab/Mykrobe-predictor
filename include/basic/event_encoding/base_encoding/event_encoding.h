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
  base_encoding/event_encoding.h
*/

#ifndef EVENT_ENCODING_H_
#define EVENT_ENCODING_H_

#include "global.h"
#include "event_encoding.h"

typedef enum
 {
    Adenine   = 0,
    Cytosine  = 1,
    Guanine   = 2,
    Thymine   = 3,
    Undefined = 4,
  } Nucleotide ;




Nucleotide char_to_binary_nucleotide(char c);
char binary_nucleotide_to_char(Nucleotide n);
Nucleotide reverse_binary_nucleotide(Nucleotide n);

boolean good_symbol(char c);

#endif /* EVENT_ENCODING_H_ */
