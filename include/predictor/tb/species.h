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
  species.h
*/
#include "dB_graph.h"
typedef enum 
 {
   Mtuberculosis = 0,
   Mafricanum = 1,
   Mcanettii = 2,
   Mmicroti = 3,
   Mpinnipedii = 4,
   Mbovis = 5,
   Mcaprae = 6,
   NTM = 7,
  } Myc_species ;


void map_species_enum_to_str(Myc_species sp, StrBuf* sbuf);
Myc_species get_species(dBGraph *db_graph,int max_gene_len, StrBuf* install_dir,
			  int ignore_first, int ignore_last);
