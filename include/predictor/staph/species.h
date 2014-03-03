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
  Scapitis = 0,
  Scaprae = 1,
  Sepidermidis = 2,
  Sequorum =3,
  Shaemolyticus = 4,
  Shominis = 5,
  Slugdunensis = 6,
  Smassiliensis = 7,
  Spettenkofer = 8,
  Spseudintermedius = 9,
  Ssaprophyticus = 10,
  Ssimiae = 11,
  Ssimulans = 12,
  Ssphgb0015 = 13,
  Sspoj82 = 14,
  Aureus = 15,
  Swarneri = 16,
  } Staph_species ;


Staph_species get_species(dBGraph *db_graph,int max_gene_len );