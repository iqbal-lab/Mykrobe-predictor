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
  known_mutations.h
*/

#include "known_mutations.h"
#include <string.h>
#include "global.h"

GeneMutationGene map_gene_name_str_to_genename(StrBuf* name)
{
  if (strcmp(name->buff, "dfrB")==0)
    {
      return dfrB;
    }
  else if (strcmp(name->buff, "fusA")==0)
    {
      return fusA;
    }
  else if (strcmp(name->buff, "rpoB")==0)
    {
      return rpoB;
    }
  else if (strcmp(name->buff, "gyrA")==0)
    {
      return gyrA;
    }
  else if (strcmp(name->buff, "grlA")==0)
    {
      return grlA;
    }
  else if (strcmp(name->buff, "grlB")==0)
    {
      return grlB;
    }
  else
    {
      die("Failed to parse gene name - got %s\n", name->buff);
    }

}
KnownMutation map_mutation_name_to_enum(StrBuf* sbuf, GeneMutationGene gene)
{
  if ( (strcmp(sbuf->buff, "L21V")==0) && (gene==dfrB) )
    {
      return dfrB_L21V;
    } 
  else 
    {
      die("Parsing error - unknown mutation %s\n", sbuf->buff);
      return NotSpecified;
    } 

}
