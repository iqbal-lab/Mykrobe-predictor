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
  if ( (strcmp(sbuf->buff, "F99Y")==0) && (gene==dfrB) )
    {
      return dfrB_F99Y;
    }
  else if ( (strcmp(sbuf->buff, "F99S")==0) && (gene==dfrB) )
    {
      return dfrB_F99S;
    } 
  else if ( (strcmp(sbuf->buff, "F99I")==0) && (gene==dfrB) )
    {
      return dfrB_F99I;
    } 
  else if ( (strcmp(sbuf->buff, "H31N")==0) && (gene==dfrB) )
    {
      return dfrB_H31N;
    } 
  else if ( (strcmp(sbuf->buff, "L41F")==0) && (gene==dfrB) )
    {
      return dfrB_L41F;
    } 
  else if ( (strcmp(sbuf->buff, "H150R")==0) && (gene==dfrB) )
    {
      return dfrB_H150R;
    } 
  else if ( (strcmp(sbuf->buff, "L21V")==0) && (gene==dfrB) )
    {
      return dfrB_L21V;
    } 
  else if ( (strcmp(sbuf->buff, "N60I")==0) && (gene==dfrB) )
    {
      return dfrB_N60I;
    } 

  else if ( (strcmp(sbuf->buff, "A655E")==0) && (gene==fusA) )
    {
      return fusA_A655E;
    } 
  else if ( (strcmp(sbuf->buff, "B434N")==0) && (gene==fusA) )
    {
      return fusA_B434N;
    } 
  else if ( (strcmp(sbuf->buff, "F441Y")==0) && (gene==fusA) )
    {
      return fusA_F441Y;
    } 
  else if ( (strcmp(sbuf->buff, "E444K")==0) && (gene==fusA) )
    {
      return fusA_E444K;
    } 
  else if ( (strcmp(sbuf->buff, "G451V")==0) && (gene==fusA) )
    {
      return fusA_G451V;
    } 
  else if ( (strcmp(sbuf->buff, "G452S")==0) && (gene==fusA) )
    {
      return fusA_G452S;
    } 
  else if ( (strcmp(sbuf->buff, "G452C")==0) && (gene==fusA) )
    {
      return fusA_G452C;
    } 
  else if ( (strcmp(sbuf->buff, "G617D")==0) && (gene==fusA) )
    {
      return fusA_G617D;
    } 
  else if ( (strcmp(sbuf->buff, "G664S")==0) && (gene==fusA) )
    {
      return fusA_G664S;
    } 
  else if ( (strcmp(sbuf->buff, "H438N")==0) && (gene==fusA) )
    {
      return fusA_H438N;
    } 
  else if ( (strcmp(sbuf->buff, "L456F")==0) && (gene==fusA) )
    {
      return fusA_L456F;
    } 
  else if ( (strcmp(sbuf->buff, "H457Y")==0) && (gene==fusA) )
    {
      return fusA_H457Y;
    } 
  else if ( (strcmp(sbuf->buff, "L461K")==0) && (gene==fusA) )
    {
      return fusA_L461K;
    } 
  else if ( (strcmp(sbuf->buff, "L461S")==0) && (gene==fusA) )
    {
      return fusA_L461S;
    } 
  else if ( (strcmp(sbuf->buff, "M453I")==0) && (gene==fusA) )
    {
      return fusA_M453I;
    } 
  else if ( (strcmp(sbuf->buff, "M651I")==0) && (gene==fusA) )
    {
      return fusA_M651I;
    } 
  else if ( (strcmp(sbuf->buff, "P114H")==0) && (gene==fusA) )
    {
      return fusA_P114H;
    } 
  else if ( (strcmp(sbuf->buff, "P404L")==0) && (gene==fusA) )
    {
      return fusA_P404L;
    } 
  else if ( (strcmp(sbuf->buff, "P404Q")==0) && (gene==fusA) )
    {
      return fusA_P404Q;
    } 
  else if ( (strcmp(sbuf->buff, "P406L")==0) && (gene==fusA) )
    {
      return fusA_P406L;
    } 
  else if ( (strcmp(sbuf->buff, "P478S")==0) && (gene==fusA) )
    {
      return fusA_P478S;
    } 
  else if ( (strcmp(sbuf->buff, "Q115L")==0) && (gene==fusA) )
    {
      return fusA_Q115L;
    } 
  else if ( (strcmp(sbuf->buff, "R464C")==0) && (gene==fusA) )
    {
      return fusA_R464C;
    } 
  else if ( (strcmp(sbuf->buff, "R464S")==0) && (gene==fusA) )
    {
      return fusA_R464S;
    } 
  else if ( (strcmp(sbuf->buff, "R464H")==0) && (gene==fusA) )
    {
      return fusA_R464H;
    } 
  else if ( (strcmp(sbuf->buff, "R659C")==0) && (gene==fusA) )
    {
      return fusA_R659C;
    } 
  else if ( (strcmp(sbuf->buff, "R659H")==0) && (gene==fusA) )
    {
      return fusA_R659H;
    } 
  else if ( (strcmp(sbuf->buff, "R659S")==0) && (gene==fusA) )
    {
      return fusA_R659S;
    } 
  else if ( (strcmp(sbuf->buff, "R659L")==0) && (gene==fusA) )
    {
      return fusA_R659L;
    } 
  else if ( (strcmp(sbuf->buff, "T385N")==0) && (gene==fusA) )
    {
      return fusA_T385N;
    } 
  else if ( (strcmp(sbuf->buff, "T436I")==0) && (gene==fusA) )
    {
      return fusA_T436I;
    } 
  else if ( (strcmp(sbuf->buff, "T656K")==0) && (gene==fusA) )
    {
      return fusA_T656K;
    } 
  else if ( (strcmp(sbuf->buff, "V607I")==0) && (gene==fusA) )
    {
      return fusA_V607I;
    } 
  else if ( (strcmp(sbuf->buff, "V90A")==0) && (gene==fusA) )
    {
      return fusA_V90A;
    } 
  else if ( (strcmp(sbuf->buff, "V90I")==0) && (gene==fusA) )
    {
      return fusA_V90I;
    } 
  else if ( (strcmp(sbuf->buff, "L456F")==0) && (gene==fusA) )
    {
      return fusA_L456F;
    } 
  else if ( (strcmp(sbuf->buff, "A376V")==0) && (gene==fusA) )
    {
      return fusA_A376V;
    } 
  else if ( (strcmp(sbuf->buff, "A655P")==0) && (gene==fusA) )
    {
      return fusA_A655P;
    } 
  else if ( (strcmp(sbuf->buff, "D463G")==0) && (gene==fusA) )
    {
      return fusA_D463G;
    } 
  else if ( (strcmp(sbuf->buff, "D189G")==0) && (gene==fusA) )
    {
      return fusA_D189G;
    } 
  else if ( (strcmp(sbuf->buff, "A70V")==0) && (gene==fusA) )
    {
      return fusA_A70V;
    } 
  else if ( (strcmp(sbuf->buff, "A160V")==0) && (gene==fusA) )
    {
      return fusA_A160V;
    } 
  else if ( (strcmp(sbuf->buff, "A71V")==0) && (gene==fusA) )
    {
      return fusA_A71V;
    } 
  else if ( (strcmp(sbuf->buff, "D189V")==0) && (gene==fusA) )
    {
      return fusA_D189V;
    } 
  else if ( (strcmp(sbuf->buff, "L430S")==0) && (gene==fusA) )
    {
      return fusA_L430S;
    } 
  else if ( (strcmp(sbuf->buff, "E444V")==0) && (gene==fusA) )
    {
      return fusA_E444V;
    } 
  else if ( (strcmp(sbuf->buff, "F652S")==0) && (gene==fusA) )
    {
      return fusA_F652S;
    } 
  else if ( (strcmp(sbuf->buff, "Y654N")==0) && (gene==fusA) )
    {
      return fusA_Y654N;
    } 
  else if ( (strcmp(sbuf->buff, "L461F")==0) && (gene==fusA) )
    {
      return fusA_L461F;
    } 
  else if ( (strcmp(sbuf->buff, "M161I")==0) && (gene==fusA) )
    {
      return fusA_M161I;
    } 
  else if ( (strcmp(sbuf->buff, "A67T")==0) && (gene==fusA) )
    {
      return fusA_A67T;
    } 
  else if ( (strcmp(sbuf->buff, "S416F")==0) && (gene==fusA) )
    {
      return fusA_S416F;
    } 
  else if ( (strcmp(sbuf->buff, "G556S")==0) && (gene==fusA) )
    {
      return fusA_G556S;
    } 
  else if ( (strcmp(sbuf->buff, "H557Y")==0) && (gene==fusA) )
    {
      return fusA_H557Y;
    } 
  else if ( (strcmp(sbuf->buff, "D373N")==0) && (gene==fusA) )
    {
      return fusA_D373N;
    } 
  else if ( (strcmp(sbuf->buff, "C473S")==0) && (gene==fusA) )
    {
      return fusA_C473S;
    } 
  else if ( (strcmp(sbuf->buff, "R76C")==0) && (gene==fusA) )
    {
      return fusA_R76C;
    } 
  else if ( (strcmp(sbuf->buff, "T387I")==0) && (gene==fusA) )
    {
      return fusA_T387I;
    } 
  else if ( (strcmp(sbuf->buff, "E449K")==0) && (gene==fusA) )
    {
      return fusA_E449K;
    } 
  else if ( (strcmp(sbuf->buff, "H457Q")==0) && (gene==fusA) )
    {
      return fusA_H457Q;
    } 
  else if ( (strcmp(sbuf->buff, "A655V")==0) && (gene==fusA) )
    {
      return fusA_A655V;
    } 
  else if ( (strcmp(sbuf->buff, "V90I")==0) && (gene==fusA) )
    {
      return fusA_V90I;
    } 
  else if ( (strcmp(sbuf->buff, "E233Q")==0) && (gene==fusA) )
    {
      return fusA_E233Q;
    } 

  else if ( (strcmp(sbuf->buff, "Q456K")==0) && (gene==rpoB) )
    {
      return rpoB_Q456K;
    } 
  else if ( (strcmp(sbuf->buff, "S463P")==0) && (gene==rpoB) )
    {
      return rpoB_S463P;
    } 
  else if ( (strcmp(sbuf->buff, "S464P")==0) && (gene==rpoB) )
    {
      return rpoB_S464P;
    } 
  else if ( (strcmp(sbuf->buff, "L466S")==0) && (gene==rpoB) )
    {
      return rpoB_L466S;
    } 
  else if ( (strcmp(sbuf->buff, "Q468L")==0) && (gene==rpoB) )
    {
      return rpoB_Q468L;
    } 
  else if ( (strcmp(sbuf->buff, "Q468R")==0) && (gene==rpoB) )
    {
      return rpoB_Q468R;
    } 
  else if ( (strcmp(sbuf->buff, "Q468K")==0) && (gene==rpoB) )
    {
      return rpoB_Q468K;
    } 
  else if ( (strcmp(sbuf->buff, "M470T")==0) && (gene==rpoB) )
    {
      return rpoB_M470T;
    } 
  else if ( (strcmp(sbuf->buff, "D471G")==0) && (gene==rpoB) )
    {
      return rpoB_D471G;
    } 
  else if ( (strcmp(sbuf->buff, "D471Y")==0) && (gene==rpoB) )
    {
      return rpoB_D471Y;
    } 
  else if ( (strcmp(sbuf->buff, "A473T")==0) && (gene==rpoB) )
    {
      return rpoB_A473T;
    } 
  else if ( (strcmp(sbuf->buff, "N474K")==0) && (gene==rpoB) )
    {
      return rpoB_N474K;
    } 
  else if ( (strcmp(sbuf->buff, "ins475G")==0) && (gene==rpoB) )
    {
      return rpoB_ins475G;
    } 
  else if ( (strcmp(sbuf->buff, "ins475H")==0) && (gene==rpoB) )
    {
      return rpoB_ins475H;
    } 
  else if ( (strcmp(sbuf->buff, "A477T")==0) && (gene==rpoB) )
    {
      return rpoB_A477T;
    } 
  else if ( (strcmp(sbuf->buff, "A477V")==0) && (gene==rpoB) )
    {
      return rpoB_A477V;
    } 
  else if ( (strcmp(sbuf->buff, "A477D")==0) && (gene==rpoB) )
    {
      return rpoB_A477D;
    } 
  else if ( (strcmp(sbuf->buff, "H481N")==0) && (gene==rpoB) )
    {
      return rpoB_H481N;
    } 
  else if ( (strcmp(sbuf->buff, "H481Y")==0) && (gene==rpoB) )
    {
      return rpoB_H481Y;
    } 
  else if ( (strcmp(sbuf->buff, "H481D")==0) && (gene==rpoB) )
    {
      return rpoB_H481D;
    } 
  else if ( (strcmp(sbuf->buff, "R484H")==0) && (gene==rpoB) )
    {
      return rpoB_R484H;
    } 
  else if ( (strcmp(sbuf->buff, "S486L")==0) && (gene==rpoB) )
    {
      return rpoB_S486L;
    } 
  else if ( (strcmp(sbuf->buff, "I527F")==0) && (gene==rpoB) )
    {
      return rpoB_I527F;
    } 
  else if ( (strcmp(sbuf->buff, "I527M")==0) && (gene==rpoB) )
    {
      return rpoB_I527M;
    } 
  else if ( (strcmp(sbuf->buff, "I527L")==0) && (gene==rpoB) )
    {
      return rpoB_I527L;
    } 
  else if ( (strcmp(sbuf->buff, "S529L")==0) && (gene==rpoB) )
    {
      return rpoB_S529L;
    } 
  else if ( (strcmp(sbuf->buff, "D550G")==0) && (gene==rpoB) )
    {
      return rpoB_D550G;
    } 
  else if ( (strcmp(sbuf->buff, "Q565R")==0) && (gene==rpoB) )
    {
      return rpoB_Q565R;
    } 


  else if ( (strcmp(sbuf->buff, "S84A")==0) && (gene==gyrA) )
    {
      return gyrA_S84A;
    } 
  else if ( (strcmp(sbuf->buff, "S84L")==0) && (gene==gyrA) )
    {
      return gyrA_S84L;
    } 
  else if ( (strcmp(sbuf->buff, "S85P")==0) && (gene==gyrA) )
    {
      return gyrA_S85P;
    } 
  else if ( (strcmp(sbuf->buff, "E88V")==0) && (gene==gyrA) )
    {
      return gyrA_E88V;
    } 
  else if ( (strcmp(sbuf->buff, "E88G")==0) && (gene==gyrA) )
    {
      return gyrA_E88G;
    } 
  else if ( (strcmp(sbuf->buff, "E88L")==0) && (gene==gyrA) )
    {
      return gyrA_E88L;
    } 
  else if ( (strcmp(sbuf->buff, "E88K")==0) && (gene==gyrA) )
    {
      return gyrA_E88K;
    } 
  else if ( (strcmp(sbuf->buff, "G106D")==0) && (gene==gyrA) )
    {
      return gyrA_G106D;
    } 

  else if ( (strcmp(sbuf->buff, "V41G")==0) && (gene==grlA) )
    {
      return grlA_V41G;
    } 
  else if ( (strcmp(sbuf->buff, "I45M")==0) && (gene==grlA) )
    {
      return grlA_I45M;
    } 
  else if ( (strcmp(sbuf->buff, "A48T")==0) && (gene==grlA) )
    {
      return grlA_A48T;
    } 
  else if ( (strcmp(sbuf->buff, "D79V")==0) && (gene==grlA) )
    {
      return grlA_D79V;
    } 
  else if ( (strcmp(sbuf->buff, "S80F")==0) && (gene==grlA) )
    {
      return grlA_S80F;
    } 
  else if ( (strcmp(sbuf->buff, "S80Y")==0) && (gene==grlA) )
    {
      return grlA_S80Y;
    } 
  else if ( (strcmp(sbuf->buff, "E84L")==0) && (gene==grlA) )
    {
      return grlA_E84L;
    } 
  else if ( (strcmp(sbuf->buff, "E84G")==0) && (gene==grlA) )
    {
      return grlA_E84G;
    } 
  else if ( (strcmp(sbuf->buff, "E84K")==0) && (gene==grlA) )
    {
      return grlA_E84K;
    } 
  else if ( (strcmp(sbuf->buff, "E84V")==0) && (gene==grlA) )
    {
      return grlA_E84V;
    } 
  else if ( (strcmp(sbuf->buff, "Y83N")==0) && (gene==grlA) )
    {
      return grlA_Y83N;
    } 
  else if ( (strcmp(sbuf->buff, "S108N")==0) && (gene==grlA) )
    {
      return grlA_S108N;
    } 
  else if ( (strcmp(sbuf->buff, "A116E")==0) && (gene==grlA) )
    {
      return grlA_A116E;
    } 
  else if ( (strcmp(sbuf->buff, "D432G")==0) && (gene==grlA) )
    {
      return grlA_D432G;
    } 
  else if ( (strcmp(sbuf->buff, "E422D")==0) && (gene==grlB) )
    {
      return grlB_E422D;
    } 
  else if ( (strcmp(sbuf->buff, "D432V")==0) && (gene==grlB) )
    {
      return grlB_D432V;
    } 
  else if ( (strcmp(sbuf->buff, "D432NorH")==0) && (gene==grlB) )
    {
      return grlB_D432NorH;
    } 
  else if ( (strcmp(sbuf->buff, "D443E")==0) && (gene==grlB) )
    {
      return grlB_D443E;
    } 
  else if ( (strcmp(sbuf->buff, "R444S")==0) && (gene==grlB) )
    {
      return grlB_R444S;
    } 
  else if ( (strcmp(sbuf->buff, "P451S")==0) && (gene==grlB) )
    {
      return grlB_P451S;
    } 
  else if ( (strcmp(sbuf->buff, "N470D")==0) && (gene==grlB) )
    {
      return grlB_N470D;
    } 
  else if ( (strcmp(sbuf->buff, "P585S")==0) && (gene==grlB) )
    {
      return grlB_P585S;
    } 
  else 
    {
      die("Parsing error\n");
      return NotSpecified;
    } 

}
