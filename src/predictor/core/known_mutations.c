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

#ifdef STAPH

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
  else if ( (strcmp(sbuf->buff, "H31N")==0) && (gene==dfrB) )
    {
      return dfrB_H31N;
    } 
  else if ( (strcmp(sbuf->buff, "L41F")==0) && (gene==dfrB) )
    {
      return dfrB_L41F;
    } 
  else if ( (strcmp(sbuf->buff, "N60I")==0) && (gene==dfrB) )
    {
      return dfrB_N60I;
    } 
  else if ( (strcmp(sbuf->buff, "F99Y")==0) && (gene==dfrB) )
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
  else if ( (strcmp(sbuf->buff, "H150R")==0) && (gene==dfrB) )
    {
      return dfrB_H150R;
    } 
  else if ( (strcmp(sbuf->buff, "V90I")==0) && (gene==fusA) )
    {
      return fusA_V90I;
    } 
  else if ( (strcmp(sbuf->buff, "P114H")==0) && (gene==fusA) )
    {
      return fusA_P114H;
    } 
  else if ( (strcmp(sbuf->buff, "Q115L")==0) && (gene==fusA) )
    {
      return fusA_Q115L;
    } 
  else if ( (strcmp(sbuf->buff, "T326I")==0) && (gene==fusA) )
    {
      return fusA_T326I;
    } 
  else if ( (strcmp(sbuf->buff, "E468V")==0) && (gene==fusA) )
    {
      return fusA_E468V;
    } 
  else if ( (strcmp(sbuf->buff, "A376V")==0) && (gene==fusA) )
    {
      return fusA_A376V;
    } 
  else if ( (strcmp(sbuf->buff, "T385N")==0) && (gene==fusA) )
    {
      return fusA_T385N;
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
  else if ( (strcmp(sbuf->buff, "D434N")==0) && (gene==fusA) )
    {
      return fusA_D434N;
    } 
  else if ( (strcmp(sbuf->buff, "T436I")==0) && (gene==fusA) )
    {
      return fusA_T436I;
    } 
  else if ( (strcmp(sbuf->buff, "H438N")==0) && (gene==fusA) )
    {
      return fusA_H438N;
    } 
  else if ( (strcmp(sbuf->buff, "E444K")==0) && (gene==fusA) )
    {
      return fusA_E444K;
    } 
  else if ( (strcmp(sbuf->buff, "E444V")==0) && (gene==fusA) )
    {
      return fusA_E444V;
    } 
  else if ( (strcmp(sbuf->buff, "G451V")==0) && (gene==fusA) )
    {
      return fusA_G451V;
    } 
  else if ( (strcmp(sbuf->buff, "G452C")==0) && (gene==fusA) )
    {
      return fusA_G452C;
    } 
  else if ( (strcmp(sbuf->buff, "G452S")==0) && (gene==fusA) )
    {
      return fusA_G452S;
    } 
  else if ( (strcmp(sbuf->buff, "M453I")==0) && (gene==fusA) )
    {
      return fusA_M453I;
    } 
  else if ( (strcmp(sbuf->buff, "L456F")==0) && (gene==fusA) )
    {
      return fusA_L456F;
    } 
  else if ( (strcmp(sbuf->buff, "H457Q")==0) && (gene==fusA) )
    {
      return fusA_H457Q;
    } 
  else if ( (strcmp(sbuf->buff, "H457Y")==0) && (gene==fusA) )
    {
      return fusA_H457Y;
    } 
  else if ( (strcmp(sbuf->buff, "L461F")==0) && (gene==fusA) )
    {
      return fusA_L461F;
    } 
  else if ( (strcmp(sbuf->buff, "L461K")==0) && (gene==fusA) )
    {
      return fusA_L461K;
    } 
  else if ( (strcmp(sbuf->buff, "L461S")==0) && (gene==fusA) )
    {
      return fusA_L461S;
    } 
  else if ( (strcmp(sbuf->buff, "D463G")==0) && (gene==fusA) )
    {
      return fusA_D463G;
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
  else if ( (strcmp(sbuf->buff, "P478S")==0) && (gene==fusA) )
    {
      return fusA_P478S;
    } 
  else if ( (strcmp(sbuf->buff, "G556S")==0) && (gene==fusA) )
    {
      return fusA_G556S;
    } 
  else if ( (strcmp(sbuf->buff, "G617D")==0) && (gene==fusA) )
    {
      return fusA_G617D;
    } 
  else if ( (strcmp(sbuf->buff, "M651I")==0) && (gene==fusA) )
    {
      return fusA_M651I;
    } 
  else if ( (strcmp(sbuf->buff, "F652S")==0) && (gene==fusA) )
    {
      return fusA_F652S;
    } 
  else if ( (strcmp(sbuf->buff, "Y654N")==0) && (gene==fusA) )
    {
      return fusA_Y654N;
    } 
  else if ( (strcmp(sbuf->buff, "A655E")==0) && (gene==fusA) )
    {
      return fusA_A655E;
    } 
  else if ( (strcmp(sbuf->buff, "A655P")==0) && (gene==fusA) )
    {
      return fusA_A655P;
    } 
  else if ( (strcmp(sbuf->buff, "T656K")==0) && (gene==fusA) )
    {
      return fusA_T656K;
    } 
  else if ( (strcmp(sbuf->buff, "R659C")==0) && (gene==fusA) )
    {
      return fusA_R659C;
    } 
  else if ( (strcmp(sbuf->buff, "R659H")==0) && (gene==fusA) )
    {
      return fusA_R659H;
    } 
  else if ( (strcmp(sbuf->buff, "R659L")==0) && (gene==fusA) )
    {
      return fusA_R659L;
    } 
  else if ( (strcmp(sbuf->buff, "R659S")==0) && (gene==fusA) )
    {
      return fusA_R659S;
    } 
  else if ( (strcmp(sbuf->buff, "G664S")==0) && (gene==fusA) )
    {
      return fusA_G664S;
    } 
  else if ( (strcmp(sbuf->buff, "S463P")==0) && (gene==rpoB) )
    {
      return rpoB_S463P;
    } 
  else if ( (strcmp(sbuf->buff, "S464P")==0) && (gene==rpoB) )
    {
      return rpoB_S464P;
    } 
  else if ( (strcmp(sbuf->buff, "Q468K")==0) && (gene==rpoB) )
    {
      return rpoB_Q468K;
    } 
  else if ( (strcmp(sbuf->buff, "Q468L")==0) && (gene==rpoB) )
    {
      return rpoB_Q468L;
    } 
  else if ( (strcmp(sbuf->buff, "Q468R")==0) && (gene==rpoB) )
    {
      return rpoB_Q468R;
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
  else if ( (strcmp(sbuf->buff, "A477D")==0) && (gene==rpoB) )
    {
      return rpoB_A477D;
    } 
  else if ( (strcmp(sbuf->buff, "A477V")==0) && (gene==rpoB) )
    {
      return rpoB_A477V;
    } 
  else if ( (strcmp(sbuf->buff, "H481D")==0) && (gene==rpoB) )
    {
      return rpoB_H481D;
    } 
  else if ( (strcmp(sbuf->buff, "H481N")==0) && (gene==rpoB) )
    {
      return rpoB_H481N;
    } 
  else if ( (strcmp(sbuf->buff, "H481Y")==0) && (gene==rpoB) )
    {
      return rpoB_H481Y;
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
  else if ( (strcmp(sbuf->buff, "D550G")==0) && (gene==rpoB) )
    {
      return rpoB_D550G;
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
  else if ( (strcmp(sbuf->buff, "E88K")==0) && (gene==gyrA) )
    {
      return gyrA_E88K;
    } 
  else if ( (strcmp(sbuf->buff, "S80F")==0) && (gene==grlA) )
    {
      return grlA_S80F;
    } 
  else if ( (strcmp(sbuf->buff, "S80Y")==0) && (gene==grlA) )
    {
      return grlA_S80Y;
    }
  else 
    {
      die("Parsing error - unknown mutation %s for gene %d\n", sbuf->buff, (int) gene);
      return NotSpecified;
    } 

}
#endif

#ifdef TB
GeneMutationGene map_gene_name_str_to_genename(StrBuf* name)
            {

    if (strcmp(name->buff, "rpoB")==0)
                    {
                      return rpoB ;
                    }

    else if (strcmp(name->buff, "katG")==0)
    {
      return katG;
    }
    

    else if (strcmp(name->buff, "fabG1")==0)
    {
      return fabG1;
    }
    

    else if (strcmp(name->buff, "pncA")==0)
    {
      return pncA;
    }
    

    else if (strcmp(name->buff, "embB")==0)
    {
      return embB;
    }
    

    else if (strcmp(name->buff, "rrs")==0)
    {
      return rrs;
    }
    

    else if (strcmp(name->buff, "eis")==0)
    {
      return eis;
    }
    

    else if (strcmp(name->buff, "gidB")==0)
    {
      return gidB;
    }
    

    else if (strcmp(name->buff, "rpsL")==0)
    {
      return rpsL;
    }
    

    else if (strcmp(name->buff, "gyrA")==0)
    {
      return gyrA;
    }
    
    else
    {
      die("Failed to parse gene name - got %s", name->buff);
    }

}
            
KnownMutation map_mutation_name_to_enum(StrBuf* sbuf, GeneMutationGene gene)
{

                        if ( (strcmp(sbuf->buff, "F425X")==0) && (gene==rpoB) )
                          {
                            return rpoB_F425X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "G426X")==0) && (gene==rpoB) )
                          {
                            return rpoB_G426X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "T427X")==0) && (gene==rpoB) )
                          {
                            return rpoB_T427X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "S428X")==0) && (gene==rpoB) )
                          {
                            return rpoB_S428X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "Q429X")==0) && (gene==rpoB) )
                          {
                            return rpoB_Q429X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "L430X")==0) && (gene==rpoB) )
                          {
                            return rpoB_L430X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "S431X")==0) && (gene==rpoB) )
                          {
                            return rpoB_S431X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "Q432X")==0) && (gene==rpoB) )
                          {
                            return rpoB_Q432X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "F433X")==0) && (gene==rpoB) )
                          {
                            return rpoB_F433X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "M434X")==0) && (gene==rpoB) )
                          {
                            return rpoB_M434X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "D435X")==0) && (gene==rpoB) )
                          {
                            return rpoB_D435X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "Q436X")==0) && (gene==rpoB) )
                          {
                            return rpoB_Q436X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "N437X")==0) && (gene==rpoB) )
                          {
                            return rpoB_N437X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "N438X")==0) && (gene==rpoB) )
                          {
                            return rpoB_N438X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "P439X")==0) && (gene==rpoB) )
                          {
                            return rpoB_P439X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "L440X")==0) && (gene==rpoB) )
                          {
                            return rpoB_L440X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "S441X")==0) && (gene==rpoB) )
                          {
                            return rpoB_S441X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "G442X")==0) && (gene==rpoB) )
                          {
                            return rpoB_G442X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "L443X")==0) && (gene==rpoB) )
                          {
                            return rpoB_L443X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "T444X")==0) && (gene==rpoB) )
                          {
                            return rpoB_T444X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "H445X")==0) && (gene==rpoB) )
                          {
                            return rpoB_H445X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "K446X")==0) && (gene==rpoB) )
                          {
                            return rpoB_K446X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "R447X")==0) && (gene==rpoB) )
                          {
                            return rpoB_R447X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "R448X")==0) && (gene==rpoB) )
                          {
                            return rpoB_R448X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "L449X")==0) && (gene==rpoB) )
                          {
                            return rpoB_L449X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "S450X")==0) && (gene==rpoB) )
                          {
                            return rpoB_S450X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "A451X")==0) && (gene==rpoB) )
                          {
                            return rpoB_A451X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "L452X")==0) && (gene==rpoB) )
                          {
                            return rpoB_L452X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "A1401X")==0) && (gene==rrs) )
                          {
                            return rrs_A1401X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "C1402X")==0) && (gene==rrs) )
                          {
                            return rrs_C1402X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "G1484X")==0) && (gene==rrs) )
                          {
                            return rrs_G1484X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "S315X")==0) && (gene==katG) )
                          {
                            return katG_S315X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "Tu8X")==0) && (gene==fabG1) )
                          {
                            return fabG1_Tu8X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "Cu15X")==0) && (gene==fabG1) )
                          {
                            return fabG1_Cu15X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "Au16X")==0) && (gene==fabG1) )
                          {
                            return fabG1_Au16X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "Gu17X")==0) && (gene==fabG1) )
                          {
                            return fabG1_Gu17X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "D49N")==0) && (gene==pncA) )
                          {
                            return pncA_D49N;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "D8N")==0) && (gene==pncA) )
                          {
                            return pncA_D8N;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "G162D")==0) && (gene==pncA) )
                          {
                            return pncA_G162D;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "H57D")==0) && (gene==pncA) )
                          {
                            return pncA_H57D;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "H57R")==0) && (gene==pncA) )
                          {
                            return pncA_H57R;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "H71Y")==0) && (gene==pncA) )
                          {
                            return pncA_H71Y;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "L172P")==0) && (gene==pncA) )
                          {
                            return pncA_L172P;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "Q141X")==0) && (gene==pncA) )
                          {
                            return pncA_Q141X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "V125G")==0) && (gene==pncA) )
                          {
                            return pncA_V125G;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "V21G")==0) && (gene==pncA) )
                          {
                            return pncA_V21G;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "M306X")==0) && (gene==embB) )
                          {
                            return embB_M306X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "G406D")==0) && (gene==embB) )
                          {
                            return embB_G406D;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "G406S")==0) && (gene==embB) )
                          {
                            return embB_G406S;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "Cu10T")==0) && (gene==eis) )
                          {
                            return eis_Cu10T;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "A138V")==0) && (gene==gidB) )
                          {
                            return gidB_A138V;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "A200E")==0) && (gene==gidB) )
                          {
                            return gidB_A200E;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "L91P")==0) && (gene==gidB) )
                          {
                            return gidB_L91P;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "P75L")==0) && (gene==gidB) )
                          {
                            return gidB_P75L;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "R137W")==0) && (gene==gidB) )
                          {
                            return gidB_R137W;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "V65G")==0) && (gene==gidB) )
                          {
                            return gidB_V65G;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "V88A")==0) && (gene==gidB) )
                          {
                            return gidB_V88A;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "K43R")==0) && (gene==rpsL) )
                          {
                            return rpsL_K43R;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "K88R")==0) && (gene==rpsL) )
                          {
                            return rpsL_K88R;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "H85X")==0) && (gene==gyrA) )
                          {
                            return gyrA_H85X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "P86X")==0) && (gene==gyrA) )
                          {
                            return gyrA_P86X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "H87X")==0) && (gene==gyrA) )
                          {
                            return gyrA_H87X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "G88X")==0) && (gene==gyrA) )
                          {
                            return gyrA_G88X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "D89X")==0) && (gene==gyrA) )
                          {
                            return gyrA_D89X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "A90X")==0) && (gene==gyrA) )
                          {
                            return gyrA_A90X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "S91X")==0) && (gene==gyrA) )
                          {
                            return gyrA_S91X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "I92X")==0) && (gene==gyrA) )
                          {
                            return gyrA_I92X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "Y93X")==0) && (gene==gyrA) )
                          {
                            return gyrA_Y93X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "D94X")==0) && (gene==gyrA) )
                          {
                            return gyrA_D94X;
                          }
                        
                    else 
                        {
                          die("Parsing error - unknown mutation %s", sbuf->buff);
                          return NotSpecified;
                        } 

}
#endif
