/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 *
 *  known_mutations.h
*/

#include "known_mutations.h"
#include <string.h>
#include "global.h"

#ifdef STAPH
{% include 'src/predictor/staph/map_gene_name_str_to_genename.c' %}
{% include 'src/predictor/staph/map_mutation_name_to_enum.c' %}
{% include 'src/predictor/staph/map_enum_to_mutation_name.c' %}
{% include 'src/predictor/staph/map_var_id_to_drug_resistance.c' %}
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
const char* map_enum_to_mutation_name(KnownMutation km)
{
   switch (km) 
   {
      case rpoB_F425X : return "rpoB_F425X";
      case rpoB_G426X : return "rpoB_G426X";
      case rpoB_T427X : return "rpoB_T427X";
      case rpoB_S428X : return "rpoB_S428X";
      case rpoB_Q429X : return "rpoB_Q429X";
      case rpoB_L430X : return "rpoB_L430X";
      case rpoB_S431X : return "rpoB_S431X";
      case rpoB_Q432X : return "rpoB_Q432X";
      case rpoB_F433X : return "rpoB_F433X";
      case rpoB_M434X : return "rpoB_M434X";
      case rpoB_D435X : return "rpoB_D435X";
      case rpoB_Q436X : return "rpoB_Q436X";
      case rpoB_N437X : return "rpoB_N437X";
      case rpoB_N438X : return "rpoB_N438X";
      case rpoB_P439X : return "rpoB_P439X";
      case rpoB_L440X : return "rpoB_L440X";
      case rpoB_S441X : return "rpoB_S441X";
      case rpoB_G442X : return "rpoB_G442X";
      case rpoB_L443X : return "rpoB_L443X";
      case rpoB_T444X : return "rpoB_T444X";
      case rpoB_H445X : return "rpoB_H445X";
      case rpoB_K446X : return "rpoB_K446X";
      case rpoB_R447X : return "rpoB_R447X";
      case rpoB_R448X : return "rpoB_R448X";
      case rpoB_S450X : return "rpoB_S450X";
      case rpoB_A451X : return "rpoB_A451X";
      case rpoB_L452X : return "rpoB_L452X";
      case rrs_A1401X : return "rrs_A1401X";
      case rrs_C1402X : return "rrs_C1402X";
      case rrs_G1484X : return "rrs_G1484X";
      case katG_S315X : return "katG_S315X";
      case fabG1_Tu8X : return "fabG1_Tu8X";
      case fabG1_Cu15X : return "fabG1_Cu15X";
      case fabG1_Au16X : return "fabG1_Au16X";
      case fabG1_Gu17X : return "fabG1_Gu17X";
      case pncA_D49N : return "pncA_D49N";
      case pncA_D8N : return "pncA_D8N";
      case pncA_H57D : return "pncA_H57D";
      case pncA_H57R : return "pncA_H57R";
      case pncA_H71Y : return "pncA_H71Y";
      case pncA_Q141X : return "pncA_Q141X";
      case pncA_V125G : return "pncA_V125G";
      case pncA_V21G : return "pncA_V21G";
      case embB_M306X : return "embB_M306X";
      case embB_G406D : return "embB_G406D";
      case embB_G406S : return "embB_G406S";
      case eis_Cu10T : return "eis_Cu10T";
      case rpsL_K43R : return "rpsL_K43R";
      case rpsL_K88R : return "rpsL_K88R";
      case rrs_C513X : return "rrs_C513X";
      case rrs_A514X : return "rrs_A514X";
      case rrs_G515X : return "rrs_G515X";
      case rrs_C516X : return "rrs_C516X";
      case rrs_C517X : return "rrs_C517X";
      case gyrA_H85X : return "gyrA_H85X";
      case gyrA_P86X : return "gyrA_P86X";
      case gyrA_H87X : return "gyrA_H87X";
      case gyrA_G88X : return "gyrA_G88X";
      case gyrA_D89X : return "gyrA_D89X";
      case gyrA_A90X : return "gyrA_A90X";
      case gyrA_S91X : return "gyrA_S91X";
      case gyrA_I92X : return "gyrA_I92X";
      case gyrA_Y93X : return "gyrA_Y93X";
      case gyrA_D94X : return "gyrA_D94X";
      case NotSpecified : return "NotSpecified";
   }
      return  "NotSpecified";

}

char* map_var_id_to_drug_resistance(KnownMutation km)
{
   switch (km) 
   {
      case rpoB_F425X : return "Rifampicin";
      case rpoB_G426X : return "Rifampicin";
      case rpoB_T427X : return "Rifampicin";
      case rpoB_S428X : return "Rifampicin";
      case rpoB_Q429X : return "Rifampicin";
      case rpoB_L430X : return "Rifampicin";
      case rpoB_S431X : return "Rifampicin";
      case rpoB_Q432X : return "Rifampicin";
      case rpoB_F433X : return "Rifampicin";
      case rpoB_M434X : return "Rifampicin";
      case rpoB_D435X : return "Rifampicin";
      case rpoB_Q436X : return "Rifampicin";
      case rpoB_N437X : return "Rifampicin";
      case rpoB_N438X : return "Rifampicin";
      case rpoB_P439X : return "Rifampicin";
      case rpoB_L440X : return "Rifampicin";
      case rpoB_S441X : return "Rifampicin";
      case rpoB_G442X : return "Rifampicin";
      case rpoB_L443X : return "Rifampicin";
      case rpoB_T444X : return "Rifampicin";
      case rpoB_H445X : return "Rifampicin";
      case rpoB_K446X : return "Rifampicin";
      case rpoB_R447X : return "Rifampicin";
      case rpoB_R448X : return "Rifampicin";
      case rpoB_S450X : return "Rifampicin";
      case rpoB_A451X : return "Rifampicin";
      case rpoB_L452X : return "Rifampicin";
      case rrs_A1401X : return "Kanamycin,Capreomycin,Amikacin";
      case rrs_C1402X : return "Kanamycin,Capreomycin,Amikacin";
      case rrs_G1484X : return "Kanamycin,Capreomycin,Amikacin";
      case katG_S315X : return "Isoniazid";
      case fabG1_Tu8X : return "Isoniazid";
      case fabG1_Cu15X : return "Isoniazid";
      case fabG1_Au16X : return "Isoniazid";
      case fabG1_Gu17X : return "Isoniazid";
      case pncA_D49N : return "Pyrazinamide";
      case pncA_D8N : return "Pyrazinamide";
      case pncA_H57D : return "Pyrazinamide";
      case pncA_H57R : return "Pyrazinamide";
      case pncA_H71Y : return "Pyrazinamide";
      case pncA_Q141X : return "Pyrazinamide";
      case pncA_V125G : return "Pyrazinamide";
      case pncA_V21G : return "Pyrazinamide";
      case embB_M306X : return "Ethambutol";
      case embB_G406D : return "Ethambutol";
      case embB_G406S : return "Ethambutol";
      case eis_Cu10T : return "Kanamycin";
      case rpsL_K43R : return "Streptomycin";
      case rpsL_K88R : return "Streptomycin";
      case rrs_C513X : return "Streptomycin";
      case rrs_A514X : return "Streptomycin";
      case rrs_G515X : return "Streptomycin";
      case rrs_C516X : return "Streptomycin";
      case rrs_C517X : return "Streptomycin";
      case gyrA_H85X : return "Quinolones";
      case gyrA_P86X : return "Quinolones";
      case gyrA_H87X : return "Quinolones";
      case gyrA_G88X : return "Quinolones";
      case gyrA_D89X : return "Quinolones";
      case gyrA_A90X : return "Quinolones";
      case gyrA_S91X : return "Quinolones";
      case gyrA_I92X : return "Quinolones";
      case gyrA_Y93X : return "Quinolones";
      case gyrA_D94X : return "Quinolones";
      case NotSpecified : return "NotSpecified";
   }
  return  "NotSpecified";

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
                        
  
                        else if ( (strcmp(sbuf->buff, "K43R")==0) && (gene==rpsL) )
                          {
                            return rpsL_K43R;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "K88R")==0) && (gene==rpsL) )
                          {
                            return rpsL_K88R;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "C513X")==0) && (gene==rrs) )
                          {
                            return rrs_C513X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "A514X")==0) && (gene==rrs) )
                          {
                            return rrs_A514X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "G515X")==0) && (gene==rrs) )
                          {
                            return rrs_G515X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "C516X")==0) && (gene==rrs) )
                          {
                            return rrs_C516X;
                          }
                        
  
                        else if ( (strcmp(sbuf->buff, "C517X")==0) && (gene==rrs) )
                          {
                            return rrs_C517X;
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