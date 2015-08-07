/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 *  known_mutations.h - from Gordon et al.
*/

#ifndef KNOWN_H_
#define KNOWN_H_

#include "string_buffer.h"

#ifdef STAPH
  typedef enum
  {
    
    rrs=0,
    
    rrs=1,
    
    embB=2,
    
    katG=3,
    
    fabG1=4,
    
    eis=5,
    
    rrs=6,
    
    pncA=7,
    
    gyrA=8,
    
    rpoB=9,
    
    rpsL=10,
    
    rrs=11,
    
    Unknown = 12
  }GeneMutationGene;
  #define NUM_KNOWN_GENES 13
 typedef enum
  {
    
    rrs_A1401X=0,
    
    rrs_C1402X=1,
    
    rrs_G1484X=2,
    
    rrs_A1401X=3,
    
    rrs_C1402X=4,
    
    rrs_G1484X=5,
    
    embB_M306X=6,
    
    embB_G406D=7,
    
    embB_G406S=8,
    
    katG_S315X=9,
    
    fabG1_Tu8X=10,
    
    fabG1_Cu15X=11,
    
    fabG1_Au16X=12,
    
    fabG1_Gu17X=13,
    
    eis_Cu10T=14,
    
    rrs_A1401X=15,
    
    rrs_C1402X=16,
    
    rrs_G1484X=17,
    
    pncA_H57D=18,
    
    gyrA_H85X=19,
    
    gyrA_P86X=20,
    
    gyrA_H87X=21,
    
    gyrA_G88X=22,
    
    gyrA_D89X=23,
    
    gyrA_A90X=24,
    
    gyrA_S91X=25,
    
    gyrA_I92X=26,
    
    gyrA_Y93X=27,
    
    gyrA_D94X=28,
    
    rpoB_F425X=29,
    
    rpoB_G426X=30,
    
    rpoB_T427X=31,
    
    rpoB_S428X=32,
    
    rpoB_Q429X=33,
    
    rpoB_L430X=34,
    
    rpoB_S431X=35,
    
    rpoB_Q432X=36,
    
    rpoB_F433X=37,
    
    rpoB_M434X=38,
    
    rpoB_D435X=39,
    
    rpoB_Q436X=40,
    
    rpoB_N437X=41,
    
    rpoB_N438X=42,
    
    rpoB_P439X=43,
    
    rpoB_L440X=44,
    
    rpoB_S441X=45,
    
    rpoB_G442X=46,
    
    rpoB_L443X=47,
    
    rpoB_T444X=48,
    
    rpoB_H445X=49,
    
    rpoB_K446X=50,
    
    rpoB_R447X=51,
    
    rpoB_R448X=52,
    
    rpoB_S450X=53,
    
    rpoB_A451X=54,
    
    rpoB_L452X=55,
    
    rpsL_K43R=56,
    
    rpsL_K88R=57,
    
    rrs_C513X=58,
    
    rrs_A514X=59,
    
    rrs_G515X=60,
    
    rrs_C516X=61,
    
    rrs_C517X=62,
    
    NotSpecified = 63
  } KnownMutation;
#define NUM_KNOWN_MUTATIONS 64

#endif

#ifdef GN
  typedef enum
  {
    
    rrs=0,
    
    rrs=1,
    
    embB=2,
    
    katG=3,
    
    fabG1=4,
    
    eis=5,
    
    rrs=6,
    
    pncA=7,
    
    gyrA=8,
    
    rpoB=9,
    
    rpsL=10,
    
    rrs=11,
    
    Unknown = 12
  }GeneMutationGene;
  #define NUM_KNOWN_GENES 13
 typedef enum
  {
    
    rrs_A1401X=0,
    
    rrs_C1402X=1,
    
    rrs_G1484X=2,
    
    rrs_A1401X=3,
    
    rrs_C1402X=4,
    
    rrs_G1484X=5,
    
    embB_M306X=6,
    
    embB_G406D=7,
    
    embB_G406S=8,
    
    katG_S315X=9,
    
    fabG1_Tu8X=10,
    
    fabG1_Cu15X=11,
    
    fabG1_Au16X=12,
    
    fabG1_Gu17X=13,
    
    eis_Cu10T=14,
    
    rrs_A1401X=15,
    
    rrs_C1402X=16,
    
    rrs_G1484X=17,
    
    pncA_H57D=18,
    
    gyrA_H85X=19,
    
    gyrA_P86X=20,
    
    gyrA_H87X=21,
    
    gyrA_G88X=22,
    
    gyrA_D89X=23,
    
    gyrA_A90X=24,
    
    gyrA_S91X=25,
    
    gyrA_I92X=26,
    
    gyrA_Y93X=27,
    
    gyrA_D94X=28,
    
    rpoB_F425X=29,
    
    rpoB_G426X=30,
    
    rpoB_T427X=31,
    
    rpoB_S428X=32,
    
    rpoB_Q429X=33,
    
    rpoB_L430X=34,
    
    rpoB_S431X=35,
    
    rpoB_Q432X=36,
    
    rpoB_F433X=37,
    
    rpoB_M434X=38,
    
    rpoB_D435X=39,
    
    rpoB_Q436X=40,
    
    rpoB_N437X=41,
    
    rpoB_N438X=42,
    
    rpoB_P439X=43,
    
    rpoB_L440X=44,
    
    rpoB_S441X=45,
    
    rpoB_G442X=46,
    
    rpoB_L443X=47,
    
    rpoB_T444X=48,
    
    rpoB_H445X=49,
    
    rpoB_K446X=50,
    
    rpoB_R447X=51,
    
    rpoB_R448X=52,
    
    rpoB_S450X=53,
    
    rpoB_A451X=54,
    
    rpoB_L452X=55,
    
    rpsL_K43R=56,
    
    rpsL_K88R=57,
    
    rrs_C513X=58,
    
    rrs_A514X=59,
    
    rrs_G515X=60,
    
    rrs_C516X=61,
    
    rrs_C517X=62,
    
    NotSpecified = 63
  } KnownMutation;
#define NUM_KNOWN_MUTATIONS 64

#endif

#ifdef TB
typedef enum
          {
            Unknown = 0,
            rpoB = 1,
    katG = 2,
    fabG1 = 3,
    pncA = 4,
    embB = 5,
    rrs = 6,
    eis = 7,
    gidB = 8,
    rpsL = 9,
    gyrA = 10
            
          }GeneMutationGene;
    #define NUM_KNOWN_GENES 10
          
typedef enum
    {
    rpoB_F425X = 0,
    rpoB_G426X = 1,
    rpoB_T427X = 2,
    rpoB_S428X = 3,
    rpoB_Q429X = 4,
    rpoB_L430X = 5,
    rpoB_S431X = 6,
    rpoB_Q432X = 7,
    rpoB_F433X = 8,
    rpoB_M434X = 9,
    rpoB_D435X = 10,
    rpoB_Q436X = 11,
    rpoB_N437X = 12,
    rpoB_N438X = 13,
    rpoB_P439X = 14,
    rpoB_L440X = 15,
    rpoB_S441X = 16,
    rpoB_G442X = 17,
    rpoB_L443X = 18,
    rpoB_T444X = 19,
    rpoB_H445X = 20,
    rpoB_K446X = 21,
    rpoB_R447X = 22,
    rpoB_R448X = 23,
    rpoB_S450X = 24,
    rpoB_A451X = 25,
    rpoB_L452X = 26,
    rrs_A1401X = 27,
    rrs_C1402X = 28,
    rrs_G1484X = 29,
    katG_S315X = 30,
    fabG1_Tu8X = 31,
    fabG1_Cu15X = 32,
    fabG1_Au16X = 33,
    fabG1_Gu17X = 34,
    pncA_D49N = 35,
    pncA_D8N = 36,
    pncA_H57D = 37,
    pncA_H57R = 38,
    pncA_H71Y = 39,
    pncA_Q141X = 40,
    pncA_V125G = 41,
    pncA_V21G = 42,
    embB_M306X = 43,
    embB_G406D = 44,
    embB_G406S = 45,
    eis_Cu10T = 46,
    rpsL_K43R = 47,
    rpsL_K88R = 48,
    rrs_C513X = 49,
    rrs_A514X = 50,
    rrs_G515X = 51,
    rrs_C516X = 52,
    rrs_C517X = 53,
    gyrA_H85X = 54,
    gyrA_P86X = 55,
    gyrA_H87X = 56,
    gyrA_G88X = 57,
    gyrA_D89X = 58,
    gyrA_A90X = 59,
    gyrA_S91X = 60,
    gyrA_I92X = 61,
    gyrA_Y93X = 62,
    gyrA_D94X = 63,
    NotSpecified = 64,
  } KnownMutation;
  #define NUM_KNOWN_MUTATIONS 64 

#endif






char* map_var_id_to_drug_resistance(KnownMutation km);
KnownMutation map_mutation_name_to_enum(StrBuf* sbuf, GeneMutationGene gene);
GeneMutationGene map_gene_name_str_to_genename(StrBuf* name);
const char* map_enum_to_mutation_name(KnownMutation km); 

#endif