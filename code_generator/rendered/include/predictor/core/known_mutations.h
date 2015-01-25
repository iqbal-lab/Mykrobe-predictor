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
    
    dfrB=0,
    
    fusA=1,
    
    rpoB=2,
    
    gyrA=3,
    
    grlA=4,
    
    Unknown = 5
  }GeneMutationGene;
  #define NUM_KNOWN_GENES 5
 typedef enum
  {
    
    dfrB_F99I=0,
    
    dfrB_F99S=1,
    
    dfrB_F99Y=2,
    
    dfrB_H150R=3,
    
    dfrB_H31N=4,
    
    dfrB_L21V=5,
    
    dfrB_L41F=6,
    
    dfrB_N60I=7,
    
    fusA_A655P=8,
    
    fusA_A655E=9,
    
    fusA_E444V=10,
    
    fusA_E444K=11,
    
    fusA_F652S=12,
    
    fusA_Y654N=13,
    
    fusA_G451V=14,
    
    fusA_G452C=15,
    
    fusA_G452S=16,
    
    fusA_G556S=17,
    
    fusA_G617D=18,
    
    fusA_G664S=19,
    
    fusA_H438N=20,
    
    fusA_H457Q=21,
    
    fusA_H457Y=22,
    
    fusA_L456F=23,
    
    fusA_L461F=24,
    
    fusA_A376V=25,
    
    fusA_D463G=26,
    
    fusA_L461K=27,
    
    fusA_L461S=28,
    
    fusA_M453I=29,
    
    fusA_M651I=30,
    
    fusA_P114H=31,
    
    fusA_P404L=32,
    
    fusA_P404Q=33,
    
    fusA_P406L=34,
    
    fusA_P478S=35,
    
    fusA_Q115L=36,
    
    fusA_R464C=37,
    
    fusA_R464H=38,
    
    fusA_R464S=39,
    
    fusA_R659C=40,
    
    fusA_R659H=41,
    
    fusA_R659L=42,
    
    fusA_R659S=43,
    
    fusA_T385N=44,
    
    fusA_T436I=45,
    
    fusA_T656K=46,
    
    fusA_V90I=47,
    
    fusA_D434N=48,
    
    fusA_T326I=49,
    
    fusA_E468V=50,
    
    rpoB_A477D=51,
    
    rpoB_A477V=52,
    
    rpoB_D471G=53,
    
    rpoB_D471Y=54,
    
    rpoB_D550G=55,
    
    rpoB_H481D=56,
    
    rpoB_H481N=57,
    
    rpoB_H481Y=58,
    
    rpoB_I527F=59,
    
    rpoB_ins475G=60,
    
    rpoB_ins475H=61,
    
    rpoB_M470T=62,
    
    rpoB_Q468K=63,
    
    rpoB_Q468L=64,
    
    rpoB_Q468R=65,
    
    rpoB_R484H=66,
    
    rpoB_S463P=67,
    
    rpoB_S464P=68,
    
    rpoB_S486L=69,
    
    rpoB_N474K=70,
    
    gyrA_E88K=71,
    
    gyrA_S84A=72,
    
    gyrA_S84L=73,
    
    gyrA_S85P=74,
    
    grlA_S80F=75,
    
    grlA_S80Y=76,
    
    NotSpecified = 77
  } KnownMutation;
#define NUM_KNOWN_MUTATIONS 77

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
    rpoB_L449X = 24,
    rpoB_S450X = 25,
    rpoB_A451X = 26,
    rpoB_L452X = 27,
    rrs_A1401X = 28,
    rrs_C1402X = 29,
    rrs_G1484X = 30,
    katG_S315X = 31,
    fabG1_Tu8X = 32,
    fabG1_Cu15X = 33,
    fabG1_Au16X = 34,
    fabG1_Gu17X = 35,
    pncA_D49N = 36,
    pncA_D8N = 37,
    pncA_H57D = 38,
    pncA_H57R = 39,
    pncA_H71Y = 40,
    pncA_Q141X = 41,
    pncA_V125G = 42,
    pncA_V21G = 43,
    embB_M306X = 44,
    embB_G406D = 45,
    embB_G406S = 46,
    eis_Cu10T = 47,
    rpsL_K43R = 48,
    rpsL_K88R = 49,
    rrs_C513X = 50,
    rrs_A514X = 51,
    rrs_G515X = 52,
    rrs_C516X = 53,
    rrs_C517X = 54,
    gyrA_H85X = 55,
    gyrA_P86X = 56,
    gyrA_H87X = 57,
    gyrA_G88X = 58,
    gyrA_D89X = 59,
    gyrA_A90X = 60,
    gyrA_S91X = 61,
    gyrA_I92X = 62,
    gyrA_Y93X = 63,
    gyrA_D94X = 64,
    NotSpecified = 65,
  } KnownMutation;
  #define NUM_KNOWN_MUTATIONS 65 

#endif






char* map_var_id_to_drug_resistance(KnownMutation km);
KnownMutation map_mutation_name_to_enum(StrBuf* sbuf, GeneMutationGene gene);
GeneMutationGene map_gene_name_str_to_genename(StrBuf* name);
const char* map_enum_to_mutation_name(KnownMutation km); 

#endif