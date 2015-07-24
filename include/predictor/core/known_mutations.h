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
    
    twentythreeS=0,
    
    dfrB=1,
    
    fusA=2,
    
    rpoB=3,
    
    gyrA=4,
    
    grlA=5,
    
    Unknown = 6
  }GeneMutationGene;
  #define NUM_KNOWN_GENES 6
 typedef enum
  {
    
    twentythreeS_G2576T=0,
    
    dfrB_F99I=1,
    
    dfrB_F99S=2,
    
    dfrB_F99Y=3,
    
    dfrB_H150R=4,
    
    dfrB_H31N=5,
    
    dfrB_L21V=6,
    
    dfrB_L41F=7,
    
    dfrB_N60I=8,
    
    fusA_A655P=9,
    
    fusA_A655E=10,
    
    fusA_E444V=11,
    
    fusA_E444K=12,
    
    fusA_F652S=13,
    
    fusA_Y654N=14,
    
    fusA_G451V=15,
    
    fusA_G452C=16,
    
    fusA_G452S=17,
    
    fusA_G556S=18,
    
    fusA_G617D=19,
    
    fusA_G664S=20,
    
    fusA_H438N=21,
    
    fusA_H457Q=22,
    
    fusA_H457Y=23,
    
    fusA_L456F=24,
    
    fusA_L461F=25,
    
    fusA_A376V=26,
    
    fusA_D463G=27,
    
    fusA_L461K=28,
    
    fusA_L461S=29,
    
    fusA_M453I=30,
    
    fusA_M651I=31,
    
    fusA_P114H=32,
    
    fusA_P404L=33,
    
    fusA_P404Q=34,
    
    fusA_P406L=35,
    
    fusA_P478S=36,
    
    fusA_Q115L=37,
    
    fusA_R464C=38,
    
    fusA_R464H=39,
    
    fusA_R464S=40,
    
    fusA_R659C=41,
    
    fusA_R659H=42,
    
    fusA_R659L=43,
    
    fusA_R659S=44,
    
    fusA_T385N=45,
    
    fusA_T436I=46,
    
    fusA_T656K=47,
    
    fusA_V90I=48,
    
    fusA_D434N=49,
    
    fusA_T326I=50,
    
    fusA_E468V=51,
    
    rpoB_A477D=52,
    
    rpoB_A477V=53,
    
    rpoB_D471G=54,
    
    rpoB_D471Y=55,
    
    rpoB_D550G=56,
    
    rpoB_H481D=57,
    
    rpoB_H481N=58,
    
    rpoB_H481Y=59,
    
    rpoB_I527F=60,
    
    rpoB_ins475G=61,
    
    rpoB_ins475H=62,
    
    rpoB_M470T=63,
    
    rpoB_Q468K=64,
    
    rpoB_Q468L=65,
    
    rpoB_Q468R=66,
    
    rpoB_R484H=67,
    
    rpoB_S463P=68,
    
    rpoB_S464P=69,
    
    rpoB_S486L=70,
    
    rpoB_N474K=71,
    
    gyrA_E88K=72,
    
    gyrA_S84A=73,
    
    gyrA_S84L=74,
    
    gyrA_S85P=75,
    
    grlA_S80F=76,
    
    grlA_S80Y=77,
    
    NotSpecified = 78
  } KnownMutation;
#define NUM_KNOWN_MUTATIONS 78

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
    rpoB_F425X = 0
    rpoB_G426X = 1
    rpoB_T427X = 2
    rpoB_S428X = 3
    rpoB_Q429X = 4
    rpoB_L430X = 5
    rpoB_S431X = 6
    rpoB_Q432X = 7
    rpoB_F433X = 8
    rpoB_M434X = 9
    rpoB_D435X = 10
    rpoB_Q436X = 11
    rpoB_N437X = 12
    rpoB_N438X = 13
    rpoB_P439X = 14
    rpoB_L440X = 15
    rpoB_S441X = 16
    rpoB_G442X = 17
    rpoB_L443X = 18
    rpoB_T444X = 19
    rpoB_H445X = 20
    rpoB_K446X = 21
    rpoB_R447X = 22
    rpoB_R448X = 23
    rpoB_S450X = 24
    rpoB_A451X = 25
    rpoB_L452X = 26
    rrs_A1401X = 27
    rrs_C1402X = 28
    rrs_G1484X = 29
    katG_S315X = 30
    fabG1_Tu8X = 31
    fabG1_Cu15X = 32
    fabG1_Au16X = 33
    fabG1_Gu17X = 34
    pncA_D49N = 35
    pncA_D8N = 36
    pncA_H57D = 37
    pncA_H57R = 38
    pncA_H71Y = 39
    pncA_Q141X = 40
    pncA_V125G = 41
    pncA_V21G = 42
    embB_M306X = 43
    embB_G406D = 44
    embB_G406S = 45
    eis_Cu10T = 46
    rpsL_K43R = 47
    rpsL_K88R = 48
    rrs_C513X = 49
    rrs_A514X = 50
    rrs_G515X = 51
    rrs_C516X = 52
    rrs_C517X = 53
    gyrA_H85X = 54
    gyrA_P86X = 55
    gyrA_H87X = 56
    gyrA_G88X = 57
    gyrA_D89X = 58
    gyrA_A90X = 59
    gyrA_S91X = 60
    gyrA_I92X = 61
    gyrA_Y93X = 62
    gyrA_D94X = 63
    NotSpecified = 64
  } KnownMutation;
  #define NUM_KNOWN_MUTATIONS 64

#endif






char* map_var_id_to_drug_resistance(KnownMutation km);
KnownMutation map_mutation_name_to_enum(StrBuf* sbuf, GeneMutationGene gene);
GeneMutationGene map_gene_name_str_to_genename(StrBuf* name);
const char* map_enum_to_mutation_name(KnownMutation km); 

#endif