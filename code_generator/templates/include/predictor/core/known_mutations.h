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
{% include 'include/predictor/common/gene_mutation_gene.h' %}
{% include 'include/predictor/common/known_mutation.h' %}
#endif

#ifdef GN
{% include 'include/predictor/common/gene_mutation_gene.h' %}
{% include 'include/predictor/common/known_mutation.h' %}
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