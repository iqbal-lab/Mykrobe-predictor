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
  known_mutations.h - from Gordon et al.
*/

#ifndef KNOWN_H_
#define KNOWN_H_

#include "string_buffer.h"

#ifdef STAPH
  typedef enum
  {
    dfrB = 0,
    fusA = 1,
    rpoB = 2,
    gyrA = 3,
    grlA = 4,
    grlB = 5,
    Unknown = 6,
  }GeneMutationGene;
  
  typedef enum
  {
    dfrB_L21V  = 0,
    dfrB_H31N  = 1,
    dfrB_L41F  = 2,
    dfrB_N60I  = 3,
    dfrB_F99Y  = 4,
    dfrB_F99S  = 5,
    dfrB_F99I  = 6,
    dfrB_H150R = 7,
    //start of fusA variants which only cause resistance when in combination
    fusA_F652S = 8,
    fusA_Y654N = 9,
    fusA_L456F = 10,
    fusA_L461F = 11,
    fusA_T326I = 12,
    fusA_A376V = 13,
    fusA_A655P = 14,
    fusA_D463G = 15,
    fusA_E444V = 16,
    fusA_E468V = 17,
    //end of fusA variants which only cause resistance when in combination
    
    //start of fusA variants, each enough on its own to cause resistance
    fusA_V90I  = 18,
    fusA_P114H = 19,
    fusA_Q115L = 20,
    fusA_T385N = 21,
    fusA_P404L = 22,
    fusA_P404Q = 23,
    fusA_P406L = 24,
    fusA_D434N = 25,
    fusA_T436I = 26,
    fusA_H438N = 27,
    fusA_E444K = 28,
    fusA_G451V = 29,
    fusA_G452C = 30,
    fusA_G452S = 31,
    fusA_M453I = 32,
    fusA_H457Q = 33,
    fusA_H457Y = 34,
    fusA_L461K = 35,
    fusA_L461S = 36,
    fusA_R464C = 37,
    fusA_R464S = 38,
    fusA_R464H = 39,
    fusA_P478S = 40,
    fusA_G556S = 41,
    fusA_G617D = 42,
    fusA_M651I = 43,
    fusA_A655E = 44,
    fusA_T656K = 45,
    fusA_R659C = 46,
    fusA_R659H = 47,
    fusA_R659L = 48,
    fusA_R659S = 49,
    fusA_G664S = 50,

    //start of rpoB mutations which must occur together
    rpoB_M470T = 51,
    rpoB_D471G = 52,
    //end of epistatic mutations

    //start of rpoB mutations which individually cause resistance
    rpoB_S463P = 53,
    rpoB_S464P = 54,
    rpoB_Q468K = 55,
    rpoB_Q468L = 56,
    rpoB_Q468R = 57,
    rpoB_D471Y = 58,
    rpoB_N474K = 59,
    rpoB_ins475G = 60,
    rpoB_ins475H = 61,
    rpoB_A477D = 62,
    rpoB_A477V = 63,
    rpoB_H481D = 64,
    rpoB_H481N = 65,
    rpoB_H481Y = 66,
    rpoB_R484H = 67,
    rpoB_S486L = 68,
    rpoB_I527F = 69,
    rpoB_D550G = 70,
    gyrA_S84A = 71,
    gyrA_S84L = 72,
    gyrA_S85P = 73,
    gyrA_E88K = 74,
    grlA_S80F = 75,
    grlA_S80Y = 76,
    NotSpecified = 77,
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
    rpoB_F425X_G426X = 28,
    rpoB_F425X_T427X = 29,
    rpoB_F425X_S428X = 30,
    rpoB_F425X_Q429X = 31,
    rpoB_G426X_T427X = 32,
    rpoB_G426X_S428X = 33,
    rpoB_G426X_Q429X = 34,
    rpoB_G426X_L430X = 35,
    rpoB_T427X_S428X = 36,
    rpoB_T427X_Q429X = 37,
    rpoB_T427X_L430X = 38,
    rpoB_T427X_S431X = 39,
    rpoB_S428X_Q429X = 40,
    rpoB_S428X_L430X = 41,
    rpoB_S428X_S431X = 42,
    rpoB_S428X_Q432X = 43,
    rpoB_Q429X_L430X = 44,
    rpoB_Q429X_S431X = 45,
    rpoB_Q429X_Q432X = 46,
    rpoB_Q429X_F433X = 47,
    rpoB_L430X_S431X = 48,
    rpoB_L430X_Q432X = 49,
    rpoB_L430X_F433X = 50,
    rpoB_L430X_M434X = 51,
    rpoB_S431X_Q432X = 52,
    rpoB_S431X_F433X = 53,
    rpoB_S431X_M434X = 54,
    rpoB_S431X_D435X = 55,
    rpoB_Q432X_F433X = 56,
    rpoB_Q432X_M434X = 57,
    rpoB_Q432X_D435X = 58,
    rpoB_Q432X_Q436X = 59,
    rpoB_F433X_M434X = 60,
    rpoB_F433X_D435X = 61,
    rpoB_F433X_Q436X = 62,
    rpoB_F433X_N437X = 63,
    rpoB_M434X_D435X = 64,
    rpoB_M434X_Q436X = 65,
    rpoB_M434X_N437X = 66,
    rpoB_M434X_N438X = 67,
    rpoB_D435X_Q436X = 68,
    rpoB_D435X_N437X = 69,
    rpoB_D435X_N438X = 70,
    rpoB_D435X_P439X = 71,
    rpoB_Q436X_N437X = 72,
    rpoB_Q436X_N438X = 73,
    rpoB_Q436X_P439X = 74,
    rpoB_Q436X_L440X = 75,
    rpoB_N437X_N438X = 76,
    rpoB_N437X_P439X = 77,
    rpoB_N437X_L440X = 78,
    rpoB_N437X_S441X = 79,
    rpoB_N438X_P439X = 80,
    rpoB_N438X_L440X = 81,
    rpoB_N438X_S441X = 82,
    rpoB_N438X_G442X = 83,
    rpoB_P439X_L440X = 84,
    rpoB_P439X_S441X = 85,
    rpoB_P439X_G442X = 86,
    rpoB_P439X_L443X = 87,
    rpoB_L440X_S441X = 88,
    rpoB_L440X_G442X = 89,
    rpoB_L440X_L443X = 90,
    rpoB_L440X_T444X = 91,
    rpoB_S441X_G442X = 92,
    rpoB_S441X_L443X = 93,
    rpoB_S441X_T444X = 94,
    rpoB_S441X_H445X = 95,
    rpoB_G442X_L443X = 96,
    rpoB_G442X_T444X = 97,
    rpoB_G442X_H445X = 98,
    rpoB_G442X_K446X = 99,
    rpoB_L443X_T444X = 100,
    rpoB_L443X_H445X = 101,
    rpoB_L443X_K446X = 102,
    rpoB_L443X_R447X = 103,
    rpoB_T444X_H445X = 104,
    rpoB_T444X_K446X = 105,
    rpoB_T444X_R447X = 106,
    rpoB_T444X_R448X = 107,
    rpoB_H445X_K446X = 108,
    rpoB_H445X_R447X = 109,
    rpoB_H445X_R448X = 110,
    rpoB_H445X_L449X = 111,
    rpoB_K446X_R447X = 112,
    rpoB_K446X_R448X = 113,
    rpoB_K446X_L449X = 114,
    rpoB_K446X_S450X = 115,
    rpoB_R447X_R448X = 116,
    rpoB_R447X_L449X = 117,
    rpoB_R447X_S450X = 118,
    rpoB_R447X_A451X = 119,
    rpoB_R448X_L449X = 120,
    rpoB_R448X_S450X = 121,
    rpoB_R448X_A451X = 122,
    rpoB_R448X_L452X = 123,
    rpoB_L449X_S450X = 124,
    rpoB_L449X_A451X = 125,
    rpoB_L449X_L452X = 126,
    rpoB_S450X_A451X = 127,
    rpoB_S450X_L452X = 128,
    rpoB_A451X_L452X = 129,
    rrs_A1401X = 130,
    rrs_C1402X = 131,
    rrs_G1484X = 132,
    rrs_A1401X_C1402X = 133,
    katG_S315X = 134,
    fabG1_Tu8X = 135,
    fabG1_Cu15X = 136,
    fabG1_Au16X = 137,
    embB_M306X = 138,
    gyrA_H85X = 139,
    gyrA_P86X = 140,
    gyrA_H87X = 141,
    gyrA_G88X = 142,
    gyrA_D89X = 143,
    gyrA_A90X = 144,
    gyrA_S91X = 145,
    gyrA_I92X = 146,
    gyrA_Y93X = 147,
    gyrA_D94X = 148,
    gyrA_H85X_P86X = 149,
    gyrA_H85X_H87X = 150,
    gyrA_H85X_G88X = 151,
    gyrA_H85X_D89X = 152,
    gyrA_P86X_H87X = 153,
    gyrA_P86X_G88X = 154,
    gyrA_P86X_D89X = 155,
    gyrA_P86X_A90X = 156,
    gyrA_H87X_G88X = 157,
    gyrA_H87X_D89X = 158,
    gyrA_H87X_A90X = 159,
    gyrA_H87X_S91X = 160,
    gyrA_G88X_D89X = 161,
    gyrA_G88X_A90X = 162,
    gyrA_G88X_S91X = 163,
    gyrA_G88X_I92X = 164,
    gyrA_D89X_A90X = 165,
    gyrA_D89X_S91X = 166,
    gyrA_D89X_I92X = 167,
    gyrA_D89X_Y93X = 168,
    gyrA_A90X_S91X = 169,
    gyrA_A90X_I92X = 170,
    gyrA_A90X_Y93X = 171,
    gyrA_A90X_D94X = 172,
    gyrA_S91X_I92X = 173,
    gyrA_S91X_Y93X = 174,
    gyrA_S91X_D94X = 175,
    gyrA_I92X_Y93X = 176,
    gyrA_I92X_D94X = 177,
    gyrA_Y93X_D94X = 178,
    NotSpecified = 179,
  } KnownMutation;
  #define NUM_KNOWN_MUTATIONS 179 

#endif







KnownMutation map_mutation_name_to_enum(StrBuf* sbuf, GeneMutationGene gene);
GeneMutationGene map_gene_name_str_to_genename(StrBuf* name);

#endif
