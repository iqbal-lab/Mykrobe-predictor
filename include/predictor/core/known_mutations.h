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
    fusA_A376V = 12,
    fusA_A655P = 13,
    fusA_D463G = 14,
    fusA_E444V = 15,
    //end of fusA variants which only cause resistance when in combination
    
    //start of fusA variants, each enough on its own to cause resistance
    fusA_V90I  = 16,
    fusA_P114H = 17,
    fusA_Q115L = 18,
    fusA_T326I = 19,
    fusA_T385N = 20,
    fusA_P404L = 21,
    fusA_P404Q = 22,
    fusA_P406L = 23,
    fusA_D434N = 24,
    fusA_T436I = 25,
    fusA_H438N = 26,
    fusA_E444K = 27,
    fusA_G451V = 28,
    fusA_G452C = 29,
    fusA_G452S = 30,
    fusA_M453I = 31,
    fusA_H457Q = 32,
    fusA_H457Y = 33,
    fusA_L461K = 34,
    fusA_L461S = 35,
    fusA_R464C = 36,
    fusA_R464S = 37,
    fusA_R464H = 38,
    fusA_P478S = 39,
    fusA_G556S = 40,
    fusA_G617D = 41,
    fusA_M651I = 42,
    fusA_A655E = 43,
    fusA_T656K = 44,
    fusA_R659C = 45,
    fusA_R659H = 46,
    fusA_R659L = 47,
    fusA_R659S = 48,
    fusA_G664S = 49,

    //start of rpoB mutations which must occur together
    rpoB_M470T = 50,
    rpoB_D471G = 51,
    //end of epistatic mutations

    //start of rpoB mutations which individually cause resistance
    rpoB_S463P = 52,
    rpoB_S464P = 53,
    rpoB_Q468K = 54,
    rpoB_Q468L = 55,
    rpoB_Q468R = 56,
    rpoB_D471Y = 57,
    rpoB_N474K = 58,
    rpoB_ins475G = 59,
    rpoB_ins475H = 60,
    rpoB_A477D = 61,
    rpoB_A477V = 62,
    rpoB_H481D = 63,
    rpoB_H481N = 64,
    rpoB_H481Y = 65,
    rpoB_R484H = 66,
    rpoB_S486L = 67,
    rpoB_I527F = 68,
    rpoB_D550G = 69,
    gyrA_S84A = 70,
    gyrA_S84L = 71,
    gyrA_S85P = 72,
    gyrA_E88K = 73,
    grlA_S80F = 74,
    grlA_S80Y = 75,
    NotSpecified = 76,
  } KnownMutation;

#define NUM_KNOWN_MUTATIONS 76

KnownMutation map_mutation_name_to_enum(StrBuf* sbuf, GeneMutationGene gene);
GeneMutationGene map_gene_name_str_to_genename(StrBuf* name);

#endif
