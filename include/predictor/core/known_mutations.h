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
    dfrB_H31N = 0,
    dfrB_L41F = 1,
    dfrB_N60I = 2,
    dfrB_F99Y = 3,
    dfrB_F99S = 4,
    dfrB_F99I = 5,
    fusA_V90I = 6,
    fusA_P114H = 7,
    fusA_Q115L = 8,
    fusA_A376V = 9,
    fusA_T385N = 10,
    fusA_P404L = 11,
    fusA_P404Q = 12,
    fusA_P406L = 13,
    fusA_D434N = 14,
    fusA_T436I = 15,
    fusA_H438N = 16,
    fusA_E444K = 17,
    fusA_E444V = 18,
    fusA_G451V = 19,
    fusA_G452C = 20,
    fusA_G452S = 21,
    fusA_M453I = 22,
    fusA_H457Q = 23,
    fusA_H457Y = 24,
    fusA_L461F = 25,
    fusA_L461K = 26,
    fusA_L461S = 27,
    fusA_D463G = 28,
    fusA_R464C = 29,
    fusA_R464S = 30,
    fusA_R464H = 31,
    fusA_P478S = 32,
    fusA_G556S = 33,
    fusA_G617D = 34,
    fusA_M651I = 35,
    fusA_F652S = 36,
    fusA_Y654N = 37,
    fusA_A655E = 38,
    fusA_A655P = 39,
    fusA_T656K = 40,
    fusA_R659C = 41,
    fusA_R659H = 42,
    fusA_R659L = 43,
    fusA_R659S = 44,
    fusA_G664S = 45,
    rpoB_S463P = 46,
    rpoB_S464P = 47,
    rpoB_Q468K = 48,
    rpoB_Q468L = 49,
    rpoB_Q468R = 50,
    rpoB_M470T = 51,
    rpoB_D471G = 52,
    rpoB_D471Y = 53,
    rpoB_ins475G = 54,
    rpoB_ins475H = 55,
    rpoB_A477D = 56,
    rpoB_A477V = 57,
    rpoB_H481D = 58,
    rpoB_H481N = 59,
    rpoB_H481Y = 60,
    rpoB_R484H = 61,
    rpoB_S486L = 62,
    rpoB_I527F = 63,
    rpoB_D550G = 64,
    rpoB_N747K = 65,
    gyrA_S84A = 66,
    gyrA_S84L = 67,
    gyrA_S85P = 68,
    gyrA_E88K = 69,
    grlA_S80F = 70,
    grlA_S80Y = 71,
    NotSpecified = 72,
  } KnownMutation;

#define NUM_KNOWN_MUTATIONS 72

KnownMutation map_mutation_name_to_enum(StrBuf* sbuf, GeneMutationGene gene);
GeneMutationGene map_gene_name_str_to_genename(StrBuf* name);

#endif
