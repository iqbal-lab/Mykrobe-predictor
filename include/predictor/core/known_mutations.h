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
    dfrB_L21V = 0,
    dfrB_H31N = 1,
    dfrB_L41F = 2,
    dfrB_N60I = 3,
    dfrB_F99Y = 4,
    dfrB_F99S = 5,
    dfrB_F99I = 6,
    dfrB_H150R = 7,
    fusA_A67T = 8,
    fusA_A70V = 9,
    fusA_A71V = 10,
    fusA_R76C = 11,
    fusA_V90A = 12,
    fusA_V90I = 13,
    fusA_P114H = 14,
    fusA_Q115L = 15,
    fusA_A160V = 16,
    fusA_M161I = 17,
    fusA_D189G = 18,
    fusA_D189V = 19,
    fusA_E233Q = 20,
    fusA_D373N = 21,
    fusA_A376V = 22,
    fusA_T385N = 23,
    fusA_T387I = 24,
    fusA_P404L = 25,
    fusA_P404Q = 26,
    fusA_P406L = 27,
    fusA_S416F = 28,
    fusA_L430S = 29,
    fusA_D434N = 30,
    fusA_T436I = 31,
    fusA_H438N = 32,
    fusA_F441Y = 33,
    fusA_E444K = 34,
    fusA_E444V = 35,
    fusA_E449K = 36,
    fusA_G451V = 37,
    fusA_G452C = 38,
    fusA_G452S = 39,
    fusA_M453I = 40,
    fusA_L456F = 41,
    fusA_H457Q = 42,
    fusA_H457Y = 43,
    fusA_L461F = 44,
    fusA_L461K = 45,
    fusA_L461S = 46,
    fusA_D463G = 47,
    fusA_R464C = 48,
    fusA_R464S = 49,
    fusA_R464H = 50,
    fusA_C473S = 51,
    fusA_P478S = 52,
    fusA_G556S = 53,
    fusA_H557Y = 54,
    fusA_V607I = 55,
    fusA_G617D = 56,
    fusA_M651I = 57,
    fusA_F652S = 58,
    fusA_Y654N = 59,
    fusA_A655E = 60,
    fusA_A655P = 61,
    fusA_A655V = 62,
    fusA_T656K = 63,
    fusA_R659C = 64,
    fusA_R659H = 65,
    fusA_R659L = 66,
    fusA_R659S = 67,
    fusA_G664S = 68,
    rpoB_N747K = 69,
    rpoB_S463P = 70,
    rpoB_S464P = 71,
    rpoB_L466S = 72,
    rpoB_Q468K = 73,
    rpoB_Q468L = 74,
    rpoB_Q468R = 75,
    rpoB_M470T = 76,
    rpoB_D471G = 77,
    rpoB_D471Y = 78,
    rpoB_A473T = 79,
    rpoB_N474K = 80,
    rpoB_ins475G = 81,
    rpoB_ins475H = 82,
    rpoB_A477D = 83,
    rpoB_A477T = 84,
    rpoB_A477V = 85,
    rpoB_H481D = 86,
    rpoB_H481N = 87,
    rpoB_H481Y = 88,
    rpoB_R484H = 89,
    rpoB_S486L = 90,
    rpoB_I527F = 91,
    rpoB_I527L = 92,
    rpoB_I527M = 93,
    rpoB_S529L = 94,
    rpoB_D550G = 95,
    rpoB_Q565R = 96,
    gyrA_S84A = 97,
    gyrA_S84L = 98,
    gyrA_S85P = 99,
    gyrA_E88V = 100,
    gyrA_E88G = 101,
    gyrA_E88K = 102,
    gyrA_E88L = 103,
    gyrA_G106D = 104,
    grlA_V41G = 105,
    grlA_I45M = 106,
    grlA_A48T = 107,
    grlA_D79V = 108,
    grlA_S80F = 109,
    grlA_S80Y = 110,
    grlA_Y83N = 112,
    grlA_E84G = 112,
    grlA_E84K = 113,
    grlA_E84L = 114,
    grlA_E84V = 115,
    grlA_S108N = 116,
    grlA_A116E = 117,
    grlB_E422D = 118,
    grlA_D432G = 119,
    grlB_D432V = 120,
    grlB_D432NorH = 121,
    grlB_D443E = 122,
    grlB_R444S = 123,
    grlB_P451S = 124,
    grlB_N470D = 125,
    grlB_P585S = 126,
    NotSpecified = 127,
  } KnownMutation;

#define NUM_KNOWN_MUTATIONS 128

KnownMutation map_mutation_name_to_enum(StrBuf* sbuf, GeneMutationGene gene);
GeneMutationGene map_gene_name_str_to_genename(StrBuf* name);

#endif
