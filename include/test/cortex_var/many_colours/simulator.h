/*
 * 
 * CORTEX project contacts:  
 *    M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 *    Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */
/*
  simulator.h
*/
 
#ifndef SIMULATOR_H_
#define SIMULATOR_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "global.h"
#include "element.h"
#include "db_variants.h"
#include "graph_info.h"
#include "model_selection.h"
#include "db_complex_genotyping.h"


void update_allele(dBNode** allele, int len, int colour, int covg, int read_len);
void zero_allele(dBNode** allele, int len, int colour_indiv, int colour_allele1, int colour_allele2, int colour_ref_minus_site);

//if you want to simulate a true hom, pass in var with both alleles the same
void simulator(int depth, int read_len, int kmer, double seq_err_per_base,
               int number_repetitions,  int colour_indiv,
               int colour_allele1, int colour_allele2, int colour_ref_minus_site,
               VariantBranchesAndFlanks* var,
               int len_genome_minus_site, zygosity true_gt,
               GraphAndModelInfo* model_info,
               char* fasta, char* true_ml_gt_name,
               int working_colour1, int working_colour2,
               boolean using_1and2_nets,
               dBGraph* db_graph);
               //dBNode** genome_minus_site
               //boolean are_the_two_alleles_identical
               //char* filelist_net1, char* filelist_net2
               //int working_colour_1net, int working_colour_2net

#endif /* SIMULATOR_H_ */
