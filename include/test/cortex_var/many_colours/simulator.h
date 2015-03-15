/*
 * 
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 *
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
