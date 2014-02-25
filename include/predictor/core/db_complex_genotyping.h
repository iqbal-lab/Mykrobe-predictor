/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
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
  db_complex_genotyping.h
*/

#ifndef DB_COMPLEX_GENOTYPING_H_
#define DB_COMPLEX_GENOTYPING_H_

#include "global.h"
#include "element.h"
#include "genotyping_element.h"
#include "db_variants.h"
#include "graph_info.h"
#include "model_selection.h"
#include "open_hash/little_hash_for_genotyping.h"
#include "genotyping_element.h"

extern int MIN_LLK;

typedef struct {
  Covg* mult11; //multiplicity of allele 1 nodes in allele 1
  Covg* mult12; //multiplicity of allele 1 nodes in allele 2,   etc
  Covg* mult21;
  Covg* mult22;
  int len1;//length of array of nodes for allele1 - ie mult11 and mult12
  int len2;

} MultiplicitiesAndOverlapsOfBiallelicVariant;


typedef struct {
  int max_allele_len;//this applies to working_blah, which have to be as long as the longest allele.
  GenotypingElement** working_g_e_one;
  Orientation* working_o_one;
  GenotypingElement** working_g_e_other;
  Orientation* working_o_other;
  Covg* working_array_self;
  Covg* working_array_shared;
  int working_colour1;
  int working_colour2;
  dBNode** path_nodes;
  Orientation* path_orientations;
  Nucleotide* path_labels;
  char* path_string;
  MultiplicitiesAndOverlapsOfBiallelicVariant* mobv;
  int max_sup_len;//used for the path_blah allocation
} GenotypingWorkingPackage;

GenotypingWorkingPackage* alloc_genotyping_work_package(int max_allele_len,
                                                        int max_sup_len,
                                                        int working_colour1,
                                                        int working_colour2);

void free_genotyping_work_package(GenotypingWorkingPackage* gwp);

MultiplicitiesAndOverlapsOfBiallelicVariant* alloc_MultiplicitiesAndOverlapsOfBiallelicVariant(
  int len_allele1, int len_allele2);

void dealloc_MultiplicitiesAndOverlapsOfBiallelicVariant(MultiplicitiesAndOverlapsOfBiallelicVariant* mobv);
void reset_MultiplicitiesAndOverlapsOfBiallelicVariant(MultiplicitiesAndOverlapsOfBiallelicVariant* mobv);

void initialise_multiplicities_of_allele_nodes_wrt_both_alleles(
  VariantBranchesAndFlanks* var,
  MultiplicitiesAndOverlapsOfBiallelicVariant* mult,
  boolean only_count_nodes_with_edge_in_specified_colour_func,
  Edges (*get_colour)(const dBNode*));
  //Covg (*get_covg)(const dBNode*) // unused

void improved_initialise_multiplicities_of_allele_nodes_wrt_both_alleles(
  VariantBranchesAndFlanks* var,
  MultiplicitiesAndOverlapsOfBiallelicVariant* mult,
  int working_colour1, int working_colour2);
  // unused paramaters
  //boolean only_count_nodes_with_edge_in_specified_colour_func,
  //Edges (*get_colour)(const dBNode*),
  //Covg (*get_covg)(const dBNode*),
  


char** alloc_array_and_get_files_from_list(char* filelist, int num_files_in_list);
void dealloc_array_of_files(char** array_files, int num_files_in_list);

double calc_log_likelihood_of_genotype_with_complex_alleles(
  VariantBranchesAndFlanks* var,
  char* name_of_this_genotype,
  MultiplicitiesAndOverlapsOfBiallelicVariant* var_mults,
  GraphAndModelInfo* model_info,
  int colour_indiv, 
  int colour_ref_minus_our_site, dBGraph* db_graph,
  Covg* working_array_self, Covg* working_array_shared,
  double* current_max_lik, double* current_max_but_one_lik,
  char* current_max_lik_name, char* current_max_but_one_lik_name,
  AssumptionsOnGraphCleaning assump,
  dBNode** p_nodes, Orientation* p_orientations, Nucleotide* p_labels,
  char* p_string, int max_allele_length,
  boolean using_1net, Covg (*get_covg_in_1net_of_genotype)(dBNode*), 
  boolean using_2net, Covg (*get_covg_in_2net_of_genotype)(dBNode*),
  double min_acceptable_llk);


void wipe_colour_and_load_binaries(dBGraph* db_graph, int colour, char* bin1, char* bin2);
void wipe_two_colours_and_load_two_binaries(dBGraph* db_graph,
                                            int colour1, int colour2,
                                            char* binary11, char* binary12,
                                            char* binary21, char* binary22);

//we ASSUME colours 0 to number_alleles are the various alternate alleles, loading in multicolour_bin
void calculate_max_and_max_but_one_llks_of_specified_set_of_genotypes_of_complex_site(
  int* colours_to_genotype, int num_colours_to_genotype,
  int colour_ref_minus_site, int number_alleles,
  int first_gt, int last_gt, // of all the possible gt's
  int max_allele_length, char* fasta,//one read per allele
  AssumptionsOnGraphCleaning assump,
  double* current_max_lik_array, double* current_max_but_one_lik_array,
  char** name_current_max_lik_array, char** name_current_max_but_one_lik_array,
  boolean print_all_liks_calculated,//not just the top two
  GraphAndModelInfo* model_info, dBGraph* db_graph,
  int working_colour1, int working_colour2,
  boolean using_1net, boolean using_2net,
  double min_acceptable_llk);


//this assumes we know which of the tewo alleles is the ref allele - needed for this genotyper
void calculate_llks_for_biallelic_site_using_full_model_for_one_colour_with_known_ref_allele(
  AnnotatedPutativeVariant* annovar,
  AssumptionsOnGraphCleaning assump,
  GraphAndModelInfo* model_info, 
  LittleHashTable* little_db_graph,
  dBGraph* db_graph, 
  GenotypingElement** working_g_e_one, 
  Orientation* working_o_one,
  GenotypingElement** working_g_e_other, 
  Orientation* working_o_other,
  Covg* working_array_self,
  Covg* working_array_shared,
  int len_working_arrays,
  MultiplicitiesAndOverlapsOfBiallelicVariant* mobv,
  int colour_to_genotype, int working_colour1, int working_colour2,
  dBNode** path_nodes, Orientation* path_orientations, 
  Nucleotide* path_labels, char* path_string , int max_sup_len);


double calc_log_prob_of_covg_on_chunk(double eff_D_over_R,
                                      Covg* working_array, int working_array_len);

double get_log_probability_of_covg_on_one_allele_given_second_allele_and_multiplicities(
  double hap_D_over_R, dBNode** allele, int len_allele,
  Covg* mult_this_allele_in_self, Covg* mult_this_allele_in_other,
  Covg* working_array_self, Covg* working_array_shared,
  Covg (*check_covg_in_ref_with_site_excised)(dBNode*),
  int colour_indiv);

double get_log_probability_of_covg_on_one_allele_given_second_allele_and_multiplicities_using_little_hash(
  double hap_D_over_R, GenotypingElement** allele, int len_allele,
  Covg* mult_this_allele_in_self, Covg* mult_this_allele_in_other,
  Covg* working_array_self, Covg* working_array_shared,
  Covg (*check_covg_in_ref_with_site_excised)(GenotypingElement*),
  int colour_indiv);

double calc_log_likelihood_of_genotype_with_complex_alleles_using_little_hash(
  GenotypingVariantBranchesAndFlanks* var,
  MultiplicitiesAndOverlapsOfBiallelicVariant* var_mults,
  GraphAndModelInfo* model_info,
  int colour_indiv, 
  LittleHashTable* little_db_graph, dBGraph* db_graph,
  Covg* working_array_self, Covg* working_array_shared,
  AssumptionsOnGraphCleaning assump,
  dBNode** p_nodes, Orientation* p_orientations,
  Nucleotide* p_labels, char* p_string, 
  int max_sup_len);

void wipe_little_graph(LittleHashTable* little_graph);

void get_all_full_model_genotype_log_likelihoods_at_PD_call_for_one_colour(
  AnnotatedPutativeVariant* annovar,
  AssumptionsOnGraphCleaning assump,
  GraphAndModelInfo* model_info,
  LittleHashTable* little_db_graph, 
  dBGraph* db_graph,
  GenotypingWorkingPackage* gwp, int colour_to_genotype);


double* alloc_ML_results_array(int num_samples_to_genotype);
char** alloc_ML_results_names_array(int num_samples_to_genotype);

void modify_character(char* str, int which_base, int which_mutant);

//set up of Putative Variant object and genotype
boolean initialise_putative_variant(AnnotatedPutativeVariant* annovar,
                                    GraphAndModelInfo* model_info,
                                    VariantBranchesAndFlanks* var,
                                    DiscoveryMethod caller,  
				    int kmer,
                                    AssumptionsOnGraphCleaning assump,
				    CovgArray* working_ca,
                                    boolean do_genotyping, 
				    boolean use_median);

#endif /* DB_COMPLEX_GENOTYPING_H_ */
