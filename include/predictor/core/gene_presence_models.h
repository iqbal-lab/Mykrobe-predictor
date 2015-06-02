/*
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 *  gene_presence_models.c
*/

#include "math.h"
#include "maths.h"
#include "gene_presence.h"
#include "mut_models.h"

#define MIN_GENE_CN 0.03 // Minimum copy number required to call r as r (this is based on agreement with consensus)
#define MIN_GENE_CN_ERY 0.17 // Minimum copy number required to call r as r (this is based on agreement with consensus)
#define MIN_GENE_CN_FUS 0.03 // Minimum copy number required to call r as r (this is based on agreement with consensus)
#define MIN_GENE_CN_GEN 0.04 // Minimum copy number required to call r as r (this is based on agreement with consensus)
#define MIN_GENE_CN_MEC 0.06 // Minimum copy number required to call r as r (this is based on agreement with consensus)
#define MIN_GENE_CN_MUP 0.21 // Minimum copy number required to call r as r (this is based on agreement with consensus)
#define MIN_GENE_CN_PEN 0.04 // Minimum copy number required to call r as r (this is based on agreement with consensus)
#define MIN_GENE_CN_TET 0.12 // Minimum copy number required to call r as r (this is based on agreement with consensus)
//epsilon =  pow(1-err_rate, cmd_line->kmer_size)
double get_log_posterior_major_resistant(double llk,
					 GeneInfo* gi,
					 int expected_covg,					 
					 double err_rate,
					 int min_expected);//given known diversity of genes

double get_log_posterior_minor_resistant(double llk,
					 GeneInfo* gi,
					 int expected_covg,
					 double err_rate,
					 int min_expected);//given known diversity of genes

double get_log_posterior_truly_susceptible(double llk,
					   GeneInfo* gi,
					   int min_expected);

double calculate_expected_gene_coverage_based_on_coverage(int coverage);
// lambda_g = expected_covg/mean_read_len on real allele
/*double get_log_lik_truly_resistant(GeneInfo* gi,
				   double lambda_g,
				   int kmer); */

// lambda = expected_covg/mean_read_len
 double get_log_lik_truly_susceptible(GeneInfo* gi,
				     double lambda_e,
				     int kmer); 

double get_log_lik_observed_coverage_on_gene(GeneInfo* gi,
			     double lambda_g,
			     double freq,//between 0 and 1
			     int expected_covg,
			     int kmer);

// epsilon = (1-e)^k
// lambda = expected_covg/mean_read_len
double get_gene_log_lik(Covg covg,//median covg on parts of gene that are present
			double lambda_g);

double log_poission_prob(double lambda, double k);

double get_log_lik_covg_due_to_errors(Covg covg,
				      int percent_nonzero,
				      int gene_len,
				      double lambda_e,
				      int kmer);


void choose_ml_gene_model(double llk_R, double llk_S, double llk_M,
			  Model* best_model);

//max a posteriori
void choose_map_gene_model(GeneInfo* gi,
			   double llk_R, double llk_S,  double llk_M,
			   Model* best_model, double epsilon, double err_rate, int expected_covg,
			   int min_expected_kmer_recovery_for_this_gene);

InfectionType resistotype_gene(GeneInfo* gi, double err_rate, int kmer,
			       double lambda_g, double lambda_e, double epsilon, int expected_covg,
			       Model* best_model,
			       ModelChoiceMethod choice,
			       int min_expected_kmer_recovery_for_this_gene,
             double min_gene_cn,
             boolean* genotyped_present);

