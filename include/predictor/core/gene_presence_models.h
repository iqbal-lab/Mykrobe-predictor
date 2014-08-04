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
  gene_presence_models.c
*/

#include "math.h"
#include "maths.h"
#include "gene_presence.h"
#include "mut_models.h"

//epsilon =  pow(1-err_rate, cmd_line->kmer_size)
double get_log_posterior_major_resistant(double llk,
					 GeneInfo* gi,
					 double recovery_given_sample_and_errors,
					 double err_rate,
					 int min_expected);//given known diversity of genes

double get_log_posterior_minor_resistant(double llk,
					 GeneInfo* gi,
					 double recovery_given_sample_and_errors,
					 double err_rate,
					 int min_expected);//given known diversity of genes

double get_log_posterior_truly_susceptible(double llk,
					   GeneInfo* gi,
					   double epsilon,
					   int min_expected);


// lambda_g = expected_covg/mean_read_len on real allele
/*double get_log_lik_truly_resistant(GeneInfo* gi,
				   double lambda_g,
				   int kmer); */

// lambda = expected_covg/mean_read_len
 double get_log_lik_truly_susceptible(GeneInfo* gi,
				     double lambda_e,
				     int kmer); 

// epsilon = (1-e)^k
// lambda = expected_covg/mean_read_len
double get_gene_log_lik(Covg covg,//median covg on parts of gene that are present
			double lambda_g, 
			int kmer);

double get_log_lik_covg_due_to_errors(Covg covg,
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
			       double lambda_g, double epsilon, int expected_covg,
			       Model* best_model,
			       ModelChoiceMethod choice,
			       int min_expected_kmer_recovery_for_this_gene);

