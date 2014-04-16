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
  mut_models.h
*/

#ifndef MUT_MODELS_H_
#define MUT_MODELS_H_


#include "genotyping_known.h"


#define MIN_CONFIDENCE 3

typedef enum
  {
    MaxLikelihood = 0,
    MaxAPosteriori =1,
  } ModelChoiceMethod;

typedef enum
  {
    Resistant = 0,
    Susceptible = 1,
    MixedInfection = 2, //means mixed resistant and susceptible
    Unsure = 3,
  } InfectionType;

double get_log_posterior_truly_resistant_plus_errors_on_suscep_allele(double llk,
								      ResVarInfo* rvi,
								      int max_perc_covg_on_res_allele);

double get_log_posterior_truly_susceptible_plus_errors_on_resistant_allele(double llk,
									   ResVarInfo* rvi,
									   int max_perc_covg_on_res_allele);

double get_log_posterior_of_mixed_infection(double llk,
					    ResVarInfo* rvi,
					    int max_perc_covg_on_res_allele,
					    int perc_covg_susc);


double get_log_lik_truly_resistant_plus_errors_on_suscep_allele(ResVarInfo* rvi,
								double lambda_g, double lambda_e,
								int kmer);

double get_log_lik_truly_susceptible_plus_errors_on_resistant_allele(ResVarInfo* rvi,
								     double lambda_g, double lambda_e,
								     int kmer);

double get_biallelic_log_lik(Covg covg_model_true,
			     Covg  covg_model_err,
			     double lambda_g, 
			     double lambda_e,
			     int kmer);

double get_log_lik_of_mixed_infection(ResVarInfo* rvi,
				      double lambda_g,
				      double err_rate,
				      int kmer);

typedef struct
{
  InfectionType type;
  double likelihood;//log likelihood
  double lp; //log posterior
  double conf;
} Model;

int model_cmp_loglik(const void *a, const void *b);
int model_cmp_logpost(const void *a, const void *b);

void choose_ml_model(double llk_R, double llk_S, double llk_M,
		     Model* best_model);

void choose_map_model(ResVarInfo* rvi,
		      double llk_R, double llk_S, double llk_M,
		      Model* best_model);

InfectionType resistotype(ResVarInfo* rvi, double err_rate, int kmer,
			  double lambda_g, double lambda_e,
			  Model* best_model,
			  ModelChoiceMethod choice);



#endif
