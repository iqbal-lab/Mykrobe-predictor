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
    Resistant = 0,
    Susceptible = 1,
    MixedInfection = 2, //means mixed resistant and susceptible
    Unsure = 3,
  } InfectionType;


double get_log_lik_truly_resistant_plus_errors_on_suscep_allele(ResVarInfo* rvi,
								double epsilon,
								double delta,
								double lambda,
								int kmer);

double get_log_lik_truly_susceptible_plus_errors_on_resistant_allele(ResVarInfo* rvi,
								     double epsilon,
								     double delta,
								     double lambda,
								     int kmer);

double get_biallelic_log_lik(Covg covg_model_true,
			     Covg  covg_model_err,
			     double epsilon,
			     double delta,
			     double lambda,
			     int kmer);

double get_log_lik_of_mixed_infection(ResVarInfo* rvi,
				      double epsilon,
				      double lambda,
				      double err_rate,
				      int kmer);

typedef struct
{
  InfectionType type;
  double likelihood;//log likelihood
} Model;

int model_cmp(const void *a, const void *b);

Model choose_best_model(double llk_R, double llk_S, double llk_M,
			double* confidence);

InfectionType best_model(ResVarInfo* rvi, double err_rate, int kmer,
			 int expected_covg, int mean_read_len,
			 double* confidence);

#endif
