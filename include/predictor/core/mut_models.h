/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 *  mut_models.h
*/

#ifndef MUT_MODELS_H_
#define MUT_MODELS_H_


#include "genotyping_known.h"

#ifdef STAPH
  #define MIN_CONFIDENCE 1
#endif
#ifdef TB
	#define MIN_CONFIDENCE 30
#endif  


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
								      int max_perc_covg_on_res_allele,
								      double epsilon);

double get_log_posterior_truly_susceptible_plus_errors_on_resistant_allele(double llk,
									   int max_perc_covg_on_res_allele,
									   double epsilon);

double get_log_posterior_of_mixed_infection(double llk,
					    Var* var,
					    int max_perc_covg_on_res_allele);


double get_log_lik_truly_resistant_plus_errors_on_suscep_allele(Var* var,
								double lambda_g, double lambda_e,
								int kmer);

double get_log_lik_minor_pop_resistant(Var* var,
				       double lambda_g, double lambda_e,
				       int kmer, double err_rate,
				       float min_frac_to_detect_minor_pops);

double get_log_lik_truly_susceptible_plus_errors_on_resistant_allele(Var* var,
								     double lambda_g, double lambda_e,
								     int kmer);

double get_biallelic_log_lik(Covg covg_model_true,
			     Covg  covg_model_err,
			     double lambda_g, 
			     double lambda_e,
			     int kmer);

double get_log_lik_of_mixed_infection(Var* var,
				      double lambda_g,
				      double err_rate,
				      int kmer);

typedef struct
{
  InfectionType type;
  double likelihood;//log likelihood
  double lp; //log posterior
  double conf;//in sorted array of probs - best - 2nd
} Model;

int model_cmp_loglik(const void *a, const void *b);
int model_cmp_logpost(const void *a, const void *b);

void choose_ml_model(double llk_R, double llk_S, double llk_M,
		     Model* best_model);

void choose_map_model(Var* var,
		      double llk_R, double llk_S, double llk_M,
		      Model* best_model, Model* mid_model, Model* worst_model, double epsilon);

InfectionType resistotype(Var* var, 
						  double err_rate, 
						  int kmer,
						  double lambda_g, 
						  double lambda_e, 
						  double epsilon,
						  Model* best_model,
						  ModelChoiceMethod choice,
						  float min_frac_to_detect_minor_pops);



#endif
