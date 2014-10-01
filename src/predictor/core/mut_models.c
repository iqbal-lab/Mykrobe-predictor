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
  mut_models.c
*/

#include "math.h"
#include "maths.h"
#include "mut_models.h"
#include "string.h"

// *** LOG LIKELIHOODS and POSTERIORS  FOR THE MODELS WITH CLONAL STRAIN + SEQUENCING ERRORS 


//this is for SNPS and indels, not for gene presence
//epsilon =  pow(1--err_rate, cmd_line->kmer_size)
double get_log_posterior_truly_resistant_plus_errors_on_suscep_allele(double llk,
								      int max_perc_covg_on_res_allele,
								      double epsilon)

{
  // prior probability that sample is resistant - look at covg gaps in resistant allele
  int p = max_perc_covg_on_res_allele;


  if (p>=80*epsilon)
    {
      return log(1)+llk;
    }
  else
    {
      return -99999999;
    } 
}



double get_log_posterior_truly_susceptible_plus_errors_on_resistant_allele(double llk,
									   int max_perc_covg_on_res_allele,
									   double epsilon)
{
      return log(1)+llk;;
}


double get_log_posterior_of_mixed_infection(double llk,
					    Var* var,
					    int max_perc_covg_on_res_allele)
{
  if ( (max_perc_covg_on_res_allele==100)
       && 
       (var->vob_best_sus->susceptible_allele.percent_nonzero==100) )
    {
      return llk;
    }
  /*  else if ((max_perc_covg_on_res_allele>80)
	   && 
	   (var->vob_best_sus->susceptible_allele.percent_nonzero > 80) )
    {
      return llk+log(0.5);
      }*/
  else
    {
      return -9999999;
    }
}




// epsilon = (1-e)^k
// delta = e(1-e)^(k-1)
// lambda = expected_covg/mean_read_len
double get_log_lik_truly_resistant_plus_errors_on_suscep_allele(Var* var,
								double lambda_g, double lambda_e,
								int kmer)
{
  Covg c = get_max_covg_on_any_resistant_allele(var->vob_best_res);
  //printf("Test res model, R covg=%d and S covg=%d\n", c, var->susceptible_allele.median_covg);
  return get_biallelic_log_lik(c, var->vob_best_sus->susceptible_allele.median_covg,
			       0.75*lambda_g, lambda_e, kmer);

}



//major population of resistant, assume poisson at 0.75 * expected. ideally dispersed poisson
// epsilon = (1-e)^k
// delta = e(1-e)^(k-1)
// lambda = expected_covg/mean_read_len
double get_log_lik_major_pop_resistant(Var* var,
				       double lambda_g, double lambda_e,
				       int kmer)
{
  Covg c = get_max_covg_on_any_resistant_allele(var->vob_best_res);
  //printf("Test res model, R covg=%d and S covg=%d\n", c, var->susceptible_allele.median_covg);
  return get_biallelic_log_lik(c, var->vob_best_sus->susceptible_allele.median_covg,
			       0.9*lambda_g, 0.1*lambda_g, kmer);

}



double get_log_lik_minor_pop_resistant(Var* var,
				       double lambda_g, double lambda_e,
				       int kmer, double err_rate,
               float min_frac_to_detect_minor_pops)
{
  Covg c = get_max_covg_on_any_resistant_allele(var->vob_best_res);

  double frac = min_frac_to_detect_minor_pops;
  if (err_rate > 0.1)
    {
      return -999999999;
    }
  return get_biallelic_log_lik(c, var->vob_best_sus->susceptible_allele.median_covg,
			       frac*lambda_g, (1-frac)*lambda_g, kmer);

}


// epsilon = (1-e)^k
// delta = e(1-e)^(k-1)
// lambda = expected_covg/mean_read_len
double get_log_lik_truly_susceptible_plus_errors_on_resistant_allele(Var* var,
								     double lambda_g,
								     double lambda_e,
								     int kmer)
{
  Covg c = get_max_covg_on_any_resistant_allele(var->vob_best_res);
  return get_biallelic_log_lik(var->vob_best_sus->susceptible_allele.median_covg, c, 
			       lambda_g, lambda_e, kmer);
}

//under a model where the first allele is the true one and second is seq error
// epsilon = (1-e)^k
// delta = e(1-e)^(k-1)
// lambda = expected_covg/mean_read_len
double get_biallelic_log_lik(Covg covg_model_true,//covg on allele the model says is true
			     Covg covg_model_err,
			     double lambda_g, 
			     double lambda_e,
			     int kmer)
{
  //Under this model, covg on th true allele is Poisson distributed
  //with rate at true allele   r_t = (1-e)^k  *  depth/read-length =: lambda_g
  double r_t = lambda_g;
  
  // P(covg_model_true) = exp(-r_t) * (r_t)^covg_model_true /covg_model_true!
  
  double log_lik_true_allele  
    = -r_t 
    + covg_model_true*log(r_t) 
    - log_factorial(covg_model_true);
    
  //Covg on the err allele is Poisson distributed with 
  // rate at error allele = e * (1-e)^(k-1) * (D/R) /3 =: lambda_e
  double r_e = lambda_e;
  
  double log_lik_err_allele  
    = -r_e 
    + covg_model_err*log(r_e) 
    - log_factorial(covg_model_err);
  
  return log_lik_true_allele+log_lik_err_allele;

}




// ***** MIXED INFECTIONS *********

//minor population of resistant
double get_log_lik_of_mixed_infection(Var* var,
				      double lambda_g,
				      double err_rate,
				      int kmer)
{

  Covg r = get_max_covg_on_any_resistant_allele(var->vob_best_res);
  Covg s = var->vob_best_sus->susceptible_allele.median_covg;

  if ( (r==0) || (s==0) )
    {
      return -9999999;
    } 

  //We look for sub pop with frequency >= 2*err_rate
  //take a flat prior for minor resistant allele freq
  double pr;
  if (err_rate <= 0.01)
    {
      pr=0.1;
    }
  else if (err_rate <= 0.05)
    {
      pr=0.2;
    }
  else if (err_rate <=0.1)
    {
      pr=0.3;
    }
  else
    {
      return -9999999;
    }


  // likelihood = likelihood( | freq of res allele =pr) 
  double r_r = lambda_g * pr;
  
  double log_lik_res_allele  
    = -r_r 
    + r*log(r_r) 
    - log_factorial(r);
  
  return log_lik_res_allele;
  
}

int model_cmp_loglik(const void *a, const void *b)
{
  // casting pointer types
  const Model *ia = (const Model *)a;
  const Model *ib = (const Model *)b;

  //actually log likelihood.
  if (ia->likelihood < ib->likelihood)
    {
      return -1;
    }
  else if (ia->likelihood > ib->likelihood)
    {
      return 1;
    }
  else
    {
      return 0;
    }
}

int model_cmp_logpost(const void *a, const void *b)
{
  // casting pointer types
  const Model *ia = (const Model *)a;
  const Model *ib = (const Model *)b;


  if (ia->lp < ib->lp)
    {
      return -1;
    }
  else if (ia->lp > ib->lp)
    {
      return 1;
    }
  else
    {
      return 0;
    }
}

void choose_ml_model(double llk_R, double llk_S, double llk_M,
		     Model* best_model)
{
  Model mR;
  mR.type=Resistant;
  mR.likelihood=llk_R;
  mR.lp = 0;
  mR.conf=0;

  Model mS;
  mS.type=Susceptible;
  mS.likelihood=llk_S;
  mS.lp = 0;
  mS.conf=0;

  Model mM;
  mM.type=MixedInfection;
  mM.likelihood=llk_M;
  mM.lp = 0;
  mM.conf=0;

  Model arr[3]={mR, mS, mM};
  qsort(arr, 3, sizeof(Model), model_cmp_loglik);
  best_model->conf = arr[2].likelihood-arr[1].likelihood;
  best_model->type = arr[2].type;
  best_model->likelihood = arr[2].likelihood;
  best_model->lp =0;
}



//max a posteriori
void choose_map_model(Var* var,
		      double llk_R, double llk_S, double llk_M,
		      Model* best_model, Model* mid_model, Model* worst_model, double epsilon)
{

  Model mR;
  mR.type=Resistant;
  mR.likelihood=llk_R;
  mR.lp=0;
  mR.conf=0;
  Model mS;
  mS.type=Susceptible;
  mS.likelihood=llk_S;
  mS.lp=0;
  mS.conf=0;
  Model mM;
  mM.type=MixedInfection;
  mM.likelihood=llk_M;
  mM.lp=0;
  mM.conf=0;

  int max_perc_covg_on_res = get_max_perc_covg_on_any_resistant_allele(var->vob_best_res);

  mR.lp = llk_R + get_log_posterior_truly_resistant_plus_errors_on_suscep_allele(llk_R, 
										 max_perc_covg_on_res,
										 epsilon);
  mS.lp = llk_S + get_log_posterior_truly_susceptible_plus_errors_on_resistant_allele(llk_S, 
										      max_perc_covg_on_res,
										      epsilon);
  mM.lp = llk_M + get_log_posterior_of_mixed_infection(llk_M, var,
						       max_perc_covg_on_res);

  Model arr[3]={mR, mS, mM};
  qsort(arr, 3, sizeof(Model), model_cmp_logpost);
  best_model->conf = arr[2].lp-arr[1].lp;
  best_model->type = arr[2].type;
  best_model->likelihood = arr[2].likelihood;
  best_model->lp = arr[2].lp;

  mid_model->conf = arr[1].lp-arr[0].lp;
  mid_model->type = arr[1].type;
  mid_model->likelihood = arr[1].likelihood;
  mid_model->lp = arr[1].lp;

  worst_model->conf = 0;
  worst_model->type = arr[0].type;
  worst_model->likelihood = arr[0].likelihood;
  worst_model->lp = arr[0].lp;


}


InfectionType resistotype(Var* var, 
                          double err_rate, int kmer,
                			  double lambda_g, 
                        double lambda_e, 
                        double epsilon,
                			  Model* best_model,
                			  ModelChoiceMethod choice,
                        float min_frac_to_detect_minor_pops)
{
  double llk_R = get_log_lik_truly_resistant_plus_errors_on_suscep_allele(var, 
									  lambda_g, lambda_e,
									  kmer);
  double llk_S = get_log_lik_truly_susceptible_plus_errors_on_resistant_allele(var, 
									       lambda_g, lambda_e,
									       kmer);
  //  double llk_M = get_log_lik_of_mixed_infection(var, lambda_g, err_rate, kmer);
  double llk_M = get_log_lik_minor_pop_resistant(var,lambda_g, lambda_e, kmer, err_rate,min_frac_to_detect_minor_pops);

  best_model->conf=0;
  Model worst_model;
  worst_model.type=Unsure;
  worst_model.likelihood=0;
  worst_model.lp=0;
  worst_model.conf=-1111111111;
  Model mid_model;
  mid_model.type=Unsure;
  mid_model.likelihood=0;
  mid_model.lp=0;
  mid_model.conf=-1111111111;
  if (choice==MaxLikelihood)
    {
      choose_ml_model(llk_R, llk_S, llk_M, best_model);
    }
  else
    {
      choose_map_model(var, llk_R, llk_S, llk_M, 
		       best_model, &mid_model,&worst_model, 
		       epsilon);
    }


  if (best_model->type==Susceptible)
    {
      return best_model->type;
    }
  else if (best_model->conf > MIN_CONFIDENCE)
  /* if ( ( (best_model->type==Susceptible) || (best_model->type==Resistant) )
       &&
       (best_model->conf > MIN_CONFIDENCE) 
       ) */
    {
      return best_model->type;
    }
  //so now the winning model is M or R. If S is the bottom of the 3 by MIN_CONFIDENCE, call r.
  else if ( (worst_model.type ==Susceptible) && (mid_model.conf>MIN_CONFIDENCE) )
    {
      return mid_model.type;
    }
  else
    {
      return Unsure;
    }
}

