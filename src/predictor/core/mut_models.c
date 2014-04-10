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

// *** LOG LIKELIHOODS FOR THE MODELS WITH CLONAL STRAIN + SEQUENCING ERRORS 

// epsilon = (1-e)^k
// delta = e(1-e)^(k-1)
// lambda = expected_covg/mean_read_len
double get_log_lik_truly_resistant_plus_errors_on_suscep_allele(ResVarInfo* rvi,
								double epsilon,
								double delta,
								double lambda,
								int kmer)
{
  Covg c = get_max_covg_on_any_resistant_allele(rvi);
  return get_biallelic_log_lik(c, rvi->susceptible_allele.median_covg,
			       epsilon, delta, lambda, kmer);
}


// epsilon = (1-e)^k
// delta = e(1-e)^(k-1)
// lambda = expected_covg/mean_read_len
double get_log_lik_truly_susceptible_plus_errors_on_resistant_allele(ResVarInfo* rvi,
								     double epsilon,
								     double delta,
								     double lambda,
								     int kmer)
{
  Covg c = get_max_covg_on_any_resistant_allele(rvi);
  return get_biallelic_log_lik(rvi->susceptible_allele.median_covg, c, 
			       epsilon, delta, lambda, kmer);
}

//under a model where the first allele is the true one and second is seq error
// epsilon = (1-e)^k
// delta = e(1-e)^(k-1)
// lambda = expected_covg/mean_read_len
double get_biallelic_log_lik(Covg covg_model_true,//covg on allele the model says is true
			     Covg covg_model_err,
			     double epsilon,
			     double delta,
			     double lambda,
			     int kmer)
{
  //Under this model, covg on th true allele is Poisson distributed
  //with rate at true allele   r_t = (1-e)^k  *  lambda
  double r_t = epsilon * lambda;
  
  // P(covg_model_true) = exp(-r_t) * (r_t)^covg_model_true /covg_model_true!
  
  double log_lik_true_allele  
    = -r_t 
    + covg_model_true*log(r_t) 
    - log_factorial(covg_model_true);
    
  //Covg on the err allele is Poisson distributed with 
  // rate at error allele = e * (1-e)^(k-1) * lambda/3
  double r_e = delta * lambda/3;
  
  double log_lik_err_allele  
    = -r_e 
    + covg_model_err*log(r_e) 
    - log_factorial(covg_model_err);
  
  return log_lik_true_allele+log_lik_err_allele;

}




// ***** MIXED INFECTIONS *********

//llk for model where frequency of resistant allele is the point
//estimate given by allele balance from covgs.
double get_log_lik_of_mixed_infection(ResVarInfo* rvi,
				      double epsilon,
				      double lambda,
				      double err_rate,
				      int kmer)
{
  Covg r = get_max_covg_on_any_resistant_allele(rvi);
  Covg s = rvi->susceptible_allele.median_covg;
  if (r+s==0)
    {
      return -99999999;
    }
  double estim_freq_res = (double) r/(double) (r+s);

  if (estim_freq_res<= 3*err_rate)
    {
      return -999999999; //not allowing this model unless looks likes both res and susc freq > 3* error rate
    }
  else if (estim_freq_res>= 1-3*err_rate)
    {
      return -999999999; 
    }
  
  //Under this model, covg on the susceptible allele is Poisson distributed
  //with rate at susc allele   r_s = (1-e)^k  *  lambda * s/(s+r)
  double r_s = epsilon * lambda * (1-estim_freq_res);
  
  // P(covg_model_susc) = exp(-r_s) * (r_s)^s /s!
  
  double log_lik_susc_allele  
    = -r_s 
    + s*log(r_s) 
    - log_factorial(s);
    
  //Covg on the resistant allele is Poisson distributed with 
  // rate at resistant allele = (1-e)^k * lambda * r/(r+s)
  double r_r = epsilon * lambda * estim_freq_res;
  
  double log_lik_res_allele  
    = -r_r 
    + r*log(r_r) 
    - log_factorial(r);
  
  return log_lik_susc_allele+log_lik_res_allele;
  
}

int model_cmp(const void *a, const void *b)
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

Model choose_best_model(double llk_R, double llk_S, double llk_M,
			double* confidence)
{
  Model mR;
  mR.type=Resistant;
  mR.likelihood=llk_R;
  Model mS;
  mS.type=Susceptible;
  mS.likelihood=llk_S;
  Model mM;
  mM.type=MixedInfection;
  mM.likelihood=llk_M;

  Model arr[3]={mR, mS, mM};
  qsort(arr, 3, sizeof(Model), model_cmp);
  *confidence = arr[2].likelihood-arr[1].likelihood;
  return arr[2];
}

InfectionType best_model(ResVarInfo* rvi, double err_rate, int kmer,
			 int expected_covg, int mean_read_len,
			 double* confidence)
{
  double epsilon = pow(1-err_rate, kmer);
  double delta = err_rate * pow(1-err_rate, kmer-1);
  double lambda = (double) expected_covg/(double) mean_read_len;

  double llk_R = get_log_lik_truly_resistant_plus_errors_on_suscep_allele(rvi, 
									  epsilon, delta, lambda,
									  kmer);
  double llk_S = get_log_lik_truly_susceptible_plus_errors_on_resistant_allele(rvi, 
									       epsilon, delta, lambda,
									       kmer);
  double llk_M = get_log_lik_of_mixed_infection(rvi, epsilon, lambda, err_rate, kmer);

  *confidence=0;
  Model best = choose_best_model(llk_R, llk_S, llk_M, confidence);
  if (*confidence>MIN_CONFIDENCE)
    {
      return best.type;
    }
  else
    {
      return Unsure;
    }
}

