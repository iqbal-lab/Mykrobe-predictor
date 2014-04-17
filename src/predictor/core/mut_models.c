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

// *** LOG LIKELIHOODS and POSTERIORS  FOR THE MODELS WITH CLONAL STRAIN + SEQUENCING ERRORS 


//this is for SNPS and indels, not for gene presence, where the prior would depend
//on divergence between the gene panel
//epsilon =  pow(1-err_rate, cmd_line->kmer_size)
double get_log_posterior_truly_resistant_plus_errors_on_suscep_allele(double llk,
								      ResVarInfo* rvi,
								      int max_perc_covg_on_res_allele,
								      double epsilon)

{
  // prior probability that sample is resistant - look at covg gaps in resistant allele
  int p = max_perc_covg_on_res_allele;
  
  //prob = 0.9 if p==100
  //     = 0.1 if p>80
  //     = 0 else
  if (p>=100*epsilon)
    {
      return log(0.9)+llk;
    }
  else if (p>=80*epsilon)
    {
      return log(0.1) + llk;
    }
  else
    {
      return -99999999;
    }
}



double get_log_posterior_truly_susceptible_plus_errors_on_resistant_allele(double llk,
									   ResVarInfo* rvi,
									   int max_perc_covg_on_res_allele,
									   double epsilon)
{

  int p = max_perc_covg_on_res_allele;
  
  if (p==0)
    {
      return llk;
    }
  else
    {
      return log(0.5) + llk;
    }
}

double get_log_posterior_of_mixed_infection(double llk,
					    ResVarInfo* rvi,
					    int max_perc_covg_on_res_allele)
{
  if ( (max_perc_covg_on_res_allele==100)
       && 
       (rvi->susceptible_allele.percent_nonzero==100) )
    {
      return llk;
    }
  else if ((max_perc_covg_on_res_allele>80)
	   && 
	   (rvi->susceptible_allele.percent_nonzero > 80) )
    {
      return llk+log(0.5);
    }
  else
    {
      return -9999999;
    }
}





// epsilon = (1-e)^k
// delta = e(1-e)^(k-1)
// lambda = expected_covg/mean_read_len
double get_log_lik_truly_resistant_plus_errors_on_suscep_allele(ResVarInfo* rvi,
								double lambda_g, double lambda_e,
								int kmer)
{
  Covg c = get_max_covg_on_any_resistant_allele(rvi);
  return get_biallelic_log_lik(c, rvi->susceptible_allele.median_covg,
			       lambda_g, lambda_e, kmer);

}


// epsilon = (1-e)^k
// delta = e(1-e)^(k-1)
// lambda = expected_covg/mean_read_len
double get_log_lik_truly_susceptible_plus_errors_on_resistant_allele(ResVarInfo* rvi,
								     double lambda_g,
								     double lambda_e,
								     int kmer)
{
  Covg c = get_max_covg_on_any_resistant_allele(rvi);
  return get_biallelic_log_lik(rvi->susceptible_allele.median_covg, c, 
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

//llk for model where frequency of resistant allele is the point
//estimate given by allele balance from covgs.
double get_log_lik_of_mixed_infection(ResVarInfo* rvi,
				      double lambda_g,
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
  //with rate at susc allele   r_s = (1-e)^k  *  (D/R) * s/(s+r)
  double r_s = lambda_g * (1-estim_freq_res);
  
  // P(covg_model_susc) = exp(-r_s) * (r_s)^s /s!
  
  double log_lik_susc_allele  
    = -r_s 
    + s*log(r_s) 
    - log_factorial(s);
    
  //Covg on the resistant allele is Poisson distributed with 
  // rate at resistant allele = (1-e)^k * lambda * r/(r+s)
  double r_r = lambda_g * estim_freq_res;
  
  double log_lik_res_allele  
    = -r_r 
    + r*log(r_r) 
    - log_factorial(r);
  
  return log_lik_susc_allele+log_lik_res_allele;
  
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

  //actually log likelihood.
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
  Model mS;
  mS.type=Susceptible;
  mS.likelihood=llk_S;
  Model mM;
  mM.type=MixedInfection;
  mM.likelihood=llk_M;

  Model arr[3]={mR, mS, mM};
  qsort(arr, 3, sizeof(Model), model_cmp_loglik);
  best_model->conf = arr[2].likelihood-arr[1].likelihood;
  best_model->type = arr[2].type;
  best_model->likelihood = arr[2].likelihood;
  best_model->lp =0;
}



//max a posteriori
void choose_map_model(ResVarInfo* rvi,
		      double llk_R, double llk_S, double llk_M,
		      Model* best_model, double epsilon)
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

  int max_perc_covg_on_res = get_max_perc_covg_on_any_resistant_allele(rvi);

  mR.lp = llk_R + get_log_posterior_truly_resistant_plus_errors_on_suscep_allele(llk_R, rvi,
										 max_perc_covg_on_res,
										 epsilon);
  mS.lp = llk_S + get_log_posterior_truly_susceptible_plus_errors_on_resistant_allele(llk_S, rvi,
										      max_perc_covg_on_res,
										      epsilon);
  mM.lp = llk_M + get_log_posterior_of_mixed_infection(llk_M, rvi,
						       max_perc_covg_on_res);

  Model arr[3]={mR, mS, mM};
  qsort(arr, 3, sizeof(Model), model_cmp_logpost);
  best_model->conf = arr[2].lp-arr[1].lp;
  best_model->type = arr[2].type;
  best_model->likelihood = arr[2].likelihood;
  best_model->lp = arr[2].lp;
}


InfectionType resistotype(ResVarInfo* rvi, double err_rate, int kmer,
			  double lambda_g, double lambda_e, double epsilon,
			  Model* best_model,
			  ModelChoiceMethod choice)
{
  double llk_R = get_log_lik_truly_resistant_plus_errors_on_suscep_allele(rvi, 
									  lambda_g, lambda_e,
									  kmer);
  double llk_S = get_log_lik_truly_susceptible_plus_errors_on_resistant_allele(rvi, 
									       lambda_g, lambda_e,
									       kmer);
  double llk_M = get_log_lik_of_mixed_infection(rvi, lambda_g, err_rate, kmer);

  best_model->conf=0;
  if (choice==MaxLikelihood)
    {
      choose_ml_model(llk_R, llk_S, llk_M, best_model);
    }
  else
    {
      choose_map_model(rvi, llk_R, llk_S, llk_M, best_model, epsilon);
    }

  if (best_model->conf > MIN_CONFIDENCE)
    {
      return best_model->type;
    }
  else
    {
      return Unsure;
    }
}

