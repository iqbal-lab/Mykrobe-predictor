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
  gene_presence_models.h
*/

#include "math.h"
#include "maths.h"
#include "gene_presence_models.h"

//on divergence between the gene panel
//epsilon =  pow(1-err_rate, cmd_line->kmer_size)
double get_log_posterior_truly_resistant(double llk,
					 GeneInfo* gi,
					 double loss_due_to_sample_and_errors,
					 int min_expected)//given known diversity of genes

{

  int p = gi->percent_nonzero;
  if (p>=loss_due_to_sample_and_errors* min_expected)
    {
      return log(1)+llk;
    }
  else
    {
      return -99999999;
    } 
}



double get_log_posterior_truly_susceptible(double llk,
					   GeneInfo* gi,
					   double loss_due_to_sample_and_errors,
					   int min_expected)
{

  int p = gi->percent_nonzero;
  
  if (p>=loss_due_to_sample_and_errors*min_expected)
    {
      return -99999999;
    }
  else
    {
      return log(1)+llk;;
    }
}



// lambda_g = expected_covg/mean_read_len on real allele
double get_log_lik_truly_resistant(GeneInfo* gi,
				   double lambda_g,
				   int kmer)
{

  return get_gene_log_lik(gi->median_covg_on_nonzero_nodes, 
			  lambda_g, 
			  kmer);

}


// epsilon = (1-e)^k
// delta = e(1-e)^(k-1)
// lambda = expected_covg/mean_read_len
double get_log_lik_truly_susceptible(GeneInfo* gi,
				     double lambda_e,
				     int kmer)
{
  
  return get_log_lik_covg_due_to_errors(gi->median_covg_on_nonzero_nodes,
					lambda_e, kmer);
}

// epsilon = (1-e)^k
// lambda = expected_covg/mean_read_len
double get_gene_log_lik(Covg covg,//median covg on parts of gene that are present
			double lambda_g, 
			int kmer)
{
  //Under this model, covg on th true allele is Poisson distributed
  //with rate at true allele   r_t = (1-e)^k  *  depth/read-length =: lambda_g
  double r_t = lambda_g;
  
  // P(covg) = exp(-r_t) * (r_t)^covg /covg!
  
  double log_lik_true_allele  
    = -r_t 
    + covg*log(r_t) 
    - log_factorial(covg);
    
  return log_lik_true_allele;

}

double get_log_lik_covg_due_to_errors(Covg covg,
				      double lambda_e,
				      int kmer)
{
    //Covg on the err allele is Poisson distributed with 
  // rate at error allele = e * (1-e)^(k-1) * (D/R) /3 =: lambda_e
  double r_e = lambda_e;
  
  double log_lik_err_allele  
    = -r_e 
    + covg*log(r_e) 
    - log_factorial(covg);
  
  return log_lik_err_allele;

}



void choose_ml_gene_model(double llk_R, double llk_S,
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

  Model arr[2]={mR, mS};
  qsort(arr, 2, sizeof(Model), model_cmp_loglik);
  best_model->conf = arr[1].likelihood-arr[0].likelihood;
  best_model->type = arr[1].type;
  best_model->likelihood = arr[1].likelihood;
  best_model->lp =0;
}



//max a posteriori
void choose_map_gene_model(GeneInfo* gi,
			   double llk_R, double llk_S, 
			   Model* best_model, double epsilon, int expected_covg,
			   int min_expected_kmer_recovery_for_this_gene)
{

  //total length of gaps = length of gene * exp(-expected_covg)
  // so percentage of kmers present = 1-exp(-expected_covg);
  double loss = 1-exp(-expected_covg);

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

  mR.lp 
    = llk_R 
    + get_log_posterior_truly_resistant(llk_R, 
					gi,
					loss,
					min_expected_kmer_recovery_for_this_gene);

  mS.lp 
    = llk_S 
    + get_log_posterior_truly_susceptible(llk_S, 
					  gi,
					  loss,
					  min_expected_kmer_recovery_for_this_gene);

  Model arr[2]={mR, mS};
  qsort(arr, 2, sizeof(Model), model_cmp_logpost);
  best_model->conf = arr[1].lp-arr[0].lp;
  best_model->type = arr[1].type;
  best_model->likelihood = arr[1].likelihood;
  best_model->lp = arr[1].lp;
}


InfectionType resistotype_gene(GeneInfo* gi, double err_rate, int kmer,
			       double lambda_g,  double epsilon, int expected_covg,
			       Model* best_model,
			       ModelChoiceMethod choice,
			       int min_expected_kmer_recovery_for_this_gene)
{
  double llk_R = get_log_lik_truly_resistant(gi, 
					     lambda_g, 
					     kmer);
  double llk_S = get_log_lik_truly_susceptible(gi, 
					       lambda_g, 
					       kmer);

  best_model->conf=0;
  if (choice==MaxLikelihood)
    {
      choose_ml_gene_model(llk_R, llk_S, best_model);
    }
  else
    {
      choose_map_gene_model(gi, llk_R, llk_S, 
			    best_model, epsilon, expected_covg,
			    min_expected_kmer_recovery_for_this_gene);
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

