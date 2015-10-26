/*
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 *
 *  mut_models.c
*/

#include "math.h"
#include "maths.h"
#include "mut_models.h"
#include "string.h"

// *** LOG LIKELIHOODS and POSTERIORS  FOR THE MODELS WITH CLONAL STRAIN + SEQUENCING ERRORS 


//this is for SNPS and indels, not for gene presence
//epsilon =  pow(1--err_rate, cmd_line->kmer_size)
double get_log_posterior_truly_resistant_plus(double llk,
								      int max_perc_covg_on_res_allele,
								      double epsilon)

{
  // prior probability that sample is resistant - look at covg gaps in resistant allele
  if (max_perc_covg_on_res_allele==100)
    {
      return log(1)+llk;
    }
  else
    {
      return -99999999;
    } 
}



double get_log_posterior_truly_susceptible_mut(double llk,
									   int max_perc_covg_on_res_allele,
									   double epsilon)
{
      return log(1)+llk;;
}


double get_log_posterior_of_mixed_infection(double llk,
					    Var* var,
					    int max_perc_covg_on_res_allele,
              int contamination_covg)
{
  if ( (max_perc_covg_on_res_allele==100)
       && 
       (var->vob_best_sus->susceptible_allele.percent_nonzero==100)
       &&
       (contamination_covg == 0)
      )
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




double get_log_lik_R_S_coverage(Var* var,
				       double R_covg, double S_covg,
				       int kmer)
{
  Covg c = get_max_covg_on_any_resistant_allele(var->vob_best_res);
  return get_biallelic_log_lik(c, var->vob_best_sus->susceptible_allele.median_covg,
			       R_covg, S_covg, kmer);
  
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
		      Model* best_model, Model* mid_model, Model* worst_model, double epsilon,
          int contamination_covg)
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

  mR.lp = llk_R + get_log_posterior_truly_resistant_plus(llk_R, 
										 max_perc_covg_on_res,
										 epsilon);
  mS.lp = llk_S + get_log_posterior_truly_susceptible_mut(llk_S, 
										      max_perc_covg_on_res,
										      epsilon);
  mM.lp = llk_M + get_log_posterior_of_mixed_infection(llk_M, var,
						       max_perc_covg_on_res, contamination_covg);
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

   // printf("LP of S, M, R are %f, %f and %f\n", mS.lp , mM.lp, mR.lp );

}


double freq_of_var(Var* var, int expected_covg){
  Covg c = get_max_covg_on_any_resistant_allele(var->vob_best_res);
  double freq = (double) c / expected_covg;
  return freq;
}

InfectionType resistotype(Var* var, 
                          double err_rate,
                          int kmer,
			  double lambda_g, 
			  double lambda_e, 
			  double epsilon,
        int expected_covg,
        int contamination_covg,
			  Model* best_model,
			  ModelChoiceMethod choice,
			  float min_frac_to_detect_minor_pops,
                          boolean* genotyped_present)
{
                          
                          // printf("kmer %d\n", kmer);
                          // printf("err_rate %f\n", err_rate );
                          // printf("lambda_g %f\n", lambda_g );
                          // printf("lambda_e %f\n", lambda_e );
                          // printf(" epsilon %f\n",  epsilon );
                          // printf("min_frac_to_detect_minor_pops %f\n ", min_frac_to_detect_minor_pops );

  if (expected_covg == 0 ){
    //covg must be >0
    expected_covg = 1;
  }
  double llk_M;
  double llk_S;
  double llk_R;


  // If contaminiation is present turn of mixed model and set S to be contam_covg on R with expected_covg on S
  if (contamination_covg > 0 ){
    llk_M = -99999999;
    // Is the resistant coverage due to contamination with S from target
    double llk_S_contaim = get_log_lik_R_S_coverage(var, contamination_covg, expected_covg, kmer); 
    // Is the resistant coverage due to errors with S from target
    double llk_S_error = get_log_lik_R_S_coverage(var, expected_covg * err_rate / 3, expected_covg, kmer);  
    if (llk_S_error > llk_S_contaim){
      llk_S = llk_S_error;
    }else{
      llk_S = llk_S_contaim;
    }
    // Is the resistant coverage due to target with S from errors
    double llk_R_error = get_log_lik_R_S_coverage(var, expected_covg, expected_covg * err_rate / 3, kmer);
    // Is the resistant coverage due to target with S from contamination
    double llk_R_contaim = get_log_lik_R_S_coverage(var, expected_covg, contamination_covg, kmer);
    if (llk_R_error > llk_R_contaim){
      llk_R = llk_R_error;
    }else{
      llk_R = llk_R_contaim;
    }    
  }
  else{
    llk_M = get_log_lik_R_S_coverage(var, min_frac_to_detect_minor_pops * expected_covg, (1-min_frac_to_detect_minor_pops) * expected_covg , kmer );
    llk_S = get_log_lik_R_S_coverage(var, expected_covg * err_rate / 3, expected_covg, kmer);  
    llk_R = get_log_lik_R_S_coverage(var, expected_covg, expected_covg * err_rate / 3, kmer);

  }  




   // printf("LLks of S, M, R are %f, %f and %f\n", llk_S, llk_M, llk_R);

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
		       epsilon, contamination_covg);
    }

    // printf("mid_model.conf %f best_model->conf %f \n", mid_model.conf , best_model->conf);
    // printf("best_model->type %i MIN_CONFIDENCE_r %i \n", best_model->type , MIN_CONFIDENCE_r);

  if (best_model->type != Susceptible){
    *genotyped_present = true;
  }

  if (best_model->type==Susceptible || best_model->type==Resistant)
    {
      return best_model->type;
    }
  //If the model is deciding between R and r return best model 
   else if ( worst_model.type ==Susceptible )
    {
      return best_model->type;
    }    
  // If the model is between r and S
  else if ( worst_model.type==Resistant )
    {
      // If the CONF is below the threshold return inconclusive 
      if (best_model->conf <= MIN_CONFIDENCE_r)
      {
        return Unsure;
      }
      else{
        // If the frequency is below the minimum threshold
        if (freq_of_var(var, expected_covg) < min_frac_to_detect_minor_pops ) {
          return Susceptible;
        }
        else{
          return best_model->type;
        }      
      }
  }
  else
    {
      return Unsure;
    }
}

