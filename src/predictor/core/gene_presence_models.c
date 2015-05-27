/*
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 *  gene_presence_models.h
*/

#include "math.h"
#include "maths.h"
#include "gene_presence_models.h"

//on divergence between the gene panel
//epsilon =  pow(1-err_rate, cmd_line->kmer_size)
double get_log_posterior_major_resistant(double llk,
					 GeneInfo* gi,
           int expected_covg,
					 double err_rate,
					 int min_expected)//given known diversity of genes

{
  int p = gi->percent_nonzero;

  double expected_percentage_coverage = calculate_expected_gene_coverage_based_on_coverage(expected_covg);
  double minimum_percentage_coverage_required =  expected_percentage_coverage * min_expected ;
  if ( p >= minimum_percentage_coverage_required) 
    {
      return log(1)+llk;
    }
  else
    {
      return -99999999;
    } 
}


double calculate_minmum_detectable_freq_given_error_rate(double err_rate){
  double freq;
  if (err_rate<0.01)
    {
      //i don't believe the error rate is that low, so fix
      err_rate=0.01;
    }

  if (err_rate<0.02)
    {
      freq=0.05;
    }
  else if (err_rate<0.1)
    {
      freq=0.25;
    }
    else{
      freq = 100;
    }
    return (freq);
}

double calculate_expected_gene_coverage_based_on_coverage(int coverage ){
  // Lander Waterman 
  //total length of gaps = length of gene * exp(-expected_covg)
  // so percentage of kmers present = 1-exp(-expected_covg);  
  return (1-exp(-coverage));
}
double get_log_posterior_minor_resistant(double llk,
					 GeneInfo* gi,
					 int expected_covg,
					 double err_rate,
					 int min_expected)//given known diversity of genes

{
  int p = gi->percent_nonzero;
  double freq = calculate_minmum_detectable_freq_given_error_rate(err_rate);
  // I the error rate is too high never call the minor model
  if (freq > 1){
    return -99999999;
  }
  // If the coverage is much higher than expected than the gene probably exists at multi copy number 

  double expected_percentage_coverage = calculate_expected_gene_coverage_based_on_coverage(expected_covg*freq);
  double minimum_percentage_coverage_required =  expected_percentage_coverage * min_expected ;


  //double step function. Coverage gap as might expect for this low frequency
  if ( (p>=minimum_percentage_coverage_required))
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
					   int min_expected)
{
      return log(1)+llk;;
}



// epsilon = (1-e)^k
// delta = e(1-e)^(k-1)
// lambda = expected_covg/mean_read_len
double get_log_lik_truly_susceptible(GeneInfo* gi,
				     double lambda_e,
				     int kmer)
{
  
  return get_log_lik_covg_due_to_errors(gi->median_covg_on_nonzero_nodes,
					gi->percent_nonzero,
					gi->len,
					lambda_e, kmer);
}


//use number of gaps (to see if contiguous), plus median covg
double get_log_lik_resistant(GeneInfo* gi,
			     double lambda_g,
			     double freq,//between 0 and 1
			     int expected_covg,
			     int kmer)
{
  double log_lk_coverage_on_gene = get_gene_log_lik(gi->median_covg_on_nonzero_nodes, expected_covg * freq, kmer);
  // double log_lk_number_of_gaps_in_gene_coverage = log_prob_gaps(gi, expected_covg, freq);
  // printf("freq : %f coverage lk %f\n", freq , log_lk_coverage_on_gene);
  double ret =  log_lk_coverage_on_gene; // + log_lk_number_of_gaps_in_gene_coverage;
  return ret;

}

// epsilon = (1-e)^k
// lambda = expected_covg/mean_read_len
double get_gene_log_lik(Covg covg,//median covg on parts of gene that are present
			double expected_covg, 
			int kmer)
{
  // What is the likely hood that the coverage seen is explained by a gene existing a this frequency. 
  // We're assuming here that the gene exists a CN = 1

  // We can assume that the expected coverage on the gene is the expected coverage of the sample

  // expected_coverage = NL/G = lambda
  // Pr(covg; expected_coverage)= expected_coverage^(covg) e^{-expected_coverage}}{covg!},
  // log(P(covg)) = -expected_coverage + covg * log(expected_coverage) - log(covg !)
  double log_lik_true_allele  
    = -expected_covg 
    + covg*log(expected_covg) 
    - log_factorial(covg);    
  return log_lik_true_allele;
}

// double log_prob_gaps(GeneInfo* gi, int expected_covg, double freq)
// {
//   return -expected_covg * gi->num_gaps;
// }



double get_log_lik_covg_due_to_errors(Covg covg,
				      int percent_nonzero,
				      int gene_len,
				      double lambda_e,
				      int kmer)
{
  //num positive kmers
  int p = (int) (percent_nonzero*gene_len/100);

  //rescale - number of SNP errors
  if (p>kmer)
    {
      p = p/kmer;
    }
  else
    {
      p =1;
    }

    //Covg on the err allele is Poisson distributed with 
  // rate at error allele = e * (1-e)^(k-1) * (D/R) /3 =: lambda_e
  double r_e = lambda_e;
  
  double log_lik_err_allele  
    = -r_e 
    + p*covg*log(r_e) 
    - log_factorial(covg*p);
  
  return log_lik_err_allele;

}



void choose_ml_gene_model(double llk_R, double llk_S, double llk_M,
			  Model* best_model)
{
  Model mR;
  mR.type=Resistant;
  mR.likelihood=llk_R;
  mR.lp = 0;
  mR.conf=0;
  Model mM;
  mM.type=MixedInfection;
  mM.likelihood=llk_M;
  mM.lp = 0;
  mM.conf=0;
  Model mS;
  mS.type=Susceptible;
  mS.likelihood=llk_S;
  mS.lp = 0;
  mS.conf=0;

  Model arr[3]={mR, mS, mM};
  qsort(arr, 3, sizeof(Model), model_cmp_loglik);
  best_model->conf = arr[2].likelihood-arr[1].likelihood;
  best_model->type = arr[2].type;
  best_model->likelihood = arr[2].likelihood;
  best_model->lp =0;
}



//max a posteriori
void choose_map_gene_model(GeneInfo* gi,
			   double llk_R, double llk_S,  double llk_M,
			   Model* best_model, double epsilon, 
			   double err_rate, int expected_covg,
			   int min_expected_kmer_recovery_for_this_gene)
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

  mR.lp 
    = get_log_posterior_major_resistant(llk_R, 
					gi,
          expected_covg,
					err_rate,
					min_expected_kmer_recovery_for_this_gene);

  mM.lp 
    = get_log_posterior_minor_resistant(llk_M, 
					gi,
					expected_covg,
					err_rate,
					min_expected_kmer_recovery_for_this_gene);


  mS.lp 
    = get_log_posterior_truly_susceptible(llk_S, 
					  gi,
					  min_expected_kmer_recovery_for_this_gene);

  Model arr[3]={mR, mS, mM};
  qsort(arr, 3, sizeof(Model), model_cmp_logpost);
  best_model->conf = arr[2].lp-arr[1].lp;
  best_model->type = arr[2].type;
  best_model->likelihood = arr[2].likelihood;
  best_model->lp = arr[2].lp;

   // printf("LP of S, M, R are %f, %f and %f\n", mS.lp , mM.lp, mR.lp );
}


InfectionType resistotype_gene(GeneInfo* gi, double err_rate, int kmer,
			       double lambda_g,  double lambda_e, double epsilon, int expected_covg,
			       Model* best_model,
			       ModelChoiceMethod choice,
			       int min_expected_kmer_recovery_for_this_gene)
{
  //depending on err rate, set freq
  double freq = calculate_minmum_detectable_freq_given_error_rate(err_rate);

  double llk_R = get_log_lik_resistant(gi, lambda_g, 1, expected_covg, kmer);
  double llk_M = get_log_lik_resistant(gi, lambda_g, freq, expected_covg, kmer);
  double llk_S = get_log_lik_truly_susceptible(gi, 
					       lambda_e, 
					       kmer);


   // printf("LLks of S, M, R are %f, %f and %f\n", llk_S, llk_M, llk_R);
  best_model->conf=0;
  if (choice==MaxLikelihood)
    {
      choose_ml_gene_model(llk_R, llk_S, llk_M, best_model);
    }
  else
    {
      choose_map_gene_model(gi, llk_R, llk_S, llk_M,
			    best_model, epsilon, err_rate, expected_covg,
			    min_expected_kmer_recovery_for_this_gene);
    }

  if (best_model->conf > MIN_CONFIDENCE_GENE)
    {
      return best_model->type;
    }
  else
    {
      return Unsure;
    }
}

