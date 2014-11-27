/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 *  species.h
*/
#include "dB_graph.h"

typedef enum 
 {
   Mtuberculosis = 0,
   Mafricanum = 1,
   Mbovis = 2 
  } Myc_species ;


typedef enum 
 {
  beijing = 0,
  animal = 1, 
  lineage1 = 2,
  lineage2 = 3, 
  lineage3 = 4, 
  lineage4 = 5, 
  lineage5 = 6, 
  lineage6 = 7, 
  lineage7 = 8
  } Myc_lineage ;

#define NUM_SPECIES 9
typedef enum
  {
    PureMTBC =0,
    MixedMTB =1,
    // MajorMTBAndMinorNonMTB = 1,
    // MinorMTBAndMajorNonMTB = 2,
    NonMTB = 2,
  } SampleType;

typedef struct
{
  SampleType type;
  double likelihood;//log likelihood
  double lp; //log posterior
  double conf;
  StrBuf* name_of_pure_mtbc_species;//how we interpret this depends on the model type.
  StrBuf* name_of_non_mtb_species;//how we interpret this depends on the model type.
  StrBuf* name_of_pure_mtbc_lineage;//how we interpret this depends on the model type.
  StrBuf* name_of_non_mtb_lineage;//how we interpret this depends on the model type.
} SampleModel;

SampleModel* alloc_and_init_sample_model();
void free_sample_model(SampleModel* sm);

void map_lineage_enum_to_str(Myc_lineage sp, StrBuf* sbuf);
void map_species_enum_to_str(Myc_species sp, StrBuf* sbuf);

Myc_species map_lineage_enum_to_species_enum(Myc_lineage sp);

SampleType get_species_model(dBGraph *db_graph,int max_branch_len, StrBuf* install_dir,
           double lambda_g_err, double lambda_e_err, 
           double err_rate, int expected_covg,
           int ignore_first, int ignore_last,
           SampleModel* best_model);


void get_stats_pure_MTBC(int expected_covg, double err_rate, 
			 double lambda_g_err,double lambda_e,
			 int* arr_perc_covg, Covg* arr_median,int* arr_tkmers, 
			 int* arr_tkmers_snps, int* arr_tkmers_mobile,
			 double* arr_prop_snps, double* arr_prop_mobile,
			 int kmer_size,
			 SampleModel* sm,
			 Myc_lineage sp,
			 boolean found_ANY_myc_evidence);

void get_stats_mix_mtbc(int expected_covg, double err_rate, 
           double lambda_g_err,
           double* arr_perc_covg, double* arr_median,
           double frac_myc,
           SampleModel* sm);


void get_stats_non_MTB(int expected_covg, double err_rate, double lambda_e,
		       int* arr_perc_covg, Covg* arr_median, int* arr_tkmers, int kmer_size,
       SampleModel* sm);

Myc_lineage get_best_hit(int* arr_perc_cov, 
			 Covg* arr_median, 
			 boolean* found, 
			 boolean exclude_sp, 
			 Myc_lineage sp);

