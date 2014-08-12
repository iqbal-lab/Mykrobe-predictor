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
  species.h
*/
#include "dB_graph.h"
typedef enum 
 {
   Mtuberculosis = 0,
   Mafricanum = 1,
   Mbovis = 2 
 } Myc_species ;

#define NUM_SPECIES 3
typedef enum
  {
    PureMTB =0,
    // MajorMTBAndMinorNonMTB = 1,
    // MinorMTBAndMajorNonMTB = 2,
    NonMTB = 3,
  } SampleType;

typedef struct
{
  SampleType type;
  double likelihood;//log likelihood
  double lp; //log posterior
  double conf;
  StrBuf* name_of_non_mtb_species;//how we interpret this depends on the model type.
} SampleModel;

SampleModel* alloc_and_init_sample_model();
void free_sample_model(SampleModel* sm);

void map_species_enum_to_str(Myc_species sp, StrBuf* sbuf);
    
SampleType get_species_model(dBGraph *db_graph,int max_branch_len, StrBuf* install_dir,
           double lambda_g_err, double lambda_e_err, 
           double err_rate, int expected_covg,
           int ignore_first, int ignore_last,
           SampleModel* best_model);


void get_stats_pure_mtb(int expected_covg, double err_rate, 
         double lambda_g_err,double lambda_e,
         double* arr_perc_covg, double* arr_median,int* arr_tkmers, 
         int* arr_tkmers_snps, int* arr_tkmers_mobile,
         double* arr_prop_snps, double* arr_prop_mobile,
         int kmer_size,
         SampleModel* sm);

void get_stats_mix_mtb_and_non_mtb(int expected_covg, double err_rate, 
           double lambda_g_err,
           double* arr_perc_covg, double* arr_median,
           double frac_myc,
           SampleModel* sm);


void get_stats_non_staph(int expected_covg, double err_rate, double lambda_e,
       double* arr_perc_covg, double* arr_median, int* arr_tkmers, int kmer_size,
       SampleModel* sm);
