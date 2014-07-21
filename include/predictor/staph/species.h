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
  Scapitis = 0,
  Scaprae = 1,
  Sepidermidis = 2,
  Sequorum =3,
  Shaemolyticus = 4,
  Shominis = 5,
  Slugdunensis = 6,
  Smassiliensis = 7,
  Spettenkofer = 8,
  Spseudintermedius = 9,
  Ssaprophyticus = 10,
  Ssimiae = 11,
  Ssimulans = 12,
  Ssphgb0015 = 13,
  Sspoj82 = 14,
  Aureus = 15,
  Swarneri = 16,
  } Staph_species ;

#define NUM_SPECIES 17

typedef enum
  {
    PureStaphAureus =0,
    MajorStaphAureusAndMinorNonCoag = 1,
    MinorStaphAureusAndMajorNonCoag = 2,
    NonStaphylococcal = 3,
  } SampleType;

typedef struct
{
  SampleType type;
  double likelihood;//log likelihood
  double lp; //log posterior
  double conf;
  StrBuf* name_of_non_aureus_species;//how we interpret this depends on the model type.
} SampleModel;

SampleModel* alloc_and_init_sample_model();
void free_sample_model(SampleModel* sm);

void map_species_enum_to_str(Staph_species sp, StrBuf* sbuf);
		
SampleType get_species_model(dBGraph *db_graph,int max_branch_len, StrBuf* install_dir,
			     double lambda_g_err, double lambda_e_err, 
			     double err_rate, int expected_covg,
			     int ignore_first, int ignore_last,
			     SampleModel* best_model);


void get_stats_pure_aureus(int expected_covg, double err_rate, 
			   double lambda_g_err,double lambda_e,
			   double* arr_perc_covg, double* arr_median,int* arr_tkmers,
			   SampleModel* sm);

void get_stats_mix_aureus_and_CONG(int expected_covg, double err_rate, 
				   double lambda_g_err,
				   double* arr_perc_covg, double* arr_median,
				   double frac_aureus,
				   SampleModel* sm);


void get_stats_non_staph(int expected_covg, double err_rate, double lambda_e,
			 double* arr_perc_covg, double* arr_median, int* arr_tkmers,
			 SampleModel* sm);



