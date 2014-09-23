/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 * **********************************************************************
 *
 * This file is part of Mykrobe.
 *
 * Mykrobe is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mykrobe is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mykrobe.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */
/*
  species.c
*/


// system headers
#include <stdlib.h>
#include <limits.h>

// third party headers
#include <string_buffer.h>

#include "build.h" 
#include "maths.h" 
#include "element.h"
#include "seq.h"
#include "open_hash/hash_table.h"
#include "dB_graph.h"
#include "species.h"
#include "gene_presence.h"
#include "genotyping_known.h"

SampleModel* alloc_and_init_sample_model()
{
  SampleModel* sm = (SampleModel*) malloc(sizeof(SampleModel));
  if (sm==NULL)
    {
      die("Cannot malloc tiny object - your computer must be about to die\n");
    }
  sm->type=PureStaphAureus;
  sm->likelihood = 0;
  sm->lp=0;
  sm->conf=0;
  sm->name_of_non_aureus_species = strbuf_new();
  return sm;
}

void free_sample_model(SampleModel* sm)
{
  strbuf_free(sm->name_of_non_aureus_species);
  free(sm);
}
void map_species_enum_to_str(Staph_species sp, StrBuf* sbuf)
{
  if (sp==Scapitis)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"S.capitis");
    }
  else if (sp==Scaprae)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"S.caprae");
    }
  else if (sp==Sepidermidis)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"S.epidermidis");
    }
  else if (sp==Sequorum)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"S.equorum");
    }
  else if (sp==Shaemolyticus)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"S.haemolyticus");
    }
  else if (sp==Shominis)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"S.hominis");
    }
  else if (sp==Slugdunensis)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"S.lugdunensis");
    }
  else if (sp==Smassiliensis)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"S.massiliensis");
    }
  else if (sp==Spettenkofer)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"S.pettenkofer");
    }
  else if (sp==Spseudintermedius)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"S.pseudintermedius");
    }
  else if (sp==Ssaprophyticus)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"S.saprophyticus");
    }
  else if (sp==Ssimiae)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"S.simiae");
    }
  else if (sp==Ssimulans)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"S.simulans");
    }
  else if (sp==Ssphgb0015)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"Staphylococcus sp. HGB0015");
    }
  else if (sp==Sspoj82)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"Staphylococcus sp. OJ82");
    }
  else if (sp==Aureus)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"S.aureus");
    }
  else if (sp==Swarneri)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"S.warneri");
    }
  else
    {
      die("Coding error - I would expect the compiler to prevent assigning a bad enum value - got %d\n", (int) sp);
    }
  
}



int sample_model_cmp_logpost(const void *a, const void *b)
{
  // casting pointer types
  const SampleModel *ia = *(const SampleModel **)a;
  const SampleModel *ib = *(const SampleModel **)b;


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

SampleType get_species_model(dBGraph *db_graph,int max_branch_len, StrBuf* install_dir,
			     double lambda_g_err, double lambda_e_err, double err_rate,
			     int expected_covg,
			     int ignore_first, int ignore_last,
			     SampleModel* best_model)

{
  // Define the paths to the possible species
  StrBuf* species_file_paths[17];
  species_file_paths[0] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[0], "data/staph/species/Scapitis_unique_branches.fasta");
  species_file_paths[1] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[1], "data/staph/species/Scaprae_unique_branches.fasta");
  species_file_paths[2] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[2], "data/staph/species/Sepidermidis_unique_branches.fasta");
  species_file_paths[3] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[3], "data/staph/species/Sequorum_unique_branches.fasta");
  species_file_paths[4] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[4], "data/staph/species/Shaemolyticus_unique_branches.fasta");
  species_file_paths[5] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[5], "data/staph/species/Shominis_unique_branches.fasta");
  species_file_paths[6] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[6], "data/staph/species/Slugdunensis_unique_branches.fasta");
  species_file_paths[7] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[7], "data/staph/species/Smassiliensis_unique_branches.fasta");
  species_file_paths[8] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[8], "data/staph/species/Spettenkofer_unique_branches.fasta");
  species_file_paths[9] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[9], "data/staph/species/Spseudintermedius_unique_branches.fasta");
  species_file_paths[10] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[10], "data/staph/species/Ssaprophyticus_unique_branches.fasta");
  species_file_paths[11] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[11], "data/staph/species/Ssimiae_unique_branches.fasta");
  species_file_paths[12] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[12], "data/staph/species/Ssimulans_unique_branches.fasta");
  species_file_paths[13] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[13], "data/staph/species/S_sp_hgb0015_unique_branches.fasta");
  species_file_paths[14] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[14], "data/staph/species/S_sp_oj82_unique_branches.fasta");
  species_file_paths[15] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[15], "data/staph/species/staph_unique_branches.fasta");
  species_file_paths[16] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[16], "data/staph/species/S_warneri_unique_branches.fasta");

  StrBuf* non_staph_fa = strbuf_create(install_dir->buff);
  strbuf_append_str(non_staph_fa, "data/staph/species/cat.fasta");

  int i;
  double pcov[17]; // for storing the percentage coverage of each reference
  double pcov_cat=0;//for checking if covg on catalase
  double mcov[17]; //median covg
  double mcov_cat=0;
  int tkmers[17];//total kmers in the unique branches
  int tkmers_cat=0;
  double tot_pos_kmers;
  double tot_kmers;
  double med;

  int number_of_reads;


  AlleleInfo* ai = alloc_allele_info();
  FILE* fp;

  //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence");
  }
  alloc_sequence(seq,max_branch_len,LINE_MAX);

  //We are going to load all the bases into a single sliding window 
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      die("Failed to malloc kmer sliding window");
    }


  CovgArray* working_ca = alloc_and_init_covg_array(max_branch_len);
  dBNode** array_nodes = (dBNode**) malloc(sizeof(dBNode*)*max_branch_len);
  Orientation* array_or =(Orientation*)  malloc(sizeof(Orientation)*max_branch_len);
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_branch_len-db_graph->kmer_size+1));
  if (kmer_window->kmer==NULL)
    {
      die("Failed to malloc kmer_window->kmer");
    }
  kmer_window->nkmers=0;
  //create file readers
  int file_reader_fasta(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;
    if (new_entry == false){
      offset = db_graph->kmer_size;
    }
    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);
    
    return ret;
  }
  
  
  
  if ( (array_nodes==NULL) || (array_or==NULL))
    {
      die("Cannot alloc array of nodes or of orientations");
    }
  
  for (i = 0; i < 17; i++)
    {
      
      fp = fopen(species_file_paths[i]->buff, "r");
      if (fp==NULL)
	{
	  die("Cannot open this file - %s", species_file_paths[i]->buff);
	}
      
      
      // while the entry is valid iterate through the fasta file
      number_of_reads = 0;
      int num_kmers=0;

      tot_pos_kmers = 0;
      tot_kmers=0;
      med=0;
      do {
	
	num_kmers= get_next_single_allele_info(fp, db_graph, ai,
					       true,
					       seq, kmer_window,
					       &file_reader_fasta,
					       array_nodes, array_or, 
					       working_ca, max_branch_len,
					       ignore_first, ignore_last);

	number_of_reads = number_of_reads + 1;
	int pos_kmers = num_kmers * (double) (ai->percent_nonzero)/100;
	//calculate a running pseudo median, before you update the tots
	if  (tot_kmers+num_kmers>0)
	  {

	    if ( (pos_kmers< 0.6 * num_kmers) //ignore repeat kmers, which might give high covg to a small fraction of the contig
		 && ( ai->median_covg_on_nonzero_nodes > 0.5* expected_covg) )
	      {
		ai->median_covg=0;
		ai->percent_nonzero=0;
		pos_kmers=0;
	      }

	    med = (med*tot_kmers + (double)ai->median_covg * num_kmers)/(tot_kmers+num_kmers);
	    tot_kmers += num_kmers;

	    tot_pos_kmers += pos_kmers;
	  }

      } while ( num_kmers>0);
      if ( (number_of_reads>0) && (tot_kmers>0) )
	{
	  pcov[i] = tot_pos_kmers/tot_kmers;
	  mcov[i] = med;
	  tkmers[i] = tot_kmers;
	}
      else
	{
	  pcov[i]=0;
	  mcov[i]=0;
	  tkmers[i]=0;
	}      
      fclose(fp);
    }


  //now get the data on catalase, which will be used for the non-staph model

  fp = fopen(non_staph_fa->buff, "r");
  if (fp==NULL)
    {
      die("Cannot open this file - %s", non_staph_fa->buff);
    }
      
  // while the entry is valid iterate through the fasta file
  number_of_reads = 0;
  int num_kmers=0;
  
  do {
    
    num_kmers= get_next_single_allele_info(fp, db_graph, ai,
					   true,
					   seq, kmer_window,
					   &file_reader_fasta,
					   array_nodes, array_or, 
					   working_ca, max_branch_len,
					   ignore_first, ignore_last);
    
    number_of_reads = number_of_reads + 1;
    //    int pos_kmers = num_kmers * (double) (ai->percent_nonzero)/100;
    //calculate a running pseudo median, before you update the tots
    if  (tot_kmers+num_kmers>0)
      {
	
	/*
	if ( (pos_kmers< 0.6 * num_kmers) //ignore repeat kmers, which might give high covg to a small fraction of the contig
	     && ( ai->median_covg_on_nonzero_nodes > 0.5* expected_covg) )
	  {
	    ai->median_covg=0;
	    ai->percent_nonzero=0;
	    pos_kmers=0;
	    }*/
	
	if (ai->percent_nonzero>pcov_cat)
	  {
	    pcov_cat=ai->percent_nonzero;;
	    mcov_cat=ai->median_covg_on_nonzero_nodes;
	  }
      }
    
  } while ( num_kmers>0);
  fclose(fp);

  
  
  free_allele_info(ai);
  free(array_nodes);
  free(array_or);
  free_covg_array(working_ca);
  free(kmer_window->kmer);
  free(kmer_window);
  free_sequence(&seq);

  for (i=0; i<17; i++)
    {
      strbuf_free(species_file_paths[i]);
    }


  
  SampleModel* M_pure_sa=alloc_and_init_sample_model();
  SampleModel* M_maj_sa=alloc_and_init_sample_model();
  SampleModel* M_min_sa=alloc_and_init_sample_model();
  SampleModel* M_non_staph=alloc_and_init_sample_model();

  get_stats_pure_aureus(expected_covg, err_rate,
			lambda_g_err, lambda_e_err,
			pcov, mcov, tkmers,
			db_graph->kmer_size,
			M_pure_sa); 
  /*  get_stats_mix_aureus_and_CONG(expected_covg, err_rate,
				lambda_g_err,
				pcov, mcov, tkmers, 
				0.99, M_pure_sa);  */
  get_stats_mix_aureus_and_CONG(expected_covg, err_rate,
				lambda_g_err,
				pcov, mcov, tkmers, 
				0.95, M_maj_sa);
  get_stats_mix_aureus_and_CONG(expected_covg, err_rate,
				lambda_g_err,
				pcov, mcov, tkmers,
				0.05, M_min_sa);

  get_stats_non_staph(expected_covg, err_rate,lambda_e_err,
		      pcov_cat, mcov_cat, tkmers_cat, 
		      db_graph->kmer_size, M_non_staph);

  SampleModel* marray[4] = {M_pure_sa, M_maj_sa, M_min_sa, M_non_staph};
  qsort(marray, 4, sizeof(SampleModel*), sample_model_cmp_logpost);
  best_model->conf = marray[3]->lp - marray[2]->lp;
  best_model->type = marray[3]->type;
  best_model->likelihood = marray[3]->likelihood;
  best_model->lp = marray[3]->lp;
  strbuf_reset(best_model->name_of_non_aureus_species);
  strbuf_append_str(best_model->name_of_non_aureus_species, 
		    marray[3]->name_of_non_aureus_species->buff);
  //cleanup
  free_sample_model(M_pure_sa);
  free_sample_model(M_maj_sa);
  free_sample_model(M_min_sa);
  free_sample_model(M_non_staph);

  return best_model->type;
  
}


Staph_species get_best_hit(double* arr_perc_cov, 
			   double* arr_median, 
			   boolean* found, 
			   boolean exclude_aureus)
{
  int i;
  double prod=0;
  int curr=-1;
  for (i=0; i<NUM_SPECIES; i++)
    {
      if ( (exclude_aureus==true) && ((Staph_species)i==Aureus))
	{
	  continue;
	}
      //      if (arr_perc_cov[i] * arr_median[i]>prod)
      if (arr_perc_cov[i] > prod)
	{
	  //prod = arr_perc_cov[i]* arr_median[i];
	  prod =arr_perc_cov[i];
	  curr=i;
	}
    }
  if (curr==-1)
    {
      *found=false;
      return Aureus;
    }
  else
    {
      *found=true;
      return (Staph_species) curr;
    }
}

void get_stats_pure_aureus(int expected_covg, double err_rate, 
			   double lambda_g_err,double lambda_e,
			   double* arr_perc_covg, double* arr_median, int* arr_tkmers, 
			   int kmer_size,
			   SampleModel* sm)

{

  //which Staph non-aureus is best hit
  boolean found=true;
  boolean exclude_aureus=true;
  int best = get_best_hit(arr_perc_covg, arr_median, &found, exclude_aureus);


  //now deal with aureus
  double recovery_expected = 1-exp(-expected_covg);
  
  double lpr=0;

  if (arr_perc_covg[Aureus] > 0.75*recovery_expected)
    {
      lpr=log(1);
    }
  else if (arr_perc_covg[Aureus] > 0.5*recovery_expected)
    {
      lpr=log(0.005);
    }
  else if (arr_perc_covg[Aureus] > 0.3*recovery_expected)
    {
      lpr=log(0.0001);
    }
  else if (arr_perc_covg[Aureus] > 0.1*recovery_expected)
    {
      lpr=-1000;
    }
  else
    {
      lpr=-99999999;
    }


  double llk = -lambda_g_err 
    + arr_median[Aureus]*log(lambda_g_err) 
    - log_factorial(arr_median[Aureus]);


  //  llk = - (double) arr_tkmers[Aureus] * ((double) (100-arr_perc_covg[Aureus]/100) * recovery_expected; //prob of a gap of that length

  //now we need to account for coverage on non-aureus
  double numk = arr_median[best];

  //now - in this model we expect errors from Aureus to give covg on cong.
  double exp_extra_cov =  (arr_median[Aureus] * err_rate) ;
  if (numk> exp_extra_cov)
    {
      numk -= exp_extra_cov;
    }
  else
    {
      numk=0;
    }
  //  numk=1; //arr_tkmers[best];//debug
  double llke =  -lambda_e 
    + numk *log(lambda_e)
    -log_factorial(numk);

  double lpre=0;

  lpre=0;
  
  sm->likelihood = llk+llke;
  sm->lp= sm->likelihood +lpr+lpre;
  map_species_enum_to_str(Aureus, sm->name_of_non_aureus_species);
  sm->type=PureStaphAureus;
}



//CONG=coag neg
void get_stats_mix_aureus_and_CONG(int expected_covg, double err_rate, double lambda_g_err, 
				   double* arr_perc_covg, double* arr_median, int* arr_tkmers,
				   double frac_aureus,
				   SampleModel* sm)
{
  if ( (frac_aureus<0.001) || (frac_aureus>0.999) )
    {
      die("Do not call get_stats_mix_aureus_and_CONG when frac==0 or 1 -= programming error\n");
    }
  boolean found=true;
  boolean exclude_aureus=true;
  int best = get_best_hit(arr_perc_covg, arr_median, &found, exclude_aureus);
  if (found==true)
    {
      map_species_enum_to_str((Staph_species) best, sm->name_of_non_aureus_species);
    }
  else
    {
      strbuf_append_str(sm->name_of_non_aureus_species, "Failed to find coag-neg - this is default text only - if you see it, it is a bug\n");
    }

  if (found==false)
    {
      sm->likelihood=-9999999;
      sm->lp=-99999999;
    }
  else
    {
      double aureus_recovery_expected = 1-exp(-frac_aureus*expected_covg);
      double lambda_aureus = lambda_g_err*frac_aureus;
      double numk=arr_median[Aureus];
      if (frac_aureus<0.5)
	{
	  //get a bit more covg from errors on the cong (major pop)
	    aureus_recovery_expected 
	    += (1-frac_aureus)*err_rate/(3*(1-err_rate));
	    lambda_aureus += lambda_g_err*(1-frac_aureus)*err_rate/(3*(1-err_rate)); 

	    //aome of the covg on Aureus is due to errors on the other, dont penalise
	    double exp_extra_cov =  (arr_median[best] * err_rate) ;
	    if (numk> exp_extra_cov)
	      {
		numk -= exp_extra_cov;
	      }
	    else
	      {
		numk=0;
	      }

	}


      double aureus_lpr=0;

      if (arr_perc_covg[Aureus] > 0.75*aureus_recovery_expected)
	{
	  aureus_lpr=log(1);
	}
      else if (arr_perc_covg[Aureus] > 0.5*aureus_recovery_expected)
        {
          aureus_lpr=log(0.005);
        }
      else if (arr_perc_covg[Aureus] > 0.3*aureus_recovery_expected)
        {
          aureus_lpr=log(0.0001);
        }
      else if (arr_perc_covg[Aureus] > 0.1*aureus_recovery_expected)
        {
          aureus_lpr=-1000;
        }
      else
	{
	  aureus_lpr=-999999;
	}

      double llk_aureus = -lambda_aureus + numk*log(lambda_aureus) - log_factorial(numk);


      //now do the same for the other population
      double cong_recovery_expected = 1-exp(-(1-frac_aureus)*expected_covg);
      double lambda_cong = lambda_g_err*(1-frac_aureus);
      numk=arr_median[best];
      if (frac_aureus>0.5)
	{
	    cong_recovery_expected 
	    += frac_aureus*err_rate/(3*(1-err_rate));
	    lambda_cong += lambda_g_err*frac_aureus*err_rate/(3*(1-err_rate));

	    double exp_extra_cov =  (arr_median[Aureus] * err_rate) ;
	    if (numk> exp_extra_cov)
	      {
		numk -= exp_extra_cov;
	      }
	    else
	      {
		numk=0;
	      }
	}



      double cong_lpr=0;
      //want to avoid calling CONG just because of small number of repeat kmers
      
      if (arr_perc_covg[best] > 0.75*cong_recovery_expected)
	{
	  cong_lpr=0;
	}
      else if (arr_perc_covg[best] > 0.5*cong_recovery_expected)
	{
	  cong_lpr=log(0.005);
	}
      else if (arr_perc_covg[best] > 0.3*cong_recovery_expected)
	{
	  cong_lpr=log(0.0001);
	}
      else if (arr_perc_covg[best] > 0.1*cong_recovery_expected)
	{
	  cong_lpr=-1000;
	}
      else
	{
	  cong_lpr=-999999;
	}

      double llk_cong = -lambda_cong + numk*log(lambda_cong)-log_factorial(numk);


      sm->likelihood = llk_aureus+llk_cong;
      sm->lp =sm->likelihood+aureus_lpr+cong_lpr;;

    }
  if (frac_aureus>0.5)
    {
      sm->type=MajorStaphAureusAndMinorNonCoag;
    }
  else
    {
      sm->type = MinorStaphAureusAndMajorNonCoag;
    }
}


void get_stats_non_staph(int expected_covg, double err_rate, double lambda_e,
			 double perc_covg_cat, double median_cat, int tkmers_cat,
			 int kmer_size,
			 SampleModel* sm)
{
  // printf("Got %f percent of cat, median covg %f\n", perc_covg_cat, median_cat);
  strbuf_append_str(sm->name_of_non_aureus_species, "Non-staphylococcal");
  
  //any coverage on catalase gene must be error
  double llk = -lambda_e 
    +  median_cat*log(lambda_e) 
    -log_factorial(median_cat);
  double lpr=0;
      
  if (perc_covg_cat>50)
    {
      lpr=-9999999;
    }
  else if (perc_covg_cat>20)
    {
      lpr=-1000;
    }
  else
    {
      lpr=log(1);
    }
  
  sm->likelihood=llk;
  sm->lp = llk+lpr;
  sm->type = NonStaphylococcal;
}

