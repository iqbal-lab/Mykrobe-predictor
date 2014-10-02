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
  sm->type=PureMTBC;
  sm->likelihood = 0;
  sm->lp=0;
  sm->conf=0;
  sm->name_of_non_mtb_species = strbuf_new();
  sm->name_of_pure_mtbc_species = strbuf_new();
  sm->name_of_non_mtb_lineage = strbuf_new();
  sm->name_of_pure_mtbc_lineage = strbuf_new();
  return sm;
}

void free_sample_model(SampleModel* sm)
{
  strbuf_free(sm->name_of_non_mtb_species);
  strbuf_free(sm->name_of_pure_mtbc_species);
  free(sm);
}
void map_species_enum_to_str(Myc_species sp, StrBuf* sbuf)
{
  if (sp==Mtuberculosis)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"M. Tuberculosis");
    }
  else if (sp==Mafricanum)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"M. Africanum");
    }
    else if (sp==Mbovis)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"M. Bovis");      
    }

  else
    {
      die("Coding error - I would expect the compiler to prevent assigning a bad enum value - we get %d\n", sp);
    }
  
}

void map_lineage_enum_to_str(Myc_lineage sp, StrBuf* sbuf)
{
  if (sp==beijing)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"Beijing/East Asia");
    }
  else if (sp==animal)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"Animal");
    }
  else if (sp==lineage1)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"Indian Ocean");
    }
  else if (sp==lineage2)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"East Africa / Indian ocean");
    }
  else if (sp==lineage3)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"Delhi/Central Asia");
    }
  else if (sp==lineage4)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"European/American");
    }
  else if (sp==lineage5)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"West Africa 1");
    }
  else if (sp==lineage6)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"West Africa 2");
    }
  else if (sp==lineage7)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"Ethiopian");
    }
  else
    {
      die("Coding error - I would expect the compiler to prevent assigning a bad enum value\n");
    }
  
}

Myc_species map_lineage_enum_to_species_enum(Myc_lineage sp)
{
  Myc_species species_enum;
  if (sp==beijing)
    {
      species_enum = Mtuberculosis;
    }
  else if (sp==animal)
    {
      species_enum = Mbovis;
    }
  else if (sp==lineage1)
    {
      species_enum = Mtuberculosis;
    }
  else if (sp==lineage2)
    {
      species_enum = Mtuberculosis;
    }
  else if (sp==lineage3)
    {
      species_enum = Mtuberculosis;
    }
  else if (sp==lineage4)
    {
      species_enum = Mtuberculosis;
    }
  else if (sp==lineage5)
    {
      species_enum = Mafricanum;
    }
  else if (sp==lineage6)
    {
      species_enum = Mafricanum;
    }
  else if (sp==lineage7)
    {
      species_enum = Mafricanum;
    }
  else
    {
      die("Coding error - I would expect the compiler to prevent assigning a bad enum value - get val %d for lineage\n", sp);
    }
    return species_enum;
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
  StrBuf* species_file_paths[NUM_SPECIES];
  // species_file_paths[0] = strbuf_create(install_dir->buff);
  // strbuf_append_str(species_file_paths[0], "data/tb/species/M.tuberculosis.fa");
  // species_file_paths[1] = strbuf_create(install_dir->buff);
  // strbuf_append_str(species_file_paths[1], "data/tb/species/M.africanum.fa");
  // species_file_paths[2] = strbuf_create(install_dir->buff);
  // strbuf_append_str(species_file_paths[2], "data/tb/species/M.bovis.fa");
  species_file_paths[0] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[0], "data/tb/species/beijing_sublineage.fa");
  species_file_paths[1] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[1], "data/tb/species/animal_lineage.fa");
  species_file_paths[2] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[2], "data/tb/species/lineage_1.fa");
  species_file_paths[3] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[3], "data/tb/species/lineage_2.fa");
  species_file_paths[4] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[4], "data/tb/species/lineage_3.fa");
  species_file_paths[5] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[5], "data/tb/species/lineage_4.fa");
  species_file_paths[6] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[6], "data/tb/species/lineage_5.fa");
  species_file_paths[7] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[7], "data/tb/species/lineage_6.fa");
  species_file_paths[8] = strbuf_create(install_dir->buff);
  strbuf_append_str(species_file_paths[8], "data/tb/species/lineage_7.fa");

  int i;
  int pcov[NUM_SPECIES]; // for storing the percentage coverage of each reference
  Covg mcov[NUM_SPECIES]; //median covg
  int tkmers[NUM_SPECIES];//total kmers in the unique branches
  int tkmers_snps[NUM_SPECIES];//total kmers in the unique branches which are SNPs
  int tkmers_mobile[NUM_SPECIES];
  double p_snps[NUM_SPECIES];//what propn of the reads which are SNP length, have >0 covg
  double p_mobile[NUM_SPECIES];
  int tot_pos_kmers;;
  int tot_kmers;
  Covg med;
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
  
  for (i = 0; i < NUM_SPECIES; i++)
    {
      
      
      fp = fopen(species_file_paths[i]->buff, "r");
      if (fp==NULL)
  {
    die("Cannot open this file - %s", species_file_paths[i]->buff);
  }
      
      
      // while the entry is valid iterate through the fasta file
      number_of_reads = 0;
      int num_kmers=0;
      int tot_snps=0;
      int tot_mobile=0;
      int tot_snps_pos=0;
      int tot_mobile_pos=0;

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

      if ( (pos_kmers< 0.15 * num_kmers) //ignore repeat kmers
     && ( ai->median_covg_on_nonzero_nodes > 3*err_rate * expected_covg) )
        {
    ai->median_covg=0;
    ai->percent_nonzero=0;
    pos_kmers=0;
        }

      med = (med*tot_kmers + (double)ai->median_covg * num_kmers)/(tot_kmers+num_kmers);
      tot_kmers += num_kmers;

      tot_pos_kmers += pos_kmers;

      if (num_kmers<=db_graph->kmer_size)
	{
	  if (pos_kmers>0)
	    {
	      tot_snps_pos++;
	    }
	  tot_snps++;
	}
      else
	{
	  if (pos_kmers>0)
	    {
	      tot_mobile_pos++;
	    }
	  tot_mobile++;
	}
      
      
    }
  
      } while ( num_kmers>0);
      if ( (number_of_reads>0) && (tot_kmers>0) )
	{
	  
	  pcov[i] = (int) 100*tot_pos_kmers/tot_kmers;
	  mcov[i] = med;
	  tkmers[i] = tot_pos_kmers;
	  tkmers_snps[i]=tot_snps;
	  p_snps[i]=(double)tot_snps_pos/ (double)tot_snps;
	  p_mobile[i]=(double)tot_mobile_pos/ (double)tot_mobile;
	  tkmers_mobile[i]=tot_mobile_pos;
	}
      else
	{
	  pcov[i]=0;
	  mcov[i]=0;
	  tkmers[i]=0;
	  tkmers_snps[i]=0;
	  p_snps[i]=0;
	  p_mobile[i]=0;
	  tkmers_mobile[i]=0;
	}
    }
  
  
  free_allele_info(ai);
  free(array_nodes);
  free(array_or);
  free_covg_array(working_ca);
  free(kmer_window->kmer);
  free(kmer_window);
  free_sequence(&seq);

  for (i=0; i<NUM_SPECIES; i++)
    {
      strbuf_free(species_file_paths[i]);
    }

  boolean found_ANY_myc=true;
  boolean exclude_sp=false;
  Myc_lineage best = get_best_hit(pcov, mcov, &found_ANY_myc, exclude_sp, 9999);

  SampleModel* M_pure_MTBC=alloc_and_init_sample_model();
  // SampleModel* M_pure_Mtuberculosis=alloc_and_init_sample_model();
  // SampleModel* M_pure_Mafricanum=alloc_and_init_sample_model();
  // SampleModel* M_pure_Mbovis=alloc_and_init_sample_model();
    

  // Mixed Models
  // SampleModel* M_maj_mixture=alloc_and_init_sample_model();
  // SampleModel* M_min_mixture=alloc_and_init_sample_model();
  
  SampleModel* M_non_MTB=alloc_and_init_sample_model();

  get_stats_pure_MTBC(expected_covg, err_rate,
		      lambda_g_err, lambda_e_err,
		      pcov, mcov, tkmers, tkmers_snps, tkmers_mobile,
		      p_snps, p_mobile, db_graph->kmer_size,
		      M_pure_MTBC, best, found_ANY_myc );
  // get_stats_pure_MTBC(expected_covg, err_rate,
  //     lambda_g_err, lambda_e_err,
  //     pcov, mcov, tkmers, tkmers_snps, tkmers_mobile,
  //     p_snps, p_mobile, db_graph->kmer_size,
  //     M_pure_Mafricanum ,Mafricanum );
  // get_stats_pure_MTBC(expected_covg, err_rate,
  //     lambda_g_err, lambda_e_err,
  //     pcov, mcov, tkmers, tkmers_snps, tkmers_mobile,
  //     p_snps, p_mobile, db_graph->kmer_size,
  //     M_pure_Mbovis,Mbovis);

  // get_stats_mix_mtbc(expected_covg, err_rate,
  //       lambda_g_err,
  //       pcov, mcov, 
  //       0.9, M_maj_mixture);
  // get_stats_mix_mtbc(expected_covg, err_rate,
  //       lambda_g_err,
  //       pcov, mcov, 
  //       0.1, M_min_mixture);
  
  get_stats_non_MTB(expected_covg, err_rate,lambda_e_err,
		    pcov, mcov, tkmers, db_graph->kmer_size, M_non_MTB);

  SampleModel* marray[2] = {M_pure_MTBC, M_non_MTB};
  printf("M_pure_MTBC %f\n", M_pure_MTBC->likelihood);
  printf("M_non_MTB %f\n", M_non_MTB->likelihood);
  qsort(marray, 2, sizeof(SampleModel*), sample_model_cmp_logpost);
  best_model->conf = marray[1]->lp - marray[0]->lp;
  best_model->type = marray[1]->type;
  best_model->likelihood = marray[1]->likelihood;
  best_model->lp = marray[1]->lp;
  strbuf_reset(best_model->name_of_non_mtb_species);
  strbuf_reset(best_model->name_of_pure_mtbc_species);

  strbuf_append_str(best_model->name_of_non_mtb_species, 
        marray[1]->name_of_non_mtb_species->buff);
  strbuf_append_str(best_model->name_of_pure_mtbc_species, 
        marray[1]->name_of_pure_mtbc_species->buff);
  
  strbuf_append_str(best_model->name_of_non_mtb_lineage, 
        marray[1]->name_of_non_mtb_lineage->buff);
  strbuf_append_str(best_model->name_of_pure_mtbc_lineage, 
        marray[1]->name_of_pure_mtbc_lineage->buff);
  //cleanup
  // free_sample_model(M_pure_Mtuberculosis);
  // free_sample_model(M_pure_Mafricanum);
  // free_sample_model(M_pure_Mbovis);
  free_sample_model(M_pure_MTBC);

  // free_sample_model(M_maj_mixture);
  // free_sample_model(M_min_mixture);

  free_sample_model(M_non_MTB);

  return best_model->type;
  
}


Myc_lineage get_best_hit(int* arr_perc_cov, 
			 Covg* arr_median, 
			 boolean* found, 
			 boolean exclude_sp, 
			 Myc_lineage sp)
{
  int i;
  int best_perc_cov_so_far=0;
  int best_median_cov_so_far=0;

  int curr=-1;
  for (i=0; i<NUM_SPECIES; i++)
    {
      if ( (exclude_sp==true) && ((Myc_lineage)i==sp))
	{
	  continue;
	}
      if (arr_perc_cov[i] >= best_perc_cov_so_far)
	{
    best_perc_cov_so_far = arr_perc_cov[i];
    // Only update if the median coverage has also improved
    if (arr_median[i] > best_median_cov_so_far){
      best_median_cov_so_far = arr_median[i];
      curr=i;
    } 
	}
    printf("Lineage : %i Per_cov :  %i median_covg : %i \n",i,arr_perc_cov[i],arr_median[i]);
    }
  if (curr==-1)
    {
      *found=false;
      return sp;
    }
  else
    {
      printf("Best Lineage : %i \n",curr);
      *found=true;
      return (Myc_lineage) curr;
    }
}

void get_stats_pure_MTBC(int expected_covg, double err_rate, 
			 double lambda_g_err,double lambda_e,
			 int* arr_perc_covg, Covg* arr_median, int* arr_tkmers, 
			 int* arr_tkmers_snps, int* arr_tkmers_mobile, 
			 double* arr_prop_snps, double* arr_prop_mobile,
			 int kmer_size,
			 SampleModel* sm,
			 Myc_lineage sp,//sp is the best of ALL lineages. We can assume it is one of our 9
			 boolean found_ANY_myc_evidence)

{

  if (found_ANY_myc_evidence==false)
    {
      sm->likelihood = -999999;
      sm->lp= -999999;
      sm->type=PureMTBC;
      sm->conf=9999999;
      return;
    }
  //which non-MTB
  boolean found=true;
  boolean exclude_sp=true;
  Myc_lineage best = get_best_hit(arr_perc_covg, arr_median, &found, exclude_sp, sp);//best is the best of the rest
  boolean no_evidence_for_ANY_alternate=false;
  if (found==false)
    {
      no_evidence_for_ANY_alternate=true;
    }


  //now deal with sp
  int recovery_expected = 100*(1-exp(-expected_covg));
  double lpr=0;

  if (arr_perc_covg[sp] > 0.75*recovery_expected)
    {
      lpr=0;
    }
  else if (arr_perc_covg[sp] > 0.5*recovery_expected)
    {
      lpr=-1000;
    }
  else if (arr_perc_covg[sp] > 0.1*recovery_expected)
    {
      lpr=-10000;
    }
  else
    {
      lpr=-99999999;
    }
    printf("lpr: %f \n", lpr);


  double llk = -lambda_g_err 
    + arr_median[sp]*log(lambda_g_err) 
    - log_factorial(arr_median[sp]);


  //  llk = - (double) arr_tkmers[sp] * ((double) (100-arr_perc_covg[sp]/100) * recovery_expected; //prob of a gap of that length

  //now we need to account for coverage on non-aureus
  //must be due to
  // a) SNP error (single base errors)s
  // b) plasmids/mobile elements. What % of our unique-to-species panel will be mobile? Say max 40%

  double llke=0;
  if (no_evidence_for_ANY_alternate==false)
    {
      int numk;
      if (arr_tkmers_snps[best]>kmer_size)
	{
	  numk=(int) (arr_tkmers_snps[best]/kmer_size);
	}
      else
	{
	  numk=1;
	}
      numk += arr_tkmers_mobile[best];
      
      
      llke +=  -lambda_e 
	+ numk*arr_median[best]*log(lambda_e)
	-log_factorial(numk*arr_median[best]);
    }
  double lpre=0;

  sm->likelihood = llk+llke;
  sm->lp= sm->likelihood +lpr+lpre;
  map_lineage_enum_to_str(sp, sm->name_of_pure_mtbc_lineage);
  Myc_species species_enum = map_lineage_enum_to_species_enum(sp);
  map_species_enum_to_str(species_enum, sm->name_of_pure_mtbc_species);
  sm->type=PureMTBC;
}


/*
void get_stats_mix_mtbc(int expected_covg, double err_rate, double lambda_g_err,
           double* arr_perc_covg, double* arr_median, 
           double frac_MTB,
           SampleModel* sm)
{
  if ( (frac_MTB<0.001) || (frac_MTB>0.999) )
    {
      die("Do not call get_stats_mix_mtb_and_non_mtb when frac==0 or 1 -= programming error\n");
    }
  boolean found=true;
  // Get the best 2 species, if it's a mixture it will be between these two. 
  boolean exclude=false;
  int best = get_best_hit(arr_perc_covg, arr_median, &found, exclude,9999);
  // Get the best hit of the remainder
  int second_best = get_best_hit(arr_perc_covg,arr_median, &found, true,best);

  // The naming of these is confusing but for now its just a way of storing both species.
  map_lineage_enum_to_str((Myc_lineage) best, sm->name_of_pure_mtbc_species);
  map_lineage_enum_to_str((Myc_lineage) second_best, sm->name_of_non_mtb_species);  

  if (found==false)
    {
      sm->likelihood=-9999999;
      sm->lp=-99999999;
    }
  else
    {
      double MTB_recovery_expected = frac_MTB*(1-exp(-expected_covg));

      double MTB_lpr;
      if (arr_perc_covg[best] > 0.75*MTB_recovery_expected)
  {
    MTB_lpr=log(1);
  }
      else if (arr_perc_covg[best] > 0.5*MTB_recovery_expected)
        {
          MTB_lpr=log(0.5);
        }
      else
  {
    MTB_lpr=-999999;
  }
      double llk_aureus = -frac_MTB*lambda_g_err + arr_median[best]*log(frac_MTB*lambda_g_err) - log_factorial(arr_median[best]);
      //now do the same for the other population
      double cong_recovery_expected = (1-frac_MTB)*(1-exp(-expected_covg));

      double cong_lpr;
      //want to avoid calling CONG just because of small number of repeat kmers
      printf("Got %f perc covg of the minor and we expect %f\n", arr_perc_covg[second_best], cong_recovery_expected );
      //      if (arr_perc_covg[second_best] < 0.3)
      //{
    //cong_lpr=-999999;
    //}
      if (arr_perc_covg[second_best] > 0.9*cong_recovery_expected)
  {
    cong_lpr=0;
  }
      else if (arr_perc_covg[second_best] > 0.75*cong_recovery_expected)
  {
    cong_lpr=log(0.1);
  }
      else
  {
    cong_lpr=-999999;
  }
      double llk_cong = -(1-frac_MTB)*lambda_g_err 
  + arr_median[second_best]*log((1-frac_MTB)*lambda_g_err) 
  - log_factorial(arr_median[second_best]);


      sm->likelihood = llk_aureus+llk_cong;
      sm->lp =sm->likelihood+MTB_lpr+cong_lpr;;

    }

    sm->type=MixedMTB;

}
*/

void get_stats_non_MTB(int expected_covg, double err_rate, double lambda_e,
		       int* arr_perc_covg, Covg* arr_median, int* arr_tkmers,
		       int kmer_size,
		       SampleModel* sm)
{
  strbuf_append_str(sm->name_of_non_mtb_species, "Non-MTBC");
  boolean found=true;
  boolean exclude_mtb=false;
  //do ANY MTB get a decent hit?
  int best = get_best_hit(arr_perc_covg, arr_median, &found, exclude_mtb, 9999);

  if (found==false)
    {
      //found no evidence of any MTB at all
      sm->likelihood=-lambda_e;
      sm->lp=sm->likelihood;
    }
  else
    {
      //all coverage must be errors
      int numk;
      if (arr_tkmers[best]>kmer_size)
	{
	  numk=arr_tkmers[best]-kmer_size;
	}
      else
	{
	  numk=1;
	}
      
      //      double t = numk*arr_median[best]*arr_perc_covg[best]/100;
      double llk = -lambda_e 
	+  numk*arr_median[best]*log(lambda_e) 
	-log_factorial(numk*arr_median[best]);
      double lpr;
      if ( (arr_perc_covg[best]>0.1) && (arr_median[best]>0.1*expected_covg) )
	{
	  lpr=-9999999;
	}
      else if (arr_perc_covg[best]>0.05)
	{
	  lpr=-log(100);
	}
      else
	{
	  lpr=log(1);
	}
      sm->likelihood=llk;
      sm->lp = llk+lpr;
    }  
  sm->type = NonMTB;
}


