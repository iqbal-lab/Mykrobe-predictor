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





void get_coverage_on_best_catalayse_gene(dBGraph *db_graph,int max_branch_len,
                                        StrBuf* install_dir,int ignore_first, 
                                        int ignore_last,
                                        int* percentage_cov_cat,
                                        Covg* median_cov_cat)
{
  FILE* fp;
  AlleleInfo* ai = alloc_allele_info();
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

  // Does the data support the existance of catalayse?
  StrBuf* catalayse_fasta_file = strbuf_create(install_dir->buff);
  strbuf_append_str(catalayse_fasta_file, "data/staph/species/cat.fasta");
  fp = fopen(catalayse_fasta_file->buff, "r");
  if (fp==NULL)
    {
      die("Cannot open this file - %s", catalayse_fasta_file->buff);
    }
  // while the entry is valid iterate through the fasta file

  int num_kmers=0;
  int per_cov_cat=0;
  Covg med_cov_cat=0;
  do {
    num_kmers= get_next_single_allele_info(fp, db_graph, ai,
             true,
             seq, kmer_window,
             &file_reader_fasta,
             array_nodes, array_or, 
             working_ca, max_branch_len,
             ignore_first, ignore_last);
    if (ai->percent_nonzero >= per_cov_cat)
    {
      if (ai->median_covg_on_nonzero_nodes > med_cov_cat)
      {
        per_cov_cat=ai->percent_nonzero;;
        med_cov_cat=ai->median_covg_on_nonzero_nodes;
      }
    }
  } while ( num_kmers>0);
  *percentage_cov_cat = per_cov_cat;
  *median_cov_cat = med_cov_cat ;
  fclose(fp); 

}

boolean catalayse_exists_in_sample(dBGraph *db_graph,int max_branch_len,
                                StrBuf* install_dir,int ignore_first, 
                                int ignore_last)
{   
    int percentage_cov_cat=0;
    Covg median_cov_cat=0;
    get_coverage_on_best_catalayse_gene(db_graph,max_branch_len,
                                        install_dir,ignore_first, 
                                        ignore_last,
                                        &percentage_cov_cat,
                                        &median_cov_cat);
    if (percentage_cov_cat > 50)
    {
      return true;
    }   
    else
    {
      return false;
    }

}

boolean sample_is_staph(boolean has_catalayse)
{
  // Check if catalayse exists
  if (has_catalayse)
  {
    return true;
  }
  else
  {
    return false;
  }
}


void load_all_species_panel_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir )
{
  panel_file_paths[Saureus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Saureus], "data/staph/species/Saureus.fasta" );
  panel_file_paths[Sepidermidis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Sepidermidis], "data/staph/species/Sepidermidis.fasta" );
  panel_file_paths[Shaemolyticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shaemolyticus], "data/staph/species/Shaemolyticus.fasta" );
  panel_file_paths[Sother] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Sother], "data/staph/species/Sother.fasta" );
}

void get_coverage_on_panels(int* percentage_coverage,int* median_coverage,
                            StrBuf** panel_file_paths,
                            int max_branch_len, dBGraph *db_graph,
                            int ignore_first, int ignore_last )

{
  int i;
  FILE* fp;
  AlleleInfo* ai = alloc_allele_info();
  int number_of_reads = 0;
  int number_of_covered_reads =0; 
  int num_kmers=0;
  Covg tot_pos_kmers;
  Covg tot_kmers;
  Covg med;  
  //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL)
  {
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
  int file_reader_fasta(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry)
  {
    long long ret;
    int offset = 0;
    if (new_entry == false)
    {
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
    // Reset counters
    tot_pos_kmers = 0; 
    tot_kmers = 0; 
    med = 0;
    number_of_covered_reads = 0;
      fp = fopen(panel_file_paths[i]->buff, "r");
      if (fp==NULL)
      {
        die("Cannot open this file - %s", panel_file_paths[i]->buff);
      } 
      // while the entry is valid iterate through the fasta file
      do 
      {
        num_kmers= get_next_single_allele_info(fp, db_graph, ai,
                                               true,
                                               seq, kmer_window,
                                               &file_reader_fasta,
                                               array_nodes, array_or, 
                                               working_ca, max_branch_len,
                                               ignore_first, ignore_last);

        number_of_reads = number_of_reads + 1;
        int pos_kmers = num_kmers * (double) (ai->percent_nonzero)/100;
        if (ai->percent_nonzero > 0 ){
          number_of_covered_reads = number_of_covered_reads +1; 
        }
        //calculate a running pseudo median, before you update the tots
        if  (tot_kmers+num_kmers>0)
          {
              med += ai->median_covg ; 
              tot_kmers += num_kmers;
              tot_pos_kmers += pos_kmers;
          }
      }while ( num_kmers > 0);

    if ( (number_of_reads > 0) && (tot_kmers>0) )
    {
      percentage_coverage[i] = (int) (100 * tot_pos_kmers) / tot_kmers;
      median_coverage[i] = (double) (med) / number_of_covered_reads ;
    }
        else
    {
      percentage_coverage[i]=0;
      median_coverage[i]=0;
    } 
  fclose(fp);
  }
}



boolean is_percentage_coverage_above_threshold(int per_cov,int threshold)
{
  if (per_cov >= threshold)
  {
      return true;
  }
  else
  {
    return false;
  }
}


void find_which_panels_are_present(int* percentage_coverage,boolean* present, 
                                  int* num_panels)
{
  boolean is_aureus_present = is_percentage_coverage_above_threshold(percentage_coverage[Saureus],90);
  boolean is_epi_present = is_percentage_coverage_above_threshold(percentage_coverage[Sepidermidis],30);
  boolean is_haem_present = is_percentage_coverage_above_threshold(percentage_coverage[Shaemolyticus],30);
  boolean is_sother_present = is_percentage_coverage_above_threshold(percentage_coverage[Sother],10);  
  present[Saureus] = is_aureus_present;
  present[Sepidermidis] = is_epi_present;
  present[Shaemolyticus] = is_haem_present;
  present[Sother] = is_sother_present && !(is_epi_present || is_haem_present) ;
  int i;
  for (i=0; i<NUM_SPECIES; i++)
  {
    if (present[i])
    {
      *num_panels = *num_panels +1;
    }
  }
}

boolean coverage_exists_on_aureus_and_at_least_one_other_panel(boolean* present)
{
  // get the coverage on all our panels
  // do we have > 10% coverage on more than one panel

  boolean is_coag_neg_present = false;
  if (present[Sepidermidis] || present[Shaemolyticus] || present[Sother])
  {
    is_coag_neg_present = true;
  }
  if (present[Saureus] && is_coag_neg_present)
  {
      return true;
  }
  else
  {
    return false;
  }
}

boolean sample_is_mixed(boolean* present)
{
  boolean is_mixed = coverage_exists_on_aureus_and_at_least_one_other_panel(present);
  return (is_mixed);
}


SampleType get_sample_type(boolean has_catalayse, boolean* present)
{
  SampleType sample_type;
  if (sample_is_staph(has_catalayse) ) 
  {
      if (sample_is_mixed(present))
      {
          sample_type=MixedStaph;
      }
      else
      {
          sample_type=PureStaph;
      }
      
  }
  else
  {
      sample_type = NonStaphylococcal;
  }
  
  return sample_type;  
}



SpeciesInfo* get_species_info(dBGraph *db_graph,int max_branch_len, 
                            StrBuf* install_dir,int expected_covg,
                            int ignore_first,int ignore_last)

{
  SpeciesInfo* species_info=(SpeciesInfo *)malloc(sizeof(SpeciesInfo)); 

  StrBuf* panel_file_paths[NUM_SPECIES];
  load_all_species_panel_file_paths(panel_file_paths,install_dir);
  int percentage_coverage[NUM_SPECIES]; // for storing the percentage coverage of each reference
  int median_coverage[NUM_SPECIES]; //median covg
  boolean present[NUM_SPECIES];
  get_coverage_on_panels(percentage_coverage,median_coverage,
                          panel_file_paths,max_branch_len,db_graph,
                          ignore_first,ignore_last);
  boolean has_catalayse =  catalayse_exists_in_sample(db_graph,max_branch_len,
                                install_dir,ignore_first, 
                                ignore_last);
  int num_panels =0 ;
  find_which_panels_are_present(percentage_coverage,present,&num_panels);
  SampleType sample_type = get_sample_type(has_catalayse,present);

  species_info->sample_type = sample_type;
  species_info->num_species = num_panels;

  memcpy (species_info->present, present, sizeof(present));
  memcpy (species_info->percentage_coverage, percentage_coverage, sizeof(percentage_coverage));
  memcpy (species_info->median_coverage, median_coverage, sizeof(median_coverage));
  return species_info;
}

void map_species_enum_to_str(Staph_species staph_species, StrBuf* sbuf)
{
  if (staph_species==Saureus)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"S. aureus");
    }
  else if (staph_species==Sepidermidis)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"S. epidermidis");
    }
  else if (staph_species==Shaemolyticus)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"S. haemolyticus");
    }
  else if (staph_species==Sother)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"Unidentified Staphylococcus");
    }        
  else
    {
      die("Coding error - I would expect the compiler to prevent assigning a bad enum value - we get %d\n", staph_species);
    }
  
}


char* get_ith_species_name(SpeciesInfo* species_info, int i)
{
  Staph_species species=Saureus;
  StrBuf* pure_species_name = strbuf_new(); 
  int j; 
  int species_index = 0;
  if (i > species_info->num_species -1 ){
    die("We only have %i species, we can't find %i \n",  species_info->num_species,i+1);
  }
  else
  {
    for (j=0; j<NUM_SPECIES; j++)
    {
      if (species_info->present[j])
      {
        if (species_index == i)
        {
            // We at the required species
            species = j;
        }
        species_index = species_index + 1;
      }
    }    
  }
  map_species_enum_to_str(species, pure_species_name);
  return pure_species_name->buff;
}


int get_ith_species_coverage(SpeciesInfo* species_info,int i)
{
  int covg=0;
  int j;
  int species_index=0;
  if (i > species_info->num_species -1 ){
    die("We only have %i species, we can't find %i \n",  species_info->num_species,i+1);
  }
  else
  {
    for (j=0; j<NUM_SPECIES; j++)
    {
      if (species_info->present[j])
      {
        if (species_index == i)
        {
            // We at the required species
            covg = species_info->median_coverage[j];
        }
        species_index = species_index + 1;
      }
    }    
  }
  return covg;

}

char* get_pure_species_name(SpeciesInfo* species_info)
{
  return get_ith_species_name(species_info, 0 );
}

int get_pure_species_coverage(SpeciesInfo* species_info)
{
  return get_ith_species_coverage(species_info, 0 );
}



// Staph_species get_best_hit(int* arr_perc_cov, 
//         Covg* arr_median, 
//         boolean* found, 
//         boolean exclude_aureus)
// {
//   int i;
//   int prod=0;
//   int curr=-1;
//   for (i=0; i<NUM_SPECIES; i++)
//     {
//       if ( (exclude_aureus==true) && ((Staph_species)i==Aureus))
//  {
//    continue;
//  }
//       //      if (arr_perc_cov[i] * arr_median[i]>prod)
//       if (arr_perc_cov[i] > prod)
//  {
//    //prod = arr_perc_cov[i]* arr_median[i];
//    prod =arr_perc_cov[i];
//    curr=i;
//  }
//     }
//   if (curr==-1)
//     {
//       *found=false;
//       return Aureus;
//     }
//   else
//     {
//       *found=true;
//       return (Staph_species) curr;
//     }
// }

