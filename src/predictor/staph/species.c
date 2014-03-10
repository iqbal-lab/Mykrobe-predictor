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
  species.c
*/


// system headers
#include <stdlib.h>
#include <limits.h>

// third party headers
#include <string_buffer.h>

#include "build.h" 
#include "element.h"
#include "seq.h"
#include "open_hash/hash_table.h"
#include "dB_graph.h"
#include "species.h"
#include "gene_presence.h"
Staph_species get_species(dBGraph *db_graph,int max_gene_len )
{
  // Define the paths to the possible species
  const char *species_file_paths[17];
  species_file_paths[0] = "../data/species/Scapitis_unique_branches.fasta";
  species_file_paths[1] = "../data/species/Scaprae_unique_branches.fasta";
  species_file_paths[2] = "../data/species/Sepidermidis_unique_branches.fasta";
  species_file_paths[3] = "../data/species/Sequorum_unique_branches.fasta";
  species_file_paths[4] = "../data/species/Shaemolyticus_unique_branches.fasta";
  species_file_paths[5] = "../data/species/Shominis_unique_branches.fasta";
  species_file_paths[6] = "../data/species/Slugdunensis_unique_branches.fasta";
  species_file_paths[7] = "../data/species/Smassiliensis_unique_branches.fasta";
  species_file_paths[8] = "../data/species/Spettenkofer_unique_branches.fasta";
  species_file_paths[9] = "../data/species/Spseudintermedius_unique_branches.fasta";
  species_file_paths[10] = "../data/species/Ssaprophyticus_unique_branches.fasta";
  species_file_paths[11] = "../data/species/Ssimiae_unique_branches.fasta";
  species_file_paths[12] = "../data/species/Ssimulans_unique_branches.fasta";
  species_file_paths[13] = "../data/species/S_sp_hgb0015_unique_branches.fasta";
  species_file_paths[14] = "../data/species/S_sp_oj82_unique_branches.fasta";
  species_file_paths[15] = "../data/species/staph_unique_branches.fasta";
  species_file_paths[16] = "../data/species/S_warneri_unique_branches.fasta";

  int i;char* tmpName;
  double pcov[17]; // for storing the percentage coverage of each reference
  int sumpcov;
  int number_of_reads;


  GeneInfo* gi = alloc_and_init_gene_info();
  FILE* fp;

  //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence");
  }
  alloc_sequence(seq,max_gene_len,LINE_MAX);

  //We are going to load all the bases into a single sliding window 
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      die("Failed to malloc kmer sliding window");
    }

    //  int max_gene_len = 5000;
  CovgArray* working_ca = alloc_and_init_covg_array(max_gene_len);
  dBNode** array_nodes = (dBNode**) malloc(sizeof(dBNode*)*max_gene_len);
  Orientation* array_or =(Orientation*)  malloc(sizeof(Orientation)*max_gene_len);
  for (i = 0; i < 17; i++)
  {


    fp = fopen(species_file_paths[i], "r");
    if (fp==NULL)
      {
        die("Cannot open this file");
      }
    
    kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_gene_len-db_graph->kmer_size+1));
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
    
    // while the entry is valid iterate through the fasta file
    number_of_reads = 0;
    sumpcov = 0;
    do {
        get_next_gene_info(fp, db_graph, gi,
         seq, kmer_window,
         &file_reader_fasta,
         array_nodes, array_or, 
         working_ca, max_gene_len);
         tmpName = gi->strbuf->buff; 
         number_of_reads = number_of_reads + 1;
         sumpcov = sumpcov + gi->percent_nonzero;

    } while ( strlen(tmpName) != 0);
    pcov[i] = sumpcov / number_of_reads;

  }

  free_gene_info(gi);
  free(array_nodes);
  free(array_or);
  free_covg_array(working_ca);
  free(kmer_window->kmer);
  free(kmer_window);
  free_sequence(&seq);

  // Look at the max of the pcov
  int c,location;
  double maximum;
  for (c = 0; c < 17; c++)
    {
      if (pcov[c] > maximum)
      {
         maximum  = pcov[c];
         location = c;
      }
    }
  printf("Maximum element is present at location %i and it's value is %f.\n", location, maximum);

  Staph_species species_out =  location;
  return species_out;

}
