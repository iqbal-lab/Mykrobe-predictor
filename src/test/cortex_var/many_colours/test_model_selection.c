/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo 
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */
/*
  test_model_selection.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <CUnit.h>

// cortex_var headers
#include "element.h"
#include "file_reader.h"
#include "dB_graph.h"
#include "dB_graph_population.h"
#include "cmd_line.h"
#include "graph_info.h"
#include "db_differentiation.h"

/*
// get_log_bayesfactor_varmodel_over_repeatmodel is an empty function now
void test_get_log_bayesfactor_varmodel_over_repeatmodel()
{
  if(NUMBER_OF_COLOURS < 100)
  {
    warn("This test requires compilation for support of >=100 colours - skipping\n");
    return;
  }

  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  // Load a single fasta, containing a single bubble. Then run several
  // iterations of a local function that resets all the coverage values on the
  // two branches in many colours to different scenarios.

  // Scenarios I want to test are:
  // 1. All colours have both alleles with equal coverage - is repeat, so
  //    negative bayes factor
  // 2. Allele has 50% frequency, so 50% of people are hets, and 25 of each type
  //    of hom - call as a variant


  void set_coverage_on_bubble(int br1cov, int br2cov, VariantBranchesAndFlanks* var, int colour)
  {
    int i;
    for (i=0; i<var->len_one_allele; i++)
      {
	db_node_set_coverage(var->one_allele[i], colour, br1cov);
      }
    for (i=0; i<var->len_other_allele; i++)
      {
	db_node_set_coverage(var->other_allele[i], colour, br2cov);
      }
  } 


  void clean_up(dBGraph* db_graph)
  {
    void wipe_node(dBNode* node)
    {
      int i;
      for (i=0; i<100; i++)
	{
	  node->coverage[i]=0;
	}
    }
    hash_table_traverse(&wipe_node, db_graph);
  }




  // First set up the hash/graph
  int kmer_size = 13;
  int number_of_bits = 10;
  int bucket_size = 10;
  
  dBGraph *db_graph = hash_table_new(number_of_bits, bucket_size, 10, kmer_size);

  // ===========================================================================
  // Then load a single bubble into the graph
  // ===========================================================================

  // Read sequence
  int fq_quality_cutoff = 20;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_read = 0, seq_loaded = 0;

  load_se_seq_data_into_graph_colour(
    "../data/test/pop_graph/example1_for_testing_genotyping.allele1.fa",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph,
    &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  load_se_seq_data_into_graph_colour(
    "../data/test/pop_graph/example1_for_testing_genotyping.allele2.fa",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph,
    &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  // Now read each allele into an array of nodes, so we can use them

  int file_reader(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;

    if (new_entry!= true){
      die("new_entry has to be true for ead_next_variant_from_full_flank_file");
    }

    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);
    
    return ret;
  }
  
  //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------
  int max_read_length = 30;
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence");
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);
  
  //We are going to load all the bases into a single sliding window 
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      die("Failed to malloc kmer sliding window in "
          "db_graph_make_reference_path_based_sv_calls. Exit.\n");
    }
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-db_graph->kmer_size-1));
  if (kmer_window->kmer==NULL)
    {
      die("Failed to malloc kmer_window->kmer in "
          "db_graph_make_reference_path_based_sv_calls. Exit.\n");
    }
  kmer_window->nkmers=0;
  

  int length_of_arrays=200;
  dBNode**     br1_path        = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays); 
  Orientation* br1_or          = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays); 
  dBNode**     br2_path        = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays); 
  Orientation* br2_or          = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays); 
  
  if ( (br1_path==NULL) || (br2_path==NULL) || (br1_or==NULL) || (br2_or==NULL) )
    {
      die("Unable to alloc array for genotyping test\n");
    }
  
  //end of intialisation 

  FILE* br1_fptr = fopen("../data/test/pop_graph/example1_for_testing_genotyping.allele1.fa", "r");
  FILE* br2_fptr = fopen("../data/test/pop_graph/example1_for_testing_genotyping.allele2.fa", "r");
  if ( (br1_fptr==NULL) || (br2_fptr==NULL) )
    {
      die("Cannot open one of:\n"
"../data/test/pop_graph/example1_for_testing_genotyping.allele1.fa and\n"
"../data/test/pop_graph/example1_for_testing_genotyping.allele2.fa\n");
    }
  
  boolean f_entry=true;
  int br1len = align_next_read_to_graph_and_return_node_array(br1_fptr, max_read_length, br1_path, br1_or,  true, &f_entry, file_reader,
							       seq, kmer_window, db_graph, 0);
  CU_ASSERT(f_entry==true);
  f_entry=true;
  int br2len = align_next_read_to_graph_and_return_node_array(br2_fptr, max_read_length, br2_path, br2_or,  true, &f_entry, file_reader,
							       seq, kmer_window, db_graph, 0);
  CU_ASSERT(f_entry==true);
  fclose(br1_fptr);
  fclose(br2_fptr);
  VariantBranchesAndFlanks* var=alloc_VariantBranchesAndFlanks_object(max_read_length,max_read_length,max_read_length,max_read_length,kmer_size);
  var->one_allele     = br1_path;
  var->len_one_allele = br1len;
  var->other_allele     = br2_path;
  var->len_other_allele = br2len;

  int i;
  GraphInfo* ginfo = graph_info_alloc_and_init();

  for (i=0; i<100; i++)
    {
      int total_seq = 4000;
      int mean_read_len = 24; // these choices explained below
      graph_info_set_seq(ginfo, i, total_seq);
      graph_info_set_mean_readlen(ginfo, i, mean_read_len);
    }

  GraphAndModelInfo model_info;
  long long genome_len=100;
  double mu = 0.8; //param of geometric describing repeat copy num
  //double err = 0.01;
  int ref_colour=-1;//no reference colour
  initialise_model_info(&model_info, ginfo, genome_len, mu, ref_colour, NUMBER_OF_COLOURS*2, EachColourADiploidSample, AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice);
  AnnotatedPutativeVariant annovar;



  clean_up(db_graph);


  //Now we are ready to start

  // 1. Repeat. Every haplotype has two copies. Lets say read length is 24, depth is 40x, k is 13. Then expected depth on a single
  //    allele is 20*(24-13+1)/24 = 10.
  //    Set up so all colours have both alleles with coverage 20

  //   Genome length? Total covg? Fix them so total covg/genome length = 40
  // say genome length=100, sequence covg = 4000;

  for (i=0; i<100; i++)
    {
      set_coverage_on_bubble(20, 20, var, i);
    }
  
  initialise_putative_variant(&annovar, &model_info, var, BubbleCaller, kmer_size, AssumeUncleaned, NULL, NULL, NULL, true);

  
  double ret = get_log_bayesfactor_varmodel_over_repeatmodel(&annovar, &model_info);

  CU_ASSERT(ret< log(100));//called as a repeat if prob var not  100x  greater than porb repeat
  clean_up(db_graph);


  // 2.  Allele has 50% frequency, so 50% of people are hets, and 25 of each type of hom - call as a variant
  //     We tell it read_length was 12, depth was 20, so expected depth at hom site is 10.
  //     Genome length? Total covg? Fix them so total covg/genome length = 20
  //     say genome length=100, sequence covg = 2000;  

  //   Colours 0 to 49 will be het. Give both alleles covg 5 for these
  //   Colours 50 to 74 will be hom_one, covg 10 on first allele and 0 on the other
  //   Colours 75-99 will be hom_other, covg 10 on second allele and 0 on first.

  for (i=0; i<50; i++)
    {
      set_coverage_on_bubble(5,5, var, i);
    }
  for (i=50; i<75; i++)
    {
      set_coverage_on_bubble(10,0, var, i);
    }
  for (i=75; i<100; i++)
    {
      set_coverage_on_bubble(0,10, var, i);
    }


  initialise_putative_variant(&annovar, &model_info, var, BubbleCaller,kmer_size, AssumeUncleaned, NULL, NULL, NULL, true);

  ret = get_log_bayesfactor_varmodel_over_repeatmodel(&annovar, &model_info);

  CU_ASSERT(ret> log(100));//called as a variant

  for (i=0; i<50; i++)
    {
      CU_ASSERT(annovar.genotype[i]==het); 
    }
  for (i=50; i<75; i++)
    {
      CU_ASSERT(annovar.genotype[i]==hom_one); 
    }
  for (i=75; i<100; i++)
    {
      CU_ASSERT(annovar.genotype[i]==hom_other); 
    }
  clean_up(db_graph);







  // 3.  Allele1 has 95% frequency, so 10% of people are hets, 90% of people are hom_one, and 0.2% of people are hom_other
  //     which we treat as zero.

  //   Colours 0 to 89 will be hom_one. covg 10 on first allele and 0 on the other
  //   Colours 90 to 99 will be het, Give both alleles covg 5 for these

  for (i=0; i<90; i++)
    {
      set_coverage_on_bubble(10,0, var, i);
    }
  for (i=90; i<100; i++)
    {
      set_coverage_on_bubble(5,5, var, i);
    }


  initialise_putative_variant(&annovar, &model_info, var, BubbleCaller,kmer_size, AssumeUncleaned, NULL, NULL, NULL, true);
  
  ret = get_log_bayesfactor_varmodel_over_repeatmodel(&annovar, &model_info);
  CU_ASSERT(ret> log(100));//called as a variant
  for (i=0; i<90; i++)
    {
      CU_ASSERT(annovar.genotype[i]==hom_one); 
    }
  for (i=90; i<100; i++)
    {
      CU_ASSERT(annovar.genotype[i]==het); 
    }
  clean_up(db_graph);





  // 4.  Allele1 has 95% frequency, so 10% of people are hets, 90% of people are hom_one, and 0.2% of people are hom_other
  //     which we treat as ONE PERSON <<<< this is the only difference with the last test


  //   Colours 0 to 89 will be hom_one. Give both alleles covg 5 for these
  //   Colours 90 to 98 will be het, covg 10 on first allele and 0 on the other
  //   Colour 99 will be hom_other, given covg 10 on second allele and 0 on first
  for (i=0; i<90; i++)
    {
      set_coverage_on_bubble(10,0, var, i);
    }
  for (i=90; i<99; i++)
    {
      set_coverage_on_bubble(5,5, var, i);
    }
  set_coverage_on_bubble(0,10, var, 99);


  initialise_putative_variant(&annovar, &model_info, var, BubbleCaller,kmer_size, AssumeUncleaned, NULL, NULL, NULL, true);

  
  ret = get_log_bayesfactor_varmodel_over_repeatmodel(&annovar, &model_info);
  CU_ASSERT(ret> log(100) );//called as a variant
  for (i=0; i<90; i++)
    {
      CU_ASSERT(annovar.genotype[i]==hom_one); 
    }
  for (i=90; i<99; i++)
    {
      CU_ASSERT(annovar.genotype[i]==het); 
    }
  CU_ASSERT(annovar.genotype[i]==hom_other); 




  // 5.  Repeat, everybody has some covg of both alleles
  //   

  for (i=0; i<20; i++)
    {
      set_coverage_on_bubble(20,20, var, i);
    }
  for (i=21; i<40; i++)
    {
      set_coverage_on_bubble(16,11, var, i);
    }
  for (i=41; i<60; i++)
    {
      set_coverage_on_bubble(17,18, var, i);
    }
  for (i=61; i<80; i++)
    {
      set_coverage_on_bubble(20,20, var, i);
    }
  for (i=81; i<90; i++)
    {
      set_coverage_on_bubble(12,20, var, i);
    }
  for (i=91; i<100; i++)
    {
      set_coverage_on_bubble(2,5, var, i);
    }


  initialise_putative_variant(&annovar, &model_info, var, BubbleCaller,kmer_size, AssumeUncleaned, NULL, NULL, NULL, true);


  ret = get_log_bayesfactor_varmodel_over_repeatmodel(&annovar, &model_info);
  CU_ASSERT(ret< log(100) );//called as a repeat

  clean_up(db_graph);



  // 6. Suppose we have a repeat that occurs a huge number of times, such as in a reptrotransposon

  for (i=0; i<100; i++)
    {
      set_coverage_on_bubble(2000,2000, var, i);
    }


  initialise_putative_variant(&annovar, &model_info, var, BubbleCaller,kmer_size, AssumeUncleaned, NULL, NULL, NULL, true);


  ret = get_log_bayesfactor_varmodel_over_repeatmodel(&annovar, &model_info);

  CU_ASSERT(ret< log(100) );//called as a repeat





  clean_up(db_graph);



 

  //cleanup
  free_sequence(&seq);
  free_VariantBranchesAndFlanks_object(var);
  free(kmer_window->kmer);
  free(kmer_window);
  free(br1_path);
  free(br2_path);
  free(br1_or);
  free(br2_or);
  hash_table_free(&db_graph);
  graph_info_free(ginfo);
}
*/
