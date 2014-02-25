/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 * 
 * CORTEX project contacts:  
 *    M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 *    Z. Iqbal (zam@well.ox.ac.uk)
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
  test_db_complex_genotyping.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <CUnit.h>
#include <Basic.h>

// cortex_var headers
#include "test_db_complex_genotyping.h"
#include "db_complex_genotyping.h"
#include "element.h"
#include "file_reader.h"
#include "simulator.h"
#include "dB_graph_population.h"


void test_initialise_multiplicities_of_allele_nodes_wrt_both_alleles()
{

  if (NUMBER_OF_COLOURS<6)
    {
      die("Must compile for >=6 colours. recompile\n");
    }
  
  int working_colour1 = 2;
  int working_colour2 = 3;
  

  dBNode* e1 = new_element();
  e1->individual_edges[0]=1;
  dBNode* e2 = new_element();
  e2->individual_edges[0]=1;
  dBNode* e3 = new_element();
  e3->individual_edges[0]=1;
  dBNode* e4 = new_element();
  e4->individual_edges[0]=1;
  dBNode* e5 = new_element();
  e5->individual_edges[0]=1;
  dBNode* e6 = new_element();
  e6->individual_edges[0]=1;
  dBNode* e7 = new_element();
  e7->individual_edges[0]=1;
  dBNode* e8 = new_element();
  e8->individual_edges[0]=1;
  dBNode* e9 = new_element();
  e9->individual_edges[0]=1;
  dBNode* e10 = new_element();
  e10->individual_edges[0]=1;


  //first test - two branches of equal length, each with globally unique nodes
  dBNode* array1[4]={e1,e2,e3,e4};
  dBNode* array2[4]={e5,e6,e7,e8};
  VariantBranchesAndFlanks var;
  set_variant_branches_and_flanks(&var,
				  NULL, NULL, 0,
				  array1, NULL, 4,
				  array2, NULL, 4,
				  NULL, NULL, 0,
				  unknown);

  
  Covg mult11[4]={0,0,0,0};
  Covg mult22[4]={0,0,0,0};
  Covg mult12[4]={0,0,0,0};
  Covg mult21[4]={0,0,0,0};

  // No allele loinger than 10 in this test
  //MultiplicitiesAndOverlapsOfBiallelicVariant* mobv
  //  = alloc_MultiplicitiesAndOverlapsOfBiallelicVariant(10,10);
  //reset_MultiplicitiesAndOverlapsOfBiallelicVariant(mobv);

  MultiplicitiesAndOverlapsOfBiallelicVariant mobv;
  mobv.mult11=mult11;
  mobv.mult22=mult22;
  mobv.mult12=mult12;
  mobv.mult21=mult21;
  mobv.len1=4;
  mobv.len2=4;

  improved_initialise_multiplicities_of_allele_nodes_wrt_both_alleles(
    &var, &mobv,
    //true, &element_get_colour_union_of_all_colours,
    //&element_get_covg_union_of_all_covgs,
    working_colour1, working_colour2);

  int p;
  for (p=0; p<3; p++)
    {
      CU_ASSERT(mobv.mult11[p]==1);
      if (mobv.mult11[p] != 1)
	{
	  warn("Failed for p=%d. mult11[p] =%d\n", p, mobv.mult11[p]);
	}
      CU_ASSERT(mobv.mult22[p]==1);
      CU_ASSERT(mobv.mult12[p]==0);
      CU_ASSERT(mobv.mult21[p]==0);
    }


  //second test - two branches of equal length, one branch containing a repeated node
  array1[1]=e1;//so array1[0]==array1[1]
  reset_MultiplicitiesAndOverlapsOfBiallelicVariant(&mobv);
  improved_initialise_multiplicities_of_allele_nodes_wrt_both_alleles(
    &var, &mobv,
    //true, &element_get_colour_union_of_all_colours,
    //&element_get_covg_union_of_all_covgs,
    working_colour1, working_colour2);

  CU_ASSERT(mobv.mult11[0]==2);
  CU_ASSERT(mobv.mult11[1]==2);
  CU_ASSERT(mobv.mult22[0]==1);
  CU_ASSERT(mobv.mult22[1]==1);
  CU_ASSERT(mobv.mult12[0]==0);
  CU_ASSERT(mobv.mult12[1]==0);
  CU_ASSERT(mobv.mult21[0]==0);
  CU_ASSERT(mobv.mult21[1]==0);

  for (p=2; p<4; p++)
    {
      CU_ASSERT(mobv.mult11[p]==1);
      if (mobv.mult11[p] != 1)
	{
	  warn("Failed for p=%d. mult11[p] =%d\n", p, mobv.mult11[p]);
	}
      CU_ASSERT(mobv.mult22[p]==1);
      CU_ASSERT(mobv.mult12[p]==0);
      CU_ASSERT(mobv.mult21[p]==0);
    }
  

  //third test - two branches, first shorter than second, with some nodes repeated within and between branches
  array1[1]=e2;//reset the second element in array1 to what it was before the last test
  //recall array1[4]={e1,e2,e3,e4}
  dBNode* array3[7]={e4,e1,e9,e10,e9,e9,e2};


  set_variant_branches_and_flanks(&var,
				  NULL, NULL, 0,
				  array1, NULL, 4,
				  array3, NULL, 7,
				  NULL, NULL, 0,
				  unknown);


  Covg mult3_11[4]={0,0,0,0};
  Covg mult3_22[7]={0,0,0,0,0,0,0};
  Covg mult3_12[4]={0,0,0,0};
  Covg mult3_21[7]={0,0,0,0,0,0,0};
  mobv.mult11=mult3_11;
  mobv.mult22=mult3_22;
  mobv.mult12=mult3_12;
  mobv.mult21=mult3_21;
  mobv.len1 = 4;
  mobv.len2 = 7;

  improved_initialise_multiplicities_of_allele_nodes_wrt_both_alleles(
    &var, &mobv,
    //true, &element_get_colour_union_of_all_colours,
    //&element_get_covg_union_of_all_covgs,
    working_colour1, working_colour2);



  CU_ASSERT(mobv.mult11[0]==1);
  CU_ASSERT(mobv.mult11[1]==1);
  CU_ASSERT(mobv.mult11[2]==1);
  CU_ASSERT(mobv.mult11[3]==1);

  CU_ASSERT(mobv.mult22[0]==1);
  CU_ASSERT(mobv.mult22[1]==1);
  CU_ASSERT(mobv.mult22[2]==3);
  CU_ASSERT(mobv.mult22[3]==1);
  CU_ASSERT(mobv.mult22[4]==3);
  CU_ASSERT(mobv.mult22[5]==3);
  CU_ASSERT(mobv.mult22[6]==1);

  CU_ASSERT(mobv.mult12[0]==1);
  CU_ASSERT(mobv.mult12[1]==1);
  CU_ASSERT(mobv.mult12[2]==0);
  CU_ASSERT(mobv.mult12[3]==1);

  CU_ASSERT(mobv.mult21[0]==1);
  CU_ASSERT(mobv.mult21[1]==1);
  CU_ASSERT(mobv.mult21[2]==0);
  CU_ASSERT(mobv.mult21[3]==0);
  CU_ASSERT(mobv.mult21[4]==0);
  CU_ASSERT(mobv.mult21[5]==0);
  CU_ASSERT(mobv.mult21[6]==1);


  // fourth test - two branches, first shorter than second, with some nodes repeated within and between branches
  array1[1]=e2;//reset the second element in array1 to what it was before the last test
  //recall array1[4]={e1,e2,e3,e4}
  dBNode* array4[4]={e1,e2,e3,e4};


  set_variant_branches_and_flanks(&var,
				  NULL, NULL, 0,
				  array1, NULL, 4,
				  array4, NULL, 4,
				  NULL, NULL, 0,
				  unknown);

  
  Covg mult4_11[4]={0,0,0,0};
  Covg mult4_22[7]={0,0,0,0};
  Covg mult4_12[4]={0,0,0,0};
  Covg mult4_21[7]={0,0,0,0};
  mobv.mult11=mult4_11;
  mobv.mult22=mult4_22;
  mobv.mult12=mult4_12;
  mobv.mult21=mult4_21;
  mobv.len1 = 4;
  mobv.len2 = 4;

  improved_initialise_multiplicities_of_allele_nodes_wrt_both_alleles(
    &var, &mobv,
    //true, &element_get_colour_union_of_all_colours,
    //&element_get_covg_union_of_all_covgs,
    working_colour1, working_colour2);



  CU_ASSERT(mobv.mult11[0]==1);
  CU_ASSERT(mobv.mult11[1]==1);
  CU_ASSERT(mobv.mult11[2]==1);
  CU_ASSERT(mobv.mult11[3]==1);

  CU_ASSERT(mobv.mult22[0]==1);
  CU_ASSERT(mobv.mult22[1]==1);
  CU_ASSERT(mobv.mult22[2]==1);
  CU_ASSERT(mobv.mult22[3]==1);

  CU_ASSERT(mobv.mult12[0]==1);
  CU_ASSERT(mobv.mult12[1]==1);
  CU_ASSERT(mobv.mult12[2]==1);
  CU_ASSERT(mobv.mult12[3]==1);

  CU_ASSERT(mobv.mult21[0]==1);
  CU_ASSERT(mobv.mult21[1]==1);
  CU_ASSERT(mobv.mult21[2]==1);
  CU_ASSERT(mobv.mult21[3]==1);


  free(e1);
  free(e2);
  free(e3);
  free(e4);
  free(e5);
  free(e6);
  free(e7);
  free(e8);
  free(e9);
  free(e10);

}

//you give it the fasta files for each allele (expect them to be the 1net fasta for both alleles, or the 2net fasta for both alleles)
void build_and_save_temp_binaries(char* filelist_binaries, 
				  char* first_allele_net_fasta, char* second_allele_net_fasta, char* stub,
				  int kmer, int number_of_bits, int bucket_size)
{
  int max_retries = 10;
  
  dBGraph * db_graph = hash_table_new(number_of_bits, bucket_size,
                                      max_retries, kmer);

  // Read sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour = 0;

  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_read = 0, seq_loaded = 0;

  load_se_seq_data_into_graph_colour(
    first_allele_net_fasta,
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph,
    &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  char bin1[200];
  bin1[0]='\0';
  sprintf(bin1, "../data/tempfiles_can_be_deleted/%s_allele1_temp.ctx", stub);
  db_graph_dump_single_colour_binary_of_colour0(bin1, &db_node_check_status_not_pruned,
                                                db_graph, NULL, BINVERSION);
  db_graph_wipe_colour(0, db_graph);

  load_se_seq_data_into_graph_colour(
    second_allele_net_fasta,
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph,
    &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  char bin2[100];
  bin2[0]='\0';
  sprintf(bin2, "../data/tempfiles_can_be_deleted/%s_allele2_temp.ctx", stub);
  db_graph_dump_single_colour_binary_of_colour0(bin2, &db_node_check_status_not_pruned,
                                                db_graph, NULL, BINVERSION);
  db_graph_wipe_colour(0, db_graph);

  FILE* fp = fopen(filelist_binaries, "w");
  if (fp==NULL)
    {
      die("Cannot open %s\n", filelist_binaries);
    }
  fprintf(fp, "%s\n%s\n", bin1, bin2);
  fclose(fp);


  hash_table_free(&db_graph);
  
}


void utility_func_test_complex_genotyping_given_two_alleles(
  char* first_allele_name, char* second_allele_name,
  char* fasta_allele1, char* fasta_allele2, char* fasta_genome_minus_site, char* fasta_alleles_then_genome,
  char* fasta_allele1_then_allele2, int read_len, int kmer, int genome_size, int number_of_bits, int bucket_size,
  int number_repeats_of_sim, boolean try_multiple_seq_errors,
  boolean using_nets, 
  char* first_allele_1net_fasta, char* first_allele_2net_fasta,
  char* second_allele_1net_fasta, char* second_allele_2net_fasta )
{
  if (NUMBER_OF_COLOURS<6)
    {
      die("Need >=6 colours for test_calc_log_likelihood_of_genotype_with_complex_alleles - recompile\n");
    }
  int colour_allele1 = 0;
  int colour_allele2 = 1;
  int colour_ref_minus_site = 2;
  int colour_indiv = 3;
  int working_colour1 = 4;
  int working_colour2 = 5;
  // These are not used by simulator anymore
  //int working_colour_net1 = 6;
  //int working_colour_net2 = 7;

  //first set up the hash/graph
  int kmer_size = kmer;
  //int number_of_bits = 15;
  //int bucket_size = 100;
  int max_retries = 10;

  if(kmer_size > NUMBER_OF_BITFIELDS_IN_BINARY_KMER*32 ||
     kmer_size < (NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1)*32)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }

  dBGraph * db_graph = hash_table_new(number_of_bits, bucket_size,
                                      max_retries, kmer_size);

  // Load them into colour indiv so it has the right edges. Coverages will be
  // fixed up later, and edges in other colours dont matter for this test

  // Read sequence
  int fq_quality_cutoff = 0;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;

  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_read = 0, seq_loaded = 0;

  load_se_seq_data_into_graph_colour(
    fasta_allele1,
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    colour_indiv, db_graph,
    &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  load_se_seq_data_into_graph_colour(
    fasta_allele2,
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    colour_indiv, db_graph,
    &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  // Now we don't want to load the whole of the rest of the genome - just to
  // annotate these nodes with whether they touch the rest of the genome,
  // create local temporary new hash, load rest of genome, dump binary, then
  // reload those nodes that are ALREADY in our allele1/2 hash

  dBGraph * temp_db_graph = hash_table_new(number_of_bits, bucket_size,
                                           max_retries, kmer_size);

  load_se_seq_data_into_graph_colour(
    fasta_genome_minus_site,
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    0, temp_db_graph,
    &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  GraphInfo* temp_db_graph_info = graph_info_alloc_and_init();

  graph_info_set_seq(temp_db_graph_info, 0, 1);//unnecessary - never used
  graph_info_set_mean_readlen(temp_db_graph_info, 0, 1);//unnecessary - never used
  db_graph_dump_single_colour_binary_of_specified_colour("../data/tempfiles_can_be_deleted/ref_minus_genome.ctx", &db_node_condition_always_true,temp_db_graph,temp_db_graph_info,0, BINVERSION);
  hash_table_free(&temp_db_graph);
  graph_info_free(temp_db_graph_info);

  GraphInfo* bin_ginfo=graph_info_alloc_and_init();
  int clean_colour = 0;
  boolean only_load_kmers_already_in_graph = true;
  load_single_colour_binary_data_from_filename_into_graph("../data/tempfiles_can_be_deleted/ref_minus_genome.ctx", 
							  db_graph, bin_ginfo, false, 
							  colour_ref_minus_site, only_load_kmers_already_in_graph, clean_colour, 
							  false);
  graph_info_free(bin_ginfo);

  int max_allele_length = 90000;

  int depth;
  double seq_err_per_base;

  //for all of these fix read length and kmer
  //read_len=50;
  //kmer = 31;
  //number_repeats_of_sim=100;
  

  //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence");
  }
  alloc_sequence(seq,max_allele_length,MAX_READ_NAME_LEN);
  
  //We are going to load all the bases into a single sliding window 
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      die("Failed to malloc kmer sliding window in db_graph_make_reference_path_based_sv_calls. Exit.\n");
    }
  
  
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_allele_length-db_graph->kmer_size-1));
  if (kmer_window->kmer==NULL)
    {
      die("Failed to malloc kmer_window->kmer in db_graph_make_reference_path_based_sv_calls. Exit.\n");
    }
  kmer_window->nkmers=0;
  
  
  //we are going to need to hold ALL the allele paths in memory at the same time
  int number_paths=3;//actiually only 2 alleles, but also one genome-minus-site will be stored in this
  dBNode*** array_of_node_arrays = (dBNode***) malloc( sizeof(dBNode**) * number_paths );
  Orientation** array_of_or_arrays = (Orientation**) malloc( sizeof(Orientation*) * number_paths );
  int* lengths_of_alleles = (int*) malloc(sizeof(int) * number_paths);
  char** array_of_allele_names = (char**) malloc( sizeof(char*) * number_paths );
  
  if ( (array_of_node_arrays==NULL) || (array_of_or_arrays==NULL) || (lengths_of_alleles==NULL) || (array_of_allele_names==NULL) )
    {
      die("Cannot alloc arrays of arrays in test_db_complex_genotyping.c\n");
    }
  
  int i;
  for (i=0; i<number_paths; i++)
    {
      array_of_node_arrays[i] = (dBNode**) malloc(sizeof(dBNode*) * max_allele_length);
      array_of_or_arrays[i]   = (Orientation*) malloc( sizeof(Orientation) * max_allele_length);
      array_of_allele_names[i]= (char*) malloc(sizeof(char) * 200 );
      
      if ( (array_of_node_arrays[i]==NULL) || (array_of_or_arrays[i]==NULL) || (array_of_allele_names==NULL) )
	{
	  die("Cannot alloc the %d -th node and or array in "
        "print_log_liks_of_specified_set_of_genotypes_of_complex_site", i);
	}
      lengths_of_alleles[i]=0;
      array_of_allele_names[i][0]='\0';
    }



  //  MultiplicitiesAndOverlapsOfBiallelicVariant* mobv = alloc_MultiplicitiesAndOverlapsOfBiallelicVariant(max_allele_length, max_allele_length);
  //if (mobv==NULL)
  //  {
  //    die("Failed to alloc mobv\n");
  //  }
  int* working_array_self = (int*) malloc(sizeof(int) * max_allele_length);
  int* working_array_shared = (int*) malloc(sizeof(int) * max_allele_length);
  
// if ( (mobv==NULL)||(working_array_self==NULL) || (working_array_shared==NULL))
  if ( (working_array_self==NULL) || (working_array_shared==NULL))
    {
      die("Cannot alloc all the arrays in "
          "calculate_max_and_max_but_one_llks_of_specified_set_of_genotypes_of_complex_site."
          "Give up and exit.");
    }
  

  //create file reader
  int file_reader(FILE * fp, Sequence * seq, int max_allele_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;
    if (new_entry == false){
      die("new_entry must be true in hsi test function");
    }
    ret =  read_sequence_from_fasta(fp,seq,max_allele_length,new_entry,full_entry,offset);
    
    return ret;
  }
  
  //end of initialisation 



  FILE* fp = fopen(fasta_alleles_then_genome, "r");
  if (fp==NULL)
    {
      die("Unable to open %s. Exit", fasta_alleles_then_genome);
    }



  int j;
  for (j=0; j<number_paths; j++)//will read allele1, then allele2, then the genome-minus-site into these arrays
    {
      boolean f_entry=true;
      int num_kmers = align_next_read_to_graph_and_return_node_array(fp, max_allele_length, 
								     array_of_node_arrays[j],
								     array_of_or_arrays[j], 
								     false, &f_entry, file_reader, seq, kmer_window, db_graph, -1);
      if (!f_entry)
	{
	  die("Read too long to be read in one go - in test - error\n");
	}
    
      strcat(array_of_allele_names[j], seq->name);
      lengths_of_alleles[j]=num_kmers;
    }
  fclose(fp);
  
  VariantBranchesAndFlanks var;
  int depths[]={20};
  int num_depths=1;

  int p;
  GraphInfo* ginfo=graph_info_alloc_and_init();
  GraphAndModelInfo model_info;

  float repeat_geometric_param_mu = 0.8;//not used in this
  // int genome_size = 554;//554 is length of one allele + rest of reference
  int num_chroms_in_expt=2;
  initialise_model_info(&model_info, ginfo, genome_size, repeat_geometric_param_mu, 
			-1, num_chroms_in_expt, EachColourADiploidSample, AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice);
  zygosity true_genotype;
  for (p=0; p<num_depths; p++)
    {
      depth=depths[p];
      //printf("\n********************************* New set of tests - depth = %d\n", depth);
      printf("First - tests where correct answer is het\n");
      true_genotype=het;

      // set up the variant so it has two different alleles, one on each branch
      set_variant_branches_but_flanks_to_null(&var, 
					      array_of_node_arrays[0], 
					      array_of_or_arrays[0], 
					      lengths_of_alleles[0],
					      array_of_node_arrays[1], 
					      array_of_or_arrays[1],
					      lengths_of_alleles[1],
					      unknown);
      //time_t now;
      //time(&now);
      //printf("%s", ctime(&now));

      graph_info_set_seq(model_info.ginfo, colour_indiv, genome_size*depth);
      graph_info_set_mean_readlen(model_info.ginfo, colour_indiv, read_len);

      char true_gt[100]="";
      strcat(true_gt,first_allele_name);
      strcat(true_gt,"/");
      strcat(true_gt,second_allele_name);
      // run for various sequencing error rates
      //printf("Zam Seq err 0.01\n");
      seq_err_per_base=0.01;
      //model_info.seq_error_rate_per_base=seq_err_per_base;
      model_info.ginfo->seq_err[colour_indiv]=seq_err_per_base; 
      //time(&now);
      //printf("%s", ctime(&now));


      //one binary per allele
      char filelist_1net_binaries[]="../data/tempfiles_can_be_deleted/filelist_hom1_1net.ctxlist";
      char filelist_2net_binaries[]="../data/tempfiles_can_be_deleted/filelist_hom1_2net.ctxlist";



      if (using_nets==true)
	{
	  build_and_save_temp_binaries(filelist_1net_binaries, 
				       first_allele_1net_fasta, second_allele_1net_fasta, "net1", kmer, number_of_bits,  bucket_size);
	  build_and_save_temp_binaries(filelist_2net_binaries, 
				       first_allele_2net_fasta, second_allele_2net_fasta, "net2", kmer, number_of_bits,  bucket_size);
	}



      //these with the [2] are the genome-minus-site  fasta_allele1_then_allele2
      // array_of_node_arrays[2], 
      simulator(depth, read_len, kmer, seq_err_per_base, number_repeats_of_sim, 
                colour_indiv, colour_allele1, colour_allele2,
                colour_ref_minus_site, &var, lengths_of_alleles[2],
                true_genotype, &model_info, fasta_allele1_then_allele2,
                true_gt, working_colour1, working_colour2, using_nets, db_graph);
      //filelist_1net_binaries, filelist_2net_binaries

      if (try_multiple_seq_errors==true)
	{
	  printf("Zam seq err 0.02\n");
	  seq_err_per_base=0.02;
	  //model_info.seq_error_rate_per_base=seq_err_per_base;
	  model_info.ginfo->seq_err[colour_indiv]=seq_err_per_base; 
	  
	  simulator(depth, read_len, kmer, seq_err_per_base, number_repeats_of_sim,
              colour_indiv, colour_allele1, colour_allele2,
              colour_ref_minus_site, &var, lengths_of_alleles[2],
              true_genotype, &model_info, fasta_allele1_then_allele2,
              true_gt, working_colour1, working_colour2, using_nets, db_graph);
	  //filelist_1net_binaries, filelist_2net_binaries

	  printf("Zam seq err 0.001\n");
	  seq_err_per_base=0.001;
	  //	  model_info.seq_error_rate_per_base=seq_err_per_base;
	  model_info.ginfo->seq_err[colour_indiv]=seq_err_per_base; 
	  
	  simulator(depth, read_len, kmer, seq_err_per_base, number_repeats_of_sim,
              colour_indiv, colour_allele1, colour_allele2,
              colour_ref_minus_site, &var, lengths_of_alleles[2],
              true_genotype, &model_info, fasta_allele1_then_allele2,
              true_gt, working_colour1, working_colour2, using_nets, db_graph);
    //filelist_1net_binaries, filelist_2net_binaries
	}
      

      // now re-run for a true hom

      printf("Now for true hom\n");
      true_genotype=hom_one;
      set_variant_branches_but_flanks_to_null(&var, 
					      array_of_node_arrays[0], 
					      array_of_or_arrays[0], 
					      lengths_of_alleles[0],
					      array_of_node_arrays[0], 
					      array_of_or_arrays[0],
					      lengths_of_alleles[0],
					      unknown);


      

      char true_gt2[100]="";
      strcat(true_gt2,first_allele_name);
      strcat(true_gt2,"/");
      strcat(true_gt2,first_allele_name);

      // run for various sequencing error rates
      printf("Zam seq err 0.001\n");
      seq_err_per_base=0.01;
      //model_info.seq_error_rate_per_base=seq_err_per_base;
      model_info.ginfo->seq_err[colour_indiv]=seq_err_per_base; 

      simulator(depth, read_len, kmer, seq_err_per_base, number_repeats_of_sim,
                colour_indiv, colour_allele1, colour_allele2,
                colour_ref_minus_site, &var, lengths_of_alleles[2],
		            true_genotype, &model_info, fasta_allele1_then_allele2,
                true_gt2, working_colour1, working_colour2, using_nets, db_graph);
      //filelist_1net_binaries, filelist_1net_binaries

      if (try_multiple_seq_errors==true)
	{      
	  printf("Zam seq err 0.02\n");
	  seq_err_per_base=0.02;
	  //model_info.seq_error_rate_per_base=seq_err_per_base;
	  model_info.ginfo->seq_err[colour_indiv]=seq_err_per_base; 
	  
	  simulator(depth, read_len, kmer, seq_err_per_base, number_repeats_of_sim,
              colour_indiv, colour_allele1, colour_allele2,
              colour_ref_minus_site, &var, lengths_of_alleles[2],
              true_genotype,  &model_info, fasta_allele1_then_allele2,
              true_gt2, working_colour1, working_colour2, using_nets, db_graph);
    //filelist_1net_binaries, filelist_1net_binaries

	  printf("Zam seq err 0.001\n");
	  seq_err_per_base=0.001;
	  //  model_info.seq_error_rate_per_base=seq_err_per_base;
	  model_info.ginfo->seq_err[colour_indiv]=seq_err_per_base; 
	  
	  simulator(depth, read_len, kmer, seq_err_per_base, number_repeats_of_sim,
              colour_indiv, colour_allele1, colour_allele2,
              colour_ref_minus_site, &var, lengths_of_alleles[2],
              true_genotype,  &model_info, fasta_allele1_then_allele2,
              true_gt2, working_colour1, working_colour2, using_nets, db_graph);
    //filelist_1net_binaries, filelist_1net_binaries
	}

    }


  for (i=0; i<number_paths; i++)
    {
      free(array_of_node_arrays[i]);
      free(array_of_or_arrays[i]);
      free(array_of_allele_names[i]);
    }
  free(working_array_self);
  free(working_array_shared);
  free_sequence(&seq);
  free(kmer_window->kmer);
  free(kmer_window);
  free(array_of_node_arrays);
  free(array_of_or_arrays);
  free(lengths_of_alleles);
  free(array_of_allele_names);
  hash_table_free(&db_graph);
  graph_info_free(ginfo);
}


void test_calc_log_likelihood_of_genotype_with_complex_alleles1()
{

  // load a simple case. 2 alleles which do not overlap at all, or with the rest of the genome
  // load the alleles, and the rest of the genome into the hash, then run the simulator, and check the right genotype
  // happens
  utility_func_test_complex_genotyping_given_two_alleles(
    "allele1", "allele2","../data/test/pop_graph/variations/complex_genotyping/simple_allele1.fa",
    "../data/test/pop_graph/variations/complex_genotyping/simple_allele2.fa",
    "../data/test/pop_graph/variations/complex_genotyping/simple_rest_of_genome.fa",
    "../data/test/pop_graph/variations/complex_genotyping/simple_alleles_then_ref_minus_site.fa",
    "../data/test/pop_graph/variations/complex_genotyping/simple_alleles.fa", 50,31,554,15,100,10, true,
    false, NULL, NULL, NULL, NULL);

}




void test_calc_log_likelihood_of_genotype_with_complex_alleles2()
{

  // load the alleles, and the rest of the genome into the hash, then run the simulator, and check the right genotype
  // happens
  utility_func_test_complex_genotyping_given_two_alleles(
    "hlab_070201_extended", "hlab_550104_extended", "../data/test/pop_graph/variations/complex_genotyping/hlab_070201.fa",
    "../data/test/pop_graph/variations/complex_genotyping/hlab_550104.fa",
    "../data/test/pop_graph/variations/complex_genotyping/chr6_minus_hlab_excerpt.fa",
    "../data/test/pop_graph/variations/complex_genotyping/hlab_two_alleles_then_ref.fa",
    "../data/test/pop_graph/variations/complex_genotyping/both_hlab_alleles.fa", 50,31,89978,15,100,10, false,//changes 24 to 15
    false, NULL, NULL, NULL, NULL);
}


void test_calc_log_likelihood_of_genotype_with_complex_alleles3()
{

  // load a simple case. 2 alleles which do not overlap at all, or with the rest of the genome
  // load the alleles, and the rest of the genome into the hash, then run the simulator, and check the right genotype
  // happens
  utility_func_test_complex_genotyping_given_two_alleles(
    "hlab_070201_extended", "hlab_550104_extended", "../data/test/pop_graph/variations/complex_genotyping/hlab_070201.fa",
    "../data/test/pop_graph/variations/complex_genotyping/hlab_550104.fa",
    "../data/test/pop_graph/variations/complex_genotyping/chr6_minus_hlab_excerpt.fa",
    "../data/test/pop_graph/variations/complex_genotyping/hlab_two_alleles_then_ref.fa",
    "../data/test/pop_graph/variations/complex_genotyping/both_hlab_alleles.fa", 50,31,89978,15,100,10, false,//change 24 to 15
    true, 
    "../data/test/pop_graph/variations/complex_genotyping/hlab_070201_extended.single_errors.fa", 
    "../data/test/pop_graph/variations/complex_genotyping/hlab_070201_extended.double_errors.fa", 
    "../data/test/pop_graph/variations/complex_genotyping/hlab_550104_extended.single_errors.fa", 
    "../data/test/pop_graph/variations/complex_genotyping/hlab_550104_extended.double_errors.fa");
}

/*
This test is based on my genotyping of a bad call - it's not even a bubble, 
i had a bug where the two branches do not rejoin with the right orientation.
The point is this - I want to know that when I genotype a site with the coverges below, 
with the graph_info taken from the chimp binary with which I made the calls - I want to make sure
you (=I) do not genotype lots of colours as het when they have no covg on branch1  at all!!

>var_827009_5p_flank length:1031 average_coverage:17.77 min_coverage:5 max_coverage:40 fst_coverage:15 fst_kmer:CTCTCTTGGTTGCTAATTTGTAGCAGTATGA fst_r:C fst_f:C lst_coverage:25 lst_kmer:CCTGAAGGATACAAACCCTATGATTTATCAT lst_r:A lst_f:GT 
CTCTCTTGGTTGCTAATTTGTAGCAGTATGACTCAAATGCACATAATGATAAAAATGACAAGGAGTTGGGAGTTGCAGTAGGCAAAAGTAGAACTAGGAAAAGTAATCATAGTAAAATTAAATAGGAATTCCAAATAATATAAATGTTTTTTACATTTGTGTTCAAGAAATATAAAATAAGAAAAAGGCTAAACTGTTCTTTCTAGCCTCCAAACACCATCAACCAAATTTTACTTTTATAATTTATATTTGAAAAAATGAAATTCGACCATATCTTCCCCTTAATGAGGAGTGTGTGTGTGTATGTATATCAATACACCTTGTTTTCCCCTTCTTTAGTTTCATGATTCTACCCTCAAATTAAGGTATTATTTTATTTAACAATAAATATTTATTGAGTGCCTGCCATGTTCCAGGTATTGTTCTGGTATTTGACAGACATAGGTAAGCAACACAGGAAAAGATGGCTTCCTTCTTGCACTCTTATCAGAGAAGAAAAAAAAATATCATAAAACCAATAATCATAATATCATAAAAATATTTTATGTTAAAAGGCTATTAATAATATGTAACATTTTAGAATTCTACATGGGCATAATACCATGAGGGAACTTGACCTCTTACTATATTTTATGCTTTATCTTTGAAAAAGTATACATGCCAATTGTTAAGTTCTAGGAAATCATATCACATTAATTCTCTTTGATAGCTGCAAGATTCAGATTGTTCCCATTTCAAGCCTCCTTTGGATTTCTAATTGGCTTGGCTATGCATTAACATGTAGAAAACTGTGGTCATAAATTTGTATTCAGCAAGCATTGGTAGCATGGTTCAGGCAGGTGAAAATGTATAAAGACATAAGAAAATGTGCTATGCAGGCCATTTGGAGAAATACGTCCTTGAAGATTTTTTCAAGTGGTCCCTGTGAATACAATGCATACTTTAATCTCATTTTAGTATAACTTCCTTGTTTAGATTTGTTTCCACCTGTCTGTGAATTCCTGAAGGATACAAACCCTATGATTTATCAT
>branch_827009_1 length:22 average_coverage: 1.38 min_coverage:1 max_coverage:2 fst_coverage:25 fst_kmer:CCTGAAGGATACAAACCCTATGATTTATCAT fst_r:A fst_f:GT lst_coverage:21 lst_kmer:ATTTATCATGGTATACAACATATACACTAAC lst_r:CT lst_f:C 
GGTATACAACATATACACTAAC
>branch_827009_2 length:1346 average_coverage:19.46 min_coverage:7 max_coverage:37 fst_coverage:25 fst_kmer:CCTGAAGGATACAAACCCTATGATTTATCAT fst_r:A fst_f:GT lst_coverage:21 lst_kmer:GTTAGTGTATATGTTGTATACCATGATAAAT lst_r:C lst_f:CT 
TAAACAGGCTGTCAGTTCCAAGACCAAGAGTGGCCAGATTAACTACCAATGAGAAGAAACCCTATCCTAGAGTCTGAACTTATTAATCTGTGGTCTGGGCTTTTGAGCCTTGACCCAATTGATGTCAGAAACCCTAACCAGTAAGACTTGGTCATATGCCTAGAAATGTTTCAAAGATTGCTCTTCCACACCTAGATTATTTCCAATATATTATTGGACAGCTTGTGAAACTAATCCATTTTCACTGACATGACTGGGAAATTCAGTCATAAAAAGTCTACGCAAAGATCATGTAGATATAGTTTTAAAAGGAGTAGGATAAGAAAAATGGCAGAAGGATAAGAAAAATAAGATAGGCTTCGAATGGTCTTTGTAAAGAGGTTCTATAACTATTGTTTAGTTTTTCTTGTATTAATTTCTTCCATTCATTAACAATTCAACAGTGGTTTCTATTTCTTTTTCTGCCTATTTTGTTTTGTTTTGTTTCGTCTTCTCTCTTTCTCTCTGGGATCTTCCTTCTGTTTTTATAAGAGTCTTCATTTTGTCTAAGGAATTGGCTTACAAAGATGTCCTCAATTTCTCCCATTTCATTCTCTAGTATACAGCTTGCCTTCCCAACCTTTATCACTGAAATTTGTCCGAGGTGGAAAAATAATACTTTCCTTATGGCATAAGATGAACAAGTTCTCTGACCTAAAGACTGTTGTAGTCAAGAATGTCTTTACATGCTTGGTCTTGCGAGTTCTGAGAATGCTGTCATTGAGTATGCTTGGCTGATACTCGCCCCTAACTACCATGCTTGTTGCTAACTGCTGGTACTCCACTAGGTGTACCATTGTGCCAATATGTATCTTTGTTCCAATATGTATCTTTAATTTACTAAGACTTTATTTCTGTTTAGAATGATTAAATTATAGATTTTGATTTTAGACCAATATTTTAAATTTATAAAGATAAAAAGTTTTAGGAAGCTCTTGCTGCCATTTTTCTTATCTTTCTCTTTTTAAAACTATGTCTATATGATCTTTGGCGTAGGCTTTTTATGTTTGAATTTCCCAATCATGTCAGTAAAAATGGATTAGTCTCACAAGATATCCAGTAATGTATTTGAAATAACCTACATGTGGAAGAACATTCTTTGAAATATTTCTAGGACTATGAGCAAGGCTTGCTGGTTAGGGTTTCTGATATCAGTTGGGGCAATTCTCAAAAGCCCAGACCATAGATCAATAAGCTCAGACTGTAGTATGGATCTAAAAGGCAGAGTGTATTGCCAGTTTTGGCTTCTGGTGGAAAATCCAAGCTGATAACTGAGGTTAGTGTATATGTTGTATACCATGATAAAT
>var_827009_3p_flank length:0 average_coverage:19.65 min_coverage:0 max_coverage:0 fst_coverage:21 fst_kmer:GTTAGTGTATATGTTGTATACCATGATAAAT fst_r:C fst_f:CT lst_coverage:21 lst_kmer:GTTAGTGTATATGTTGTATACCATGATAAAT lst_r:C lst_f:CT 


branch1 coverages
Covg in Colour 0:
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
Covg in Colour 1:
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
Covg in Colour 2:
6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
Covg in Colour 3:
3 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 
Covg in Colour 4:
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
Covg in Colour 5:
4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
Covg in Colour 6:
2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
Covg in Colour 7:
3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
Covg in Colour 8:
2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
Covg in Colour 9:
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
Covg in Colour 10:
4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
branch2 coverages
Covg in Colour 0:
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
Covg in Colour 1:
0 0 0 0 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 4 3 3 3 3 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 2 4 4 4 4 4 4 4 4 4 4 4 3 3 3 3 3 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 3 3 5 5 5 4 4 4 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 5 5 5 5 5 5 5 5 5 3 3 3 3 3 3 3 4 4 4 4 4 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 1 1 1 2 3 3 3 4 4 4 4 5 5 5 5 5 4 4 4 5 5 5 5 5 4 3 3 3 2 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 3 3 4 4 4 4 4 4 3 3 3 3 3 3 3 3 3 3 3 2 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 2 2 3 3 3 3 3 3 3 3 2 2 2 3 3 3 2 2 3 3 3 3 3 2 2 2 2 2 2 2 2 3 3 3 2 2 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 2 2 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 0 0 0 0 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 0 0 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 3 3 3 4 4 4 4 5 5 5 4 4 4 3 3 3 3 3 3 3 3 2 2 2 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 2 3 3 4 5 6 6 6 6 8 8 8 8 7 7 6 6 6 6 6 6 5 4 4 3 3 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 4 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 3 3 3 3 3 1 1 1 1 1 1 2 2 2 2 1 1 1 1 1 1 2 2 3 4 4 4 4 4 4 4 4 3 3 3 3 3 3 3 3 3 3 2 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 1 1 2 2 2 2 3 3 3 3 3 3 3 4 6 7 7 7 6 7 7 7 7 6 6 6 6 5 5 5 5 5 5 5 4 2 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 1 1 1 2 2 2 2 2 2 2 2 2 3 3 3 3 2 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 2 2 2 3 3 3 3 3 4 3 3 3 3 3 3 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
Covg in Colour 2:
6 6 6 6 4 3 3 3 2 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 2 3 3 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 1 2 2 2 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 1 1 2 2 3 3 3 3 3 3 3 4 3 3 3 3 3 3 3 3 3 3 3 2 3 2 2 2 2 2 2 2 1 1 1 1 1 1 1 2 2 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 0 0 0 0 0 0 0 0 0 0 0 1 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 4 2 2 2 1 1 1 1 1 1 1 1 3 3 3 4 4 4 4 4 3 3 3 4 5 5 5 5 5 5 5 5 5 3 3 3 2 2 2 2 2 2 2 2 1 0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 4 4 4 4 3 3 3 2 2 2 2 2 2 2 2 2 2 1 1 2 2 1 1 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 1 1 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 2 2 2 2 2 2 3 5 5 6 6 6 6 6 6 5 5 5 5 5 5 4 4 5 5 5 5 4 2 2 1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 1 2 2 2 2 2 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 3 3 3 3 3 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 3 3 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 3 3 4 3 3 3 3 3 3 2 2 2 2 3 3 4 4 4 5 5 5 4 4 3 3 3 3 3 3 3 4 5 5 4 4 4 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 1 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 2 2 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 5 5 4 4 4 3 3 3 3 3 3 3 3 2 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 2 2 2 2 2 1 1 2 2 3 3 3 3 3 3 3 3 3 3 3 3 4 4 5 5 5 5 5 4 4 3 3 3 3 3 3 3 4 4 4 4 4 2 2 1 1 1 1 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5 5 5 5 4 4 4 3 3 3 3 4 3 3 3 3 4 3 3 3 3 3 3 3 3 3 3 4 3 3 3 3 2 2 2 2 2 1 1 2 2 2 2 2 3 3 3 3 2 2 2 2 2 2 2 3 3 3 3 3 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 1 1 1 2 2 3 3 3 4 4 4 4 4 4 4 4 4 3 3 3 3 3 3 4 4 4 3 4 4 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 3 2 3 3 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 1 1 1 1 1 1 1 2 2 4 4 4 4 5 5 5 5 4 4 4 4 4 4 4 4 4 5 5 4 4 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 3 5 4 4 4 5 5 5 6 6 6 6 6 7 8 8 10 10 10 10 9 8 6 6 6 6 5 5 5 4 4 4 4 4 3 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 
Covg in Colour 3:
3 2 2 2 2 2 1 2 2 2 2 2 1 1 1 1 0 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 0 0 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 2 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 2 2 2 2 2 1 2 3 3 3 4 4 4 4 4 4 4 4 4 6 6 5 6 6 6 6 6 5 4 4 4 3 3 3 3 3 3 3 3 3 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 3 3 1 1 2 2 2 2 2 2 2 2 2 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 1 1 2 2 2 2 3 2 2 2 2 2 2 2 2 2 2 2 2 3 3 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 2 2 2 2 2 3 3 3 2 2 2 2 2 3 3 3 3 3 3 3 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 2 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 3 3 3 4 3 1 2 2 2 2 2 2 3 3 3 4 3 3 3 3 3 3 3 3 2 2 2 1 1 1 1 2 2 1 1 1 1 1 2 2 2 2 3 3 3 3 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 0 0 0 1 1 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 2 2 1 1 1 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 2 2 2 2 3 4 4 4 4 4 4 4 4 4 4 4 4 3 3 3 3 2 2 2 2 1 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 3 3 3 3 3 3 2 2 2 3 3 3 3 3 3 3 3 3 3 4 4 3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 3 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 2 2 2 3 3 3 2 2 2 2 3 3 3 3 3 3 3 3 3 3 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 2 2 3 3 3 3 4 4 4 3 3 3 4 4 4 5 4 5 5 5 5 5 5 4 3 3 3 2 2 2 2 2 2 2 2 2 1 1 0 0 0 2 2 2 2 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 4 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 2 2 3 2 3 3 4 4 4 4 5 5 5 5 5 6 6 6 6 6 6 6 6 5 5 4 5 4 5 5 5 4 4 4 4 4 3 3 3 3 3 3 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 4 4 4 4 5 5 5 5 5 5 5 5 5 4 4 4 4 4 4 4 4 1 1 1 1 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3 3 3 3 
Covg in Colour 4:
1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 2 2 2 2 3 3 3 3 3 3 4 4 4 4 3 3 3 3 3 3 3 2 2 2 2 1 1 1 1 1 2 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 2 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 6 6 6 5 5 4 3 2 2 2 2 2 2 2 2 3 3 3 3 3 3 2 2 2 2 2 2 3 4 4 4 4 4 4 4 5 4 4 4 4 4 4 3 4 4 4 4 4 3 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3 5 5 4 4 5 5 5 5 5 5 5 5 5 3 3 3 3 4 4 4 4 2 2 3 3 2 2 2 2 2 2 2 2 2 3 3 3 3 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 3 1 1 1 1 1 1 1 1 1 1 1 2 2 2 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 3 3 4 3 3 3 2 2 2 2 1 1 2 2 3 3 3 3 3 3 3 3 3 2 2 3 3 4 4 4 4 4 4 3 3 3 3 3 3 3 3 5 5 5 5 5 4 4 4 4 4 4 4 4 4 4 3 3 3 3 3 3 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 3 3 3 3 2 2 2 2 2 3 4 4 4 4 4 4 4 4 4 4 5 4 4 4 4 4 4 4 4 4 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 1 1 0 0 0 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 2 2 2 2 2 2 2 3 3 1 1 1 1 2 3 3 2 2 2 2 2 2 2 2 2 2 3 4 4 4 4 4 4 5 5 4 4 4 4 5 5 6 6 6 6 6 6 5 4 4 4 4 4 4 3 2 2 2 2 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 2 2 1 1 1 1 1 1 1 1 2 2 2 3 3 4 4 4 4 4 4 3 3 3 3 3 4 4 4 4 4 3 4 4 3 3 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 4 4 4 4 3 3 2 2 2 2 2 3 2 2 2 2 2 2 2 2 3 2 2 2 2 2 2 4 4 4 4 4 3 3 3 3 3 3 3 3 4 3 3 3 3 3 3 3 1 1 1 1 3 4 4 4 4 4 4 4 4 3 3 3 3 3 3 3 3 5 5 5 5 4 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 2 4 4 4 3 3 3 3 3 3 3 3 3 3 2 4 5 5 5 5 5 5 3 3 3 3 3 3 3 3 4 4 5 7 7 7 5 4 4 4 4 4 4 5 5 5 5 5 5 4 4 4 4 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 4 5 5 6 6 6 6 6 5 5 5 5 5 5 5 5 5 4 4 6 6 4 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 1 1 2 2 2 2 2 3 3 3 3 4 4 4 4 4 4 4 4 4 4 3 3 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 2 2 2 2 2 2 2 2 2 3 3 3 3 2 2 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 1 1 1 2 2 2 1 1 1 1 1 2 2 2 2 2 2 
Covg in Colour 5:
4 4 4 4 4 4 4 4 4 4 4 4 2 3 4 3 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 1 2 2 2 3 3 4 4 4 5 5 5 5 5 6 6 6 7 8 8 8 7 6 6 6 5 5 4 4 4 3 4 4 4 4 3 3 3 2 1 1 1 1 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 2 2 2 2 2 2 3 4 4 4 4 4 4 4 4 4 4 4 3 3 4 4 4 4 4 4 4 3 2 2 2 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 4 4 4 4 3 3 2 2 3 3 3 3 2 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 2 2 2 2 3 3 3 3 3 3 3 4 4 4 4 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 2 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 2 2 2 2 2 2 2 2 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 2 2 2 2 2 3 3 2 2 2 2 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 3 3 3 2 2 2 2 3 4 4 4 5 5 5 5 5 5 5 5 5 3 3 6 6 6 5 5 6 6 5 6 6 5 5 5 5 5 5 5 5 5 5 5 3 3 3 3 3 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 3 4 4 5 5 5 5 6 6 6 6 6 6 6 6 5 5 5 5 5 5 4 3 3 2 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 1 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 3 3 3 2 2 2 2 2 2 2 2 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 1 2 3 3 4 5 5 5 5 4 4 4 4 4 5 5 5 5 5 5 5 5 4 3 3 2 1 1 1 2 3 3 3 2 2 3 3 3 3 3 3 6 6 6 6 6 6 7 7 7 6 5 5 5 6 6 5 5 5 5 5 5 2 2 2 2 2 2 1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 2 2 2 3 4 4 4 4 4 4 4 4 5 4 4 4 4 4 4 4 3 3 3 3 4 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 2 2 1 1 1 1 1 1 1 0 0 0 0 0 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 
Covg in Colour 6:
2 2 2 2 2 2 1 1 1 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 3 3 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 3 3 3 3 2 1 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 2 2 2 2 3 3 2 2 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 2 2 2 2 2 2 3 3 3 3 4 4 4 4 4 4 4 4 4 3 3 2 2 2 2 2 2 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 2 2 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 2 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 3 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 3 3 4 4 3 3 3 3 4 4 4 4 4 4 4 4 4 4 5 5 3 3 3 2 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 3 3 2 2 2 2 2 2 2 2 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 2 2 2 2 2 2 2 2 2 2 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 6 6 6 5 4 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 2 2 2 2 1 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 5 5 5 6 6 6 6 6 7 7 6 6 6 6 6 6 6 6 6 6 6 4 4 4 3 3 3 2 2 1 1 1 2 2 2 2 2 4 4 4 4 3 3 3 3 3 3 3 3 3 3 3 3 2 3 3 3 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 2 2 2 2 2 2 3 3 4 4 4 4 4 5 5 5 5 5 4 4 4 4 4 4 4 4 4 3 4 3 3 3 3 3 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 3 3 3 4 5 5 5 5 4 4 3 3 3 3 3 3 3 3 3 3 3 1 1 1 0 1 1 1 1 1 1 1 1 1 2 2 3 3 3 3 3 3 4 4 4 4 3 3 3 3 3 4 4 4 4 3 3 2 1 1 1 1 1 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 2 2 3 3 2 2 2 2 2 3 3 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 2 2 2 1 2 2 3 3 3 3 3 4 4 4 4 4 4 4 4 3 3 2 2 2 2 
Covg in Colour 7:
3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 3 3 3 4 4 4 4 4 4 4 4 4 4 4 6 6 6 5 4 5 5 4 4 5 4 4 4 3 3 3 3 3 3 3 3 2 2 2 2 2 1 1 1 1 1 1 1 1 2 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 3 3 3 3 3 2 1 2 2 2 2 2 2 2 2 2 2 2 2 3 3 2 2 2 2 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 3 3 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 2 2 2 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 2 2 2 2 1 1 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 3 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 3 3 2 2 2 2 2 2 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 2 3 2 2 2 2 1 1 1 1 1 2 3 3 3 3 3 4 4 4 4 4 3 4 4 4 4 4 4 5 5 5 4 3 4 4 4 5 4 4 4 4 4 4 3 3 3 3 3 4 3 3 3 4 4 3 3 3 2 2 2 2 2 2 2 2 2 2 3 3 2 2 2 2 1 1 1 1 1 2 4 4 4 4 4 4 4 4 5 4 5 5 6 6 6 6 6 6 6 6 6 4 4 5 5 5 5 5 5 4 5 4 4 3 3 3 3 4 4 4 4 3 3 3 2 2 3 3 3 3 3 2 2 2 3 3 3 4 3 3 3 3 3 3 3 3 3 2 2 1 1 2 2 3 4 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 3 3 3 3 2 2 2 2 2 2 2 2 1 2 2 2 2 2 4 4 4 3 3 3 3 4 4 4 4 4 4 6 6 6 5 4 4 4 4 2 2 2 2 2 2 2 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 2 2 2 2 2 3 4 5 5 5 4 4 4 5 5 6 6 6 6 6 6 6 6 6 6 6 5 4 3 3 3 3 4 4 4 4 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 3 3 2 2 2 2 1 1 1 1 2 2 2 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 3 3 3 4 4 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 1 0 1 3 3 3 3 3 4 4 4 4 5 5 5 5 5 5 6 6 6 6 6 5 3 3 3 3 3 3 3 3 3 2 2 2 2 2 3 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 4 5 5 5 5 4 4 4 4 4 4 5 5 5 5 5 5 3 3 3 3 2 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 2 2 4 4 4 4 5 5 5 5 5 5 5 5 5 4 4 4 4 4 4 4 4 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 3 3 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 2 2 2 2 2 2 2 2 3 4 4 4 4 4 4 4 4 4 4 3 3 2 2 2 2 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 2 3 3 3 3 3 3 3 4 3 3 4 4 4 4 4 4 4 4 4 4 3 2 2 3 3 3 3 3 2 2 2 1 1 1 1 1 2 2 2 2 2 2 4 4 3 5 5 5 5 5 6 6 7 7 7 7 9 8 8 8 8 8 8 6 6 6 5 5 5 6 6 5 5 4 4 4 4 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 6 6 7 7 7 7 7 7 7 7 6 8 7 8 8 7 7 6 6 6 6 3 3 3 3 3 3 3 3 3 3 3 1 1 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 1 1 1 1 2 2 2 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 3 3 3 3 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 2 2 2 2 2 2 3 3 3 4 4 4 4 4 4 4 3 3 3 3 3 2 2 2 2 2 2 1 1 1 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 1 1 1 2 2 2 2 2 2 3 2 2 2 2 2 1 1 1 1 1 1 2 2 2 1 1 
Covg in Colour 8:
2 4 4 3 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 0 0 0 0 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 4 4 4 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 4 4 3 3 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 2 3 3 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 2 2 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 1 1 1 2 2 2 2 2 4 4 4 3 3 3 3 3 3 3 6 6 7 7 7 7 7 7 8 8 8 4 4 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 3 3 3 4 4 4 4 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 3 3 4 4 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 5 5 5 5 6 6 6 4 4 4 5 5 6 6 5 5 5 5 7 8 8 7 7 7 7 7 7 7 6 7 7 6 6 6 6 6 8 8 8 6 5 5 5 5 5 5 3 3 3 3 3 3 3 5 5 5 5 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 4 4 3 2 2 2 2 2 2 2 2 2 2 3 3 2 2 2 2 2 2 2 2 2 2 2 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 3 3 3 3 3 3 2 2 2 2 2 2 2 3 4 4 4 4 4 4 4 2 2 2 2 2 3 2 2 2 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 2 2 3 3 3 3 3 5 5 8 8 8 8 8 8 8 9 9 9 8 7 7 7 7 7 7 7 7 6 6 4 4 4 5 5 4 4 1 1 1 4 4 4 4 4 4 4 4 4 4 3 3 3 3 3 3 3 3 3 3 3 0 0 0 1 1 1 1 2 2 2 2 2 4 4 4 4 5 5 5 5 5 5 5 5 4 4 4 4 3 3 3 3 3 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 2 4 4 4 4 4 3 3 3 3 3 3 3 4 4 4 4 4 3 3 3 3 1 1 1 1 1 1 1 1 1 2 2 3 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 4 4 3 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 0 0 0 0 0 0 0 0 0 0 0 1 1 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 5 6 2 2 2 2 2 2 3 3 3 4 4 4 4 4 4 4 3 3 3 3 2 2 3 3 3 3 3 2 2 3 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 5 5 3 3 6 6 6 5 5 5 5 5 5 5 5 4 4 4 4 4 4 3 3 3 3 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 2 2 2 3 3 3 3 3 5 5 
Covg in Colour 9:
0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 1 1 1 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 2 2 2 2 2 1 1 1 2 2 2 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2 2 2 3 2 2 2 2 2 2 2 2 2 2 1 1 1 2 2 2 2 3 3 3 2 2 2 2 2 3 3 3 3 3 3 5 5 7 6 6 6 6 5 5 5 5 5 5 4 5 5 5 5 5 5 5 3 3 2 2 2 2 2 2 2 2 2 2 2 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 1 1 2 2 2 2 2 3 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 1 1 1 1 2 3 3 3 3 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 3 4 4 5 5 5 5 5 5 5 5 3 3 3 3 3 3 3 3 3 3 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 0 2 2 2 3 3 3 3 3 3 2 2 2 2 4 4 4 5 5 5 5 5 3 3 4 4 4 4 4 4 4 4 4 4 4 2 2 2 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 2 2 2 2 2 2 2 2 1 1 1 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 1 1 1 2 2 2 2 3 2 3 3 3 3 3 3 3 3 3 4 4 4 4 4 5 4 4 4 4 3 3 2 2 2 2 2 2 2 2 2 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 2 2 2 3 3 3 3 4 4 4 4 4 5 5 5 5 5 4 4 4 5 4 5 4 4 4 4 4 3 3 3 3 3 2 2 2 2 2 2 2 2 1 1 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 3 4 4 4 4 3 4 4 4 5 4 4 4 4 4 4 7 7 8 8 8 7 6 6 6 6 6 5 5 5 4 4 4 4 4 4 4 1 1 0 0 0 0 0 0 0 0 0 1 1 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 1 1 1 1 2 2 2 3 3 3 2 2 2 2 2 3 3 3 3 3 3 2 2 2 3 3 3 3 2 2 2 3 3 4 4 4 3 3 4 4 4 4 4 4 4 3 3 3 3 3 3 3 2 2 1 1 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 2 2 1 1 1 2 3 4 3 3 2 3 3 3 3 3 3 3 3 3 4 5 5 5 5 5 4 4 3 3 3 3 2 2 2 3 3 3 3 3 4 3 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 4 4 4 4 3 3 3 3 3 3 3 2 2 2 2 1 1 1 2 2 2 1 1 1 1 2 3 3 3 3 4 5 5 5 5 4 4 4 4 5 5 5 5 5 5 5 4 3 3 3 3 2 1 1 1 1 1 1 1 1 0 0 0 1 1 1 1 1 2 2 2 2 2 2 2 2 3 4 5 5 5 5 5 5 4 4 4 4 5 3 3 4 4 4 4 4 4 3 2 1 1 1 1 1 1 1 1 1 2 2 2 2 1 1 1 1 1 1 1 2 3 3 3 4 4 4 4 5 5 4 4 4 4 4 5 5 5 5 5 5 4 3 2 3 2 2 2 3 5 5 5 5 5 5 6 5 5 5 5 5 6 6 6 6 5 5 6 6 5 3 3 3 3 3 3 2 2 3 3 3 3 2 2 2 2 2 2 1 1 1 1 1 2 2 3 3 4 4 3 4 4 4 4 4 4 4 4 4 5 5 5 5 5 4 4 3 3 2 2 2 1 1 1 1 1 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 3 3 2 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 4 4 4 4 4 4 3 3 4 4 4 3 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 2 2 2 3 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 3 2 4 4 3 3 3 3 
Covg in Colour 10:
4 4 4 3 3 3 4 4 4 3 3 3 3 4 4 4 4 4 3 3 3 3 3 3 3 3 3 2 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3 4 4 4 4 4 4 3 3 3 3 3 3 3 1 2 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 3 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 2 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 4 4 4 4 3 3 3 3 4 4 3 3 3 3 3 3 3 2 2 2 2 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 3 3 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 3 3 3 3 3 3 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 3 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 2 3 3 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 3 4 5 5 5 4 4 4 3 2 3 3 3 3 3 3 3 3 3 3 3 3 2 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 1 0 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 1 1 1 1 1 0 0 0 0 0 0 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 6 5 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 4 4 4 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 2 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 3 3 3 2 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 4 3 3 3 4 4 4 4 5 5 6 5 5 5 5 6 7 7 6 6 6 5 5 5 5 4 4 4 5 4 4 3 3 3 3 3 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 3 4 4 5 6 6 6 6 6 7 7 6 6 6 6 5 5 5 5 5 5 5 4 4 3 2 2 2 2 2 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 3 3 4 4 4 3 3 3 3 3 3 3 2 2 2 2 2 2 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 4 4 4 4 4 5 4 4 4 4 3 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 3 4 4 4 3 3 3 3 3 3 3 3 3 4 4 4 3 4 4 5 6 5 5 5 5 5 6 6 7 7 7 7 7 7 6 6 6 6 5 5 4 3 3 2 2 2 2 1 1 0 0 0 0 0 0 1 1 1 2 2 2 2 2 2 2 3 4 4 4 4 4 4 4 4 4 4 3 3 3 2 2 2 2 2 2 2 1 0 0 0 0 0 0 0 1 2 2 2 2 3 3 3 3 3 4 4 4 3 3 3 3 3 4 4 4 3 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 1 1 1 2 2 2 2 2 2 2 


 */

void regression_test_1_single_bubble_call_one_allele_shorter_than_k_one_very_long()
{
  if (NUMBER_OF_COLOURS<11)
  {
    warn("This test requires >=11 colours, skipping\n");
    return;
  }

  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1)
  {
    warn("Test not configured for NUMBER_OF_BITFIELDS_IN_BINARY_KMER > 1\n");
    return;
  }


  VariantBranchesAndFlanks var;
  dBNode* branch1[22];
  dBNode* branch2[134];
  int i,j;
  for (i=0; i<22; i++)
    {
      branch1[i]=new_element();
    }
  for (i=0; i<134; i++)
    {
      branch2[i]=new_element();
    }
  //int widths[] = { [0 ... 9] = 1, [10 ... 99] = 2, [100] = 3 };
  int br1_covg_0[]={[0]=1, [1 ... 21]=0};
  for (j=0; j<22; j++)
    {
      branch1[j]->coverage[0]=br1_covg_0[j];
    }
  int br1_covg_1[]={[0 ... 21]=0};
  for (j=0; j<22; j++)
    {
      branch1[j]->coverage[1]=br1_covg_1[j];
    }
  int br1_covg_2[]={[0]=6, [1 ... 21]=0};
  for (j=0; j<22; j++)
    {
      branch1[j]->coverage[2]=br1_covg_2[j];
    }

  int br1_covg_3[]={3,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,1};
  for (j=0; j<22; j++)
    {
      branch1[j]->coverage[3]=br1_covg_3[j];
    }

  int br1_covg_4[]={[0]=1, [1 ... 21]=0};
  for (j=0; j<22; j++)
    {
      branch1[j]->coverage[4]=br1_covg_4[j];
    }


  int br1_covg_5[]={[0]=4, [1 ... 21]=0};
  for (j=0; j<22; j++)
    {
      branch1[j]->coverage[5]=br1_covg_5[j];
    }

  int br1_covg_6[]={[0]=2, [1 ... 21]=0};
  for (j=0; j<22; j++)
    {
      branch1[j]->coverage[6]=br1_covg_6[j];
    }

  int br1_covg_7[]={[0]=3, [1 ... 21]=0};
  for (j=0; j<22; j++)
    {
      branch1[j]->coverage[7]=br1_covg_7[j];
    }

  int br1_covg_8[]={[0]=2, [1 ... 21]=0};
  for (j=0; j<22; j++)
    {
      branch1[j]->coverage[8]=br1_covg_8[j];
    }
  //colour 9 all zeroes, no need to do anything
  
  int br1_covg_10[]={[0]=4, [1 ... 21]=0};
  for (j=0; j<22; j++)
    {
      branch1[j]->coverage[10]=br1_covg_10[j];
    }


  //branch2 coverages
  int br2_covg_0[]={[0 ... 133]=1};
  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[0]=br2_covg_0[j];
    }

  int br2_covg_1[]={0,0,0,0,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,4,3,3,3,3,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,2,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4};
  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[1]=br2_covg_1[j];
    }


  int br2_covg_2[]={6,6,6,6,4,3,3,3,2,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1};
  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[2]=br2_covg_2[j];
    }


  int br2_covg_3[]={3,2,2,2,2,2,1,2,2,2,2,2,1,1,1,1,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,0,0,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[3]=br2_covg_3[j];
    }



  int br2_covg_4[]={1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,2,2,2,2,3,3,3,3,3,3,4,4,4,4,3,3,3,3,3,3,3,2,2,2,2,1,1,1,1,1,2,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1};

  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[4]=br2_covg_4[j];
    }


  int br2_covg_5[]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3,2,2,2,2,3,4,4,4,5,5,5,5,5,5,5,5,5,3,3,6,6,6,5,5,6,6,5,6,6,5,5,5,5,5,5,5,5,5,5,5,3,3,3,3,3,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0};

  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[5]=br2_covg_5[j];
    }


  int br2_covg_6[]={2,2,2,2,2,2,1,1,1,2,2,2,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[6]=br2_covg_6[j];
    }



  int br2_covg_7[]={3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,3,3,3,4,4,4,4,4,4,4,4,4,4,4,6,6,6,5,4,5,5,4,4,5,4,4,4,3,3,3,3,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1,1,1,2,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,3,3,3,3,3,2,1,2,2,2,2,2,2,2,2,2,2,2,2,3,3,2,2,2,2,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2};

  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[7]=br2_covg_7[j];
    }


  int br2_covg_8[]={2,4,4,3,3,3,3,3,3,3,3,3,3,3,3,1,1,1,1,1,1,1,0,0,0,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,4,4,4,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,4,4,3,3,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,2,2,2};


  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[8]=br2_covg_8[j];
    }


  int br2_covg_9[]={0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,1,1,1,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,2,2,2,2,2,1,1,1,2,2,2,2,2,2,2,2,1,1,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2};

  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[9]=br2_covg_9[j];
    }



  int br2_covg_10[]={4,4,4,3,3,3,4,4,4,3,3,3,3,4,4,4,4,4,3,3,3,3,3,3,3,3,3,2,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,4,4,4,4,4,4,3,3,3,3,3,3,3,1,2,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,2,1,1,1,1,1,2,2,2,3,2,2,2,2,2,2,4,4};


  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[10]=br2_covg_10[j];
    }


  /*
    output from making the chimp calls:

SUMMARY:
Colour:	MeanReadLen	TotalSeq
0	0	2912760135
1	50	18286122352
2	50	16816361244
3	50	18039181209
4	50	15879192506
5	50	17729089947
6	50	15750659112
7	50	26196361173
8	52	20202087523
9	50	18907785783
10	50	16870486574
****************************************

   */

  GraphInfo* ginfo = graph_info_alloc_and_init();

  graph_info_set_seq(ginfo, 0, 2912760135);
  graph_info_set_seq(ginfo, 1, 18286122352 );
  graph_info_set_seq(ginfo, 2, 16816361244);
  graph_info_set_seq(ginfo, 3, 18039181209);
  graph_info_set_seq(ginfo, 4, 15879192506);
  graph_info_set_seq(ginfo, 5, 17729089947);
  graph_info_set_seq(ginfo, 6, 15750659112);
  graph_info_set_seq(ginfo, 7, 26196361173);
  graph_info_set_seq(ginfo, 8, 20202087523);
  graph_info_set_seq(ginfo, 9, 18907785783);
  graph_info_set_seq(ginfo, 10, 16870486574);
  graph_info_set_mean_readlen(ginfo, 0, 0);
  graph_info_set_mean_readlen(ginfo, 1, 50);
  graph_info_set_mean_readlen(ginfo, 2, 50);
  graph_info_set_mean_readlen(ginfo, 3, 50);
  graph_info_set_mean_readlen(ginfo, 4, 50);
  graph_info_set_mean_readlen(ginfo, 5, 50);
  graph_info_set_mean_readlen(ginfo, 6, 50);
  graph_info_set_mean_readlen(ginfo, 7, 50);
  graph_info_set_mean_readlen(ginfo, 8, 52);
  graph_info_set_mean_readlen(ginfo, 9, 50);
  graph_info_set_mean_readlen(ginfo, 10, 50);
 

  GraphAndModelInfo model_info;
  double mu=0.8;
  //double seq_err_rate_per_base=0.01;
  int ref_colour=0;
  int num_chroms=20; 
  long long genome_len = 3000000000;
  initialise_model_info(&model_info, ginfo, genome_len, mu, //seq_err_rate_per_base, 
			ref_colour, num_chroms, EachColourADiploidSampleExceptTheRefColour, AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice);

  var.one_allele       = branch1;
  var.len_one_allele   = 22;
  var.other_allele     = branch2;
  var.len_other_allele = 134;

  CovgArray* working_ca = alloc_and_init_covg_array(200);

  //should not make a difference if you measure coverage using median or not, should still give same gt
  boolean count_covg_using_jumps[2];
  count_covg_using_jumps[0]=false;
  count_covg_using_jumps[1]=true;
  int c;
  for (c=0; c<2; c++)
    {
      AnnotatedPutativeVariant annovar;
      initialise_putative_variant(&annovar, &model_info, 
				  &var, BubbleCaller, 31, 
				  AssumeUncleaned,
				  working_ca, true, count_covg_using_jumps[c] );
      
      
      // Since none of the colours except colour 3 has any coverage AT ALL on branch1, I simply
      // cannot accept a genotype call which is het or hom_one for those colours
      
      CU_ASSERT(annovar.genotype[1]!=hom_one);
      CU_ASSERT(annovar.genotype[2]!=hom_one);
      //leaving out colour 3
      CU_ASSERT(annovar.genotype[4]!=hom_one);
      CU_ASSERT(annovar.genotype[5]!=hom_one);
      CU_ASSERT(annovar.genotype[6]!=hom_one);
      CU_ASSERT(annovar.genotype[7]!=hom_one);
      CU_ASSERT(annovar.genotype[8]!=hom_one);
      CU_ASSERT(annovar.genotype[9]!=hom_one);
      CU_ASSERT(annovar.genotype[10]!=hom_one);
      
      //hom one is least likely
      CU_ASSERT(annovar.gen_log_lh[10].log_lh[0]<annovar.gen_log_lh[10].log_lh[1]);
      CU_ASSERT(annovar.gen_log_lh[10].log_lh[0]<annovar.gen_log_lh[10].log_lh[2]);
      
      CU_ASSERT(annovar.gen_log_lh[1].log_lh[0]<annovar.gen_log_lh[1].log_lh[1]);
      CU_ASSERT(annovar.gen_log_lh[1].log_lh[0]<annovar.gen_log_lh[1].log_lh[2]);
      
      CU_ASSERT(annovar.gen_log_lh[2].log_lh[0]<annovar.gen_log_lh[2].log_lh[1]);
      CU_ASSERT(annovar.gen_log_lh[2].log_lh[0]<annovar.gen_log_lh[2].log_lh[2]);
      
      CU_ASSERT(annovar.gen_log_lh[3].log_lh[0]<annovar.gen_log_lh[3].log_lh[1]);
      CU_ASSERT(annovar.gen_log_lh[3].log_lh[0]<annovar.gen_log_lh[3].log_lh[2]);
      
      CU_ASSERT(annovar.gen_log_lh[4].log_lh[0]<annovar.gen_log_lh[4].log_lh[1]);
      CU_ASSERT(annovar.gen_log_lh[4].log_lh[0]<annovar.gen_log_lh[4].log_lh[2]);
      
      CU_ASSERT(annovar.gen_log_lh[5].log_lh[0]<annovar.gen_log_lh[5].log_lh[1]);
      CU_ASSERT(annovar.gen_log_lh[5].log_lh[0]<annovar.gen_log_lh[5].log_lh[2]);
      
      CU_ASSERT(annovar.gen_log_lh[6].log_lh[0]<annovar.gen_log_lh[6].log_lh[1]);
      CU_ASSERT(annovar.gen_log_lh[6].log_lh[0]<annovar.gen_log_lh[6].log_lh[2]);
      
      CU_ASSERT(annovar.gen_log_lh[7].log_lh[0]<annovar.gen_log_lh[7].log_lh[1]);
      CU_ASSERT(annovar.gen_log_lh[7].log_lh[0]<annovar.gen_log_lh[7].log_lh[2]);
      
      CU_ASSERT(annovar.gen_log_lh[8].log_lh[0]<annovar.gen_log_lh[8].log_lh[1]);
      CU_ASSERT(annovar.gen_log_lh[8].log_lh[0]<annovar.gen_log_lh[8].log_lh[2]);
      
      CU_ASSERT(annovar.gen_log_lh[9].log_lh[0]<annovar.gen_log_lh[9].log_lh[1]);
      CU_ASSERT(annovar.gen_log_lh[9].log_lh[0]<annovar.gen_log_lh[9].log_lh[2]);
      
    }
      
      
  //cleanup
  free_covg_array(working_ca);

  for (i=0; i<22; i++)
    {
      free(branch1[i]);
    }
  for (i=0; i<134; i++)
    {
      free(branch2[i]);
    }
  graph_info_free(ginfo);
}



void regression_test_2_genotyping_of_PD_SNP_call()
{
  if(NUMBER_OF_BITFIELDS_IN_BINARY_KMER != 2)
  {
    warn("Null test - compile for k=55\n");
    return;
  }

  int kmer_size = 55;

  if(NUMBER_OF_COLOURS < 2)
  {
    warn("This test requires >= 2 colours, skipping\n");
    return;
  }

  int file_reader_fasta(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;
    if (new_entry == false){
      //offset = kmer_size;
      die("new_entry must be true in hsi test function");
    }
    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);
    
    return ret;
  }



  //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------

  int max_read_length=400;

  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence");
  }
  alloc_sequence(seq,max_read_length,max_read_length);
  
  //We are going to load all the bases into a single sliding window 
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      die("Failed to malloc kmer sliding window in "
          "align_list_of_fastaq_to_graph_and_print_coverages_in_all_colours. "
          "Exit.\n");
    }
  

  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-kmer_size-1));
  if (kmer_window->kmer==NULL)
    {
      die("Failed to malloc kmer_window->kmer in "
          "align_list_of_fastaq_to_graph_and_print_coverages_in_all_colours. "
          "Exit.\n");
    }
  kmer_window->nkmers=0;
  
  
  //end of intialisation 


  /*

One of the variants called on chr22 on NA12878

>var_18_5p_flank length:165 average_coverage: 5.00 min_coverage:3 max_coverage:11 mode_coverage: 4 percent_nodes_with_modal_covg: 36.00 percent_novel:  0.00 fst_coverage:11 fst_kmer:AATACTATTTTTAGGCCAGGCATGGTGGCTCATACCTGTAATCCCAACATTTTGG fst_r:A fst_f:G lst_coverage:6 lst_kmer:TAGTGAGACCTTTCCTTTATTAAATAAATAAATAGGATGGGCACTGTGGCTCATA lst_r:T lst_f:T 
AATACTATTTTTAGGCCAGGCATGGTGGCTCATACCTGTAATCCCAACATTTTGGGAGGCCAAGTTTGGAGACTCATTTGAGTCCAGGAGTTGACCAGCCTGGGCAACATAGTGAGACCTTTCCTTTATTAAATAAATAAATAGGATGGGCACTGTGGCTCATAT
>var_18_trusted_branch length:56 average_coverage: 1.00 min_coverage:0 max_coverage:19 mode_coverage: 0 percent_nodes_with_modal_covg: 92.00 percent_novel:  0.00 fst_coverage:0 fst_kmer:GTGAGACCTTTCCTTTATTAAATAAATAAATAGGATGGGCACTGTGGCTCATATC fst_r: fst_f: lst_coverage:8 lst_kmer:TGTAATCCCAGCATTTTGGGTTGCCAAGGCAGGAGGCTTGCTTGAGCCCAGGAGT lst_r:C lst_f:T covgs of trusted not variant nodes:  0 0 0 0 0 0 0 0 16 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 19 0 0 0 0 0 0 0 0 16 0 0 0 0 0 number_of_such_nodes: 55
CTGTAATCCCAGCATTTTGGGTTGCCAAGGCAGGAGGCTTGCTTGAGCCCAGGAGT
>var_18_variant_branch length:56 average_coverage: 6.00 min_coverage:4 max_coverage:9 mode_coverage: 7 percent_nodes_with_modal_covg: 46.00 percent_novel: 98.00 fst_coverage:6 fst_kmer:AGTGAGACCTTTCCTTTATTAAATAAATAAATAGGATGGGCACTGTGGCTCATAT fst_r:A fst_f:G lst_coverage:8 lst_kmer:TGTAATCCCAGCATTTTGGGTTGCCAAGGCAGGAGGCTTGCTTGAGCCCAGGAGT lst_r:C lst_f:T covgs of variant not trusted nodes:  7 7 5 8 5 8 8 7 8 8 7 7 7 5 7 4 7 6 8 7 9 7 9 7 7 7 7 7 6 7 8 7 6 9 7 7 7 7 8 5 8 7 4 7 8 6 7 7 7 8 6 8 6 5 8 number_of_such_nodes: 55
GTGTAATCCCAGCATTTTGGGTTGCCAAGGCAGGAGGCTTGCTTGAGCCCAGGAGT
>var_18_3p_flank length:45 average_coverage: 4.00 min_coverage:1 max_coverage:8 mode_coverage: 2 percent_nodes_with_modal_covg: 35.00 percent_novel:  0.00 fst_coverage:8 fst_kmer:TGTAATCCCAGCATTTTGGGTTGCCAAGGCAGGAGGCTTGCTTGAGCCCAGGAGT fst_r:C fst_f:T lst_coverage:1 lst_kmer:GCCCAGGAGTTTGAGACCAGCCTGGACAGCATAGCAAGACTCCATCTCTACAAAT lst_r:T lst_f: 
TTGAGACCAGCCTGGACAGCATAGCAAGACTCCATCTCTACAAAT
var_18 - extra information

branch1 coverages
Mult in  hum ref:
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
Covg in indiv:
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16 16 19 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 
branch2 coverages
Mult in  hum ref:
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
Covg in indiv:
6 6 6 7 7 7 7 7 7 7 8 8 9 9 9 8 8 8 8 8 8 7 5 6 6 6 5 5 4 4 5 5 6 7 8 8 8 7 7 7 7 7 7 7 7 7 7 7 7 8 7 7 7 7 7 8 
  */




  // First set up the hash/graph

  int number_of_bits = 18;
  int bucket_size = 100;
  int max_retries = 10;

  dBGraph * db_graph = hash_table_new(number_of_bits, bucket_size,
                                      max_retries, kmer_size);

  // I load this into the graph so the kmers are there, but then I am going to
  // just create a var object with the covgs I want

  // Read sequence
  int fq_quality_cutoff = 20;
  int homopolymer_cutoff = 0;
  boolean remove_duplicates_se = false;
  char ascii_fq_offset = 33;
  int into_colour;

  unsigned long long bad_reads = 0, dup_reads = 0;
  unsigned long long seq_read = 0, seq_loaded = 0;

  // Load into colour 0
  into_colour = 0;
  load_se_seq_data_into_graph_colour(
    "../data/test/pop_graph/variations/complex_genotyping/pd_example1_branch1andflanks.fa",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  // Load into colour 1
  into_colour = 1;
  load_se_seq_data_into_graph_colour(
    "../data/test/pop_graph/variations/complex_genotyping/pd_example1_branch2andflanks.fa",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  // For second var/test

  // Load into colour 0
  into_colour = 0;
  load_se_seq_data_into_graph_colour(
    "../data/test/pop_graph/variations/complex_genotyping/pd_example2_branch1_and_flanks.fa",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  // Load into colour 1
  into_colour = 1;
  load_se_seq_data_into_graph_colour(
    "../data/test/pop_graph/variations/complex_genotyping/pd_example2_branch2_and_flanks.fa",
    fq_quality_cutoff, homopolymer_cutoff,
    remove_duplicates_se, ascii_fq_offset,
    into_colour, db_graph, &bad_reads, &dup_reads, &seq_read, &seq_loaded,
    NULL, 0, &subsample_null);

  VariantBranchesAndFlanks var;
  dBNode** branch1 = (dBNode**) malloc(sizeof(dBNode*) * max_read_length);
  dBNode** branch2 = (dBNode**) malloc(sizeof(dBNode*) * max_read_length);
  Orientation* branch1_o = (Orientation*) malloc(sizeof(Orientation)*max_read_length);
  Orientation* branch2_o = (Orientation*) malloc(sizeof(Orientation)*max_read_length);

  if (  (branch1==NULL) || (branch2==NULL) || (branch1_o==NULL) || (branch2_o==NULL) )
    {
      die("Unable to malloc branch arrays in test\n");
    }


  FILE* fp = fopen("../data/test/pop_graph/variations/complex_genotyping/pd_example1_just_alleles.fa", "r");
  if (fp==NULL)
    {
      die("Unable to open test file ../data/test/pop_graph/variations/complex_genotyping/pd_example1_just_alleles.fa\n");
    }
  boolean full_entry=true;
  align_next_read_to_graph_and_return_node_array(fp, max_read_length, branch1, branch1_o, false, &full_entry, &file_reader_fasta, 
								      seq, kmer_window, db_graph, 0);

  align_next_read_to_graph_and_return_node_array(fp, max_read_length, branch2, branch2_o, false, &full_entry, &file_reader_fasta, 
								      seq, kmer_window, db_graph, 0);

  fclose(fp);
  int j;

  int br1_covg_0[]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

  for (j=0; j<56; j++)
    {
      branch1[j]->coverage[0]=br1_covg_0[j];
    }
  int br1_covg_1[]= {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,16,19,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8};

  for (j=0; j<56; j++)
    {
      branch1[j]->coverage[1]=br1_covg_1[j];
    }



  int br2_covg_0[]={1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  for (j=0; j<56; j++)
    {
      branch2[j]->coverage[0]=br2_covg_0[j];
    }
  int br2_covg_1[]= {6,6,6,7,7,7,7,7,7,7,8,8,9,9,9,8,8,8,8,8,8,7,5,6,6,6,5,5,4,4,5,5,6,7,8,8,8,7,7,7,7,7,7,7,7,7,7,7,7,8,7,7,7,7,7,8 };

  for (j=0; j<56; j++)
    {
      branch2[j]->coverage[1]=br2_covg_1[j];
    }
  
  /*

The above site was called by the PD caller on NA12878 with k=55 - this is the graph info:
Colour 0 = reference


****************************************
  SUMMARY:
  Colour: MeanReadLen     TotalSeq
0       100     200
1       90      73500000000
****************************************
*/

  GraphInfo* ginfo=graph_info_alloc_and_init();

  graph_info_set_seq(ginfo, 0, 200);
  graph_info_set_seq(ginfo, 1, 73500000000 );
  graph_info_set_mean_readlen(ginfo, 0, 100);
  graph_info_set_mean_readlen(ginfo, 1, 90);

  for (j=2; j<NUMBER_OF_COLOURS; j++)
    { 
      graph_info_set_seq(ginfo, j, 0 );
      graph_info_set_mean_readlen(ginfo, j, 0);
    }

  GraphAndModelInfo model_info;
  double mu=0.8;
  // double seq_err_rate_per_base=0.01;
  int ref_colour=0;
  int num_chroms=2; 
  long long genome_len = 3000000000;
  initialise_model_info(&model_info, ginfo, genome_len, mu, //seq_err_rate_per_base, 
			ref_colour, num_chroms, EachColourADiploidSampleExceptTheRefColour, AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice);

  var.one_allele       = branch1;
  var.one_allele_or    =branch1_o;
  var.len_one_allele   = 56;
  var.other_allele     = branch2;
  var.other_allele_or    =branch2_o;
  var.len_other_allele = 56;
  var.which = first;
  CovgArray* working_ca = alloc_and_init_covg_array(200);

  AnnotatedPutativeVariant annovar;


  //counting covg by counting jumps
  initialise_putative_variant(&annovar, &model_info, 
			      &var, SimplePathDivergenceCaller, 31,
			      AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice, 
			      working_ca, true, false );
  CU_ASSERT(annovar.genotype[1]==hom_other);

  //using median
  initialise_putative_variant(&annovar, &model_info, 
			      &var, SimplePathDivergenceCaller, 31,
			      AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice,
			      working_ca, true, true );
  CU_ASSERT(annovar.genotype[1]==hom_other);



  //******************
  /// second test

  fp = fopen("../data/test/pop_graph/variations/complex_genotyping/pd_example2_just_alleles.fa", "r");
  if (fp==NULL)
    {
      die("Unable to open test file ../data/test/pop_graph/variations/complex_genotyping/pd_example2_just_alleles.fa\n");
    }
  align_next_read_to_graph_and_return_node_array(fp, max_read_length, branch1, branch1_o, false, &full_entry,&file_reader_fasta, 
								      seq, kmer_window, db_graph, 0);

  align_next_read_to_graph_and_return_node_array(fp, max_read_length, branch2, branch2_o, false, &full_entry,&file_reader_fasta, 
								      seq, kmer_window, db_graph, 0);


  fclose(fp);
  //second meaning second test
  int second_br1_covg_0[]={1,1,1,1,1,1,1,1,1,1,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

  for (j=0; j<56; j++)
    {
      branch1[j]->coverage[0]=second_br1_covg_0[j];
    }

  //the 6 on the end basically is there to make sure I ignore the star/end nodes! Makes a big difference
  int second_br1_covg_1[]= {10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6};

  for (j=0; j<56; j++)
    {
      branch1[j]->coverage[1]=second_br1_covg_1[j];
    }



  //  int second_br2_covg_0[]={1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  // for DEBUG
  int second_br2_covg_0[]={1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10};

  for (j=0; j<56; j++)
    {
      branch2[j]->coverage[0]=second_br2_covg_0[j];
    }
  int second_br2_covg_1[]= {10,10,11,11,11,12,12,12,11,11,11,11,11,10,11,11,11,11,10,10,10,10,10,10,10,10,10,10,9,8,7,6,5,5,5,5,6,6,5,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,6,6};

  for (j=0; j<56; j++)
    {
      branch2[j]->coverage[1]=second_br2_covg_1[j];
    }
  


  AnnotatedPutativeVariant annovar2;


  initialise_putative_variant(&annovar2, &model_info, 
			      &var, SimplePathDivergenceCaller, 31,
			      AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice, 
			      working_ca, true, false );
  CU_ASSERT(annovar2.genotype[1]==hom_other);//1==colour


  //using median
  initialise_putative_variant(&annovar2, &model_info, 
			      &var, SimplePathDivergenceCaller, 31,
			      AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice, 
			      working_ca, true, true );
  CU_ASSERT(annovar2.genotype[1]==hom_other);//1==colour








  hash_table_free(&db_graph);
  free(kmer_window->kmer);
  free(kmer_window);
  free_sequence(&seq);
  free(branch1);
  free(branch2);
  free(branch1_o);
  free(branch2_o);
  graph_info_free(ginfo);



}
