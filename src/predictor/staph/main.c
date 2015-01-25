/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 *
  main.c 
*/

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <string_buffer.h>

// Mykrobe headers
#include "element.h"
#include "file_reader.h"
#include "build.h"
#include "cmd_line.h"
#include "graph_info.h"
#include "db_differentiation.h"
#include "maths.h"
#include "gene_presence.h"
#include "genotyping_known.h"
#include "antibiotics.h"
#include "species.h"
#include "json.h"

void timestamp();

int main(int argc, char **argv)
{
  setvbuf(stdout, NULL, _IOLBF, 0);




  CmdLine* cmd_line = cmd_line_alloc();
  if (cmd_line==NULL)
    {
      return -1;
    }
    // VERSION_STR is passed from the makefile -- usually last commit hash

  parse_cmdline(cmd_line, argc,argv,sizeof(Element));

  if (cmd_line->format==Stdout){
    printf("myKrobe.predictor for Staphylococcus, version %d.%d.%d.%d"VERSION_STR"\n",
           VERSION, SUBVERSION, SUBSUBVERSION, SUBSUBSUBVERSION);  
  }

  dBGraph * db_graph = NULL;



  boolean (*subsample_function)();
  
  //local func
  boolean subsample_as_specified()
  {
    double ran = drand48();
    if (ran <= cmd_line->subsample_propn)
      {
	return true;
      }
    return false;
  }
  //end of local func
  
  if (cmd_line->subsample==true)
    {
      subsample_function = &subsample_as_specified;
    }
  else
    {
      subsample_function = &subsample_null;
    }
  





  //  int lim = cmd_line->max_expected_sup_len;
    /* CovgArray* working_ca_for_median=alloc_and_init_covg_array(lim);//will die if fails to alloc
  if (working_ca_for_median==NULL)
    {
      return -1;
    }*/

  //Create the de Bruijn graph/hash table
  int max_retries=15;
  db_graph = hash_table_new(cmd_line->mem_height,
			    cmd_line->mem_width,
			    max_retries, 
			    cmd_line->kmer_size);
  if (db_graph==NULL)
    {
      return -1;
    }

  //some setup
  int file_reader_fasta(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;
    if (new_entry == false){
      offset = db_graph->kmer_size;
      //die("new_entry must be true in hsi test function");
    }
    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);
    return ret;
  }

  ReadingUtils* ru = alloc_reading_utils(MAX_LEN_GENE, db_graph->kmer_size);
  
  //init the dbNodes, mostly to keep valgrind happy
  int i;
  BinaryKmer b;
  binary_kmer_initialise_to_zero(&b);
  dBNode dummy_node;
  element_initialise(&dummy_node, &b, cmd_line->kmer_size);
  for (i=0; i<MAX_LEN_GENE; i++)
    {
      ru->array_nodes[i]=&dummy_node;
    }


  VarOnBackground* tmp_vob = alloc_and_init_var_on_background();
  GeneInfo* tmp_gi = alloc_and_init_gene_info();
  AntibioticInfo* abi = alloc_antibiotic_info();
  // Alloc an array to store called variants and gene presence
  CalledVariant* called_variants = alloc_and_init_called_variant_array();
  CalledGene* called_genes = alloc_and_init_called_genes_array();

  if ( (ru==NULL) || (tmp_vob==NULL) || (abi==NULL) || (tmp_gi==NULL) || (called_variants==NULL) || (called_genes==NULL) )

    {
      return -1;
    }
  

  if (cmd_line->format==Stdout)
    {
      printf("** Start time\n");
      timestamp();
      printf("** Sample:\n%s\n", cmd_line->id->buff);
    }
 
  int into_colour=0;
  boolean only_load_pre_existing_kmers=false;
  uint64_t bp_loaded=0;

  if (cmd_line->method==WGAssemblyThenGenotyping)
    {
      only_load_pre_existing_kmers=false;
    }
  else if (cmd_line->method==InSilicoOligos)
    {
      StrBuf* sk = strbuf_new();
      strbuf_append_str(sk, cmd_line->install_dir->buff);
      strbuf_append_str(sk, "data/skeleton_binary/staph/skeleton.k15.ctx");
      if (access(sk->buff,F_OK)!=0)
	{
	  StrBuf* skeleton_flist = strbuf_new();
	  strbuf_append_str(skeleton_flist, 
			    cmd_line->install_dir->buff);
	  strbuf_append_str(skeleton_flist, 
			    "data/skeleton_binary/staph/list_speciesbranches_genes_and_muts");
	  uint64_t dummy=0;
	  boolean is_rem=true;
	  timestamp();
	  build_unclean_graph(db_graph, 
			      skeleton_flist,
			      true,
			      cmd_line->kmer_size,
			      NULL, 0,
			      NULL, 0,
			      false,
			      into_colour,
			      &subsample_null,
			      false, &dummy, 0, &is_rem);

	  //dump binary so can reuse

	  db_graph_dump_binary(sk->buff, 
			       &db_node_condition_always_true,
			       db_graph,
			       NULL,
			       BINVERSION);
	  strbuf_free(skeleton_flist);
	  set_all_coverages_to_zero(db_graph, 0);
	}
      else
	{
	  if (cmd_line->format==Stdout)
	    {
	      timestamp();
	    }
	  int num=0;
	  GraphInfo* ginfo=graph_info_alloc_and_init();//will exit it fails to alloc.
	  load_multicolour_binary_from_filename_into_graph(sk->buff, db_graph, ginfo,&num);
	  graph_info_free(ginfo);
	  set_all_coverages_to_zero(db_graph, 0);
    if (cmd_line->format==Stdout){
    }
	}
      strbuf_free(sk);
      
      only_load_pre_existing_kmers=true;
    }
  else
    {
      die("For now --method only allowed to take InSilicoOligos or WGAssemblyThenGenotyping\n");
    }


  //only need this for progress
  uint64_t total_reads = 0;
  uint64_t count_so_far=0;
  if (cmd_line->progress==true)
    {
      //timestamp();
      total_reads=count_all_reads(cmd_line->seq_path, cmd_line->input_list);
      //timestamp();
    }
  boolean progressbar_remainder=true;
  if (cmd_line->format==Stdout){
    timestamp();
  }

  bp_loaded = build_unclean_graph(db_graph, 
				  cmd_line->seq_path,
				  cmd_line->input_list,
				  cmd_line->kmer_size,
				  cmd_line->readlen_distrib,
				  cmd_line->readlen_distrib_size,
				  cmd_line->kmer_covg_array, 
				  cmd_line->len_kmer_covg_array,
				  only_load_pre_existing_kmers,
				  into_colour, subsample_function,
				  cmd_line->progress, &count_so_far, total_reads,
				  &progressbar_remainder);
  if (  (cmd_line->progress==true) && (progressbar_remainder==true) )
    {
      printf("Progress %" PRIu64 "/%" PRIu64 "\n", total_reads, total_reads);
    }
  fflush(stdout);
  if (bp_loaded==0)
    {
      printf("No data\n");
      return 1;

    }
  if (cmd_line->format==Stdout){
    timestamp();
  }
  
  unsigned long mean_read_length = calculate_mean_uint64_t(cmd_line->readlen_distrib,
							   cmd_line->readlen_distrib_size);

  double err_rate = estimate_err_rate(cmd_line->seq_path, cmd_line->input_list);
  if (err_rate<0.01)
    {
      err_rate=0.01;
    }

  // if minor_pop_frac is not specified then set the default based on the error rate. 
  if (cmd_line->min_frac_to_detect_minor_pops == 0.0){
    if (err_rate<0.02)
    {
      cmd_line->min_frac_to_detect_minor_pops = 0.1;
    }
  else
    {
      cmd_line->min_frac_to_detect_minor_pops = 0.25;
    }
  }

  //given the error rate and other params, we can estimate expected depth of covg, and read-arrival rate
  // lambda_g = Depth/read_len _g means lambda on the true genome
  double lambda_g_err_free = ((double) bp_loaded/(double)(cmd_line->genome_size)) / (double) mean_read_length ; 

  int expected_depth 
    = (int) ( pow(1-err_rate, cmd_line->kmer_size)  
	      * (mean_read_length-cmd_line->kmer_size+1)
	      * lambda_g_err_free );
  //printf("Expected depth %d\n", expected_depth);
  //  clean_graph(db_graph, cmd_line->kmer_covg_array, cmd_line->len_kmer_covg_array,
  //   expected_depth, cmd_line->max_expected_sup_len); 
  
  
  //calculate expected read-arrival rates on true and error alleles
  double lambda_g_err = 
    lambda_g_err_free
    * pow(1-err_rate, cmd_line->kmer_size) ;

  //rate of arrival of reads on a sequencing error allele
  double lambda_e_err = 
    lambda_g_err_free
    * err_rate/3
    * pow(1-err_rate, cmd_line->kmer_size-1);
  
  StrBuf* tmp_name = strbuf_new();
  // SampleModel* species_mod = alloc_and_init_sample_model();
  SpeciesInfo* species_info = get_species_info(db_graph, 10000, cmd_line->install_dir,
                                  expected_depth,1,1);

  if (species_info->sample_type == MixedStaph)
    {
      strbuf_append_str(tmp_name, "Staphylococcus Mixed ");
      // strbuf_append_str(tmp_name, species_mod->name_of_non_aureus_species->buff);
    }
  // else if (species_info->sample_type == MinorStaphAureusAndMajorNonCoag)
  //   {
  //     strbuf_append_str(tmp_name, species_mod->name_of_non_aureus_species->buff);
  //   }
  else if (species_info->sample_type == NonStaphylococcal)
    {
      strbuf_append_str(tmp_name, "Non-staphylococcal");
    }
  else
    {
      strbuf_append_str(tmp_name, "S.aureus");
    }

  if (cmd_line->format==Stdout)
    {
      printf("** Species\n");
      if (species_info->sample_type != PureStaph)
	{
	  printf("%s\n No AMR predictions given.\n** End time\n", tmp_name->buff);
	  timestamp();
	  // free_sample_model(species_mod);

	  //cleanup
	  strbuf_free(tmp_name);
	  free_antibiotic_info(abi);
	  free_var_on_background(tmp_vob);
	  free_gene_info(tmp_gi);
	  free_reading_utils(ru);
	  
	  cmd_line_free(cmd_line);
	  hash_table_free(&db_graph);
	  
	  return 0;
	}
      else
	{
	  printf("%s\n", tmp_name->buff);
	  timestamp();
	  // free_sample_model(species_mod);

	  printf("** Antimicrobial susceptibility predictions\n");
	}
    }
  else//JSON
    {
      print_json_start();
      print_json_called_variant_item("expected_depth",expected_depth,false);
      print_json_called_variant_item("mean_read_length",mean_read_length,false);
      print_json_species_start();
      if (species_info->sample_type == PureStaph)
	    {
	     print_json_called_variant_item(get_pure_species_name(species_info), get_pure_species_coverage(species_info), true);
	   }
      else if (species_info->sample_type == MixedStaph) 
  	{
      int i;
      boolean last = false;
      for (i=0; i < species_info->num_species; i++)
      {
        if (i == species_info->num_species-1){
          last = true;
        }
        print_json_called_variant_item(get_ith_species_name(species_info,i), get_ith_species_coverage(species_info,i), last);
      }
  	}
    else
  	{
  	  print_json_item("Non-staphylococcal", "Major", true);
  	}
      print_json_species_end(); 

      if (species_info->sample_type == NonStaphylococcal)
	{
	  print_json_susceptibility_start(); 
	  print_json_susceptibility_end();
	  print_json_virulence_start();
	  print_json_virulence_end();
	  print_json_end();

	  //cleanup
	  strbuf_free(tmp_name);
	  free_antibiotic_info(abi);
	  free_var_on_background(tmp_vob);
	  free_gene_info(tmp_gi);
	  free_reading_utils(ru);
	  
	  cmd_line_free(cmd_line);
	  hash_table_free(&db_graph);
	  
	  return 0;
	}
    }
  
  //assumption is num_bases_around_mut_in_fasta is at least 30, to support all k<=31.
  //if k=31, we want to ignore 1 kmer at start and end
  //if k=29, we want to ignore 3 kmers at start and end.. etc


  //if we get here, is s. aureus

  if (cmd_line->format==JSON)
    {
      print_json_susceptibility_start();
    }
  int ignore = cmd_line->num_bases_around_mut_in_fasta - cmd_line->kmer_size +2;  
  boolean output_last=false;
  
  
  boolean any_erm_present=false;
  print_erythromycin_susceptibility(db_graph, &file_reader_fasta, ru, tmp_vob, tmp_gi, abi,
            &is_erythromycin_susceptible, tmp_name, cmd_line->install_dir,
            ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate, cmd_line, output_last,     
            &any_erm_present,
            called_variants,called_genes);
  
  
  
  
  
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_vob, tmp_gi, abi,
          &is_gentamicin_susceptible, tmp_name, cmd_line->install_dir,
          ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate, cmd_line, output_last,
                      called_variants,called_genes);   
  
  
  
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_vob, tmp_gi, abi,
          &is_biocides_susceptible, tmp_name, cmd_line->install_dir,
          ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate, cmd_line, output_last,
                      called_variants,called_genes);   
  
  
  
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_vob, tmp_gi, abi,
          &is_chloramphenicol_susceptible, tmp_name, cmd_line->install_dir,
          ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate, cmd_line, output_last,
                      called_variants,called_genes);   
  
  
  
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_vob, tmp_gi, abi,
          &is_linezolid_susceptible, tmp_name, cmd_line->install_dir,
          ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate, cmd_line, output_last,
                      called_variants,called_genes);   
  
  
  
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_vob, tmp_gi, abi,
          &is_penicillin_susceptible, tmp_name, cmd_line->install_dir,
          ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate, cmd_line, output_last,
                      called_variants,called_genes);   
  
  
  
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_vob, tmp_gi, abi,
          &is_mupirocin_susceptible, tmp_name, cmd_line->install_dir,
          ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate, cmd_line, output_last,
                      called_variants,called_genes);   
  
  
  
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_vob, tmp_gi, abi,
          &is_spectinomycin_susceptible, tmp_name, cmd_line->install_dir,
          ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate, cmd_line, output_last,
                      called_variants,called_genes);   
  
  
  
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_vob, tmp_gi, abi,
          &is_trimethoprim_susceptible, tmp_name, cmd_line->install_dir,
          ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate, cmd_line, output_last,
                      called_variants,called_genes);   
  
  
  
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_vob, tmp_gi, abi,
          &is_methicillin_susceptible, tmp_name, cmd_line->install_dir,
          ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate, cmd_line, output_last,
                      called_variants,called_genes);   
  
  
  
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_vob, tmp_gi, abi,
          &is_tetracycline_susceptible, tmp_name, cmd_line->install_dir,
          ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate, cmd_line, output_last,
                      called_variants,called_genes);   
  
  
  
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_vob, tmp_gi, abi,
          &is_streptothricin_susceptible, tmp_name, cmd_line->install_dir,
          ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate, cmd_line, output_last,
                      called_variants,called_genes);   
  
  
  
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_vob, tmp_gi, abi,
          &is_vancomycin_susceptible, tmp_name, cmd_line->install_dir,
          ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate, cmd_line, output_last,
                      called_variants,called_genes);   
  
  
  
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_vob, tmp_gi, abi,
          &is_fusidicacid_susceptible, tmp_name, cmd_line->install_dir,
          ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate, cmd_line, output_last,
                      called_variants,called_genes);   
  
  
  
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_vob, tmp_gi, abi,
          &is_rifampicin_susceptible, tmp_name, cmd_line->install_dir,
          ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate, cmd_line, output_last,
                      called_variants,called_genes);   
  
  
  
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_vob, tmp_gi, abi,
          &is_ciprofloxacin_susceptible, tmp_name, cmd_line->install_dir,
          ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate, cmd_line, output_last,
                      called_variants,called_genes);   
  
  
  output_last=true;
  print_clindamycin_susceptibility(db_graph, &file_reader_fasta, ru, tmp_vob, tmp_gi, abi,
           &is_clindamycin_susceptible, tmp_name, 
           any_erm_present,cmd_line->install_dir,
           ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate, cmd_line, output_last,
           called_variants,called_genes);

   
  if (cmd_line->format==JSON)
    {
      print_json_susceptibility_end();
    }
  print_called_variants(called_variants,cmd_line->format,false);// false here means that this is NOT the last element of the JSON output
  print_called_genes(called_genes,cmd_line->format);
  if (cmd_line->format==Stdout)
    {
      printf("** Virulence markers\n");
    }
  print_pvl_presence(db_graph, &file_reader_fasta, ru,  tmp_gi, 
		     &is_pvl_positive, 
		     cmd_line->install_dir, cmd_line->format); 

  if (cmd_line->format==Stdout)
    {
      timestamp();
    }
  else
    {
      print_json_end();
      printf("\n");
    }


  if ( (cmd_line ->output_supernodes==true) && (cmd_line->format==Stdout) )
    {
      printf("Print contigs\n");
      db_graph_print_supernodes_defined_by_func_of_colours(cmd_line->contig_file->buff, "", 
							   cmd_line->max_expected_sup_len,
							   db_graph, 
							   &element_get_colour_union_of_all_colours, 
							   &element_get_covg_union_of_all_covgs, 
							   &print_no_extra_supernode_info);
      printf("Completed printing contigs\n");
    }


  //cleanup
  strbuf_free(tmp_name);
  free_antibiotic_info(abi);
  free_var_on_background(tmp_vob);
  free_gene_info(tmp_gi);
  free_reading_utils(ru);
  
  cmd_line_free(cmd_line);
  hash_table_free(&db_graph);
  return 0;
}



void timestamp(){
 time_t ltime;
 ltime = time(NULL);
 printf("%s",asctime(localtime(&ltime)));
 fflush(stdout);
}