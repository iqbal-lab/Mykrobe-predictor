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
  main.c 
*/

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <string_buffer.h>

// cortex_var headers
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

void timestamp();

int main(int argc, char **argv)
{

  // VERSION_STR is passed from the makefile -- usually last commit hash
  printf("myKrobe.predictor for Staphylococcus, version %d.%d.%d.%d"VERSION_STR"\n",
         VERSION, SUBVERSION, SUBSUBVERSION, SUBSUBSUBVERSION);

  CmdLine* cmd_line = cmd_line_alloc();
  if (cmd_line==NULL)
    {
      return -1;
    }
  
  parse_cmdline(cmd_line, argc,argv,sizeof(Element));
  printf("Kmer size %d\n", cmd_line->kmer_size);
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




  int lim = cmd_line->max_expected_sup_len;
  CovgArray* working_ca_for_median=alloc_and_init_covg_array(lim);//will die if fails to alloc
  if (working_ca_for_median==NULL)
    {
      return -1;
    }

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
  ResVarInfo* tmp_rvi = alloc_and_init_res_var_info();
  GeneInfo* tmp_gi = alloc_and_init_gene_info();
  AntibioticInfo* abi = alloc_antibiotic_info();
  if ( (ru==NULL) || (tmp_rvi==NULL) || (abi==NULL) || (tmp_gi==NULL) )
    {
      return -1;
    }
  

  printf("** Start time\n");
  timestamp();
  printf("** Sample:\n%s\n", cmd_line->id->buff);

 
  int into_colour=0;
  boolean only_load_pre_existing_kmers=false;
  //  GraphInfo* gi = graph_info_alloc_and_init();
  uint64_t bp_loaded=0;

  if (cmd_line->method==WGAssemblyThenGenotyping)
    {
      only_load_pre_existing_kmers=false;
    }
  else if (cmd_line->method==InSilicoOligos)
    {
      StrBuf* skeleton_flist = strbuf_new();
      strbuf_append_str(skeleton_flist, 
			cmd_line->install_dir->buff);
      strbuf_append_str(skeleton_flist, 
			"data/skeleton_binary/list_speciesbranches_genes_and_muts");
      build_unclean_graph(db_graph, 
			  skeleton_flist,
			  true,
			  cmd_line->kmer_size,
			  NULL, 0,
			  NULL, 0,
			  false,
			  into_colour,
			  &subsample_null);
      set_all_coverages_to_zero(db_graph, 0);
      strbuf_free(skeleton_flist);
      only_load_pre_existing_kmers=true;
    }
  else
    {
      die("For now --method only allowed to take InSilicoOligos or WGAssemblyThenGenotyping\n");
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
				  into_colour, subsample_function);

  if (bp_loaded==0)
    {
      printf("No data\n");
      return 1;
    }
  
  unsigned long mean_read_length = calculate_mean_uint64_t(cmd_line->readlen_distrib,
							   cmd_line->readlen_distrib_size);

  double err_rate = estimate_err_rate(cmd_line->seq_path, cmd_line->input_list);

  //given the error rate and other params, we can estimate expected depth of covg, and read-arrival rate
  // lambda_g = Depth/read_len _g means lambda on the true genome
  double lambda_g_err_free = ((double) bp_loaded/(double)(cmd_line->genome_size)) / (double) mean_read_length ; 
  
  int expected_depth 
    = (int) ( pow(1-err_rate, cmd_line->kmer_size)  
	      * (mean_read_length-cmd_line->kmer_size+1)
	      * lambda_g_err_free );
  
  printf("Expected covg\t%d\n", expected_depth);
  clean_graph(db_graph, cmd_line->kmer_covg_array, cmd_line->len_kmer_covg_array,
  	      expected_depth, cmd_line->max_expected_sup_len);
  
  
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
  Staph_species sp = get_species(db_graph, 10000, cmd_line->install_dir,
				 1,1);
  map_species_enum_to_str(sp,tmp_name);
  printf("** Species\n%s\n", tmp_name->buff);
  if (sp != Aureus)
    {
      printf("** No AMR predictions for coag-negative\n** End time\n");
      timestamp();
      return 1;
    }
  else
    {
      printf("** Antimicrobial susceptibility predictions\n");
    }
  
  //assumption is num_bases_around_mut_in_fasta is at least 30, to support all k<=31.
  //if k=31, we want to ignore 1 kmer at start and end
  //if k=29, we want to ignore 3 kmers at start and end.. etc


  int ignore = cmd_line->num_bases_around_mut_in_fasta - cmd_line->kmer_size +2;  

  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_rvi, tmp_gi, abi,
				  &is_gentamycin_susceptible, tmp_name, cmd_line->install_dir,
				  ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate);  
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_rvi, tmp_gi, abi,
				  &is_penicillin_susceptible, tmp_name, cmd_line->install_dir,
				  ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate); 
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_rvi, tmp_gi, abi,
				  &is_trimethoprim_susceptible, tmp_name, cmd_line->install_dir,
				  ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate); 
    boolean any_erm_present=false;
  print_erythromycin_susceptibility(db_graph, &file_reader_fasta, ru, tmp_rvi, tmp_gi, abi,
				    &is_erythromycin_susceptible, tmp_name, cmd_line->install_dir,
				    ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate,
				    &any_erm_present);

  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_rvi, tmp_gi, abi,
				  &is_methicillin_susceptible, tmp_name, cmd_line->install_dir,
				  ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate);
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_rvi, tmp_gi, abi,
				  &is_fusidic_acid_susceptible, tmp_name, cmd_line->install_dir,
				  ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate);
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_rvi, tmp_gi, abi,
				  &is_ciprofloxacin_susceptible, tmp_name, cmd_line->install_dir,
				  ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate);
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_rvi, tmp_gi, abi,
				  &is_rifampicin_susceptible, tmp_name, cmd_line->install_dir,
				  ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate);
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_rvi, tmp_gi, abi,
				  &is_tetracycline_susceptible, tmp_name, cmd_line->install_dir,
				  ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate);
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_rvi, tmp_gi, abi,
				  &is_vancomycin_susceptible, tmp_name, cmd_line->install_dir,
				  ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate);
  print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_rvi, tmp_gi, abi,
				  &is_mupirocin_susceptible, tmp_name, cmd_line->install_dir,
				  ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate);
  print_clindamycin_susceptibility(db_graph, &file_reader_fasta, ru, tmp_rvi, tmp_gi, abi,
				   &is_clindamycin_susceptible, tmp_name, 
				   any_erm_present,cmd_line->install_dir,
				   ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate);
  printf("** Virulence markers\n");
  print_pvl_presence(db_graph, &file_reader_fasta, ru,  tmp_gi, 
		     &is_pvl_positive, cmd_line->install_dir); 

  timestamp();

  if (cmd_line ->output_supernodes==true)
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
  free_res_var_info(tmp_rvi);
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
