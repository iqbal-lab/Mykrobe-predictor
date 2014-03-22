/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * This file is part of myKrobe
 *
 * myKrobe is free software: you can redistribute it and/or modify
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
 * along with myKrobe.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */
/*
  cmd_line.c
*/

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <err.h>
#include <errno.h>

// myKrobe headers
#include "cmd_line.h"
#include "file_reader.h"
#include "file_format.h"

#define MAX_LINE 500

// Supress a warning from the compiler about an unused variable
#define UNUSED(x) ((void)(x))


const char* usage=
"myKrobe predictorper\n"\
"   [--help] \t\t\t\t\t\t\t=\t This help screen.\n\n" \
"   [--list FILENAME] \t\t\t\t\t=\t List of fastq or bam. Cannot use --list and --file\n" \
"   [--file FILENAME] \t\t\t\t\t=\t Single fastq or bam. Cannot use --file and --list\n" \
"   [--sample_id STRING] \t\t\t\t\t=\t Identifier for sample under test\n" \
"   [--method STRING] \t\t\t\t\t=\t Default is WGAssemblyThenGenotyping. Or can have InSilicoOligos\n" \
"   [--install_dir PATH] \t\t\t\t\t=\t myKrobe.predictor needs to use config files that come in the install, so you need to specify the full path to your install\n\n" ;

int default_opts(CmdLine * c)
{
  strbuf_reset(c->seq_path);
  strbuf_reset(c->id);
  strbuf_append_str(c->id, "UnknownSample");
  strbuf_reset(c->install_dir);
  c->genome_size = 2800000;
  c->num_bases_around_mut_in_fasta=30;//our antibiotic fasta have 30 bases before/after the mutation
  c->kmer_size = 15;
  c->mem_width = 100;
  c->mem_height= 19;
  c->max_expected_sup_len=50000;
  c->method=WGAssemblyThenGenotyping;//other options InSilicoOligos and WGAssemblyAndTranslation
  c->input_file=false;
  c->input_list=false;
  return 1;
}


CmdLine* cmd_line_alloc()
{
  CmdLine* cmd = (CmdLine*) malloc(sizeof(CmdLine));
  if (cmd==NULL)
    {
      die("Out of memory before we even start Cortex, cannot alloc space to hold the commandline variables! Abort\n");
    }
  cmd->seq_path = strbuf_new();
  cmd->install_dir = strbuf_new();
  cmd->id = strbuf_new();
  int max_expected_read_len = 500;//illumina
  cmd->readlen_distrib_size = max_expected_read_len + 1;
  cmd->readlen_distrib
    = (uint64_t*) calloc(cmd->readlen_distrib_size, sizeof(uint64_t));
    
  if(cmd->readlen_distrib == NULL)
    {
      die("Unable to malloc array to hold readlen distirbution!Exit.\n");
    }


  int max_expected_mean_kmer_covg = 500;
  cmd->len_kmer_covg_array = max_expected_mean_kmer_covg+1;
  cmd->kmer_covg_array
    = (uint64_t*) calloc(cmd->len_kmer_covg_array, sizeof(uint64_t));
    
  if(cmd->kmer_covg_array== NULL)
    {
      die("Unable to malloc array to hold kmer_covgs! Exit.\n");
    }
  
 return cmd;
}

void cmd_line_free(CmdLine* cmd)
{
  strbuf_free(cmd->seq_path);
  strbuf_free(cmd->install_dir);
  strbuf_free(cmd->id);
  free(cmd->readlen_distrib);
  free(cmd);
}


//inner loop called by the parse_cmdline function, which returns various error codes.
int parse_cmdline_inner_loop(int argc, char* argv[], int unit_size, CmdLine* cmdline_ptr, char* error_string)
{
  int opt;
  int longopt_index;
  
  static struct option long_options[] = {
    {"help",no_argument,NULL,'h'},
    {"list",required_argument, NULL, 'l'},
    {"file",required_argument, NULL, 'f'},
    {"method", required_argument, NULL, 'm'},
    {"sample_id", required_argument, NULL, 's'},
    {"oligo_bin", required_argument, NULL, 'b'},
    {"install_dir", required_argument, NULL, 'i'},
    {0,0,0,0}	
  };
  

  //do not change this! Only really matters for testing, but getopt_long uses 
  //variables which are not local to this function, so when we run tests, and call this function
  // repeatedly, those variable values are carried across. The following line resets this.
  optind=1;
  
 
  opt = getopt_long(argc, argv, "hf:l:m:s:", long_options, &longopt_index);

  while ((opt) > 0) {
	       
    //Parse the default options
    switch(opt) {

    case 'h':
      {
	printf("***********************\n");
	printf("myKrobe.predictor for Staphylococcus\n");
	printf("***********************\n");

	printf("%s",usage);
	exit(0);
	break;
      }

    case 'i':
      {
	StrBuf* tmp = strbuf_create(optarg);
	strbuf_add_slash_on_end(tmp);
	if (0 != access(tmp->buff, F_OK)) {
	  if (ENOENT == errno) {
	    die("You have specified with --install_dir, a directory which does not exist (%s)\n", optarg);
	  }
	  if (ENOTDIR == errno) {
	    // not a directory
	    die("You have specified with --install_dir, a directory which is not a directory (%s)\n", optarg);
	  }
	}	
	strbuf_append_str(cmdline_ptr->install_dir, tmp->buff);
	strbuf_append_str(tmp, "data/staph/antibiotics/penicillin.fa");
	if (access(tmp->buff,F_OK)!=0)
	  {
	    die("You have specified with --install_dir, a directory which does not seem to be the install directory of myKrobe.predictor. Cannot find %s\n", tmp->buff);
	  }
	
	break;
      }
    case 'l':
      {
	if (access(optarg,F_OK)==0) 
	  {
	    char* full_path = realpath(optarg,NULL);
	    strbuf_append_str(cmdline_ptr->seq_path,full_path );
	    free(full_path);
	    cmdline_ptr->input_list=true;
	  }
	else
	  {
	    errx(1,"Cannot open file %s",optarg);
	    return -1;
	  }
	break;
      }
    case 'f':
      {
	if (access(optarg,F_OK)==0) 
	  {
	    char* full_path = realpath(optarg,NULL);
	    strbuf_append_str(cmdline_ptr->seq_path, full_path);
	    free(full_path);
	    cmdline_ptr->input_file=true;
	  }
	else
	  {
	    errx(1,"Cannot open file %s",optarg);
	    return -1;
	  }
	break;
      }
    case 'm':
      {
	if (strcmp(optarg, "WGAssemblyThenGenotyping")==0)
	  {
	    cmdline_ptr->method=WGAssemblyThenGenotyping;
	  }
	else if (strcmp(optarg, "InSilicoOligos")==0)
	  {
	    cmdline_ptr->method=InSilicoOligos;
	    /* cmdline_ptr->mem_height=19;
	       cmdline_ptr->mem_width=100; */
	  }
	else if (strcmp(optarg, "WGAssemblyAndTranslation")==0)
	  {
	    cmdline_ptr->method=WGAssemblyAndTranslation;
	  }
	else
	  {
	    errx(1, "--method requires an argument, which must be one of WGAssemblyThenGenotyping, InSilicoOligos, WGAssemblyAndTranslation\n");
	  }
	break;
      }
    case 's':
      {
	strbuf_reset(cmdline_ptr->id);
	strbuf_append_str(cmdline_ptr->id, optarg);
	break;
      }
    default:
      {
	errx(1, "Unknown option %c\n", opt);
	return -1;
      }      

    }
    opt = getopt_long(argc, argv, "hf:l:m:s:", long_options, &longopt_index);
    
  }   
  
  return 0;
}
  

void parse_cmdline(CmdLine* cmd_line, 
		   int argc, 
		   char* argv[], 
		   int unit_size) 
{	
  int i;
  printf("Command: ");
  for(i=0;i<argc;i++){
    printf("%s ",argv[i]);
  }
  printf("\n");

  default_opts(cmd_line);

  char error_string[LEN_ERROR_STRING]="";
  int err = parse_cmdline_inner_loop(argc, argv, unit_size, cmd_line, error_string);
  
  if (err==-1) 
    {
      die("Error in cmd-line input: %s\n", error_string);
    }

  check_cmdline(cmd_line, error_string);
}

int check_cmdline(CmdLine* cmd_line, char* error_string)
{
  if ( (cmd_line->input_file==true) && (cmd_line->input_list==true) )
    {
      die("You cannot specify both --list and --file\n");
    }
  else if ((cmd_line->input_file==false) && (cmd_line->input_list==false) )
    {
      die("You must specify one of  --list (list of FASTQ or BAM files to load)  or  --file (single FASTQ or BAM). In both cases, these are assumed to come from one sample\n");
    }
  if ( 
      (cmd_line->method!=WGAssemblyThenGenotyping)
      &&
      (cmd_line->method!=InSilicoOligos)
       )
    {
      die("--method only takes argument WGAssemblyThenGenotyping or InSilicoOligos\n");
    }
  if (strcmp(cmd_line->install_dir->buff, "")==0)
    {
      die("You must specify --install_dir\n");
    }
  return 0;
}


