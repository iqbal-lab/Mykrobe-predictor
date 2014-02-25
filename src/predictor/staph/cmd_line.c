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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <err.h>

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
"   [--list FILENAME] \t\t\t\t\t=\t List of fastq.\n" ;

int default_opts(CmdLine * c)
{
  strbuf_reset(c->list_of_fastq);
  c->genome_size = 2800000;
  c->kmer_size = 31;
  c->mem_width = 100;
  c->mem_height= 19;
  c->max_expected_sup_len=50000;
  return 1;
}


CmdLine* cmd_line_alloc()
{
  CmdLine* cmd = (CmdLine*) malloc(sizeof(CmdLine));
  if (cmd==NULL)
    {
      die("Out of memory before we even start Cortex, cannot alloc space to hold the commandline variables! Abort\n");
    }
  cmd->list_of_fastq = strbuf_new();

  int max_expected_read_len = 500;//illumina
  cmd->readlen_distrib_size = max_expected_read_len + 1;
  cmd->readlen_distrib
    = (uint64_t*) calloc(cmd->readlen_distrib_size, sizeof(uint64_t));
    
  if(cmd->readlen_distrib == NULL)
    {
      die("Unable to malloc array to hold readlen distirbution!Exit.\n");
    }


  int max_expected_mean_kmer_covg = 1000;
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
  strbuf_free(cmd->list_of_fastq);
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
    {"list",required_argument, NULL, 'f'},
    {0,0,0,0}	
  };
  

  //do not change this! Only really matters for testing, but getopt_long uses 
  //variables which are not local to this function, so when we run tests, and call this function
  // repeatedly, those variable values are carried across. The following line resets this.
  optind=1;
  
 
  opt = getopt_long(argc, argv, "hf:", long_options, &longopt_index);

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

    case 'f':
      {
	if (access(optarg,F_OK)==0) 
	  {
	    strbuf_append_str(cmdline_ptr->list_of_fastq, optarg);
	  }
	else
	  {
	    errx(1,"Cannot open file %s",optarg);
	    return -1;
	  }
	break;
      }
    default:
      {
	errx(1, "Unknown option %c\n", opt);
	return -1;
      }      

    }
    opt = getopt_long(argc, argv, "hf:", long_options, &longopt_index);
    
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
  if (strcmp(cmd_line->list_of_fastq->buff, "")==0)
    {
      die("You must speoify --list\n");
    }

  return 0;
}


