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


int main(int argc, char **argv)
{
  // VERSION_STR is passed from the makefile -- usually last commit hash
  printf("Starting myKrobe.predictor for Staphylococcus, version %d.%d.%d.%d"VERSION_STR"\n",
         VERSION, SUBVERSION, SUBSUBVERSION, SUBSUBSUBVERSION);

  CmdLine* cmd_line = cmd_line_alloc();
  if (cmd_line==NULL)
    {
      return -1;
    }
  
  parse_cmdline(cmd_line, argc,argv,sizeof(Element));

  int hash_key_bits, bucket_size;
  dBGraph * db_graph = NULL;
  short kmer_size;



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

  uint64_t bp_loaded = build_unclean_graph(db_graph, 
					   cmd_line->list_of_fastq, 
					   cmd_line->kmer_size,
					   cmd_line->readlen_distrib,
					   cmd_line->readlen_distrib_size,
					   cmd_line->kmer_covg_array, 
					   cmd_line->len_kmer_covg_array);

  unsigned long mean_read_length = calculate_mean_uint64_t(cmd_line->readlen_distrib,
							   cmd_line->readlen_distrib_size);


  int expected_depth = (mean_read_length-cmd_line->kmer_size+1)*(bp_loaded/cmd_line->genome_size) / mean_read_length;


  clean_graph(db_graph, cmd_line->kmer_covg_array, cmd_line->len_kmer_covg_array,
	      expected_depth, cmd_line->max_expected_sup_len);

  cmd_line_free(cmd_line);
  hash_table_free(&db_graph);
  return 0;
}



void timestamp(){
 time_t ltime;
 ltime = time(NULL);
 printf("\n-----\n%s",asctime(localtime(&ltime)));
 fflush(stdout);
}
