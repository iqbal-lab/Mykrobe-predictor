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
  cmd_line.h
*/


#ifndef CMD_LINE_H_
#define CMD_LINE_H_

#include <stdio.h>

#include "global.h"

#define MAX_FILENAME_LEN 1000
#define LEN_ERROR_STRING 400


typedef struct
{
  StrBuf* list_of_fastq;
  long long genome_size;
  uint16_t kmer_size;
  uint64_t* readlen_distrib;
  uint64_t readlen_distrib_size;
  uint64_t* kmer_covg_array;
  int len_kmer_covg_array;
  int bucket_size;
  int number_of_buckets_bits;
  int max_expected_sup_len;
  int mem_height;
  int mem_width;
} CmdLine;


CmdLine* cmd_line_alloc();
void cmd_line_free(CmdLine* cmd);

int parse_cmdline_inner_loop(int argc, char* argv[], int unit_size, CmdLine* cmdline_ptr, char* error_string);
int check_cmdline(CmdLine* cmd_ptr, char* error_string);
void parse_cmdline(CmdLine* cmd_line, int argc, char* argv[],int unit_size); 
int default_opts(CmdLine *);

#endif /* CMD_LINE_H_ */
