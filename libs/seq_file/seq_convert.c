/*
 seq_convert.c
 project: seq_file
 url: https://github.com/noporpoise/seq_file
 author: Isaac Turner <turner.isaac@gmail.com>
 Copyright (C) 21 June 2012

 To build type 'make'

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "string_buffer.h"
#include "seq_file.h"

void print_usage()
{
  printf("usage: seq_convert <in> <out> [line-wrap]\n"
"  Output file must end with: fa|fq|fa.gz|fq.gz|sam|bam\n"
"  If line-wrap value is given, lines are wrapped.\n");
  exit(EXIT_FAILURE);
}

int main(int argc, char** argv)
{
  if(argc < 3 || argc > 4)
  {
    print_usage();
  }

  char* in_path = argv[1];
  char* out_path = argv[2];

  unsigned long line_wrap = 0;
  
  if(argc == 4)
  {
    char *line_wrap_str = argv[3];
    char *endptr;
    line_wrap = strtoul(line_wrap_str, &endptr, 10);

    if((unsigned)(endptr-line_wrap_str) != strlen(line_wrap_str))
    {
      print_usage();
    }
  }

  SeqFileType out_file_type = SEQ_UNKNOWN;
  char out_zipped = 0;

  seq_guess_filetype_from_path(out_path, &out_file_type, &out_zipped);

  if(out_file_type == SEQ_UNKNOWN)
  {
    fprintf(stderr, "%s:%i: Sorry, I cannot identify the output file's format "
                    "from its path [file: %s]\n", __FILE__, __LINE__, out_path);
    exit(EXIT_FAILURE);
  }

  SeqFile* in_file = seq_file_open(in_path);
  SeqFile* out_file = seq_file_open_write(out_path, out_file_type,
                                          out_zipped, line_wrap);

  if(in_file == NULL)
  {
    fprintf(stderr, "%s:%i: Couldn't open input file: %s\n",
            __FILE__, __LINE__, in_path);
    exit(EXIT_FAILURE);
  }

  if(out_file == NULL)
  {
    fprintf(stderr, "%s:%i: Couldn't open output file: %s\n",
            __FILE__, __LINE__, out_path);
    exit(EXIT_FAILURE);
  }

  printf(" In : %s [%s]\n", in_path, seq_file_get_type_str(in_file));
  printf(" Out: %s [%s]\n", out_path, seq_file_get_type_str(out_file));

  // Start converting
  size_t bytes_written = 0;
  char c[2] = ".";

  if(out_file_type == SEQ_PLAIN)
  {
    // Example reading in an entire entry at a time using seq_read_all_bases()
    StrBuf *bases = strbuf_new();

    while(seq_next_read(in_file))
    {
      while(seq_read_all_bases(in_file, bases))
      {
        bytes_written += seq_file_write_seq(out_file, bases->buff);
      }
    }
  }
  else
  {
    // Example reading in a char at a time using seq_read_base()
    while(seq_next_read(in_file))
    {
      const char* read_name = seq_get_read_name(in_file);
      bytes_written += seq_file_write_name(out_file, read_name);

      unsigned long seq_length = 0;

      while(seq_read_base(in_file, c))
      {
        seq_length++;

        if(!(bytes_written += seq_file_write_seq(out_file, c)))
        {
          fprintf(stderr, "%s:%i: Couldn't write base to file "
                          "[file: %s; line: %lu]\n",
                  __FILE__, __LINE__,
                  seq_get_path(out_file), seq_curr_line_number(in_file));

          exit(EXIT_FAILURE);
        }
      }
    
      if(seq_has_quality_scores(out_file))
      {
        size_t bytes_written_before_qual = bytes_written;

        while(seq_read_qual(in_file, c))
        {
          if(!(bytes_written += seq_file_write_qual(out_file, c)))
          {
            fprintf(stderr, "%s:%i: Couldn't write quality score to file "
                            "[file: %s; line: %lu]\n",
                    __FILE__, __LINE__, seq_get_path(out_file),
                    seq_curr_line_number(in_file));

            exit(EXIT_FAILURE);
          }
        }

        if(bytes_written == bytes_written_before_qual)
        {
          // No quality scores were read - fill in
          unsigned long i;
          *c = '?';
          for(i = 0; i < seq_length; i++)
            bytes_written += seq_file_write_qual(out_file, c);
        }
      }
    }
  }

  unsigned long seq_total_bases_read = seq_total_bases_passed(in_file);
  unsigned long total_entries = seq_get_read_index(in_file);

  seq_file_close(in_file);
  bytes_written += seq_file_close(out_file);

  printf("%lu entries read\n", total_entries);
  printf("%lu bases read\n", seq_total_bases_read);
  printf("%lu bytes written\n", bytes_written);
  printf("Done. \n");

  return EXIT_SUCCESS;
}
