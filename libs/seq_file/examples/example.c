/*
 examples/example.c
 project: seq_file
 url: https://github.com/noporpoise/seq_file
 author: Isaac Turner <turner.isaac@gmail.com>
 Copyright (C) 13-Dec-2013
 
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

#include "seq_file.h"

int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    fprintf(stderr, "Usage: example <inputfile>\n");
    fprintf(stderr, "  Example code reads in a sequence file and prints as FASTQ\n");
    exit(EXIT_FAILURE);
  }

  const char *path = argv[1];

  SeqFile *sf = seq_file_open(path);

  if(sf == NULL)
  {
    // Error opening file
    fprintf(stderr, "Error: cannot read seq file '%s'\n", path);
    exit(EXIT_FAILURE);
  }

  char c, q;

  while(seq_next_read(sf))
  {
    printf("@%s\n", seq_get_read_name(sf));

    while(seq_read_base(sf, &c))
      putc(c, stdout);

    printf("\n+\n");

    if(seq_has_quality_scores(sf))
    {
      // Input file has quality scores, read and print them
      while(seq_read_qual(sf, &q))
        fputc(q, stdout);
    }
    else
    {
      // Print enough '.' characters in place of quality scores
      int i, len = seq_get_length(sf);
      for(i = 0; i < len; i++)
      {
        putc('.', stdout);
      }
    }

    putc('\n', stdout);
  }

  seq_file_close(sf);
}
