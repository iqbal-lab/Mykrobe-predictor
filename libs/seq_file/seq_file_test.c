/*
 seq_file_test.c
 project: seq_file
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

#include "seq_file.h"

/*
void test_seq_reader(char* file_path)
{
  Seq* seq = seq_create();
  SequenceFile* file = seq_file_open(file_path);

  printf(" Filetype: %s\n", seq_file_get_type_str(file));

  int i;

  for(i = 0; seq_file_read(file, seq); i++)
  {
    printf("%i>%s\n%s\n%s\n", i, seq->name->buff,
           seq->seq->buff, seq->qual->buff);
  }

  seq_file_close(file);
  seq_free(seq);
}
*/

void test_read_all_bases(char *file_path)
{
  SeqFile *sf = seq_file_open(file_path);

  if(sf == NULL)
  {
    fprintf(stderr, "%s:%i: Couldn't open file: %s\n",
            __FILE__, __LINE__, file_path);
    return;
  }

  while(seq_next_read(sf))
  {
    printf("read: %s\n", seq_get_read_name(sf));

    printf("  ");
    char c;
    while(seq_read_base(sf, &c))
      printf("%c", c);

    printf("\n");
  }
}

void test_read_names(char *file_path)
{
  SeqFile *sf = seq_file_open(file_path);

  if(sf == NULL)
  {
    fprintf(stderr, "%s:%i: Couldn't open file %s\n",
            __FILE__, __LINE__, file_path);
    return;
  }

  while(seq_next_read(sf))
  {
    printf("read: %s\n", seq_get_read_name(sf));
  }

  printf("%lu bases read\n", seq_total_bases_passed(sf));
  printf("%lu bases skipped\n", seq_total_bases_skipped(sf));
}

void test_get_type(char* file_path)
{
  SeqFile* file = seq_file_open(file_path);

  if(file == NULL)
  {
    fprintf(stderr, "%s:%i: Cannot open file: %s\n",
            __FILE__, __LINE__, file_path);
    return;
  }

  printf(" Filetype: %s\n", seq_file_get_type_str(file));
  seq_file_close(file);
}

int main(int argc, char** argv)
{
  if(argc != 2)
  {
    printf("usage: seq_file_test <in.fa|fq|sam|bam>\n");
    return -1;
  }

  // Print file path
  char* file_path = argv[1];
  printf(" Filepath: %s\n", file_path);

  // Test get type
  //test_get_type(file_path);

  test_read_names(file_path);
  
  //test_read_all_bases(file_path);

  return EXIT_SUCCESS;
}
