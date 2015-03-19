/*
 seq_fasta.c
 project: seq_file
 url: https://github.com/noporpoise/seq_file
 author: Isaac Turner <turner.isaac@gmail.com>
 Copyright (C) 20-June-2012
 
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

#include "seq_fasta.h"

char seq_next_read_fasta(SeqFile *sf)
{
  if(sf->read_line_start)
  {
    // Read name
    seq_readline(sf->entry_name, sf);
    strbuf_chomp(sf->entry_name);

    sf->line_number++;
    sf->read_line_start = 0;

    return 1;
  }
  else
  {
    int c;

    // Look for line starting with >
    do
    {
      // Read until the end of the line
      while((c = seq_getc(sf)) != -1 && c != '\n' && c != '\r')
      {
        sf->total_bases_skipped++;
      }

      if(c == -1)
        return 0;
      else
        sf->line_number++; // Must have read a new line

      // Read through end of line chars
      while((c = seq_getc(sf)) != -1 && (c == '\n' || c == '\r'))
      {
        sf->line_number++;
      }

      if(c == -1)
      {
        return 0;
      }
      else if(c != '>')
      {
        sf->total_bases_skipped++;
      }
    }
    while(c != '>');

    // Read name
    seq_readline(sf->entry_name, sf);
    strbuf_chomp(sf->entry_name);
    sf->line_number++;

    return 1;
  }
}

char seq_read_base_fasta(SeqFile *sf, char *c)
{
  int next;
  
  while((next = seq_getc(sf)) != -1 && (next == '\n' || next == '\r'))
    sf->line_number++;

  if(next == -1)
  {
    return 0;
  }
  else if(next == '>')
  {
    sf->read_line_start = 1;
    return 0;
  }
  else
  {
    *c = (char)next;
    return 1;
  }
}

char seq_read_all_bases_fasta(SeqFile *sf, StrBuf *sbuf)
{
  int c;

  while((c = seq_getc(sf)) != -1 && c != '>')
  {
    if(c != '\r' && c != '\n')
    {
      strbuf_append_char(sbuf, (char)c);
      seq_readline(sbuf, sf);
      strbuf_chomp(sbuf);
    }

    sf->line_number++;
  }

  if(c == '>')
    sf->read_line_start = 1;

  return 1;
}

unsigned long seq_file_write_name_fasta(SeqFile *sf, const char *name)
{
  size_t num_bytes_printed = 0;

  if(sf->write_state == WS_BEGIN)
  {
    num_bytes_printed += seq_puts(sf, ">");
  }
  else
  {
    num_bytes_printed += seq_puts(sf, "\n>");
    sf->line_number++;
  }

  num_bytes_printed += seq_puts(sf, name);

  return num_bytes_printed;
}

size_t seq_file_write_seq_fasta(SeqFile *sf, const char *seq, size_t str_len)
{
  if(sf->write_state == WS_BEGIN)
  {
    fprintf(stderr, "%s:%i: writing in the wrong order (seq) [path: %s]\n",
            __FILE__, __LINE__, sf->path);
    exit(EXIT_FAILURE);
  }

  size_t num_bytes_printed = 0;

  if(sf->write_state == WS_NAME)
  {
    num_bytes_printed += seq_puts(sf, "\n");
    sf->line_number++;
  }

  num_bytes_printed += _write(sf, seq, str_len);

  return num_bytes_printed;
}

size_t seq_file_close_write_fasta(SeqFile *sf)
{
  size_t num_bytes_printed = 0;

  if(sf->write_state == WS_NAME || sf->write_state == WS_SEQ)
  {
    num_bytes_printed += seq_puts(sf, "\n");
    sf->line_number++;
  }

  return num_bytes_printed;
}
