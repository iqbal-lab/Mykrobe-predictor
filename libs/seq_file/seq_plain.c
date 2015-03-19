/*
 seq_plain.c
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

#include "seq_plain.h"

char seq_next_read_plain(SeqFile *sf)
{
  int c;

  if(!sf->read_line_start)
  {
    // Skip the rest of line
    while((c = seq_getc(sf)) != -1 && c != '\n' && c != '\r')
    {
      sf->total_bases_skipped++;
    }

    sf->read_line_start = 1;
  }

  // check we can get another base
  if((c = seq_getc(sf)) == -1)
  {
    return 0;
  }

  seq_ungetc(c, sf);

  return 1;
}

char seq_read_base_plain(SeqFile *sf, char *cc)
{
  int c = seq_getc(sf);

  if(c == -1)
  {
    return 0;
  }
  else if(c == '\r' || c == '\n')
  {
    sf->line_number++;
    sf->read_line_start = 1;
    return 0;
  }
  else
  {
    *cc = (char)c;
    return 1;
  }
}

char seq_read_all_bases_plain(SeqFile *sf, StrBuf *sbuf)
{
  seq_readline(sbuf, sf);

  // chomp returns the number of charactes removed
  if(strbuf_chomp(sbuf) > 0)
  {
    sf->read_line_start = 1;
    sf->line_number++;
  }

  return 1;
}

size_t seq_file_write_seq_plain(SeqFile *sf, const char *seq)
{
  size_t num_bytes_printed = 0;

  if(sf->write_state != WS_BEGIN)
  {
    num_bytes_printed += seq_puts(sf, "\n");
    sf->line_number++;
  }

  num_bytes_printed += seq_puts(sf, seq);

  return num_bytes_printed;
}

size_t seq_file_close_write_plain(SeqFile *sf)
{
  size_t num_bytes_printed = 0;

  if(sf->write_state == WS_SEQ)
  {
    num_bytes_printed += seq_puts(sf, "\n");
    sf->line_number++;
  }

  return num_bytes_printed;
}

