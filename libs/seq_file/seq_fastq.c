/*
 seq_fastq.c
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

#include "seq_fastq.h"

void _seq_read_fastq_sequence(SeqFile *sf)
{
  strbuf_reset(sf->bases_buff);

  int c;

  while((c = seq_getc(sf)) != -1 && c != '+')
  {
    if(c != '\r' && c != '\n')
    {
      strbuf_append_char(sf->bases_buff, (char)c);
      seq_readline(sf->bases_buff, sf);
      strbuf_chomp(sf->bases_buff);
    }

    sf->line_number++;
  }

  if(c == -1)
  {
    fprintf(stderr, "%s:%i: Missing '+' in FASTQ [file: %s; line: %lu]\n",
            __FILE__, __LINE__, sf->path, sf->line_number);
  }

  // Read to end of separator line
  if(c != '\r' && c != '\n')
  {
    seq_skip_line(sf);
    sf->line_number++;
  }
}

char seq_next_read_fastq(SeqFile *sf)
{
  if(sf->read_line_start)
  {
    // Read name
    seq_readline(sf->entry_name, sf);
    strbuf_chomp(sf->entry_name);
    sf->line_number++;

    // Read whole sequence
    _seq_read_fastq_sequence(sf);

    sf->read_line_start = 0;
    return 1;
  }
  else
  {
    int c;

    // Count bases not read in
    sf->total_bases_skipped += (strbuf_len(sf->bases_buff) - sf->entry_offset);

    // Skip over remaining quality values
    while(sf->entry_offset_qual < strbuf_len(sf->bases_buff))
    {
      if((c = seq_getc(sf)) == -1)
        return 0;

      if(c != '\r' && c != '\n')
        sf->entry_offset_qual++;
      else
        sf->line_number++;
    }

    // Skip newlines
    while((c = seq_getc(sf)) != -1 && (c == '\n' || c == '\r'))
      sf->line_number++;

    if(c == -1)
      return 0;

    if(c != '@')
    {
      fprintf(stderr, "%s:%i: FASTQ header does not begin with '@' (%c) "
                      "[file: %s; line: %lu]\n",
              __FILE__, __LINE__, c, sf->path, sf->line_number);
      return 0;
    }

    // Read name
    seq_readline(sf->entry_name, sf);
    strbuf_chomp(sf->entry_name);
    sf->line_number++;

    // Read whole sequence
    _seq_read_fastq_sequence(sf);

    return 1;
  }
}

char seq_read_base_fastq(SeqFile *sf, char *c)
{
  if(sf->entry_offset < strbuf_len(sf->bases_buff))
  {
    *c = strbuf_get_char(sf->bases_buff, sf->entry_offset);
    return 1;
  }
  else
  {
    return 0;
  }
}

char seq_read_qual_fastq(SeqFile *sf, char *c)
{
  if(sf->entry_offset_qual >= strbuf_len(sf->bases_buff))
    return 0;

  int next;

  while((next = seq_getc(sf)) != -1 && (next == '\n' || next == '\r'))
    sf->line_number++;

  if(next == -1)
  {
    fprintf(stderr, "%s:%i: FASTQ file ended without finishing quality scores "
                    "[file: %s; line: %lu]\n",
            __FILE__, __LINE__, sf->path, sf->line_number);
    return 0;
  }

  *c = (char)next;
  return 1;
}

char seq_read_all_bases_fastq(SeqFile *sf, StrBuf *sbuf)
{
  // Copy from buffer
  t_buf_pos len = sf->bases_buff->len - sf->entry_offset;
  strbuf_copy(sbuf, 0, sf->bases_buff, sf->entry_offset, len);

  return 1;
}

char seq_read_all_quals_fastq(SeqFile *sf, StrBuf *sbuf)
{
  if(sf->entry_offset_qual >= strbuf_len(sf->bases_buff))
    return 0;

  // Expect the same number of quality scores as bases
  t_buf_pos expected_len = strbuf_len(sf->bases_buff) -
                           sf->entry_offset_qual;

  int next = -1;
  t_buf_pos i;

  for(i = 0; i < expected_len && (next = seq_getc(sf)) != -1; i++)
  {
    if(next != '\r' && next != '\n')
    {
      strbuf_append_char(sbuf, (char)next);
    }
    else
      sf->line_number++;
  }

  if(next == -1)
  {
    fprintf(stderr, "%s:%i: FASTQ file ended without finishing quality "
                    "scores (FASTQ) [file: %s; line: %lu]\n",
            __FILE__, __LINE__, sf->path, sf->line_number);
  }

  return 1;
}

unsigned long seq_file_write_name_fastq(SeqFile *sf, const char *name)
{
  unsigned long num_bytes_printed = 0;

  if(sf->write_state == WS_BEGIN)
  {
    num_bytes_printed += seq_puts(sf, "@");
  }
  else if(sf->write_state == WS_NAME)
  {
    num_bytes_printed += seq_puts(sf, "\n\n+\n\n@");
  }
  else if(sf->write_state == WS_QUAL)
  {
    num_bytes_printed += seq_puts(sf, "\n@");
    sf->line_number++;
  }
  else if(sf->write_state == WS_SEQ)
  {
    fprintf(stderr, "%s:%i: writing in the wrong order (name) "
                    "[path: %s; line: %lu]\n",
            __FILE__, __LINE__, sf->path, sf->line_number);
    exit(EXIT_FAILURE);
  }

  num_bytes_printed += seq_puts(sf, name);

  return num_bytes_printed;
}

size_t seq_file_write_seq_fastq(SeqFile *sf, const char *seq, size_t str_len)
{
  if(sf->write_state == WS_BEGIN || sf->write_state == WS_QUAL)
  {
    fprintf(stderr, "%s:%i: writing in the wrong order (seq) [path: %s; line: %lu]\n",
            __FILE__, __LINE__, sf->path, sf->line_number);
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

size_t seq_file_write_qual_fastq(SeqFile *sf, const char *qual)
{
  if(sf->write_state == WS_BEGIN || sf->write_state == WS_NAME)
  {
    fprintf(stderr, "%s:%i: Writing in the wrong order (qual) [path: %s; line: %lu]\n",
            __FILE__, __LINE__, sf->path, sf->line_number);
    exit(EXIT_FAILURE);
  }

  size_t str_len = strlen(qual);

  if(str_len == 0)
    return 0;

  size_t num_bytes_printed = 0;

  if(sf->write_state == WS_SEQ)
  {
    num_bytes_printed += seq_puts(sf, "\n+\n");
    sf->line_number += 2;
  }

  num_bytes_printed += _write(sf, qual, str_len);
  sf->write_state = WS_QUAL;

  return num_bytes_printed;
}

size_t seq_file_close_write_fastq(SeqFile *sf)
{
  size_t num_bytes_printed = 0;

  if(sf->write_state == WS_NAME)
  {
    num_bytes_printed += seq_puts(sf, "\n\n+\n\n");
    sf->line_number += 4;
  }
  else if(sf->write_state == WS_SEQ)
  {
    num_bytes_printed += seq_puts(sf, "\n+\n\n");
    sf->line_number += 3;
  }
  else if(sf->write_state == WS_QUAL)
  {
    num_bytes_printed += seq_puts(sf, "\n");
    sf->line_number += 1;
  }

  return num_bytes_printed;
}
