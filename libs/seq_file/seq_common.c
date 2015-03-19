/*
 seq_common.c
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

#include "seq_common.h"

size_t _write_wrapped(SeqFile *sf, const char *str, size_t str_len)
{
  size_t num_bytes_printed = 0;

  if(sf->curr_line_length == sf->line_wrap)
  {
    sf->curr_line_length = 0;
    num_bytes_printed += seq_puts(sf, "\n");
    sf->line_number++;
  }
  
  if(sf->curr_line_length + str_len <= sf->line_wrap)
  {
    // Doesn't go over a single line
    sf->curr_line_length += str_len;
    num_bytes_printed += seq_puts(sf, str);
    return num_bytes_printed;
  }

  size_t bytes_to_print = sf->line_wrap - sf->curr_line_length;

  num_bytes_printed += seq_write(sf, str, bytes_to_print);
  num_bytes_printed += seq_puts(sf, "\n");
  sf->line_number++;

  size_t offset;

  for(offset = bytes_to_print; offset < str_len; offset += sf->line_wrap)
  {
    bytes_to_print = MIN(str_len - offset, sf->line_wrap);
    num_bytes_printed += (size_t)seq_write(sf, str + offset, bytes_to_print);

    if(bytes_to_print < sf->line_wrap)
    {
      num_bytes_printed += seq_puts(sf, "\n");
      sf->line_number++;
    }
  }

  sf->curr_line_length = bytes_to_print;

  return num_bytes_printed;
}

