/*
 seq_plain.h
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

#ifndef SEQ_PLAIN_HEADER_SEEN
#define SEQ_PLAIN_HEADER_SEEN

#include "seq_common.h"

char seq_next_read_plain(SeqFile *sf);
char seq_read_base_plain(SeqFile *sf, char *c);
char seq_read_all_bases_plain(SeqFile *sf, StrBuf *sbuf);

size_t seq_file_write_seq_plain(SeqFile *sf, const char *seq);

size_t seq_file_close_write_plain(SeqFile *sf);

#endif
