/*
 seq_fastq.h
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

#ifndef SEQ_FASTQ_HEADER_SEEN
#define SEQ_FASTQ_HEADER_SEEN

#include "seq_common.h"

char seq_next_read_fastq(SeqFile *sf);

char seq_read_base_fastq(SeqFile *sf, char *c);
char seq_read_qual_fastq(SeqFile *sf, char *c);

char seq_read_all_bases_fastq(SeqFile *sf, StrBuf *sbuf);
char seq_read_all_quals_fastq(SeqFile *sf, StrBuf *sbuf);

size_t seq_file_write_name_fastq(SeqFile *sf, const char *name);
size_t seq_file_write_seq_fastq(SeqFile *sf, const char *seq, size_t str_len);
size_t seq_file_write_qual_fastq(SeqFile *sf, const char *seq);

size_t seq_file_close_write_fastq(SeqFile *sf);

#endif
