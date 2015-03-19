/*
 seq_common.h
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

#ifndef SEQ_COMMON_HEADER_SEEN
#define SEQ_COMMON_HEADER_SEEN

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include "sam.h"

#include "seq_file.h"

// Printed nothing, then name, then some sequence, then some qualities
// ... then ready to print a name again
typedef enum WriteState
  {WS_READ_ONLY, WS_BEGIN, WS_NAME, WS_SEQ, WS_QUAL} WriteState;

struct SeqFile
{
  const char *path;

  // for reading FASTA/FASTQ/plain
  gzFile gz_file;

  // For reading sam/bams
  samFile *sam_file;
  bam1_t *bam;
  bam_hdr_t *sam_header;

  enum SeqFileType file_type;

  // have we seen a '>' at the start of a line in a fasta file?
  // or a '@' in a fastq?
  // For 'plain' format files this is used to store the first char per entry
  char read_line_start;

  // name, index and bases-read/offset of current entry
  StrBuf *entry_name;
  unsigned long entry_index;
  
  unsigned long entry_offset, entry_offset_qual;

  // Whether an entry has been read in
  char entry_read, entry_read_qual;

  // Buffer for reading in bases in FASTQ files
  StrBuf *bases_buff;

  // Total bases read/written - initially 0
  unsigned long total_bases_passed;
  // Total bases skipped (not read through API) in file so far
  unsigned long total_bases_skipped;

  unsigned long line_number;

  /* Writing preferences */

  // Output plain file for writing if not gzipping output
  // (newer zlib allows you to do this with gzFile, but mac version is outdated)
  FILE *plain_file;

  // 0 if no wrap, otherwise max bases per line
  unsigned long line_wrap, curr_line_length;

  // State of writing
  WriteState write_state;
};

// Write output MACROs
// wrapper for fputs/gzputs
#define seq_puts(f,str) (size_t) \
((f)->plain_file != NULL ? (size_t)fputs((str), (f)->plain_file) \
                         : (size_t)gzputs((f)->gz_file, (str)))

// wrapper for fwrite/gzwrite
#define seq_write(f,str,len) (size_t) \
((f)->plain_file != NULL \
  ? (size_t)fwrite((str), sizeof(char), (size_t)(len), (f)->plain_file) \
  : (size_t)gzwrite((f)->gz_file, (str), (unsigned int)(len)))

#define seq_getc(seq) ((seq)->plain_file != NULL ? fgetc((seq)->plain_file) \
                                                 : gzgetc((seq)->gz_file))

#define seq_ungetc(c,seq) ((seq)->plain_file != NULL \
  ? ungetc((c),(seq)->plain_file) \
  : gzungetc((c),(seq)->gz_file))

#define seq_readline(sbuf,seq) ((seq)->plain_file != NULL \
  ? strbuf_readline((sbuf), (seq)->plain_file) \
  : strbuf_gzreadline((sbuf), (seq)->gz_file))

#define seq_skip_line(seq) ((seq)->plain_file != NULL \
  ? strbuf_skip_line((seq)->plain_file) \
  : strbuf_gzskip_line((seq)->gz_file))

#define MIN(x,y) ((x) <= (y) ? (x) : (y))

#define is_base_char(x) ((x) == 'a' || (x) == 'A' || \
                         (x) == 'c' || (x) == 'C' || \
                         (x) == 'g' || (x) == 'G' || \
                         (x) == 't' || (x) == 'T' || \
                         (x) == 'n' || (x) == 'N')


#define _write(s,str,len) \
((sf)->line_wrap == 0 ? seq_puts((sf), (str)) : _write_wrapped((sf),(str),(len)))

size_t _write_wrapped(SeqFile *sf, const char *str, size_t str_len);

#endif
