/*
 seq_sam.c
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

#include "hts.h"

#include "seq_sam.h"

// Array for complementing bases read from BAM/SAM files
static int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 9, 14,
                                     1, 6, 5, 13, 3, 11, 7, 15 };

char seq_next_read_sam(SeqFile *sf)
{
  if(sf->entry_read)
  {
    // Count skipped bases
    sf->total_bases_skipped += (unsigned long)sf->bam->core.l_qseq -
                               sf->entry_offset;
  }

  if(sam_read1(sf->sam_file, sf->sam_header, sf->bam) < 0)
    return 0;

  // Get name
  strbuf_append_str(sf->entry_name, bam_get_qname(sf->bam));

  return 1;
}

char seq_read_base_sam(SeqFile *sf, char *c)
{
  unsigned long query_len = (unsigned long)sf->bam->core.l_qseq;

  if(sf->entry_offset >= query_len)
    return 0;

  uint8_t *seq = bam_get_seq(sf->bam);

  if(bam_is_rev(sf->bam))
  {
    unsigned long index = query_len - sf->entry_offset - 1;
    int8_t b = bam_seqi(seq, index);
    *c = seq_nt16_str[seq_comp_table[b]];
  }
  else
  {
    int8_t b = bam_seqi(seq, sf->entry_offset);
    *c = seq_nt16_str[b];
  }

  return 1;
}

char seq_read_qual_sam(SeqFile *sf, char *c)
{
  unsigned long query_len = (unsigned long)sf->bam->core.l_qseq;

  if(sf->entry_offset_qual >= query_len)
    return 0;

  uint8_t *seq = bam_get_qual(sf->bam);
  unsigned long index;

  if(bam_is_rev(sf->bam))
    index = query_len - sf->entry_offset_qual - 1;
  else
    index = sf->entry_offset_qual;

  *c = 33 + seq[index];

  return 1;
}

char seq_read_all_bases_sam(SeqFile *sf, StrBuf *sbuf)
{
  unsigned long qlen = (unsigned long)sf->bam->core.l_qseq;

  if(sf->entry_offset >= qlen)
    return 0;

  // Get reverse
  char is_reversed = bam_is_rev(sf->bam);

  strbuf_ensure_capacity(sbuf, qlen - sf->entry_offset);

  uint8_t *seq = bam_get_seq(sf->bam);

  // read in and reverse complement (if needed)
  unsigned long i;
  for(i = sf->entry_offset; i < qlen; i++)
  {
    unsigned long index = (is_reversed ? qlen - i - 1 : i);
    int8_t b = bam_seqi(seq, index);

    char c = seq_nt16_str[is_reversed ? seq_comp_table[b] : b];
    strbuf_append_char(sbuf, c);
  }

  return 1;
}

char seq_read_all_quals_sam(SeqFile *sf, StrBuf *sbuf)
{
  unsigned long qlen = (unsigned long)sf->bam->core.l_qseq;

  if(sf->entry_offset_qual >= qlen)
    return 0;

  // Get reverse
  char is_reversed = bam_is_rev(sf->bam);

  strbuf_ensure_capacity(sbuf, qlen - sf->entry_offset);

  uint8_t *seq = bam_get_qual(sf->bam);

  // read in and reverse complement (if needed)
  unsigned long i;
  for(i = sf->entry_offset; i < qlen; i++)
  {
    char c = 33 + seq[is_reversed ? qlen - i - 1 : i];
    strbuf_append_char(sbuf, c);
  }

  return 1;
}


char seq_read_all_bases_and_quals_sam(SeqFile *sf, StrBuf *sbuf_seq, StrBuf *sbuf_qual)
{
  unsigned long qlen = (unsigned long)sf->bam->core.l_qseq;

  if(sf->entry_offset >= qlen)
    return 0;

  // Get reverse
  char is_reversed = bam_is_rev(sf->bam);

  strbuf_ensure_capacity(sbuf_seq, qlen - sf->entry_offset);
  strbuf_ensure_capacity(sbuf_qual, qlen - sf->entry_offset);

  uint8_t *seq = bam_get_seq(sf->bam);
  uint8_t *qual = bam_get_qual(sf->bam);

  // read in and reverse complement (if needed)
  unsigned long i;
  for(i = sf->entry_offset; i < qlen; i++)
  {
    unsigned long index = (is_reversed ? qlen - i - 1 : i);
    int8_t b = bam_seqi(seq, index);
    char c = seq_nt16_str[is_reversed ? seq_comp_table[b] : b];
    strbuf_append_char(sbuf_seq, c);

    char cq = 33 + qual[is_reversed ? qlen - i - 1 : i];
    strbuf_append_char(sbuf_qual, cq);
  }

  return 1;
}
