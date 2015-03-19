/*
 seq_file.h
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

#ifndef SEQ_HEADER_SEEN
#define SEQ_HEADER_SEEN

#include "string_buffer.h"

typedef struct SeqFile SeqFile;
typedef enum SeqFileType SeqFileType;

enum SeqFileType
{
  SEQ_UNKNOWN = 0, SEQ_FASTA = 1, SEQ_FASTQ = 2, SEQ_PLAIN = 3,
  SEQ_SAM = 4, SEQ_BAM = 5,
};

// Open for reading
SeqFile* seq_file_open(const char* path);

// Open a file assuming a given filetype
//SeqFile* seq_file_open_filetype(const char* file_path,
//                                SeqFileType file_type);

// If open for writing: writes newline to file and returns 1 on success
// otherwise return 0
size_t seq_file_close(SeqFile *sf);

// Guess a filetype from path
void seq_guess_filetype_from_path(const char* path,
                                  SeqFileType* file_type,
                                  char* zipped);

// Various methods for getting file type
SeqFileType seq_file_get_type(const SeqFile *sf);
const char* seq_file_get_type_str(const SeqFile  *sf);
const char* seq_file_type_str(SeqFileType file_type, char zipped);

// Get path
const char* seq_get_path(const SeqFile *sf);

// Get min and max quality values in the first `num` quality scores
// Returns -1 on error, 0 if no quality scores or no reads, 1 on success
int seq_estimate_qual_limits(const char *path, int num, int *min, int *max);

// Get the number of bases read/written so far
unsigned long seq_total_bases_passed(const SeqFile *sf);

// Get the total bases skipped (not read through API) in file so far
unsigned long seq_total_bases_skipped(const SeqFile *sf);

// Get current line number
unsigned long seq_curr_line_number(const SeqFile *sf);

// Returns 1 if file has quality scores, 0 otherwise
// Note: quality scores may still all be set to 'null' e.g. ??? or ### etc.
char seq_has_quality_scores(const SeqFile *sf);

// Returns 1 on success 0 if no more to read
char seq_next_read(SeqFile *sf);

// Get the name of the next read
const char* seq_get_read_name(SeqFile *sf);

// Get this read index -- 0 if no reads read yet,
// otherwise read number of prev read entry
unsigned long seq_get_read_index(SeqFile *sf);

// Get the distance into this read that we have read
unsigned long seq_get_bases_read(SeqFile *sf);
// Get the distance into this read's quality scores that we have read
unsigned long seq_get_quals_read(SeqFile *sf);

// If seq_next_read() returned 1 and seq_read_base() is now returning 0,
// seq_get_length() will now report the correct read length
unsigned long seq_get_length(SeqFile *sf);

// Read a single base from the current read
// Returns 1 on success, 0 if no more quality scores or run out of bases
char seq_read_base(SeqFile *sf, char *c);

// Read a single quality score from the current read
// Returns 1 on success, 0 if no more quality scores or run out of bases
char seq_read_qual(SeqFile *sf, char *c);

// str must be at least k+1 bytes long
// returns 1 on success, 0 otherwise
char seq_read_k_bases(SeqFile *sf, char* str, int k);
char seq_read_k_quals(SeqFile *sf, char* str, int k);

// returns 1 on success, 0 otherwise
char seq_read_all_bases(SeqFile *sf, StrBuf *sbuf);
char seq_read_all_quals(SeqFile *sf, StrBuf *sbuf);
char seq_read_all_bases_and_quals(SeqFile *sf, StrBuf *sbuf_seq, StrBuf *sbuf_qual);

// Write to a file.  Any of name, seq, quals may be NULL
// Returns the number of bytes written or 0 on failure

SeqFile* seq_file_open_write(const char* file_path, SeqFileType file_type,
                             char gzip, unsigned long line_wrap);

// Returns 1 if open for writing, 0 otherwise
char seq_is_open_for_write(const SeqFile *sf);

size_t seq_file_write_name(SeqFile *sf, const char *name);
size_t seq_file_write_seq(SeqFile *sf, const char *seq);
size_t seq_file_write_qual(SeqFile *sf, const char *qual);

#endif
