seq_file
========
C Library for reading multiple bioinformatics sequence file formats  
https://github.com/noporpoise/seq_file  
Isaac Turner turner.isaac@gmail.com  
4 July 2012, GPLv3  

    *NOTE*: new_api/ will soon replace the code and current API. This implements
    buffered and unbuffered reading which provides a significant speed increase
    for reading gzip files on older systems. It also drops the dependency
    for string_buffer and is implemented entirely in two header files.  It does
    not provide all the same features, please get in touch if there's anything
    you'd like to see added to new_api. -Isaac

About
=====

The aim is to provide a C library that allows programs to transparently read
sequence data from multiple file formats, without having to worry about the
format.  Pass a file to seq_file_open() and then read sequences using
seq_next_read() without having to worry about what the file format is.  

Currently supports:
* SAM & BAM
* FASTA (& gzipped fasta)
* FASTQ (& gzipped fastq)
* 'plain' format (one sequence per line [.txt]) (& gzipped plain [.txt.gz])

Also included is a simple example program (seq_convert) that converts between
file formats, for example:

    $ seq_convert in.fq.gz out.fa
    $ seq_convert in.fa out.fq.gz
    $ seq_convert in.bam out.fa.gz
    $ seq_convert in.sam out.txt

It can also 'wrap lines' in the output:

    $ seq_convert in.fq.gz out.fa 80


Example Code
============

Example code to read a file and print as a FASTA file using seq_file.
For more example code see examples/example.c

    SeqFile *sf = seq_file_open(path);

    if(sf == NULL)
    {
      // Error opening file
      fprintf(stderr, "Error: cannot read seq file '%s'\n", path);
      exit(EXIT_FAILURE);
    }

    char c, q;

    while(seq_next_read(sf))
    {
      printf(">%s\n", seq_get_read_name(sf));

      while(seq_read_base(sf, &c))
        putc(c, stdout);

      putc('\n', stdout);
    }

    seq_file_close(sf);

Build
=====

Seq_file requires:

* htslib [https://github.com/samtools/htslib]
* string_buffer [https://github.com/noporpoise/string_buffer]

It also requires zlib, which should already be installed.  

To build the test code and the program seq_convert:

    make STRING_BUF_PATH=path/to/string_buffer/ HTS_PATH=path/to/htslib/

Sometimes the linker can't find your libz.a file (zlib), so you may need to try:

    make STRING_BUF_PATH=path/to/string_buffer/ HTS_PATH=path/to/htslib/ ZLIB_PATH=/dir/with/libz/in/

To call from your own programs, use the following in your Makefile etc.

    LIBS=-lseqfile -lhts -lstrbuf -lz
    INCS=-I$(PATH_TO_SAMTOOLS) -I$(PATH_TO_STRING_BUFFER) -I$(PATH_TO_SEQ_FILE) \
         -L$(PATH_TO_SAMTOOLS) -L$(PATH_TO_STRING_BUFFER) -L$(PATH_TO_SEQ_FILE)
    gcc $(INCS) <your files / args etc.> $(LIBS)


Functions
=========

Opening / Closing
-----------------

Open a file for reading.

    SeqFile* seq_file_open(const char* path);

Open a file for writing.

    SeqFile* seq_file_open_write(const char* file_path, const SeqFileType file_type,
                                 const char gzip, const unsigned long line_wrap);

Close a file. If open for writing: flushes changes and returns number of bytes
written.  

    size_t seq_file_close(SeqFile *sf);

Reading
-------

Get the next sequence read.  Returns 1 on success 0 if no more to read

    char seq_next_read(SeqFile *sf);

Get the name of the current read

    const char* seq_get_read_name(SeqFile *sf);

Get this index of the current read -- starts from 0

    unsigned long seq_get_read_index(SeqFile *sf);

Get the distance into this read that we have read

    unsigned long seq_get_base_offset(SeqFile *sf);

Get the distance into this read's quality scores that we have read

    unsigned long seq_get_qual_offset(SeqFile *sf);

If `seq_next_read()` returned 1 and `seq_read_base()` is now returning 0 (i.e we
have read all the bases), `seq_get_length()` will now report the correct read
length

    unsigned long seq_get_length(SeqFile *sf);

Read a single base from the current read.
Returns 1 on success, 0 if no more bases

    char seq_read_base(SeqFile *sf, char *c);

Read a single quality score from the current read.
Returns 1 on success, 0 if no more quality scores

    char seq_read_qual(SeqFile *sf, char *c);

Read `k` bases/quality scores of the current read into `str`.  `str` must be at
least k+1 bytes long. Returns 1 on success, 0 otherwise

    char seq_read_k_bases(SeqFile *sf, char* str, int k);
    char seq_read_k_quals(SeqFile *sf, char* str, int k);

Read all remaining bases/quality scores from the current read into a
string_buffer. Returns 1 on success, 0 otherwise.

    char seq_read_all_bases(SeqFile *sf, StrBuf *sbuf);
    char seq_read_all_quals(SeqFile *sf, StrBuf *sbuf);

Writing
-------

Check if a file is open for writing. Returns 1 if open for writing, 0 otherwise.

    char seq_is_open_for_write(const SeqFile *sf);

Write various values to the file.  Must be called in order of: name, seq, qual.
Any of `name`, `seq`, `quals` may be NULL.  (These API calls may change).  

    size_t seq_file_write_name(SeqFile *sf, const char *name);
    size_t seq_file_write_seq(SeqFile *sf, const char *seq);
    size_t seq_file_write_qual(SeqFile *sf, const char *qual);

Get file details
----------------

Guess a filetype from path

    void seq_guess_filetype_from_path(const char* path,
                                      SeqFileType* file_type,
                                      char* zipped);

Various methods for getting file type

    SeqFileType seq_file_get_type(const SeqFile *sf);
    const char* seq_file_get_type_str(const SeqFile  *sf);
    const char* seq_file_type_str(const SeqFileType file_type, const char zipped);

Get file path

    const char* seq_get_path(const SeqFile *sf);

Get the number of bases read/written so far

    unsigned long seq_num_bases_passed(const SeqFile *sf);

Whether or not this file has quality scores.  Returns 1 if file has quality
scores, 0 otherwise.  Note: quality scores may still all be set to 'null' e.g.
??? or ### etc.

    char seq_has_quality_scores(const SeqFile *sf);


Get min and max quality values in the first 500 quality scores of a file.
Returns -1 on error, 0 if no quality scores or no reads, 1 on success.

    char seq_estimate_qual_limits(const char *path, char *min, char *max);

Development
===========

Please get in touch if you find bugs, have questions or can suggest features.  
Isaac <turner.isaac@gmail.com>

To run some basic tests:

    for f in test/*; do echo $f >&2; ./seq_file_test $f; done > test/tests.out

TODO:
seq_file
 * Add support for [file] to [.bam|.sam] -- unmapped
 * re-write output writing code
 * Retrieve mate pair data etc from a sam/bam, pairs of fastq/fasta
 * Add pair-end support to seq_convert
 * Add support for sra?

Proposed new API to deal with mate pairs:

    typedef struct SeqRead;

    SeqRead* seq_file_create_read();
    seq_file_free_read(SeqRead *sr);

    // 1)
    char seq_file_next_read(SeqFile *sf, SeqRead *sr);
    char seq_file_next_read_mp(SeqFile *sf, SeqRead *sr1, SeqRead *sr2);

    // OR 2)
    char seq_file_next_template(SeqFile *sf, SeqTemplate *st);
    char seq_file_next_read(SeqTemplate st*, SeqRead* sr);

    // Then we can use the following to read
    char seq_read_base(SeqRead *sr, char *c);
    char seq_read_qual(SeqRead *sr, char *c);
    char seq_read_k_bases(SeqRead *sr, char* str, int k);
    char seq_read_k_quals(SeqRead *sr, char* str, int k);
    char seq_read_all_bases(SeqRead *sr, StrBuf *sbuf);
    char seq_read_all_quals(SeqRead *sr, StrBuf *sbuf);

    // Get position etc
    unsigned long seq_get_base_offset(SeqRead *sr);
    unsigned long seq_get_qual_offset(SeqRead *sr);
    unsigned long seq_get_length(SeqRead *sr);

License
=======

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
 
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.
 
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
