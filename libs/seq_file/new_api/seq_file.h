
#ifndef _SEQ_FILE_HEADER
#define _SEQ_FILE_HEADER

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <zlib.h>

#include "hts.h"
#include "sam.h"
#include "buffered_input.h"

typedef struct seq_file_t seq_file_t;
typedef struct read_t read_t;

struct seq_file_t
{
  char *path;
  FILE *f_file;
  gzFile gz_file;
  samFile *s_file;
  bam1_t *bam;
  bam_hdr_t *bam_header;
  int (*read)(seq_file_t *sf, read_t *r);
  buffer_t in;
  int headc, nextc;
};

struct read_t
{
  buffer_t name, seq, qual;
};

#define seq_is_bam(sf) ((sf)->s_file != NULL && (sf)->s_file->is_bin)
#define seq_is_sam(sf) ((sf)->s_file != NULL && !(sf)->s_file->is_bin)
#define seq_use_gzip(sf) ((sf)->gz_file != NULL)

// The following require a read to have been read successfully using seq_read
#define seq_is_fasta(sf) ((sf)->s_file == NULL && (sf)->headc == '>')
#define seq_is_fastq(sf) ((sf)->s_file == NULL && (sf)->headc == '@')
#define seq_is_plain(sf) ((sf)->s_file == NULL && (sf)->headc != '@' && (sf)->headc != '>')

#define seq_get_path(sf) ((sf)->path)

// return 1 on success, 0 on eof, -1 if partially read / syntax error
#define seq_read(sf,r) (sf)->read(sf,r)

// fh could be sam,FASTA,FASTQ,txt (+gzip)
// file could be sam,bam,FASTA,FASTQ,txt (+gzip)

#define DEFAULT_BUFSIZE (1<<20)

static const int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 9, 14,
                                           1, 6, 5, 13, 3, 11, 7, 15 };

static inline int sread_s(seq_file_t *sf, read_t *read)
{
  read->name.end = read->seq.end = read->qual.end = 0;
  buffer_terminate(read->name);
  buffer_terminate(read->seq);
  buffer_terminate(read->qual);

  if(sam_read1(sf->s_file, sf->bam_header, sf->bam) < 0) return 0;

  char *str = bam_get_qname(sf->bam);
  buffer_append_str(&read->name, str);

  size_t qlen = sf->bam->core.l_qseq;
  buffer_ensure_capacity(&read->seq, qlen);
  buffer_ensure_capacity(&read->qual, qlen);
  uint8_t *bamseq = bam_get_seq(sf->bam);
  uint8_t *bamqual = bam_get_qual(sf->bam);

  size_t i, j;
  if(bam_is_rev(sf->bam))
  {
    for(i = 0, j = qlen - 1; i < qlen; i++, j--)
    {
      int8_t b = bam_seqi(bamseq, j);
      read->seq.b[i] = seq_nt16_str[seq_comp_table[b]];
      read->qual.b[i] = 33 + bamqual[j];
    }
  }
  else
  {
    for(i = 0, j = qlen - 1; i < qlen; i++, j--)
    {
      int8_t b = bam_seqi(bamseq, i);
      read->seq.b[i] = seq_nt16_str[b];
      read->qual.b[i] = 33 + bamqual[i];
    }
  }

  read->seq.end = qlen;
  read->qual.end = qlen;
  return 1;
}

#define _func_read(fname,__getc,__readline) \
  static inline int fname(seq_file_t *sf, read_t *read)                        \
  {                                                                            \
    read->name.end = read->seq.end = read->qual.end = 0;                       \
    buffer_terminate(read->name);                                              \
    buffer_terminate(read->seq);                                               \
    buffer_terminate(read->qual);                                              \
                                                                               \
    int c;                                                                     \
    sf->headc = (sf->nextc != -1 ? sf->nextc : __getc(sf));                    \
    sf->nextc = -1;                                                            \
    if(sf->headc == -1) return 0;                                              \
    else if(sf->headc == '@')                                                  \
    {                                                                          \
      if(__readline(sf, read->name) == 0) return -1;                           \
      buffer_chomp(read->name);                                                \
                                                                               \
      while((c = __getc(sf)) != '+') {                                         \
        if(c == -1) return -1;                                                 \
        if(c != '\r' && c != '\n') {                                           \
          buffer_append_char(&read->seq,c);                                    \
          if(__readline(sf, read->seq) == 0) return -1;                        \
          buffer_chomp(read->seq);                                             \
        }                                                                      \
      }                                                                        \
      while((c = __getc(sf)) != -1 && c != '\n');                              \
      if(c == -1) return -1;                                                   \
      do {                                                                     \
        if(__readline(sf,read->qual) > 0) buffer_chomp(read->qual);            \
        else return 1;                                                         \
      } while(read->qual.end < read->seq.end);                                 \
    }                                                                          \
    else if(sf->headc == '>')                                                  \
    {                                                                          \
      if(__readline(sf, read->name) == 0) return -1;                           \
      buffer_chomp(read->name);                                                \
                                                                               \
      while((c = __getc(sf)) != '>') {                                         \
        if(c == -1) return 1;                                                  \
        if(c != '\r' && c != '\n') {                                           \
          buffer_append_char(&read->seq,c);                                    \
          size_t r = __readline(sf, read->seq);                                \
          buffer_chomp(read->seq);                                             \
          if(r == 0) return 1;                                                 \
        }                                                                      \
      }                                                                        \
                                                                               \
      sf->nextc = c;                                                           \
    }                                                                          \
    else {                                                                     \
      buffer_append_char(&read->seq, sf->headc);                               \
      __readline(sf, read->seq);                                               \
      buffer_chomp(read->seq);                                                 \
    }                                                                          \
    return 1;                                                                  \
  }

// perform reading on seq_file_t
#define _sf_gzgetc(sf)              gzgetc(sf->gz_file)
#define _sf_gzgetc_buf(sf)          gzgetc_buf(sf->gz_file,&sf->in)
#define _sf_fgetc(sf)               fgetc(sf->f_file)
#define _sf_fgetc_buf(sf)           fgetc_buf(sf->f_file,&sf->in)

#define _sf_gzreadline(sf,buf)      gzreadline(sf->gz_file,&buf)
#define _sf_gzreadline_buf(sf,buf)  gzreadline_buf(sf->gz_file,&sf->in,&buf)
#define _sf_freadline(sf,buf)       freadline(sf->f_file,&buf)
#define _sf_freadline_buf(sf,buf)   freadline_buf(sf->f_file,&sf->in,&buf)

_func_read(sread_gz_buf, _sf_gzgetc_buf, _sf_gzreadline_buf)
_func_read(sread_gz,     _sf_gzgetc,     _sf_gzreadline)
_func_read(sread_f_buf,  _sf_fgetc_buf,  _sf_freadline_buf)
_func_read(sread_f,      _sf_fgetc,      _sf_freadline)

// Create and destroy read structs
static inline read_t* seq_read_alloc()
{
  read_t *r = calloc(1, sizeof(read_t));
  buffer_init(&r->name, 512);
  buffer_init(&r->seq, DEFAULT_BUFSIZE);
  buffer_init(&r->qual, DEFAULT_BUFSIZE);
  return r;
}

static inline void seq_read_destroy(read_t *r)
{
  free(r->name.b);
  free(r->seq.b);
  free(r->qual.b);
  free(r);
}

#define seq_file_init(sf) do { \
    sf->gz_file = NULL;                                                        \
    sf->f_file = NULL;                                                         \
    sf->s_file = NULL;                                                         \
    sf->headc = sf->nextc = -1;                                                \
    sf->in.size = sf->in.begin = sf->in.end = 0;                               \
    sf->path = sf->in.b = NULL;                                                \
  } while(0)

// I have removed gzbuffer for now since it causes linking errors on systems
// that are not set up properly.  Feel free to uncomment and remove empty
// definition
#if defined(ZLIB_VERNUM) && ZLIB_VERNUM >= 0x1240
//#define SET_ZLIB_BUFFER(sf,s) gzbuffer((sf)->gz_file, DEFAULT_BUFSIZE)
#define SET_ZLIB_BUFFER(sf,s)
#else
#define SET_ZLIB_BUFFER(sf,s)
#endif

#define seq_setup(sf,use_zlib,buf_size) do \
  {                                                                            \
    sf->in.size = buf_size;                                                    \
    if(sf->in.size > 0 && (sf->in.b = (char*)malloc(sf->in.size)) == NULL) {   \
      free(sf); return NULL;                                                   \
    }                                                                          \
    if((sf)->in.size == 0 && use_zlib) { SET_ZLIB_BUFFER(sf,DEFAULT_BUFSIZE); }\
    if(use_zlib) { sf->read = sf->in.size > 0 ? sread_gz_buf : sread_gz; }     \
    else         { sf->read = sf->in.size > 0 ? sread_f_buf  : sread_f; }      \
  } while(0)

// Can only be used for plain,fasta,fastq (+gzip) formats
static inline seq_file_t* seq_open_fh(FILE *fh, char use_zlib, size_t buf_size)
{
  seq_file_t *sf = (seq_file_t*)calloc(1, sizeof(seq_file_t));
  seq_file_init(sf);
  if(buf_size == 0) sf->f_file = fh;
  else if((sf->gz_file = gzdopen(fileno(fh), "r")) == NULL) {
    free(sf);
    return NULL;
  }
  seq_setup(sf, use_zlib, buf_size);
  sf->path = strdup("-");
  return sf;
}

static inline seq_file_t* seq_open2(const char *p, char sam,
                                    char use_zlib, size_t buf_size)
{
  seq_file_t *sf = (seq_file_t*)calloc(1, sizeof(seq_file_t));
  seq_file_init(sf);

  if(sam == 1 || sam == 2)
  {
    if((sf->s_file = sam_open(p, sam == 1 ? "rs" : "rb", 0)) == NULL) {
      free(sf);
      return NULL;
    }
    sf->bam = bam_init1();
    sf->bam_header = sam_hdr_read(sf->s_file);
    sf->read = sread_s;
  }
  else
  {
    if((use_zlib && ((sf->gz_file = gzopen(p, "r")) == NULL)) ||
       (!use_zlib && ((sf->f_file = fopen(p, "r")) == NULL)))
    {
      free(sf);
      return NULL;
    }
    seq_setup(sf, use_zlib, buf_size);
  }
  sf->path = strdup(p);
  return sf;
}

// Guess file type from file path or contents
#define EXT_ARRLEN 21
typedef enum { IS_ERROR, IS_UNKNOWN,
               IS_SEQ, IS_SEQ_GZIP,
               IS_SAM, IS_BAM } filetype_t;

static inline filetype_t _guess_filetype_from_path(const char *path)
{
  size_t plen = strlen(path);
  const char* exts[EXT_ARRLEN]
    = {".fa", ".fasta", ".fq", ".fastq", ".txt", ".fz",
       ".faz", ".fagz", ".fa.gz", ".fa.gzip", ".fasta.gzip", ".fqz", ".fqgz",\
       ".fq.gz", ".fq.gzip", ".fastq.gzip", ".txtgz", ".txt.gz", ".txt.gzip",\
       ".sam", ".bam"};
  const filetype_t types[EXT_ARRLEN]
    = {IS_SEQ, IS_SEQ, IS_SEQ, IS_SEQ, IS_SEQ, IS_SEQ_GZIP,
       IS_SEQ_GZIP, IS_SEQ_GZIP, IS_SEQ_GZIP, IS_SEQ_GZIP, IS_SEQ_GZIP,
       IS_SEQ_GZIP, IS_SEQ_GZIP, IS_SEQ_GZIP, IS_SEQ_GZIP, IS_SEQ_GZIP,
       IS_SEQ_GZIP, IS_SEQ_GZIP, IS_SEQ_GZIP, IS_SAM, IS_BAM};

  size_t extlens[EXT_ARRLEN];
  size_t i;
  for(i = 0; i < EXT_ARRLEN; i++)
    extlens[i] = strlen(exts[i]);

  for(i = 0; i < EXT_ARRLEN; i++)
    if(extlens[i] <= plen && strcasecmp(path+plen-extlens[i], exts[i]) == 0)
      return types[i];

  return IS_UNKNOWN;
}

static inline filetype_t _guess_filetype_from_content(const char *path)
{
  char buf[20];
  gzFile gz = gzopen(path, "r");
  if(gz == NULL) return IS_ERROR;
  int read = gzread(gz,buf,20);
  gzclose(gz);

  if(read >= 4 && strncmp(buf, "@HD\t", 4) == 0) return IS_SAM;
  if(read >= 3 && strncmp(buf, "BAM", 3) == 0) return IS_BAM;
  if(strchr(buf,'\t') != NULL) return IS_SAM;
  return strchr("@>acgtn", tolower(buf[0])) != NULL ? IS_SEQ : IS_UNKNOWN;
}

static inline seq_file_t* seq_open(const char *p)
{
  if(strcmp(p,"-") == 0) return seq_open_fh(stdin, 0, 0);
  filetype_t type = _guess_filetype_from_path(p);
  if(type == IS_UNKNOWN) {
    type = _guess_filetype_from_content(p);
    if(type == IS_UNKNOWN) {
      fprintf(stderr, "Unknown filetype: %s\n", p);
      fprintf(stderr, "Email turner.isaac@gmail.com to report a bug\n");
      exit(EXIT_FAILURE);
    }
    else if(type == IS_ERROR) return NULL;
  }
  char sam_bam = 0;
  char zipped = 1;
  size_t buf_size = (size_t)DEFAULT_BUFSIZE;
  if(type == IS_SAM) sam_bam = 1;
  else if(type == IS_BAM) sam_bam = 2;
  else if(type == IS_SEQ) zipped = buf_size = 0;
  return seq_open2(p, sam_bam, zipped, buf_size);
}

static inline void seq_close(seq_file_t *sf)
{
  if(sf->f_file != NULL) fclose(sf->f_file);
  else if(sf->gz_file != NULL) gzclose(sf->gz_file);
  else if(sf->s_file != NULL)
  {
    sam_close(sf->s_file);
    free(sf->bam);
    free(sf->bam_header);
  }
  if(sf->in.size != 0) free(sf->in.b);
}

// Get min and max quality values in the first `num` bases of a file.
// Returns 0 if no qual scores, 1 on success, -1 if read error
static inline int seq_get_qual_limits(const char *path, size_t num,
                                      int *minptr, int *maxptr)
{
  seq_file_t *sf = seq_open(path);
  if(sf == NULL) return -1;
  read_t *r = seq_read_alloc();
  int i, limit, min = INT_MAX, max = 0;
  size_t count = 0;

  while(count < num && seq_read(sf,r))
  {
    limit = (r->qual.end < num-count ? r->qual.end : num-count);
    for(i = 0; i < limit; i++)
    {
      char q = r->qual.b[i];
      if(q > max) max = q;
      if(q < min) min = q;
    }
    count += limit;
  }

  seq_close(sf);

  if(count > 0) {
    *minptr = min;
    *maxptr = max;
  }
  return (count > 0);
}

static inline char _seq_read_looks_valid(read_t *r, const char *alphabet)
{
  size_t i;
  if(r->qual.end != 0) {
    if(r->qual.end != r->seq.end) return 0;
    for(i = 0; i < r->seq.end; i++) {
      char b = tolower(r->seq.b[i]);
      char q = r->qual.b[i];
      if(strchr(alphabet, b) == NULL) return 0;
      if(q < 33 || q > 104) return 0;
    }
  }
  else {
    for(i = 0; i < r->seq.end; i++) {
      char b = tolower(r->seq.b[i]);
      if(strchr(alphabet, b) == NULL) return 0;
    }
  }
  return 1;
}

#define seq_read_looks_valid_dna(r) _seq_read_looks_valid(r,"acgtn")
#define seq_read_looks_valid_rna(r) _seq_read_looks_valid(r,"acgun")
#define seq_read_looks_valid_protein(r) _seq_read_looks_valid(r,"acdefghiklmnopqrstuvwy")

#define SETUP_SEQ_FILE()

// read_t* seq_read_alloc()
// seq_read_destroy(read_t* r)

// seq_open(path)
// seq_open2(path,sam_bam,use_gzip,buffer_size)
// seq_open_fh(fh,use_gzip,buffer_size)
// seq_close(seq_file_t *sf)

#endif
