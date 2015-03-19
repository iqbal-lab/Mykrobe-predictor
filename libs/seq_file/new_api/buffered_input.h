
#ifndef _BUFFERED_INPUT_HEADER
#define _BUFFERED_INPUT_HEADER

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

typedef struct
{
  char *b;
  // begin is index of first char (unless begin >= end)
  // end is index of \0
  // size should be >= end+1 to allow for \0
  // (end-size) is the number of bytes in buffer
  size_t begin, end, size;
} buffer_t;

// Buffer functions
#define ROUNDUP2POW(x) (0x1 << (64 - __builtin_clzl(x)))

static inline char buffer_init(buffer_t *b, size_t s)
{
  if((b->b = (char*)malloc(sizeof(char)*s)) == NULL) return 0;
  b->size = s;
  b->begin = b->end = 0;
  return 1;
}

static inline buffer_t* buffer_alloc(size_t s)
{
  buffer_t *b = (buffer_t*)malloc(sizeof(buffer_t));
  if(b == NULL) return NULL;
  else if(buffer_init(b,s)) return b;
  free(b); return NULL; /* couldn't malloc */
}

#define buffer_destroy(buf) do{free((buf)->b); free(buf); } while(0)

// size_t s is the number of bytes you want to be able to store
// the actual buffer is created with s+1 bytes to allow for the \0
static inline void buffer_ensure_capacity(buffer_t *buf, size_t s)
{
  if(buf->size < ++s) {
    buf->size = ROUNDUP2POW(s);
    buf->b = realloc(buf->b, buf->size);
  }
}

static inline void buffer_append_str(buffer_t *buf, char *str)
{
  size_t len = buf->end + strlen(str);
  buffer_ensure_capacity(buf, len);
  memcpy(buf->b+buf->end, str, len);
  buf->b[buf->end = len] = 0;
}

static inline void buffer_append_char(buffer_t *buf, char c)
{
  buffer_ensure_capacity(buf, buf->end+1);
  buf->b[buf->end++] = c;
  buf->b[buf->end] = '\0';
}

#define buffer_terminate(buf) (buf.b[buf.end] = 0)

#define buffer_chomp(buf) do { \
    if(buf.end > 0 && buf.b[buf.end-1] == '\n') {                              \
      buf.end--;                                                               \
      if(buf.end > 0 && buf.b[buf.end-1] == '\r') buf.end--;                   \
      buf.b[buf.end] = 0;                                                      \
    }                                                                          \
  } while(0)

/* 
Unbuffered

fgetc(f)
gzgetc(gz)
gzread2(gz,buf,len)
fread2(f,buf,len)
gzgets2(gz,buf,len)
fgets2(f,buf,len)
gzreadline(gz,out)
freadline(f,out)
*/

// Define read for gzFile and FILE (unbuffered)
#define gzread2(gz,buf,len) gzread(gz,buf,len)
#define fread2(f,buf,len) fread(buf,sizeof(char),len,file)

#define gzgets2(gz,buf,len) gzgets(gz,buf,len)
#define fgets2(f,buf,len) fgets(buf,len,f)

// fgetc(f), gzgetc(gz) are already good to go

// Define readline for gzFile and FILE (unbuffered)
#define _func_readline(name,type_t,__gets) \
  static inline size_t name(type_t file, buffer_t *buf)                        \
  {                                                                            \
    buffer_ensure_capacity(buf, buf->end+1);                                   \
    size_t n, total_read = 0;                                                  \
    while(__gets(file, buf->b+buf->end, buf->size-buf->end))                   \
    {                                                                          \
      n = strlen(buf->b+buf->end);                                             \
      buf->end += n; total_read += n;                                          \
      if(buf->b[buf->end-1] == '\n') return total_read;                        \
      else buf->b = realloc(buf->b, buf->size <<= 1);                          \
    }                                                                          \
    return total_read;                                                         \
  }

_func_readline(gzreadline,gzFile,gzgets2)
_func_readline(freadline,FILE*,fgets2)

/* Buffered */

/*
fgetc_buf(f,in)
gzgetc_buf(gz,in)
gzreadline_buf(gz,in,out)
freadline_buf(f,in,out)
*/

#define READ_BUFFER(file,in,__read) do { \
    int input = __read(file,in->b,in->size);                                   \
    (in)->end = input < 0 ? 0 : input;                                         \
    (in)->begin = 0;                                                           \
  } while(0)

// Define getc for gzFile and FILE (buffered)
#define _func_getc_buf(fname,type_t,__read)                                    \
  static inline int fname(type_t file, buffer_t *in)                           \
  {                                                                            \
    if(in->begin >= in->end) {                                                 \
      READ_BUFFER(file,in,__read);                                             \
      return in->end == 0 ? -1 : in->b[in->begin++];                           \
    }                                                                          \
    return in->b[in->begin++];                                                 \
  }

_func_getc_buf(gzgetc_buf,gzFile,gzread2)
_func_getc_buf(fgetc_buf,FILE*,fread2)

// Define readline for gzFile and FILE (buffered)
#define _func_readline_buf(fname,type_t,__read) \
  static inline size_t fname(type_t file, buffer_t *in, buffer_t *buf)         \
  {                                                                            \
    if(in->begin >= in->end) { READ_BUFFER(file,in,__read); }                  \
    size_t new_buf_len, total_read = 0;                                        \
    while(in->end != 0)                                                        \
    {                                                                          \
      size_t offset = in->begin;                                               \
      while(offset < in->end) { offset++; if(in->b[offset-1] == '\n') break; } \
      offset -= in->begin;                                                     \
      new_buf_len = buf->end+offset;                                           \
      buffer_ensure_capacity(buf, new_buf_len);                                \
      memcpy(buf->b+buf->end, in->b+in->begin, offset);                        \
      buf->end += offset;                                                      \
      in->begin += offset;                                                     \
      total_read += offset;                                                    \
      if(in->begin < in->end) break;                                           \
      READ_BUFFER(file,in,__read);                                             \
    }                                                                          \
    buf->b[buf->end] = 0;                                                      \
    return total_read;                                                         \
  }

_func_readline_buf(gzreadline_buf,gzFile,gzread2)
_func_readline_buf(freadline_buf,FILE*,fread2)

#define BUFFERED_INPUT_SETUP()

#endif
