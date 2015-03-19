/*
 string_buffer.c
 project: string_buffer
 url: https://github.com/noporpoise/StringBuffer
 author: Isaac Turner <turner.isaac@gmail.com>

 Copyright (c) 2011, Isaac Turner
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN_SIZE 16

#include <stdlib.h>
#include <string.h>
#include <ctype.h> // toupper() and tolower()

#include "string_buffer.h"

/******************************/
/*  Constructors/Destructors  */
/******************************/

// Return default length string buffer
StrBuf* strbuf_new()
{
  return strbuf_init(256);
}

StrBuf* strbuf_init(t_buf_pos size)
{
  t_buf_pos new_size = size < MIN_SIZE ? MIN_SIZE : size;

  StrBuf* sbuf = (StrBuf*) malloc(sizeof(StrBuf));

  sbuf->buff = (char*) malloc(new_size);
  sbuf->len = 0;
  sbuf->size = new_size;

  if(sbuf->buff == NULL)
  {
    fprintf(stderr, "Error: StrBuf couldn't be created with %lui bytes.",
            new_size);
    exit(-1);
  }

  return sbuf;
}

StrBuf* strbuf_create(const char* str)
{
  t_buf_pos str_len = strlen(str);

  StrBuf* sbuf = strbuf_init(str_len+1);

  strcpy(sbuf->buff, str);
  sbuf->buff[str_len] = '\0';

  sbuf->len = str_len;

  return sbuf;
}

void strbuf_reset(StrBuf* sbuf)
{
  sbuf->len = 0;
  sbuf->buff[0] = '\0';
}

void strbuf_free(StrBuf* sbuf)
{
  free(sbuf->buff);
  free(sbuf);
}

// Free sbuf struct, but retain and return the char array
char* strbuf_free_get_str(StrBuf* sbuf)
{
  char *buff = sbuf->buff;
  free(sbuf->buff);
  return buff;
}

StrBuf* strbuf_clone(const StrBuf* sbuf)
{
  // One byte for the string end / null char \0
  StrBuf* sbuf_cpy = strbuf_init(sbuf->len+1);
  
  strcpy(sbuf_cpy->buff, sbuf->buff);
  sbuf_cpy->buff[sbuf->len] = '\0';

  sbuf_cpy->len = sbuf->len;
  
  return sbuf_cpy;
}

// Get a copy of this StrBuf as a char array
// Returns NULL if not enough memory
char* strbuf_as_str(const StrBuf* sbuf)
{
  char* cpy = (char*) malloc(sbuf->len+1);

  if(cpy == NULL)
  {
    return NULL;
  }

  memcpy(cpy, sbuf->buff, sbuf->len);
  cpy[sbuf->len] = '\0';

  return cpy;
}

// Get string length
t_buf_pos strbuf_len(const StrBuf *sbuf)
{
  return sbuf->len;
}

// Get buffer length
t_buf_pos strbuf_size(const StrBuf *sbuf)
{
  return sbuf->size;
}

void strbuf_update_len(StrBuf *sbuf)
{
  sbuf->len = strlen(sbuf->buff);
}

// Get / set characters

char strbuf_get_char(const StrBuf *sbuf, t_buf_pos index)
{
  // Bounds checking
  if(index >= sbuf->len)
  {
    fprintf(stderr, "StrBuf OutOfBounds Error: "
                    "strbuf_get_char(index: %lu) [strlen: %lu]\n",
            (unsigned long)index, (unsigned long)sbuf->len);

    return -1;
  }

  return sbuf->buff[index];
}

void strbuf_set_char(StrBuf *sbuf, t_buf_pos index, char c)
{
  // Bounds checking
  if(index > sbuf->len)
  {
    fprintf(stderr, "StrBuf OutOfBounds Error: "
                    "strbuf_set_char(index: %lu, %c) [strlen: %lu]\n",
            (unsigned long)index, c, (unsigned long)sbuf->len);

    return;
  }
  else if(index == sbuf->len)
  {
    // Extend
    strbuf_ensure_capacity(sbuf, sbuf->len + 1);
    sbuf->buff[sbuf->len++] = c;
    sbuf->buff[sbuf->len] = '\0';
  }
  else
  {
    sbuf->buff[index] = c;
  }
}

// Set string buffer to contain a given string
// The string can be a string within the given string buffer
void strbuf_set(StrBuf *sbuf, const char *str)
{
  size_t len = strlen(str);
  strbuf_ensure_capacity(sbuf, len+1);

  // Use memmove to allow overlapping strings
  memmove(sbuf->buff, str, sizeof(char) * len);

  sbuf->buff[len] = '\0';
  sbuf->len = len;
}

/******************************/
/*  Resize Buffer Functions   */
/******************************/

char strbuf_resize(StrBuf *sbuf, t_buf_pos new_size)
{
  char *new_buff = realloc(sbuf->buff, new_size);

  if(new_buff == NULL)
  {
    return 0;
  }

  sbuf->buff = new_buff;
  sbuf->size = new_size;

  if(sbuf->len+1 >= sbuf->size)
  {
    // Buffer was shrunk - add null byte
    sbuf->len = sbuf->size-1;
    sbuf->buff[sbuf->len] = '\0';
  }

  return 1;
}

void strbuf_resize_vital(StrBuf *sbuf, t_buf_pos new_size)
{
  if(!strbuf_resize(sbuf, new_size))
  {
    fprintf(stderr, "Error: StrBuf couldn't be given more memory.  "
                    "Requested %lui bytes.  StrBuf begins '%s...'",
            new_size, strbuf_substr(sbuf, 0, 5));
    
    free(sbuf->buff);
    
    exit(EXIT_FAILURE);
  }
}

// Ensure capacity for len characters plus '\0' character
void strbuf_ensure_capacity(StrBuf *sbuf, t_buf_pos len)
{
  if(sbuf->size > len+1)
  {
    return;
  }

  // Need to resize
  t_buf_pos new_size = 2*sbuf->size;

  while(len+1 >= new_size)
  {
    new_size = 2*new_size;
  }

  strbuf_resize_vital(sbuf, new_size);
}

void strbuf_shrink(StrBuf *sbuf, t_buf_pos new_len)
{
  sbuf->len = new_len;
  sbuf->buff[new_len] = '\0';
}

/********************/
/* String functions */
/********************/

void strbuf_append_str(StrBuf* sbuf, const char* txt)
{
  size_t str_len = strlen(txt);
  strbuf_append_strn(sbuf, txt, str_len);
}

void strbuf_append_strn(StrBuf* sbuf, const char* txt, t_buf_pos len)
{
  // plus 1 for '\0'
  strbuf_ensure_capacity(sbuf, sbuf->len + len);

  strncpy(sbuf->buff+sbuf->len, txt, len);
  sbuf->len += len;

  sbuf->buff[sbuf->len] = '\0';
}

void strbuf_append_char(StrBuf* sbuf, char c)
{
  // Adding one character
  strbuf_ensure_capacity(sbuf, sbuf->len + 1);

  sbuf->buff[(sbuf->len)++] = c;
  sbuf->buff[sbuf->len] = '\0';
}

// Copy a StrBuf to the end of this StrBuf
void strbuf_append_buff(StrBuf* dst, StrBuf* src)
{
  strbuf_ensure_capacity(dst, dst->len + src->len);

  memcpy(dst->buff + dst->len, src->buff, src->len);

  dst->len += src->len;
  dst->buff[dst->len] = '\0';
}


// Remove \r and \n characters from the end of this StrBuf
// Returns the number of characters removed
t_buf_pos strbuf_chomp(StrBuf *sbuf)
{
  t_buf_pos original_length = sbuf->len;

  while(sbuf->len >= 1)
  {
    char last_char = sbuf->buff[sbuf->len-1];

    if(last_char == '\n' || last_char == '\r')
    {
      sbuf->buff[--(sbuf->len)] = '\0';
    }
    else
    {
      break;
    }
  }

  return (original_length - sbuf->len);
}

// Reverse a string
void strbuf_reverse(StrBuf *sbuf)
{
  t_buf_pos half_way = sbuf->len >> 1;
  t_buf_pos i;

  for (i = 0; i < half_way; i++)
  {
    char tmp = sbuf->buff[sbuf->len - 1 - i];
    sbuf->buff[sbuf->len - 1 - i] = sbuf->buff[i];
    sbuf->buff[i] = tmp;
  }
}

// Reverse a string region
void strbuf_reverse_region(StrBuf *sbuf, t_buf_pos start, t_buf_pos length)
{
  t_buf_pos left = start;
  t_buf_pos right = start + length - 1;

  // bounds check
  if(right >= sbuf->len)
  {
    fprintf(stderr, "StrBuf OutOfBounds Error: "
                    "strbuf_reverse_region(start: %lui, len: %lui) "
                    "[strlen: %lu]\n",
            (unsigned long)start, (unsigned long)length,
            (unsigned long)sbuf->len);
    return;
  }

  for(; left < right; left++, right--)
  {
    char tmp = sbuf->buff[left];
    sbuf->buff[left] = sbuf->buff[right];
    sbuf->buff[right] = tmp;
  }
}

char* strbuf_substr(StrBuf *sbuf, t_buf_pos start, t_buf_pos len)
{
  // Bounds checking
  //commented out by Zam  if(start + len >= sbuf->len)
  if(start + len > sbuf->len)
  {
    fprintf(stderr, "StrBuf OutOfBounds Error: "
                    "strbuf_substr(start: %lui, len: %lui) [strlen: %lu]\n",
            (unsigned long)start, (unsigned long)len, (unsigned long)sbuf->len);

    return NULL;
  }

  char* new_string = (char*) malloc(len+1);
  strncpy(new_string, sbuf->buff + start, len);
  new_string[len] = '\0';

  return new_string;
}

void strbuf_substr_prealloced(StrBuf *sbuf, t_buf_pos start, t_buf_pos len, char* outstr)
{
  strncpy(outstr, sbuf->buff + start, len);
  outstr[len] = '\0';
}

void strbuf_to_uppercase(StrBuf *sbuf)
{
  char* pos;
  char* end = sbuf->buff + sbuf->len;

  for(pos = sbuf->buff; pos < end; pos++)
  {
    *pos = toupper(*pos);
  }
}

void strbuf_to_lowercase(StrBuf *sbuf)
{
  char* pos;
  char* end = sbuf->buff + sbuf->len;

  for(pos = sbuf->buff; pos < end; pos++)
  {
    *pos = tolower(*pos);
  }
}

// Get
void strbuf_copy(StrBuf* dst, t_buf_pos dst_pos,
                 const StrBuf* src, t_buf_pos src_pos,
                 t_buf_pos len)
{
  if(dst_pos > dst->len)
  {
    fprintf(stderr, "StrBuf OutOfBounds Error: "
                    "strbuf_copy [dst; pos: %lu; len: %lu; strlen: %lu]",
            (unsigned long)dst_pos, (unsigned long)len, (unsigned long)dst->len);

    return;
  }
  else if(src_pos + len > src->len)
  {
    fprintf(stderr, "StrBuf OutOfBounds Error: "
                    "strbuf_copy [src; pos: %lu; len: %lu; strlen: %lu]",
            (unsigned long)src_pos, (unsigned long)len, (unsigned long)src->len);

    return;
  }

  strbuf_overwrite_str(dst, dst_pos, src->buff+src_pos, len);
}

void strbuf_overwrite_str(StrBuf* dst, t_buf_pos ddst_pos,
                          const char* src, t_buf_pos len)
{
  t_buf_pos dst_pos = ddst_pos;

  if(src == NULL || len == 0)
  {
    return;
  }
  else if(dst_pos > dst->len)
  {
    // Insert position cannot be greater than current string length
    fprintf(stderr, "StrBuf OutOfBounds Error: "
                    "strbuf_str_copy(index: %lu) [strlen: %lu]",
            (unsigned long)dst_pos, (unsigned long)dst->len);

    return;
  }

  // Check if dest buffer can handle string
  strbuf_ensure_capacity(dst, dst_pos + len);

  // memmove instead of strncpy, as it can handle overlapping regions
  memmove(dst->buff+dst_pos, src, (size_t)len);

  if(dst_pos + len > dst->len)
  {
    // Extended string - add '\0' char
    dst->len = dst_pos + len;
    dst->buff[dst->len] = '\0';
  }
}

void strbuf_insert_strn(StrBuf* dst, t_buf_pos dst_pos,
                        const char* src, t_buf_pos len)
{
  if(src == NULL || len == 0)
  {
    return;
  }
  else if(dst_pos > dst->len)
  {
    // Insert position cannot be greater than current string length
    fprintf(stderr, "StrBuf OutOfBounds Error: "
                    "strbuf_insert_strn(index: %lu) [strlen: %lu]",
            (unsigned long)dst_pos, (unsigned long)dst->len);

    return;
  }

  // Check if dest buffer can handle string plus \0
  strbuf_ensure_capacity(dst, dst->len + len);

  // dst_pos could be at the end (== dst->len)
  if(dst_pos < dst->len)
  {
    // Shift some characters up
    memmove(dst->buff + dst_pos + len,
            dst->buff + dst_pos,
            (size_t)(dst->len - dst_pos));
  }

  // Insert
  memmove(dst->buff + dst_pos, src, (size_t)len);

  // Update size
  dst->len = dst->len + len;
  dst->buff[dst->len] = '\0';
}

void strbuf_insert_str(StrBuf* dst, t_buf_pos dst_pos, const char* src)
{
  strbuf_insert_strn(dst, dst_pos, src, strlen(src));
}

void strbuf_insert(StrBuf* dst, t_buf_pos dst_pos,
                   const StrBuf* src, t_buf_pos src_pos,
                   t_buf_pos len)
{
  strbuf_insert_strn(dst, dst_pos, src->buff+src_pos, len);
}

void strbuf_insert_char(StrBuf* dst, t_buf_pos dst_pos, char c)
{
  strbuf_insert_strn(dst, dst_pos, &c, 1);
}

/**********************/
/* Printing functions */
/**********************/

// Print to stdout. Returns number of bytes printed
int strbuf_puts(StrBuf* sbuf)
{
  return puts(sbuf->buff);
}

int strbuf_fputs(StrBuf* sbuf, FILE* out)
{
  return fputs(sbuf->buff, out);
}

size_t strbuf_fwrite(StrBuf* sbuf, t_buf_pos pos, t_buf_pos len, FILE* out)
{
  return fwrite(sbuf->buff + pos, sizeof(char), len, out);
}

// gz versions
int strbuf_gzputs(StrBuf* sbuf, gzFile gzout)
{
  return gzputs(gzout, sbuf->buff);
}

int strbuf_gzwrite(StrBuf* sbuf, t_buf_pos pos, t_buf_pos len, gzFile gzout)
{
  return gzwrite(gzout, sbuf->buff + pos, (unsigned)len);
}

/**************************/
/*         sprintf        */
/**************************/

int strbuf_vsprintf(StrBuf *sbuf, t_buf_pos pos, const char* fmt, va_list argptr)
{
  // Bounds check
  if(pos > sbuf->len)
  {
    fprintf(stderr, "StrBuf OutOfBounds Error: "
                    "strbuf_vsprintf(index: %lu) [strlen: %lu]",
            (unsigned long)pos, (unsigned long)sbuf->len);

    return -1;
  }

  // Length of remaining buffer
  size_t buf_len = (size_t)(sbuf->size - pos);

  if(buf_len == 0)
  {
    strbuf_resize(sbuf, 2*sbuf->size);
  }

  // Make a copy of the list of args incase we need to resize buff and try again
  va_list argptr_cpy;
  va_copy(argptr_cpy, argptr);

  int num_chars = vsnprintf(sbuf->buff+pos, buf_len, fmt, argptr);

  // num_chars is the number of chars that would be written (not including '\0')
  // num_chars < 0 => failure
  if((t_buf_pos)(num_chars + 1) >= buf_len)
  {
    strbuf_ensure_capacity(sbuf, pos+num_chars);

    // now use the argptr copy we made earlier
    // Don't need to use vsnprintf now, vsprintf will do since we know it'll fit
    num_chars = vsprintf(sbuf->buff+pos, fmt, argptr_cpy);

    va_end(argptr_cpy);
  }

  // Don't need to NUL terminate, vsprintf/vnsprintf does that for us
  if(num_chars < 0)
  {
    // Errors occurred - report, and make sure string is terminated
    fprintf(stderr, "Warning: strbuf_sprintf something went wrong..\n");
    sbuf->buff[sbuf->len] = '\0';
  }
  else
  {
    // Update length
    sbuf->len = pos + num_chars;
  }

  return num_chars;
}

// Appends sprintf
int strbuf_sprintf(StrBuf *sbuf, const char* fmt, ...)
{
  va_list argptr;
  va_start(argptr, fmt);
  int num_chars = strbuf_vsprintf(sbuf, sbuf->len, fmt, argptr);
  va_end(argptr);

  return num_chars;
}

int strbuf_sprintf_at(StrBuf *sbuf, t_buf_pos pos, const char* fmt, ...)
{
  // Bounds check
  if(pos > sbuf->len)
  {
    fprintf(stderr, "StrBuf OutOfBounds Error: "
                    "strbuf_sprintf_at(index: %lu) [strlen: %lu]",
            (unsigned long)pos, (unsigned long)sbuf->len);

    return -1;
  }

  va_list argptr;
  va_start(argptr, fmt);
  int num_chars = strbuf_vsprintf(sbuf, pos, fmt, argptr);
  va_end(argptr);

  return num_chars;
}

// Does not prematurely end the string if you sprintf within the string
// (vs at the end)
int strbuf_sprintf_noterm(StrBuf *sbuf, t_buf_pos pos, const char* fmt, ...)
{
  // Bounds check
  if(pos > sbuf->len)
  {
    fprintf(stderr, "StrBuf OutOfBounds Error: "
                    "strbuf_sprintf_noterm(index: %lu) [strlen: %lu]",
            (unsigned long)pos, (unsigned long)sbuf->len);

    return -1;
  }

  va_list argptr;
  va_start(argptr, fmt);

  int num_chars = vsnprintf(NULL, 0, fmt, argptr);
  
  va_end(argptr);
  
  // Save overwritten char
  char last_char;
  
  if(pos + num_chars < sbuf->len)
  {
    last_char = sbuf->buff[pos + num_chars];
  }
  else
  {
    last_char = '\0';
  }

  va_start(argptr, fmt);

  num_chars = strbuf_vsprintf(sbuf, pos, fmt, argptr);

  va_end(argptr);
  
  // Re-instate overwritten character
  sbuf->buff[pos+num_chars] = last_char;

  return num_chars;
}


/*****************/
/* File handling */
/*****************/

#define _func_read(name,type,func) \
t_buf_pos name(StrBuf *sbuf, type file, t_buf_pos len)                         \
{                                                                              \
  if(len == 0)                                                                 \
  {                                                                            \
    return 0;                                                                  \
  }                                                                            \
                                                                               \
  strbuf_ensure_capacity(sbuf, sbuf->len + len);                               \
                                                                               \
  if(func)                                                                     \
  {                                                                            \
    return 0;                                                                  \
  }                                                                            \
                                                                               \
  t_buf_pos num_of_chars_read = (t_buf_pos)strlen(sbuf->buff + sbuf->len);     \
  sbuf->len += num_of_chars_read;                                              \
  return num_of_chars_read;                                                    \
}

_func_read(strbuf_gzread, gzFile,
           gzgets(file, (char*)(sbuf->buff + sbuf->len), (unsigned)(len+1)) == Z_NULL)

_func_read(strbuf_read, FILE*,
           fgets((char*)(sbuf->buff + sbuf->len), (unsigned)(len+1), file) == NULL)

#define _func_readline(name,type,func) \
t_buf_pos name(StrBuf *sbuf, type file)                                        \
{                                                                              \
  t_buf_pos init_str_len = sbuf->len;                                          \
  strbuf_ensure_capacity(sbuf, sbuf->len+1);                                   \
  while(func)                                                                  \
  {                                                                            \
    t_buf_pos num_of_chars_read = (t_buf_pos)strlen(sbuf->buff + sbuf->len);   \
    char* last_char = (char*)(sbuf->buff + sbuf->len + num_of_chars_read - 1); \
    sbuf->len += num_of_chars_read;                                            \
    if(*last_char == '\n' || *last_char == '\r')                               \
    {                                                                          \
      return sbuf->len - init_str_len;                                         \
    }                                                                          \
    else                                                                       \
    {                                                                          \
      strbuf_resize_vital(sbuf, 2*sbuf->size);                                 \
    }                                                                          \
  }                                                                            \
  return (sbuf->len - init_str_len);                                           \
}

// read gzFile
// returns number of characters read
// or 0 if EOF
_func_readline(strbuf_gzreadline, gzFile,
               gzgets(file, (char*)(sbuf->buff + sbuf->len),
                      (int)(sbuf->size - sbuf->len)) != Z_NULL)

// read FILE
// returns number of characters read
// or 0 if EOF
_func_readline(strbuf_readline, FILE*,
               fgets((char*)(sbuf->buff + sbuf->len),
                     (int)(sbuf->size - sbuf->len), file) != NULL)

// read FILE
// returns number of characters read
// or 0 if EOF
t_buf_pos strbuf_reset_readline(StrBuf *sbuf, FILE *file)
{
  strbuf_reset(sbuf);
  return strbuf_readline(sbuf, file);
}


// read gzFile
// returns number of characters read
// or 0 if EOF
t_buf_pos strbuf_reset_gzreadline(StrBuf *sbuf, gzFile gz_file)
{
  strbuf_reset(sbuf);
  return strbuf_gzreadline(sbuf, gz_file);
}


// These two functions are the same apart from calling fgetc / gzgetc
// (strbuf_skip_line, strbuf_gzskip_line)
t_buf_pos strbuf_skip_line(FILE *file)
{
  char c;
  t_buf_pos count = 0;
  
  while((c = fgetc(file)) != -1)
  {
    count++;
    
    if(c == '\n' || c == '\r')
    {
      break;
    }
  }

  return count;
}

t_buf_pos strbuf_gzskip_line(gzFile gz_file)
{
  char c;
  t_buf_pos count = 0;
  
  while((c = gzgetc(gz_file)) != -1)
  {
    count++;
    
    if(c == '\n' || c == '\r')
    {
      break;
    }
  }

  return count;
}

//
// String functions
//

// Trim whitespace characters from the start and end of a string
void strbuf_trim(StrBuf *sbuf)
{
  if(sbuf->len == 0)
    return;

  // Trim end first
  while(sbuf->len > 0 && isspace(sbuf->buff[sbuf->len-1]))
    sbuf->len--;

  sbuf->buff[sbuf->len] = '\0';

  if(sbuf->len == 0)
    return;

  t_buf_pos start = 0;

  while(start < sbuf->len && isspace(sbuf->buff[start]))
    start++;

  if(start != 0)
  {
    sbuf->len -= start;
    memmove(sbuf->buff, sbuf->buff+start, sbuf->len);
    sbuf->buff[sbuf->len] = '\0';
  }
}

// Trim the characters listed in `list` from the left of `sbuf`
// `list` is a null-terminated string of characters
void strbuf_ltrim(StrBuf *sbuf, char* list)
{
  t_buf_pos start = 0;

  while(start < sbuf->len && strchr(list, sbuf->buff[start]) != NULL)
    start++;

  if(start != 0)
  {
    sbuf->len -= start;
    memmove(sbuf->buff, sbuf->buff+start, sbuf->len);
    sbuf->buff[sbuf->len] = '\0';
  }
}

// Trim the characters listed in `list` from the right of `sbuf`
// `list` is a null-terminated string of characters
void strbuf_rtrim(StrBuf *sbuf, char* list)
{
  if(sbuf->len == 0)
    return;

  while(sbuf->len > 0 && strchr(list, sbuf->buff[sbuf->len-1]) != NULL)
    sbuf->len--;

  sbuf->buff[sbuf->len] = '\0';
}

/**************************/
/* Other String Functions */
/**************************/

char string_is_all_whitespace(const char* s)
{
  int i;

  for(i = 0; s[i] != '\0'; i++)
  {
    if(!isspace(s[i]))
    {
      return 0;
    }
  }

  return 1;
}

char* string_next_nonwhitespace(const char* s)
{
  while(*s != '\0')
  {
    if(!isspace(*s))
    {
      return (char*)s;
    }

    s++;
  }

  return NULL;
}

// Strip whitespace the the start and end of a string.  
// Strips whitepace from the end of the string with \0, and returns pointer to
// first non-whitespace character
char* string_trim(char* str)
{
  // Work backwards
  size_t len = strlen(str);

  while(len > 0 && isspace(*(str+len-1)))
  {
    len--;
  }

  *(str+len) = '\0';

  // Work forwards
  while(isspace(*str)) // don't need start < len because will hit \0
  {
    str++;
  }

  return str;
}

// Removes \r and \n from the ends of a string and returns the new length
size_t string_chomp(char* str)
{
  size_t len = strlen(str);

  while(len > 0 && (str[len-1] == '\r' || str[len-1] == '\n'))
  {
    len--;
  }

  str[len] = '\0';

  return len;
}

// Returns count
size_t string_count_char(const char* str, int c)
{
  size_t count = 0;
  const char *tmp = str;

  while((tmp = strchr(tmp, c)) != NULL)
  {
    tmp++;
    count++;
  }

  return count;
}

// Returns the number of strings resulting from the split
long string_split(const char* split, const char* txt, char*** result)
{
  size_t split_len = strlen(split);
  size_t txt_len = strlen(txt);

  // result is temporarily held here
  char** arr;

  if(split_len == 0)
  {
    // Special case
    if(txt_len == 0)
    {
      *result = NULL;
      return 0;
    }
    else
    {
      arr = (char**) malloc(txt_len * sizeof(char*));
    
      t_buf_pos i;

      for(i = 0; i < txt_len; i++)
      {
        arr[i] = (char*) malloc(2 * sizeof(char));
        arr[i][0] = txt[i];
        arr[i][1] = '\0';
      }

      *result = arr;
      return txt_len;
    }
  }
  
  const char* find = txt;
  long count = 1; // must have at least one item

  while((find = strstr(find, split)) != NULL)
  {
    //printf("Found1: '%s'\n", find);
    count++;
    find += split_len;
  }

  // Create return array
  arr = (char**) malloc(count * sizeof(char*));
  
  count = 0;
  const char* last_position = txt;

  while((find = strstr(last_position, split)) != NULL)
  {
    long str_len = find - last_position;

    arr[count] = (char*) malloc((str_len+1) * sizeof(char));
    strncpy(arr[count], last_position, str_len);
    arr[count][str_len] = '\0';
    
    count++;
    last_position = find + split_len;
  }

  // Copy last item
  long str_len = txt + txt_len - last_position;
  arr[count] = (char*) malloc((str_len+1) * sizeof(char));

  if(count == 0)
  {
    strcpy(arr[count], txt);
  }
  else
  {
    strncpy(arr[count], last_position, str_len);
  }

  arr[count][str_len] = '\0';
  count++;
  
  *result = arr;
  
  return count;
}
