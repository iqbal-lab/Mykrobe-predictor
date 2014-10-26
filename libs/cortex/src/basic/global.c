/*
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@tgac.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * The MIT License (MIT)
 * Copyright (c) 2009-2014 <Z. Iqbal and M. Caccamo>
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * **********************************************************************
 */

#include <stdlib.h>
#include <stdarg.h> // needed for va_list
#include <stdio.h>
#include <string.h>
#include "global.h"
#include <string_buffer.h>
#include "binary_kmer.h"
#include <limits.h>

const Covg COVG_MAX = UINT_MAX;
//const Covg COVG_MAX = INT_MAX

//want to avoid overflow issues,
//see discussion here:http://stackoverflow.com/questions/1020188/fast-average-without-division
Covg mean_of_covgs(Covg a, Covg b)
{
  return (a&b) + ((a^b)>>1);
}

Covg sum_covgs(Covg a, Covg b)
{
  if(COVG_MAX - b >= a)
    {
      return a+b;
    }
    else
    {
      return  COVG_MAX;
    }
}

boolean test_file_existence(char* file)
{
  FILE* fp = fopen(file, "r");
  if(fp == NULL)
  {
    return false;
  }
  else
  {
    fclose(fp);
    return true;
  }
}


// integer comparison: returns negative if a < b
//                                    0 if a == b
//                             positive if a > b
int int_cmp(const void *a, const void *b)
{
  // casting pointer types
  const int *ia = (const int *)a;
  const int *ib = (const int *)b;

  //return (*ia  - *ib); Hmm - this could overflow if you subtract a negative from a positive.
  return (*ia > *ib) - (*ia < *ib);//this won't overflow
}


int Covg_cmp(const void *a, const void *b)
{
  // casting pointer types
  const Covg *ca = (const Covg *)a;
  const Covg *cb = (const Covg *)b;

  return (*ca > *cb) - (*ca < *cb);
}

int float_cmp(const void* a, const void* b)
{
  float fa = *(float*) a;
  float fb = *(float*) b;
  return (fa > fb) - (fa < fb);
}

int double_cmp(const void* a, const void* b)
{
  double fa = *(double*) a;
  double fb = *(double*) b;
  return (fa > fb) - (fa < fb);
}

int long_double_cmp(const void* a, const void* b)
{
  long double fa = *(long double*) a;
  long double fb = *(long double*) b;
  return (fa > fb) - (fa < fb);
}

void set_string_to_null(char* str, int len)
{
  memset(str, 0, sizeof(char)*len);
}

void set_int_array(int* arr, int len, int val)
{
  int  i;
  for (i=0; i<len; i++)
    {
      arr[i]=val;
    }
}

void set_uint64_t_array(uint64_t* arr, int len, uint64_t val)
{
  int  i;
  for (i=0; i<len; i++)
    {
      arr[i]=val;
    }
}

void die(const char* fmt, ...)
{
  fflush(stdout);

  // Print error
  fprintf(stderr, "Error: ");

  va_list argptr;
  va_start(argptr, fmt);
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);

  // Check if we need to print a newline
  if(*(fmt+strlen(fmt)-1) != '\n')
  {
    fprintf(stderr, "\n");
  }

  exit(EXIT_FAILURE);
}

void warn(const char* fmt, ...)
{
  fflush(stdout);

  // Print warning
  fprintf(stderr, "Warning: ");

  va_list argptr;
  va_start(argptr, fmt);
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);

  // Check if we need to print a newline
  if(*(fmt+strlen(fmt)-1) != '\n')
  {
    fprintf(stderr, "\n");
  }
}

// Placeholder for message() -- a function for standard output
void message(const char* fmt, ...)
{
  va_list argptr;
  va_start(argptr, fmt);
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);
}


void strbuf_remove_all_whitespace(StrBuf* sbuf)
{
  char* str = strbuf_as_str(sbuf);
  strbuf_reset(sbuf);
  int i;

  for(i = 0; str[i]!='\0'; i++)
  {
    if(!isspace(str[i]))
    {
      strbuf_append_char(sbuf, str[i]);
    }
  }
  free(str);
}

void strbuf_search_replace(StrBuf* sbuf, char find, char repl)
{
  uint32_t i;
  for (i=0; i<sbuf->len; i++)
    {
      if (sbuf->buff[i]==find)
	{
	  sbuf->buff[i]=repl;
	}
    }
}

//returns -1 on failure, else the index of the position where the first match is
int strbuf_find_first(StrBuf* sbuf, char find)
{
  uint32_t i;
  for (i=0; i<sbuf->len; i++)
    {
      if (sbuf->buff[i]==find)
	{
	  return i;
	}
    }
  return -1;
}

//useful if represents a path to a directory
void strbuf_add_slash_on_end(StrBuf* sbuf)
{
  if (sbuf->buff[sbuf->len]!='/')
    {
      strbuf_append_str(sbuf, "/");
    }
}

void strbuf_rev_comp(StrBuf* sb)
{
  strbuf_reverse(sb);
  int i;
  for (i=0; i<sb->len; i++)
    {
      char r = reverse_char_nucleotide(sb->buff[i]);
      strbuf_set_char(sb, i, r);
    }
}

//backported by Zam
// Bounds check when reading a range (start+len < strlen is valid)
static void _bounds_check_read_range(const StrBuf *sbuf, size_t start, size_t len,
                                     const char* file, int line, const char* func)
{
  if(start + len > sbuf->len)
    {
      fprintf(stderr, "%s:%i:%s() - out of bounds error "
	      "[start: %zu; length: %zu; strlen: %zu; buf:%.*s%s]\n",
	      file, line, func, start, len, sbuf->len,
	      (int)MIN(5, sbuf->len), sbuf->buff, sbuf->len > 5 ? "..." : "");
      exit(1);
    }
}

// backported by Zam
void strbuf_delete(StrBuf *sbuf, size_t pos, size_t len)
{
  _bounds_check_read_range(sbuf, pos, len, __FILE__, __LINE__, "strbuf_delete");
  memmove(sbuf->buff+pos, sbuf->buff+pos+len, sbuf->len-pos-len);
  sbuf->len -= len;
  sbuf->buff[sbuf->len] = '\0';
}

