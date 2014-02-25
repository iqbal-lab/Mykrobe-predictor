/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
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

