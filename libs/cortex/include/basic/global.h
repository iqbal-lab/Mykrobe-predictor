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
/*
  global.h
*/

#ifndef GLOBAL_H_
#define GLOBAL_H_

#ifndef _WIN32
  #include <err.h>
#endif

#include <stdint.h>
#include "string_buffer.h"

#ifdef __mingw__
  #define __fopen_read_override_mingw(a,b) fopen(a,b)
#else
  #define __fopen_read_override_mingw(a,b) fopen(a,"r")
#endif

typedef uint32_t Covg;
// COVG_MAX is defined as UINT_MAX in global.c
extern const Covg COVG_MAX;

Covg mean_of_covgs(Covg a, Covg b);
Covg sum_covgs(Covg a, Covg b);

typedef signed char boolean;
#ifndef true
#define true 1
#define false 0
#endif


typedef enum
  {
    _True=0,
    _False=1,
    _Inconclusive=2,
  } Troolean;

typedef enum
{
  forward = 0,
  reverse = 1
} Orientation;


#define MAX_READ_NAME_LEN 300
#define VERSION 0
#define SUBVERSION 0
#define SUBSUBVERSION 0
#define SUBSUBSUBVERSION 1
boolean DEBUG;

#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define MIN(x,y) ((x) <= (y) ? (x) : (y))

boolean test_file_existence(char* file);

int int_cmp(const void *a, const void *b);
int Covg_cmp(const void *a, const void *b);
int float_cmp(const void* a, const void* b);
int double_cmp(const void* a, const void* b);
int long_double_cmp(const void* a, const void* b);

void set_string_to_null(char* str, int len);

void set_int_array(int*arr, int len, int val);

void set_uint64_t_array(uint64_t* arr, int len, uint64_t val);

void die(const char* fmt, ...)
  __attribute__ ((format(printf, 1, 2)))
  __attribute__ ((noreturn));

void warn(const char* fmt, ...)
  __attribute__ ((format(printf, 1, 2)));

// Placeholder for message() -- a function for standard output
void message(const char* fmt, ...)
  __attribute__ ((format(printf, 1, 2)));

void strbuf_remove_all_whitespace(StrBuf* sbuf);
void strbuf_search_replace(StrBuf* sbuf, char find, char repl);
int strbuf_find_first(StrBuf* sbuf, char find);
void strbuf_add_slash_on_end(StrBuf* sbuf);
void strbuf_rev_comp(StrBuf* sb);
// backported by Zam
void strbuf_delete(StrBuf *sbuf, size_t pos, size_t len);

#endif /* GLOBAL_H_ */
