/*
 * Copyright 2009-2013 Zamin Iqbal and Mario Caccamo
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
/*
  global.h
*/

#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <err.h>
#include <stdint.h>
#include "string_buffer.h"

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
