/*
 sb_test.c
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include "string_buffer.h"

void _test_trim(const char* str)
{
  size_t len = strlen(str);
  char* str_cpy = (char*) malloc(len+1);
  strcpy(str_cpy, str);
  char* new_str = string_trim(str_cpy);

  printf("trim('%s'): '%s'\n", str, new_str);
  free(str_cpy);
}

void test_trim()
{
  printf("Test trim:\n");
  _test_trim("  \t asdf asdf \n ");
  _test_trim("a");
  _test_trim("");
}

void test_all_whitespace()
{
  printf("Test string_is_all_whitespace:\n");
  char* str = "  \tasdf";
  printf("string_is_all_whitespace('%s'): %i\n", str, string_is_all_whitespace(str));
  str = "  \t ";
  printf("string_is_all_whitespace('%s'): %i\n", str, string_is_all_whitespace(str));
}

void _test_split(char* split, char* txt)
{
  char** results;
  
  printf("split '%s' by '%s': (", txt, split);
  
  long count = string_split(split, txt, &results);

  if(count > 0)
  {
    printf("'%s'", results[0]);
    free(results[0]);
  
    int i;
    for(i = 1; i < count; i++)
    {
      printf(", '%s'", results[i]);
      free(results[i]);
    }

    free(results);
  }

  printf(")\n");
}

void test_split()
{
  _test_split("/", "a/b");
  _test_split("/", "/");
  _test_split("/", "/b");
  _test_split("/", "a/");
  _test_split("/", "asdf");
  _test_split("/", "");
  _test_split("", "asdf");
  _test_split("", "");
}

void test_add_char()
{
  StrBuf* sbuf = strbuf_init(100);
  
  strbuf_append_char(sbuf, 'a');
  strbuf_append_char(sbuf, 'b');
  printf("'%s' (length: %lu)\n", sbuf->buff, sbuf->len);
  strbuf_free(sbuf);
}

void test_sprintf()
{
  printf("printf:\n");
  StrBuf* sbuf = strbuf_init(100);
  
  strbuf_sprintf(sbuf, "hi ello");
  printf("'%s' (length: %lu)\n", sbuf->buff, sbuf->len);
  
  strbuf_sprintf_noterm(sbuf, 0, "woot %i %s;", 12, "byebye");
  printf("'%s' (length: %lu)\n", sbuf->buff, sbuf->len);

  char *a = "wooo-%s-xx";
  char *b = "hihi";

  strbuf_sprintf_at(sbuf, strbuf_len(sbuf), a, b);
  printf("'%s' (length: %lu)\n", sbuf->buff, sbuf->len);

  strbuf_reset(sbuf);
  strbuf_resize(sbuf, 10);
  strbuf_append_str(sbuf, "asdfasdf");

  strbuf_sprintf(sbuf, a, b);
  printf("'%s' (length: %lu)\n", sbuf->buff, sbuf->len);
  strbuf_free(sbuf);
}

void test_sscanf()
{
  StrBuf *sbuf = strbuf_new();
  char *input = "I'm sorry Dave I can't do that";
  
  strbuf_ensure_capacity(sbuf, strlen(input));
  sscanf(input, "I'm sorry %s I can't do that", sbuf->buff);
  strbuf_update_len(sbuf);

  printf("Name: '%s'\n", sbuf->buff);

  strbuf_free(sbuf);
}

int main(int argc, char* argv[])
{
  if(argc != 1)
  {
    int i;
    printf("Unused arguments: %s", argv[1]);
    for(i = 2; i < argc; i++) {
      printf(", %s", argv[i]);
    }
    printf("\n");
    exit(EXIT_FAILURE);
  }

  test_sscanf();

  // StringBuffer functions
  test_add_char();
  test_sprintf();
  
  // char array functions
  test_all_whitespace();
  test_split();
  test_trim();
  
  return EXIT_SUCCESS;
}
