C String Buffer
===============
Library code for handling strings and reading from files  
project: string_buffer  
url: https://github.com/noporpoise/StringBuffer  
author: Isaac Turner <turner.isaac@gmail.com>  

About
=====

A string buffer library for C. Only has zlib as a dependency. Compiles with gcc
and clang.

Features:
- copying, inserting, appending, substring, chomp, trim
- reverse region, convert to upper/lower case
- sprintf into string buffer
- read a line at a time from a file
- write to a file
- gzip file support

To build the test code:

    $ make
    $ ./strbuf_test

Calling
=======

To use in your code, include the following arguments in your gcc command:

    gcc ... -I$(STRING_BUF_PATH) -L$(STRING_BUF_PATH) ... -lstrbuf -lz

and include in your source code:

    include "string_buffer.h"

Example Code
============

    StrBuf* myBuff = strbuf_new()

    // Read from a file:

    gzFile fgz = gzopen("path/here.txt.gz")

    while(strbuf_gzreadline(myBuff, fgz))
    {
      // Do something with the line

      // e.g. chomp (remove newline)
      strbuf_chomp(myBuff)

      // e.g. print to stdout
      strbuf_puts(myBuff)

      // Print out a newline
      putc("\n")

      // Reset StrBuf so you're not just concatenating all the lines in memory
      strbuf_reset(myBuff)
    }

    // Close up
    gzclose(fgz)

    strbuf_free(myBuff)


String buffers can still be used as input to standard str functions by accessing
the char* in the StrBuf struct. e.g.:

Get the position of the first 'a' in a StrBuf

    char* ptr = strchr(strbuf->buff, 'a')
    int pos = (ptr == NULL ? -1 : ptr - strbuf->buff)

Test if the StrBuf contains 'hello' from index 12

    if(strncasecmp(strbuf->buff+12, "hello", 5) == 0)
      puts("world!\n")


Functions
=========

    struct StrBuf
    {
      char *buff;
      t_buf_pos len; // length of the string
      t_buf_pos size; // buffer size - includes '\0' (size is always >= len+1)
    };

Creators, destructors etc.
--------------------------

Constructors

    StrBuf* strbuf_new()
    StrBuf* strbuf_init(const t_buf_pos size)
    StrBuf* strbuf_create(const char* str)

Destructors

    void strbuf_free(StrBuf* sbuf)

Free sbuf struct, but retain and return the char array

    char* strbuf_free_get_str(StrBuf* sbuf)

Clone a string buffer (including content)

    StrBuf* strbuf_clone(const StrBuf* sbuf)

Get a copy of this StrBuf as a char array.
Returns NULL if not enough memory.

    char* strbuf_as_str(const StrBuf* sbuf)

Clear the content of an existing StrBuf (sets size to 0)

    void strbuf_reset(StrBuf* sbuf)

Get number of characters in buffer

    t_buf_pos strbuf_len(const StrBuf* sbuf)

Get current capacity

    t_buf_pos strbuf_size(const StrBuf* sbuf)

Resizing
--------

Ensure capacity for len characters plus '\0' character - exits on FAILURE

    void strbuf_ensure_capacity(StrBuf *sbuf, const t_buf_pos len)

*More specialised -- used less frequently:*

reallocs to exact memory specified - return 1 on success 0 on failure

    char strbuf_resize(StrBuf *sbuf, const t_buf_pos new_size)

convenience function: prints error and exits with EXIT_FAILURE if it fails

    void strbuf_resize_vital(StrBuf *sbuf, const t_buf_pos new_size)

Shorten string without reallocating memory

    void strbuf_shrink(StrBuf *sbuf, const t_buf_pos new_len)

Useful functions
----------------

get/set chars

    char strbuf_get_char(const StrBuf *sbuf, const t_buf_pos index)
    void strbuf_set_char(StrBuf *sbuf, const t_buf_pos index, const char c)

Set string buffer to contain a given string
The string can be a string within the given string buffer

    void strbuf_set(StrBuf *sbuf, const char *str)

Add a character to the end of this StrBuf

    void strbuf_append_char(StrBuf* sbuf, const char txt)
Copy a StrBuf to the end of this StrBuf

    void strbuf_append_buff(StrBuf* dst, StrBuf* src)
Copy a character array to the end of this StrBuf

    void strbuf_append_str(StrBuf* sbuf, const char* txt)
Copy N characters from a character array to the end of this StrBuf

    void strbuf_append_strn(StrBuf* sbuf, const char* txt, const t_buf_pos len)

Remove \r and \n characters from the end of this StrBuf.
Returns the number of characters removed

    t_buf_pos strbuf_chomp(StrBuf *sbuf)

Reverse a string

    void strbuf_reverse(StrBuf *sbuf)

Reverse a string region

    void strbuf_reverse_region(StrBuf *sbuf, t_buf_pos start, t_buf_pos length)

Get a substring as a new null terminated char array
(remember to free the returned char* after you're done with it!)

    char* strbuf_substr(StrBuf *sbuf, const t_buf_pos start, const t_buf_pos len)

Change to upper or lower case

    void strbuf_to_uppercase(StrBuf *sbuf)
    void strbuf_to_lowercase(StrBuf *sbuf)

Copy a string to this StrBuf, overwriting any existing characters

    void strbuf_copy(StrBuf* dst, const t_buf_pos dst_pos,
                     const StrBuf* src, const t_buf_pos src_pos,
                     const t_buf_pos len)

Overwrite a portion of an StrBuf with a new string
Note: dst_pos + len can be longer the the current dst StrBuf

    void strbuf_overwrite_str(StrBuf* dst, const t_buf_pos dst_pos,
                              const char* src, const t_buf_pos len)

Insert -- copy to a StrBuf, shifting any existing characters along to the right

    void strbuf_insert(StrBuf* dst, const t_buf_pos dst_pos,
                       const StrBuf* src, const t_buf_pos src_pos,
                       const t_buf_pos len)

Insert a from a `char*`

    void strbuf_insert_str(StrBuf* dst, const t_buf_pos dst_pos,
                           const char* src, const t_buf_pos len)

Insert a single char

    void strbuf_insert_char(StrBuf* dst, const t_buf_pos dst_pos, const char c)

Printing to streams
-------------------

Print to stdout. Returns number of bytes printed

    int strbuf_puts(StrBuf* sbuf)

Print to FILE stream. Returns number of bytes printed

    int strbuf_fputs(StrBuf* sbuf, FILE* out)

    size_t strbuf_fwrite(StrBuf* sbuf, const t_buf_pos pos, const t_buf_pos len,
                         FILE* out)

    int strbuf_gzputs(StrBuf* sbuf, gzFile gzout)

    int strbuf_gzwrite(StrBuf* sbuf, const t_buf_pos pos, const t_buf_pos len,
                       gzFile gzout)

Formatted strings (sprintf)
---------------------------

sprintf to a StrBuf (adds string terminator after sprint)

    int strbuf_sprintf(StrBuf *sbuf, const char* fmt, ...)
    int strbuf_sprintf_at(StrBuf *sbuf, const t_buf_pos pos, const char* fmt, ...)
    int strbuf_vsprintf(StrBuf *sbuf, const t_buf_pos pos,
                        const char* fmt, va_list argptr)

sprintf without terminating character.
Does not prematurely end the string if you sprintf within the string
(terminates string if sprintf to the end)

    int strbuf_sprintf_noterm(StrBuf *sbuf, const t_buf_pos pos,
                              const char* fmt, ...)

Reading files
-------------

Reading a FILE

    t_buf_pos strbuf_reset_readline(StrBuf *sbuf, FILE *file)
    t_buf_pos strbuf_readline(StrBuf *sbuf, FILE *gz_file)

Reading a gzFile

    t_buf_pos strbuf_reset_gzreadline(StrBuf *sbuf, gzFile gz_file)
    t_buf_pos strbuf_gzreadline(StrBuf *sbuf, gzFile gz_file)

Skip a line and return how many characters were skipped

    t_buf_pos strbuf_skip_line(FILE *file)
    t_buf_pos strbuf_gzskip_line(gzFile gz_file)

Read a line but no more than len bytes

    t_buf_pos strbuf_read(StrBuf *sbuf, FILE *file, t_buf_pos len)
    t_buf_pos strbuf_gzread(StrBuf *sbuf, gzFile file, t_buf_pos len)

String functions
----------------

Trim whitespace characters from the start and end of a string

    void strbuf_trim(StrBuf *sbuf)

Trim the characters listed in `list` from the left of `sbuf`.
`list` is a null-terminated string of characters

    void strbuf_ltrim(StrBuf *sbuf, char* list)

Trim the characters listed in `list` from the right of `sbuf`.
`list` is a null-terminated string of characters

    void strbuf_rtrim(StrBuf *sbuf, char* list)

Use with sscanf
---------------

To read strings into a string buffer using `sscanf`, first you must ensure the
buffer is big enough, and afterwards you must ensure the length is stored
correctly.

    StrBuf *sbuf = strbuf_new();
    char *input = "I'm sorry Dave I can't do that";
    
    strbuf_ensure_capacity(sbuf, strlen(input));
    sscanf(input, "I'm sorry %s I can't do that", sbuf->buff);
    strbuf_update_len(sbuf);

    printf("Name: '%s'\n", sbuf->buff);

    strbuf_free(sbuf);

Other string functions
----------------------

These work on `char*` not `StrBuf`, but they're here because they're useful. 

    char string_is_all_whitespace(const char* s)
    char* string_next_nonwhitespace(const char* s)
    char* string_trim(char* str)
    size_t string_chomp(char* str)
    size_t string_count_char(const char* str, const int c)
    long string_split(const char* split, const char* txt, char*** result)


License
=======

    Copyright (c) 2011-2, Isaac Turner
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

Development
===========

Short term goals: none -- please suggest some!

I like to hear about how you're using it, what bugs you've found and what
features you'd like to see!  Contact me: Isaac Turner <turner.isaac@gmail>
