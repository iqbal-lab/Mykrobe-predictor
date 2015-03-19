/*
 string_buffer.h
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
 DIRECT, INDIRECT, INCIDEStrBufNTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef STRING_BUFFER_FILE_SEEN
#define STRING_BUFFER_FILE_SEEN

#include <stdio.h> // needed for FILE
#include <zlib.h> // needed for gzFile
#include <stdarg.h> // needed for va_list

typedef unsigned long t_buf_pos;
typedef struct StrBuf StrBuf;

struct StrBuf
{
  char *buff;
  t_buf_pos len; // length of the string
  t_buf_pos size; // buffer size - includes '\0' (size >= len+1)
};

//
// Creation, reset, free and memory expansion
//

// Constructors
StrBuf* strbuf_new();
StrBuf* strbuf_init(t_buf_pos size);
StrBuf* strbuf_create(const char* str);

// Destructors
void strbuf_free(StrBuf* sbuf);
// Free sbuf struct, but retain and return the char array
char* strbuf_free_get_str(StrBuf* sbuf);

// Clone a buffer (including content)
StrBuf* strbuf_clone(const StrBuf* sbuf);

// Get a copy of this StrBuf as a char array
// Returns NULL if not enough memory
char* strbuf_as_str(const StrBuf* sbuf);

// Clear the content of an existing StrBuf (sets size to 0)
void strbuf_reset(StrBuf* sbuf);

// Get number of characters in buffer
t_buf_pos strbuf_len(const StrBuf* sbuf);

// Get current capacity
t_buf_pos strbuf_size(const StrBuf* sbuf);

// If you alter the buffer, call strbuf_update_len to correct the struct
void strbuf_update_len(StrBuf * sbuf);

//
// Resizing
//

// Ensure capacity for len characters plus '\0' character - exits on FAILURE
void strbuf_ensure_capacity(StrBuf *sbuf, t_buf_pos len);

/* More focused -- less used */

// reallocs to exact memory specified - return 1 on success 0 on failure
char strbuf_resize(StrBuf *sbuf, t_buf_pos new_size);

// convenience function: prints error and exits with EXIT_FAILURE if it fails
void strbuf_resize_vital(StrBuf *sbuf, t_buf_pos new_size);

// Shorten string without reallocating memory
void strbuf_shrink(StrBuf *sbuf, t_buf_pos new_len);

//
// Useful String functions
//

// get/set chars
char strbuf_get_char(const StrBuf *sbuf, t_buf_pos index);
void strbuf_set_char(StrBuf *sbuf, t_buf_pos index, char c);

// Set string buffer to contain a given string
void strbuf_set(StrBuf *sbuf, const char *str);

// Add a character to the end of this StrBuf
void strbuf_append_char(StrBuf* sbuf, char c);
// Copy a StrBuf to the end of this StrBuf
void strbuf_append_buff(StrBuf* dst, StrBuf* src);
// Copy a character array to the end of this StrBuf
void strbuf_append_str(StrBuf* sbuf, const char* txt);
// Copy N characters from a character array to the end of this StrBuf
void strbuf_append_strn(StrBuf* sbuf, const char* txt, t_buf_pos len);

// Remove \r and \n characters from the end of this StrBuf
// Returns the number of characters removed
t_buf_pos strbuf_chomp(StrBuf *sbuf);

// Reverse a string
void strbuf_reverse(StrBuf *sbuf);

// Reverse a string region
void strbuf_reverse_region(StrBuf *sbuf, t_buf_pos start, t_buf_pos length);

// Get a substring as a new null terminated char array
// (remember to free the returned char* after you're done with it!)
char* strbuf_substr(StrBuf *sbuf, t_buf_pos start, t_buf_pos len);

//added by Zam - as above, but passing in a preallocated string.
//caller's job to ensure the atring you pass in is long enough
void strbuf_substr_prealloced(StrBuf *sbuf, t_buf_pos start, t_buf_pos len, char* outstr);

// Change to upper or lower case
void strbuf_to_uppercase(StrBuf *sbuf);
void strbuf_to_lowercase(StrBuf *sbuf);

// Copy a string to this StrBuf, overwriting any existing characters
void strbuf_copy(StrBuf* dst, t_buf_pos dst_pos,
                 const StrBuf* src, t_buf_pos src_pos,
                 t_buf_pos len);

// Overwrite a portion of an StrBuf with a new string
// Note: dst_pos + len can be longer the the current dst StrBuf
void strbuf_overwrite_str(StrBuf* dst, t_buf_pos dst_pos,
                          const char* src, t_buf_pos len);

// Insert: copy to a StrBuf, shifting any existing characters along
void strbuf_insert(StrBuf* dst, t_buf_pos dst_pos,
                   const StrBuf* src, t_buf_pos src_pos,
                   t_buf_pos len);

// Insert a string
void strbuf_insert_strn(StrBuf* dst, t_buf_pos dst_pos,
                        const char* src, t_buf_pos len);

void strbuf_insert_str(StrBuf* dst, t_buf_pos ddst_pos, const char* src);

// Insert a single char
void strbuf_insert_char(StrBuf* dst, t_buf_pos dst_pos, char c);

//
// Print to stream
//

// Print to stdout. Returns number of bytes printed
int strbuf_puts(StrBuf* sbuf);

// Print to FILE stream. Returns number of bytes printed
int strbuf_fputs(StrBuf* sbuf, FILE* out);
int strbuf_gzputs(StrBuf* sbuf, gzFile gzout);

size_t strbuf_fwrite(StrBuf* sbuf, t_buf_pos pos, t_buf_pos len, FILE* out);
int strbuf_gzwrite(StrBuf* sbuf, t_buf_pos pos, t_buf_pos len, gzFile gzout);

//
// sprintf
//

// sprintf to a StrBuf (adds string terminator after sprint)
int strbuf_sprintf(StrBuf *sbuf, const char* fmt, ...)
  __attribute__ ((format(printf, 2, 3)));

int strbuf_sprintf_at(StrBuf *sbuf, t_buf_pos pos, const char* fmt, ...)
  __attribute__ ((format(printf, 3, 4)));

int strbuf_vsprintf(StrBuf *sbuf, t_buf_pos pos, const char* fmt, va_list argptr)
  __attribute__ ((format(printf, 3, 0)));

// sprintf without terminating character
// Does not prematurely end the string if you sprintf within the string
// (terminates string if sprintf to the end)
int strbuf_sprintf_noterm(StrBuf *sbuf, t_buf_pos pos, const char* fmt, ...)
  __attribute__ ((format(printf, 3, 4)));

//
// Reading files
//

// Reading a FILE
t_buf_pos strbuf_reset_readline(StrBuf *sbuf, FILE *file);
t_buf_pos strbuf_readline(StrBuf *sbuf, FILE *file);

// Reading a gzFile
t_buf_pos strbuf_reset_gzreadline(StrBuf *sbuf, gzFile gz_file);
t_buf_pos strbuf_gzreadline(StrBuf *sbuf, gzFile gz_file);

// Skip a line and return how many characters were skipped
t_buf_pos strbuf_skip_line(FILE *file);
t_buf_pos strbuf_gzskip_line(gzFile gz_file);

// Read a line but no more than len bytes
t_buf_pos strbuf_read(StrBuf *sbuf, FILE *file, t_buf_pos len);
t_buf_pos strbuf_gzread(StrBuf *sbuf, gzFile gz_file, t_buf_pos len);

//
// String functions
//

// Trim whitespace characters from the start and end of a string
void strbuf_trim(StrBuf *sbuf);

// Trim the characters listed in `list` from the left of `sbuf`
// `list` is a null-terminated string of characters
void strbuf_ltrim(StrBuf *sbuf, char* list);

// Trim the characters listed in `list` from the right of `sbuf`
// `list` is a null-terminated string of characters
void strbuf_rtrim(StrBuf *sbuf, char* list);

/**************************/
/* Other String functions */
/**************************/

char string_is_all_whitespace(const char* s);
char* string_next_nonwhitespace(const char* s);
char* string_trim(char* str);
size_t string_chomp(char* str);
size_t string_count_char(const char* str, int c);
long string_split(const char* split, const char* txt, char*** result);

#endif
