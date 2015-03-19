ifndef CC
  CC = gcc
endif

CFLAGS := -Wall -Wextra
LIBFLAGS := -L. -lstrbuf -lz

ifdef DEBUG
	CFLAGS := $(CFLAGS) -DDEBUG=1 --debug -g -ggdb
else
	CFLAGS := $(CFLAGS) -O3
endif

all: clean
	$(CC) $(CFLAGS) -c string_buffer.c -o string_buffer.o
	ar -csru libstrbuf.a string_buffer.o
	$(CC) $(CFLAGS) strbuf_test.c -o strbuf_test $(LIBFLAGS)

clean:
	rm -rf string_buffer.o libstrbuf.a strbuf_test string_buffer.dSYM string_buffer.greg
