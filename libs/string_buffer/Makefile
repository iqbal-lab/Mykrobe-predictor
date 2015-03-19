UNAME=$(shell uname)
ifeq ($(UNAME),Darwin)
	MAC=1
else
	ifeq ($(UNAME),CYGWIN_NT-6.1)
	WIN=1
	endif
endif

ifndef CC
  CC = gcc
endif

ifdef WIN
	CFLAGS := -m64 -D__mingw__ -Wall -Wc++-compat  -Wall -Wextra -I/usr/include
	LIBFLAGS := -L. -L/usr/lib -lstrbuf -lz
	CC = gcc
else
	CFLAGS := -Wall -Wextra
	LIBFLAGS := -L. -lstrbuf -lz
endif


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
