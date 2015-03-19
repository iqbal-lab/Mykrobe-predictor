ifndef CC
	CC=gcc
endif

MINGW=1

CFLAGS= -m64 -D__mingw__ -Wall -Wc++-compat -O3 -Wall -Wextra -I .. -I $(HTSLIB)/htslib/ -L $(HTSLIB)/htslib/


LINKING=-lhts -lpthread -lz

ifdef MINGW
	LINKING := $(LINKING) -lws2_32
endif

all: btest seqcmdline seqstat

btest: htslib buffered_input.c ../buffered_input.h
	$(CC) -O3 -I .. -o btest buffered_input.c -lz
seqcmdline: htslib seq_cmdline.c ../seq_file.h ../buffered_input.h
	$(CC) -O3 $(CFLAGS) -o seqcmdline seq_cmdline.c $(LINKING)
seqstat: htslib seq_stat.c ../seq_file.h ../buffered_input.h
	$(CC) $(CFLAGS) -o seqstat seq_stat.c $(LINKING)

htslib:
	if [[ '$(HTSLIB)' == '' ]]; \
	then echo "Error: Please pass HTSLIB=... with path to htslib"; exit 1; fi

clean:
	rm -rf btest seqcmdline seqstat

.PHONY: all clean htslib

