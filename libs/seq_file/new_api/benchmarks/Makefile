ifndef CC
	CC=gcc
endif

MINGW=1

CFLAGS= -m64 -D__mingw__ -Wall -Wc++-compat -O3 -Wall -Wextra -I .. -I $(HTSLIB)/htslib/ -L $(HTSLIB)/htslib/


LINKING=-lhts -lpthread -lz

ifdef MINGW
	LINKING := $(LINKING) -lws2_32
endif

all: seqtest ktest

seqtest: htslib seq_test.c ../seq_file.h ../buffered_input.h
	$(CC) -O3 $(CFLAGS) -o seqtest seq_test.c $(LINKING)
ktest: htslib ktest.c kseq.h
	$(CC) -O3 -o ktest ktest.c -lz

htslib:
	if [[ '$(HTSLIB)' == '' ]]; \
	then echo "Error: Please pass HTSLIB=... with path to htslib"; exit 1; fi

clean:
	rm -rf seqtest ktest

.PHONY: all clean htslib
