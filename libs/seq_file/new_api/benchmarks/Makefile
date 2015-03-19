ifndef CC
	CC=gcc
endif

CFLAGS=-Wall -Wextra -I .. -I $(HTSLIB)/htslib/ -L $(HTSLIB)/htslib/
LINKING=-lhts -lpthread -lz

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
