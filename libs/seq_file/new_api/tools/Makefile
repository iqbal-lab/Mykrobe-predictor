ifndef CC
	CC=gcc
endif

CFLAGS=-O3 -Wall -Wextra -I .. -I $(HTSLIB)/htslib/ -L $(HTSLIB)/htslib/
LINKING=-lhts -lpthread -lz

all: ../bin/facat ../bin/fqcat ../bin/seqcat

../bin/facat: htslib ../bin seq_cat.c ../seq_file.h ../buffered_input.h
	$(CC) -DFASTA=1 $(CFLAGS) -o ../bin/facat seq_cat.c $(LINKING)
../bin/fqcat: htslib ../bin seq_cat.c ../seq_file.h ../buffered_input.h
	$(CC) -DFASTQ=1 $(CFLAGS) -o ../bin/fqcat seq_cat.c $(LINKING)
../bin/seqcat: htslib ../bin seq_cat.c ../seq_file.h ../buffered_input.h
	$(CC) -DPLAIN=1 $(CFLAGS) -o ../bin/seqcat seq_cat.c $(LINKING)

htslib:
	if [[ '$(HTSLIB)' == '' ]]; \
	then echo "Error: Please pass HTSLIB=... with path to htslib"; exit 1; fi

../bin:
	mkdir -p ../bin

clean:
	rm -rf ../bin/*

.PHONY: all clean htslib
