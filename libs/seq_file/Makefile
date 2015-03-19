ifndef CC
  CC = gcc
endif

ifdef DEBUG
	CFLAGS := -DDEBUG=1 --debug -g
else
	CFLAGS := -O3
endif

CFLAGS := $(CFLAGS) -Wall -Wextra -I $(HTS_PATH)/htslib/ -I $(STRING_BUF_PATH)

LIB_FLAGS = -L $(HTS_PATH)/htslib/ -L $(STRING_BUF_PATH) -lstrbuf -lhts -lpthread -lz -lm

OBJS = seq_file.o seq_common.o seq_fasta.o seq_fastq.o seq_plain.o seq_sam.o

all: htslib string_buffer clean $(OBJS)
	ar -csru libseqfile.a $(OBJS)
	$(CC) -o seq_convert $(CFLAGS) seq_convert.c libseqfile.a $(LIB_FLAGS)
	$(CC) -o seq_file_test $(CFLAGS) seq_file_test.c libseqfile.a $(LIB_FLAGS)
	cd new_api; make HTSLIB=$(shell readlink -f $(HTS_PATH))

htslib:
	if [[ '$(HTS_PATH)' == '' ]]; \
	then echo "Error: Please pass HTS_PATH=... with path to htslib dir"; exit 1; fi

string_buffer:
	if [[ '$(STRING_BUF_PATH)' == '' ]]; \
	then echo "Error: Please pass STRING_BUF_PATH=... with path to string_buffer dir"; exit 1; fi

clean:
	rm -rf $(OBJS) libseqfile.a seq_convert seq_file_test \
	       seq_file.dSYM seq_file.greg seq_convert.dSYM seq_file_test.dSYM
	cd new_api; make clean

.PHONY: all clean htslib string_buffer

%.o : %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@
