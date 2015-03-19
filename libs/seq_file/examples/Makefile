ifndef CC
  CC = gcc
endif

ifndef STRING_BUF_PATH
	STRING_BUF_PATH=$(HOME)/c/libs/string_buffer
endif

ifndef HTS_PATH
	HTS_PATH=$(HOME)/bioinf/htslib/
endif

HTS_PATH:=$(HTS_PATH)/htslib/

INC_ARGS=-I../ -L../ \
         -I $(STRING_BUF_PATH) -L $(STRING_BUF_PATH) \
         -I $(HTS_PATH) -L $(HTS_PATH)

LIB_ARGS=-lseqfile -lstrbuf -lhts -lz

all:
	$(CC) $(LIB_ARGS) -o example example.c $(INC_ARGS)
