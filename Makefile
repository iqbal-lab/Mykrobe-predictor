CC = gcc
BIN = bin
BITFIELDS = 1
NUM_COLS = 1

# Test if running on a mac
UNAME=$(shell uname)
ifeq ($(UNAME),Darwin)
	MAC = 1
endif


# Library paths
IDIR_STRS = libs/string_buffer
IDIR_SEQ = libs/seq_file
IDIR_BAM = libs/htslib/htslib

# Main program includes
IDIR_BASIC = libs/cortex/include/basic
IDIR_BASE_ENCODING = ${IDIR_BASIC}/event_encoding/base_encoding
IDIR_HASH = libs/cortex/include/hash_table
IDIR_CORE = libs/cortex/include/core
IDIR_PREDICTOR_STAPH = include/predictor/staph
IDIR_PREDICTOR_TB = include/predictor/tb
IDIR_PREDICTOR_CORE = include/predictor/core


# Test code includes
IDIR_BASIC_TESTS = libs/cortex/include/test/basic
IDIR_HASH_TABLE_TESTS = libs/cortex/include/test/hash_table
IDIR_PREDICTOR_TESTS = include/test/myKrobe/predictor


IDIR_CUNIT = /home/zam/dev/hg/CUnit/CUnit-2.1-0/CUnit/Headers
LDIR_CUNIT = /home/zam/bin/lib
#IDIR_CUNIT = /home/phelimb/local/include/CUnit
#LDIR_CUNIT = /home/phelimb/local/lib


ifdef MAC
	MACFLAG = -fnested-functions
endif

ARCH = -m64

ifdef 32_BITS
	ARCH =
endif


ifndef STAPH
	STAPH=0
endif
ifndef TB
  TB=0
endif


ifeq ($(STAPH),1)
	ifeq ($(TB),1)
	else
		BINNAME="Mykrobe.predictor.staph"
	endif
endif
ifeq ($(STAPH),0)
	ifeq ($(TB),0)
	else
		BINNAME="Mykrobe.predictor.tb"
	endif
endif



# Comment out this line to turn off adding the commit version
# (it already checks if hg is installed)
#VERSION_STR=$(shell if [ `command -v hg` ]; then echo ' (commit' `hg id --num --id`')'; else echo; fi)

OPT := $(ARCH) -Wall $(MACFLAG) -DVERSION_STR='"$(VERSION_STR)"' \
       -DNUMBER_OF_BITFIELDS_IN_BINARY_KMER=$(BITFIELDS) \
       -DNUMBER_OF_COLOURS=$(NUM_COLS)

ifeq ($(STAPH),1)
	OPT := $(OPT) -DSTAPH=$(STAPH)
else
	OPT := $(OPT) -DTB=$(TB)
endif


ifdef DEBUG
	OPT := -O0 -g $(OPT)
else
	OPT := -O3 $(OPT)
endif

LIBLIST = -lseqfile -lstrbuf -lhts -lpthread -lz -lm
TEST_LIBLIST = -lcunit -lncurses $(LIBLIST)

# Add -L/usr/local/lib/ to satisfy some systems that struggle to link libz
LIBINCS = -L/usr/local/lib  -I$(IDIR_BAM) \
          -I$(IDIR_SEQ) -I$(IDIR_STRS) \
          -L$(IDIR_BAM) -L$(IDIR_SEQ) -L$(IDIR_STRS)

TEST_LIBINCS = -I$(IDIR_CUNIT) -L$(LDIR_CUNIT) $(LIBINCS)

CFLAGS_BASIC      = -I$(IDIR_BASIC) -I$(IDIR_BASE_ENCODING) $(LIBINCS)
CFLAGS_PREDICTOR_CORE = -I$(IDIR_CORE) -I$(IDIR_PREDICTOR_CORE) -I$(IDIR_BASIC) -I$(IDIR_HASH)  -I$(IDIR_BASE_ENCODING) $(LIBINCS)
ifeq ($(STAPH),1)
	CFLAGS_PREDICTOR = -I$(IDIR_PREDICTOR_CORE) -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_PREDICTOR_STAPH) -I$(IDIR_BASE_ENCODING) $(LIBINCS)
else
	CFLAGS_PREDICTOR = -I$(IDIR_PREDICTOR_CORE) -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_PREDICTOR_TB) -I$(IDIR_BASE_ENCODING) $(LIBINCS)
endif


CFLAGS_BASIC_TESTS      = -I$(IDIR_BASIC_TESTS) -I$(IDIR_BASIC) -I$(IDIR_BASE_ENCODING) $(TEST_LIBINCS)
CFLAGS_HASH_TABLE_TESTS = -I$(IDIR_HASH) -I$(IDIR_HASH_TABLE_TESTS) -I$(IDIR_CUNIT) -I$(IDIR_BASIC) $(TEST_LIBINCS) -I$(IDIR_BASE_ENCODING) -I$(IDIR_PREDICTOR_CORE) -I$(IDIR_PREDICTOR)
CFLAGS_PREDICTOR_TESTS = -I$(IDIR_PREDICTOR_TESTS) -I$(IDIR_BASIC) -I$(IDIR_BASE_ENCODING) -I$(IDIR_HASH) -I$(IDIR_PREDICTOR) -I$(IDIR_PREDICTOR_CORE) $(TEST_LIBINCS)

PREDICTOR_OBJ = src/obj/predictor/global.o src/obj/predictor/main.o src/obj/predictor/antibiotics.o src/obj/predictor/binary_kmer.o src/obj/predictor/element.o src/obj/predictor/seq.o src/obj/predictor/hash_value.o src/obj/predictor/hash_table.o src/obj/predictor/build.o src/obj/predictor/dB_graph_supernode.o src/obj/predictor/graph_info.o  src/obj/predictor/dB_graph.o src/obj/predictor/db_variants.o src/obj/predictor/cmd_line.o src/obj/predictor/event_encoding.o src/obj/predictor/db_differentiation.o src/obj/predictor/genotyping_known.o src/obj/predictor/known_mutations.o src/obj/predictor/gene_presence.o src/obj/predictor/species.o src/obj/predictor/maths.o src/obj/predictor/file_reader.o src/obj/predictor/mut_models.o src/obj/predictor/gene_presence_models.o  src/obj/predictor/json.o

BASIC_TESTS_OBJ = src/obj/basic/binary_kmer.o src/obj/basic/global.o src/obj/basic/seq.o src/obj/test/basic/test_binary_kmer.o src/obj/test/basic/test_seq.o src/obj/test/basic/run_basic_tests.o src/obj/basic/event_encoding.o

HASH_TABLE_TESTS_OBJ = src/obj/basic/global.o src/obj/test/hash_table/run_hash_table_tests.o src/obj/predictor/element.o src/obj/predictor/hash_value.o src/obj/predictor/hash_table.o src/obj/test/hash_table/test_hash.o src/obj/basic/binary_kmer.o  src/obj/basic/seq.o src/obj/basic/event_encoding.o

PREDICTOR_TESTS_OBJ = src/obj/test/predictor/run_predictor_tests.o src/obj/test/predictor/test_build.o src/obj/test/predictor/test_genotyping_known.o src/obj/test/predictor/test_gene_presence.o src/obj/test/predictor/test_species_prediction.o src/obj/predictor/global.o src/obj/predictor/binary_kmer.o src/obj/predictor/element.o src/obj/predictor/seq.o src/obj/predictor/hash_value.o src/obj/predictor/hash_table.o src/obj/predictor/build.o src/obj/predictor/dB_graph_supernode.o  src/obj/predictor/dB_graph.o src/obj/predictor/db_variants.o src/obj/predictor/event_encoding.o src/obj/predictor/db_differentiation.o src/obj/predictor/maths.o src/obj/predictor/file_reader.o src/obj/predictor/genotyping_known.o src/obj/predictor/known_mutations.o src/obj/predictor/gene_presence.o src/obj/predictor/species.o src/obj/predictor/antibiotics.o src/obj/predictor/graph_info.o src/obj/predictor/mut_models.o src/obj/predictor/gene_presence_models.o  src/obj/predictor/json.o


predictor : remove_objects $(PREDICTOR_OBJ)
	mkdir -p $(BIN); $(CC) $(CFLAGS_PREDICTOR) $(OPT) $(OPT_COLS) -o $(BIN)/$(BINNAME) $(PREDICTOR_OBJ) $(LIBLIST)

run_basic_tests : remove_objects $(BASIC_TESTS_OBJ)
	mkdir -p $(BIN);  $(CC) $(CFLAGS_BASIC_TESTS) $(OPT) -o $(BIN)/run_basic_tests $(BASIC_TESTS_OBJ) $(TEST_LIBLIST)

run_hash_table_tests : remove_objects $(HASH_TABLE_TESTS_OBJ)
	mkdir -p $(BIN);  $(CC) $(CFLAGS_HASH_TABLE_TESTS) $(OPT) -o $(BIN)/run_hash_table_tests $(HASH_TABLE_TESTS_OBJ) $(TEST_LIBLIST)

run_predictor_tests : remove_objects $(PREDICTOR_TESTS_OBJ)
	mkdir -p $(BIN);  $(CC) $(CFLAGS_PREDICTOR_TESTS) $(OPT) -o $(BIN)/run_predictor_tests $(PREDICTOR_TESTS_OBJ) $(TEST_LIBLIST)


.PHONY : clean
clean :
	rm -rf $(BIN)/*
	rm -rf src/obj

remove_objects:
	rm -rf src/obj/*

#pattern rules

src/obj/predictor/%.o : src/predictor/core/%.c include/predictor/core/%.h
	mkdir -p src/obj/predictor; $(CC) $(CFLAGS_PREDICTOR_CORE) $(CFLAGS_PREDICTOR) $(OPT) -c $< -o $@
src/obj/predictor/%.o : libs/cortex/src/core/%.c libs/cortex/include/core/%.h
	mkdir -p src/obj/predictor; $(CC) $(CFLAGS_PREDICTOR_CORE) $(CFLAGS_PREDICTOR) $(OPT) -c $< -o $@


ifeq ($(STAPH),1)
src/obj/predictor/%.o : src/predictor/staph/%.c include/predictor/staph/%.h
	mkdir -p src/obj/predictor; $(CC) $(CFLAGS_PREDICTOR_CORE) $(CFLAGS_PREDICTOR) $(OPT) -c $< -o $@
else
src/obj/predictor/%.o : src/predictor/tb/%.c include/predictor/tb/%.h
	mkdir -p src/obj/predictor; $(CC) $(CFLAGS_PREDICTOR_CORE) $(CFLAGS_PREDICTOR) $(OPT) -c $< -o $@
endif


src/obj/basic/%.o : libs/cortex/src/basic/%.c libs/cortex/include/basic/%.h
	mkdir -p src/obj/basic/; $(CC) $(CFLAGS_BASIC) $(OPT) -c $< -o $@

src/obj/basic/%.o : libs/cortex/src/basic/event_encoding/base_encoding/%.c  libs/cortex/include/basic/event_encoding/base_encoding/%.h
	mkdir -p src/obj/basic/; $(CC) $(CFLAGS_BASIC) $(OPT) -c $< -o $@

src/obj/test/basic/%.o : libs/cortex/src/test/basic/%.c libs/cortex/include/test/basic/%.h
	mkdir -p src/obj/test/basic; $(CC) $(CFLAGS_BASIC_TESTS) $(OPT) -c $< -o $@

src/obj/test/hash_table/open_hash/%.o : libs/cortex/src/hash_table/open_hash/%.c libs/cortex/include/hash_table/open_hash/%.h
	mkdir -p src/obj/test/hash_table/open_hash; $(CC) $(CFLAGS_HASH_TABLE_TESTS) $(OPT) -c $< -o $@

src/obj/test/hash_table/hash_key/bob_jenkins/%.o : libs/cortex/src/hash_table/hash_key/bob_jenkins/%.c libs/cortex/include/hash_table/hash_key/bob_jenkins/%.h
	mkdir -p src/obj/test/hash_table/hash_key/bob_jenkins; $(CC) $(CFLAGS_HASH_TABLE_TESTS) $(OPT) -c $< -o $@

src/obj/test/hash_table/%.o : libs/cortex/src/test/hash_table/%.c libs/cortex/include/test/hash_table/%.h
	mkdir -p src/obj/test/hash_table; $(CC) $(CFLAGS_HASH_TABLE_TESTS) $(OPT) -c $< -o $@

src/obj/test/predictor/%.o : src/test/myKrobe/predictor/%.c include/test/myKrobe/predictor/%.h
	mkdir -p src/obj/test/predictor; $(CC) $(CFLAGS_PREDICTOR_TESTS) $(OPT) -c $< -o $@


src/obj/predictor/%.o : libs/cortex/src/basic/event_encoding/base_encoding/%.c libs/cortex/include/basic/event_encoding/base_encoding/%.h
	mkdir -p src/obj/predictor; $(CC) $(CFLAGS_PREDICTOR_CORE) $(CFLAGS_PREDICTOR) $(OPT) -c $< -o $@

src/obj/predictor/hash_table/open_hash/%.o : libs/cortex/src/hash_table/open_hash/%.c libs/cortex/include/hash_table/open_hash/%.h
	mkdir -p src/obj/predictor/hash_table/open_hash; $(CC) $(CFLAGS_PREDICTOR) $(OPT) -c $< -o $@

src/obj/predictor/hash_table/hash_key/bob_jenkins/%.o : libs/cortex/src/hash_table/hash_key/bob_jenkins/%.c libs/cortex/include/hash_table/hash_key/bob_jenkins/%.h
	mkdir -p src/obj/predictor/hash_table/hash_key/bob_jenkins; $(CC) $(CFLAGS_PREDICTOR) $(OPT) -c $< -o $@


src/obj/predictor/%.o : libs/cortex/src/basic/%.c libs/cortex/include/basic/%.h
	mkdir -p src/obj/predictor; $(CC) $(CFLAGS_PREDICTOR_CORE) $(CFLAGS_PREDICTOR) $(OPT) -c $< -o $@

src/obj/predictor/%.o : libs/cortex/src/hash_table/hash_key/bob_jenkins/%.c libs/cortex/include/hash_table/%.h
	mkdir -p src/obj/predictor; $(CC) $(CFLAGS_PREDICTOR_CORE) $(CFLAGS_PREDICTOR) $(OPT) -c $< -o $@

src/obj/predictor/%.o : libs/cortex/src/hash_table/open_hash/%.c libs/cortex/include/hash_table/open_hash/%.h
	mkdir -p src/obj/predictor; $(CC) $(CFLAGS_PREDICTOR_CORE) $(CFLAGS_PREDICTOR) $(OPT) -c $< -o $@

src/obj/test/predictor/%.o : src/test/predictor/%.c include/test/predictor/%.h
	mkdir -p src/obj/test/predictor; $(CC) $(CFLAGS_PREDICTOR_TESTS) $(CFLAGS_PREDICTOR_CORE) $(OPT) -c $< -o $@

