ifndef CC
	CC=gcc
endif

HTSLIB_ABS=$(shell readlink -f $(HTSLIB))

all: tools benchmarks dev

tools: htslib
	cd tools; make HTSLIB=$(HTSLIB_ABS)
benchmarks: htslib
	cd benchmarks; make HTSLIB=$(HTSLIB_ABS)
dev: htslib
	cd dev; make HTSLIB=$(HTSLIB_ABS)

htslib:
	if [[ '$(HTSLIB)' == '' ]]; \
	then echo "Error: Please pass HTSLIB=... with path to htslib"; exit 1; fi

clean:
	rm -rf bin/*
	cd benchmarks; make clean
	cd dev; make clean

.PHONY: all clean tools benchmarks dev htslib
