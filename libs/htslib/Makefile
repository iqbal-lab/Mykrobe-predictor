UNAME=$(shell uname)
ifeq ($(UNAME),Darwin)
	MAC=1
else
	ifeq ($(UNAME),CYGWIN_NT-6.1)
	WIN=1
	endif
endif

CC=			gcc
ifdef WIN
	CFLAGS=		-m64 -D__mingw__ -g -Wall -Wc++-compat -O2
	LIBLIST = -Lhtslib -lhts -lpthread -lz -lm -lws2_32
else
	CFLAGS=		-g -Wall -Wc++-compat -O2
	LIBLIST = -Lhtslib -lhts -lpthread -lz -lm
endif
 

DFLAGS=
OBJS=		main.o samview.o vcfview.o bamidx.o bcfidx.o bamshuf.o bam2fq.o tabix.o \
			abreak.o bam2bed.o vcfcheck.o vcfisec.o vcfmerge.o
INCLUDES=	-Ihtslib
PROG=		htscmd

.SUFFIXES:.c .o
.PHONY:all lib

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

lib:
		cd htslib; $(MAKE) CC="$(CC)" CFLAGS="$(CFLAGS)" libhts.a || exit 1; cd ..

htscmd:lib $(OBJS)
		$(CC) $(CFLAGS) -o $@ $(OBJS) $(LIBLIST)


clean:
		rm -fr gmon.out *.o a.out *.dSYM *~ $(PROG); cd htslib; $(MAKE) clean; cd ..
