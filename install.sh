#!/bin/bash

if [ ! -d libs ]
then
  echo "Directory libs does not exist" 1>&2
  exit 0
fi

cd libs

cd htslib
make clean
make
cd ..

cd gsl-1.15
make clean
./configure
make
cd ..

cd string_buffer
make clean
make
cd ..

cd seq_file
make clean
make STRING_BUF_PATH=../string_buffer HTS_PATH=../htslib
cd ..

cd ..

cd scripts/analyse_variants/needleman_wunsch/
make clean
make
cd ..

cd vcf-hack
make clean
make
cd ..
cd ..
cd ..
