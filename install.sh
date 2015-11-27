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

cd string_buffer
make clean
make
cd ..

cd seq_file
make clean
#make STRING_BUF_PATH="../string_buffer" HTS_PATH="../htslib"
make STRING_BUF_PATH="../string_buffer" HTS_PATH="../htslib"


cd ..

cd ..

./make.sh