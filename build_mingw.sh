#!/bin/bash


cp Makefile.mingw Makefile



export CFLAGS="-fPIC -pie -I/usr/include -O3"
export CPPFLAGS="$CFLAGS"

export LDFLAGS="-rdynamic -L/usr/lib"


make STAPH=1 predictor
make TB=1 predictor

