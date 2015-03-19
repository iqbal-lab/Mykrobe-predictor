#!/bin/bash


cp Makefile.master Makefile



make STAPH=1 predictor
make TB=1 predictor

