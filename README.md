myKrobe-predictor
=================



------ Compilation ---------

To compile for Staphylococcus

make STAPH=1 predictor

creates an executable bin/Mykrobe.predictor.staph

To compile for Tuberculosis

make TB=1 predictor

creates an executable bin/Mykrobe.predictor.tb


### On Mac ###

In order to install on Mac you may need to install gcc via homebrew or macports. 

Running 
$ gcc -v 

should return something like:

gcc version 4.9.2 (Homebrew gcc49 4.9.2_1)

not 

Apple LLVM version 6.0 (clang-600.0.57)

---- Usage -----




Standard usage

Mykrobe.predictor.staph --file some_file.fastq --install_dir /path/to/myKrobe-predictor

Add --progress
if you want to see progress.
