Mykrobe-predictor
=================

This repo uses submodules so to download it, type the following

git clone --recursive https://github.com/iqbal-lab/Mykrobe-predictor.git

Then run the following command to compile libraries. 

sh install.sh 

### Compilation ###

To compile for S.aureus, run the following command 

make STAPH=1 predictor

to create this executable 

bin/Mykrobe.predictor.staph

To compile for M.tuberculosis, run the following command 

make TB=1 predictor

to create this executable 

bin/Mykrobe.predictor.tb


To download and compile on Mac

In order to install on Mac you may need to install gcc via homebrew or macports. Running 
$ gcc -v 
should return something like:
gcc version 4.9.2 (Homebrew gcc49 4.9.2_1)
not 
Apple LLVM version 6.0 (clang-600.0.57)

### Usage ###

Standard usage

Mykrobe.predictor.staph --file some_file.fastq --install_dir /path/to/Mykrobe-predictor

Finally, there are GUI versions of Mykrobe-predictor for Windows and Mac OS X, which you can download from Releases


### Paper, citation ###
We have a preprint of the paper describing Mykrobe predictor here:





.
