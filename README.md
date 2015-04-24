Mykrobe-predictor
=================

This repo uses submodules so to download it, type the following

git clone --recursive https://github.com/iqbal-lab/Mykrobe-predictor.git

Then run the following command to compile libraries. 

sh install.sh 

### Compilation ###

To compile for S.aureus, the following command creates an executable bin/Mykrobe.predictor.staph

make STAPH=1 predictor

To compile for M.tuberculosis, the following command creates an executable bin/Mykrobe.predictor.tb

make TB=1 predictor

See below for compiling on a Mac

### Usage ###

Standard usage

Mykrobe.predictor.staph --file some_file.fastq --install_dir /path/to/myKrobe-predictor

Finally, there are GUI versions of Mykrobe-predictor for Windows and Mac OS X, which you can download from Releases


### Paper, citation ###
We have a preprint of the paper describing Mykrobe predictor here:


--- To download and compile on Mac ---

In order to install on Mac you may need to install gcc via homebrew or macports. 

Running 
$ gcc -v 

should return something like:

gcc version 4.9.2 (Homebrew gcc49 4.9.2_1)

not 

Apple LLVM version 6.0 (clang-600.0.57)


.
