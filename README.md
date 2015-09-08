Mykrobe-predictor
=================

This repo uses submodules so to download it, type the following

git clone --recursive https://github.com/iqbal-lab/Mykrobe-predictor.git

Then run the following command to compile libraries. 

sh install.sh 

### Compilation ###

To compile for S.aureus, run the following command 

**make STAPH=1 predictor**

to create this executable 

bin/Mykrobe.predictor.staph

To compile for M.tuberculosis, run the following command 

**make TB=1 predictor**

to create this executable 

bin/Mykrobe.predictor.tb


[In order to install on Mac you may need to install gcc via homebrew or macports. Running  "gcc -v" 
should return something like:
"gcc version 4.9.2 (Homebrew gcc49 4.9.2_1)"
not 
"Apple LLVM version 6.0 (clang-600.0.57)"

### Usage ###

Standard usage

Mykrobe.predictor.staph --file some_file.fastq --install_dir /path/to/Mykrobe-predictor

Finally, there are GUI versions of Mykrobe-predictor for Windows and Mac OS X, which you can download from Releases


### Paper, citation ###
We have a preprint of the paper describing Mykrobe predictor here:
http://biorxiv.org/content/early/2015/04/26/018564
Please cite us if you use Mykrobe predictor in a publication

### Extending Mykrobe (modifying panel or supporting a new species) ### 
Please see [README_EXTENDING_MYKROBE](https://github.com/iqbal-lab/Mykrobe-predictor/blob/extending/README_EXTENDING_MYKROBE.md).





.
