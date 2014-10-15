myKrobe-predictor
=================



------ Compilation ---------

To compile for Staphylococcus

make STAPH=1 predictor

creates an executable bin/Mykrobe.predictor.staph

To compile for Tuberculosis

make TB=1 predictor

creates an executable bin/Mykrobe.predictor.tb




---- Usage -----




Two typical uses

1. Mykrobe.predictor.staph --file some_file.fastq --install_dir /path/to/myKrobe-predictor

or, to build the full genome and then draw inferences from that

2. Mykrobe.predictor.tb --file some_file.fastq --method WgAssemblyAndGenotyping --install_dir /path/to/myKrobe-predictor

Add 
--format JSON
if you want output in JSON format with a bit more detail.

Add --progress
if you want to see progress.
