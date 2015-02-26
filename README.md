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




Standard usage

Mykrobe.predictor.staph --file some_file.fastq --install_dir /path/to/myKrobe-predictor

Add --progress
if you want to see progress.
