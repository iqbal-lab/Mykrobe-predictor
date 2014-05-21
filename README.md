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

or, to run in fixed memory using a skeleton binary

2. Mykrobe.predictor.tb --file some_file.fastq --method InSilicoOligos --install_dir /path/to/myKrobe-predictor
