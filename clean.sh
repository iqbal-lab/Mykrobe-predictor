rm data/skeleton_binary/staph/skeleton.k15.ctx data/skeleton_binary/tb/skeleton.k15.ctx
cd data/skeleton_binary/tb/
ls ../../tb/antibiotics/*.fa -1 > list_speciesbranches_genes_and_muts
ls ../../tb/species/*.fa -1 >> list_speciesbranches_genes_and_muts
ls ../../tb/antibiotics/*.fasta -1 >> list_speciesbranches_genes_and_muts
ls ../../tb/species/*.fasta -1 >> list_speciesbranches_genes_and_muts


cd ../../../
cd data/skeleton_binary/staph/
ls ../../staph/antibiotics/*.fa -1 > list_speciesbranches_genes_and_muts
ls ../../staph/species/*.fa -1 >> list_speciesbranches_genes_and_muts
ls ../../staph/antibiotics/*.fasta -1 >> list_speciesbranches_genes_and_muts
ls ../../staph/species/*.fasta -1 >> list_speciesbranches_genes_and_muts
cd ../../../

cd code_generator
~/banyanvenv/bin/python main.py
cd ..
cp code_generator/rendered/src/predictor/staph/* src/predictor/staph/
cp code_generator/rendered/include/predictor/staph/* include/predictor/staph
cp code_generator/rendered/src/predictor/core/* src/predictor/core/
cp code_generator/rendered/include/predictor/core/* include/predictor/core

make TB=1 predictor && make STAPH=1 predictor

bin/Mykrobe.predictor.staph --install_dir ~/git/myKrobe-predictor/ --file data/staph/species/Saureus.fasta & 
bin/Mykrobe.predictor.tb --install_dir ~/git/myKrobe-predictor/ --file data/tb/species/MTBC.fa
