set -e
make clean
rm -f data/skeleton_binary/*/skeleton.k15.ctx
mkdir -p data/skeleton_binary/{staph,tb}
mkdir -p src/predictor/{staph,tb}
cd code_generator
python staph.py
cd ..
set +e
cp code_generator/rendered/include/predictor/core/* include/predictor/core/
cp code_generator/rendered/include/predictor/staph/* include/predictor/staph/
cp code_generator/rendered/src/predictor/core/* src/predictor/core/
cp code_generator/rendered/src/predictor/staph/* src/predictor/staph/

cd data/skeleton_binary/staph/
ls ../../../data/staph/*/*.fa > list_speciesbranches_genes_and_muts
ls ../../../data/staph/*/*.fasta >> list_speciesbranches_genes_and_muts
ls ../../../data/staph/*/*/*.fasta >> list_speciesbranches_genes_and_muts
ls ../../../data/staph/*/*/*.fasta >> list_speciesbranches_genes_and_muts
cd ../../../

make STAPH=1 predictor
cd bin
./Mykrobe.predictor.staph  --file ../data/staph/phylo/species/Saureus.fasta --install_dir ../
cd ..


cd code_generator
python tb.py
cd ..

cd data/skeleton_binary/tb/
ls ../../tb/antibiotics/*.fa > list_speciesbranches_genes_and_muts
ls ../../tb/virulence/*.fa >> list_speciesbranches_genes_and_muts
ls ../../tb/phylo/*/*.fa >> list_speciesbranches_genes_and_muts
cd ../../../


cp code_generator/rendered/include/predictor/core/* include/predictor/core/
cp code_generator/rendered/include/predictor/tb/* include/predictor/tb/
cp code_generator/rendered/src/predictor/core/* src/predictor/core/
cp code_generator/rendered/src/predictor/tb/* src/predictor/tb/

make TB=1 predictor 
cd bin
./Mykrobe.predictor.tb  --file ../data/tb/phylo/phylo_group/MTBC.fa --install_dir ../
cd ..
