sh render.sh
rm data/skeleton_binary/*/skeleton.k15.ctx

cd data/skeleton_binary/staph/
ls ../../../data/staph/*/*.fa > list_speciesbranches_genes_and_muts
ls ../../../data/staph/*/*.fasta >> list_speciesbranches_genes_and_muts

cd ../../../
cd data/skeleton_binary/tb/
ls ../../../data/tb/*/*.fa > list_speciesbranches_genes_and_muts
ls ../../../data/tb/*/*.fasta >> list_speciesbranches_genes_and_muts
cd ../../../

make TB=1 predictor 
make STAPH=1 predictor

cd bin
./Mykrobe.predictor.tb  --file ../data/tb/species/MTBC.fa --install_dir ../
./Mykrobe.predictor.staph  --file ../data/staph/species/Saureus.fasta --install_dir ../
