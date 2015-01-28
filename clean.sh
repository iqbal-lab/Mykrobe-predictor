rm data/skeleton_binary/staph/skeleton.k15.ctx data/skeleton_binary/tb/skeleton.k15.ctx
cd data/skeleton_binary/tb/
ls ../../tb/antibiotics/*.fa -1 > list_speciesbranches_genes_and_muts
ls ../../tb/species/*.fa -1 >> list_speciesbranches_genes_and_muts
cd ../../../
cd data/skeleton_binary/staph/
ls ../../staph/antibiotics/*.fa -1 > list_speciesbranches_genes_and_muts
ls ../../staph/species/*.fa -1 >> list_speciesbranches_genes_and_muts
ls ../../staph/antibiotics/*.fasta -1 >> list_speciesbranches_genes_and_muts
ls ../../staph/species/*.fasta -1 >> list_speciesbranches_genes_and_muts
cd ../../../