Simply demo for Atlas - a database of genetic variation in bateria.  

# Usage

Add new variants to atlas

./add-new-variants-to-database.py sample.vcf --db_name :db_name --kmer :kmer_size

# AMR prediction

atlas amr -s C00003204 -1 C00003204.bam

# Make panels

./bin/make-variant-panel.py -f staph-panel.txt -g data/BX571856.1.gb data/BX571856.1.fasta




