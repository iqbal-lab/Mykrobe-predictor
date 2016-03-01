Master : [![Build Status](https://travis-ci.com/Phelimb/atlas.svg?token=zS56Z2pmznVQKhUTxqcq&branch=master)](https://travis-ci.com/Phelimb/atlas)

Dev : [![Build Status](https://travis-ci.com/Phelimb/atlas.svg?token=zS56Z2pmznVQKhUTxqcq&branch=dev)](https://travis-ci.com/Phelimb/atlas)

## Installation

git clone **--recursive** https://github.com/Phelimb/atlas.git

### Install requirements mccortex and mongodb

	cd atlas
	cd mccortex
	make	
	export PATH=$PATH:$(pwd)/bin
	cd ..

Install mongodb

	Follow instructions at https://docs.mongodb.org/manual/installation/


## Usage

Add new variants to atlas (requires mongod running in background)

	atlas add sample.vcf --db_name :db_name --kmer :kmer_size

### Make panels

	atlas make-probes -f example-data/staph-panel.txt -g BX571856.1.gb BX571856.1.fasta

	atlas make-probes -f example-data/tb-walker-2015-panel.txt -g data/NC_000962.3.gb data/NC_000962.3.fasta

### Dump variant probes (requires mongod running in background)

	atlas dump-probes data/NC_000962.3.fasta > panel_tb_k31.fasta

### Genotype using these variants

	atlas genotype panel_tb_k31.fasta 31 -s 10564-01 -1 /data2/users/phelim/data/tb/atlas/fastq/10564-01/10564-01.fastq.gz

	atlas genotype tb-amr-probes.fasta 31 -s 10564-01 -1 /data2/users/phelim/data/tb/atlas/fastq/10564-01/10564-01.fastq.gz

### AMR prediction

	atlas predict -s C00003204 -1 C00003204.bam


