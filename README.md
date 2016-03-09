Master : [![Build Status](https://travis-ci.com/Phelimb/atlas.svg?token=zS56Z2pmznVQKhUTxqcq&branch=master)](https://travis-ci.com/Phelimb/atlas)

Dev : [![Build Status](https://travis-ci.com/Phelimb/atlas.svg?token=zS56Z2pmznVQKhUTxqcq&branch=dev)](https://travis-ci.com/Phelimb/atlas)

## Installation

git clone **--recursive** https://github.com/Phelimb/atlas.git

### Install requirements mccortex

	cd atlas
	cd mccortex
	make	
	export PATH=$PATH:$(pwd)/bin
	cd ..

### Install Mykrobe predictor
	
### Install Mykrobe predictor with virtualenv (recommended but optional)

	1. [Install virtualenv](https://virtualenv.readthedocs.org/en/latest/installation.html)
	2. Create virtualenv 'virtualenv venv'
	3. Activate the virtualenv 'source venv/bin/activate'. You can deactivate at anytime by typing deactivate. 

	python setup.py install


## Usage

	mykrobe --help

	usage: mykrobe [-h] [--version] {predict,genotype} ...

	optional arguments:
	  -h, --help          show this help message and exit
	  --version           mykrobe version

	[sub-commands]:
	  {predict,genotype}
	    predict           Predict the sample's antibiogram
	    genotype          Genotype a sample

	----------------------------------------------------------------

	mykrobe predict --help

	usage: mykrobe predict [-h] [-k kmer] [--tmp TMP]
	                       [--skeleton_dir SKELETON_DIR] [-q] [--panel panel]
	                       [--force]
	                       sample seq [seq ...] species

	positional arguments:
	  sample                sample id
	  seq                   sequence files (fastq or bam)
	  species               species

	optional arguments:
	  -h, --help            show this help message and exit
	  -k kmer, --kmer kmer  kmer length (default:21)
	  --tmp TMP             tmp directory (default: /tmp/)
	  --skeleton_dir SKELETON_DIR
	                        directory for skeleton binaries
	  -q, --quiet           do not output warnings to stderr
	  --panel panel         variant panel (default:bradley-2015)	

### AMR prediction

	mykrobe predict tb_sample_id tb_sequence.bam tb

	mykrobe predict staph_sample_id staph_sequence.bam tb

### Output

Output is in JSON format. To convert to a less verbose tabular format use [json_to_tsv](scripts/json_to_tsv.py).

{
    "sample_id": {
        "susceptibility": {
            "Rifampicin": {
                "predict": "S"
            },
            ...
            "Streptomycin": {
                "predict": "S"
            }
        "phylogenetics": {
            "lineage": {
                "Unknown": {
                    "percent_coverage": -1,
                    "median_depth": -1
                }
            },
			...
            "species": {
                "Mycobacterium_tuberculosis": {
                    "percent_coverage": 98.0,
                    "median_depth": 53
                }
            }
        },  
        "typed_variants": {
            "rpoB_N438S-AAC761118AGT": {
                "info": {
                    "contamination_depths": [],
                    "coverage": {
                        "alternate": {
                            "percent_coverage": 47.62,
                            "median_depth": 0.0,
                            "min_depth": 47.0
                        },
                        "reference": {
                            "percent_coverage": 100.0,
                            "median_depth": 49.0,
                            "min_depth": 44.0
                        }
                    },
                    "expected_depths": [
                        56.0
                    ]
                },
                "_cls": "Call.VariantCall",
                "genotype": [
                    0,
                    0
                ],
                "genotype_likelihoods": [
                    -4.25684443365591,
                    -99999999.0,
                    -99999999.0
                ]
            },   ...               
        },			

### Change the panel for resistance prediction (TB only)
	
	mykrobe predict tb_sample_id tb_sequence.bam tb --panel walker-2015

	Walker, Timothy M., et al. "Whole-genome sequencing for prediction of Mycobacterium tuberculosis drug susceptibility and resistance: a retrospective cohort study." The Lancet Infectious Diseases 15.10 (2015): 1193-1202.

	mykrobe predict tb_sample_id tb_sequence.bam tb --panel bradley-2015

	Bradley, Phelim, et al. "Rapid antibiotic-resistance predictions from genome sequence data for Staphylococcus aureus and Mycobacterium tuberculosis." Nature communications 6 (2015).




Add new variants to atlas (requires mongod running in background)

	mykrobe add sample.vcf --db_name :db_name --kmer :kmer_size

### Make panels

	mykrobe make-probes -f example-data/staph-panel.txt -g BX571856.1.gb BX571856.1.fasta

	mykrobe make-probes -f example-data/tb-walker-2015-panel.txt -g data/NC_000962.3.gb data/NC_000962.3.fasta

### Dump variant probes (requires mongod running in background)

	mykrobe dump-probes data/NC_000962.3.fasta > panel_tb_k31.fasta

### Genotype using these variants

	mykrobe genotype panel_tb_k31.fasta 31 -s 10564-01 -1 /data2/users/phelim/data/tb/atlas/fastq/10564-01/10564-01.fastq.gz

	mykrobe genotype tb-amr-probes.fasta 31 -s 10564-01 -1 /data2/users/phelim/data/tb/atlas/fastq/10564-01/10564-01.fastq.gz

### Genotype a panel of genes

	mykrobe genotype ~/git/atlas-core/data/panels/staph-amr-genes.fasta 31 -s C00001283 -1 /data2/users/phelim/data/staph/atlas/bams/C00001283.bam


## Extending Mykrobe

Install mongodb

	Follow instructions at https://docs.mongodb.org/manual/installation/


