Master : [![Build Status](https://travis-ci.com/Phelimb/atlas.svg?token=zS56Z2pmznVQKhUTxqcq&branch=master)](https://travis-ci.com/Phelimb/atlas)

Dev : [![Build Status](https://travis-ci.com/Phelimb/atlas.svg?token=zS56Z2pmznVQKhUTxqcq&branch=dev)](https://travis-ci.com/Phelimb/atlas)

## Installation

git clone **--recursive** git@github.com:Phelimb/atlas.git

### Install requirements mccortex
**NB: You must install the version of mccortex that comes with this repostitory**


	cd atlas
	cd mccortex
	make	
	export PATH=$PATH:$(pwd)/bin
	cd ..


## Install Mykrobe predictor

### Download probes

	./scripts/download-probes.sh
	
### Install Mykrobe predictor with virtualenv (recommended but optional)

#### Install virtualenv

	https://virtualenv.readthedocs.org/en/latest/installation.html

#### Create virtualenv 

	virtualenv venv

#### Activate the virtualenv

	source venv/bin/activate

You can deactivate at anytime by typing 'deactivate'. 


#### Install Mykrobe predictor


	python setup.py install


## Usage

	mykrobe predict predict --help
	usage: mykrobe predict [-h] [-k kmer] [--tmp TMP]
	                       [--skeleton_dir SKELETON_DIR]
	                       [--mccortex31_path MCCORTEX31_PATH] [-q]
	                       [--panel panel] [--force]
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
	  --mccortex31_path MCCORTEX31_PATH (default:mccortex31)
	                        Path to mccortex31
	  -q, --quiet           do not output warnings to stderr
	  --panel panel         variant panel (default:bradley-2015)
	  --force

### AMR prediction

	mykrobe predict tb_sample_id tb_sequence.bam/fq tb

	mykrobe predict staph_sample_id staph_sequence.bam/fq staph

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

> Walker, Timothy M., et al. "Whole-genome sequencing for prediction of Mycobacterium tuberculosis drug susceptibility and resistance: a retrospective cohort study." The Lancet Infectious Diseases 15.10 (2015): 1193-1202.

	mykrobe predict tb_sample_id tb_sequence.bam tb --panel bradley-2015

> Bradley, Phelim, et al. "Rapid antibiotic-resistance predictions from genome sequence data for Staphylococcus aureus and Mycobacterium tuberculosis." Nature communications 6 (2015).

### Paper, citation 

> [Bradley, Phelim, et al. "Rapid antibiotic-resistance predictions from genome sequence data for Staphylococcus aureus and Mycobacterium tuberculosis."Nature communications 6 (2015).](http://www.nature.com/ncomms/2015/151221/ncomms10063/full/ncomms10063.html)

Please cite us if you use Mykrobe predictor in a publication

All analysis in this paper was done with release [v0.1.3-beta](https://github.com/iqbal-lab/Mykrobe-predictor/releases/tag/v0.1.3-beta).


### Common issues

mccortex fails to make. 

Likely problem: Submodules have not been pulled with the repo. 

Solution : Run 
	
	git pull && git submodule init && git submodule update
	cd mccortex && git submodule init && git submodule update
	make

