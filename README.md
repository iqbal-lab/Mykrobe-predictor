Master : [![Build Status](https://travis-ci.org/Phelimb/Mykrobe-predictor.svg?branch=master)](https://travis-ci.org/Phelimb/Mykrobe-predictor)

Dev : [![Build Status](https://travis-ci.org/Phelimb/Mykrobe-predictor.svg?branch=m2)](https://travis-ci.org/Phelimb/Mykrobe-predictor)

## Installation

## Python
	
	(sudo) pip install git+https://github.com/Phelimb/atlas@f90622e5a029cd3fa28695aab5ca4137002c3e86 ## Required dependancy. 
	(sudo) pip install mykrobe

We recommend that you use a virtualenv to install mykrobe. If you haven't used virtualenv before please read a guide [here](virutualenvREADME.md).

## Docker 

	docker run phelimb/mykrobe_predictor mykrobe --help

## Usage

	mykrobe predict --help
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

	mykrobe predict tb_sample_id tb -1 tb_sequence.bam/fq

	mykrobe predict staph_sample_id staph -1 staph_sequence.bam/fq

e.g.

	mykrobe predict ERR117639 /download/ena/ERR117639*.gz tb

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
	
	mykrobe predict tb_sample_id  tb --panel walker-2015 -1 tb_sequence.bam

> Walker, Timothy M., et al. "Whole-genome sequencing for prediction of Mycobacterium tuberculosis drug susceptibility and resistance: a retrospective cohort study." The Lancet Infectious Diseases 15.10 (2015): 1193-1202.

	mykrobe predict tb_sample_id  tb --panel bradley-2015 -1 tb_sequence.bam

> Bradley, Phelim, et al. "Rapid antibiotic-resistance predictions from genome sequence data for Staphylococcus aureus and Mycobacterium tuberculosis." Nature communications 6 (2015).

### Choosing a panel for TB

In short, bradley-2015 is more specific, walker-2015 more sensitive. You can see a comparision of the results on ~3000 MTUB sequence [here](example-data/walker-vs-bradley.md)

### Paper, citation 

> [Bradley, Phelim, et al. "Rapid antibiotic-resistance predictions from genome sequence data for Staphylococcus aureus and Mycobacterium tuberculosis."Nature communications 6 (2015).](http://www.nature.com/ncomms/2015/151221/ncomms10063/full/ncomms10063.html)

Please cite us if you use Mykrobe predictor in a publication

All analysis in this paper was done with release [v0.1.3-beta](https://github.com/iqbal-lab/Mykrobe-predictor/releases/tag/v0.1.3-beta).



### Genotype on a catalog 

	mykrobe genotype [-h] [-k kmer] [--tmp TMP] [--keep_tmp]
	                        [--skeleton_dir SKELETON_DIR]
	                        [--mccortex31_path MCCORTEX31_PATH] [-t THREADS]
	                        [--expected_depth EXPECTED_DEPTH] [-1 seq [seq ...]]
	                        [-c ctx] [-f] [-q] [--ignore_filtered IGNORE_FILTERED]
	                        sample probe_set

	positional arguments:
	  sample                sample id
	  probe_set             probe_set

	optional arguments:
	  -h, --help            show this help message and exit
	  -1 seq [seq ...], --seq seq [seq ...]
	                        sequence files (fasta,fastq,bam)
	  -c ctx, --ctx ctx     cortex graph binary	                        	  
	  -k kmer, --kmer kmer  kmer length (default:21)
	  --tmp TMP             tmp directory (default: /tmp/)
	  --keep_tmp            Dont remove tmp files
	  --skeleton_dir SKELETON_DIR
	                        directory for skeleton binaries
	  --mccortex31_path MCCORTEX31_PATH
	                        Path to mccortex31
	  -t THREADS, --threads THREADS
	                        threads
	  --expected_depth EXPECTED_DEPTH
	                        expected depth
	  -f, --force           force
	  -q, --quiet           do not output warnings to stderr

	  e.g. 

	   head example-data/staph-amr-bradley_2015.fasta

		>mecA?name=mecA&version=1
		ATGAATATAGTTGAAAATGAAATATGTATAAGA...ATAAAAGGACTTATAAAGATTGA
		>mecA?name=mecA&version=2
		ATGAATATAGTTGAAAATGAAATATGTATAAGA...TGAAGATTTGCCAGAACATGAAT	   
		>fusA?name=fusA&version=1
		ATGAATATAGTTGAAAATGAAATATGTATAAGA...TGAAGATTTGCCAGAACATGAAT	 	  

	   mykrobe genotype sample_id example-data/staph-amr-bradley_2015.fasta -1 seq.fq 

		

	{
	    "sample_id": {
	        "files": [
	            "seq.fq "
	        ],
	        "kmer": 21,
	        "sequence_calls": {
	            "mecA": {
	                "info": {
	                    "copy_number": 0.0,
	                    "contamination_depths": [],
	                    "coverage": {
	                        "percent_coverage": 0.0,
	                        "median_depth": 0.0,
	                        "min_non_zero_depth": 0.0
	                    },
	                    "expected_depths": [
	                        1
	                    ]
	                },
	                "_cls": "Call.SequenceCall",
	                "genotype": [
	                    0,
	                    0
	                ],
	                "genotype_likelihoods": [
	                    -0.001,
	                    -99999999.0,
	                    -99999999.0
	                ]
	            },
	            "fusA": {
	                "info": {
	                    "copy_number": 1.0276923076923077,
	                    "contamination_depths": [],
	                    "version": "10",
	                    "coverage": {
	                        "percent_coverage": 100.0,
	                        "median_depth": 167.0,
	                        "min_non_zero_depth": 116.0
	                    },
	                    "expected_depths": [
	                        162.5
	                    ]
	                },
	                "_cls": "Call.SequenceCall",
	                "genotype": [
	                    1,
	                    1
	                ],
	                "genotype_likelihoods": [
	                    -994.7978064088725,
	                    -349.45246450237215,
	                    -10.95808091830304
	                ]
	            },	            	   
	        ....
	    }
	}

### Comparing results

	./scripts/compare.py 

	usage: compare.py [-h] [--analysis {summary,table,diff}]
	                  [--ana1 ANA1 [ANA1 ...]] [--ana2 ANA2 [ANA2 ...]]
	                  [--format {short,long}]
	                  truth



e.g. Get a summary of a commit vs the truth.  

	Truth JSON is in the form: 

	{
	    "TRL0079551-S12": {
	        "susceptibility": {
	            "Rifampicin": {
	                "predict": "S"
	            },
	            "Capreomycin": {
	                "predict": "S"
	            },
	            ...
	            "Quinolones": {
	                "predict": "S"
	            }
	        }
	    },


	}	        		

	compare.py --analysis summary truth.json --ana1 *.json

Output 

	Ana1
	Drug    Total   FN(R)   FP(S)   VME     ME      sensitivity     specificity     PPV     NPV
	Rifampicin      99      1 (24)  1 (75)  4.2%       1.3%       95.8%      98.7%      95.8%      98.7% 
	Capreomycin     14      0 (4)   0 (10)  0.0%       0.0%       100.0%     100.0%     100.0%     100.0% 
	all     539     23 (128)        14 (411)        18.0%      3.4%       82.0%      96.6%      88.2%      94.5% 

	...	


Other options for --analysis are `diff` and `table` which will report a table showing samples where two analyses differed and a full table of comparisions respectively. 

e.g. 

	sample  drug    truth   ana1    ana2
	10205-03        Quinolones      NA      R       S
	10091-01        Isoniazid      S      R       S


### Common issues

* mccortex install fails

mccortex should be installed automatically with the atlas dependancy. If for some reason it doesn't you can run a manual install:

	git clone --recursive https://github.com/iqbal-lab/Mykrobe-predictor.git
	cd Mykrobe-predictor
	cd mccortex
	make    
	export PATH=$PATH:$(pwd)/bin
	cd ..

	

