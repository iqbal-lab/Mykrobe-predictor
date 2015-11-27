Mykrobe-predictor
=================

This repo uses submodules so to download it, type the following

git clone --recursive https://github.com/iqbal-lab/Mykrobe-predictor.git

Then run the following command to compile libraries. 

sh install.sh 

pip install -r code_generator/requirements.txt
### Compilation ###

To compile for S.aureus, run the following command 

**make STAPH=1 predictor**

to create this executable 

bin/Mykrobe.predictor.staph

To compile for M.tuberculosis, run the following command 

**make TB=1 predictor**

to create this executable 

bin/Mykrobe.predictor.tb


[In order to install on Mac you may need to install gcc via homebrew or macports. Running  "gcc -v" 
should return something like:
"gcc version 4.9.2 (Homebrew gcc49 4.9.2_1)"
not 
"Apple LLVM version 6.0 (clang-600.0.57)"

### Usage ###

Standard usage

Mykrobe.predictor.staph --file some_file.fastq --install_dir /path/to/Mykrobe-predictor

Finally, there are GUI versions of Mykrobe-predictor for Windows and Mac OS X, which you can download from Releases

#### Force resistance predictions if species is not target species. 

By default, if we can't find S. aureus or MTBC in the data we don't show resistance predictions. 

You can force resistance predictions with --force flag

Mykrobe.predictor.staph --file some_non_staph.fastq --install_dir /path/to/Mykrobe-predictor --force


### Output ### 

The output of Mykrobe is in JSON format. An exemplar output might looks like this:

	{
		"expected_depth": "77",
		"mean_read_length": "48",
		"phylogenetics": {
			"phylo_group": {
				"Staphylococcus aureus": "82"
			},
			"species": {
				"S. aureus": "82"
			},
			"lineage": {
				"N/A": "-1"
			}
		},
		"susceptibility" :{
			"Gentamicin": "S",
			"Penicillin": "R",
			"Methicillin": "R",
			"Trimethoprim": "S",
			"Erythromycin": "R",
			"FusidicAcid": "S",
			"Ciprofloxacin": "R",
			"Rifampicin": "S",
			"Tetracycline": "S",
			"Vancomycin": "S",
			"Mupirocin": "S",
			"Clindamycin": "R(inducible)"
		},
		"called_variants" :{
			"gyrA_S84L" :{
				"R_per_cov": "100",
				"S_per_cov": "0",
				"R_median_cov": "101",
				"S_median_cov": "0",
				"conf": "1092",
			"induced_resistance": "Ciprofloxacin"
			}
		},
		"called_genes" :{
			"blaZ" :{
				"per_cov": "74",
				"median_cov": "82",
				"conf": "42298",
			"induced_resistance": "Penicillin"
			},
			"ermC" :{
				"per_cov": "99",
				"median_cov": "834",
				"conf": "2496",
			"induced_resistance": "Erythromycin"
			},
			"mecA" :{
				"per_cov": "99",
				"median_cov": "96",
				"conf": "285",
			"induced_resistance": "Methicillin"
			}
		},
		"virulence_toxins" :{
			"PVL": "negative"
		}
	}

The phylogenetics section gives results for the species identification results. In this case we confidently call S. aureus (at ~82x coverage) with no evidence of contamination from other Staphylococcal species. 

The "susceptibility" section gives the DST results. The possible options are "S", "r", "R" for "Susceptible", "minor resistant" and "Resistant" respectively for all drugs with the exception of Clindamycin. The options for Clindamycin susceptibiliy prediction are: "S", "R(constitutive)", "R(inducible)", "r(constitutive)" and "r(inducible)". The "constitutive" is implied for all the other drugs and just means that there's a direct relationship between the gene and resistance (the normal case). However, Clindamycin has a special case. When any of the ERM* genes are present the isolate will become resistant to Clindamycin if first treated with Erythromycin (hence induced) but not otherwise [http://www.ncbi.nlm.nih.gov/pubmed/16380772](http://www.ncbi.nlm.nih.gov/pubmed/16380772). (constitutive) and (inducible) distinguish between these cases. 

The "called_genes" and "called_variants" provide the evidence (if any) for the susceptibility calls. 

"R_per_cov" - Percentage of kmers for the resistant allele recovered. 

"S_per_cov" - Percentage of kmers for the susceptible (wild-type) allele recovered. 

"R_median_cov" - Median coverage on all the recovered kmers on the resistant allele. 

"S_median_cov" - Median coverage on all the recovered kmers on the susceptible allele. 


"per_cov" - Percentage of kmers recovered for the gene. 

"median_cov" - Median coverage on all the recovered kmers for the gene. 


"induced_resistance" - shows the drug resistance(s) normally associated with this genetic element. 

### Paper, citation ###
We have a preprint of the paper describing Mykrobe predictor here:
http://biorxiv.org/content/early/2015/04/26/018564
Please cite us if you use Mykrobe predictor in a publication

All analysis in this paper was done with release [v0.1.3-beta](https://github.com/iqbal-lab/Mykrobe-predictor/releases/tag/v0.1.3-beta).

### Extending Mykrobe (modifying panel or supporting a new species) ### 
Please see [README_EXTENDING_MYKROBE](https://github.com/iqbal-lab/Mykrobe-predictor/blob/extending/README_EXTENDING_MYKROBE.md).





.
