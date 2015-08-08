## Extending Mykrobe

To make extending a little easier than editing the optimised C code we offer a code generator (in python) to allow you to extend Mykrobe a higher level. Please see our extending branch for development of this code.  

Below is a quick tutorial of how to use the code generator to extend Mykrobe to new species or to use it with different panels

Let's say we want to extend Mykrobe to a new species Foo barius. 

First, create a new file in /code_generator

foo_bar.py and fill it with the following template:

	import os
	from base import CodeGenerator
	from base import DrugCodeGenerator
	from base import SpeciesCodeGenerator

	species = "fbarius"
	class FooBarCodeGenerator(CodeGenerator):

	    def __init__(self):
	        self.species = species
	        super(FooBarCodeGenerator, self).__init__()  
	        
	    @property 
	    def drugs(self):
	        return [StaphDrugCodeGenerator(drug) for drug in self.drug_names ]

	class FooBarDrugCodeGenerator(DrugCodeGenerator):
	    
	    def __init__(self,name):
	        self.species = species
	        super(FooBarDrugCodeGenerator, self).__init__(name)  
	 
	            
	cg = FooBarCodeGenerator()
	cg.render_and_write_all()

	spc = SpeciesCodeGenerator(species)
	spc.render_and_write_phylo()

we also need to create a new folder in /code_generator/data called fbarius

mkdir code_generator/data/fbarius

and populate the directory with gene_to_drug.json which gives a mapping between genes and the drugs the presence of these genes induce resistance to. And virulence_genes.csv which list the names of any genes you want to type but do not nessecarily induce AMR. 

	gene_to_drug.json

	{
		"abcd": ["antimicrobial1"],
		"efgh": ["antimicrobial2", "antimicrobial3"]
	}

The code generator assumes that any gene mentioned about has a fasta file stored in data/fbarius/antibiotics/:gene_name.fa or data/fbarius/virulence/:gene_name.fa with the relevant sequence. You can have more than one version of a gene in the fasta file. 

The code generator uses the Jinga template engine to render the C code. By default we search for a template in a fbarius subdirectory and falls back to the file in /common/ if one is not found. 

To generate the Mykrobe.predictor.fbarius binary run the following commands. 

	cd code_generator
	python fbarius.py
	cd ..

	cd data/skeleton_binary/fbarius/
	ls ../../fbarius/antibiotics/*.fa > list_speciesbranches_genes_and_muts
	ls ../../fbarius/virulence/*.fa >> list_speciesbranches_genes_and_muts
	ls ../../fbarius/phylo/*/*.fa >> list_speciesbranches_genes_and_muts
	cd ../../../


	cp code_generator/rendered/include/predictor/core/* include/predictor/core/
	cp code_generator/rendered/include/predictor/fbarius/* include/predictor/fbarius/
	cp code_generator/rendered/src/predictor/core/* src/predictor/core/
	cp code_generator/rendered/src/predictor/fbarius/* src/predictor/fbarius/

	Modify the Makefile to include:

	IDIR_PREDICTOR_FBARIUS = include/predictor/fbarius
	src/obj/predictor/%.o : src/predictor/fbarius/%.c include/predictor/fbarius/%.h
		mkdir -p src/obj/predictor; $(CC) $(CFLAGS_PREDICTOR_CORE) $(CFLAGS_PREDICTOR) $(OPT) -c $< -o $@	

	and run:

	make FBARIUS=1 predictor 


If we want to extend the FBARIUS app to also look for mutations associated with resistance we need to do the following:

Add the following property to our FbariusCodeGenerator class:

    @property 
    def mutation_induced_drug_names(self):
        return ['drug1','drug2']      


So your fbarius.py file will look like the following:

	import os
	from base import CodeGenerator
	from base import DrugCodeGenerator
	from base import SpeciesCodeGenerator

	species = "fbarius"
	class FooBarCodeGenerator(CodeGenerator):

	    def __init__(self):
	        self.species = species
	        super(FooBarCodeGenerator, self).__init__()  
	        
	    @property 
	    def drugs(self):
	        return [StaphDrugCodeGenerator(drug) for drug in self.drug_names ]

	    @property 
	    def mutation_induced_drug_names(self):
	        return ['drug1','drug2']    	        

	class FooBarDrugCodeGenerator(DrugCodeGenerator):
	    
	    def __init__(self,name):
	        self.species = species
	        super(FooBarDrugCodeGenerator, self).__init__(name)  
	 
	            
	cg = FooBarCodeGenerator()
	cg.render_and_write_all()

	spc = SpeciesCodeGenerator(species)
	spc.render_and_write_phylo()


This makes the code generator look for files in data/fbarius/drug1.fa and data/fbarius/drug2.fa that have the following structure. The ref_* reads are the wildtype allele, alt_* are the alternate alleles that induce resistance. sub-N is the number of alternate alleles associated with this reference allele. alt-N gives the 1-based index of the alternate. These are followed in both cases by -abc indicating the names of the gene. 


	>ref_A999X_sub-3-abc
	ACGTTCCCGGGCCTTGTACACACCGCCCGTCACGTCATGAAAGTCGGTAACACCCGAAGC
	CAG
	>alt_A999X_alt-1-abc
	ACGTTCCCGGGCCTTGTACACACCGCCCGTCTCGTCATGAAAGTCGGTAACACCCGAAGC
	CAG
	>alt_A999X_alt-2-abc
	ACGTTCCCGGGCCTTGTACACACCGCCCGTCCCGTCATGAAAGTCGGTAACACCCGAAGC
	CAG
	>alt_A999X_alt-3-abc
	ACGTTCCCGGGCCTTGTACACACCGCCCGTCGCGTCATGAAAGTCGGTAACACCCGAAGC
	CAG
	>ref_T111X_sub-3-abc
	CGTTCCCGGGCCTTGTACACACCGCCCGTCACGTCATGAAAGTCGGTAACACCCGAAGCC
	AGT
	>alt_T111X_alt-1-abc
	CGTTCCCGGGCCTTGTACACACCGCCCGTCAAGTCATGAAAGTCGGTAACACCCGAAGCC
	AGT
	>alt_T111X_alt-2-abc
	CGTTCCCGGGCCTTGTACACACCGCCCGTCATGTCATGAAAGTCGGTAACACCCGAAGCC
	AGT
	>alt_T111X_alt-3-abc
	CGTTCCCGGGCCTTGTACACACCGCCCGTCAGGTCATGAAAGTCGGTAACACCCGAAGCC
	AGT


This is a very simple example and we'll be updating this tutorial shortly to include a more real-world example. However, all of the code generator is writing in an object orientated manner. So, if you want to overwrite a base method you can just inherit from the base class and overwrite the method. To extend the default C files you can use the Jinga templating package to do so. 
