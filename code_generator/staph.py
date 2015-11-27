#! /usr/bin/env python
import os
from base import unique,flatten
from base import CodeGenerator
from base import DrugCodeGenerator
from base import SpeciesCodeGenerator

species = "staph"
class StaphCodeGenerator(CodeGenerator):

    def __init__(self):
        self.species = species
        self.genome_size = 2800000
        self.kmer_size = 15
        self.species_long_name = "S. aureus"        
        super(StaphCodeGenerator, self).__init__()  
        
    @property 
    def drug_names(self):
        drug_names = unique(self.gene_induced_drug_names + self.mutation_induced_drug_names) 
        drug_names.insert(0, drug_names.pop(drug_names.index("Erythromycin")))
        return drug_names  

    @property 
    def drugs(self):
        return [StaphDrugCodeGenerator(drug) for drug in self.drug_names ]

    @property 
    def mutation_induced_drug_names(self):
        return ['Rifampicin','Ciprofloxacin']                  
                

class StaphDrugCodeGenerator(DrugCodeGenerator):
    
    def __init__(self,name):
        self.species = species
        super(StaphDrugCodeGenerator, self).__init__(name)  

    @property 
    def has_epistatic_muts(self):
        if self.name in ['FusidicAcid','Rifampicin']:
            return True
        else:
            return False  
            
cg = StaphCodeGenerator()
cg.render_and_write_all()

spc = SpeciesCodeGenerator(species)
spc.render_and_write_phylo()
