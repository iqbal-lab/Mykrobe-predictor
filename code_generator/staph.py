import os
from utils import unique,flatten
from base import CodeGenerator
from base import DrugCodeGenerator
from base import GeneCodeGenerator

class StaphCodeGenerator(CodeGenerator):

    def __init__(self):
        self.species = "staph"
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
        self.species = "staph"
        super(StaphDrugCodeGenerator, self).__init__(name)  

    @property 
    def has_epistatic_muts(self):
        if self.name in ['FusidicAcid','Rifampicin']:
            return True
        else:
            return False  

    @property
    def genes_resistance_induced_by(self):
        genes = []
        for gene,drugs in self.gene_enum_to_drug_name.iteritems():
            if self.name in drugs:
                genes.append(StaphGene(gene))
        return genes                             
        

class StaphGene(GeneCodeGenerator):
    
    def __init__(self,name):
        self.species = "staph" 
        super(StaphGene, self).__init__(name=name)                    
                               



cg = StaphCodeGenerator()
cg.render_and_write_all()
