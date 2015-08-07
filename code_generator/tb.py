import os
from base import CodeGenerator
from base import DrugCodeGenerator
from base import SpeciesCodeGenerator

species = "tb"
class TBCodeGenerator(CodeGenerator):

    def __init__(self):
        self.species = species
        super(TBCodeGenerator, self).__init__()  
        
    @property 
    def drugs(self):
        return [TBDrugCodeGenerator(drug) for drug in self.drug_names ]

    @property 
    def mutation_induced_drug_names(self):
        return ["amikacin","capreomycin","ethambutol","isoniazid",
                "kanamycin","pyrazinamide","quinolones","rifampicin",
                "streptomycin"]                 
                

class TBDrugCodeGenerator(DrugCodeGenerator):
    
    def __init__(self,name):
        self.species = species
        super(TBDrugCodeGenerator, self).__init__(name)  
 
            
cg = TBCodeGenerator()
cg.render_and_write_all()

spc = SpeciesCodeGenerator(species)
spc.render_and_write_phylo()