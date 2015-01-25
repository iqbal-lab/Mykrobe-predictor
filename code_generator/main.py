import jinja2
import os
import json
from utils import unique,flatten
from mutations import MutationFasta

class CodeGenerator(object):

    def __init__(self):
        self.template_dir = os.path.join(os.path.dirname(__file__), 'templates')
        self.jinja_env = jinja2.Environment(loader = jinja2.FileSystemLoader(self.template_dir),
                               autoescape = False)        
       
    def render_str(self, template, **params):
        try:
            t = self.jinja_env.get_template(template)
        except AttributeError:
            return None
        selfer = self
        return t.render(params,selfer = selfer)

    def load_gene_enum_to_drug_name(self):
        with open('data/%s/gene_to_drug.json'  % self.species ,'r') as infile:
            d =  json.load(infile)
        return d          



class StaphCodeGenerator(CodeGenerator):

    def __init__(self):
        self.species = "staph"
        self.gene_enum_to_drug_name = self.load_gene_enum_to_drug_name()
        self.template_dir = os.path.join(os.path.dirname(__file__), 'templates/' )
        self.render_dir = os.path.join(os.path.dirname(__file__), 'rendered/' )

        self.jinja_env = jinja2.Environment(loader = jinja2.FileSystemLoader(self.template_dir),
                               autoescape = False)
    @property 
    def genes(self):
        return unique(self.gene_enum_to_drug_name.keys() + self.virulence_genes)

    @property
    def virulence_genes(self):
        return ['luk']

    @property 
    def drugs(self):

        return [StaphDrugCodeGenerator(drug) for drug in self.drug_names ]

    @property 
    def drug_names(self):
        drug_names = unique(self.gene_induced_drug_names + self.mutation_induced_drug_names) 
        drug_names.insert(0, drug_names.pop(drug_names.index("Erythromycin")))
        return drug_names                        

    @property 
    def gene_induced_drug_names(self):
        return unique(flatten(self.gene_enum_to_drug_name.values()))

    @property 
    def mutation_induced_drug_names(self):
        return ['Rifampicin','Ciprofloxacin']

    def render_antibiotics_src(self):
        return self.render_str('src/predictor/staph/antibiotics.c')

    def render_antibiotics_include(self):
        return self.render_str('include/predictor/staph/antibiotics.h')        

    def render_and_write_antibiotics(self):
        st = self.render_antibiotics_src()
        with open(os.path.join(self.render_dir,'src/predictor/staph/antibiotics.c'),'w') as outfile:
            outfile.write(st)
        st = self.render_antibiotics_include()
        with open(os.path.join(self.render_dir,'include/predictor/staph/antibiotics.h'),'w') as outfile:
            outfile.write(st) 

    def render_known_mutations_src(self):
        return self.render_str('src/predictor/core/known_mutations.c')

    def render_known_mutations_include(self):
        return self.render_str('include/predictor/core/known_mutations.h')        

    def render_and_write_known_mutations(self):
        st = self.render_known_mutations_src()
        with open(os.path.join(self.render_dir,'src/predictor/core/known_mutations.c'),'w') as outfile:
            outfile.write(st)
        st = self.render_known_mutations_include()
        with open(os.path.join(self.render_dir,'include/predictor/core/known_mutations.h'),'w') as outfile:
            outfile.write(st) 

    def render_gene_presence_src(self):
        return self.render_str('src/predictor/core/gene_presence.c')

    def render_gene_presence_include(self):
        return self.render_str('include/predictor/core/gene_presence.h')        

    def render_and_write_gene_presence(self):
        st = self.render_gene_presence_src()
        with open(os.path.join(self.render_dir,'src/predictor/core/gene_presence.c'),'w') as outfile:
            outfile.write(st)
        st = self.render_gene_presence_include()
        with open(os.path.join(self.render_dir,'include/predictor/core/gene_presence.h'),'w') as outfile:
            outfile.write(st) 

    def render_main_src(self):
        return self.render_str('src/predictor/staph/main.c')

    def render_and_write_main(self):
        st = self.render_main_src()
        with open(os.path.join(self.render_dir,'src/predictor/staph/main.c'),'w') as outfile:
            outfile.write(st)    

    def render_and_write_all(self):
        self.render_and_write_antibiotics() 
        self.render_and_write_known_mutations() 
        self.render_and_write_gene_presence()
        self.render_and_write_main()
        

    @property 
    def all_mutations(self):
        muts = []
        for drug in self.drugs:
            muts.extend(drug.mut_list)
        return muts

    @property 
    def all_mut_genes(self):
        mut_genes = []
        for drug in self.drugs:
            mut_genes.extend(drug.mut_gene_list)
        return mut_genes        




class DrugCodeGenerator(CodeGenerator):

    def __init__(self,name):
        self.name = name
        self.gene_enum_to_drug_name = self.load_gene_enum_to_drug_name()
        self.mutation_fasta = MutationFasta('../data/staph/antibiotics/%s' % self.name.lower() )

    def __str__(self):
        return self.name   

    @property
    def genes_resistance_induced_by(self):
        genes = []
        for gene,drugs in self.gene_enum_to_drug_name.iteritems():
            if self.name in drugs:
                genes.append(StaphGene(gene))
        return genes 

    @property 
    def has_epistatic_muts(self):
        if self.name in ['FusidicAcid','Rifampicin']:
            return True
        else:
            return False

    @property 
    def epistaic_file_path(self):
        return os.path.join(self.template_dir,'epistatic',self.name)

    @property 
    def num_mutations(self):
        return self.mutation_fasta.num_mutations

    @property 
    def first_mut(self):
        return self.mutation_fasta.first_mut

    @property 
    def last_mut(self):
        return self.mutation_fasta.last_mut  

    @property 
    def mut_list(self):
        return self.mutation_fasta.mut_list 

    @property 
    def mut_gene_list(self):
        return self.mutation_fasta.mut_gene_list                    

    @property 
    def num_genes(self):
        return len(self.genes_resistance_induced_by)                   

class StaphDrugCodeGenerator(DrugCodeGenerator):
    
    def __init__(self,name):
        self.species = "staph"
        super(StaphDrugCodeGenerator, self).__init__(name)            
        

class GeneCodeGenerator(CodeGenerator):

    def __init__(self,name):
        self.name = name 
        self.gene_enum_to_drug_name = self.load_gene_enum_to_drug_name()

    def __str__(self):
        return self.name     

        
class StaphGene(GeneCodeGenerator):
    
    def __init__(self,name):
        self.species = "staph" 
        super(StaphGene, self).__init__(name=name)                    
                               


cg = StaphCodeGenerator()
cg.render_and_write_all()
# print cg.render_map_gene_to_drug_resistance()
# print cg.render_antibiotics()
# print cg.is_drug_susceptible(drug="Gentamicin")
# print cg.is_drug_susceptible(drug="FusidicAcid")