import os
import sys
try:
    import jinja2
except:
    print "Missing python package. Please run 'pip install -r code_generator/requriements.txt'"
    sys.exit()
import json
import csv
import string
import glob
import logging
from jinja2.exceptions import TemplateNotFound



from mutations import MutationFasta

def make_safe_string(s):
    valid_chars = "%s%s" % (string.ascii_letters, string.digits)
    return ''.join(c for c in s if c in valid_chars)

def unique(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]

def flatten(l):
    return [item for sublist in l for item in sublist]


class CodeGenerator(object):

    def __init__(self):
        self.gene_enum_to_drug_name = self.load_gene_enum_to_drug_name()
        self.virulence_genes = unique(self.load_virulence_genes())
        self.template_dir = os.path.join(os.path.dirname(__file__), 'templates/' )
        self.render_dir = os.path.join(os.path.dirname(__file__), 'rendered/' )

        self.jinja_env = jinja2.Environment(loader = jinja2.FileSystemLoader(self.template_dir),
                               autoescape = False,
                               extensions=['jinja2.ext.with_'])  

    @property 
    def genes(self):
        return unique(self.gene_enum_to_drug_name.keys() + self.virulence_genes)

    @property 
    def drugs(self):
        raise NotImplementedError("Implemented in child")

    @property 
    def mutation_induced_drug_names(self):
        return []   

    @property 
    def drug_names(self):
        drug_names = unique(self.gene_induced_drug_names + self.mutation_induced_drug_names) 
        return drug_names                        

    @property 
    def gene_induced_drug_names(self):
        return unique(flatten(self.gene_enum_to_drug_name.values()))
       
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

    def load_virulence_genes(self):
        with open('data/%s/virulence_genes.csv'  % self.species ,'r') as infile:
            reader = csv.reader(infile)
            row = reader.next()
        return row

    def render_antibiotics_src(self):
        try:
            return self.render_str('src/predictor/%s/antibiotics.c' % self.species)        
        except TemplateNotFound:
            logging.warning("Using default antibiotics.c template")
            return self.render_str('src/predictor/common/antibiotics.c')        

    def render_cmd_line_include(self):
        try:
            return self.render_str('include/predictor/%s/cmd_line.h' % self.species)        
        except TemplateNotFound:
            logging.warning("Using default cmd_line.h template")
            return self.render_str('include/predictor/common/cmd_line.h')

    def render_cmd_line_src(self):
        try:
            return self.render_str('src/predictor/%s/cmd_line.c' % self.species)        
        except TemplateNotFound:
            logging.warning("Using default cmd_line.c template")
            return self.render_str('src/predictor/common/cmd_line.c')

    def render_antibiotics_include(self):
        try:
            return self.render_str('include/predictor/%s/antibiotics.h' % self.species)
        except TemplateNotFound:
            logging.warning("Using default antibiotics.h template")    
            return self.render_str('include/predictor/common/antibiotics.h')        
                    

    def write_to_file(self, st, filepath):
        dir_out = os.path.join(self.render_dir,os.path.dirname(filepath))
        if not os.path.exists(dir_out):
            os.makedirs(dir_out)
        with open(os.path.join(self.render_dir, filepath) ,'w') as outfile:
            outfile.write(st)                

    def render_and_write_antibiotics(self):
        st = self.render_antibiotics_src()
        self.write_to_file(st, 'src/predictor/%s/antibiotics.c' % self.species)
        st = self.render_antibiotics_include()
        self.write_to_file(st, 'include/predictor/%s/antibiotics.h' % self.species)      

    def render_and_write_cmd_line(self):
        st = self.render_cmd_line_src()
        self.write_to_file(st, 'src/predictor/%s/cmd_line.c' % self.species)
        st = self.render_cmd_line_include()
        self.write_to_file(st, 'include/predictor/%s/cmd_line.h' % self.species)

    def render_known_mutations_src(self):
        return self.render_str('src/predictor/core/known_mutations.c')

    def render_known_mutations_include(self):
        return self.render_str('include/predictor/core/known_mutations.h')        

    def render_and_write_known_mutations(self):
        st = self.render_known_mutations_src()
        self.write_to_file(st, 'src/predictor/core/known_mutations.c')            
        st = self.render_known_mutations_include()
        self.write_to_file(st, 'include/predictor/core/known_mutations.h')

    def render_gene_presence_src(self):
        return self.render_str('src/predictor/core/gene_presence.c')

    def render_gene_presence_include(self):
        return self.render_str('include/predictor/core/gene_presence.h')        

    def render_and_write_gene_presence(self):
        st = self.render_gene_presence_src()
        self.write_to_file(st, 'src/predictor/core/gene_presence.c')            
        st = self.render_gene_presence_include()
        self.write_to_file(st, 'include/predictor/core/gene_presence.h')            

    def render_main_src(self):
        try:
            return self.render_str('src/predictor/%s/main.c' % self.species )
        except TemplateNotFound:
            logging.warning("Using default main.c template")
            return self.render_str('src/predictor/common/main.c' )


    def render_and_write_main(self):
        st = self.render_main_src()
        self.write_to_file(st, 'src/predictor/%s/main.c' % self.species)            

    def render_and_write_all(self):
        self.render_and_write_antibiotics() 
        self.render_and_write_cmd_line() 
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

    def __init__(self, name):
        self.name = name
        self.gene_enum_to_drug_name = self.load_gene_enum_to_drug_name()
        self.mutation_fasta = MutationFasta('../data/%s/antibiotics/%s' % (self.species, self.name.lower() )  )

    def __str__(self):
        return self.name   

    @property
    def genes_resistance_induced_by(self):
        genes = []
        for gene,drugs in self.gene_enum_to_drug_name.iteritems():
            if self.name in drugs:
                genes.append(Gene(gene, self.species))
        return genes             

    @property 
    def has_epistatic_muts(self):
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


class GeneCodeGenerator(CodeGenerator):

    def __init__(self,name):
        self.name = name 
        self.gene_enum_to_drug_name = self.load_gene_enum_to_drug_name()

    def __str__(self):
        return self.name   

class Gene(GeneCodeGenerator):
    
    def __init__(self,name, species):
        self.species = species
        super(Gene, self).__init__(name=name) 


class SpeciesCodeGenerator(CodeGenerator):

    def __init__(self, species):
        self.species = species
        super(SpeciesCodeGenerator, self).__init__()
        self.taxon_coverage_threshold_dict = self._load_taxon_coverage_thresholds()

    def render_species_src(self):
        try:
            return self.render_str('src/predictor/%s/phylo/species.c' % self.species )
        except TemplateNotFound:
            logging.warning("Using default species.c template")
            return self.render_str('src/predictor/common/phylo/species.c' )

    def render_species_include(self):
        try:
            return self.render_str('include/predictor/%s/phylo/species.h' % self.species )
        except TemplateNotFound:
            logging.warning("Using default species.h template")
            return self.render_str('include/predictor/common/phylo/species.h' )
                                                
    def render_and_write_phylo(self):
        st = self.render_species_src()
        with open(os.path.join(self.render_dir,'src/predictor/%s/species.c' % self.species),'w') as outfile:
            outfile.write(st)
        st = self.render_species_include()
        with open(os.path.join(self.render_dir,'include/predictor/%s/species.h' % self.species),'w') as outfile:
            outfile.write(st) 

    @property
    def phylo_groups(self):
        phylo_groups = glob.glob('../data/%s/phylo/*'  %  (self.species))
        return [PhyloGroup(self.species, os.path.basename(f), self.taxon_coverage_threshold_dict) for f in phylo_groups]

    def _load_taxon_coverage_thresholds(self):
        try:
            with open('data/%s/taxon_coverage_threshold.json'  % self.species ,'r') as infile:
                d =  json.load(infile) 
        except IOError:
            d = {}
        return d  

class PhyloGroup(object):

    def __init__(self, species,  name, taxon_coverage_threshold_dict):
        self.name = name
        self.species = species
        self.enum = name.title()
        self.taxon_coverage_threshold_dict = taxon_coverage_threshold_dict
        self._load_basename_to_taxon()

    @property 
    def taxons(self):
        fasta_list = glob.glob('../data/%s/phylo/%s/*fasta'  %  (self.species, self.name))
        fasta_list.extend(glob.glob('../data/%s/phylo/%s/*.fa'  %  (self.species, self.name)))  
        return [Taxon(os.path.basename(f), self.name, self.taxon_coverage_threshold_dict, basename_to_taxon_name = self.basename_to_taxon_name) for f in fasta_list]

    def _load_basename_to_taxon(self):
        try:
            with open('data/%s/basename_to_taxon_name.json'  % self.species ,'r') as infile:
                self.basename_to_taxon_name =  json.load(infile)
        except:
            logging.warning("Cannot find data/%s/basename_to_taxon_name.json " % self.species)
            self.basename_to_taxon_name = {}        

class Taxon(object):

    def __init__(self, filename, phylo_group, taxon_coverage_threshold_dict = {}, basename_to_taxon_name = {}):
        self.filename = filename
        self.phylo_group = phylo_group
        self.basename_to_taxon_name = basename_to_taxon_name
        self.basename = make_safe_string(filename.split('.')[0])
        self.enum = self.basename.title() 
        self.threshold = taxon_coverage_threshold_dict.get(self.enum, self.default_threshold)
        

    @property 
    def name(self):
        return self.basename_to_taxon_name.get(self.basename, self.basename)

    @property
    def default_threshold(self):
        if self.phylo_group == "lineage":
            return 90
        elif self.phylo_group == "species":
            return 30
        else:
            return 50
    
