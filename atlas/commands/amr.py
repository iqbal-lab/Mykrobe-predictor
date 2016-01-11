## Read the kmer counts into a hash
from atlas.utils import check_args
from atlas.typing import Genotyper
from atlas.pheno import TBPredictor
from atlas.pheno import StaphPredictor
from atlas.metagenomics import SpeciesPredictor
from pprint import pprint
def run(parser, args):
    args = parser.parse_args()
    check_args(args)  
    

    # # Genotype
    # gt = Genotyper(args, panel = "panel-tb-21")
    # gt.run()
    # # 
    # TBPredictor(typed_variants = gt.variant_covgs,
    # 			called_genes = gt.gene_presence_covgs,
    # 			sample = args.sample).run()

    q = args.quiet
    args.quiet = True
    ## Run Cortex
    gt = Genotyper(args, panels = ["Coagneg", "Staphaureus", "Saureus", "Sepidermidis", 
                                   "Shaemolyticus", "Sother","staph_amr_genes",
                                   "staph-amr-mutations"], 
                                   depths = [100],
                                   verbose = False)
    gt.run()
    args.quiet = q
    # Detect species
    species_predictor = SpeciesPredictor(phylo_group_covgs = gt.covgs["phylo_group"],
                    species_covgs = gt.covgs["species"],
                    lineage_covgs = gt.covgs.get("lineage", {}),
                    base_json = gt.out_json[args.sample])
    species_predictor.run()
    StaphPredictor(typed_variants = gt.variant_covgs,
                   called_genes = gt.gene_presence_covgs,
                   base_json = gt.out_json[args.sample]).run()    
