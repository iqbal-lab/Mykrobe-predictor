## Read the kmer counts into a hash
from atlas.utils import check_args
from atlas.typing import CoverageParser
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


    ## Run Cortex
    cp = CoverageParser(args, panels = ["Coagneg", "Staphaureus", "Saureus", "Sepidermidis", 
                                   "Shaemolyticus", "Sother","staph_amr_genes",
                                   "staph-amr-mutations"], 
                                   verbose = False)
    cp.run()
    
    # Detect species
    species_predictor = SpeciesPredictor(phylo_group_covgs = cp.covgs["phylo_group"],
                    species_covgs = cp.covgs["species"],
                    lineage_covgs = cp.covgs.get("lineage", {}),
                    base_json = cp.out_json[args.sample])
    species_predictor.run()

    ## Genotype
    q = args.quiet
    args.quiet = True    
    gt = Genotyper(args, depths = [species_predictor.out_json["phylogenetics"]["species"]["Saureus"]["median_depth"]],
                variant_covgs = cp.covgs["variant"],
                gene_presence_covgs = cp.covgs["presence"],
                verbose = False)
    gt.run()
    args.quiet = q

    ## AMR prediction
    staph_predictor = StaphPredictor(typed_variants = gt.variant_covgs,
                   called_genes = gt.gene_presence_covgs,
                   base_json = gt.out_json[args.sample])
    staph_predictor.run()    
