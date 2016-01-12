## Read the kmer counts into a hash
from atlas.utils import check_args
from atlas.typing import CoverageParser
from atlas.typing import Genotyper
from atlas.pheno import TBPredictor
from atlas.pheno import StaphPredictor
from atlas.metagenomics import AMRSpeciesPredictor
from pprint import pprint
def run(parser, args):
    args = parser.parse_args()
    check_args(args)  
    
    ## Run Cortex
    cp = CoverageParser(args, panels = ["MTBC", "NTM", "Coagneg",
                                        "Staphaureus",
                                        "Saureus", "Sepidermidis", 
                                       "Shaemolyticus",
                                        "Sother","staph_amr_genes",
                                       "staph-amr-mutations",
                                       "panel-tb-21",
                                       "abscessus",
                                        "africanum",
                                        "aromaticivorans",
                                        "avium",
                                        "bovis",
                                        "branderi",
                                        "caprae",
                                        "chelonae",
                                        "chlorophenolicum",
                                        "chubuense",
                                        "colombiense",
                                        "crocinum",
                                        "flavescens",
                                        "fluoranthenivorans",
                                        "fortuitum",
                                        "gilvum",
                                        "gordonae",
                                        "hodleri",
                                        "interjectum",
                                        "intracellulare",
                                        "kansasii",
                                        "lentiflavum",
                                        "leprae",
                                        "malmoense",
                                        "marinum",
                                        "mucogenicum",
                                        "pallens",
                                        "peregrinum",
                                        "phage",
                                        "pyrenivorans",
                                        "rhodesiae",
                                        "rufum",
                                        "rutilum",
                                        "scrofulaceum",
                                        "senegalense",
                                        "smegmatis",
                                        "sphagni",
                                        "szulgai",
                                        "triplex",
                                        "tuberculosis",
                                        "tusciae",
                                        "ulcerans",
                                        "vaccae",
                                        "xenopi"], 
                                       verbose = False,
                                       panel_name = "amr")
    cp.run()
    
    # Detect species
    species_predictor = AMRSpeciesPredictor(phylo_group_covgs = cp.covgs["phylo_group"],
                    species_covgs = cp.covgs["species"],
                    lineage_covgs = cp.covgs.get("lineage", {}),
                    base_json = cp.out_json[args.sample])
    species_predictor.run()



    # ## AMR prediction
    if species_predictor.is_saureus_present():
        depths = [species_predictor.out_json["phylogenetics"]["phylo_group"]["Staphaureus"]["median_depth"]]
        Predictor = StaphPredictor
    elif species_predictor.is_mtbc_present():
        depths = [species_predictor.out_json["phylogenetics"]["phylo_group"]["MTBC"]["median_depth"]]
        Predictor = TBPredictor
    ## Genotype
    q = args.quiet
    args.quiet = True    
    gt = Genotyper(args, depths = depths,
                variant_covgs = cp.covgs["variant"],
                gene_presence_covgs = cp.covgs["presence"],
                verbose = False, 
                base_json = cp.out_json,
                contamination_depths = species_predictor.contamination_depths())
    gt.run()
    args.quiet = q
    predictor = Predictor(typed_variants = gt.variant_covgs,
                   called_genes = gt.gene_presence_covgs,
                   base_json = gt.out_json[args.sample])
    predictor.run()    
