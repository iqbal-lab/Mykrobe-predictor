## Read the kmer counts into a hash
from atlas.utils import check_args
from atlas.typing import CoverageParser
from atlas.typing import Genotyper
from atlas.pheno import TBPredictor
from atlas.pheno import StaphPredictor
from atlas.pheno import GramNegPredictor
from atlas.metagenomics import AMRSpeciesPredictor
from pprint import pprint

STAPH_PANELS = ["Coagneg",
                "Staphaureus",
                "Saureus",
                "Sepidermidis", 
                "Shaemolyticus",
                "Sother",
                "staph_amr_genes",
                "staph-amr-mutations"]
GN_PANELS = ["gn-amr-genes","Escherichia_coli", "Klebsiella_pneumoniae","gn-amr-genes-extended"]
TB_PANELS = ["MTBC", "NTM",
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
            "xenopi"]

def run(parser, args):
    args = parser.parse_args()
    check_args(args)  
    if not args.species:
        panels = TB_PANELS + GN_PANELS + STAPH_PANELS
        panel_name = "tb-gn-staph-amr"
    elif args.species == "staph":
        panels = STAPH_PANELS
        panel_name = "staph-amr"
    elif args.species == "tb":
        panels = TB_PANELS
        panel_name = "tb-amr"
    elif args.species == "gn":
        panels = GN_PANELS
        panel_name = "gn-amr"                
    ## Run Cortex
    cp = CoverageParser(args, panels = panels, 
                                       verbose = False,
                                       panel_name = panel_name)
    cp.run()
    
    # Detect species
    species_predictor = AMRSpeciesPredictor(phylo_group_covgs = cp.covgs.get("phylo_group", {}),
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
    elif species_predictor.is_gram_neg_present():
        Predictor = GramNegPredictor
        try:
            depths = [species_predictor.out_json["phylogenetics"]["species"]["Klebsiella_pneumoniae"]["median_depth"]]
        except KeyError:
            depths = [species_predictor.out_json["phylogenetics"]["species"]["Escherichia_coli"]["median_depth"]]

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
