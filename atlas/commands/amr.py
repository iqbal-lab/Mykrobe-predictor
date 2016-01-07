## Read the kmer counts into a hash
from atlas.utils import check_args
from atlas.typing import Genotyper
from atlas.pheno import TBPredictor
from atlas.pheno import StaphPredictor

def run(parser, args):
    args = parser.parse_args()
    check_args(args)  
    # Detect species

    # # Genotype
    # gt = Genotyper(args, panel = "panel-tb-21")
    # gt.run()
    # # 
    # TBPredictor(typed_variants = gt.variant_covgs,
    # 			called_genes = gt.gene_presence_covgs,
    # 			sample = args.sample).run()


    # Genotype
    q = args.quiet
    args.quiet = True
    gt = Genotyper(args, panels = ["staph_amr_genes", "staph-amr-mutations"])
    gt.run()
    args.quiet = q
    # 
    StaphPredictor(typed_variants = gt.variant_covgs,
                   called_genes = gt.gene_presence_covgs,
                   sample = args.sample, 
                   base_json = gt.out_json).run()    
