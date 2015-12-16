## Read the kmer counts into a hash
from atlas.utils import check_args
from atlas.typing import Genotyper
from atlas.pheno import TBPredictor

def run(parser, args):
    args = parser.parse_args()
    check_args(args)  
    # Detect species

    # Genotype
    gt = Genotyper(args, panel = "panel-tb-21")
    gt.run()
    # 
    TBPredictor(typed_variants = gt.variant_covgs,
    			called_genes = gt.gene_presence_covgs,
    			sample = args.sample).run()
