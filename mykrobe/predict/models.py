## Defines Mykrobe predictor output
from mongoengine import Document
from mongoengine import DictField

from mykrobe.utils import unique
class MykrobePredictorSusceptibilityResult(Document):

    susceptibility = DictField()

    @classmethod
    def create(cls, susceptibility):
        return cls(susceptibility = susceptibility)

    def to_dict(self):
    	return self.to_mongo().to_dict()

    def __eq__(self,other):
        return self.susceptibility == other.susceptibility

    def diff(self, other):
        diff = {}
        ## Compares the antibiogram of two predictor results
        drugs = unique(self.drugs + other.drugs)
        for drug in drugs:
            predict1,predict2 = self.susceptibility.get(drug, {"predict": "NA"}).get("predict"), other.susceptibility.get(drug,{"predict": "NA"}).get("predict")
            if predict1 != predict2:
                diff[drug] = {}
                diff[drug]["predict"] = (predict1,predict2)
        return diff

    @property
    def drugs(self):
        return self.susceptibility.keys()



    # def compare(self, truth):
    #     """
    #     truth = {
    #             "Rifampicin": "R",
    #             "Ethambutol":"S"
    #             ...
    #             }
    #     """
    #     self.comparision = {}
    #     self.truth = truth
    #     for drug in self.drug_list:
    #         try:
    #             predict = self.called_drugs.get(drug, None).call
    #         except AttributeError:
    #             raise AttributeError("No susceptibility predictions %s" % drug)
    #         truth = self.truth.get(drug, None)
    #         indiv_compare = self.compare_calls(predict, truth)
    #         self.comparision[drug] = {
    #             'predict': predict,
    #             'truth': truth,
    #             'comparision': indiv_compare
    #         }
    #     return self.comparision

    # def compare_calls(self, predict, truth):
    #     return compare_calls(predict, truth)

    # @property
    # def called_variants(self):
    #     return [self.CalledVariant(k, v) for k, v in self.results_json.get(
    #         'called_variants').iteritems()]

    # class CalledVariant(object):

    #     def __init__(self, var, called_variants_dict):
    #         """
    #         "katG_S315X" :{
    #             "R_per_cov": "100",
    #             "S_per_cov": "100",
    #             "R_median_cov": "51",
    #             "S_median_cov": "1",
    #             "conf": "124",
    #             "induced_resistance": "Isoniazid"
    #         }"""
    #         self.var = var
    #         self.R_per_cov = int(called_variants_dict.get("R_per_cov"))
    #         self.S_per_cov = int(called_variants_dict.get("S_per_cov"))
    #         self.R_median_cov = int(called_variants_dict.get("R_median_cov"))
    #         self.S_median_cov = int(called_variants_dict.get("S_median_cov"))
    #         self.conf = int(called_variants_dict.get("conf"))
    #         self.induced_resistance = called_variants_dict.get(
    #             "induced_resistance")

    # @property
    # def called_genes(self):
    #     try:
    #         return [self.CalledGene(k, v) for k, v in self.results_json.get(
    #             'called_genes').iteritems()]
    #     except AttributeError:
    #         return []

    # class CalledGene(object):

    #     def __init__(self, gene, called_genes_dict):
    #         """
    #         "blaZ" :{
    #                     "per_cov": "85",
    #                     "median_cov": "11",
    #                     "conf": "6054",
    #             "induced_resistance": "Penicillin"
    #             }"""
    #         self.gene = gene
    #         self.median_cov = int(called_genes_dict.get("median_cov"))
    #         self.conf = int(called_genes_dict.get("conf"))
    #         self.induced_resistance = called_genes_dict.get(
    #             "induced_resistance")

    # class CalledDrug(object):

    #     def __init__(self, mykrobe_result, drug, call):
    #         self.drug = drug
    #         if self.drug in ['CIP', "MOX", "OFX"]:
    #             self.mykrobe_drug = "Quinolones"
    #         else:
    #             self.mykrobe_drug = self.drug
    #         self.call = call[0]
    #         self.mykrobe_result = mykrobe_result

    #     def __str__(self):
    #         return "%s:%s" % (self.drug, self.call)

    #     @property
    #     def called_variants(self):
    #         """Return the called variants objects associated that have induced this drug"""
    #         return [
    #             cv for cv in self.mykrobe_result.called_variants if self.mykrobe_drug in cv.induced_resistance]

    #     @property
    #     def called_genes(self):
    #         """Return the called gene objects associated that have induced this drug"""
    #         return [
    #             gene for gene in self.mykrobe_result.called_genes if self.mykrobe_drug in gene.induced_resistance]

    #     @property
    #     def is_induced_by_gene(self):
    #         return len(self.called_genes) >= 1

    #     @property
    #     def is_induced_by_variant(self):
    #         return len(self.called_variants) >= 1

    #     @property
    #     def called_gene_with_max_conf(self):
    #         try:
    #             return sorted(
    #                 self.called_genes, key=lambda x: x.conf, reverse=True)[0]
    #         except IndexError:
    #             return None

    #     @property
    #     def called_variant_with_max_conf(self):
    #         try:
    #             return sorted(
    #                 self.called_variants, key=lambda x: x.conf, reverse=True)[0]
    #         except IndexError:
    #             return None

    #     @property
    #     def max_conf(self):
    #         if self.is_induced_by_variant:
    #             try:
    #                 return self.called_variant_with_max_conf.conf
    #             except AttributeError:
    #                 return None
    #         elif self.is_induced_by_gene:
    #             try:
    #                 return self.called_gene_with_max_conf.conf
    #             except AttributeError:
    #                 return None

    #     @property
    #     def called_variant_names(self):
    #         return ",".join(cv.var for cv in self.called_variants)

    #     @property
    #     def comparision(self):
    #         try:
    #             return self.mykrobe_result.comparision.get(
    #                 self.drug).get('comparision')
    #         except AttributeError:
    #             raise AttributeError(
    #                 "You must run compare(truth) before accessing drug comparisions")

    #     @property
    #     def R_median_cov(self):
    #         if self.is_induced_by_variant:
    #             try:
    #                 return self.called_variant_with_max_conf.R_median_cov
    #             except AttributeError:
    #                 raise AttributeError(
    #                     "S/N calls don't have coverage information")
    #         elif self.is_induced_by_gene:
    #             try:
    #                 return self.called_gene_with_max_conf.median_cov
    #             except AttributeError:
    #                 raise AttributeError(
    #                     "S/N calls don't have coverage information")

    #     @property
    #     def S_median_cov(self):
    #         try:
    #             return self.called_variant_with_max_conf.S_median_cov
    #         except AttributeError:
    #             raise AttributeError(
    #                 "S/N calls don't have coverage information")

    #     @property
    #     def R_per_cov(self):
    #         try:
    #             return self.called_variant_with_max_conf.R_per_cov
    #         except AttributeError:
    #             raise AttributeError(
    #                 "S/N calls don't have coverage information")

    #     @property
    #     def S_per_cov(self):
    #         try:
    #             return self.called_variant_with_max_conf.S_per_cov
    #         except AttributeError:
    #             raise AttributeError(
    #                 "S/N calls don't have coverage information")

    # @property
    # def phylo_group(self):
    #     return self.PhyloGroup(self)

    # @property
    # def lineage(self):
    #     return self.Lineage(self)

    # @property
    # def major_phylo_group(self):
    #     return self.phylo_group.major_name

    # @property
    # def major_species(self):
    #     return self.species.major_name

    # @property
    # def major_lineage(self):
    #     return self.lineage.major_name 

    # @property
    # def is_mixed(self):
    #     return len(self.phylo_group.names) > 1                       

    # class Phylo(object):

    #     def __init__(self, mykrobe_result):
    #         self.mykrobe_result = mykrobe_result

    #     def __str__(self):
    #         return ",".join(self.names)

    #     @property
    #     def dict(self):
    #         return self.mykrobe_result.results_json.get(
    #             'phylogenetics').get(self.group)

    #     @property
    #     def names(self):
    #         return self.dict.keys()

    #     @property
    #     def median_covg(self):
    #         return [int(i) for i in self.dict.values()]

    #     @property
    #     def major_index(self):
    #         return self.median_covg.index(max(self.median_covg))

    #     @property
    #     def major_name(self):
    #         return self.names[self.major_index]

    # class PhyloGroup(Phylo):

    #     def __init__(self, mykrobe_result):
    #         super(MykrobeResult.PhyloGroup, self).__init__(mykrobe_result)
    #         self.group = "phylo_group"

    # class Species(Phylo):

    #     def __init__(self, mykrobe_result):
    #         super(MykrobeResult.Species, self).__init__(mykrobe_result)
    #         self.group = "species"

    # class Lineage(Phylo):

    #     def __init__(self, mykrobe_result):
    #         super(MykrobeResult.Lineage, self).__init__(mykrobe_result)
    #         self.group = "lineage"
