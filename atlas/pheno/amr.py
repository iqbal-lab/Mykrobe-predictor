import os
import json
from atlas.utils import unique
from atlas.utils import flatten

from atlas.typing import SequenceCoverage
from atlas.typing import TypedVariant

from pprint import pprint


DEFAULT_MIN_GENE_CN = 0.03
DEFAULT_MIN_VARIANT_CN = 0.1


class BasePredictor(object):

    def __init__(self, typed_variants, called_genes, base_json={}):
        self.typed_variants = typed_variants
        self.called_genes = called_genes
        self.drugs = self._get_drug_list_from_variant_to_resistance_drug()
        self.resistance_prediction = self._create_initial_resistance_prediction()
        self.out_json = base_json
        self._coveage_threshold = {
            "ermA": 0.19,
            "ermB": 0.19,
            "ermC": 0.19,
            "ermT": 0.19,
            "ermY": 0.19,
            "fusA": 0.03,
            "aacAaphD": 0.04,
            "mecA": 0.06,
            "mupA": 0.21,
            "blaZ": 0.04,
            "tetK": 0.13
        }

    def _create_initial_resistance_prediction(self):
        self.resistance_predictions = dict((k, "I") for k in self.drugs)

    def _get_drug_list_from_variant_to_resistance_drug(self):
        return unique(
            flatten(
                self.variant_or_gene_name_to_resistance_drug.values()))

    def predict_antibiogram(self):
        for name, variants in self.typed_variants.items():
            for variant in variants:
                self._update_resistance_prediction(variant)
        for name, gene in self.called_genes.items():
            self._update_resistance_prediction(gene)

    def _update_resistance_prediction(self, variant_or_gene):
        name = self._get_name(variant_or_gene)
        drugs = self._get_drugs(name)

        resistance_prediction = self._resistance_prediction(variant_or_gene)
        for drug in drugs:
            variant_or_gene.add_induced_resistance(drug)
            current_resistance_prediction = self.resistance_predictions.get(
                drug)
            assert resistance_prediction is not None
            if current_resistance_prediction in ["I", "N"]:
                self.resistance_predictions[drug] = resistance_prediction
            elif current_resistance_prediction == "S":
                if resistance_prediction in ["r", "R"]:
                    self.resistance_predictions[drug] = resistance_prediction
            elif current_resistance_prediction == "r":
                if resistance_prediction == "R":
                    self.resistance_predictions[drug] = resistance_prediction

    def _get_name(self, variant_or_gene):
        if variant_or_gene.alt_name:
            name = variant_or_gene.alt_name
        else:
            name = variant_or_gene.name
        return name

    def _get_drugs(self, name, lower=False):
        if lower:
            name = name.lower()
        try:
            drugs = self.variant_or_gene_name_to_resistance_drug[name]
        except KeyError:
            try:
                drugs = self.variant_or_gene_name_to_resistance_drug[
                    name.split("-")[0]]
            except KeyError:
                talt_name = list(name)
                talt_name[-1] = "X"
                try:
                    drugs = self.variant_or_gene_name_to_resistance_drug[
                        "".join(talt_name)]
                except KeyError:
                    drugs = []
                    if not lower:
                        return self._get_drugs(name, lower=True)
                    else:
                        pass
                        # print ("Warning:NoEntry for %s" % name)

        return drugs

    def _resistance_prediction(self, variant_or_gene):
        if variant_or_gene.gt == "1/1":
            if self._coverage_greater_than_threshold(variant_or_gene):
                return "R"
            else:
                return "S"
        elif variant_or_gene.gt == "0/1":
            if self._coverage_greater_than_threshold(variant_or_gene):
                return "r"
            else:
                return "S"
        elif variant_or_gene.gt == "0/0":
            return "S"
        else:
            return "I"

    def _coverage_greater_than_threshold(self, variant_or_gene):
        if isinstance(variant_or_gene, SequenceCoverage):
            return variant_or_gene.copy_number > self._coveage_threshold.get(
                variant_or_gene.name,
                DEFAULT_MIN_GENE_CN)
        elif isinstance(variant_or_gene, TypedVariant):
            return variant_or_gene.copy_number > self._coveage_threshold.get(
                variant_or_gene.name,
                DEFAULT_MIN_VARIANT_CN)
        else:
            raise TypeError(
                "Must be either SequenceCoverage or TypedVariant object")

    def run(self):
        self.predict_antibiogram()
        self.out_json["susceptibility"] = self.resistance_predictions


def load_json(f):
    with open(f, 'r') as infile:
        return json.load(infile)


class TBPredictor(BasePredictor):

    def __init__(self, typed_variants, called_genes, base_json={}):
        self.data_dir = os.path.abspath(
            os.path.join(
                os.path.dirname(__file__),
                '../../data/predict/tb/'))
        self.variant_or_gene_name_to_resistance_drug = load_json(
            os.path.join(
                self.data_dir,
                "variant_to_resistance_drug.json"))
        super(
            TBPredictor,
            self).__init__(
            typed_variants,
            called_genes,
            base_json)


class StaphPredictor(BasePredictor):

    def __init__(self, typed_variants, called_genes, base_json={}):
        self.data_dir = os.path.abspath(
            os.path.join(
                os.path.dirname(__file__),
                '../../data/predict/staph/'))
        self.variant_or_gene_name_to_resistance_drug = load_json(
            os.path.join(
                self.data_dir,
                "variant_to_resistance_drug.json"))
        super(
            StaphPredictor,
            self).__init__(
            typed_variants,
            called_genes,
            base_json)


class GramNegPredictor(BasePredictor):

    def __init__(self, typed_variants, called_genes, base_json={}):
        self.data_dir = os.path.abspath(
            os.path.join(
                os.path.dirname(__file__),
                '../../data/predict/gn/'))
        self.variant_or_gene_name_to_resistance_drug = load_json(
            os.path.join(
                self.data_dir,
                "variant_to_resistance_drug.json"))
        super(
            GramNegPredictor,
            self).__init__(
            typed_variants,
            called_genes,
            base_json)
