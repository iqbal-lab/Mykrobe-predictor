from atlas.typing.typer.base import Typer


class PresenceTyper(Typer):

    "Initiated with expected depths and contamination depths"
    
    def __init__(self, depths, contamination_depths = []):
        super(PresenceTyper, self).__init__(depths, contamination_depths)

    def genotype(self, sequence_coverage):
        "Takes a single SequenceCoverage object (or child) and returns genotype"
        if sequence_coverage.percent_coverage > sequence_coverage.percent_coverage_threshold:
            sequence_coverage.set_genotype("1/1")
        else:
            sequence_coverage.set_genotype("0/0")
        return sequence_coverage

class GeneCollectionTyper(Typer):

    """Types a collection of genes returning only the most likely version 
        in the collection"""

    def __init__(self, depths, contamination_depths = []):
        super(GeneCollectionTyper, self).__init__(depths, contamination_depths)
        self.presence_typer = PresenceTyper(depths, contamination_depths)

    def genotype(self, sequence_coverage_collection):
        """Types a collection of genes returning the most likely gene version 
            in the collection with it's genotype"""
        best_version = self.get_best_version(sequence_coverage_collection.values())
        return self.presence_typer.genotype(best_version)

    def get_best_version(self, sequence_coverages):
        sequences.sort(key=lambda x: x.percent_coverage, reverse=True)
        current_best_gene = sequences[0]
        for gene in sequences[1:]:
            if gene.percent_coverage < current_best_gene.percent_coverage:
                return current_best_gene
            else:
                if gene.median_depth > current_best_gene.median_depth:
                    current_best_gene = gene
        return current_best_gene

