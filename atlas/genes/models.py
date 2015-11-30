from Bio import SeqIO
from Bio.Data import CodonTable 
import itertools

from atlas.vcf2db import split_var_name

def flatten(l):
    return [item for sublist in l for item in sublist]

class Region(object):

    def __init__(self, reference, start, end, forward = True):
        self.reference = reference
        self.start = start
        self.end = end
        self.forward = forward

    @property
    def strand(self):
        if self.forward:
            return "forward"
        else:
            return "reverse"

    @property
    def seq(self):
        if self.forward:
            return self.reference[self.start - 1:self.end]
        else:
            return self.reference[self.start - 1:self.end].reverse_complement()

    def get_reference_position(self, pos):
        if pos < 0 and self.forward:
            return self.start + pos
        elif pos < 0 and not self.forward:
            ## Upstream of a gene on the reverse stand
            return self.end - pos            
        elif pos > 0:
            return self.start + pos - 1
        else:
            raise ValueError("Positions are 1-based")

class Gene(Region):

    def __init__(self, name, reference, start, end, forward = True):
        super(self.__class__, self).__init__(reference, start, end, forward)
        self.name = name
        self.translation_table = 11

    @property
    def prot(self):
        return self.seq.translate(table = self.translation_table).rstrip("*")

    def get_codon(self, pos):
        if pos > len(self.prot):
            raise ValueError("There are only %s aminoacids in this gene" % len(self.prot))
        else:
            return self.seq[(3 *(pos - 1)) : pos * 3]

    def __str__(self):
        return "Gene:%s" % self.name

    def __repr__(self):
        return "Gene:%s" % self.name        

class GeneAminoAcidChangeToDNAVariants():

    def __init__(self, reference, genbank):
        self.reference = self._parse_reference(reference)
        self.genbank = self._parse_genbank(genbank)
        self.backward_codon_table = self._make_backward_codon_table()

    def _parse_reference(self, reference):
        with open(reference, "r") as infile:
            return list(SeqIO.parse(infile, "fasta"))[0].seq

    def _parse_genbank(self, genbank):
        d = {}
        with open(genbank, 'r') as infile:
            for feat in SeqIO.read(infile, "genbank").features:
                if feat.type == "gene":
                    name = feat.qualifiers.get("gene",feat.qualifiers.get("db_xref"))[0]
                    ## SeqIO converts to 0-based
                    d[name] = Gene(name, self.reference, start = feat.location.start + 1, 
                        end = feat.location.end + 1, forward = bool(feat.location.strand))
        return d

    def _make_backward_codon_table(self):
        table = {}
        standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
        codons = self._generate_all_possible_codons()
        for codon in codons:
            if codon not in standard_table.stop_codons:
                try:
                    table[standard_table.forward_table[codon]].append(codon)
                except:
                    table[standard_table.forward_table[codon]] = [codon]
        return table

    def _generate_all_possible_codons(self):
        return ["".join(subset) for subset in itertools.product(["A", "T", "C", "G"], repeat = 3)]

    def get_alts(self, amino_acid):
        if amino_acid == "X":
            return flatten(self.backward_codon_table.values())
        else:
            return self.backward_codon_table[amino_acid]

    def get_location(self, gene, pos):
        dna_pos = (3 * (pos - 1)) + 1#, (pos * 3) + 1
        # return [gene.get_reference_position(p) for p in dna_pos]
        return gene.get_reference_position(dna_pos)

    def get_variant_names(self, gene, mutation):
        ref, start, alt = split_var_name(mutation)
        gene = self.get_gene(gene)
        assert gene.prot[start - 1] == ref
        ref_codon = gene.get_codon(start)
        alt_codons = self.get_alts(alt)
        if ref_codon in alt_codons: alt_codons.remove(ref_codon)
        location = self.get_location(gene, start)
        alternative = "/".join(alt_codons)
        names = ["".join(["".join(ref_codon), str(location), "".join(alt_codon)]) for alt_codon in alt_codons]
        return names

    def get_gene(self, gene):
        return self.genbank[gene]





