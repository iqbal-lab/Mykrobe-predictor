from Bio import SeqIO
from os import path

class MutationFasta(object):

	def __init__(self,filepath):
		self.filepath = filepath
		self.mut_list = []
		self.mut_gene_list = []
		self.read_mut_file()

	def read_mut_file(self):
		filepath = "%s.fa" % self.filepath
		if path.exists(filepath):
			handle = open(filepath, "r")
			for record in SeqIO.parse(handle, "fasta") :
				mut = self.mut_enum(record.id)
				if not mut in self.mut_list:
					self.mut_list.append(mut)
				gene = self.gene_enum(record.id)
				if not gene in self.mut_gene_list:
					self.mut_gene_list.append(gene)

	def mut_enum(self,record_id):
		mut_code = record_id.split('_')[1]
		gene = self.gene_enum(record_id)
		return "_".join([gene,mut_code])

	def gene_enum(self,record_id):
		return record_id.split('-')[-1]

	@property
	def num_mutations(self):
		return len(self.mut_list)

	@property 
	def first_mut(self):
		return self.mut_list[0]

	@property 
	def last_mut(self):
		return self.mut_list[-1]	

