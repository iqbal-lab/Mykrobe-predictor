from Bio import SeqIO
from Bio.Seq import Seq
from copy import copy 
class AlleleGenerator(object):
	"""docstring for PanelGenerator"""
	def __init__(self, reference_filepath, kmer = 31):
		self.reference_filepath = reference_filepath
		self.kmer = kmer
		self.ref = []
		self.ref_length = 0
		self._read_reference()

	def _read_reference(self):
		for record in SeqIO.parse(self.reference_filepath, 'fasta'):
			self.ref += list(record.seq)
		self.ref_length = len(self.ref)

	def create(self, ref, pos, alt):
		## Position should be 1 based
		if pos < 1:
			pos = self.ref_length + pos
		assert pos > 0
		index = pos - 1
		if self.ref[index] != ref:
			raise ValueError("Cannot create alleles as ref at pos %i is not %s (it's %s) are you sure you're using one-based co-ordinates?" % (pos, ref, self.ref[index]))
		i, start_index, end_index = self._get_start_end(pos)
		reference = self.ref[start_index:end_index]
		alternate = copy(reference)
		alternate[i] = alt
		return Panel(reference, pos, alternate)	


	def _get_start_end(self, pos):
		start_index = pos - self.kmer
		end_index =  pos + self.kmer +1
		i = self.kmer - 1
		if start_index < 0:
			diff = abs(start_index)
			start_index = 0 
			end_index += diff
			i -= diff
		elif end_index > self.ref_length:
			diff = abs(end_index - self.ref_length)
			end_index = self.ref_length
			start_index -= diff
			i += diff
		return (i, start_index, end_index)


class Panel(object):

	def __init__(self, ref, pos, alt):
		self.ref = "".join(ref)
		self.pos = pos
		self.alt = "".join(alt)

		