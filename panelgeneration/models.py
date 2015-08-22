from Bio import SeqIO
from Bio.Seq import Seq
from copy import copy 
import itertools
from collections import Counter

def unique(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]

class Variant(object):

	def __init__(self, ref, pos, alt):
		self.ref = ref 
		self.pos = pos
		self.alt = alt

	def __str__(self):
		return "".join([self.ref, str(self.pos), self.alt])

	def __repr__(self):
		return "".join([self.ref, str(self.pos), self.alt])

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

	def create(self, v, context = []):
		## Position should be 1 based
		assert v.pos > 0
		i, start_index, end_index = self._get_start_end(v.pos)
		self._check_valid_variant(v)
		wild_type_reference = copy(self.ref[start_index:end_index])
		backgrounds = self._generate_all_backgrounds_using_context(v, i, wild_type_reference, context)
		alternates = self._generate_alternates_on_all_backgrounds(i, v,  backgrounds)
		return Panel(wild_type_reference, v.pos, alternates)

	def _generate_alternates_on_all_backgrounds(self, i, v, backgrounds):
		alternates = []
		for background in backgrounds:
			alternate = copy(background)
			assert alternate[i] == v.ref	
			alternate[i] = v.alt
			alternates.append(alternate)
		return alternates
		
	def _check_valid_variant(self, v):
		index = v.pos - 1
		if self.ref[index] != v.ref:
			raise ValueError("""Cannot create alleles as ref at pos %i is not %s (it's %s) are you sure you're using one-based co-ordinates?
							 """ % (v.pos, v.ref, self.ref[index]))			

	def _generate_all_backgrounds_using_context(self, v, i, wild_type_reference, context = []):
		backgrounds = [wild_type_reference]
		if context:
			contexts = self._create_multiple_contexts(context)
			for context in contexts:
				context_combinations = self._get_combinations_of_backgrounds(context)
				for variants in context_combinations:
					new_background = copy(wild_type_reference)
					for variant in variants:
						j = i + variant.pos - v.pos
						assert new_background[j] == variant.ref						
						new_background[j] = variant.alt
					backgrounds.append(new_background)
		return backgrounds

	def _create_multiple_contexts(self, context):
		new_contexts = self._recursive_context_creator([context])
		return new_contexts

	def _recursive_context_creator(self, contexts):
		valid = True
		for i,context in enumerate(contexts):
			position_counts = Counter([v.pos for v in context])
			if not all([c == 1 for c in position_counts.values()]):
				valid = False
		if valid:
			return contexts
		else:
			for i,context in enumerate(contexts):
				position_counts = Counter([v.pos for v in context])
				if not all([c == 1 for c in position_counts.values()]):
					position = position_counts.most_common(1)[0][0]
					repeated_vars = []
					common_vars = []
					for var in context:
						if var.pos == position:
							repeated_vars.append(var)
						else:
							common_vars.append(var)
					new_contexts = []
					for var in repeated_vars:
						context = [var] + common_vars
						new_contexts.append(context)
					contexts.pop(i)
					contexts += new_contexts
					self._recursive_context_creator(contexts)
		return contexts

	def _get_combinations_of_backgrounds(self, context):
		combination_context = []
		for L in range(1, len(context)+1):
		  for subset in itertools.combinations(context, L):
		    combination_context.append(subset)	
		return combination_context

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

	def __init__(self, ref, pos, alts):
		self.ref = "".join(ref)
		self.pos = pos
		self.alts = unique(["".join(alt) for alt in alts])

		