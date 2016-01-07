import os

class Panel(object):

	def __init__(self, name):
		self.name = name

	def __str__(self):
		return self.name

	def __repr__(self):
		return self.name

	@property
	def filepath(self):
		return os.path.abspath("data/panels/%s.fasta" % self.name)  