class Typer(object):

	def __init__(self, depths, contamination_depths = [], error_rate = 0.05):
		self.depths = depths
		self.contamination_depths = contamination_depths
		self.error_rate = error_rate

	def type(self, l):
		raise NotImplemented("Implemented in sub class")
