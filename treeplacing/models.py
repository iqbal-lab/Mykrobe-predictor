from vcf2db import Variant
from vcf2db import VariantSet
from vcf2db import GenotypedVariant
from vcf2db import CallSet

class Placer(object):

    """Placer"""
    def __init__(self, root):
        super(Placer, self).__init__()
        self.root = root

    def place(self, sample):
    	gvs = GenotypedVariant.objects(call_set = CallSet.objects.get(name = sample)).distinct('name')
    	return self.root.search(variants = gvs)

# class Tree(dict):
#     """Tree is defined by a dict of nodes"""
#     def __init__(self):
#         super(Tree, self).__init__()

class Node(object):
    """docstring for Node"""
    def __init__(self,  children = []):
        super(Node, self).__init__()
        self.parent = None         
        self.children = children ## List of nodes
        for child in self.children:
        	child.add_parent(self)

    def add_parent(self, parent):
    	self.parent = parent

    def other_child(self, node):
    	for child in self.children:
    		if child != node:
    			return child

    @property
    def samples(self):
        samples = []
        for child in self.children:
            samples.extend(child.samples)
        return samples ## List of sample below node in tree

    @property 
    def num_samples(self):
        return len(self.samples)

    @property
    def is_leaf(self):
        return False        

    @property
    def is_node(self):
        return True  

    @property
    def phylo_snps(self):
    	ingroup = VariantSet.objects(name__in = self.samples)
    	if self.parent:
    		outgroup = VariantSet.objects(name__in = self.parent.other_child(self).samples)
    	else:
    		outgroup = []
    	non_unique_variant_names = Variant.objects(variant_set__nin = ingroup , variant_set__in = outgroup).distinct('name')
    	phylo_snps = Variant.objects(variant_set__in = ingroup, name__nin = non_unique_variant_names ).distinct('name')
        return phylo_snps

    def search(self, variants):
    	assert self.children[0].parent is not None
    	assert self.children[1].parent is not None
    	overlap = []
    	overlap = (float(len(set(self.children[0].phylo_snps) & set(variants)) )/ len(set(self.children[0].phylo_snps) | set(variants) ),
    	           float(len(set(self.children[1].phylo_snps) & set(variants))) / len( set(self.children[1].phylo_snps) | set(variants)  ) )
    	print overlap, self.children[0], self.children[1]
    	# print overlap, len(variants), len(self.children[0].phylo_snps)
    	if overlap[0] > overlap[1]:
    		return self.children[0].search(variants)
    	elif overlap[1] > overlap[0]:
    		return self.children[1].search(variants)
    	else:
    		return self.samples

    def __repr__(self):
    	return "Node : %s " % ",".join(self.samples)

class Leaf(Node):

    def __init__(self, sample):
        super(Leaf, self).__init__()
        self.sample = sample 

    @property 
    def samples(self):
        return [self.sample]

    @property
    def is_leaf(self):
        return True

    @property
    def is_node(self):
        return False 

    def search(self, variants):
    	return self.sample

    def __repr__(self):
    	return "Leaf : %s " % self.sample
    	



    
        
        