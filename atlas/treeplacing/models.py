# import sys
# import os
# from os import path
# from atlas.schema import Variant
# from atlas.schema import VariantSet
# from atlas.schema import VariantCallSet

# sys.path.append(path.abspath("../"))


# class Placer(object):

#     """Placer"""

#     def __init__(self, root):
#         super(Placer, self).__init__()
#         self.root = root

#     def place(self, sample, verbose=False):
#         gvs = TypedVariant.objects(
#             call_set=CallSet.objects.get(
#                 name=sample)).distinct('name')
#         return self.root.search(variants=gvs, verbose=verbose)

# # class Tree(dict):
# #     """Tree is defined by a dict of nodes"""
# #     def __init__(self):
# #         super(Tree, self).__init__()


# def lazyprop(fn):
#     attr_name = '_lazy_' + fn.__name__

#     @property
#     def _lazyprop(self):
#         if not hasattr(self, attr_name):
#             setattr(self, attr_name, fn(self))
#         return getattr(self, attr_name)
#     return _lazyprop


# class Node(object):

#     """docstring for Node"""

#     def __init__(self, children=[]):
#         super(Node, self).__init__()
#         self.parent = None
#         self.children = children  # List of nodes
#         for child in self.children:
#             child.add_parent(self)
#         if self.is_node:
#             self.phylo_snps

#     def add_parent(self, parent):
#         self.parent = parent

#     def other_child(self, node):
#         for child in self.children:
#             if child != node:
#                 return child

#     @property
#     def samples(self):
#         samples = []
#         for child in self.children:
#             samples.extend(child.samples)
#         return samples  # List of sample below node in tree

#     @property
#     def num_samples(self):
#         return len(self.samples)

#     @property
#     def is_leaf(self):
#         return False

#     @property
#     def is_node(self):
#         return True

#     @property
#     def phylo_snps(self):
#         out_dict = {}
#         ingroup = VariantSet.objects(name__in=self.samples)
#         number_of_ingroup_samples = float(ingroup.count())

#         if self.parent:
#             # outgroup = VariantSet.objects(id__nin = [vs.id for vs in ingroup])
#             outgroup = VariantSet.objects(
#                 name__in=self.parent.other_child(self).samples)
#             number_of_outgroup_samples = float(outgroup.count())
#         else:
#             outgroup = []
#             number_of_outgroup_samples = 0

#         phylo_snp_names = Variant.objects(
#             variant_set__in=ingroup).distinct('name')
#         for name in phylo_snp_names:
#             count_ingroup = Variant.objects(
#                 name=name,
#                 variant_set__in=ingroup).count()
#             ingroup_freq = float(count_ingroup) / number_of_ingroup_samples
#             if number_of_outgroup_samples != 0:
#                 count_outgroup = Variant.objects(
#                     name=name,
#                     variant_set__in=outgroup).count()
#                 outgroup_freq = float(
#                     count_outgroup) / number_of_outgroup_samples
#             else:
#                 outgroup_freq = 0
#             out_dict[name] = ingroup_freq - outgroup_freq
#         return out_dict

#     def search(self, variants, verbose=False):
#         assert self.children[0].parent is not None
#         assert self.children[1].parent is not None
#         overlap = []
#         # Get the overlapping SNPS
#         variant_set = set(variants)
#         l0 = list(set(self.children[0].phylo_snps.keys()) & variant_set)
#         l1 = list(set(self.children[1].phylo_snps.keys()) & variant_set)
#         count0 = 0
#         count1 = 0
#         for k in l0:
#             count0 += self.children[0].phylo_snps[k]
#         for k in l1:
#             count1 += self.children[1].phylo_snps[k]
#         overlap = (count0, count1)
#         if verbose:
#             print (self.children[0], self.children[1], overlap)
#         if overlap[0] > overlap[1]:
#             return self.children[0].search(variants)
#         elif overlap[1] > overlap[0]:
#             return self.children[1].search(variants)
#         else:
#             return self.samples

#     def __repr__(self):
#         return "Node : %s " % ",".join(self.samples)


# class Leaf(Node):

#     def __init__(self, sample):

#         super(Leaf, self).__init__()
#         self.sample = sample

#     @property
#     def samples(self):
#         return [self.sample]

#     @property
#     def is_leaf(self):
#         return True

#     @property
#     def is_node(self):
#         return False

#     def search(self, variants):
#         return self.sample

#     def __repr__(self):
#         return "Leaf : %s " % self.sample
