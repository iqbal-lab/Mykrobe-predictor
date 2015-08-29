from models import Tree 
from models import Node 
from models import Leaf 
from nose.tools import assert_raises

class BaseTest():

    def setUp(self):
        pass
        
    def teardown(self):
        pass

class TestNodes(BaseTest):

	def test_single_node_no_children(self):
		node = Node()
		assert node.children == []
		assert node.num_samples == 0
		assert node.is_leaf == False


	def test_node_triplet(self):
		node1 = Leaf(sample = 'C1')
		node2 = Leaf(sample = 'C2')
		root = Node(children = [node1, node2]) 
		assert root.num_samples == 2
		assert root.samples == ['C1', 'C2']

class TestMultiNode(TestNodes):

	def __init__(self):
		self.l1 = Leaf(sample = 'C1') 
		self.l2 = Leaf(sample = 'C2') 
		self.l3 = Leaf(sample = 'C3') 
		self.l4 = Leaf(sample = 'C4') 
		self.l5 = Leaf(sample = 'C5') 		

		self.node1 = Node(children = [l1, l2])
		self.node2 = Node(children = [l4, l5])
		self.node3 = Node(children = [node2, l3])
		self.root = Node(children = [node1, node3])


	def test_multi_node(self):
		assert self.root.num_samples == 5
		assert sorted(self.root.samples) == ['C1', 'C2', 'C3', 'C4', 'C5']