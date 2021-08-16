#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Usage: simulate_synteny.py -s <INT> -k <INT> -g <INT> -a <FLT> -l <INT>

  [Options]
    -s, --simulations <INT>                     Number of simulations
    -k, --karyotype <INT>                       Number of initial chromosomes in simulations
    -g, --genes <INT>                           Number of genes in simulations
    -a, --rate <FLT>                            Rate of rearrangements per unit branch length
    -l, --leaves <INT>                          Leaves in tree
    -h, --help                                  Show this message

"""

import sys
from docopt import docopt
import collections
import copy
import numpy
from scipy import stats
import random
import ete3

sys.path.append('/ceph/users/amackintosh/software/syngraph/source/')
import syngraph as sg


random.seed(44)
numpy.random.seed(44)


args = docopt(__doc__)
s_arg = int(args['--simulations'])
k_arg = int(args['--karyotype'])
g_arg = int(args['--genes'])
a_arg = float(args['--rate'])
l_arg = int(args['--leaves'])


def generate_random_tree(leaves):
	tree = ete3.Tree()
	tree.populate(size=leaves, random_branches=True, branch_range=(0.01, 10))
	for idx, node in enumerate(tree.traverse()):
		if not node.is_leaf():
				node.name = "n%s" % idx
	print(tree.get_ascii(), "\n")
	return tree

def init_genome(k_arg, g_arg):
	ios = [set() for i in range(k_arg)]
	gene_name = 1
	for gene in range(g_arg):
		idx = random.sample(range(k_arg), 1)[0]
		ios[idx].add(gene_name)
		gene_name += 1
	if min([len(i) for i in ios]) == 0:
		sys.exit("[X] Generated a genome with empty chromosomes")
	return ios

def dice_roll(branch_length, rate):
	rearrangement_count = stats.poisson.rvs(branch_length*rate, size=1)[0]
	rearrangements = random.choices(["fission", "fusion"], weights=None, k=rearrangement_count)
	return rearrangements

def implement_fission(genome):
	chromo = random.choices(genome, weights=[len(chrom) for chrom in genome], k=1)[0]
	if len(chromo) == 1:
		return genome
	cut = random.sample(range(1, len(chromo)), 1)[0]
	new_chrom = set()
	for i in range(cut):
		sampled_gene = random.sample(chromo, 1)[0]
		new_chrom.add(sampled_gene)
		chromo.remove(sampled_gene)
	genome.append(new_chrom)
	return genome

def implement_fusion(genome):
	if len(genome) == 1:
		return genome
	chroms = random.sample(genome, 2)
	new_chrom = chroms[0].union(chroms[1])
	genome.remove(chroms[0])
	genome.remove(chroms[1])
	genome.append(new_chrom)
	return genome


def tree_traversal(tree, k_arg, g_arg, a_arg):
	genome_dict = collections.defaultdict(list)
	rearrangement_log = []
	for node in tree.traverse(strategy="preorder"):
		if not node.is_leaf():
			if node.is_root():
				genome_dict[node.name] = init_genome(k_arg, g_arg)
			for child in node.get_children():
				branch_length = node.get_distance(child)
				rearrangements = dice_roll(branch_length, a_arg)
				genome_dict[child.name] = copy.deepcopy(genome_dict[node.name])
				for rearrangement in rearrangements:
					if rearrangement == "fission":
						rearrangement_log.append([node.name, child.name, branch_length, "fission"])
						genome_dict[child.name] = implement_fission(genome_dict[child.name])
					elif rearrangement == "fusion":
						genome_dict[child.name] = implement_fusion(genome_dict[child.name])
						rearrangement_log.append([node.name, child.name, branch_length, "fusion"])
	return genome_dict, rearrangement_log

def syngraph_from_dict(genome_dict, tree):
	syngraph = sg.Syngraph()
	for tree_node in genome_dict:
		if tree_node in [leaf.name for leaf in tree.get_leaves()]:
			syngraph.graph['taxa'].add(tree_node)
			chromosome_number = 1
			for chromosome in genome_dict[tree_node]:
				for marker in chromosome:
					if marker not in syngraph:
						syngraph.add_node(marker, taxa=set(), seqs_by_taxon={})
					syngraph.nodes[marker]['taxa'].add(tree_node)
					syngraph.nodes[marker]['seqs_by_taxon'][tree_node] = tree_node + "_" + str(chromosome_number)
				chromosome_number += 1
	return syngraph

class ParameterObj():
    def __init__(self, simulated_syngraph, tree):
        self.syngraph = simulated_syngraph
        self.tree = tree
        self.minimum = 1
        self.model = 2


for i in range(0, s_arg):
	tree = generate_random_tree(l_arg)
	genome_dict, rearrangement_log = tree_traversal(tree, k_arg, g_arg, a_arg)
	simulated_syngraph = syngraph_from_dict(genome_dict, tree)
	parameterObj = ParameterObj(simulated_syngraph, tree)
	solved_syngraph, log = sg.tree_traversal(simulated_syngraph, parameterObj)