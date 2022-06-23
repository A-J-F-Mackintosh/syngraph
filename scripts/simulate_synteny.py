#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Usage: simulate_synteny.py -s <INT> -k <INT> -g <INT> -a <FLT> -r <STR> -m <INT> -l <INT>

  [Options]
    -s, --simulations <INT>                     Number of simulations
    -k, --karyotype <INT>                       Number of initial chromosomes in simulations
    -g, --genes <INT>                           Number of genes in simulations
    -a, --rearrangements <INT>                  Number of rearrangements placed on the tree
    -r, --ratio <STR>                           Ratio of fission,fusion,translocation, e.g. 2,2,1
    -m, --model <INT>                           Model to infer back rearrangements with, 2 or 3
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
a_arg = int(args['--rearrangements'])
r_arg = [int(i) for i in (args['--ratio']).split(",")]
m_arg = int(args['--model'])
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

def throw_rearrangements(rearrangements, parent_child, branch_lengths):
	return random.choices(parent_child, weights=branch_lengths, k=rearrangements)

def implement_fission(genome):
	chromo = random.choices(genome, weights=[len(chrom) for chrom in genome], k=1)[0]
	if len(chromo) == 1:
		return genome, False
	cut = random.sample(range(1, len(chromo)), 1)[0]
	new_chrom = set()
	for i in range(cut):
		sampled_gene = random.sample(chromo, 1)[0]
		new_chrom.add(sampled_gene)
		chromo.remove(sampled_gene)
	genome.append(new_chrom)
	return genome, True

def implement_fusion(genome):
	if len(genome) == 1:
		return genome, False
	chroms = random.sample(genome, 2)
	new_chrom = chroms[0].union(chroms[1])
	genome.remove(chroms[0])
	genome.remove(chroms[1])
	genome.append(new_chrom)
	return genome, True

def implement_translocation(genome):
	#print("##", len(genome), sum([len(j) for j in genome]), genome)
	if len(genome) == 1:
		return genome, False
	chrom_A, chrom_B = random.sample(genome, 2)
	if len(chrom_A) == 1 or len(chrom_B) == 1:
		return genome, False
	chrom_A_bp = random.sample(range(1, len(chrom_A)), 1)[0]
	chrom_B_bp = random.sample(range(1, len(chrom_B)), 1)[0]
	new_chrom = set()
	for i in range(chrom_A_bp):
		sampled_gene = random.sample(chrom_A, 1)[0]
		new_chrom.add(sampled_gene)
		chrom_A.remove(sampled_gene)
	for i in range(chrom_B_bp):
		sampled_gene = random.sample(chrom_B, 1)[0]
		new_chrom.add(sampled_gene)
		chrom_B.remove(sampled_gene)
	genome.append(chrom_A.union(chrom_B))
	genome.remove(chrom_A)
	genome.remove(chrom_B)
	genome.append(new_chrom)
	#print("###", len(genome), sum([len(j) for j in genome]), genome)
	return genome, True

def tree_traversal(tree, k_arg, g_arg, a_arg, r_arg):
	genome_dict = collections.defaultdict(list)
	rearrangement_log = []
	parent_child = []
	branch_lengths = []
	for node in tree.traverse(strategy="preorder"):
		if not node.is_leaf():
			if node.is_root():
				genome_dict[node.name] = init_genome(k_arg, g_arg)
			for child in node.get_children():
				branch_length = node.get_distance(child)
				parent_child.append((node.name, child.name))
				branch_lengths.append(branch_length)
	pre_rearrangement_log = throw_rearrangements(a_arg, parent_child, branch_lengths)
	for node in tree.traverse(strategy="preorder"):
		if not node.is_leaf():
			for child in node.get_children():
				genome_dict[child.name] = copy.deepcopy(genome_dict[node.name])
				for entry in pre_rearrangement_log:
					if (node.name, child.name) == entry:
						rearrangement = random.choices(["fission", "fusion", "translocation"], weights=r_arg, k=1)[0]
						if rearrangement == "fission":
							genome_dict[child.name], rearranged = implement_fission(genome_dict[child.name])
							if not node.is_root() and rearranged:
								rearrangement_log.append([node.name, child.name, "fission"])
						elif rearrangement == "fusion":
							genome_dict[child.name], rearranged = implement_fusion(genome_dict[child.name])
							if not node.is_root() and rearranged:
								rearrangement_log.append([node.name, child.name, "fusion"])
						elif rearrangement == "translocation":
							genome_dict[child.name], rearranged = implement_translocation(genome_dict[child.name])
							if not node.is_root() and rearranged:
								rearrangement_log.append([node.name, child.name, "translocation"])
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

def syngraph_with_ancestors_from_dict(genome_dict, tree):
	syngraph = sg.Syngraph()
	for node in tree.traverse(strategy="preorder"):
		if not node.is_root():
			tree_node = node.name
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
    def __init__(self, simulated_syngraph, tree, model):
        self.syngraph = simulated_syngraph
        self.tree = tree
        self.minimum = 1
        self.model = model

def compare_genomes(simulated_syngraph, solved_syngraph):
	for taxon in simulated_syngraph.graph['taxa']:
		taxon_simulated_genome = collections.defaultdict(set)
		taxon_solved_genome = collections.defaultdict(set)
		for graph_node_id in simulated_syngraph.nodes():
			taxon_simulated_genome[simulated_syngraph.nodes[graph_node_id]['seqs_by_taxon'][taxon]].add(graph_node_id)
			taxon_solved_genome[solved_syngraph.nodes[graph_node_id]['seqs_by_taxon'][taxon]].add(graph_node_id)
		for chromosome in taxon_simulated_genome.values():
			if chromosome in taxon_solved_genome.values():
				pass
			else:
				return 0
	return 1

def compare_rearrangements(rearrangement_log, inferred_log):
	simulated_rearrangements = []
	inferred_rearrangements = []
	for entry in rearrangement_log:
		simulated_rearrangements.append(entry)
	for entry in inferred_log:
		if entry == ['parent', 'child', 'event', 'multiplicity', 'markers']:
			pass
		else:
			for i in range(0, int(entry[3])):
				inferred_rearrangements.append(entry[0:3])
	simulated_rearrangements.sort()
	inferred_rearrangements.sort()
	if simulated_rearrangements == inferred_rearrangements:
		return 1
	else:
		return 0


total_sims = 0
correctly_inferred_genomes = 0
correctly_inferred_histories = 0
for i in range(0, s_arg):
	tree = generate_random_tree(l_arg)
	genome_dict, rearrangement_log = tree_traversal(tree, k_arg, g_arg, a_arg, r_arg)
	simulated_syngraph = syngraph_from_dict(genome_dict, tree)
	simulated_syngraph_with_ancestors = syngraph_with_ancestors_from_dict(genome_dict, tree)
	parameterObj = ParameterObj(simulated_syngraph, tree, m_arg)
	solved_syngraph, inferred_log = sg.tree_traversal(simulated_syngraph, parameterObj)
	#print("simulated:", rearrangement_log)
	#print("inferred:", inferred_log)
	total_sims += 1
	correctly_inferred_genomes += compare_genomes(simulated_syngraph_with_ancestors, solved_syngraph)
	print("genomes:", correctly_inferred_genomes, "/", total_sims)
	correctly_inferred_histories += compare_rearrangements(rearrangement_log, inferred_log)
	print("histories:", correctly_inferred_histories, "/", total_sims)
