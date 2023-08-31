#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Usage: simulate_synteny.py -s <INT> -k <INT> -g <INT> -a <FLT> -r <STR> -m <INT> -i <FLT> -e <FLT> -l <INT> [-z <STR> -M <INT> -h]

  [Options]
    -s, --simulations <INT>                     Number of simulations
    -k, --karyotype <INT>                       Number of initial chromosomes in simulations
    -g, --genes <INT>                           Number of genes in simulations
    -a, --rearrangements <INT>                  Number of rearrangements placed on the tree
    -r, --ratio <STR>                           Ratio of fission,fusion,translocation, e.g. 2,2,1
    -m, --model <INT>                           Model to infer back rearrangements with, 2 or 3
    -i, --missingness <FLT>                     Proportion of markers missing from each genome
    -e, --error <FLT>                           Proportion of markers assigned to the wrong chromosome
    -l, --leaves <INT>                          Leaves in tree
    -z, --anc_inference <STR>                   Infer ancestral chromosomes approximately (quick) or more accurately (slow) [default: quick]
    -M, --min_set <INT>                         Minimum synteny-set size [default: 1]
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
from dendropy.simulate import treesim
import itertools

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
i_arg = float(args['--missingness'])
e_arg = float(args['--error'])
l_arg = int(args['--leaves'])

ai_arg = args['--anc_inference']
ms_arg = int(args['--min_set'])


def generate_random_tree(leaves):
	simed_tree = treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5, num_extant_tips=leaves+1, 
		repeat_until_success=True) # simulate initial tree with leaves + 1
	newick_string = simed_tree.as_string(schema='newick')[5:-1] # convert to newick
	tree = ete3.Tree(newick_string) # load into ete3
	leaves = []
	for idx, node in enumerate(tree.traverse()): # collect leaf names
		if node.is_leaf():
				leaves.append(node.name)
	for combo in itertools.combinations(leaves, 2): # find one of the two leaves with an external branch length of 0
		up_length, down_length = sg.get_branch_distances(tree, combo[0], combo[1])
		if up_length + down_length == 0:
			leaf_to_remove = combo[0]
			break
	simed_tree.prune_taxa_with_labels([leaf_to_remove]) # remove that leaf
	# now we have a tree with the correct number of leaves
	# this way of making trees replaces the n + 1 speciation event with sampling
	newick_string = simed_tree.as_string(schema='newick')[5:-1]
	tree = ete3.Tree(newick_string)
	print(newick_string)
	tree = ete3.Tree(newick_string)
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

def branch_traversal(k_arg, g_arg, a_arg, r_arg):
	rearrangement_log = []
	genome_dict = collections.defaultdict(list)
	genome_dict["A"] = init_genome(k_arg, g_arg)
	genome_dict["B"] = copy.deepcopy(genome_dict["A"])
	for rearrangement in range(0, a_arg):
		rearrangement = random.choices(["fission", "fusion", "translocation"], weights=r_arg, k=1)[0]
		if rearrangement == "fission":
			genome_dict["B"], rearranged = implement_fission(genome_dict["B"])
			if rearranged:
				rearrangement_log.append(["A", "B", "fission"])
		elif rearrangement == "fusion":
			genome_dict["B"], rearranged = implement_fusion(genome_dict["B"])
			if rearranged:
				rearrangement_log.append(["A", "B", "fusion"])
		elif rearrangement == "translocation":
			genome_dict["B"], rearranged = implement_translocation(genome_dict["B"])
			if rearranged:
				rearrangement_log.append(["A", "B", "translocation"])
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

def syngraph_from_dict_for2(genome_dict):
	syngraph = sg.Syngraph()
	for tree_node in ["A", "B"]:
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
    def __init__(self, simulated_syngraph, tree, model, ancinf, min_set):
        self.syngraph = simulated_syngraph
        self.tree = tree
        self.minimum = min_set
        self.model = model
        self.ancinf = ancinf

def compare_genomes(simulated_syngraph, solved_syngraph, taxon, g_arg):
	markers_in_perfect_LGs = 0
	# lets say an approx LG contains 90% of markers from a single true LG and no more than 5% of markers from any other true LG
	markers_in_approx_LGs = 0
	taxon_simulated_genome = collections.defaultdict(set) # record each simulated chromosome as a set of markers
	taxon_solved_genome = collections.defaultdict(set) # record each inferred chromosome as a set of markers
	for graph_node_id in simulated_syngraph.nodes(): # populate the sets
		if taxon in simulated_syngraph.nodes[graph_node_id]['seqs_by_taxon'].keys():
			taxon_simulated_genome[simulated_syngraph.nodes[graph_node_id]['seqs_by_taxon'][taxon]].add(graph_node_id)
		if graph_node_id in solved_syngraph.nodes():
			if taxon in solved_syngraph.nodes[graph_node_id]['seqs_by_taxon'].keys():
				taxon_solved_genome[solved_syngraph.nodes[graph_node_id]['seqs_by_taxon'][taxon]].add(graph_node_id)
	for chromosome in taxon_simulated_genome.values(): # do we find a one - one match between a simulated and inferred LG?
		if chromosome in taxon_solved_genome.values():
			markers_in_perfect_LGs += len(chromosome)
	for solved_chromosome in taxon_solved_genome.values(): # do we find a fuzzy match?
		approx_criterea = 0
		for chromosome in taxon_simulated_genome.values():
			if len(solved_chromosome.intersection(chromosome)) / len(chromosome) >= 0.9:
				approx_criterea += 1
		for chromosome in taxon_simulated_genome.values():
			if 0.05 < len(solved_chromosome.intersection(chromosome)) / len(chromosome) < 0.9:
				approx_criterea = 0
		if approx_criterea == 1:
			markers_in_approx_LGs += len(solved_chromosome)
	return markers_in_perfect_LGs / g_arg, markers_in_approx_LGs / g_arg

def compare_rearrangements(rearrangement_log, inferred_log, l_arg):
	correct_branches = (2 * l_arg) - 4 # assume all branches are perfect to begin with
	total_branches = (2 * l_arg) - 4 # it is 4 rather than 2 because we don't infer rearrangements on branches beneath the root
	simulated_rearrangements = collections.defaultdict(list)
	inferred_rearrangements = collections.defaultdict(list)
	for entry in rearrangement_log:
		if entry[0] != "n0": # collect all branches with simmed rearrangements (expect those beneath the root)
			simulated_rearrangements[entry[0] + "_" + entry[1]].append(entry[2])
			simulated_rearrangements[entry[0] + "_" + entry[1]] = sorted(simulated_rearrangements[entry[0] + "_" + entry[1]])
	for entry in inferred_log:
		if entry == ['#parent', 'child', 'event', 'multiplicity', 'ref_seqs']:
			pass
		else:
			for i in range(0, int(entry[3])): # collect inferred rearrangements
				inferred_rearrangements[entry[0] + "_" + entry[1]].append(entry[2])
				inferred_rearrangements[entry[0] + "_" + entry[1]] = sorted(inferred_rearrangements[entry[0] + "_" + entry[1]])
	for key in simulated_rearrangements: # for every branch with simmed rearrangements, -1 from correct if inferred incorrectly
		if simulated_rearrangements[key] != inferred_rearrangements[key]:
			correct_branches -= 1
	# what have we calculated?
	# the proportion of branches with correctly inferred rearrangements **counts**
	return correct_branches / total_branches

def get_deepest_node(tree):
	n1_count = 0
	n2_count = 0
	for node in tree.traverse(strategy="preorder"):
		if node.name in ["n1", "n2"]:
			for desc in node.iter_descendants():
				if node.name == "n1":
					n1_count += 1
				elif node.name == "n2":
					n2_count += 1
	if n1_count >= n2_count:
		return "n1"
	else:
		return "n2"

def add_missingness(genome_dict, missingness, genes, tree):
	for node in tree.traverse(strategy="preorder"):
		if node.is_leaf():
			missing_markers = random.sample(range(1, 1 + genes), k= int(missingness * genes))
			for chromosome in genome_dict[node.name]:
				for m_marker in missing_markers:
					if m_marker in chromosome:
						chromosome.remove(m_marker)
			while set() in genome_dict[node.name]:
				genome_dict[node.name].remove(set())
	return genome_dict

def add_error(genome_dict, error, genes, tree):
	for node in tree.traverse(strategy="preorder"):
		if node.is_leaf():
			error_markers = random.sample(range(1, 1 + genes), k= int(error * genes))
			for chromosome in genome_dict[node.name]:
				for e_marker in error_markers:
					if e_marker in chromosome:
						chromosome.remove(e_marker)
			for e_marker in error_markers: # can return marker to orginal chrom, which is realistic
				new_chromosome = random.sample(genome_dict[node.name], k=1)[0]
				new_chromosome.add(e_marker)		
			while set() in genome_dict[node.name]:
				genome_dict[node.name].remove(set())
	return genome_dict


total_sims = 0

if l_arg <= 1:
	sys.exit("[X] Cannot simulate rearrangements with 0 or 1 leaves")

elif l_arg == 2:
	for i in range(0, s_arg):
		genome_dict, rearrangement_log = branch_traversal(k_arg, g_arg, a_arg, r_arg)
		simulated_syngraph = syngraph_from_dict_for2(genome_dict)
		LMSs, unassignable_markers = sg.get_LMSs(simulated_syngraph, ["A", "B"], 1)
		ios_A = sg.compact_synteny_1(simulated_syngraph, LMSs, "A")
		ios_B = sg.compact_synteny_1(simulated_syngraph, LMSs, "B")
		if m_arg == 2:
			inferred_log = sg.ffsd(sg.compact_synteny_2(ios_B, ios_A), [], 0)
		elif m_arg == 3:
			inferred_log = sg.ferretti(sg.compact_synteny_2(ios_B, ios_A), [], 0)
		inferred_rearrangements = 0
		for rearrangement in inferred_log:
			inferred_rearrangements += rearrangement[3]
		print(a_arg, inferred_rearrangements)

else:
	ALG_accuracy = 0
	approx_ALGs_accuracy = 0
	deep_ALG_accuracy = 0
	deep_approx_ALGs_accuracy = 0
	rearrangement_accuracy = 0
	for i in range(0, s_arg):
		tree = generate_random_tree(l_arg)
		deepest_node = get_deepest_node(tree)
		genome_dict, rearrangement_log = tree_traversal(tree, k_arg, g_arg, a_arg, r_arg)
		genome_dict = add_missingness(genome_dict, i_arg, g_arg, tree)
		genome_dict = add_error(genome_dict, e_arg, g_arg, tree)
		simulated_syngraph = syngraph_from_dict(genome_dict, tree)
		simulated_syngraph_with_ancestors = syngraph_with_ancestors_from_dict(genome_dict, tree)
		parameterObj = ParameterObj(simulated_syngraph, tree, m_arg, ai_arg, ms_arg)
		solved_syngraph, inferred_log = sg.tree_traversal(simulated_syngraph, parameterObj)
		total_sims += 1
		total_nodes = 0
		raw_ALG_accuracy = 0
		raw_approx_ALGs_accuracy = 0
		for node in tree.traverse(strategy="preorder"):
			if node.is_root():
				pass
			elif node.is_leaf():
				pass
			else:
				if node.name == deepest_node:
					deep_perfect_genome, deep_imperfect_genome = compare_genomes(simulated_syngraph_with_ancestors, solved_syngraph, 
						deepest_node, g_arg)
				perfect_genome, imperfect_genome = compare_genomes(simulated_syngraph_with_ancestors, solved_syngraph, 
					node.name, g_arg)
				raw_ALG_accuracy += perfect_genome
				raw_approx_ALGs_accuracy += imperfect_genome
				total_nodes += 1

		deep_ALG_accuracy += deep_perfect_genome
		deep_approx_ALGs_accuracy += deep_imperfect_genome

		ALG_accuracy += (raw_ALG_accuracy / total_nodes)
		approx_ALGs_accuracy += (raw_approx_ALGs_accuracy / total_nodes)

		print("ALG accuracy:", round(ALG_accuracy, 3), "/", total_sims)
		print("approx-ALG accuracy:", round(approx_ALGs_accuracy, 3), "/", total_sims)

		print("deep ALG accuracy:", round(deep_ALG_accuracy, 3), "/", total_sims)
		print("deep approx-ALG accuracy:", round(deep_approx_ALGs_accuracy, 3), "/", total_sims)

		rearrangement_accuracy += compare_rearrangements(rearrangement_log, inferred_log, l_arg)
		print("rearrangements accuracy:", round(rearrangement_accuracy, 3), "/", total_sims)

