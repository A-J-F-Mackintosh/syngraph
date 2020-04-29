"""

Usage: syngraph synulate [-k <INT> -g <INT> -t <NWK> -i <INT> -h]

  [Options]
    -k, --init_karyotype <INT>                  Number of initial chromosomes [default: 10]
    -g, --genes_per_chromosome <INT>            Number of initial genes per chromosome [default: 50]
    -t, --tree <NWK>							Species tree in Newick format with branch lengths in generations
    -i, --inversion_rate <FLT>                  Per generation probability of an inversion happening [default: 0]
    -h, --help                                  Print this message

"""

import os
import sys
from docopt import docopt
import ete3
import matplotlib
import matplotlib.cm as cm
import pandas as pd
from tqdm import tqdm
import collections
import networkx as nx
from timeit import default_timer as timer
import random
import numpy as np
pd.set_option('display.max_rows', None)



def init_graph(k, g):
	k = int(k)
	g = int(g)
	#make an empty graph
	G_init = nx.Graph()
	#populate graph with edges
	for i in range(0, k):
		for j in range(1, g):
			G_init.add_edge(str((i*g)+j), str((i*g)+j+1))
	return G_init

def simulate(graph_init, tree, inversion_rate):
	# load in args
	inv_rate = float(inversion_rate)
	t = ete3.Tree(tree)
	# generate node numbers and store branch info in dicts
	branch_parentage = {} # the child node is the key and the parental node is the value
	branch_lengths = {} # the child node is the key and the length of the branch is the value
	branch_inversions = {} # the child node is the key and the number of inversion events is the value
	node_number = 0
	for tree_node in t.traverse(strategy='preorder'):
		if tree_node.name == "":
			node_number += 1
			tree_node.name = str(node_number)
		for child in tree_node.get_children():
			if child.name == "":
				node_number += 1
				child.name = str(node_number)
			branch_parentage[child.name] = tree_node.name
			branch_lengths[child.name] = child.get_distance(tree_node)
	print(t.get_ascii(), "\n")
	# record inversion events given the inversion_rate arg and the branch lengths
	total_inversions = 0
	for key in branch_lengths:
		branch_inversions[key] = np.random.poisson(lam=(int(branch_lengths[key])*inv_rate))
		total_inversions += branch_inversions[key]
	print("Total inversions: {}\n".format(total_inversions))
	# now iterate over branches and alter graphs given the inversions happening on that graph
	# first need to make a dict to store graph info
	graphs = {"1" : graph_init} # keys are nodes and values are graphs
	for tree_node in t.iter_descendants('preorder'):
		# get graph of parent
		parent_graph = graphs[branch_parentage[tree_node.name]]
		print("Inversions on the branch leading to node {}:".format(tree_node.name))
		current_graph = inversion(parent_graph, branch_inversions[tree_node.name])
		graphs[tree_node.name] = current_graph
		print("\n")
	return(graphs, t)

def inversion(graph, reps):
	reps_so_far = 0
	while reps_so_far < int(reps):
		# we assume that each chromosome has an equal chance of rearrangement
		chromosomes = []
		for chromosome in nx.connected_components(graph):
			chromosomes.append(chromosome)
		# select a chromosome to rearrange
		r_chromosome = random.choice(chromosomes)
		# select two genes to be breakpoints
		r_genes = random.sample(r_chromosome, 2)
		print("Inversion between marker {} and {}".format(r_genes[0], r_genes[1]))
		# get nodes captured in the inversion
		inversion_nodes_list = [node for node in nx.shortest_simple_paths(graph, r_genes[0], r_genes[1])][0]
		inversion_nodes_set = set(inversion_nodes_list)
		# record new graph edges
		new_edges = []
		# loop over all edges in the graph
		for n1, n2 in graph.edges():
			# difference can be both nodes -> therefore not in the inversion 
		    diff = {n1,n2}.difference(inversion_nodes_set)
		    if len(diff) == 2:
		        new_edges.append(diff)
		    # or one node -> therefore one node is right next to the inversion
		    elif len(diff) == 1:
		    	# n is a node next to the inversion, graph[n] are the nodes that n connects to
		    	# inversion_nodes_list[0] and inversion_nodes_list[-1] are the ends of the inversion
		        n = list(diff)[0]
		        if inversion_nodes_list[0] in graph[n] or n in graph[inversion_nodes_list[0]]:
		            new_edges.append((n, inversion_nodes_list[-1]))
		        if inversion_nodes_list[-1] in graph[n] or n in graph[inversion_nodes_list[-1]]:
		            new_edges.append((n, inversion_nodes_list[0]))
		    # or no nodes -> therefore fully within the inversion
		    else:
		        new_edges.append((n1, n2))
		# make a graph from new_edges
		reps_so_far += 1
		new_graph = nx.Graph()
		new_graph.add_edges_from(new_edges)
		graph = new_graph
	return(graph)

def reconstruct(graphs, tree):
	# make a single syngraph containing all edges of leaf graphs
	syngraph = nx.Graph()
	for tree_node in tree.traverse(strategy='postorder'):
		if tree_node.is_leaf():
			for u, v in graphs[tree_node.name].edges():
				if syngraph.has_edge(u, v):
					syngraph[u][v]['taxa'].add(tree_node.name)
				else:
					syngraph.add_edge(u, v, taxa=set())
					syngraph[u][v]['taxa'].add(tree_node.name)
	# make a recon graph to store reconstructed edges
	recongraph = nx.Graph()
	# now identify the ancestral tree node we want to reconstruct
	# this can be defined as the node with the greatest number of children that is not the root
	max_leaves = 0
	for tree_node in tree.iter_descendants('preorder'):
		if len(tree_node.get_leaves()) > max_leaves:
			anc_tree_node = tree_node.name
			max_leaves = len(tree_node.get_leaves())
	# and we can get the anc_tree_node edges from the graphs dict
	anc_tree_node_graph = graphs[anc_tree_node]
	# okie doke, time to reconstruct
	# iterate over all graph nodes
	for graph_node_id in syngraph.nodes:
		# get edges associated a graph node
		edges_verbose = syngraph.edges([graph_node_id], data=True)
		# iterate up through the tree, storing reconstructions in a dict
		reconstruct_dict = {}
		for tree_node in tree.traverse(strategy='postorder'):
			if not tree_node.is_leaf():
				children = tree_node.get_children()
				# collect child states 
				child_counter = 0
				child_state_1 = set()
				child_state_2 = set()
				for child in children:
					child_counter += 1
					# if child is a leaf then collect edge state from edges_verbose
					if len(child.get_leaves()) == 1:
						for edge in edges_verbose:
							if child.name in edge[2]['taxa']:
								if child_counter == 1:
									child_state_1.add(edge[0:2])
								elif child_counter == 2:
									child_state_2.add(edge[0:2])
					# if child is an internal node then look in the reconstruct_dict for edge_state
					else:
						if child_counter == 1:
							child_state_1 = reconstruct_dict[child.name]
						elif child_counter == 2:
							child_state_2 = reconstruct_dict[child.name]
				# if a child has only one edge then add a terminal edge for reconstruction
				for child_state_n in child_state_1, child_state_2:
					if len(child_state_n) == 1:
						child_state_n.add((graph_node_id, "Terminal_%s" % graph_node_id))
				# reconstruct current tree_node from children
				if len(child_state_1) == 2 and len(child_state_2) == 2:
					reconstruct_dict[tree_node.name] = child_state_1.union(child_state_2)
				else:
					if len(child_state_1.intersection(child_state_2)) >= 2:
						reconstruct_dict[tree_node.name] = child_state_1.intersection(child_state_2)
					else:
						reconstruct_dict[tree_node.name] = child_state_1.union(child_state_2)
		# iterate down through the tree, updating the reconstructions_dict
		for tree_node in tree.traverse(strategy='preorder'):
			if not tree_node.is_leaf():
				tree_node_state = reconstruct_dict[tree_node.name]
				children = tree_node.get_children()
				child_counter = 0
				for child in children:
					if not child.is_leaf():
						child_counter += 1
						if child_counter == 1:
							child_state_1 = reconstruct_dict[child.name]
							if len(tree_node_state.intersection(child_state_1)) >= 2:
								reconstruct_dict[child.name] = tree_node_state.intersection(child_state_1)
							elif child_counter == 2:
								child_state_2 =  reconstruct_dict[child.name]
								if len(tree_node_state.intersection(child_state_2)) >= 2:
									reconstruct_dict[child.name] = tree_node_state.intersection(child_state_2)
		# add edges to recon_graph
		for e in reconstruct_dict[anc_tree_node]:
			if "Terminal_" in e[0] or "Terminal_" in e[1]:
				pass
			else:
				recongraph.add_edge(e[0], e[1])
	# finally we compare the recongraph to the anc_tree_node graph
	return(anc_tree_node_graph, recongraph)

def compare(true, estimate):
	# sensitivity - what percent of edges in true are also contained in estimate
	total_t_edges = 0
	shared_edges = 0
	for true_e in true.edges():
		total_t_edges += 1
		for estimated_e in estimate.edges():
			if set(true_e) == set(estimated_e):
				shared_edges += 1
	sensitivity = (shared_edges / total_t_edges) * 100
	print(sensitivity)
	# specificty - what percent of edges in estimate are also contained in true
	total_e_edges = 0
	for estimated_e in estimate.edges():
		total_e_edges += 1
	specificty = (shared_edges / total_e_edges) * 100
	print(specificty)

# output each event of each lineage - prob selected, event, subgraph, nodes
# intervals and genes can be seperate - start with just genes
# fragment assemblies to a target N50 to test reconstruction algorithm
# can also incorporate othology errors and deletions/missingness
# first gen should be a speciation event to ensure no inversions before anc node



args = docopt(__doc__)
print(args, "\n")

print("Initialising graph with {} chromosomes and {} genes per chromosome".format(args['--init_karyotype'], args['--genes_per_chromosome']))
G_init = init_graph(args['--init_karyotype'], args['--genes_per_chromosome'])
print("Simulating chromosome evolution")
graphs_for_reconstruction, tree_for_reconstruction = simulate(G_init, args['--tree'], args['--inversion_rate'])
print("Reconstrucing synteny and gene order for the most ancestral tree node that is not the root")
true_graph, estimated_graph = reconstruct(graphs_for_reconstruction, tree_for_reconstruction)
compare(true_graph, estimated_graph)




