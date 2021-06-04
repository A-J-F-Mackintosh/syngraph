import sys
import itertools
import more_itertools
import matplotlib
import matplotlib.cm as cm
import pandas as pd
pd.set_option('display.max_rows', None)
from tqdm import tqdm
import networkx as nx
import collections
import functools
from scipy import stats
from operator import attrgetter
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
import copy

'''
[To Do]

- Alex:
    - test for consistency of taxon-names in tree/filenames in parameterObj
'''


def load_markerObjs(parameterObj):
    '''
    marker := orthogroup/busco
    locus := protein i.e. position of "marker" in genome of taxon

    - independent of "status", all markers/loci should be parsed into markerObjs
    - sort them out when loading into syngraph (syngraph has to know about allowed statuses, pass as args upon __init__)

    '''
    if parameterObj.label == 'busco':
        tmp_markerObjs = []
        status_allowed = {'Complete'}
        if parameterObj.fragmented:
            status_allowed.add('Fragmented') 
        for label, infiles in parameterObj.infiles_by_label.items(): # this might be needed later when OrthoFinder parsing is implemented...
            for infile in infiles:
                taxon_markerObjs = []
                taxon = infile.split("/")[-1].split(".")[0]
                df = pd.read_csv(infile, 
                        sep='\t', 
                        names=['name', 'status', 'seq', 'start', 'end', 'sign'], 
                        dtype={'name': str, 'status': str , 'seq': str, 'start': float, 'end': float, 'sign': str}
                        ).sort_values(['seq', 'start'], ascending=[True, True])
                for idx, (name, status, seq, start, end, sign) in enumerate(df.values.tolist()):
                    if status in status_allowed:
                        if parameterObj.sign:
                            if sign == "+":
                                coords_by_name = {"%s_head" % name: start, "%s_tail" % name: end}
                            elif sign == "-":
                                coords_by_name = {"%s_head" % name: end, "%s_tail" % name: start}
                            else:
                                sys.exit("[X] Sign must be '+' or '-' in line %r in file %r, not %r" % (idx, infile, sign))
                            for _name, _coord in coords_by_name.items():
                                markerObj= MarkerObj(name=_name, status=status, taxon=taxon, seq=taxon+"_"+seq, coord=_coord)
                                taxon_markerObjs.append(markerObj)
                        else:
                            markerObj = MarkerObj(name=name, status=status, taxon=taxon, seq=taxon+"_"+seq, coord=start)
                            taxon_markerObjs.append(markerObj)
                taxon_markerObjs = sorted(taxon_markerObjs, key=attrgetter('seq', 'coord'))
                tmp_markerObjs.extend(taxon_markerObjs)
            if parameterObj.missing == True:
                return tmp_markerObjs
            else:
                markerObjs = []
                counter = collections.Counter([tmp.name for tmp in tmp_markerObjs])
                for markerObj in tmp_markerObjs:
                    if counter[markerObj.name] == len(infiles):
                        markerObjs.append(markerObj)
                return markerObjs
    elif parameterObj.label == 'ortho':
        '''
        not dealing yet with:
        - missing
        - duplicates
        '''
        orthogroup_df = pd.read_csv(parameterObj.orthogroups, sep='\t')
        header = orthogroup_df.columns.values.tolist()
        marker_by_locus_by_taxon = collections.defaultdict(dict)
        locus_by_taxon_by_marker = collections.defaultdict(dict)
        locus_by_type_by_taxon = collections.defaultdict(lambda: collections.defaultdict(set))
        taxa = []
        for row in orthogroup_df.values.tolist():
            marker = row[0]
            for taxon, loci in zip(header[1:], row[1:]):
                if isinstance(loci, str):
                    for locus in loci.split(", "):                        
                        marker_by_locus_by_taxon[taxon][locus] = marker
                        locus_by_taxon_by_marker[marker][taxon] = locus
                        locus_by_type_by_taxon[taxon]['tsv'].add(locus)
        markerObjs = []
        locus_types = ['bed', 'tsv', 'both', 'not_bed', 'not_tsv']
        for label, infiles in parameterObj.infiles_by_label.items():
            for infile in tqdm(infiles, total=len(infiles), desc="[%] ", ncols=100):
                taxon = infile.split("/")[-1].split(".")[0]
                taxa.append(taxon)
                bed_df = pd.read_csv(infile, sep='\t', 
                    names=['seq', 'start', 'end', 'desc'], 
                    dtype={'seq': str, 'start': int , 'end': int, 'desc': str})
                for idx, (seq, start, end, locus) in enumerate(bed_df.values.tolist()):
                    locus_by_type_by_taxon[taxon]['bed'].add(locus)
                    marker = marker_by_locus_by_taxon[taxon].get(locus, None)
                    if not marker is None:
                        if len(locus_by_taxon_by_marker[marker][taxon]) > 1:
                            status = 'Duplicted'
                        elif len(locus_by_taxon_by_marker[marker][taxon]) == 1:
                            status = 'Single'
                        else:
                            sys.exit("[X] %s" % (idx, seq, start, end, locus))
                        markerObj = MarkerObj(name=marker, desc=locus, status=status, taxon=taxon, seq=taxon+"_"+seq, coord=start)
                        markerObjs.append(markerObj)
                locus_by_type_by_taxon[taxon]['both'] = locus_by_type_by_taxon[taxon]['tsv'].intersection(locus_by_type_by_taxon[taxon]['bed'])
                locus_by_type_by_taxon[taxon]['not_bed'] = locus_by_type_by_taxon[taxon]['tsv'] - locus_by_type_by_taxon[taxon]['bed']
                locus_by_type_by_taxon[taxon]['not_tsv'] = locus_by_type_by_taxon[taxon]['bed'] - locus_by_type_by_taxon[taxon]['tsv']
        # Sanity check whether all loci were found in both TSV and BED
        for taxon in taxa:
            if not len(locus_by_type_by_taxon[taxon]['tsv']) == len(locus_by_type_by_taxon[taxon]['both']):
                print("[-] Not all loci in Orthogroups file were found in BED file: %s" % taxon)
                for locus_type in locus_types:
                    print("[-]\t%s: %s [%s...]" % (locus_type, len(locus_by_type_by_taxon[taxon][locus_type]), ",".join(list(locus_by_type_by_taxon[taxon][locus_type])[0:3])))
        return markerObjs

# def fitch(states_by_taxon, number_of_states, tree): # states_by_taxon_node should be a dict, with keys as taxon and states as sets
#     states_by_tree_node = {}
#     for tree_node in tree.traverse(strategy='postorder'):
#         if tree_node.name in states_by_taxon:
#             states_by_tree_node[tree_node.name] = states_by_taxon[tree_node.name]
#         elif not tree_node.is_leaf():
#             intersection = set.intersection(*[states_by_tree_node.get(child_node.name, {False}) for child_node in tree_node.get_children()])
#             if len(intersection) >= number_of_states:
#                 states_by_tree_node[tree_node.name] = intersection
#             else:
#                 states_by_tree_node[tree_node.name] = set.union(
#                     *[states_by_tree_node.get(child_node.name, {False}) for child_node in tree_node.get_children()])
#         else:
#             pass
#     for tree_node in tree.traverse(strategy='levelorder'):
#         if not tree_node.is_root():
#             parent_tree_node = tree_node.up
#             intersection = states_by_tree_node.get(parent_tree_node.name, {False}).intersection(states_by_tree_node.get(tree_node.name, {False}))
#             if len(intersection) >= number_of_states:
#                 states_by_tree_node[tree_node.name] = intersection
#     return(states_by_tree_node)

# def reconstruct_syngraphs_by_tree_node(syngraph, tree, algorithm='fitch'):
#     '''
#     - input: syngraph, tree
#     - output: novel graphs with fitch edges for each internal tree node
#     '''
#     if algorithm == 'fitch':
#         edges_by_tree_node_by_graph_node = collections.defaultdict(dict) # nested dict, graph_node --> taxon --> edges
#         edges_by_tree_node = collections.defaultdict(list)
#         taxa_by_tree_node = collections.defaultdict(set)
#         for graph_node_id in syngraph.nodes:
#             edge_sets_by_taxon = syngraph.get_target_edge_sets_by_taxon(graph_node_id)
#             #print('edge_sets_by_taxon', edge_sets_by_taxon)
#             edges_by_tree_node_by_graph_node[graph_node_id] = fitch(edge_sets_by_taxon, 2, tree)
#         for graph_node_id, _edges_by_tree_node in edges_by_tree_node_by_graph_node.items():
#             for tree_node, edges in _edges_by_tree_node.items():
#                 #print("tree_node", tree_node)
#                 #print("edges", edges)
#                 for (u, v) in edges:
#                     taxa_under_node = set([node.name for node in (tree&tree_node).iter_leaves()])
#                     taxa_by_tree_node[tree_node].update(taxa_under_node)
#                     edge_taxa = []
#                     for taxon in taxa_under_node:
#                         if frozenset([u, v]) in _edges_by_tree_node[taxon]:
#                             edge_taxa.append(taxon)
#                     edges_by_tree_node[tree_node].append((u, v, {'taxa': edge_taxa}))
#         syngraph_by_tree_node = {}
#         for tree_node, edges in edges_by_tree_node.items():
#             syngraph_by_tree_node[tree_node] = Syngraph()
#             syngraph_by_tree_node[tree_node].from_edges(edges, taxa=taxa_by_tree_node[tree_node])
#             syngraph_by_tree_node[tree_node].show_recon_metrics(True, tree_node)
#             #syngraph_by_tree_node[tree_node].plot(outprefix="node_%s" % tree_node)
#     return syngraph_by_tree_node

# def plot_histogram(x, out_f):
#     fig, ax = plt.subplots(figsize=(14, 5))
#     hist, bins = np.histogram(x, bins=50)
#     width = 0.7 * (bins[1] - bins[0])
#     center = (bins[:-1] + bins[1:]) / 2
#     ax.bar(center, hist, align='center', width=width)
#     fig.savefig('%s.png' % out_f, format="png")

# def get_hex_colours_by_taxon(taxa, cmap='Spectral'):
#     return {taxon: matplotlib.colors.rgb2hex(cm.get_cmap(cmap, len(taxa))(i)[:3]) for i, taxon in enumerate(sorted(taxa))}

#############################################################################################
################################################################################:)###########
### below are functions for implementing fusion+fission models ##############################
################################################################################:(###########
#############################################################################################

# a function for getting the closest outgroup to a tree node
# outgroup must be in 'available_taxa' and cannot be a descenant of the tree_node or the tree_nodes itself.
# think about using the outgroup with the shortest rearrangement distance
def get_closest_outgroup(tree, tree_node, available_taxa):
    closest_taxon_so_far = None
    closest_distance_so_far = float("inf")
    for some_tree_node in tree.search_nodes(): # traverse in some order?
        if some_tree_node.name in available_taxa:
            if not some_tree_node.name in [descendant.name for descendant in tree_node.get_descendants()]:
                if not some_tree_node.name == tree_node.name:
                    if tree_node.get_distance(some_tree_node.name) < closest_distance_so_far:
                        closest_distance_so_far = tree_node.get_distance(some_tree_node.name)
                        closest_taxon_so_far = some_tree_node.name
    return closest_taxon_so_far

# a function which gives branch lengths by orientation, i.e. how much up and how much down
# the median/target will be A, the child/outgroup/parent will be B
def get_branch_distances(tree, tree_node_A, tree_node_B):
    mrca = tree.get_common_ancestor(tree_node_A, tree_node_B)
    up_distance = mrca.get_distance(tree_node_A)
    down_distance = mrca.get_distance(tree_node_B)
    return up_distance, down_distance

def compact_synteny_1(syngraph, LMSs, taxon):
    querychrom2index = collections.defaultdict(set)
    for graph_node_id in syngraph.nodes():
        for LMS in LMSs:
            if graph_node_id in LMSs[LMS]:
                querychrom2index[syngraph.nodes[graph_node_id]['seqs_by_taxon'][taxon]].add(LMS)
    return querychrom2index

def compact_synteny_2(taxon_A, taxon_B):
    querychrom2index = collections.defaultdict(set)
    for B_chrom in taxon_B:
        for index in taxon_B[B_chrom]:
            for A_chrom in taxon_A:
                if index in taxon_A[A_chrom]:
                    querychrom2index[B_chrom].add(A_chrom)
    return querychrom2index

def label_genome(tree_node, genome):
    # chroms have no names when the genome is represented as a list of sets so am giving them a number
    labelled_genome = {}
    chrom_number = 0
    for chrom in genome:
        chrom_number += 1
        labelled_genome[tree_node + "_" + str(chrom_number)] = chrom
    return labelled_genome

def check_for_fusions(instance_of_synteny, fusion_log):
    # can be sped up
    for combo in itertools.combinations(instance_of_synteny.keys(), 2):
        if instance_of_synteny[combo[0]].intersection(instance_of_synteny[combo[1]]):
            instance_of_synteny[combo[0]+ "_" + combo[1]] = instance_of_synteny[combo[0]].union(instance_of_synteny[combo[1]])
            fusion_log.append((combo[0], combo[1]))
            del instance_of_synteny[combo[0]]
            del instance_of_synteny[combo[1]]
            return check_for_fusions(instance_of_synteny, fusion_log)
    return instance_of_synteny, fusion_log

def check_for_fissions(instance_of_synteny):
    # log something like chrom_name and number of fissions?
    fission_log = []
    for chrom in instance_of_synteny:
        indices = len(instance_of_synteny[chrom])
        if indices > 1:
            fission_log.append(tuple([chrom, indices]))
    return fission_log

def ffsd(instance_of_synteny):
    #print("#####", instance_of_synteny)
    instance_of_synteny, fusion_log = check_for_fusions(instance_of_synteny, [])
    fission_log = check_for_fissions(instance_of_synteny)
    return fusion_log, fission_log

def get_LMSs(syngraph, list_of_taxa, minimum):
    # could be a numpy solution...
    Linked_marker_sets = collections.defaultdict(set)
    Unassignable_markers = set()
    for graph_node_id in syngraph.nodes():
        target_seqs_by_taxon_set = set()
        for taxon in list_of_taxa:
            if taxon in syngraph.nodes[graph_node_id]['seqs_by_taxon']:
                target_seqs_by_taxon_set.add(syngraph.nodes[graph_node_id]['seqs_by_taxon'][taxon])
        # are seqs_by_taxon unique to each taxon? If not the above will not work
        target_seqs_by_taxon_set = frozenset(target_seqs_by_taxon_set)
        if len(target_seqs_by_taxon_set) == len(list_of_taxa):
            Linked_marker_sets[target_seqs_by_taxon_set].add(graph_node_id)
        else:
            Unassignable_markers.add(graph_node_id)
    Filtered_linked_marker_sets = {}
    Filtered_LMS_count = 0
    for LMS in Linked_marker_sets:
        if len(Linked_marker_sets[LMS]) >= minimum:
            Filtered_linked_marker_sets["LMS_" + str(Filtered_LMS_count)] = Linked_marker_sets[LMS]
            Filtered_LMS_count += 1
        else:
            for graph_node_id in Linked_marker_sets[LMS]:
                Unassignable_markers.add(graph_node_id)
    return Filtered_linked_marker_sets, Unassignable_markers

# bubble will add to the set until nothing else can be added
def bubble(i, connected_components):
    bubbled = False
    for j in range(len(connected_components)):
        if connected_components[j] != set() and i != j:
            if len(connected_components[i].intersection(connected_components[j])) > 0:
                connected_components[i] = connected_components[i].union(connected_components[j])
                connected_components[j] = set()
                bubbled = True
    if bubbled:
        return bubble(i, connected_components)
    else:
        return connected_components

def generate_connected_components(ios_1, ios_2, ios_3):
    # get unique chromosomes / i.e. LMS relationships
    connected_components = list(set([frozenset(x) for x in list(ios_1.values()) + list(ios_3.values()) + list(ios_2.values())]))
    # convert back to sets
    connected_components = [set(component) for component in connected_components]
    # for each LMs relationship, bubble
    for i in range(len(connected_components)):
        if connected_components[i] != set():
            connected_components = bubble(i, connected_components)
    # get rid of empty sets then return
    connected_components = [component for component in connected_components if component != set()]
    return connected_components

def generate_parsimony_genome(ios_1, ios_2, ios_3, LMSs):
    parsimony_genome = []
    synteny_counts = collections.defaultdict(int)
    observed_chromosomes = list(ios_1.values()) + list(ios_3.values()) + list(ios_2.values())
    for combo in itertools.combinations(LMSs.keys(), 2):
        for observed_chromosome in observed_chromosomes:
            if combo[0] in observed_chromosome and combo[1] in observed_chromosome:
                synteny_counts[frozenset({combo[0], combo[1]})] += 1
    for combo_key in synteny_counts.keys():
        if synteny_counts[combo_key] == 2:
            parsimony_genome.append(set(combo_key))
    for i in range(len(parsimony_genome)):
        if parsimony_genome[i] != set():
            parsimony_genome = bubble(i, parsimony_genome)
    # add in shared LMSs
    for LMS in LMSs.keys():
        already_in = False
        for chrom in parsimony_genome:
            if LMS in chrom:
                already_in = True
        if not already_in:
            parsimony_genome.append({LMS})
    # get rid of empty sets then return
    parsimony_genome = [component for component in parsimony_genome if component != set()]
    return parsimony_genome

def evaluate_genome_with_parsimony(branch_lengths, ios_1, ios_2, ios_3, possible_median, 
    best_total, best_max_branch_rate, best_genome, best_fusion_log, best_fission_log):
    total_fissions = 0
    total_fusions = 0
    max_branch_rate = 0
    fusion_log = [[], [], []]
    fission_log = [[], [], []]
    for ios, i in zip([ios_1, ios_2, ios_3], range(0, 3)):
        fusion_log[i], fission_log[i] = ffsd(compact_synteny_2(ios, possible_median))
        fusions = len(fusion_log[i])
        fissions = sum(entry[1]-1 for entry in fission_log[i])
        branch_rate = (fusions + fissions) / (branch_lengths["up"][i] + branch_lengths["down"][i])
        if branch_rate > max_branch_rate:
            max_branch_rate = branch_rate
        total_fissions += fissions
        total_fusions += fusions
        if best_total < total_fissions + total_fissions:
            return best_genome, best_total, best_max_branch_rate, best_fusion_log, best_fission_log
    if best_total > total_fissions + total_fusions or \
    best_total == total_fissions + total_fusions and max_branch_rate < best_max_branch_rate:
        best_total = total_fissions + total_fusions
        best_max_branch_rate = max_branch_rate
        best_genome = possible_median
        best_fusion_log = fusion_log
        best_fission_log = fission_log
    return best_genome, best_total, best_max_branch_rate, best_fusion_log, best_fission_log

def exhaustive_solver(tree_node, connected_components, ios_1, ios_2, ios_3, branch_lengths):
    # if there are no 'tangles' then there are no events
    if max([len(connected_component) for connected_component in connected_components]) == 1:
        return connected_components
    # else, generate all possible candidate genomes from the connected components
    else:
        possible_medians = []
        partition_dict = collections.defaultdict(list)
        for connected_component in connected_components:
            for par in more_itertools.set_partitions(connected_component):
                par = [set(pa) for pa in par]
                partition_dict[frozenset(list(connected_component))].append(par)
        for combo in itertools.product(*partition_dict.values()):
            a_possible_median = []
            for lists in combo:
                for sets in lists:
                    a_possible_median.append(sets)
            possible_medians.append(a_possible_median)
        # now evaluate candidate genomes
        best_total = float("inf")
        best_max_branch_rate = float("inf") # use this to solve ties
        best_genome = []
        best_fusion_log = [[], [], []]
        best_fission_log = [[], [], []]
        for possible_median in possible_medians:
            possible_median = label_genome(tree_node, possible_median)
            best_genome, best_total, best_max_branch_rate, best_fusion_log, best_fission_log = \
            evaluate_genome_with_parsimony(branch_lengths, ios_1, ios_2, ios_3, possible_median, best_total, 
                best_max_branch_rate, best_genome, best_fusion_log, best_fission_log)
        return best_genome, [best_fusion_log[0], best_fusion_log[1]], [best_fission_log[0], best_fission_log[1]]

# permutations by fusion should only fuse LMSs found in the same connected component
def permute_genome_by_fusion(genome, connected_components):
    if len(genome) == len(connected_components):
        return genome
    else:
        temp_genome = copy.deepcopy(genome)
        sampled_chroms = random.sample(temp_genome, 2)
        successful_fusion = False
        for connected_component in connected_components:
            if list(sampled_chroms[0])[0] in connected_component and list(sampled_chroms[1])[0] in connected_component:
                successful_fusion = True
                for chrom in sampled_chroms:
                    temp_genome.remove(chrom)
                temp_genome.append(sampled_chroms[0].union(sampled_chroms[1]))
                return temp_genome
        if not successful_fusion:
            return genome

def permute_genome_by_fission(genome):
    if max([len(chrom) for chrom in genome]) == 1:
        return genome
    else:
        temp_genome = copy.deepcopy(genome)
        samplable_chroms = [chrom for chrom in temp_genome if len(chrom) > 1]
        sampled_chrom = random.sample(samplable_chroms, 1)
        temp_genome.remove(sampled_chrom[0])
        sampled_int = random.sample(range(1, len(sampled_chrom[0])), 1)[0]
        fission_product_1 = set()
        fission_product_2 = set()
        for index in random.sample(sampled_chrom[0], sampled_int):
            fission_product_1.add(index)
        for index in sampled_chrom[0]:
            if index not in fission_product_1:
                fission_product_2.add(index)
        temp_genome.append(fission_product_1)
        temp_genome.append(fission_product_2)
        return temp_genome

def heuristic_solver(tree_node, connected_components, ios_1, ios_2, ios_3, parsimony_genome, branch_lengths):
    print("[+] There are many possible genomes at this node so syngraph will undertake a heuristic search.", \
        "To avoid this, consider increasing the -m parameter.")
    print("[+] Searching ...")
    best_total = float("inf")
    best_max_branch_rate = float("inf")
    best_genome = []
    best_fusion_log = [[], [], []]
    best_fission_log = [[], [], []]
    starting_genome = copy.deepcopy(parsimony_genome)
    for rw_param in [[9, 1000], [7, 1000], [5, 1000], [3, 1000], [1, 1000]]:
        rw_length = rw_param[0]
        rw_iterations = rw_param[1]
        for iteration in range(0, rw_iterations):
            rw_genome = copy.deepcopy(starting_genome)
            if iteration == 0:
                pass
            else:
                for step in range(0, rw_length):
                    coin_flip = random.sample(["fusion", "fission"], 1)[0]
                    if coin_flip == "fusion":
                        rw_genome = permute_genome_by_fusion(rw_genome, connected_components)
                    elif coin_flip == "fission":
                        rw_genome = permute_genome_by_fission(rw_genome)
            rw_genome = label_genome(tree_node, rw_genome)
            # evaluate permuted genomes
            best_genome, best_total, best_max_branch_rate, best_fusion_log, best_fission_log = \
            evaluate_genome_with_parsimony(branch_lengths, ios_1, ios_2, ios_3, rw_genome, best_total, best_max_branch_rate, 
                best_genome, best_fusion_log, best_fission_log)
                #print(len(rw_genome), len(best_genome), best_total)
        starting_genome = [value for value in best_genome.values()]
    return best_genome, [best_fusion_log[0], best_fusion_log[1]], [best_fission_log[0], best_fission_log[1]]

def write_in_unassigned(tree_node, syngraph, taxa, LMSs, unassignable_markers):
    # if there are bugs in this function then this would be bad
    # how will this function behave when there are multiple traversals?
    reassigned_markers = 0
    LMS_sbt_chrom = {}
    for LMS in LMSs:
        l_sbt = set()
        for taxon in taxa:
            if taxon in syngraph.nodes[list(LMSs[LMS])[0]]['seqs_by_taxon']:
                l_sbt.add(syngraph.nodes[list(LMSs[LMS])[0]]['seqs_by_taxon'][taxon])
        l_chrom = syngraph.nodes[list(LMSs[LMS])[0]]['seqs_by_taxon'][tree_node]
        LMS_sbt_chrom[frozenset(l_sbt)] = l_chrom
    for unassignable in unassignable_markers:
        u_chroms = set()
        u_sbt = set()
        for taxon in taxa:
            if taxon in syngraph.nodes[unassignable]['seqs_by_taxon']:
                u_sbt.add(syngraph.nodes[unassignable]['seqs_by_taxon'][taxon])
        for l_sbt in LMS_sbt_chrom:
            matches = 0
            for seq in u_sbt:
                if seq in l_sbt:
                    matches += 1
            if matches == 2:
                u_chroms.add(LMS_sbt_chrom[l_sbt])
        if len(u_chroms) == 1:
            syngraph.nodes[unassignable]['taxa'].add(tree_node)
            syngraph.nodes[unassignable]['seqs_by_taxon'][tree_node] = list(u_chroms)[0]
            reassigned_markers += 1
    print("[=] Assigned {} of the markers not within any LMS".format(reassigned_markers))
    return syngraph

# name function variables, e.g. tree_node = None and tree_node = my_tree_node when you call

def median_genome(tree_node, working_syngraph, taxa, branch_lengths, minimum):
    # get LMSs >= minimum
    # LMSs < minimum are dealt with later and don't contribute to events
    # remaining LMSs represent a fissioned median genome
    LMSs, unassignable_markers = get_LMSs(working_syngraph, taxa, minimum)
    print("[=] Generated {} LMSs containing {} markers".format(len(LMSs.keys()), sum([len(LMSs[LMS]) for LMS in LMSs])))
    print("[=] A total of {} markers are not assigned to an LMS".format(len(unassignable_markers)))
    # given compact synteny of each extant genome to the fissioned median, generate connected components of LMSs
    # rewrite as instance_of_synteny
    ios_1 = compact_synteny_1(working_syngraph, LMSs, taxa[0])
    ios_2 = compact_synteny_1(working_syngraph, LMSs, taxa[1])
    ios_3 = compact_synteny_1(working_syngraph, LMSs, taxa[2])
    connected_components = generate_connected_components(ios_1, ios_2, ios_3)
    parsimony_genome = generate_parsimony_genome(ios_1, ios_2, ios_3, LMSs)
    print("[=] Generated {} connected components".format(len(connected_components)))
    # now check whether, given the connected components, the best genome can be found exhaustively
    # and then call the appropriate solving function
    bell_numbers = [1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975, 678570, 4213597]
    if max([len(connected_component) for connected_component in connected_components]) > 9:
        solved_connected_components, fusion_log, fission_log = heuristic_solver(tree_node, connected_components, ios_1, ios_2, 
            ios_3, parsimony_genome, branch_lengths)
    else:
        possible_medians = 1
        for connected_component in connected_components:
            possible_medians *= bell_numbers[len(connected_component)-1]
        if possible_medians > 50000:
            solved_connected_components, fusion_log, fission_log = heuristic_solver(tree_node, connected_components, ios_1, 
                ios_2, ios_3, parsimony_genome, branch_lengths)
        else:
            solved_connected_components, fusion_log, fission_log = exhaustive_solver(tree_node, connected_components, ios_1, 
                ios_2, ios_3, branch_lengths)
    print("[=] Found a median genome with {} chromosomes".format(len(solved_connected_components)))
    # from solved_connected_components, add this new ancestral genome to a syngraph
    working_syngraph.graph['taxa'].add(tree_node)
    for chrom in solved_connected_components:
        for LMS in solved_connected_components[chrom]:
            for graph_node_id in LMSs[LMS]:
                working_syngraph.nodes[graph_node_id]['taxa'].add(tree_node)
                working_syngraph.nodes[graph_node_id]['seqs_by_taxon'][tree_node] = chrom
    # finally, write in LMSs/markers that were too small or missing from a taxon, but can be assigned by parismony
    ### last part to fix
    ### then proof read
    working_syngraph = write_in_unassigned(tree_node, working_syngraph, taxa, LMSs, unassignable_markers)
    print("")
    return working_syngraph, fusion_log, fission_log

def tree_traversal(syngraph, params):
    # write a log
    log = ["#parent\tchild\tfusions\tfissions"]
    # copy syngraph
    traversal_0_syngraph = copy.deepcopy(syngraph)
    # define which taxa are extant and so can be sampled from the start
    available_taxa = set()
    for leaf in params.tree.get_leaves():
        available_taxa.add(leaf.name)
    print("[+] Starting first traversal ...")
    print("[+] ========================================================================")
    for tree_node in params.tree.traverse(strategy='postorder'):
        if not tree_node.is_leaf() and not tree_node.is_root():
            child_1 = tree_node.get_children()[0].name
            child_2 = tree_node.get_children()[1].name
            outgroup = get_closest_outgroup(params.tree, tree_node, available_taxa)
            branch_lengths = collections.defaultdict(list)
            for taxon in child_1, child_2, outgroup:
                up_distance, down_distance = get_branch_distances(params.tree, tree_node, taxon)
                branch_lengths["up"].append(up_distance)
                branch_lengths["down"].append(down_distance)
            print("[+] Inferring median genome for {} using data from {}, {}, and {} ...". format(tree_node.name, child_1, 
                child_2, outgroup))
            traversal_0_syngraph, fusion_log, fission_log = median_genome(tree_node.name, traversal_0_syngraph, [child_1, child_2, 
                outgroup], branch_lengths, params.minimum)
            available_taxa.add(tree_node.name)
            for i in 0, 1:
                log.append("{}\t{}\t{}\t{}".format(tree_node.name, [child_1, child_2][i], fusion_log[i], fission_log[i]))
    print("[=] ========================================================================")
    return traversal_0_syngraph, log


#############################################################################################
#############################################################################################
#############################################################################################

class Syngraph(nx.Graph):

    def __init__(self, name='', **attr):
        nx.Graph.__init__(self, name=name, taxa=set(), **attr)
        
    def __repr__(self):
       return "Syngraph(name=%r, taxa=%r, ...)" % (self.name, self.graph['taxa']) 

    def __eq__(self, other):
        if isinstance(other, Syngraph):
            return (self.nodes == other.nodes and 
                    self.edges == other.edges and 
                    self.graph == other.graph)
        return False

    def save(self, parameterObj, check_consistency=True):
        graph_file = '%s.pickle' % parameterObj.outprefix
        nx.write_gpickle(self, graph_file)
        if check_consistency:
            syngraph_loaded = Syngraph()
            syngraph_loaded.from_file(graph_file)
            # syngraph_loaded.add_node('1')             # this will fuck up saving/loading consistency
            if self == syngraph_loaded:
                return graph_file
            else:
                print("[-] Consistency check failed ...")
                print("     [Nodes]", self.nodes == syngraph_loaded.nodes)
                print("     [Edges]", self.edges == syngraph_loaded.edges)
                print("     [Graph]", self.graph == syngraph_loaded.graph)
                sys.exit("[X] Saving/Loading consistency is compromised.")
        else:
            return graph_file

    def get_taxon_syngraph(self, taxon=None):
        if not taxon is None:
            edges = []
            for u, v, taxa in self.edges(data='taxa'):
                if taxon in set(taxa):
                    edges.append((u, v, {'taxa': set([taxon])}))
            syngraph = Syngraph()
            syngraph.from_edges(edges, taxa=set([taxon]))
            for graph_node_id in self.nodes:
                if graph_node_id not in list(syngraph):
                    syngraph.add_node(graph_node_id)
            return syngraph

    def from_file(self, graph_file):
        g = nx.read_gpickle(graph_file)
        self.add_nodes_from(g.nodes(data=True))
        self.add_edges_from(g.edges(data=True))
        self.graph = g.graph

    def from_edges(self, edges=[], taxa=set()):
        '''
        edges: is list of 3-tuples (u,v, {'taxa': taxa})
        taxa: is set of all taxa 
        '''
        self.add_edges_from(edges)
        self.graph['taxa'] = set(taxa)

    def from_markerObjs(self, markerObjs):
        prev_markerObj = MarkerObj(None)
        marker_ids_by_seq_id_by_taxon = collections.defaultdict(functools.partial(collections.defaultdict, list))
        for markerObj in markerObjs:
            marker_ids_by_seq_id_by_taxon[markerObj.taxon][markerObj.seq].append(markerObj.name)
            self.graph['taxa'].add(markerObj.taxon)                                   
            if not markerObj.name in self:
                # add node if it does not exist yet                                              
                self.add_node(markerObj.name, taxa=set(), terminal=set(), seqs_by_taxon={})             
            # add taxon to node (useful if '--missing')
            self.nodes[markerObj.name]['taxa'].add(markerObj.taxon)
            self.nodes[markerObj.name]['seqs_by_taxon'][markerObj.taxon] = markerObj.seq
            # check whether same seq/taxon
            if markerObj.is_syntenic(prev_markerObj):
                # add edge if it does not exist
                if not self.has_edge(prev_markerObj.name, markerObj.name):
                    self.add_edge(prev_markerObj.name, markerObj.name, taxa={})
                # add distance to edge
                self[prev_markerObj.name][markerObj.name]['taxa'][markerObj.taxon] = prev_markerObj.distance(markerObj)    
            else:
                # first marker on seq
                if not prev_markerObj.name is None:
                    # deal with terminal marker on prev seq
                    self.nodes[prev_markerObj.name]['terminal'].add(prev_markerObj.taxon)
                self.nodes[markerObj.name]['terminal'].add(markerObj.taxon)
            prev_markerObj = markerObj
        self.nodes[markerObj.name]['terminal'].add(markerObj.taxon)
        self.graph['marker_ids_by_seq_id_by_taxon'] = marker_ids_by_seq_id_by_taxon

    def get_target_edge_sets_by_taxon(self, graph_node_id):
        target_edge_sets_by_taxon = collections.defaultdict(set)
        for u, v, taxa in self.edges([graph_node_id], data='taxa'):
            for taxon in taxa:
                target_edge_sets_by_taxon[taxon].add(frozenset((u, v)))
        for taxon, target_edge_set in target_edge_sets_by_taxon.items():
            if len(target_edge_set) == 1:
                target_edge_sets_by_taxon[taxon].add(frozenset((graph_node_id, '%s.terminal' % (graph_node_id))))
        for taxon in self.nodes[graph_node_id]['seqs_by_taxon']:
            if taxon not in target_edge_sets_by_taxon:
                target_edge_sets_by_taxon[taxon].add(frozenset((graph_node_id, '%s.terminal' % (graph_node_id))))
                target_edge_sets_by_taxon[taxon].add(frozenset((graph_node_id, '%s.also_terminal' % (graph_node_id))))
        return target_edge_sets_by_taxon

    def show_metrics(self):
        taxon_count = len(self.graph['taxa'])
        node_total_count = nx.number_of_nodes(self)
        edge_total_count = nx.number_of_edges(self)
        connected_component_count = nx.number_connected_components(self)
        # taxa_per_edge = {}
        # edges_per_node = {}
        # for i in range(1, len(self.graph['taxa'])+1):
        #     taxa_per_edge[i] = sum([1 for (_, __, taxa) in self.edges.data('taxa') if len(taxa) == i])
        # for j in range(1, (len(self.graph['taxa'])+1)*2):
        #     edges_per_node[j] = 0
        # for graph_node_id in self.nodes:
        #     neighbours = self.degree(graph_node_id)
        #     edges_per_node[neighbours] += 1
        print("[=] ====================================")
        print("[=] Taxa = %s" % taxon_count)
        print("[=] Nodes (Markers) = %s" % node_total_count)
        print("[=] Distinct Edges (Adjacencies) = %s" % edge_total_count) 
        # print("[=]   distinct = %s" % edge_total_count)
        # print("[=]   taxa per edge\tcount")
        # for key in taxa_per_edge:
        #     print("[=]   {}\t{}".format(key, taxa_per_edge[key]))
        # print("[=]   edges per node\tcount")
        # for key in edges_per_node:
        #     print("[=]   {}\t{}".format(key, edges_per_node[key]))
        print("[=] Subgraphs (connected components) = %s" % connected_component_count)
        print("[=] ====================================")

    # maybe this could be combined with the above?
    def show_recon_metrics(self, gene_order, name):
        connected_component_count = nx.number_connected_components(self)
        cc_lengths = []
        for component in nx.connected_components(self):
            cc_lengths.append(len(component))
        cc_lengths.sort()
        if gene_order == True:
            edges_per_node = {}
            for j in range(1, 99):
                edges_per_node[j] = 0
            for graph_node_id in self.nodes:
                neighbours = self.degree(graph_node_id)
                edges_per_node[neighbours] += 1
            resolved = 0
            unresolved = 0
            for key in edges_per_node:
                if int(key) > 2:
                    unresolved += edges_per_node[key]
                else:
                    resolved += edges_per_node[key]
        print("[=] ====================================")
        print(name)
        if gene_order == True:
            print("[=] Resolved edges = %s" % (resolved / (resolved + unresolved)))
        print("[=] Subgraphs (connected components) = %s" % connected_component_count)
        print("[=] Subgraph lengths = %s" % cc_lengths)
        print("[=] ====================================")

    def taxon_jaccard(self):
        """
        for a syngraph
            - extract taxon graphs
            - do pairwise comparisons between taxon graphs
            - print jaccard index of edges
        """
        taxon_graphs = {}
        taxon_jaccard = collections.defaultdict(dict)
        for taxon in self.graph['taxa']:
            taxon_graphs[taxon] = self.get_taxon_syngraph(taxon=taxon)
        for taxon_A in taxon_graphs:
            for taxon_B in taxon_graphs:
                intersection_L = len(set(taxon_graphs[taxon_A].edges()).intersection(set(taxon_graphs[taxon_B].edges())))
                union_L = len(set(taxon_graphs[taxon_A].edges()).union(set(taxon_graphs[taxon_B].edges())))
                jaccard = intersection_L / union_L
                taxon_jaccard[taxon_A][taxon_B] = jaccard
        df = pd.DataFrame.from_dict(taxon_jaccard)
        print(df.round(decimals=3))
       
    def plot(self, outprefix, cmap='Set2', as_multigraph=True):
        taxon_count = len(self.graph['taxa'])
        colour_by_taxon = get_hex_colours_by_taxon(self.graph['taxa'], cmap=cmap)
        if as_multigraph:
            multi_graph = nx.MultiGraph() # create new multigraph
            for node, data in self.nodes.data(True):
                multi_graph.add_node(node, label=node, color='black', width=0.5, shape='rect')
                if 'terminal' in data:
                    for taxon in data['terminal']:
                        new_terminal_node = "%s_%s_terminal" % (taxon, node)
                        multi_graph.add_node(new_terminal_node, color=colour_by_taxon[taxon], label='', width=0.5, shape='rect')
                        multi_graph.add_edge(new_terminal_node, node, penwidth=0.5)
            for u, v, taxa in self.edges.data('taxa'):
                for taxon in taxa:
                    multi_graph.add_edge(u, v, color=colour_by_taxon[taxon], style='solid', penwidth=4)
            plot_graph = multi_graph
        else:
            plot_graph = nx.Graph()
            for node, data in self.nodes.data(True):
                plot_graph.add_node(node, label='', width=0.5, shape='point')
                if 'terminal' in data:
                    for taxon in data['terminal']:
                        new_terminal_node = "%s_%s_terminal" % (taxon, node) 
                        plot_graph.add_node(new_terminal_node, color=colour_by_taxon[taxon], label='', width=0.5, shape='point')
                        plot_graph.add_edge(new_terminal_node, node, penwidth=0.5)
            for u, v, taxa in self.edges.data('taxa'):
                edge_style = 'solid'
                if len(taxa) == 1:
                    edge_style = 'dashed'
                plot_graph.add_edge(u, v, color='grey', style=edge_style, penwidth=4+(len(taxa)/taxon_count))
        # add legend
        for taxon, colour in colour_by_taxon.items():
            plot_graph.add_node(taxon, color=colour_by_taxon[taxon], label=taxon, width=0.5, style='filled', shape='circle')
        # layout
        if nx.number_of_nodes(self) < 50:
            layout='dot'
        else:
            layout='neato'
        # generate graphviz graph
        G = nx.drawing.nx_agraph.to_agraph(plot_graph)
        out_f = "%s.pdf" % outprefix

        G.draw(out_f, prog=layout, args=" -Grankdir=BT")
        return out_f

class MarkerObj():
    def __init__(self, name=None, desc=None, status=None, taxon=None, seq=None, coord=float("nan")):
        self.name = name
        self.desc = desc if desc is not None else name
        self.status = status
        self.taxon = taxon
        self.seq = seq
        self.coord = coord

    def __repr__(self):
        return "MarkerObj(name=%r, desc=%r, status=%r, taxon=%r, seq=%r, coord=%f)" % (
            self.name, self.desc, self.status, self.taxon, self.seq, self.coord) 

    def __eq__(self, other):
        if isinstance(other, MarkerObj):
            return (self.name == other.name)
        return False

    def is_syntenic(self, other):
        if isinstance(other, MarkerObj):
            if self.taxon == other.taxon and self.seq == other.seq:
                return True
        return False

    def distance(self, other):
        if self.is_syntenic(other):
            if self.coord == other.coord:
                return 0.0
            else:
                return float(max(self.coord, other.coord) - min(self.coord, other.coord))                
        return float("nan")