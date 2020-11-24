import sys
import itertools
import matplotlib
import matplotlib.cm as cm
import pandas as pd
from tqdm import tqdm
import networkx as nx
import collections
pd.set_option('display.max_rows', None)
import functools
import statistics
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
                                markerObj= MarkerObj(name=_name, status=status, taxon=taxon, seq=seq, coord=_coord)
                                taxon_markerObjs.append(markerObj)
                        else:
                            markerObj = MarkerObj(name=name, status=status, taxon=taxon, seq=seq, coord=start)
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
                        markerObj = MarkerObj(name=marker, desc=locus, status=status, taxon=taxon, seq=seq, coord=start)
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

def fitch(states_by_taxon, number_of_states, tree): # states_by_taxon_node should be a dict, with keys as taxon and states as sets
    states_by_tree_node = {}
    for tree_node in tree.traverse(strategy='postorder'):
        if tree_node.name in states_by_taxon:
            states_by_tree_node[tree_node.name] = states_by_taxon[tree_node.name]
        elif not tree_node.is_leaf():
            intersection = set.intersection(*[states_by_tree_node.get(child_node.name, {False}) for child_node in tree_node.get_children()])
            if len(intersection) >= number_of_states:
                states_by_tree_node[tree_node.name] = intersection
            else:
                states_by_tree_node[tree_node.name] = set.union(
                    *[states_by_tree_node.get(child_node.name, {False}) for child_node in tree_node.get_children()])
        else:
            pass
    for tree_node in tree.traverse(strategy='levelorder'):
        if not tree_node.is_root():
            parent_tree_node = tree_node.up
            intersection = states_by_tree_node.get(parent_tree_node.name, {False}).intersection(states_by_tree_node.get(tree_node.name, {False}))
            if len(intersection) >= number_of_states:
                states_by_tree_node[tree_node.name] = intersection
    return(states_by_tree_node)

def reconstruct_syngraphs_by_tree_node(syngraph, tree, algorithm='fitch'):
    '''
    - input: syngraph, tree
    - output: novel graphs with fitch edges for each internal tree node
    '''
    if algorithm == 'fitch':
        edges_by_tree_node_by_graph_node = collections.defaultdict(dict) # nested dict, graph_node --> taxon --> edges
        edges_by_tree_node = collections.defaultdict(list)
        taxa_by_tree_node = collections.defaultdict(set)
        for graph_node_id in syngraph.nodes:
            edge_sets_by_taxon = syngraph.get_target_edge_sets_by_taxon(graph_node_id)
            #print('edge_sets_by_taxon', edge_sets_by_taxon)
            edges_by_tree_node_by_graph_node[graph_node_id] = fitch(edge_sets_by_taxon, 2, tree)
        for graph_node_id, _edges_by_tree_node in edges_by_tree_node_by_graph_node.items():
            for tree_node, edges in _edges_by_tree_node.items():
                #print("tree_node", tree_node)
                #print("edges", edges)
                for (u, v) in edges:
                    taxa_under_node = set([node.name for node in (tree&tree_node).iter_leaves()])
                    taxa_by_tree_node[tree_node].update(taxa_under_node)
                    edge_taxa = []
                    for taxon in taxa_under_node:
                        if frozenset([u, v]) in _edges_by_tree_node[taxon]:
                            edge_taxa.append(taxon)
                    edges_by_tree_node[tree_node].append((u, v, {'taxa': edge_taxa}))
        syngraph_by_tree_node = {}
        for tree_node, edges in edges_by_tree_node.items():
            syngraph_by_tree_node[tree_node] = Syngraph()
            syngraph_by_tree_node[tree_node].from_edges(edges, taxa=taxa_by_tree_node[tree_node])
            syngraph_by_tree_node[tree_node].show_recon_metrics(True, tree_node)
            #syngraph_by_tree_node[tree_node].plot(outprefix="node_%s" % tree_node)
    return syngraph_by_tree_node

def plot_histogram(x, out_f):
    fig, ax = plt.subplots(figsize=(14, 5))
    hist, bins = np.histogram(x, bins=50)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    ax.bar(center, hist, align='center', width=width)
    fig.savefig('%s.png' % out_f, format="png")

def find_LinkedMarkerSets(syngraph):
    '''
    input: syngraph
    output: linked marker sets
    '''
    LinkedMarkerSets = [] # list of lists
    for graph_node_id in syngraph.nodes():
        if len(LinkedMarkerSets) == 0:
            LinkedMarkerSets.append([graph_node_id])
        else:
            match = False
            for LMS in LinkedMarkerSets:
                if syngraph.nodes[graph_node_id]['seqs_by_taxon'] == syngraph.nodes[LMS[0]]['seqs_by_taxon']:
                    LMS.append(graph_node_id)
                    match = True
            if match == False:
                LinkedMarkerSets.append([graph_node_id])
    LinkedMarkerSets.sort(reverse=True, key=len)
    # make an ID dict
    LinkedMarkerSet_IDs = {} # dict with a ID key for each FSBlock
    ID = 0
    for LMS in LinkedMarkerSets:
        ID += 1
        LinkedMarkerSet_IDs[ID] = LMS
    Total_IDs = ID
    print("[=] Found {} linked marker sets".format(Total_IDs))
    print("[=] Largest linked marker set has a length of {}".format(len(LinkedMarkerSets[0])))
    print("[=] Smallest linked marker set has a length of {}".format(len(LinkedMarkerSets[-1])))
    for LMS in LinkedMarkerSets:
        print(len(LMS))
    return LinkedMarkerSets, LinkedMarkerSet_IDs


def reconstruct_linkage_groups_for_each_tree_node(syngraph, linked_marker_sets, linked_marker_set_IDs, tree, algorithm='fitch'):

    '''
    - input: syngraph, tree
    - output: for each tree node, nodes assigned to LGs
    ''' 
    LG_store = collections.defaultdict(lambda: collections.defaultdict(dict)) # put inferred LGs in here later
    # calculate coocurence matrix for each taxon and produce a graph of the matrix, a C_graph
    coocurence_by_taxon_by_LMS_ID = collections.defaultdict(lambda: collections.defaultdict(lambda: {False})) # nested dict, set_of_2_LMS --> taxon --> True/False
    if algorithm == 'fitch':
        edges_by_tree_node = collections.defaultdict(list)
        print("[+] Collecting coocurence data from taxa")
        for taxon in syngraph.graph['taxa']:
            for ID_comparison in itertools.combinations(range(1, (len(linked_marker_sets)+1)), 2):
                a_graph_node = linked_marker_set_IDs[ID_comparison[0]][0]
                another_graph_node = linked_marker_set_IDs[ID_comparison[1]][0]
                if syngraph.nodes[a_graph_node]['seqs_by_taxon'][taxon] == syngraph.nodes[another_graph_node]['seqs_by_taxon'][taxon]:
                    coocurence_by_taxon_by_LMS_ID[frozenset(ID_comparison)][taxon] = {True}
        # reconstruct coocurence matrix for internal nodes
        print("[+] Estimating coocurence data for ancestral genomes")
        for ID_comparison in coocurence_by_taxon_by_LMS_ID:
            coocurence_by_taxon_by_LMS_ID[ID_comparison] = fitch(coocurence_by_taxon_by_LMS_ID[ID_comparison], 1, tree)
            u, v = tuple(ID_comparison)
            for tree_node, coocurence in coocurence_by_taxon_by_LMS_ID[ID_comparison].items():
                if coocurence == {True}:
                    edges_by_tree_node[tree_node].append((u, v, {'taxa': tree_node}))
        # construct a graph where edges represent coocurence
        print("[+] Building coocurence graphs")
        C_graphs_by_tree_node = {}
        for tree_node, edges in edges_by_tree_node.items():
            C_graphs_by_tree_node[tree_node] = Syngraph()
            C_graphs_by_tree_node[tree_node].from_edges(edges, taxa=set(tree_node))
            # add LMS that connect to nothing as island nodes
            for graph_node_id in range(1, len(linked_marker_sets)+1):
                if graph_node_id not in C_graphs_by_tree_node[tree_node].nodes():
                    C_graphs_by_tree_node[tree_node].add_node(graph_node_id)
    # apply heuristic algorithm for pruning out LGs from C
    print("[+] Building linkage group graphs")
    LG_graphs_by_tree_node = {}
    for tree_node in C_graphs_by_tree_node:
        LG_graphs_by_tree_node[tree_node] = nx.DiGraph()
        # iterate over nodes in C graph, make a directed edge for each node's greatest size connection
        # edges are only made from a smaller to a greater sized node, or between equally sized nodes
        # if there is a tie in size, then multiple edges will be made
        for graph_node_id in C_graphs_by_tree_node[tree_node]:
            neighbours = list(C_graphs_by_tree_node[tree_node].neighbors(graph_node_id))
            biggest_neighbour = []
            if len(neighbours) > 0:
                for neighbour in neighbours:
                    if len(linked_marker_set_IDs[neighbour]) >= len(linked_marker_set_IDs[graph_node_id]):
                        if len(biggest_neighbour) == 0:
                            biggest_neighbour = [neighbour]
                        else:
                            if len(linked_marker_set_IDs[neighbour]) == len(linked_marker_set_IDs[biggest_neighbour[0]]):
                                biggest_neighbour.append(neighbour)
                            elif len(linked_marker_set_IDs[neighbour]) > len(linked_marker_set_IDs[biggest_neighbour[0]]):
                                biggest_neighbour = [neighbour]
                if len(biggest_neighbour) > 0:
                    for big_neighbour in biggest_neighbour:
                        LG_graphs_by_tree_node[tree_node].add_edge(graph_node_id, big_neighbour)
                else:
                    LG_graphs_by_tree_node[tree_node].add_node(graph_node_id)
            else:
                LG_graphs_by_tree_node[tree_node].add_node(graph_node_id)
        # use destinations in the graph to assign nodes to LGs
        # a destination is either a node with no successors/neighbors of greater size
        LG_store[tree_node][frozenset(['unassignables'])] = []
        destinations = []
        for ID in range(1, len(linked_marker_sets)+1):
            if len(linked_marker_set_IDs[ID]) > 1: # here make a rule that destinations must be 
                is_destination = True
                for descendant in nx.descendants(LG_graphs_by_tree_node[tree_node], ID):
                    if len(linked_marker_set_IDs[descendant]) > len(linked_marker_set_IDs[ID]):
                        is_destination = False
                if is_destination == True:
                    if len(destinations) == 0:
                        destinations.append([ID])
                    else:
                        in_cycle = False
                        for descendant in nx.descendants(LG_graphs_by_tree_node[tree_node], ID):
                            for destination in destinations:
                                if descendant in destination:
                                    destination.append(ID)
                                    in_cycle = True
                        if in_cycle == False:
                            destinations.append([ID])
        # destinations (terminals and cycles) have now been defined
        # loop through IDs to find their destination
        # if they are a destination then assign them to that LG
        # if they have 1 destination then they should be assigned to that destination's LG
        # if they have >1 destinations then they cannot be assigned
        for destination in destinations:
            LG_store[tree_node][frozenset(destination)] = []
        for ID in range(1, len(linked_marker_sets)+1):
            is_destination = False
            for destination in LG_store[tree_node]:
                if ID in destination:
                    LG_store[tree_node][destination].append(ID)
                    is_destination = True
            if is_destination == False:
                descendants = nx.descendants(LG_graphs_by_tree_node[tree_node], ID)
                destination_of_ID = set()
                for descendant in descendants:
                    for destination in LG_store[tree_node]:
                        if descendant in destination:
                            destination_of_ID.add(destination)
                if len(destination_of_ID) == 1:
                    for dest in destination_of_ID:
                        LG_store[tree_node][dest].append(ID)
                else:
                    LG_store[tree_node][frozenset(['unassignables'])].append(ID)
    # all done, return the stored LGs
    return LG_store

def analyse_chromosomal_units(syngraph, LG_store, linked_marker_set_IDs, mrca):
    ''''
    input : LG store, linked_marker_set_IDs, a target node for defining units
    output : tsv for each species - chrom - which unit(s) they contain, plus information about units
    '''
    print("[+] Using the following internal node to define chromosomal units \n{}\n".format(mrca.get_ascii(show_internal=True)))
    # extract units using specified mrca
    chromosomal_units = {} # key is unit name, e.g. unit_1 ordered by size, value is list of markers
    temp_units = {} # these will be stored here with the wrong names
    temp_unit_names = 0
    for LG in LG_store[mrca.name]:
        if LG != frozenset(['unassignables']):
            temp_unit_names += 1
            temp_units[str(temp_unit_names)] = []
            for ID in LG_store[mrca.name][LG]:
                for marker in linked_marker_set_IDs[ID]:
                    temp_units[str(temp_unit_names)].append(marker)
    # now give names based on length and store in chromosomal_units
    chromosomal_units_names = 0
    for unit in sorted(temp_units, key=lambda unit: len(temp_units[unit]), reverse=True):
        chromosomal_units_names += 1
        chromosomal_units["unit_"+str(chromosomal_units_names)] = temp_units[unit]
    print("[=] Estimated {} ancestral chromosome units".format(chromosomal_units_names))
    print("[=] {} markers could not be assigned to any unit".format(len(LG_store[mrca.name][frozenset(['unassignables'])])))
    print("[+] Calculating and then outputting how extant chromosomes relate to ancestral units\n")
    for tree_node in LG_store:
        if tree_node in syngraph.graph['taxa']:
            temporary_results_dict = collections.defaultdict(lambda: collections.defaultdict(int))
            for LG in LG_store[tree_node]:
                for ID in LG_store[tree_node][LG]:
                    for marker in linked_marker_set_IDs[ID]:
                        for unit in chromosomal_units:
                            if marker in chromosomal_units[unit]:
                                #print("{}\t{}\t{}\t{}".format(tree_node, marker, syngraph.nodes[marker]['seqs_by_taxon'][tree_node], unit))
                                taxon_chromo = syngraph.nodes[marker]['seqs_by_taxon'][tree_node]
                                temporary_results_dict[taxon_chromo][unit] += 1
            for chromo in sorted(temporary_results_dict):
                results_string = "NA"
                for unit in sorted(temporary_results_dict[chromo]):
                    percent_found = round(((temporary_results_dict[chromo][unit] / len(chromosomal_units[unit])) * 100), 3)
                    if results_string == "NA":
                        results_string = str(percent_found) + "%_" + unit
                    else:
                        results_string = results_string + "\t" + str(percent_found) + "%_" + unit
                print("{}\t{}\t{}".format(tree_node, chromo, results_string))


def get_hex_colours_by_taxon(taxa, cmap='Spectral'):
    return {taxon: matplotlib.colors.rgb2hex(cm.get_cmap(cmap, len(taxa))(i)[:3]) for i, taxon in enumerate(sorted(taxa))}

#############################################################################################
################################################################################:)###########
### below are functions for implementing fusion+fission models ##############################
################################################################################:(###########
#############################################################################################

def compact_synteny(syngraph, ref_taxon, query_taxon, minimum, mode):
    querychrom2index = collections.defaultdict(set)
    if mode == "LMS_syngraph":
        for graph_node_id in syngraph.nodes():
            for LMS in ref_taxon:
                if graph_node_id in ref_taxon[LMS]:
                    querychrom2index[syngraph.nodes[graph_node_id]['seqs_by_taxon'][query_taxon]].add(LMS)
    elif mode == "index_index": # list of sets
        for querychrom in query_taxon:
            for index in querychrom:
                for refchrom in ref_taxon:
                    if index in refchrom:
                        querychrom2index[frozenset(querychrom)].add(frozenset(refchrom))
    return(querychrom2index)

def ffsd(instance_of_synteny):
    def check_for_fusions(instance_of_synteny, fusions_so_far):
        fusions = fusions_so_far
        new_fusions = 0
        for combo in itertools.combinations(instance_of_synteny.keys(), 2):
            if instance_of_synteny[combo[0]].intersection(instance_of_synteny[combo[1]]):
                instance_of_synteny[combo[0]] = instance_of_synteny[combo[0]].union(instance_of_synteny[combo[1]])
                instance_of_synteny[combo[1]] = set()
                fusions += 1
                new_fusions += 1
        if new_fusions > 0:
            return check_for_fusions(instance_of_synteny, fusions)
        else:
            return(fusions)
    def check_for_fissions(instance_of_synteny, fissions_so_far):
        fissions = fissions_so_far
        for querychrom in instance_of_synteny:
            indices = len(instance_of_synteny[querychrom])
            if indices > 1:
                fissions += (indices - 1)
        return(fissions)
    total_fusions = check_for_fusions(instance_of_synteny, 0)
    total_fissions = check_for_fissions(instance_of_synteny, 0)
    return(total_fusions, total_fissions)

def get_LMS_triplets(syngraph, ref_taxon, query_taxon, query2_taxon, minimum):
    LinkedMarkerSets = collections.defaultdict(set)
    for graph_node_id in syngraph.nodes():
        triplet_seqs_by_taxon = set()
        for taxon in ref_taxon, query_taxon, query2_taxon:
            if taxon in syngraph.nodes[graph_node_id]['seqs_by_taxon']:
                triplet_seqs_by_taxon.add(syngraph.nodes[graph_node_id]['seqs_by_taxon'][taxon])
        # are seqs_by_taxon unique to each taxon? If not the above will not work.
        triplet_seqs_by_taxon = frozenset(list(triplet_seqs_by_taxon))
        # need to think about this function given missingness, below is a temp fix
        if len(triplet_seqs_by_taxon) == 3:
            LinkedMarkerSets[triplet_seqs_by_taxon].add(graph_node_id)
    Filtered_LinkedMarkerSets = {}
    Filtered_LMS_count = 1
    for LMS in LinkedMarkerSets:
        if len(LinkedMarkerSets[LMS]) >= minimum:
            Filtered_LinkedMarkerSets["LMS_" + str(Filtered_LMS_count)] = LinkedMarkerSets[LMS]
            Filtered_LMS_count += 1
    return(Filtered_LinkedMarkerSets)

def generate_connected_components(ios_ref, ios_query, ios_query2, iteration, connected_components):
    if iteration == 0:
        for ios in ios_ref, ios_query, ios_query2:
            for chrom in ios.values():
                connected_components.append(chrom)
    another_iteration = False
    for combo in itertools.combinations(range(0, len(connected_components)), 2):
        if connected_components[combo[0]].intersection(connected_components[combo[1]]):
            connected_components[combo[0]] = connected_components[combo[0]].union(connected_components[combo[1]])
            connected_components[combo[1]] = set()
            another_iteration = True
    if another_iteration == True:
        return(generate_connected_components(ios_ref, ios_query, ios_query2, iteration+1, connected_components))
    else:
        connected_components = [component for component in connected_components if component != set()]
        return(connected_components)

def solve_connected_component(connected_component, ios_ref, ios_query, ios_query2):
    if len(connected_component) == 1:
        connected_component = [connected_component]
        return(0, 0, connected_component)
    else:
        best_fissions = float("inf")
        best_fusions = float("inf")
        best_total = float("inf")
        best_partition = []
        for ios in ios_ref, ios_query, ios_query2:
            temp_chrom = []
            for chrom in ios.values():
                if len(set(connected_component).intersection(chrom)) > 0:
                    temp_chrom.append(chrom)
            if ios == ios_ref:
                ios_ref_f = temp_chrom
            if ios == ios_query:
                ios_query_f = temp_chrom
            if ios == ios_query2:
                ios_query2_f = temp_chrom

        if len(connected_component) > 12:
            print("[+] Too many (>12) LMSs are in a connected component so syngraph will undertake a heuristic search. To avoid this try increasing the -m parameter.")
            print("[+] Searching genome land ...")
            def permute_par_by_fusion(par):
                temp_par = copy.deepcopy(par)
                if len(temp_par) == 1:
                    return(temp_par)
                else:
                    sampled_chroms = random.sample(temp_par, 2)
                    for chrom in sampled_chroms:
                        temp_par.remove(chrom)
                    temp_par.append(sampled_chroms[0].union(sampled_chroms[1]))
                    return(temp_par)
            def permute_par_by_fission(par):
                temp_par = copy.deepcopy(par)
                if max([len(chrom) for chrom in temp_par]) == 1:
                    return(temp_par)
                else:
                    sampled_chrom = random.sample(temp_par, 1)
                    if len(sampled_chrom[0]) < 2:
                        return permute_par_by_fission(temp_par)
                    else:
                        temp_par.remove(sampled_chrom[0])
                        sampled_int = random.sample(range(1, len(sampled_chrom[0])), 1)[0]
                        fission_product_1 = set()
                        fission_product_2 = set()
                        for index in range(0, sampled_int):
                            fission_product_1.add(random.sample(sampled_chrom[0], 1)[0])
                        for index in sampled_chrom[0]:
                            if index not in fission_product_1:
                                fission_product_2.add(index)
                        temp_par.append(fission_product_1)
                        temp_par.append(fission_product_2)
                        return(temp_par)
            for starting_par in ios_ref_f, ios_query_f, ios_query2_f:
                # five random walk rounds
                # each starts from the previous round's best point
                # a round consists of N walks
                # the length of the walk is reduced each round
                for rw_param in [[25, 500], [10, 400], [5, 300], [3, 200], [1, 100]]:
                    rw_length = rw_param[0]
                    rw_iterations = rw_param[1]
                    for iteration in range(0, rw_iterations):
                        rw_par = copy.deepcopy(starting_par)
                        for step in range(0, rw_length):
                            coin_flip = random.sample(["fusion", "fission"], 1)[0]
                            if coin_flip == "fusion":
                                rw_par = permute_par_by_fusion(rw_par)
                            elif coin_flip == "fission":
                                rw_par = permute_par_by_fission(rw_par)
                        total_fusions = 0
                        total_fissions = 0
                        for ios in ios_ref_f, ios_query_f, ios_query2_f:
                            fusions, fissions = ffsd(compact_synteny("null", ios, rw_par, "null", "index_index"))
                            total_fusions += fusions
                            total_fissions += fissions
                        if total_fusions + total_fissions < best_total:
                            best_fissions = total_fissions
                            best_fusions = total_fusions
                            best_total = total_fusions + total_fissions
                            best_partition = rw_par                   
                    starting_par = best_partition
        else:
            # function is from https://stackoverflow.com/questions/19368375/set-partitions-in-python
            # this will exhaustively generate all possible partitions
            def partition(cc):
                if len(cc) == 1:
                    yield [cc]
                    return
                first = cc[0]
                for smaller in partition(cc[1:]):
                    for n, subset in enumerate(smaller):
                        yield smaller[:n] + [[first] + subset]  + smaller[n+1:] 
                    yield [[first]] + smaller
            connected_component_listed = list(connected_component)
            for par in partition(connected_component_listed):
                par = [set(index) for index in par]
                total_fusions = 0
                total_fissions = 0
                for ios in ios_ref_f, ios_query_f, ios_query2_f:
                    fusions, fissions = ffsd(compact_synteny("null", ios, par, "null", "index_index"))
                    total_fusions += fusions
                    total_fissions += fissions
                if total_fusions + total_fissions < best_total:
                    best_fissions = total_fissions
                    best_fusions = total_fusions
                    best_total = total_fusions + total_fissions
                    best_partition = par        
        return(best_fissions, best_fusions, best_partition)


def median_genome(tree_node, syngraph, ref_taxon, query_taxon, query2_taxon, minimum):
    # get LMSs >= minimum
    # LMSs < minimum are dealt with later and don't contribute to events
    # remaining LMSs represent a fissioned median genome
    LMSs = get_LMS_triplets(syngraph, ref_taxon, query_taxon, query2_taxon, minimum)
    print("[=] Generated {} LMSs containing {} markers".format(len(LMSs.keys()), sum([len(LMSs[LMS]) for LMS in LMSs])))
    # given compact synteny of each extant genome to the fissioned median, generate connected components of LMSs
    ios_ref = compact_synteny(syngraph, LMSs, ref_taxon, minimum, "LMS_syngraph")
    ios_query = compact_synteny(syngraph, LMSs, query_taxon, minimum, "LMS_syngraph")
    ios_query2 = compact_synteny(syngraph, LMSs, query2_taxon, minimum, "LMS_syngraph")
    connected_components = []
    solved_connected_components = []
    total_fissions = 0
    total_fusions = 0
    connected_components = generate_connected_components(ios_ref, ios_query, ios_query2, 0, connected_components)
    print("[=] Generated {} connected components".format(len(connected_components)))
    # for each connected component return the configuration which minimises ffsd
    for connected_component in connected_components:
        fissions, fusions, solved_connected_component = solve_connected_component(connected_component, ios_ref, ios_query, ios_query2)
        total_fissions += fissions
        total_fusions += fusions
        [solved_connected_components.append(chrom) for chrom in solved_connected_component]
    print("[=] Found a median genome with {} chromosomes that only requires:".format(len(solved_connected_components)))
    print("[=]\t{}\tfissions".format(total_fissions))
    print("[=]\t{}\tfusions".format(total_fusions))
    # from solved_connected_components, add this new ancestral genome to the syngraph
    # also write in LMSs that were too small
    syngraph.graph['taxa'].add(tree_node)
    for LMS in LMSs:
        for graph_node_id in LMSs[LMS]:
            syngraph.nodes[graph_node_id]['taxa'].add(tree_node)
            syngraph.nodes[graph_node_id]['seqs_by_taxon'][tree_node] = tree_node + "_" + LMS
    return syngraph







    

#############################################################################################

class Syngraph(nx.Graph):
    '''
    [ToDo]
    - write documentation for Syngraph about which data is stored where
    - create test-dataset 'butterfly_hard': with P. napi
    - create test-dataset 'mammals' : chromosomal BUSCOs (10: cow, possum, rat, mouse, human, chimps...)
    
    '''
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