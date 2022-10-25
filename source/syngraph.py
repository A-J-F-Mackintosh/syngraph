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


#############################################################################################
################# function for loading a graph from marker objects ##########################
#############################################################################################

def load_markerObjs(parameterObj):
    '''
    marker := orthogroup/busco
    locus := protein i.e. position of "marker" in genome of taxon

    '''
    tmp_markerObjs = []
    for infile in parameterObj.infiles:
        taxon_markerObjs = []
        taxon = infile.split("/")[-1].split(".")[0]
        df = pd.read_csv(infile, 
                sep='\t', 
                names=['name', 'seq', 'start', 'end'], 
                dtype={'name': str, 'seq': str, 'start': int, 'end': int}
                ).sort_values(['seq', 'start'], ascending=[True, True])
        for idx, (name, seq, start, end) in enumerate(df.values.tolist()):
            markerObj = MarkerObj(name=name, taxon=taxon, seq=taxon+"_"+seq, start=start, end=end)
            taxon_markerObjs.append(markerObj)
        taxon_markerObjs = sorted(taxon_markerObjs, key=attrgetter('seq', 'start'))
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


#############################################################################################
### below are functions that infer and other modules rely on ################################
#############################################################################################

# a function for getting the closest outgroup to a tree node
# outgroup must be in 'available_taxa' and cannot be a descenant of the tree_node or the tree_nodes itself
def get_closest_outgroup(tree, tree_node, available_taxa): # should this be branch or toplogical distance?
    closest_taxon_so_far = None
    closest_distance_so_far = float("inf")
    for some_tree_node in tree.search_nodes():
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

# write the genome of a taxon in terms of LMSs
def compact_synteny_1(syngraph, LMSs, taxon):
    chrom2LMS = GenomeObj()
    for graph_node_id in syngraph.nodes():
        for LMS in LMSs:
            if graph_node_id in LMSs[LMS]:
                chrom2LMS.labelled_CWAL[syngraph.nodes[graph_node_id]['seqs_by_taxon'][taxon]].add(LMS)
    chrom2LMS.CWAL = [chromosome for chromosome in chrom2LMS.labelled_CWAL.values()]
    return chrom2LMS

# returns a genome written in terms of another taxon's genome
# also include a mapping of LMS --> other taxon's chromosomes
def compact_synteny_2(taxon_A, taxon_B):
    new_taxon_B = copy.deepcopy(taxon_B)
    for B_chrom in taxon_B.labelled_CWAL:
        for LMS in taxon_B.labelled_CWAL[B_chrom]:
            for A_chrom in taxon_A.labelled_CWAL:
                if LMS in taxon_A.labelled_CWAL[A_chrom]:
                    new_taxon_B.labelled_CWAOC[B_chrom].add(A_chrom)
                    new_taxon_B.LMS_OC[LMS] = A_chrom
    return new_taxon_B

# chroms are not labelled in some GenomeObjs
# this function can give them a label (labelled_CWAL)
# or relabels them after permutation
def label_genome(tree_node, genome, label_count):
    label_state = copy.deepcopy(label_count)
    labelled_genome = copy.deepcopy(genome)
    labelled_genome.labelled_CWAL = collections.defaultdict(set)
    for chrom in labelled_genome.CWAL:
        label_state += 1
        labelled_genome.labelled_CWAL[tree_node + "_" + str(label_state)] = chrom
    return labelled_genome

def get_union_of_chrom_dicts(chrom_dict_A, chrom_dict_B):
    union_of_chrom_dicts = collections.defaultdict(set)
    for chrom_dict in chrom_dict_A, chrom_dict_B:
        for key in chrom_dict:
            union_of_chrom_dicts[key] = union_of_chrom_dicts[key].union(chrom_dict[key])
    return union_of_chrom_dicts

def get_all_LMSs_of_a_chrom(instance_of_synteny, chrom):
    LMSs = set()
    for element in instance_of_synteny[chrom]:
        LMSs.add(element)
    return LMSs

def check_for_fusions(instance_of_synteny, fusion_log, i):
    # can be sped up
    for combo in itertools.combinations(instance_of_synteny.labelled_CWAOC, 2):
        if len(instance_of_synteny.labelled_CWAOC[combo[0]].intersection(instance_of_synteny.labelled_CWAOC[combo[1]])) > 0:
            instance_of_synteny.labelled_CWAOC[combo[0]+ "_" + combo[1].split("_", 1)[1]] = \
            instance_of_synteny.labelled_CWAOC[combo[0]].union(instance_of_synteny.labelled_CWAOC[combo[1]])
            instance_of_synteny.labelled_CWAL[combo[0]+ "_" + combo[1].split("_", 1)[1]] = \
            instance_of_synteny.labelled_CWAL[combo[0]].union(instance_of_synteny.labelled_CWAL[combo[1]])
            fusion_log.append(["parent_node", i, "fusion", 1, 
                list([get_all_LMSs_of_a_chrom(instance_of_synteny.labelled_CWAL, combo[0]), 
                    get_all_LMSs_of_a_chrom(instance_of_synteny.labelled_CWAL, combo[1])])])
            del instance_of_synteny.labelled_CWAOC[combo[0]]
            del instance_of_synteny.labelled_CWAOC[combo[1]]
            del instance_of_synteny.labelled_CWAL[combo[0]]
            del instance_of_synteny.labelled_CWAL[combo[1]]
            return check_for_fusions(instance_of_synteny, fusion_log, i)
    return instance_of_synteny, fusion_log

def check_for_fissions(instance_of_synteny, fission_log, i):
    for chrom in instance_of_synteny.labelled_CWAOC:
        indices = len(instance_of_synteny.labelled_CWAOC[chrom])
        if indices > 1:
            fission_log.append(["parent_node", i, "fission", indices-1, 
                get_all_LMSs_of_a_chrom(instance_of_synteny.labelled_CWAL, chrom)])
    return fission_log

def ffsd(instance_of_synteny, rearrangement_log, i):
    instance_of_synteny, rearrangement_log = check_for_fusions(instance_of_synteny, rearrangement_log, i)
    rearrangement_log = check_for_fissions(instance_of_synteny, rearrangement_log, i)
    return rearrangement_log

def get_LMSs(syngraph, list_of_taxa, minimum):
    # could be a numpy solution...
    linked_marker_sets = collections.defaultdict(set)
    unassignable_markers = set()
    for graph_node_id in syngraph.nodes():
        target_seqs_by_taxon_set = set()
        for taxon in list_of_taxa:
            if taxon in syngraph.nodes[graph_node_id]['seqs_by_taxon']:
                target_seqs_by_taxon_set.add(syngraph.nodes[graph_node_id]['seqs_by_taxon'][taxon])
        target_seqs_by_taxon_set = frozenset(target_seqs_by_taxon_set)
        if len(target_seqs_by_taxon_set) == len(list_of_taxa):
            linked_marker_sets[target_seqs_by_taxon_set].add(graph_node_id)
        else:
            unassignable_markers.add(graph_node_id)
    filtered_linked_marker_sets = {}
    filtered_LMS_count = 0
    for LMS in linked_marker_sets:
        if len(linked_marker_sets[LMS]) >= minimum:
            filtered_linked_marker_sets["LMS_" + str(filtered_LMS_count)] = linked_marker_sets[LMS]
            filtered_LMS_count += 1
        else:
            for graph_node_id in linked_marker_sets[LMS]:
                unassignable_markers.add(graph_node_id)
    return filtered_linked_marker_sets, unassignable_markers

# bubble will combine intersecting sets recursively
# used for generating connected components
def bubble(i, some_dict):
    bubbled = False
    for j in range(len(some_dict)):
        if some_dict[j] != set() and i != j:
            if len(some_dict[i].intersection(some_dict[j])) > 0:
                some_dict[i] = some_dict[i].union(some_dict[j])
                some_dict[j] = set()
                bubbled = True
    if bubbled:
        return bubble(i, some_dict)
    else:
        return some_dict

def generate_connected_components(ios_1, ios_2, ios_3):
    connected_components = GenomeObj()
    # get unique chromosomes i.e. LMS relationships
    connected_components.CWAL = list(set([frozenset(x) for x in list(ios_1.CWAL) + list(ios_3.CWAL) + list(ios_2.CWAL)]))
    # convert back to sets
    connected_components.CWAL = [set(component) for component in connected_components.CWAL]
    # for each LMs relationship, bubble
    for i in range(len(connected_components.CWAL)):
        if len(connected_components.CWAL[i]) > 0:
            connected_components.CWAL = bubble(i, connected_components.CWAL)
    # get rid of empty sets then return
    connected_components.CWAL = [component for component in connected_components.CWAL if component != set()]
    #print(connected_components.CWAL)
    return connected_components

def generate_parsimony_genome(ios_1, ios_2, ios_3, LMSs):
    parsimony_genome = GenomeObj()
    synteny_counts = collections.defaultdict(int)
    observed_chromosomes = ios_1.CWAL + ios_3.CWAL + ios_2.CWAL
    # find LMS-LMS edges found in at least two genomes
    for combo in itertools.combinations(LMSs.keys(), 2):
        for observed_chromosome in observed_chromosomes:
            if combo[0] in observed_chromosome and combo[1] in observed_chromosome:
                synteny_counts[frozenset({combo[0], combo[1]})] += 1
    for combo_key in synteny_counts.keys():
        if synteny_counts[combo_key] == 2:
            parsimony_genome.CWAL.append(set(combo_key))
    # bubble to combine intersecting sets
    for i in range(len(parsimony_genome.CWAL)):
        if len(parsimony_genome.CWAL[i]) > 0:
            parsimony_genome.CWAL = bubble(i, parsimony_genome.CWAL)
    # add in other LMSs
    for LMS in LMSs.keys():
        already_in = False
        for chrom in parsimony_genome.CWAL:
            if LMS in chrom:
                already_in = True
        if not already_in:
            parsimony_genome.CWAL.append({LMS})
    # get rid of empty sets then return
    parsimony_genome.CWAL = [component for component in parsimony_genome.CWAL if len(component) > 0]
    return parsimony_genome

def prune_genome(ios, connected_component): # both ios and cc should be GenomeObjs
    pruned_ios = copy.deepcopy(ios)
    pruned_chroms = []
    for chrom in pruned_ios.CWAL:
        if chrom.intersection(connected_component.CWAL[0]):
            pruned_chroms.append(chrom)
    pruned_ios.CWAL = pruned_chroms
    return pruned_ios

def evaluate_genome_with_parsimony(branch_lengths, ios_1, ios_2, ios_3, possible_median, 
    best_total, best_max_branch_rate, best_genome, best_rearrangement_log, model):
    total = 0
    max_branch_rate = 0
    rearrangement_log = []
    # for each known genome in the triplet, what rearrangements does this solution require?
    for ios, i in zip([ios_1, ios_2, ios_3], range(0, 3)):
        # translocation+fission+fusion model
        if model == 3:
            rearrangement_log = ferretti(compact_synteny_2(ios, possible_median), rearrangement_log, i)
        # fission+fusion model
        elif model == 2:
            rearrangement_log  = ffsd(compact_synteny_2(ios, possible_median), rearrangement_log, i)
        # count rearrangements
        for rearrangement in rearrangement_log:
            if rearrangement[1] == i and rearrangement[2] != "translocation":
                total += rearrangement[3]
            elif rearrangement[1] == i and rearrangement[2] == "translocation":
                total += rearrangement[3] * 1.5
        # calculate branch rate, can use this to settle ties
        branch_rate = total / (branch_lengths["up"][i] + branch_lengths["down"][i])
        if branch_rate > max_branch_rate:
            max_branch_rate = branch_rate
        # if after two branches the total is already higher than the best so far then we can stop the search
        if best_total < total:
            return best_genome, best_total, best_max_branch_rate, best_rearrangement_log
    # if better than we have found so far, save the details
    if best_total > total or best_total == total and max_branch_rate < best_max_branch_rate:
        best_total = total
        best_max_branch_rate = max_branch_rate
        best_genome = possible_median
        best_rearrangement_log = rearrangement_log
    return best_genome, best_total, best_max_branch_rate, best_rearrangement_log

def parsimony_genome_solver(tree_node, parsimony_genome, ios_1, ios_2, ios_3, branch_lengths, model, label_count):
    best_rearrangement_log = []
    best_total = float("inf")
    best_max_branch_rate = float("inf")
    best_genome = []
    labelled_parsimony_genome = label_genome(tree_node, parsimony_genome, label_count)
    best_genome, best_total, best_max_branch_rate, best_rearrangement_log  = \
    evaluate_genome_with_parsimony(branch_lengths, ios_1, ios_2, ios_3, labelled_parsimony_genome, best_total, 
        best_max_branch_rate, best_genome, best_rearrangement_log, model)
    # only return the rearrangemnts on branches leading to children
    return best_genome, \
    [rearrangement for rearrangement in best_rearrangement_log if rearrangement[1] == 0 or rearrangement[1] == 1]

def exhaustive_solver(tree_node, single_component, ios_1, ios_2, ios_3, branch_lengths, model, label_count):
    best_rearrangement_log = []
    # if there are no 'tangles' then there are no events
    if len(single_component.CWAL[0]) == 1:
        labelled_single_component = label_genome(tree_node, single_component, label_count)
        return labelled_single_component, []
    # else, generate all possible candidate solutions from the single component
    # i.e. all set partitions
    else:
        possible_medians = []
        for par in more_itertools.set_partitions(single_component.CWAL[0]):
            par = [set(pa) for pa in par]
            possible_median = GenomeObj()
            possible_median.CWAL = par
            possible_medians.append(possible_median)
        # now evaluate candidate genomes
        best_total = float("inf")
        best_max_branch_rate = float("inf") # use this to solve ties
        best_genome = []
        for possible_median in possible_medians:
            possible_median = label_genome(tree_node, possible_median, label_count)
            best_genome, best_total, best_max_branch_rate, best_rearrangement_log  = \
            evaluate_genome_with_parsimony(branch_lengths, ios_1, ios_2, ios_3, possible_median, best_total, 
                best_max_branch_rate, best_genome, best_rearrangement_log, model)
        # only return the rearrangemnts on branches leading to children
        return best_genome, \
        [rearrangement for rearrangement in best_rearrangement_log if rearrangement[1] == 0 or rearrangement[1] == 1]

def permute_genome_by_fusion(genome):
    if len(genome.CWAL) == 1: # can't fuse if only 1 chromosome
        return genome
    else:
        temp_genome = copy.deepcopy(genome)
        sampled_chroms = random.sample(temp_genome.CWAL, 2)
        temp_genome.CWAL.append(sampled_chroms[0].union(sampled_chroms[1]))
        for chrom in sampled_chroms:
            temp_genome.CWAL.remove(chrom)
        return temp_genome

def permute_genome_by_fission(genome):
    if max([len(chrom) for chrom in genome.CWAL]) == 1: # can't fission single LMS chromosomes
        return genome
    else:
        temp_genome = copy.deepcopy(genome)
        samplable_chroms = [chrom for chrom in temp_genome.CWAL if len(chrom) > 1]
        sampled_chrom = random.sample(samplable_chroms, 1)[0]
        temp_genome.CWAL.remove(sampled_chrom)
        sampled_int = random.sample(range(1, len(sampled_chrom)), 1)[0]
        fission_product_1 = set()
        fission_product_2 = set()
        for index in random.sample(sampled_chrom, sampled_int):
            fission_product_1.add(index)
        for index in sampled_chrom:
            if index not in fission_product_1:
                fission_product_2.add(index)
        temp_genome.CWAL.append(fission_product_1)
        temp_genome.CWAL.append(fission_product_2)
        return temp_genome

def heuristic_solver(tree_node, ios_1, ios_2, ios_3, parsimony_genome, branch_lengths, model, label_count):
    print("[+] There are many possible genomes at this node so syngraph will undertake a heuristic search.", \
        "To avoid this, consider increasing the -m parameter.")
    print("[+] Searching ...")
    best_total = float("inf")
    best_max_branch_rate = float("inf")
    best_genome = None
    best_rearrangement_log = []
    starting_genome = copy.deepcopy(parsimony_genome)
    # peturb the parsimony genome using a decreasing random walk
    for rw_param in [[9, 2000], [7, 2000], [5, 2000], [3, 2000], [1, 2000]]:
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
                        rw_genome = permute_genome_by_fusion(rw_genome)
                    elif coin_flip == "fission":
                        rw_genome = permute_genome_by_fission(rw_genome)
            rw_genome = label_genome(tree_node, rw_genome, label_count)
            # evaluate permuted genomes
            best_genome, best_total, best_max_branch_rate, best_rearrangement_log  = \
            evaluate_genome_with_parsimony(branch_lengths, ios_1, ios_2, ios_3, rw_genome, best_total, best_max_branch_rate, 
                best_genome, best_rearrangement_log , model)
        starting_genome = best_genome
    return best_genome, \
    [rearrangement for rearrangement in best_rearrangement_log if rearrangement[1] == 0 or rearrangement[1] == 1]

def write_in_unassigned(tree_node, syngraph, taxa, LMSs, unassignable_markers):
    # if there are bugs in this function then this would be bad
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

def edit_rearrangement_log(rearrangement_log, LMSs):
    for i in range(0, len(rearrangement_log)):
        if rearrangement_log[i][2] == "fission":
            markers = set()
            for LMS in rearrangement_log[i][4]:
                for marker in LMSs[LMS]:
                    markers.add(marker)
            rearrangement_log[i][4] = markers
        else:
            markers_0 = set()
            markers_1 = set()
            for LMS in rearrangement_log[i][4][0]:
                for marker in LMSs[LMS]:
                    markers_0.add(marker)
            for LMS in rearrangement_log[i][4][1]:
                for marker in LMSs[LMS]:
                    markers_1.add(marker)
            rearrangement_log[i][4][0] = markers_0
            rearrangement_log[i][4][1] = markers_1
    return rearrangement_log

# could be good to name function variables, e.g. tree_node = None and tree_node = my_tree_node when you call

def median_genome(tree_node, working_syngraph, taxa, branch_lengths, minimum, model):
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
    # generate connected components
    connected_components = generate_connected_components(ios_1, ios_2, ios_3)
    # generate a genome where edges with frequency >= 2 are used to resolve connected components
    parsimony_genome = generate_parsimony_genome(ios_1, ios_2, ios_3, LMSs)
    print("[=] Generated {} connected components".format(len(connected_components.CWAL)))
    solved_connected_components = GenomeObj()
    triplet_rearrangement_log = []
    label_count = 0
    # for fission+fusion we can just use the parsimony genome
    if model == 2:
        solved_connected_component, new_rearrangements = parsimony_genome_solver(tree_node, 
            parsimony_genome, ios_1, ios_2, ios_3, branch_lengths, model, label_count)
        # update solution, log, and labels
        solved_connected_components.labelled_CWAL.update(solved_connected_component.labelled_CWAL)
        triplet_rearrangement_log += new_rearrangements
        label_count += len(solved_connected_component.labelled_CWAL.keys())
    # for more complex models
    # check whether, for each connected component, the best arrangement can be found exhaustively
    # and then call the appropriate solving function
    else:
        for connected_component in connected_components.CWAL:
            single_component = GenomeObj()
            single_component.CWAL.append(connected_component)
            # pruning gives the connected component specific data
            if len(connected_component) > 10:
                # solve heuristicly
                solved_connected_component, new_rearrangements = heuristic_solver(tree_node, 
                    prune_genome(ios_1, single_component), prune_genome(ios_2, single_component), 
                    prune_genome(ios_3, single_component), prune_genome(parsimony_genome, single_component), branch_lengths, 
                    model, label_count) 
            else:
                # solve exhaustively
                solved_connected_component, new_rearrangements = exhaustive_solver(tree_node, single_component, 
                    prune_genome(ios_1, single_component), prune_genome(ios_2, single_component), 
                    prune_genome(ios_3, single_component), branch_lengths, model, label_count)
            # update solution, log, and labels
            solved_connected_components.labelled_CWAL.update(solved_connected_component.labelled_CWAL)
            triplet_rearrangement_log += new_rearrangements
            label_count += len(solved_connected_component.labelled_CWAL.keys())
    print("[=] Found a median genome with {} chromosomes".format(len(solved_connected_components.labelled_CWAL.keys())))
    # from solved_connected_components, add this new ancestral genome to a syngraph
    working_syngraph.graph['taxa'].add(tree_node)
    for chrom in solved_connected_components.labelled_CWAL.keys():
        for LMS in solved_connected_components.labelled_CWAL[chrom]:
            for graph_node_id in LMSs[LMS]:
                working_syngraph.nodes[graph_node_id]['taxa'].add(tree_node)
                working_syngraph.nodes[graph_node_id]['seqs_by_taxon'][tree_node] = chrom
    # replace LMSs in the rearrangment_log with the actual markers
    triplet_rearrangement_log = edit_rearrangement_log(triplet_rearrangement_log, LMSs)
    # finally, write in LMSs/markers that were too small or missing from a taxon, but can be assigned by parismony
    working_syngraph = write_in_unassigned(tree_node, working_syngraph, taxa, LMSs, unassignable_markers)
    return working_syngraph, triplet_rearrangement_log

def tree_traversal(syngraph, params):
    # write a log header
    log = [["#parent", "child", "event", "multiplicity", "ref_seqs"]]
    # copy syngraph
    traversal_0_syngraph = copy.deepcopy(syngraph)
    # define which taxa are extant and so can be sampled from the start
    available_taxa = set([leaf.name for leaf in params.tree.get_leaves()])
    print("[+] Starting traversal ...")
    print("[+] ========================================================================")
    # are there still unsolved nodes? If yes, then keeping solving
    nodes_left = True
    while nodes_left:
        # store info about easiest triplet to solve in a list
        best_triplet = [None, None, float("inf")]
        for tree_node in params.tree.traverse(strategy='postorder'):
            # is there info to solve this node?
            if not tree_node.is_leaf() and not tree_node.is_root() and not tree_node.name in available_taxa:
                child_1 = tree_node.get_children()[0].name
                child_2 = tree_node.get_children()[1].name
                if child_1 in available_taxa and child_2 in available_taxa:
                    outgroup = get_closest_outgroup(params.tree, tree_node, available_taxa)
                    # how easy is it to solve?
                    evaluation = evaluate_triplet([child_1, child_2, outgroup], traversal_0_syngraph, params.minimum)
                    if evaluation < best_triplet[2]:
                        best_triplet = [[child_1, child_2, outgroup], tree_node.name, evaluation]
                        if evaluation == 0:
                            break
        # if no more triplets to solve then we are done
        if best_triplet == [None, None, float("inf")]:
            nodes_left = False
        else:
            # solve the easiest triplet found
            branch_lengths = collections.defaultdict(list)
            for taxon in best_triplet[0]:
                up_distance, down_distance = get_branch_distances(params.tree, best_triplet[1], taxon)
                branch_lengths["up"].append(up_distance)
                branch_lengths["down"].append(down_distance)
            print("[+] Inferring median genome for {} using data from {}, {}, and {} ...". format(best_triplet[1], best_triplet[0][0], 
                best_triplet[0][1], best_triplet[0][2]))
            # call the solver
            traversal_0_syngraph, rearrangement_log = \
            median_genome(best_triplet[1], traversal_0_syngraph, best_triplet[0], branch_lengths, params.minimum, params.model)
            available_taxa.add(best_triplet[1])
            # add rearrangements to log
            for rearrangement in rearrangement_log:
                log.append([best_triplet[1], [best_triplet[0][0], best_triplet[0][1]][rearrangement[1]], rearrangement[2], 
                    rearrangement[3], rearrangement[4]])
            print("[=] ========================================================================")
    return traversal_0_syngraph, log

# returns the minimum fission+fusion distance between a triplet of taxa
# the smaller this value, the easier it is to find a median genome
def evaluate_triplet(taxa, syngraph, minimum):
    smallest = float("inf")
    LMSs, unassignables = get_LMSs(syngraph, taxa, minimum)
    ios_1 = compact_synteny_1(syngraph, LMSs, taxa[0])
    ios_2 = compact_synteny_1(syngraph, LMSs, taxa[1])
    ios_3 = compact_synteny_1(syngraph, LMSs, taxa[2])
    for i, combo in enumerate(itertools.combinations([ios_1, ios_2, ios_3], 2)):
        rearrangement_log = []
        rearrangement_log = ffsd(compact_synteny_2(combo[0], combo[1]), rearrangement_log, i)
        total = 0
        for rearrangement in rearrangement_log:
            total += rearrangement[3]
        if total < smallest:
            smallest = total
    return smallest

#############################################################################################
###### below are functions for a fission + fusion + translocation model #####################
#############################################################################################

def ferretti(instance_of_synteny, rearrangement_log, i):
    # ios -> chrom -> element(s) (chrom(s) mapped to in target genome) --> LMS(s)
    # remove chroms that are already solved
    instance_of_synteny = strip_chromes(instance_of_synteny)
    if len(instance_of_synteny.labelled_CWAOC.keys()) == 0:
        return rearrangement_log
    else:
        # get element multiplicy (l) and minimum syntenic element multiplicity (r_min)
        # actually not sure if we need r_min...
        l_dict = get_l(instance_of_synteny)
        # if there are elements with l=1 # should we do this first?
        if 1 in l_dict.values():
            # pick a chrom with l=1
            small_l_element = sample_small_l(l_dict, 1)
            # fission
            instance_of_synteny, rearrangement_log = implement_fission_1(instance_of_synteny, small_l_element, 
                rearrangement_log, i)
            return ferretti(instance_of_synteny, rearrangement_log, i)
        # if elements with l=2
        elif 2 in l_dict.values():
            # pick a chrom with l=2
            small_l_element = sample_small_l(l_dict, 2)
            # check if elements are in chroms with lens >=2
            # if yes, then translocate, else fuse
            l_2_status, chrom_content = check_l_2_status(instance_of_synteny, small_l_element)
            if l_2_status == "translocation":
                instance_of_synteny, rearrangement_log = implement_translocation_2(instance_of_synteny, chrom_content, 
                    small_l_element, rearrangement_log, i)
                return ferretti(instance_of_synteny, rearrangement_log, i)
            elif l_2_status == "fusion":
                instance_of_synteny, rearrangement_log = implement_fusion_2(instance_of_synteny, chrom_content, 
                    rearrangement_log, i)
                return ferretti(instance_of_synteny, rearrangement_log, i)
        else:
            # if only elements with l>2, fuse
            big_l_element = sample_big_l(l_dict)
            instance_of_synteny, rearrangement_log = implement_fusion_3(instance_of_synteny, big_l_element, rearrangement_log, i)
            return ferretti(instance_of_synteny, rearrangement_log, i)


def strip_chromes(ios):
    chromes_to_be_stripped = []
    for chrom in ios.labelled_CWAOC:
        if len(ios.labelled_CWAOC[chrom]) == 1:
            multiplicity = 0
            for chrom_2 in ios.labelled_CWAOC:
                if list(ios.labelled_CWAOC[chrom])[0] in ios.labelled_CWAOC[chrom_2]:
                    multiplicity += 1
            if multiplicity == 1: # i.e. it is a chrom that does not need to be rearranged
                chromes_to_be_stripped.append(chrom)
    for chrom in chromes_to_be_stripped:
        ios.labelled_CWAL.pop(chrom, None)
        ios.labelled_CWAOC.pop(chrom, None)
    return ios

def get_l(ios):
    l_dict = collections.defaultdict(int)
    for chrom in ios.labelled_CWAOC:
        for element in ios.labelled_CWAOC[chrom]:
            l_dict[element] += 1
    return l_dict

# def get_r_min(ios, l_dict):
#     r_min_dict = collections.defaultdict(lambda:float("inf"))
#     for chrom in ios:
#         for element in ios[chrom]:
#             for other_element in ios[chrom]:
#                 if other_element != element:
#                     if l_dict[other_element] < r_min_dict[element]:
#                         r_min_dict[element] = l_dict[other_element]
#     return r_min_dict

def sample_small_l(l_dict, l_min):
    potential_small_l = set()
    for element in l_dict:
        if l_dict[element] == l_min:
            potential_small_l.add(element)
    return random.sample(potential_small_l, 1)[0]

def implement_fission_1(ios, element, fission_log, i):
    for chrom in ios.labelled_CWAOC:
        if element in ios.labelled_CWAOC[chrom]:
            ios.labelled_CWAOC[chrom + "_a"].add(element)
            for LMS in ios.labelled_CWAL[chrom]:
                if ios.LMS_OC[LMS] == element:
                    ios.labelled_CWAL[chrom + "_a"].add(LMS)
            for other_elements in ios.labelled_CWAOC[chrom]:
                if other_elements != element:
                    ios.labelled_CWAOC[chrom + "_b"].add(other_elements)
                    for LMS in ios.labelled_CWAL[chrom]:
                        if ios.LMS_OC[LMS] == other_elements:
                            ios.labelled_CWAL[chrom + "_b"].add(LMS)
            fission_log.append(["parent_node", i, "fission", 1, get_all_LMSs_of_a_chrom(ios.labelled_CWAL, chrom)])
            ios.labelled_CWAOC.pop(chrom, None)
            ios.labelled_CWAL.pop(chrom, None)
            return ios, fission_log

def check_l_2_status(ios, element):
    chrom_content = []
    for chrom in ios.labelled_CWAOC:
        if element in ios.labelled_CWAOC[chrom]:
            chrom_content.append(chrom)
    if len(ios.labelled_CWAOC[chrom_content[0]]) >= 2 and len(ios.labelled_CWAOC[chrom_content[1]]) >= 2:
        return "translocation", chrom_content
    else:
        return "fusion", chrom_content

def implement_fusion_2(ios, chrom_content, fusion_log, i):
    ios.labelled_CWAOC[chrom_content[0] + "_" + chrom_content[1].split("_", 1)[1]] = \
    ios.labelled_CWAOC[chrom_content[0]].union(ios.labelled_CWAOC[chrom_content[1]])
    ios.labelled_CWAL[chrom_content[0] + "_" + chrom_content[1].split("_", 1)[1]] = \
    ios.labelled_CWAL[chrom_content[0]].union(ios.labelled_CWAL[chrom_content[1]])
    fusion_log.append(["parent_node", i, "fusion", 1, 
        list([get_all_LMSs_of_a_chrom(ios.labelled_CWAL, chrom_content[0]), 
            get_all_LMSs_of_a_chrom(ios.labelled_CWAL, chrom_content[1])])])
    ios.labelled_CWAOC.pop(chrom_content[0], None)
    ios.labelled_CWAOC.pop(chrom_content[1], None)
    ios.labelled_CWAL.pop(chrom_content[0], None)
    ios.labelled_CWAL.pop(chrom_content[1], None)
    return ios, fusion_log

def implement_translocation_2(ios, chrom_content, element, translocation_log, i):
    ios.labelled_CWAOC[chrom_content[0] + "_a_" + chrom_content[1].split("_", 1)[1] + "_a"].add(element)
    for index in 0, 1:
        for LMS in ios.labelled_CWAL[chrom_content[index]]:
            if ios.LMS_OC[LMS] == element:
                ios.labelled_CWAL[chrom_content[0] + "_a_" + chrom_content[1].split("_", 1)[1] + "_a"].add(LMS)
    for index in 0, 1:
        for other_elements in ios.labelled_CWAOC[chrom_content[index]]:
            if other_elements != element:
                ios.labelled_CWAOC[chrom_content[0] + "_b_" + chrom_content[1].split("_", 1)[1] + "_b"].add(other_elements)
                for LMS in ios.labelled_CWAL[chrom_content[index]]:
                    if ios.LMS_OC[LMS] == other_elements:
                        ios.labelled_CWAL[chrom_content[0] + "_b_" + chrom_content[1].split("_", 1)[1] + "_b"].add(LMS)
    translocation_log.append(["parent_node", i, "translocation", 1, 
        list([get_all_LMSs_of_a_chrom(ios.labelled_CWAL, chrom_content[0]), 
            get_all_LMSs_of_a_chrom(ios.labelled_CWAL, chrom_content[1])])])
    ios.labelled_CWAOC.pop(chrom_content[0], None)
    ios.labelled_CWAOC.pop(chrom_content[1], None)
    ios.labelled_CWAL.pop(chrom_content[0], None)
    ios.labelled_CWAL.pop(chrom_content[1], None)
    return ios, translocation_log

def sample_big_l(l_dict):
    potential_big_l = set()
    l_min = min(l_dict.values())
    for element in l_dict:
        if l_dict[element] == l_min:
            potential_big_l.add(element)
    return random.sample(potential_big_l, 1)[0]

def implement_fusion_3(ios, element, fusion_log, i):
    # currently choosing fusions with maximal intersections
    fusion_candidates = []
    for chrom in ios.labelled_CWAOC:
        if element in ios.labelled_CWAOC[chrom]:
            fusion_candidates.append(chrom)
    best_intersection = 0
    best_pair = []
    for combo in itertools.combinations(fusion_candidates, 2):
        if len(ios.labelled_CWAOC[combo[0]].intersection(ios.labelled_CWAOC[combo[1]])) > best_intersection:
                best_intersection = len(ios.labelled_CWAOC[combo[0]].intersection(ios.labelled_CWAOC[combo[1]]))
                best_pair = list(combo)
    ios.labelled_CWAOC[best_pair[0] + "_" + best_pair[1].split("_", 1)[1]] = \
    ios.labelled_CWAOC[best_pair[0]].union(ios.labelled_CWAOC[best_pair[1]])
    ios.labelled_CWAL[best_pair[0] + "_" + best_pair[1].split("_", 1)[1]] = \
    ios.labelled_CWAL[best_pair[0]].union(ios.labelled_CWAL[best_pair[1]])
    fusion_log.append(["parent_node", i, "fusion", 1, list([get_all_LMSs_of_a_chrom(ios.labelled_CWAL, best_pair[0]), 
                    get_all_LMSs_of_a_chrom(ios.labelled_CWAL, best_pair[1])])])
    ios.labelled_CWAOC.pop(best_pair[0], None)
    ios.labelled_CWAOC.pop(best_pair[1], None)
    ios.labelled_CWAL.pop(best_pair[0], None)
    ios.labelled_CWAL.pop(best_pair[1], None)
    return ios, fusion_log


# need to double check the Ferretti et al. paper to see how my algorithm differs

#############################################################################################
######## below are functions for formatting the output of infer #############################
#############################################################################################

# map chromosomes involved in rearrangements to chromosomes belonging to some taxon (tip or internal node)
def map_log(log, reference_taxon, reference_dict, syngraph, minimum):
    for i in range(1, len(log)): # start at 1 to avoid parsing the header line
        if log[i][2] == "fission":
            chroms = collections.defaultdict(int)
            for marker in log[i][4]:
                if reference_taxon:
                    if reference_taxon in syngraph.nodes[marker]['seqs_by_taxon'].keys():
                        chroms[syngraph.nodes[marker]['seqs_by_taxon'][reference_taxon]] += 1
                elif reference_dict:
                    if marker in reference_dict.keys():
                        chroms[reference_dict[marker]] += 1
            chroms = pop_bad_mappings(chroms, minimum)
            log[i][4] = list(chroms.keys())  
        else:
            for anc_chrom in 0, 1:
                chroms = collections.defaultdict(int)
                for marker in log[i][4][anc_chrom]:
                    if reference_taxon:
                        if reference_taxon in syngraph.nodes[marker]['seqs_by_taxon'].keys():
                            chroms[syngraph.nodes[marker]['seqs_by_taxon'][reference_taxon]] += 1
                    elif reference_dict:
                        if marker in reference_dict.keys():
                            chroms[reference_dict[marker]] += 1
                chroms = pop_bad_mappings(chroms, minimum)
                log[i][4][anc_chrom] = list(chroms.keys())
    return log

# get rid of mappings below minimum (error threshold)
def pop_bad_mappings(chroms, minimum):
    bad_mappings = []
    for mapping in chroms:
        if chroms[mapping] < minimum:
            bad_mappings.append(mapping)
    for bad_mapping in bad_mappings:
        chroms.pop(bad_mapping, None)
    return chroms

def clusters_by_descent(log, tree, syngraph):
    clusters = []
    for node in tree.traverse(strategy="preorder"):
        if not node.is_leaf() and not node.is_root():
            in_cluster = False
            for cluster in clusters:
                if node.name in cluster:
                    in_cluster = True
            if not in_cluster:
                clusters.append([node.name])
            for child in node.get_children():
                rearrangement = False
                for entry in log:
                    if entry[0] == node.name and entry[1] == child.name:
                        rearrangement = True
                if rearrangement:
                    clusters.append([child.name])
                else:
                    for cluster in clusters:
                        if node.name in cluster:
                            cluster.append(child.name)
    cluster_ID = 0
    formatted_clusters = []
    for cluster in clusters:
        cluster_ID += 1
        for taxon in cluster:
            formatted_clusters.append(["cluster_" + str(cluster_ID), taxon, get_karyotype(taxon, syngraph)])
    return formatted_clusters

def get_karyotype(taxon, syngraph):
    seqs = set()
    for graph_node_id in syngraph.nodes():
        if taxon in syngraph.nodes[graph_node_id]['seqs_by_taxon']:
            seqs.add(syngraph.nodes[graph_node_id]['seqs_by_taxon'][taxon])
    return len(seqs)


#############################################################################################
############################### defining classes ############################################
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

    def save(self, parameterObj, check_consistency, with_ancestors):
        if with_ancestors:
            graph_file = '%s.with_ancestors.pickle' % parameterObj.outprefix
        else:
            graph_file = '%s.pickle' % parameterObj.outprefix
        nx.write_gpickle(self, graph_file)
        if check_consistency:
            syngraph_loaded = Syngraph()
            syngraph_loaded.from_file(graph_file)
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
                self.add_node(markerObj.name, taxa=set(), terminal=set(), seqs_by_taxon={}, starts_by_taxon={}, 
                    ends_by_taxon={})             
            # add taxon to node (useful if '--missing')
            self.nodes[markerObj.name]['taxa'].add(markerObj.taxon)
            self.nodes[markerObj.name]['seqs_by_taxon'][markerObj.taxon] = markerObj.seq
            self.nodes[markerObj.name]['starts_by_taxon'][markerObj.taxon] = markerObj.start
            self.nodes[markerObj.name]['ends_by_taxon'][markerObj.taxon] = markerObj.end
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
        self.graph['marker_ids_by_seq_id_by_taxon'] = marker_ids_by_seq_id_by_taxon # do we need/want this?

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
        node_non_lonely_count = 0
        node_complete_count = 0
        for graph_node_id in self.nodes:
            if len(self.nodes[graph_node_id]['taxa']) > 1:
                node_non_lonely_count += 1
                if len(self.nodes[graph_node_id]['taxa']) == len(self.graph['taxa']):
                    node_complete_count += 1
        edge_total_count = nx.number_of_edges(self)
        connected_component_count = nx.number_connected_components(self)
        print("[=] ====================================")
        print("[=] Taxa = %s" % taxon_count)
        print("[=] Nodes (Markers) = %s" % node_total_count)
        print("[=] Nodes (Markers) shared by > 1 taxon = %s" % node_non_lonely_count)
        print("[=] Nodes (Markers) shared by all taxa = %s" % node_complete_count)
        print("[=] Distinct Edges (Adjacencies) = %s" % edge_total_count)
        print("[=] Subgraphs (connected components) = %s" % connected_component_count)
        print("[=] ====================================")


class MarkerObj():
    def __init__(self, name=None, desc=None, status=None, taxon=None, seq=None, start=None, end=None):
        self.name = name
        self.desc = desc if desc is not None else name
        self.status = status
        self.taxon = taxon
        self.seq = seq
        self.start = start
        self.end = end

    def __repr__(self):
        return "MarkerObj(name=%r, desc=%r, status=%r, taxon=%r, seq=%r, start=%d, end=%d)" % (
            self.name, self.desc, self.status, self.taxon, self.seq, self.start, self.end) 

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
            return int(max((self.start - other.end), (other.start - self.end)))
        return float("nan")

class GenomeObj():
    # CWAL = chromosomes written as LMSs
    # CWAOC = chromosomes written as other chromosomes
    # either a list of sets or, if labelled, a dict of sets
    # LMS_OC = mapping of LMS to other chromosomes
    # This could all be done with one nested dict, but I tried that and it was ugly
    def __init__(self, labelled_CWAL=None, CWAL=None, labelled_CWAOC=None, CWAOC=None, LMS_OC=None, taxon=None):
        if labelled_CWAL is None:
            labelled_CWAL = collections.defaultdict(set)
        if CWAL is None:
            CWAL = []
        if labelled_CWAOC is None:
            labelled_CWAOC = collections.defaultdict(set)
        if CWAOC is None:
            CWAOC = []
        if LMS_OC is None:
            LMS_OC = {}
        self.labelled_CWAL = labelled_CWAL
        self.CWAL = CWAL
        self.labelled_CWAOC = labelled_CWAOC
        self.CWAOC = CWAOC
        self.LMS_OC = LMS_OC
        self.taxon = taxon