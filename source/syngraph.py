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
### below are functions for implementing fusion + fission models ############################
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
    chrom2LMS = collections.defaultdict(set)
    for graph_node_id in syngraph.nodes():
        for LMS in LMSs:
            if graph_node_id in LMSs[LMS]:
                chrom2LMS[syngraph.nodes[graph_node_id]['seqs_by_taxon'][taxon]].add(LMS)
    return chrom2LMS

def compact_synteny_2(taxon_A, taxon_B):
    chrom2chrom2LMS = collections.defaultdict(lambda : collections.defaultdict(set))
    for B_chrom in taxon_B:
        for LMS in taxon_B[B_chrom]:
            for A_chrom in taxon_A:
                if LMS in taxon_A[A_chrom]:
                    chrom2chrom2LMS[B_chrom][A_chrom].add(LMS)
    return chrom2chrom2LMS

def label_genome(tree_node, genome, label_count):
    # chroms have no names when the genome is represented as a list of sets so am giving them a number
    label_state = copy.deepcopy(label_count)
    labelled_genome = {}
    for chrom in genome:
        label_state += 1
        labelled_genome[tree_node + "_" + str(label_state)] = chrom
    return labelled_genome

def check_for_fusions(instance_of_synteny, fusion_log, i):
    # can be sped up
    for combo in itertools.combinations(instance_of_synteny.keys(), 2):
        if set(instance_of_synteny[combo[0]].keys()).intersection(set(instance_of_synteny[combo[1]].keys())):
            instance_of_synteny[combo[0]+ "_" + combo[1].split("_", 1)[1]] = \
            get_union_of_chrom_dicts(instance_of_synteny[combo[0]], instance_of_synteny[combo[1]])
            fusion_log.append(["parent_node", i, "fusion", (combo[0], combo[1]), 
                list([get_all_LMSs_of_a_chrom(instance_of_synteny, combo[0]), 
                    get_all_LMSs_of_a_chrom(instance_of_synteny, combo[1])])])
            del instance_of_synteny[combo[0]]
            del instance_of_synteny[combo[1]]
            return check_for_fusions(instance_of_synteny, fusion_log, i)
    return instance_of_synteny, fusion_log

def get_union_of_chrom_dicts(chrom_dict_A, chrom_dict_B):
    union_of_chrom_dicts = collections.defaultdict(set)
    for chrom_dict in chrom_dict_A, chrom_dict_B:
        for key in chrom_dict:
            union_of_chrom_dicts[key] = union_of_chrom_dicts[key].union(chrom_dict[key])
    return union_of_chrom_dicts

def get_all_LMSs_of_a_chrom(instance_of_synteny, chrom):
    LMSs = set()
    for chromie in instance_of_synteny[chrom]:
        LMSs = LMSs.union(instance_of_synteny[chrom][chromie])
    return LMSs

def check_for_fissions(instance_of_synteny, fission_log, i):
    for chrom in instance_of_synteny:
        indices = len(instance_of_synteny[chrom].keys())
        if indices > 1:
            fission_log.append(["parent_node", i, "fission", list([chrom, indices]), get_all_LMSs_of_a_chrom(instance_of_synteny, chrom)])
    return fission_log

def ffsd(instance_of_synteny, rearrangement_log, i):
    #print("# instance_of_synteny", instance_of_synteny, "\n")
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

def prune_ios(ioses, connected_component):
    pruned_ioses = []
    for ios in ioses:
        pruned_ios = {}
        for chrom in ios:
            if ios[chrom].intersection(connected_component):
                pruned_ios[chrom] = ios[chrom]
        pruned_ioses.append(pruned_ios)
    return pruned_ioses[0], pruned_ioses[1], pruned_ioses[2]

def prune_pg(parsimony_genome, connected_component):
    pruned_pg = []
    for chrom in parsimony_genome:
        if chrom.intersection(connected_component):
            pruned_pg.append(chrom)
    return pruned_pg

def evaluate_genome_with_parsimony(branch_lengths, ios_1, ios_2, ios_3, possible_median, 
    best_total, best_max_branch_rate, best_genome, best_rearrangement_log, model):
    total = 0
    max_branch_rate = 0
    rearrangement_log = []
    for ios, i in zip([ios_1, ios_2, ios_3], range(0, 3)):
        if model == 3:
            rearrangement_log = ferretti(compact_synteny_2(ios, possible_median), rearrangement_log, i)
        elif model == 2:
            rearrangement_log  = ffsd(compact_synteny_2(ios, possible_median), rearrangement_log, i)
        for rearrangement in rearrangement_log:
            if rearrangement[1] == i:
                if rearrangement[2] == "fission":
                    total += rearrangement[3][1] - 1
                else:
                    total += 1
        branch_rate = total / (branch_lengths["up"][i] + branch_lengths["down"][i])
        if branch_rate > max_branch_rate:
            max_branch_rate = branch_rate
        if best_total < total:
            return best_genome, best_total, best_max_branch_rate, best_rearrangement_log
    if best_total > total or best_total == total and max_branch_rate < best_max_branch_rate:
        best_total = total
        best_max_branch_rate = max_branch_rate
        best_genome = possible_median
        best_rearrangement_log = rearrangement_log
    return best_genome, best_total, best_max_branch_rate, best_rearrangement_log

def exhaustive_solver(tree_node, connected_component, ios_1, ios_2, ios_3, branch_lengths, model, label_count):
    best_rearrangement_log = []
    # if there are no 'tangles' then there are no events
    if len(connected_component) == 1:
        connected_component = [connected_component]
        labelled_cc = label_genome(tree_node, connected_component, label_count)
        return labelled_cc, []
    # else, generate all possible candidate genomes from the connected components
    else:
        possible_medians = []
        for par in more_itertools.set_partitions(connected_component):
            par = [set(pa) for pa in par]
            possible_medians.append(par)
        # now evaluate candidate genomes
        best_total = float("inf")
        best_max_branch_rate = float("inf") # use this to solve ties
        best_genome = []
        for possible_median in possible_medians:
            possible_median = label_genome(tree_node, possible_median, label_count)
            best_genome, best_total, best_max_branch_rate, best_rearrangement_log  = \
            evaluate_genome_with_parsimony(branch_lengths, ios_1, ios_2, ios_3, possible_median, best_total, 
                best_max_branch_rate, best_genome, best_rearrangement_log, model)
        return best_genome, \
        [rearrangement for rearrangement in best_rearrangement_log if rearrangement[1] == 0 or rearrangement[1] == 1]

# permutations by fusion should only fuse LMSs found in the same connected component
def permute_genome_by_fusion(genome):
    if len(genome) == 1:
        return genome
    else:
        temp_genome = copy.deepcopy(genome)
        sampled_chroms = random.sample(temp_genome, 2)
        for chrom in sampled_chroms:
            temp_genome.remove(chrom)
            temp_genome.append(sampled_chroms[0].union(sampled_chroms[1]))
        return temp_genome

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

def heuristic_solver(tree_node, connected_component, ios_1, ios_2, ios_3, parsimony_genome, branch_lengths, model, label_count):
    print("[+] There are many possible genomes at this node so syngraph will undertake a heuristic search.", \
        "To avoid this, consider increasing the -m parameter.")
    print("[+] Searching ...")
    best_total = float("inf")
    best_max_branch_rate = float("inf")
    best_genome = []
    best_rearrangement_log = []
    starting_genome = prune_pg(parsimony_genome, connected_component)
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
                        rw_genome = permute_genome_by_fusion(rw_genome)
                    elif coin_flip == "fission":
                        rw_genome = permute_genome_by_fission(rw_genome)
            rw_genome = label_genome(tree_node, rw_genome, label_count)
            # evaluate permuted genomes
            best_genome, best_total, best_max_branch_rate, best_rearrangement_log  = \
            evaluate_genome_with_parsimony(branch_lengths, ios_1, ios_2, ios_3, rw_genome, best_total, best_max_branch_rate, 
                best_genome, best_rearrangement_log , model)
        starting_genome = [value for value in best_genome.values()]
    return best_genome, \
    [rearrangement for rearrangement in best_rearrangement_log if rearrangement[1] == 0 or rearrangement[1] == 1]

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

# name function variables, e.g. tree_node = None and tree_node = my_tree_node when you call

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
    connected_components = generate_connected_components(ios_1, ios_2, ios_3)
    parsimony_genome = generate_parsimony_genome(ios_1, ios_2, ios_3, LMSs)
    print("[=] Generated {} connected components".format(len(connected_components)))
    # now check whether, for each connected component, the best arrangement can be found exhaustively
    # and then call the appropriate solving function
    solved_connected_components = {}
    triplet_rearrangement_log = []
    label_count = 0
    for connected_component in connected_components:
        pruned_ios_1, pruned_ios_2, pruned_ios_3 = prune_ios([ios_1, ios_2, ios_3], connected_component)
        if len(connected_component) > 10:
            # how to use parsimony genome to solve tricky CCs?
            solved_connected_component, new_rearrangements = heuristic_solver(tree_node, connected_component, pruned_ios_1, 
                pruned_ios_2, pruned_ios_3, parsimony_genome, branch_lengths, model, label_count)
        else:
            solved_connected_component, new_rearrangements = exhaustive_solver(tree_node, connected_component, pruned_ios_1, 
                pruned_ios_2, pruned_ios_3, branch_lengths, model, label_count)
        solved_connected_components.update(solved_connected_component)
        triplet_rearrangement_log += new_rearrangements
        label_count += len(solved_connected_component)
    print("[=] Found a median genome with {} chromosomes".format(len(solved_connected_components)))
    print("solved_ccs:", solved_connected_components)
    print("rearrangement_log:", triplet_rearrangement_log)
    # from solved_connected_components, add this new ancestral genome to a syngraph
    working_syngraph.graph['taxa'].add(tree_node)
    for chrom in solved_connected_components:
        for LMS in solved_connected_components[chrom]:
            for graph_node_id in LMSs[LMS]:
                working_syngraph.nodes[graph_node_id]['taxa'].add(tree_node)
                working_syngraph.nodes[graph_node_id]['seqs_by_taxon'][tree_node] = chrom
    # replace LMSs in the rearrangment_log with the markers
    triplet_rearrangement_log = edit_rearrangement_log(triplet_rearrangement_log, LMSs) 
    # finally, write in LMSs/markers that were too small or missing from a taxon, but can be assigned by parismony
    working_syngraph = write_in_unassigned(tree_node, working_syngraph, taxa, LMSs, unassignable_markers)
    return working_syngraph, triplet_rearrangement_log

def tree_traversal(syngraph, params):
    # write a log
    log = [["parent", "child", "event", "info", "extra"]]
    # copy syngraph
    traversal_0_syngraph = copy.deepcopy(syngraph)
    # define which taxa are extant and so can be sampled from the start
    available_taxa = set()
    for leaf in params.tree.get_leaves():
        available_taxa.add(leaf.name)
    print("[+] Starting first traversal ...")
    print("[+] ========================================================================")
    nodes_left = True
    while nodes_left == True:
        best_triplet = [None, None, float("inf")]
        for tree_node in params.tree.traverse(strategy='postorder'):
            if not tree_node.is_leaf() and not tree_node.is_root() and not tree_node.name in available_taxa:
                child_1 = tree_node.get_children()[0].name
                child_2 = tree_node.get_children()[1].name
                if child_1 in available_taxa and child_2 in available_taxa:
                    outgroup = get_closest_outgroup(params.tree, tree_node, available_taxa)
                    evaluation = evaluate_triplet([child_1, child_2, outgroup], traversal_0_syngraph, params.minimum)
                    if evaluation < best_triplet[2]:
                        best_triplet = [[child_1, child_2, outgroup], tree_node.name, evaluation]
        if best_triplet == [None, None, float("inf")]:
            nodes_left = False
        else:
            branch_lengths = collections.defaultdict(list)
            for taxon in best_triplet[0]:
                up_distance, down_distance = get_branch_distances(params.tree, best_triplet[1], taxon)
                branch_lengths["up"].append(up_distance)
                branch_lengths["down"].append(down_distance)
            print("[+] Inferring median genome for {} using data from {}, {}, and {} ...". format(best_triplet[1], best_triplet[0][0], 
                best_triplet[0][1], best_triplet[0][2]))
            traversal_0_syngraph, rearrangement_log = \
            median_genome(best_triplet[1], traversal_0_syngraph, best_triplet[0], branch_lengths, params.minimum, params.model)
            available_taxa.add(best_triplet[1])
            for rearrangement in rearrangement_log:
                log.append([best_triplet[1], [best_triplet[0][0], best_triplet[0][1]][rearrangement[1]], rearrangement[2], 
                    rearrangement[3], rearrangement[4]])
            print("[=] ========================================================================")
    return traversal_0_syngraph, log


# should tree node be saved as name or node??^


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
            if rearrangement[2] == "fission":
                    total += rearrangement[3][1] - 1
            else:
                total += 1
        if total < smallest:
            smallest = total
    return smallest





#############################################################################################
###### Below are functions for a fission + fusion + translocation model #####################
#############################################################################################

def ferretti(instance_of_synteny, rearrangement_log, i):
    #print("# instance_of_synteny", instance_of_synteny, "\n")
    instance_of_synteny = strip_chromes(instance_of_synteny)
    #print("# stripped instance_of_synteny", instance_of_synteny, "\n")
    if len(instance_of_synteny) == 0:
        return rearrangement_log
    else:
        # get element multiplicy (l) and minimum syntenic element multiplicity (r_min)
        # actually not sure if we need r_min
        l_dict = get_l(instance_of_synteny)
        #print("# ldict", l_dict, "\n")
        #r_min_dict = get_r_min(instance_of_synteny, l_dict)
        # if there are elements with l=1
        if 1 in l_dict.values():
            small_l_element = sample_small_l(l_dict, 1)
            instance_of_synteny, rearrangement_log = implement_fission_1(instance_of_synteny, small_l_element, 
                rearrangement_log, i)
            return ferretti(instance_of_synteny, rearrangement_log, i)
        elif 2 in l_dict.values():
            small_l_element = sample_small_l(l_dict, 2)
            # check if both elements are in chroms of >=len(2)
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
    for chrom in ios:
        if len(ios[chrom].keys()) == 1:
            multiplicity = 0
            for chrom_2 in ios:
                if list(ios[chrom].keys())[0] in ios[chrom_2].keys():
                    multiplicity += 1
            if multiplicity == 1:
                chromes_to_be_stripped.append(chrom)
    for chrom in chromes_to_be_stripped:
        ios.pop(chrom, None)
    return ios

def get_l(ios):
    l_dict = collections.defaultdict(int)
    for chrom in ios:
        for element in ios[chrom]:
            l_dict[element] += 1
    return l_dict

def get_r_min(ios, l_dict):
    r_min_dict = collections.defaultdict(lambda:float("inf"))
    for chrom in ios:
        for element in ios[chrom]:
            for other_element in ios[chrom]:
                if other_element != element:
                    if l_dict[other_element] < r_min_dict[element]:
                        r_min_dict[element] = l_dict[other_element]
    return r_min_dict

def sample_small_l(l_dict, l_min):
    potential_small_l = set()
    for element in l_dict:
        if l_dict[element] == l_min:
            potential_small_l.add(element)
    return random.sample(potential_small_l, 1)[0]

def implement_fission_1(ios, element, fission_log, i):
    for chrom in ios:
        if element in list(ios[chrom].keys()):
            ios[chrom + "_a"][element] = ios[chrom][element]
            for other_elements in list(ios[chrom].keys()):
                if other_elements != element:
                    ios[chrom + "_b"][other_elements] = ios[chrom][other_elements]
            fission_log.append(["parent_node", i, "fission", list([chrom, len(ios[chrom].keys())]), \
                get_all_LMSs_of_a_chrom(ios, chrom)])
            ios.pop(chrom, None)
            return ios, fission_log

def check_l_2_status(ios, element):
    chrom_lengths = [0, 0]
    chrom_content = []
    for chrom in ios:
        if element in list(ios[chrom].keys()) and chrom_lengths[0] == 0:
            chrom_lengths[0] = len(ios[chrom].keys())
            chrom_content.append(chrom)
        elif element in list(ios[chrom].keys()) and chrom_lengths[0] != 0:
            chrom_lengths[1] = len(ios[chrom].keys())
            chrom_content.append(chrom)
    if chrom_lengths[0] > 1 and chrom_lengths[1] > 1:
        return "translocation", chrom_content
    else:
        return "fusion", chrom_content

def implement_fusion_2(ios, chrom_content, fusion_log, i):
    ios[chrom_content[0] + "_" + chrom_content[1].split("_", 1)[1]] = \
    get_union_of_chrom_dicts(ios[chrom_content[0]], ios[chrom_content[1]])
    fusion_log.append(["parent_node", i, "fusion", list(chrom_content), 
        list([get_all_LMSs_of_a_chrom(ios, chrom_content[0]), get_all_LMSs_of_a_chrom(ios, chrom_content[1])])])
    ios.pop(chrom_content[0], None)
    ios.pop(chrom_content[1], None)
    return ios, fusion_log

def implement_translocation_2(ios, chrom_content, element, translocation_log, i):
    ios[chrom_content[0] + "_a_" + chrom_content[1].split("_", 1)[1] + "_a"][element] = \
    ios[chrom_content[0]][element].union(ios[chrom_content[1]][element])
    for index in 0, 1:
        for other_elements in list(ios[chrom_content[index]].keys()):
            if other_elements != element:
                ios[chrom_content[0] + "_b_" + chrom_content[1].split("_", 1)[1] + "_b"][other_elements] = \
                ios[chrom_content[0] + "_b_" + chrom_content[1].split("_", 1)[1] + "_b"][other_elements].union(ios[chrom_content[index]][other_elements])
    translocation_log.append(["parent_node", i, "translocation", list(chrom_content), 
        list([get_all_LMSs_of_a_chrom(ios, chrom_content[0]), get_all_LMSs_of_a_chrom(ios, chrom_content[1])])])
    ios.pop(chrom_content[0], None)
    ios.pop(chrom_content[1], None)
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
    for chrom in ios:
        if element in list(ios[chrom].keys()):
            fusion_candidates.append(chrom)
    best_intersection = float(0)
    best_pair = []
    for combo in itertools.combinations(fusion_candidates, 2):
        if len(set(ios[combo[0]].keys()).intersection(set(ios[combo[1]].keys()))) > best_intersection:
                best_intersection = len(set(ios[combo[0]].keys()).intersection(set(ios[combo[1]].keys())))
                best_pair = list(combo)
    ios[best_pair[0] + "_" + best_pair[1].split("_", 1)[1]] = get_union_of_chrom_dicts(ios[best_pair[0]], ios[best_pair[1]])
    fusion_log.append(["parent_node", i, "fusion", list(best_pair), list([get_all_LMSs_of_a_chrom(ios, best_pair[0]), 
                    get_all_LMSs_of_a_chrom(ios, best_pair[1])])])
    ios.pop(best_pair[0], None)
    ios.pop(best_pair[1], None)
    return ios, fusion_log

# need to QC all implement_* functions
# need to double check the Ferretti et al. paper to see how my algorithm differs

#############################################################################################
#############################################################################################
#############################################################################################

def map_log(log, reference, syngraph, minimum):
    for i in range(1, len(log)): # start at 1 to avoid parsing the header line
        if log[i][2] == "fission":
            chroms = collections.defaultdict(int)
            for marker in log[i][4]:
                if reference in syngraph.nodes[marker]['seqs_by_taxon'].keys():
                    chroms[syngraph.nodes[marker]['seqs_by_taxon'][reference]] += 1
            chroms = pop_bad_mappings(chroms, minimum)
            log[i][4] = list(chroms.keys())      
        else:
            for anc_chrom in 0, 1:
                chroms = collections.defaultdict(int)
                for marker in log[i][4][anc_chrom]:
                    if reference in syngraph.nodes[marker]['seqs_by_taxon'].keys():
                        chroms[syngraph.nodes[marker]['seqs_by_taxon'][reference]] += 1
                chroms = pop_bad_mappings(chroms, minimum)
                log[i][4][anc_chrom] = list(chroms.keys())
    return log

def pop_bad_mappings(chroms, minimum):
    bad_mappings = []
    for mapping in chroms:
        if chroms[mapping] < minimum:
            bad_mappings.append(mapping)
    for bad_mapping in bad_mappings:
        chroms.pop(bad_mapping, None)
    return chroms

def clusters_by_descent(log, tree):
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
            formatted_clusters.append(["cluster_" + str(cluster_ID), taxon])
    return formatted_clusters
                    


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