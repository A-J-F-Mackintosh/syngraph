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

from operator import attrgetter
import pandas as pd

'''
[To Do]
- Dom:
    - orthofinder parsing has to be implemented
    - write LG plotting function, restructure plotting in general (legends!)
    - all other trello things

- Alex:
    - test for consistency of taxon-names in tree/filenames in parameterObj
    - Ask pablo for BUSCOs, test, let him now how to run it himself
    - Test signed
    - metrics: 
        - marker-order-length distributions by tree-node
        - placements of distinct edges of the same nodes along branches 
            - basis for event detection (?)
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

def reconstruct_linkage_groups_for_each_tree_node(syngraph, tree, algorithm='fitch'):
    '''
    - input: syngraph, tree
    - output: for each internal tree node, a list of sets with each set being markers in a linkage group
    ''' 
    edges_by_tree_node = collections.defaultdict(list)
    coocurence_by_taxon_by_marker_set = collections.defaultdict(lambda: collections.defaultdict(lambda: {False})) # nested dict, set_of_2_markers --> taxon --> True/False
    count = 0 
    if algorithm == 'fitch':
        print("[+] Collecting coocurence data from taxa")
        for taxon in tqdm(syngraph.graph['taxa'], total=len(syngraph.graph['taxa']), desc="[%] ", ncols=100):
            for seq_id, syntenic_markers in syngraph.graph['marker_ids_by_seq_id_by_taxon'][taxon].items():
                for marker_set in itertools.combinations(syntenic_markers, 2):
                    count += 1
                    coocurence_by_taxon_by_marker_set[frozenset(marker_set)][taxon] = {True}
        print("[+] Estimating coocurence data for ancestral genomes")
        for marker_set in tqdm(coocurence_by_taxon_by_marker_set, total=len(coocurence_by_taxon_by_marker_set), desc="[%] ", ncols=100):
            coocurence_by_taxon_by_marker_set[marker_set] = fitch(coocurence_by_taxon_by_marker_set[marker_set], 1, tree)
            # coocurence_by_taxon_by_marker_set is now a nested dict that tells you for a given pair of markers and an internal node
            # whether the two markers are syntenic
            u, v = tuple(marker_set)
            for taxon, coocurence in coocurence_by_taxon_by_marker_set[marker_set].items():
                if coocurence == {True}:
                    edges_by_tree_node[taxon].append((u, v, {'taxa': taxon}))
        LG_by_tree_node = {}
        for tree_node, edges in edges_by_tree_node.items():
            LG_by_tree_node[tree_node] = Syngraph()
            LG_by_tree_node[tree_node].from_edges(edges, taxa=set(tree_node))
            #LG_by_tree_node[tree_node].plot(outprefix="LG_node_%s" % tree_node)
            LG_by_tree_node[tree_node].show_recon_metrics(False, "LG_"+tree_node)

def get_hex_colours_by_taxon(taxa, cmap='Spectral'):
    return {taxon: matplotlib.colors.rgb2hex(cm.get_cmap(cmap, len(taxa))(i)[:3]) for i, taxon in enumerate(sorted(taxa))}


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
        taxa_per_edge = {}
        edges_per_node = {}
        for i in range(1, len(self.graph['taxa'])+1):
            taxa_per_edge[i] = sum([1 for (_, __, taxa) in self.edges.data('taxa') if len(taxa) == i])
        for j in range(1, (len(self.graph['taxa'])+1)*2):
            edges_per_node[j] = 0
        for graph_node_id in self.nodes:
            neighbours = self.degree(graph_node_id)
            edges_per_node[neighbours] += 1
        print("[=] ====================================")
        print("[=] Taxa = %s" % taxon_count)
        print("[=] Nodes (Markers) = %s" % node_total_count)
        print("[=] Edges (Adjacencies)") 
        print("[=]   distinct = %s" % edge_total_count)
        print("[=]   taxa per edge\tcount")
        for key in taxa_per_edge:
            print("[=]   {}\t{}".format(key, taxa_per_edge[key]))
        print("[=]   edges per node\tcount")
        for key in edges_per_node:
            print("[=]   {}\t{}".format(key, edges_per_node[key]))
        print("[=] Subgraphs (connected components) = %s" % connected_component_count)
        print("[=] ====================================")

    # maybe this could be combined with the above?
    def show_recon_metrics(self, gene_order, name):
        connected_component_count = nx.number_connected_components(self)
        chromosome_lengths = []
        for component in nx.connected_components(self):
            chromosome_lengths.append(len(component))
        chromosome_lengths.sort()
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
        print("[=] Subgraph lengths = %s" % chromosome_lengths)
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