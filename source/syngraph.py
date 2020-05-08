import sys
import ete3
import matplotlib
import matplotlib.cm as cm
import pandas as pd
import networkx as nx
import collections
pd.set_option('display.max_rows', None)

def fitch(states_by_taxon, number_of_states, tree): # states_by_taxon should be a dict, with keys as taxon and states as sets
    states_by_tree_node = {}
    for tree_node in tree.traverse(strategy='postorder'):
        if tree_node.name in states_by_taxon:
            states_by_tree_node[tree_node.name] = states_by_taxon[tree_node.name]
        else:
            intersection = set.intersection(*[states_by_tree_node[child_node.name] for child_node in tree_node.get_children()])
            if len(intersection) >= number_of_states:
                states_by_tree_node[tree_node.name] = intersection
            else:
                states_by_tree_node[tree_node.name] = set.union(
                    *[states_by_tree_node[child_node.name] for child_node in tree_node.get_children()])
    for tree_node in tree.traverse(strategy='levelorder'):
        if not tree_node.is_root():
            parent_tree_node = tree_node.up
            intersection = states_by_tree_node[parent_tree_node.name].intersection(states_by_tree_node[tree_node.name])
            if len(intersection) >= number_of_states:
                states_by_tree_node[tree_node.name] = intersection
    return(states_by_tree_node)

def reconstruct_syngraphs_for_each_tree_node(syngraph, tree, algorithm='fitch'):
    '''
    - input: syngraph, tree
    - output: novel graphs with fitch edges for each internal tree node
    '''
    gene_order_data = collections.defaultdict(dict) # nested dict, graph_node --> taxon --> edges
    if algorithm == 'fitch':
        for graph_node_id in syngraph.nodes:
            edge_sets_by_taxon = syngraph.get_target_edge_sets_by_taxon(graph_node_id)
            gene_order_data[graph_node_id] = fitch(edge_sets_by_taxon, 2, tree)
            #######
            # TBI #
            #######
            # make graphs from edge_lists_by_tree_node

def reconstruct_linkage_groups_for_each_tree_node(syngraph, tree, algorithm='fitch'):
    '''
    - input: syngraph, tree
    - output: for each internal tree node, a list of sets with each set being markers in a linkage group
    '''
    synteny_data = collections.defaultdict(dict) # nested dict, set_of_2_markers --> taxon --> True/False
    if algorithm == 'fitch':
        for i in range(0, len(list(syngraph))): # i and j form unique pairwise comparisons between markers
            for j in range(i+1, len(list(syngraph))):
                set_of_2_markers = {list(syngraph)[i], list(syngraph)[j]}
                for taxon in syngraph.graph['taxa']:
                    synteny = set() # whether markers are syntenic, needs to be a set for fitch algorithm
                    if syngraph.nodes[list(syngraph)[i]]['seqs_by_taxon'][taxon] == \
                    syngraph.nodes[list(syngraph)[j]]['seqs_by_taxon'][taxon]:
                        synteny.add("True")
                    else:
                        synteny.add("False")
                    synteny_data[frozenset(set_of_2_markers)][taxon] = synteny
                # check whether there are any matches, if not then remove markers
                synteny_match = False
                for taxon in syngraph.graph['taxa']:
                    if synteny_data[frozenset(set_of_2_markers)][taxon] == {"True"}:
                        synteny_match = True
                if synteny_match == False:
                    del synteny_data[frozenset(set_of_2_markers)]
        for set_of_2_markers in synteny_data: # loop over keys and fitch reconstruct the internal nodes
            synteny_data[set_of_2_markers] = fitch(synteny_data[frozenset(set_of_2_markers)], 1, tree)
        # synteny_data is now a nested dict that tells you for a given pair of markers and an internal node
        # whether the two markers are syntenic
    #######
    # TBI #
    #######
    # cannot be saved (yet)
    # cannot be plotted (yet)



def get_hex_colours_by_taxon(taxa, cmap='Spectral'):
    return {taxon: matplotlib.colors.rgb2hex(cm.get_cmap(cmap, len(taxa))(i)[:3]) for i, taxon in enumerate(sorted(taxa))}


#############################################################################################

class Syngraph(nx.Graph):
    '''
    [ToDo]
    - write documentation for Syngraph about which data is stored where

    - Cut_links : edge-pruning
        Caution: dom broke ['seqs_by_taxon']
        should work on any edge-attribute  (by taxon) instead
            - 'synteny': intra-sequence > inter-sequence
            - 'distance': shorter_edges > longer_edges 
            ...
    
    '''
    def __init__(self, name='', **attr):
        nx.Graph.__init__(self, name=name, taxa=set(), **attr)
        
    def __repr__(self):
       return "Syngraph(name=%r, taxa=%r, ...)" % (self.name, self.taxa) 

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

    def from_file(self, graph_file):
        g = nx.read_gpickle(graph_file)
        self.add_nodes_from(g.nodes(data=True))
        self.add_edges_from(g.edges(data=True))
        self.graph = g.graph

    def from_edges(self, edges, taxa):
        if isinstance(taxa, str):
            taxa = [taxa] 
        self.add_edges_from(edges)
        for u, v in self.edges:
            self.edges[u, v]['taxa'] = set(taxa)
        for node in self.nodes:
            self.nodes[node]['taxa'] = set(taxa)
        self.graph['taxa'] = set(taxa)

    def from_markerObjs(self, markerObjs):
        prev_markerObj = MarkerObj(None)
        for markerObj in markerObjs:
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

    def get_target_edge_sets_by_taxon(self, graph_node_id):
        target_edge_sets_by_taxon = collections.defaultdict(set)
        for u, v, taxa in self.edges([graph_node_id], data='taxa'):
            for taxon in taxa:
                target_edge_sets_by_taxon[taxon].add(frozenset((u, v)))
        for taxon, target_edge_set in target_edge_sets_by_taxon.items():
            if len(target_edge_set) == 1:
                target_edge_sets_by_taxon[taxon].add(frozenset((graph_node_id, '%s.terminal' % (graph_node_id))))
        return target_edge_sets_by_taxon

    def show_metrics(self):
        taxon_count = len(self.graph['taxa'])
        node_total_count = nx.number_of_nodes(self)
        edge_total_count = nx.number_of_edges(self)
        connected_component_count = nx.number_connected_components(self)
        edge_monomorphic_count = sum([1 for (_, __, taxa) in self.edges.data('taxa') if len(taxa) == taxon_count])
        edge_singleton_count = sum([1 for _, __, taxa in self.edges.data('taxa') if len(taxa) == 1]) 
        print("[=] ====================================")
        print("[=] Taxa = %s" % taxon_count)
        print("[=] Nodes (Markers) = %s" % node_total_count)
        print("[=] Edges (Adjacencies)") 
        print("[=]   distinct = %s" % edge_total_count)
        print("[=]   monomorphic = %s" % edge_monomorphic_count)
        print("[=]   singleton = %s" % edge_singleton_count)
        print("[=] Subgraphs (connected components) = %s" % connected_component_count)
        print("[=] ====================================")
        
    def plot(self, outprefix, cmap='Set2', as_multigraph=True):
        taxon_count = len(self.graph['taxa'])
        colour_by_taxon = get_hex_colours_by_taxon(self.graph['taxa'], cmap=cmap)
        if as_multigraph:
            multi_graph = nx.MultiGraph() # create new multigraph
            for node, data in self.nodes.data(True):
                multi_graph.add_node(node, label='', color='black', width=0.5, shape='point')
                if 'terminal' in data:
                    for taxon in data['terminal']:
                        new_terminal_node = "%s_%s_terminal" % (taxon, node)
                        multi_graph.add_node(new_terminal_node, color=colour_by_taxon[taxon], label='', width=0.5, shape='point')
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
    def __init__(self, name=None, desc=None, status=None, taxon=None, seq=None, start=float("nan"), end=float("nan"), length=float("nan")):
        self.name = name
        self.desc = desc if desc is not None else name
        self.status = status
        self.taxon = taxon
        self.seq = seq
        self.start = start
        self.end = end
        self.length = length

    def __repr__(self):
        return "MarkerObj(name=%r, desc=%r, status=%r, taxon=%r, seq=%r, start=%f, end=%f, length=%f)" % (
            self.name, self.desc, self.status, self.taxon, self.seq, self.start, self.end, self.length) 

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
            if self.start == other.start and self.end == other.end:
                return 0.0
            else:
                return float(max(self.start, other.start) - min(self.end, other.end))                
        return float("nan")