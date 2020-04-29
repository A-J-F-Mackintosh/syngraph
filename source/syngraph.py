import sys
import ete3
import matplotlib
import matplotlib.cm as cm
import pandas as pd
import networkx as nx
import collections
pd.set_option('display.max_rows', None)

def reconstruct_fitch(syngraph, tree, node):
    '''
    - input: syngraph, tree, target_tree_node
    - output: novel Graph with fitch-parsimonious edges for a given tree node
    '''
    edge_sets_by_taxon = syngraph.get_edge_sets_by_taxon() 
    edge_sets_by_tree_node = {}
    for tree_node in tree.traverse(strategy='postorder'):
        if tree_node.name in edge_sets_by_taxon:
            edge_sets_by_tree_node[tree_node.name] = edge_sets_by_taxon[tree_node.name]
        else:
            intersection = set.intersection(
                *[edge_sets_by_tree_node[child_node.name] for child_node in tree_node.get_children()]
                )
            if len(intersection) >= 2:
                edge_sets_by_tree_node[tree_node.name] = intersection
            else:
                edge_sets_by_tree_node[tree_node.name] = set.union(
                    *[edge_sets_by_tree_node[child_node.name] for child_node in tree_node.get_children()]
                    )   
    for tree_node in tree.traverse(strategy='levelorder'):
        if not tree_node.is_root():
            parent_tree_node = tree_node.up
            intersection = edge_sets_by_tree_node[parent_tree_node.name].intersection(edge_sets_by_tree_node[tree_node.name])
            if len(intersection) >= 2:
                edge_sets_by_tree_node[tree_node.name] = intersection
            else:
                edge_sets_by_tree_node[tree_node.name] = edge_sets_by_tree_node[parent_tree_node.name].union(edge_sets_by_tree_node[tree_node.name])

    edges = [(u,v) for u,v in edge_sets_by_tree_node[node.name]]
    taxa = set([leaf.name for leaf in node.get_leaves()])
    reconstructed_syngraph = Syngraph()
    reconstructed_syngraph.from_edges(edges, taxa)
    reconstructed_syngraph.show_metrics()
    reconstructed_syngraph.plot('recon.pdf')
    return reconstructed_syngraph

def get_hex_colours_by_taxon(taxa, cmap='Spectral'):
    return {taxon: matplotlib.colors.rgb2hex(cm.get_cmap(cmap, len(taxa))(i)[:3]) for i, taxon in enumerate(sorted(taxa))}

def cut_links(edges, syngraph, tree, node):
    #####################################################
    # good reconstructions should have two edges per node
    false_links = []
    if len(edges) == 2:
        pass
    # bad reconstructions have > two edges and therefore uncertainty
    else:
        for edge in edges:
            # check whether edges link nodes with the same seqs_by_taxon
            # this can't be done for edge linking to terminal nodes
            if "Terminal_" in edge[0] or "Terminal_" in edge[1]:
                pass
            else:
                if syngraph.graph.nodes[edge[0]]['seqs_by_taxon'] != syngraph.graph.nodes[edge[1]]['seqs_by_taxon']:
                    a = syngraph.graph.nodes[edge[0]]['seqs_by_taxon']
                    b = syngraph.graph.nodes[edge[1]]['seqs_by_taxon']
                    # links is a set of taxa for which these two nodes are on different scaffolds/chromosomes
                    links = {i for i, j in zip(a.keys(), b.keys()) if a[i] != b[j]}
                    ########################################################################################
                    # now reconstruct whether the ancestral genome has these two node on the same chromosome
                    # "1" = same chromosome, "0" = different chromosome
                    reconstruct_links = {}
                    target = node
                    # go up the tree
                    for tree_node in tree.traverse(strategy='postorder'):
                        if not tree_node.is_leaf():
                            # get info about tree_node
                            leaves = frozenset(tree_node.get_leaves())
                            children = tree_node.get_children()
                            # collect child states 
                            child_counter = 0
                            child_links_1 = set()
                            child_links_2 = set()
                            for child in children:
                                child_counter += 1
                                # if child is a leaf then collect edge state from c
                                if len(child.get_leaves()) == 1:
                                    taxon = str(child)[3:]
                                    if taxon in links:
                                        if child_counter == 1:
                                            child_links_1.add("0")
                                        elif child_counter == 2:
                                            child_links_2.add("0")
                                    else:
                                        if child_counter == 1:
                                            child_links_1.add("1")
                                        elif child_counter == 2:
                                            child_links_2.add("1")
                                # if child is an internal node then look in the reconstruct_links for link state
                                else:
                                    if child_counter == 1:
                                        child_links_1 = reconstruct_links[frozenset(child.get_leaves())]
                                    elif child_counter == 2:
                                        child_links_2 = reconstruct_links[frozenset(child.get_leaves())]
                            # reconstruct current link state from children
                            if len(child_links_1.intersection(child_links_2)) == 1:
                                reconstruct_links[leaves] = child_links_1.intersection(child_links_2)
                            else:
                                reconstruct_links[leaves] = child_links_1.union(child_links_2)
                    # now go down the tree
                    for tree_node in tree.traverse(strategy='preorder'):
                        if not tree_node.is_leaf():
                            tree_node_links = reconstruct_links[frozenset(tree_node.get_leaves())]
                            children = tree_node.get_children()
                            child_counter = 0
                            for child in children:
                                if not child.is_leaf():
                                    child_counter += 1
                                    if child_counter == 1:
                                        child_links_1 = reconstruct_links[frozenset(child.get_leaves())]
                                        if len(tree_node_links.intersection(child_links_1)) == 1:
                                            reconstruct_links[frozenset(child.get_leaves())] = tree_node_links.intersection(child_links_1)
                                    elif child_counter == 2:
                                        child_links_2 =  reconstruct_links[frozenset(child.get_leaves())]
                                        if len(tree_node_links.intersection(child_links_2)) == 1:
                                            reconstruct_links[frozenset(child.get_leaves())] = tree_node_links.intersection(child_links_2)
                    #############################################################################################################
                    # if the ancestral node has these graph nodes on different chromosomes, then the edge represents a false link 
                    if set("0") == reconstruct_links[frozenset(target.get_leaves())]:
                        false_links.append(edge)
    return(false_links)

#############################################################################################

class Syngraph(nx.Graph):
    def __init__(self, name='', **attr):
        nx.Graph.__init__(self, name=name, taxa=set(), terminal_nodes=set(), **attr)
        
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
                self.add_node(markerObj.name, taxa=set(), terminal=set())               
            # add taxon to node (useful if '--missing')
            self.nodes[markerObj.name]['taxa'].add(markerObj.taxon)
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
                    self.graph['terminal_nodes'].add(prev_markerObj.name)
                self.nodes[markerObj.name]['terminal'].add(markerObj.taxon)
                self.graph['terminal_nodes'].add(markerObj.name)
            prev_markerObj = markerObj
        self.graph['terminal_nodes'].add(markerObj.name)
        self.nodes[markerObj.name]['terminal'].add(markerObj.taxon)

    def get_edge_sets_by_taxon(self):
        edge_sets_by_taxon = collections.defaultdict(set)
        for u, v, taxa in self.edges.data('taxa'):
            for taxon in taxa:
                edge_sets_by_taxon[taxon].add(frozenset((u, v)))
        return edge_sets_by_taxon

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