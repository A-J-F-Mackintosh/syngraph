"""

Usage: syngraph ffsd -g <FILE> -t <NWK> [-m <INT> -o <STR> -h]

  [Options]
    -g, --syngraph <FILE>                       Syngraph file
    -t, --tree <NWK>                            Tree in Newick format
    -m, --minimum <INT>                         The minimum number of markers for a synteny relationship [default: 5]
    -o, --outprefix <STR>                       Outprefix [default: test]
    -h, --help                                  Show this screen.

"""

import sys
from docopt import docopt
import pathlib
import ete3
import copy
from timeit import default_timer as timer
from source import syngraph as sg

class ParameterObj():
    def __init__(self, args):
        self.syngraph = self._get_path(args['--syngraph'])
        self.tree = self._get_tree(args['--tree'])
        self.outprefix = args['--outprefix']
        self.minimum = int(args['--minimum'])

    def _get_path(self, infile):
        path = pathlib.Path(infile).resolve()
        if not path.exists():
            sys.exit("[X] File not found: %r" % str(infile))
        return path

    def _get_tree(self, tree_f):
        tree = ete3.Tree(str(self._get_path(tree_f)))
        for idx, node in enumerate(tree.traverse()):
            if not node.is_leaf():
                node.name = "n%s" % idx
        return tree

def main(run_params):
    try:
        main_time = timer()
        args = docopt(__doc__)
        print(args)
        print("[+] Sorting out commandline arguments ...")
        parameterObj = ParameterObj(args)
        print("[+] Creating Syngraph from file ...")
        syngraph = sg.Syngraph()
        syngraph.from_file(parameterObj.syngraph)
        print("[+] Show Syngraph metrics ...")
        syngraph.show_metrics()

        #### TO DO LIST AND QUESTIONS
        ####
        #### For each triplet, pre-calculate from-median probs of fis/fus given branch length and direction
        #### Might be easiest to evaluate overall likelihood/parismony with a distinct function, rather than on the fly
        #### Should the second traversal be recursive?
        #### How easy would including RTs be?
        #### Could adding in unassignables cause problems in later traversals?

        # copy syngraph
        traversal_0_syngraph = copy.deepcopy(syngraph)
        traversal_1_syngraph = copy.deepcopy(syngraph)

        # a function for getting the closest outgroup to a tree node
        def get_closest_outgroup(tree, tree_node, child_1, child_2, available_taxa):
            closest_taxon_so_far = "null"
            closest_distance_so_far = float("inf")
            for some_tree_node in tree.search_nodes():
                if some_tree_node.name in available_taxa:
                    if not some_tree_node.name in [descendant.name for descendant in tree_node.get_descendants()]:
                        if not some_tree_node.name == tree_node.name:
                            if tree_node.get_distance(some_tree_node.name) < closest_distance_so_far:
                                closest_distance_so_far = tree_node.get_distance(some_tree_node.name)
                                closest_taxon_so_far = some_tree_node.name
            return(closest_taxon_so_far)


        # define which taxa are extant and so can be used as outgroups from the start
        available_taxa = set()
        for leaf in parameterObj.tree.get_leaves():
            available_taxa.add(leaf.name)


        # first traversal forms triplets from a node's two children and an outgroup
        # syngraph is updated each iteration
        print("[+] Starting first traversal ...")
        print("[+] ========================================================================")

        for tree_node in parameterObj.tree.traverse(strategy='postorder'):
            if not tree_node.is_leaf() and not tree_node.is_root():
                child_1 = tree_node.get_children()[0].name
                child_2 = tree_node.get_children()[1].name
                outgroup = get_closest_outgroup(parameterObj.tree, tree_node, child_1, child_2, available_taxa)
                print("[+] Inferring median genome for {} using data from {}, {}, and {} ...". format(tree_node.name, child_1, child_2, outgroup))
                traversal_0_syngraph = sg.median_genome(tree_node.name, traversal_0_syngraph, traversal_0_syngraph, child_1, child_2, outgroup, parameterObj.minimum)
                available_taxa.add(tree_node.name)

        print("[=] ========================================================================")

        # second traversal forms triplets from the two children and the parent, as we now have the parents from the first traversal
        # if a node is a child of the root then the triplet is formed from the node's two children and an outgroup, i.e. the other child of the root
        # info is read from the first traversals syngraoh but a new syngraph is written to
        print("[+] Starting second traversal ...")
        print("[+] ========================================================================")

        for tree_node in parameterObj.tree.traverse(strategy='postorder'):
            if not tree_node.is_leaf() and not tree_node.is_root():
                child_1 = tree_node.get_children()[0].name
                child_2 = tree_node.get_children()[1].name
                if tree_node.up.is_root():
                    outgroup = get_closest_outgroup(parameterObj.tree, tree_node, child_1, child_2, available_taxa)
                    print("[+] Inferring median genome for {} using data from {}, {}, and {} ...". format(tree_node.name, child_1, child_2, outgroup))
                    traversal_1_syngraph = sg.median_genome(tree_node.name, traversal_0_syngraph, traversal_1_syngraph, child_1, child_2, outgroup, parameterObj.minimum)
                else:
                    parent = tree_node.up.name
                    print("[+] Inferring median genome for {} using data from {}, {}, and {} ...". format(tree_node.name, child_1, child_2, parent))
                    traversal_1_syngraph = sg.median_genome(tree_node.name, traversal_0_syngraph, traversal_1_syngraph, child_1, child_2, parent, parameterObj.minimum)

        print("[=] ========================================================================")

        # evaluator should iterate from the children of the root to the leaves
        print("[=] ========================================================================")

        for tree_node in parameterObj.tree.traverse(strategy='preorder'):
            if not tree_node.is_leaf() and not tree_node.is_root():
                child_1 = tree_node.get_children()[0].name
                child_2 = tree_node.get_children()[1].name
                for child in child_1, child_2:
                    #print(tree_node.name, child, tree_node.get_distance(child))
                    # need to write a get_LMS_for_pair function, or use the existing function in a hacky way
                    pass


        print("[*] Total runtime: %.3fs" % (timer() - main_time))
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()