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

        def get_closest_outgroup(tree, tree_node, child_1, child_2, available_taxa):
            closest_taxon_so_far = "null"
            closest_distance_so_far = float("inf")
            for some_tree_node in tree.search_nodes():
                if some_tree_node.name in available_taxa:
                    if not some_tree_node.name in [descendant.name for descendant in tree_node.get_descendants()]:
                        if tree_node.get_distance(some_tree_node.name) < closest_distance_so_far:
                            closest_distance_so_far = tree_node.get_distance(some_tree_node.name)
                            closest_taxon_so_far = some_tree_node.name
            return(closest_taxon_so_far)

        # first traversal forms triplets from the two children and an outgroup

        available_taxa = set()
        for leaf in parameterObj.tree.get_leaves():
            available_taxa.add(leaf.name)

        for tree_node in parameterObj.tree.traverse(strategy='postorder'):
            if not tree_node.is_leaf() and not tree_node.is_root():
                child_1 = tree_node.get_children()[0].name
                child_2 = tree_node.get_children()[1].name
                outgroup = get_closest_outgroup(parameterObj.tree, tree_node, child_1, child_2, available_taxa)
                print("[+] Inferring median genome for {} using data from {}, {}, and {} ...". format(tree_node.name, child_1, child_2, outgroup))
                syngraph = sg.median_genome(tree_node.name, syngraph, child_1, child_2, outgroup, parameterObj.minimum)
                available_taxa.add(tree_node.name)

        # second traversal forms triplets from the two children and the parent
        # can be done recursively

        for tree_node in parameterObj.tree.traverse(strategy='postorder'):
            if not tree_node.is_leaf() and not tree_node.is_root():
                # remove that tree_node from the current syngraph?
                # or something cleverer?
                child_1 = tree_node.get_children()[0].name
                child_2 = tree_node.get_children()[1].name
                parent = tree_node.up.name
                # what about for children of the root? Should probably just use an outgroup here
                pass



        print("[*] Total runtime: %.3fs" % (timer() - main_time))
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()