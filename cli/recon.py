"""

Usage: syngraph recon -g <FILE> -t <FILE> -n <STR,STR> [-o <STR> -h]

  [Options]
    -g, --syngraph <FILE>                       Syngraph file
    -t, --tree <FILE>                           Tree in Newick format
    -n, --node <STR,STR>                        Target tree node for reconstruction, specified as MRCA of taxa
                                                    e.g. -n taxon_A,taxon_B
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
        self.outprefix = args['--outprefix']
        self.tree = self._get_tree(args['--tree'])     # ete3.coretype.tree.TreeNode 
        self.node = self._get_mrca(args['--node'])     # list of nodes

    def _get_tree(self, tree_f):
        tree = ete3.Tree(str(self._get_path(tree_f)))
        for idx, node in enumerate(tree.traverse()):
            if not node.is_leaf():
                node.name = "n%s" % idx
        return tree

    def _get_mrca(self, taxa_string):
        return self.tree.get_common_ancestor(taxa_string.split(",")) # ete3.coretype.tree.TreeNode

    def _get_path(self, infile):
        path = pathlib.Path(infile).resolve()
        if not path.exists():
            sys.exit("[X] File not found: %r" % str(infile))
        return path

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
        print("[+] Reconstructing syngraphs and linkage groups at internal nodes of the following tree:\n%s" % (parameterObj.tree.get_ascii(show_internal=True)))
        reconstructed_syngraphs_by_tree_node = sg.reconstruct_syngraphs_by_tree_node(syngraph, parameterObj.tree)
        # have not tested further ...
        #reconstructed_linkage_groups = sg.reconstruct_linkage_groups_for_each_tree_node(syngraph, parameterObj.tree, algorithm='fitch')
        print("[+] Save Syngraph to file ...")
        #graph_file = reconstructed_syngraph.save(parameterObj, check_consistency=True)
        #print("[+] Saved Syngraph in %r" % graph_file)

        print("[*] Total runtime: %.3fs" % (timer() - main_time))
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()