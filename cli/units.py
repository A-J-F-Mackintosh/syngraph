"""

Usage: syngraph units -g <FILE> -t <FILE> -n <STR> [-o <STR> -h]

  [Options]
    -g, --syngraph <FILE>                       Syngraph file
    -t, --tree <FILE>                           Tree in Newick format
    -n, --mrca <STR>                            Two taxon names seperated by a comma, e.g. taxon_A,taxon_B, the mrca (internal node) of these taxa will be used to define units
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
        self.tree = self._get_tree(args['--tree'])
        self.mrca = self._get_mrca(args['--mrca'])

    def _get_tree(self, tree_f):
        tree = ete3.Tree(str(self._get_path(tree_f)))
        for idx, node in enumerate(tree.traverse()):
            if not node.is_leaf():
                node.name = "n%s" % idx
        return tree

    def _get_path(self, infile):
        path = pathlib.Path(infile).resolve()
        if not path.exists():
            sys.exit("[X] File not found: %r" % str(infile))
        return path

    def _get_mrca(self, instring):
        mrca = self.tree.get_common_ancestor(instring.split(","))
        return mrca
        

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
        print("[+] Identifying LMS (linked marker sets) within the syngraph")
        LinkedMarkerSets, LinkedMarkerSet_IDs = sg.find_LinkedMarkerSets(syngraph)
        print("[+] Reconstructing linkage groups at internal nodes of the following tree:\n%s\n" % (parameterObj.tree.get_ascii(show_internal=True)))
        LG_store = sg.reconstruct_linkage_groups_for_each_tree_node(syngraph, LinkedMarkerSets, LinkedMarkerSet_IDs, parameterObj.tree, algorithm='fitch')
        print("[+] Etimating chromosomal units and variation")
        sg.analyse_chromosomal_units(syngraph, LG_store, LinkedMarkerSet_IDs, parameterObj.mrca)
        print("[*] Total runtime: %.3fs" % (timer() - main_time))
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()