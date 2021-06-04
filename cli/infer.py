"""

Usage: syngraph infer -g <FILE> -t <NWK> [-m <INT> -o <STR> -h]

  [Options]
    -g, --syngraph <FILE>                       Syngraph file
    -t, --tree <NWK>                            Tree in Newick format
    -m, --minimum <INT>                         Minimum number of markers for a synteny relationship [default: 5]
    -o, --outprefix <STR>                       Outprefix [default: test]
    -h, --help                                  Show this message

"""

import sys
from docopt import docopt
import pathlib
import ete3
import copy
import random
from functools import partial
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
        print(tree.get_ascii())
        print("")
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
        random.seed(44)

        #### TO DO LIST AND QUESTIONS
        ####
        #### Should the traversals be recursive?
        #### How easy would including RTs be?
        #### Could adding in unassignables cause problems in later traversals?

        solved_syngraph, log = sg.tree_traversal(syngraph, parameterObj)
        with open("{}.tsv".format(parameterObj.outprefix), 'w') as fh:
            fh.write("\n".join(log))
            fh.write("\n")
        with open("{}.newick.txt".format(parameterObj.outprefix), 'w') as fh:
            fh.write(parameterObj.tree.get_ascii())
            fh.write("\n")
        print("[*] Total runtime: %.3fs" % (timer() - main_time))
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)

if __name__ == '__main__':
    main()