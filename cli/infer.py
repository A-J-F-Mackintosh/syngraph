"""

Usage: syngraph infer -g <FILE> -t <NWK> (-s <STR> | -S <STR>) [-m <INT> -r <STR> -o <STR> -h]

  [Options]
    -g, --syngraph <FILE>        Syngraph file
    -t, --tree <NWK>             Tree in Newick format
    -m, --minimum <INT>          Minimum number of markers for a synteny relationship [default: 5]
    -r, --rearrangements <INT>   Rearrangements to be modelled: 2 (fission+fusion) or 3 (fission+fusion+translocation) [default: 2]
    -s, --reference_taxon <STR>  Taxon name to map ancestral seqs to
    -S, --reference_tsv <STR>    Predefined linkage groups in tsv format (marker_ID, LG_name) to map ancestral seqs to
    -o, --outprefix <STR>        Outprefix [default: test]
    -h, --help                   Show this message

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
import pandas as pd
import collections


class ParameterObj():
    def __init__(self, args):
        self.syngraph = self._get_path(args['--syngraph'])
        self.tree = self._get_tree(args['--tree'])
        self.minimum = int(args['--minimum'])
        self.model = self._check_model(int(args['--rearrangements']))
        self.outprefix = args['--outprefix']
        if args['--reference_taxon']:
            self.reference_taxon = args['--reference_taxon']
            self.reference_tsv = None
        elif args['--reference_tsv']:
            self.reference_tsv = self._get_path(args['--reference_tsv'])
            self.reference_taxon = None

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

    def _check_model(self, model):
        if model not in [2, 3]:
            sys.exit("[X] Invalid model specified: %r" % str(model))
        return model

def check_tree_syngraph_concordance(tree, syngraph):
    for leaf in tree.get_leaves():
        if leaf.name not in syngraph.graph['taxa']:
            sys.exit("\n[X] A leaf in the tree, {}, is not in the syngraph.\n".format(leaf.name))

def check_reference_syngraph_concordance(reference, syngraph):
    if reference not in syngraph.graph['taxa']:
        print("\n[X] WARNING: The reference taxon, {}, is not in the syngraph.\n".format(reference))

def generate_reference_dict(reference_path):
    reference_dict = {}
    with open(reference_path, "r") as reference_f:
        for line in reference_f:
            marker_ID = line.rstrip().split()[0]
            ref_chrom = line.rstrip().split()[1]
            reference_dict[marker_ID] = ref_chrom
    return reference_dict

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
        check_tree_syngraph_concordance(parameterObj.tree, syngraph)
        print("[+] Show Syngraph metrics ...")
        syngraph.show_metrics()
        random.seed(44) 

        solved_syngraph, log = sg.tree_traversal(syngraph, parameterObj)
        if parameterObj.reference_taxon:
            check_reference_syngraph_concordance(parameterObj.reference_taxon, solved_syngraph)
            mapped_log = sg.map_log(log, parameterObj.reference_taxon, None, 
                solved_syngraph, parameterObj.minimum)
        else:
            reference_dict = generate_reference_dict(parameterObj.reference_tsv)
            mapped_log = sg.map_log(log, parameterObj.reference_taxon, reference_dict, 
                solved_syngraph, parameterObj.minimum)

        clusters = sg.clusters_by_descent(log, parameterObj.tree, solved_syngraph)

        mapped_log_df = pd.DataFrame(mapped_log)
        clusters_df = pd.DataFrame(clusters)
        mapped_log_tsv = mapped_log_df.to_csv(sep="\t", header=None, index=None)
        clusters_tsv = clusters_df.to_csv(sep="\t", header=None, index=None)

        with open("{}.rearrangements.tsv".format(parameterObj.outprefix), 'w') as fh:
            fh.write(mapped_log_tsv)
            fh.write("\n")
        with open("{}.clusters.tsv".format(parameterObj.outprefix), 'w') as fh:
            fh.write(clusters_tsv)
            fh.write("\n")
        with open("{}.newick.ascii".format(parameterObj.outprefix), 'w') as fh:
            fh.write(parameterObj.tree.get_ascii())
            fh.write("\n")
        with open("{}.newick.txt".format(parameterObj.outprefix), 'w') as fh:
            fh.write(parameterObj.tree.write(format=1))
            fh.write("\n")

        print("[+] Save Syngraph to file ...")
        graph_file = solved_syngraph.save(parameterObj, check_consistency=False, with_ancestors=True)
        print("[+] Saved Syngraph in %r" % graph_file)


        print("[*] Total runtime: %.3fs" % (timer() - main_time))
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)

if __name__ == '__main__':
    main()


# label branches on tree, output as ascii
# can split tsv up:
#
# event summary:
#   branch, event_type,    event_ID,
# 
# event details:   
#   event_ID,    LMS representation before/after event (backwards in time)
#
# starting point:
#   taxon, LMS representation