"""

Usage: syngraph query -g <FILE> [-s <STR> -b -o <STR> -h]

  [Options]
    -g, --syngraph <FILE>          Syngraph file
    -f, --focal_taxon <STR>        Taxon in graph to query
    -r, --reference_taxon <STR>    For each marker in focal_taxon, write which seq the marker is on in this taxon
    -o, --outprefix <STR>          Outprefix [default: test]
    -h, --help                     Show this message

"""

import os
import sys
from docopt import docopt
import pathlib
import collections
from timeit import default_timer as timer
from source import syngraph as sg
import pandas as pd

class ParameterObj():
    def __init__(self, args):
        self.syngraph = self._get_path(args['--syngraph'])
        self.focal_taxon = args['--focal_taxon']
        self.reference_taxon = args['--reference_taxon']
        self.outprefix = args['--outprefix']

    def _get_path(self, infile):
        path = pathlib.Path(infile).resolve()
        if not path.exists():
            sys.exit("[X] File not found: %r" % str(infile))
        return path

def collect_info_from_graph(syngraph, focal, reference):
    info = []
    for graph_node_id in syngraph.nodes:
        entry = [0, 0, 0, 0, 0]
        for graph_node_taxon in syngraph.nodes[graph_node_id]['taxa']:
            if focal in taxa:
                entry[0] = syngraph.nodes[graph_node_id]['seqs_by_taxon'][focal]
                # need coords!
                
    return info

def main(run_params):
    try:
        main_time = timer()
        args = docopt(__doc__)
        print("[+] Sorting out commandline arguments ...")
        parameterObj = ParameterObj(args)  
        print("[+] Creating Syngraph from file ...")
        syngraph = sg.Syngraph()
        syngraph.from_file(parameterObj.syngraph)
        print("[+] Show Syngraph metrics ...")
        syngraph.show_metrics()

        info = collect_info_from_graph(syngraph, parameterObj.focal_taxon, parameterObj.reference_taxon)
        info_df = pd.DataFrame(info)
        info_csv = info_df.to_csv(sep="\t", header=None, index=None)
        with open("{}.query.tsv".format(parameterObj.outprefix), 'w') as fh:
            fh.write(info_csv)
            fh.write("\n")

        print("[*] Total runtime: %.3fs" % (timer() - main_time))
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()