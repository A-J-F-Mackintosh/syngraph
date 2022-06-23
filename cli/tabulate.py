"""

Usage: syngraph tabulate -g <FILE> [-o <STR> -h]

  [Options]
    -g, --syngraph <FILE>          Syngraph file
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
        self.outprefix = args['--outprefix']

    def _get_path(self, infile):
        path = pathlib.Path(infile).resolve()
        if not path.exists():
            sys.exit("[X] File not found: %r" % str(infile))
        return path

def collect_info_from_graph(syngraph):
    info = []
    for graph_node_id in syngraph.nodes:
        entry = [graph_node_id]
        for taxon in syngraph.graph['taxa']:
            if taxon in syngraph.nodes[graph_node_id]['taxa']:
                entry.append(syngraph.nodes[graph_node_id]['seqs_by_taxon'][taxon])
                if taxon in syngraph.nodes[graph_node_id]['starts_by_taxon']:
                    entry.append(syngraph.nodes[graph_node_id]['starts_by_taxon'][taxon])
                else:
                    entry.append("NA")
                if taxon in syngraph.nodes[graph_node_id]['ends_by_taxon']:
                    entry.append(syngraph.nodes[graph_node_id]['ends_by_taxon'][taxon])
                else:
                    entry.append("NA")
            else:
                entry.append("NA")
                entry.append("NA")
            info.append(entry)
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

        info = collect_info_from_graph(syngraph)
        info_df = pd.DataFrame(info)
        info_csv = info_df.to_csv(sep="\t", header=None, index=None)
        with open("{}.table.tsv".format(parameterObj.outprefix), 'w') as fh:
            fh.write(info_csv)
            fh.write("\n")

        print("[*] Total runtime: %.3fs" % (timer() - main_time))
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()