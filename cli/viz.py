"""

Usage: syngraph viz -g <FILE> -t <STR> -n <INT> [-o <STR> -s <STR> -m <INT> -h]

  [Options]
    -g, --syngraph <FILE>                       Syngraph file
    -o, --outprefix <STR>                       Outprefix [default: test]
    -t, --taxa <STR>                            Comma delimited list of taxa to visualise, or 'all'
    -n, --threshold <INT>                       Minimum number of supporting taxa required for an edge to be drawn
    -s, --size <STR>                            Node size in the form 'a,b' where size = a + (b * markers) [default: 10,0.1]
    -m, --minimum <INT>                         Minimum synteny set size [default: 2]
    -h, --help                                  Show this screen.

"""

import sys
import os
from docopt import docopt
import pathlib
from timeit import default_timer as timer
from source import syngraph as sg
import collections
from graph_tool.all import *
import itertools
from matplotlib import cm

class ParameterObj():
    def __init__(self, args):
        self.syngraph = self._get_path(args['--syngraph'])
        self.outprefix = args['--outprefix']
        self.taxa = args['--taxa'].split(",")
        self.threshold = int(args['--threshold'])
        self.a_size = float(args['--size'].split(",")[0])
        self.b_size = float(args['--size'].split(",")[1])
        self.minimum = int(args['--minimum'])
        
    def _get_path(self, infile):
        path = pathlib.Path(infile).resolve()
        if not path.exists():
            sys.exit("[X] File not found: %r" % str(infile))
        return str(path)

def get_taxa(taxa, syngraph):
    if taxa == ['all']:
        return list(syngraph.graph['taxa'])
    else:
        return taxa


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
        print("[+] Collecting Syngraph info ...")
        parameterObj.taxa = get_taxa(parameterObj.taxa, syngraph)
        print(parameterObj.taxa)
        LMSs, unassingables = sg.get_LMSs(syngraph, parameterObj.taxa, parameterObj.minimum)
        LMS_plotable = {}
        taxa_plotable = {}
        for LMS in LMSs:
            LMS_plotable[LMS] = len(LMSs[LMS])
        #print(LMS_plotable)
        for taxon in parameterObj.taxa:
            genome_object = sg.compact_synteny_1(syngraph, LMSs, taxon)
            taxa_plotable[taxon] = genome_object.CWAL
        #print(taxa_plotable)
        print("[+] Plotting Syngraph ...")
        g = graph_tool.Graph(directed=False)
        g.properties[("v", "size")] = g.new_vertex_property("int")
        g.properties[("v", "name")] = g.new_vertex_property("string")
        g.properties[("e", "taxon")] = g.new_edge_property("string")
        g.properties[("e", "colour")] = g.new_edge_property("string")
        for LMS in LMS_plotable:
            v = g.add_vertex()
            g.properties[("v", "name")][v] = LMS
            g.properties[("v", "size")][v] = parameterObj.a_size + (parameterObj.b_size  * LMS_plotable[LMS])
        edge_support = collections.defaultdict(int)
        for taxon in taxa_plotable:
            print(taxon)
            for chromosome in taxa_plotable[taxon]:
                if len(chromosome) > 1:
                    for combo in itertools.combinations(chromosome, 2):
                        LMS_1 = combo[0]
                        LMS_2 = combo[1]
                        for v_1 in g.vertices():
                            if g.properties[("v", "name")][v_1] == LMS_1:
                                for v_2 in g.vertices():
                                    if g.properties[("v", "name")][v_2] == LMS_2:
                                        edge_support[frozenset([LMS_1, LMS_2])] += 1
                                        if edge_support[frozenset([LMS_1, LMS_2])] == parameterObj.threshold:
                                            e = g.add_edge(v_1, v_2)
                                            g.properties[("e", "taxon")][e] = taxon
                                            g.properties[("e", "colour")][e] = "#000000"
        g.list_properties()
        graph_draw(g, vertex_size=g.properties[("v", "size")], edge_color=g.properties[("e", "colour")], output="test.pdf", 
            edge_pen_width=0.5, vertex_color="#000000", vertex_fill_color="#fff0cb", pos=sfdp_layout(g, r=10))
        print("[+] Saved plot in %r" % parameterObj.outprefix)
        print("[*] Total runtime: %.3fs" % (timer() - main_time))
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()
