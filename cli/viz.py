"""

Usage: syngraph viz -g <FILE> -t <STR> [-o <STR> -h]

  [Options]
    -g, --syngraph <FILE>                       Syngraph file
    -o, --outprefix <STR>                       Outprefix [default: test]
    -t, --taxa <STR>                            Comma delimited list of taxa to visualise
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
        
    def _get_path(self, infile):
        path = pathlib.Path(infile).resolve()
        if not path.exists():
            sys.exit("[X] File not found: %r" % str(infile))
        return str(path)

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
        print(parameterObj.taxa)
        LMSs, unassingables = sg.get_LMSs(syngraph, parameterObj.taxa, 5)
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
        colours = ["#dfac53", "#dfac53", "#dfac53", "#000000"]
        colours_plotable = {}
        for i, taxon in enumerate(parameterObj.taxa):
            colours_plotable[taxon] = colours[i]
        g = graph_tool.Graph(directed=False)
        g.properties[("v", "size")] = g.new_vertex_property("int")
        g.properties[("v", "name")] = g.new_vertex_property("string")
        g.properties[("e", "taxon")] = g.new_edge_property("string")
        g.properties[("e", "colour")] = g.new_edge_property("string")
        for LMS in LMS_plotable:
            v = g.add_vertex()
            g.properties[("v", "name")][v] = LMS
            g.properties[("v", "size")][v] = 10 + (LMS_plotable[LMS] * 0.05)
        for taxon in taxa_plotable:
            for chromosome in taxa_plotable[taxon]:
                if len(chromosome) > 1:
                    for combo in itertools.combinations(chromosome, 2):
                        LMS_1 = combo[0]
                        LMS_2 = combo[1]
                        for v_1 in g.vertices():
                            if g.properties[("v", "name")][v_1] == LMS_1:
                                for v_2 in g.vertices():
                                    if g.properties[("v", "name")][v_2] == LMS_2:
                                        e = g.add_edge(v_1, v_2)
                                        g.properties[("e", "taxon")][e] = taxon
                                        g.properties[("e", "colour")][e] = colours_plotable[taxon]
        edge_dict = collections.defaultdict(set)
        summary_dict = collections.defaultdict(int)
        for e in g.edges():
            edge_dict[frozenset({str(e.source()), str(e.target())})].add(g.properties[("e", "taxon")][e])
        for entry in edge_dict:
            summary_dict[frozenset(edge_dict[entry])] += 1
        for entry in summary_dict:
            print(entry, summary_dict[entry])
        g.list_properties()
        graph_draw(g, vertex_size=g.properties[("v", "size")], edge_color=g.properties[("e", "colour")], output="test.png", 
            edge_pen_width=3, vertex_color="#000000", vertex_fill_color="#ffffff", pos=sfdp_layout(g))
        print("[+] Saved plot in %r" % parameterObj.outprefix)
        print("[*] Total runtime: %.3fs" % (timer() - main_time))
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()
