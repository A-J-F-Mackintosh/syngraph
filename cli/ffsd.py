"""

Usage: syngraph ffsd -g <FILE> -r <STR> -q <STR> [-m <INT> -Q <STR> -o <STR> -h]

  [Options]
    -g, --syngraph <FILE>                       Syngraph file
    -r, --reference <STR>                       Reference taxon, e.g. genus_speciesA
    -q, --query <STR>                           Query taxon, e.g. genus_speciesB
    -Q, --query2 <STR>                          A second query taxon that can be used to calculate the median genome between ref, query, and query2 [default: None]
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
        self.outprefix = args['--outprefix']
        self.reference = args['--reference']
        self.query = args['--query']
        self.query2 = args['--query2']
        self.minimum = int(args['--minimum'])

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

        if parameterObj.query2 == "None":
            instance_of_synteny = sg.compact_synteny(syngraph, parameterObj.reference, parameterObj.query, parameterObj.minimum, "syngraph_syngraph")
            sg.ffsd(instance_of_synteny)
        else:
            sg.median_genome(syngraph, parameterObj.reference, parameterObj.query, parameterObj.query2, parameterObj.minimum)

        print("[*] Total runtime: %.3fs" % (timer() - main_time))
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()