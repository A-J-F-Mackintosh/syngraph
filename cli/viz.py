"""

Usage: syngraph viz -g <FILE> [-m -o <STR> -h]

  [Options]
    -g, --syngraph <FILE>                       Syngraph file
    -o, --outprefix <STR>                       Outprefix [default: test]
    -m, --as_multigraph                         Plot as multigraph
    -h, --help                                  Show this screen.

"""

import sys
from docopt import docopt
import pathlib
from timeit import default_timer as timer
from source import syngraph as sg

class ParameterObj():
    def __init__(self, args):
        self.syngraph = self._get_path(args['--syngraph'])
        self.outprefix = args['--outprefix']
        self.as_multigraph = args['--as_multigraph']
        
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
        print("[+] Plotting Syngraph ...")
        plot_file = syngraph.plot(
            parameterObj.outprefix, 
            as_multigraph=parameterObj.as_multigraph)
        print("[+] Saved plot in %r" % plot_file)
        print("[*] Total runtime: %.3fs" % (timer() - main_time))
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()