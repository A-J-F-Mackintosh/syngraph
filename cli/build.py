"""

Usage: syngraph build -d <DIR> [-o <STR> -m -h]

  [Options]
    -d, --dir <DIR>                             Directory containing marker tsvs
    -o, --outprefix <STR>                       Outprefix [default: test]
    -m, --missing                               Allow markers that are missing in ≥1 taxa
    -h, --help                                  Show this screen.

"""

import os
import sys
from docopt import docopt
import pathlib
import collections
from timeit import default_timer as timer
from source import syngraph as sg

class ParameterObj():
    def __init__(self, args):
        print(args)
        self.directory = self._get_path(args['--dir'])
        self.infiles = self._get_infiles(args['--dir'])
        self.outprefix = args['--outprefix']
        self.missing = True if (args['--missing']) else False
        print(self.__dict__)

    def _get_infiles(self, directory):
        infiles = []
        extension = '.tsv'
        for f in os.listdir(directory):
            if f.endswith(extension):
                infiles.append(os.path.join(directory, f))
        if len(infiles) == 0:
            sys.exit("[X] No files ending in '*.tsv' in folder %s" % os.path.abspath(directory))
        return infiles

    def _get_path(self, infile):
        if infile is None:
            return infile
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
        print("[+] Creating Syngraph from markers ...")
        syngraph = sg.Syngraph()
        markerObjs = sg.load_markerObjs(parameterObj)
        syngraph.from_markerObjs(markerObjs)
        print("[+] Show Syngraph metrics ...")
        syngraph.show_metrics()
        print("[+] Save Syngraph to file ...")
        graph_file = syngraph.save(parameterObj, check_consistency=True, with_ancestors=False)
        print("[+] Saved Syngraph in %r" % graph_file)
        print("[*] Total runtime: %.3fs" % (timer() - main_time))
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()
