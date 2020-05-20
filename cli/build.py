"""

Usage: syngraph build -d <DIR> [-l <STR>] [-g <FILE>] [-m|-M] [-f|-F] [-s] [-o <STR> -h]

  [Options]
    -l, --label <STR>                           Type of input data to analyse (busco|ortho) [default: busco]
    -g, --orthogroups <FILE>                    If '-l busco', provide Orthofinder Orthogroups.tsv file 
    -d, --dir <DIR>                             Directory containing 
                                                    - BUSCO (v4.0.3+) TAXONNAME.*.full_table.tsv files, if '-l busco'
                                                    - BED files TAXONNAME.*.bed files (where TAXONNAME is in header of Orthogroups.tsv, if '-l ortho'
    -o, --outprefix <STR>                       Outprefix [default: test]
    -m, --missing                               Allow markers that are missing in ≥1 taxa (not implemented)
    -M, --ignore_missing                        Ignore markers that are missing in ≥1 taxa (default)
    -f, --fragmented                            Allow fragmented markers (default)
    -F, --ignore_fragmented                     Ignore fragmented markers 
    -s, --sign                                  Use sign information 
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
        self.orthogroups = self._get_path(args['--orthogroups'])
        self.label = args['--label']
        self.infiles_by_label = self._get_infiles_by_label(args['--dir'])
        self.outprefix = args['--outprefix']
        self.missing = not(args['--ignore_missing']) if args['--ignore_missing'] else (args['--missing'] != args['--ignore_missing'])
        self.fragmented = args['--fragmented'] if args['--fragmented'] else (args['--fragmented'] == args['--ignore_fragmented'])
        self.sign = args['--sign']
        print(self.__dict__)
        self._check()

    def _get_infiles_by_label(self, directory):
        infiles_by_label = collections.defaultdict(list)
        if self.label == 'busco':
            extension = '.tsv'
        elif self.label == 'ortho':
            extension = '.bed'
        else:
            sys.exit("[X] Extension is not supported: %r" % extension)
        for f in os.listdir(directory):
            if f.endswith(extension):
                infiles_by_label[self.label].append(os.path.join(directory, f))
        if len(infiles_by_label) == 0:
            sys.exit("[X] No files ending in '*.bed' in folder %s" % os.path.abspath(directory))
        return infiles_by_label
        
    def _check(self):
        if not self.infiles_by_label:
            sys.exit("[X] Please specify a directory containing input files (*.tsv).")
        if self.label == 'ortho' and not self.orthogroups:
            sys.exit("[X] Please specify an Orthofinder Orthogroups.tsv file.")

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
        graph_file = syngraph.save(parameterObj, check_consistency=True)
        print("[+] Saved Syngraph in %r" % graph_file)
        print("[*] Total runtime: %.3fs" % (timer() - main_time))
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()
