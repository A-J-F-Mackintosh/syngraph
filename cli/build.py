"""

Usage: syngraph build -d <DIR> [-m|-M] [-f|-F] [-s|-S] [-o <STR> -h]

  [Options]
    -d, --dir <DIR>                             Directory containing BUSCO (v4.0.3+) TAXONNAME.*.full_table.tsv files
    -o, --outprefix <STR>                       Outprefix [default: test]
    -m, --missing                               Allow markers that are missing in ≥1 taxa (not implemented)
    -M, --ignore_missing                        Ignore markers that are missing in ≥1 taxa (default)
    -f, --fragmented                            Allow fragmented markers (default)
    -F, --ignore_fragmented                     Ignore fragmented markers 
    -s, --sign                                  Use sign information (default)
    -S, --ignore_sign                           Ignore sign information
    -h, --help                                  Show this screen.

"""

import os
import sys
from docopt import docopt
from operator import attrgetter
import pathlib
import pandas as pd
import collections
from timeit import default_timer as timer
from source import syngraph as sg

def parse_markerObjs(parameterObj):
    tmp_markerObjs = []
    status_allowed = {'Complete'}
    if parameterObj.fragmented:
        status_allowed.add('Fragmented') 
    for label, infiles in parameterObj.infiles_by_label.items(): # this might be needed later when OrthoFinder parsing is implemented...
        for infile in infiles:
            taxon_markerObjs = []
            taxon = infile.split("/")[-1].split(".")[0]
            df = pd.read_csv(infile, 
                    sep='\t', 
                    usecols=[0, 1, 2, 3, 4, 5, 7], 
                    skiprows=3,
                    names=['name', 'status', 'seq', 'start', 'end', 'sign', 'length'], 
                    dtype={'name': str, 'status': str , 'seq': str, 'start': float, 'end': float, 'sign': str, 
                    'length': float}
                    ).sort_values(['seq', 'start'], ascending=[True, True])
            for name, status, seq, start, end, sign, length in df.values.tolist():
                if status in status_allowed:
                    if parameterObj.sign:
                        if sign == "+":
                            markerObj_head = sg.MarkerObj(name=name+"_head", status=status, taxon=taxon, seq=seq, coord=start)
                            markerObj_tail = sg.MarkerObj(name=name+"_tail", status=status, taxon=taxon, seq=seq, coord=end)
                        elif sign == "-":
                            markerObj_head = sg.MarkerObj(name=name+"_head", status=status, taxon=taxon, seq=seq, coord=end)
                            markerObj_tail = sg.MarkerObj(name=name+"_tail", status=status, taxon=taxon, seq=seq, coord=start)
                        taxon_markerObjs.append(markerObj_head)
                        taxon_markerObjs.append(markerObj_tail)
                    else:
                        markerObj = sg.MarkerObj(name=name, status=status, taxon=taxon, seq=seq, coord=start)
                        taxon_markerObjs.append(markerObj)
            taxon_markerObjs = sorted(taxon_markerObjs, key=attrgetter('seq', 'coord'))
            tmp_markerObjs += taxon_markerObjs
        if parameterObj.missing == True:
            return tmp_markerObjs
        else:
            markerObjs = []
            counter = collections.Counter([tmp.name for tmp in tmp_markerObjs])
            for markerObj in tmp_markerObjs:
                if counter[markerObj.name] == len(infiles):
                    markerObjs.append(markerObj)
            return markerObjs

class ParameterObj():
    def __init__(self, args):
        self.directory = self._get_path(args['--dir'])
        self.infiles_by_label = self._get_infiles_by_label(args['--dir'])
        self.outprefix = args['--outprefix']
        self.missing = not(args['--ignore_missing']) if args['--ignore_missing'] else (args['--missing'] != args['--ignore_missing'])
        self.fragmented = args['--fragmented'] if args['--fragmented'] else (args['--fragmented'] == args['--ignore_fragmented'])
        self.sign = args['--sign'] if args['--sign'] else (args['--sign'] == args['--ignore_sign'])
        self._check()

    def _get_infiles_by_label(self, directory, data='busco'):
        infiles_by_label = collections.defaultdict(list)
        if data == 'busco':
            for f in os.listdir(directory):
                if f.endswith(".tsv"):
                    infiles_by_label['busco'].append(os.path.join(directory, f))
            if len(infiles_by_label) == 0:
                sys.exit("[X] No files ending in '*.tsv' in folder %s" % os.path.abspath(directory))
        return infiles_by_label
        
    def _check(self):
        if not self.infiles_by_label:
            sys.exit("[X] Please specify a directory containing input files (*.tsv)")

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
        print("[+] Parsing markers from %r ..." % parameterObj.directory)
        markerObjs = parse_markerObjs(parameterObj)    
        print("[+] Creating Syngraph from markers ...")
        syngraph = sg.Syngraph()
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