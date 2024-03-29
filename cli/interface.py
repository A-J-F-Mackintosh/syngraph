"""

Usage: syngraph <module> [<args>...] [-D -V -h]

  [Modules]
    build               Build graph from orthology data (e.g. BUSCO *.full_table.tsv)
    infer               Model rearrangements over a tree
    tabulate            Get table of extant and ancestral genomes
    viz                 Visualise graph/data [TBI]
    
  [Options]
    -h, --help          Show this screen.
    -D, --debug         Print debug information [TBI]
    -v, --version       Show version

  [Dependencies] 
    ------------------------------------------------------------------------------
    | $ conda install -c conda-forge networkx pandas docopt tqdm ete3 pygraphviz |
    ------------------------------------------------------------------------------

"""

import sys
from docopt import docopt
from timeit import default_timer as timer

def main():
    try:
        __version__ = '0.1.0a'
        start_time = timer()
        args = docopt(__doc__, help=True, version=__version__, options_first=True)
        run_params = {
            'module': args['<module>'],
            'debug': True if '--debug' in args['<args>'] or '-D' in args['<args>'] else False,
            'version': __version__
        }
        if args['<module>'] == 'build':
            import cli.build as build
            build.main(run_params)
        elif args['<module>'] == 'infer':
            import cli.infer as infer
            infer.main(run_params)
        elif args['<module>'] == 'tabulate':
            import cli.tabulate as tabulate
            tabulate.main(run_params)
        elif args['<module>'] == 'viz':
            import cli.viz as viz
            viz.main(run_params)
        else:
            sys.exit("%r is not a syngraph module. See 'syngraph -help'." % args['<module>'])
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %i seconds!\n" % (timer() - start_time))
        sys.exit(-1)