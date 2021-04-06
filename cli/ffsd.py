"""

Usage: syngraph ffsd -g <FILE> -t <NWK> -i <STR> [-m <INT> -f <STR> -F <STR> -o <STR> -h]

  [Options]
    -g, --syngraph <FILE>                       Syngraph file
    -t, --tree <NWK>                            Tree in Newick format
    -i, --inference <STR>                       Inference method, either parsimony or likelihood
    -m, --minimum <INT>                         Minimum number of markers for a synteny relationship [default: 5]
    -f, --fissions <STR>                        Minimum, starting, and maximum rate of fissions per branch length unit, comma delimited [default: 0,0.01,1]
    -F, --fusions <STR>                         Minimum, starting, and maximum rate of fusions per branch length unit, comma delimited [default: 0,0.01,1]
    -o, --outprefix <STR>                       Outprefix [default: test]
    -h, --help                                  Show this message

"""

import sys
from docopt import docopt
import pathlib
import ete3
import copy
import nlopt
import random
from functools import partial
from timeit import default_timer as timer
from source import syngraph as sg

class ParameterObj():
    def __init__(self, args):
        self.syngraph = self._get_path(args['--syngraph'])
        self.tree = self._get_tree(args['--tree'])
        self.inference = self._get_inference_method(args['--inference'])
        self.outprefix = args['--outprefix']
        self.minimum = int(args['--minimum'])
        self.fission_rates = [float(rate) for rate in args["--fissions"].split(",")]
        self.fusion_rates = [float(rate) for rate in args["--fusions"].split(",")]

    def _get_path(self, infile):
        path = pathlib.Path(infile).resolve()
        if not path.exists():
            sys.exit("[X] File not found: %r" % str(infile))
        return path

    def _get_tree(self, tree_f):
        tree = ete3.Tree(str(self._get_path(tree_f)))
        for idx, node in enumerate(tree.traverse()):
            if not node.is_leaf():
                node.name = "n%s" % idx
        print(tree.get_ascii())
        return tree

    def _get_inference_method(self, inference_string):
        if not inference_string == "parsimony" and not inference_string == "likelihood":
            sys.exit("[X] Inference method not supported, please choose from parsimony/likelihood")
        return inference_string

def master_function(rates, grad, syngraph, params):
    solved_syngraph = sg.tree_traversal(rates, syngraph, params)
    evaluation = sg.solution_evaluator(rates, solved_syngraph, params)
    return evaluation

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
        random.seed(44)

        #### TO DO LIST AND QUESTIONS
        ####
        #### Should the traversals be recursive?
        #### How easy would including RTs be?
        #### Could adding in unassignables cause problems in later traversals?

        # depending on inference type, optimise the model or run once in parsimony mode 
        if parameterObj.inference == "likelihood":
            if len(parameterObj.fission_rates) == 1 and len(parameterObj.fusion_rates) == 1:
                master_function([parameterObj.fission_rates[0], parameterObj.fusion_rates[0]], None, syngraph, parameterObj)
            elif len(parameterObj.fission_rates) == 3 and len(parameterObj.fusion_rates) == 3:
                opt = nlopt.opt(nlopt.LN_NELDERMEAD, 2)
                opt.set_lower_bounds([parameterObj.fission_rates[0], parameterObj.fusion_rates[0]])
                opt.set_upper_bounds([parameterObj.fission_rates[2], parameterObj.fusion_rates[2]])
                opt.set_xtol_rel(0.01)
                specified_master_function = partial(master_function, syngraph = syngraph, params = parameterObj)
                opt.set_max_objective(specified_master_function)
                xopt = opt.optimize([parameterObj.fission_rates[1], parameterObj.fusion_rates[1]])
                print("[=] Optimised rates:\t{}".format(xopt))
                print("[=] Optimised likelihood:\t{}".format(opt.last_optimum_value()))
        elif parameterObj.inference == "parsimony":
            master_function([float(0), float(0)], None, syngraph, parameterObj)

        print("[*] Total runtime: %.3fs" % (timer() - main_time))
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)

if __name__ == '__main__':
    main()