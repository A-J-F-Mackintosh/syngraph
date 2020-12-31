"""

Usage: syngraph ffsd -g <FILE> -t <NWK> -i <STR> [-m <INT> -r <STR> -o <STR> -h]

  [Options]
    -g, --syngraph <FILE>                       Syngraph file
    -t, --tree <NWK>                            Tree in Newick format
    -i, --inference <STR>                       Inference method of choice, either parsimony or likelihood
    -m, --minimum <INT>                         Minimum number of markers for a synteny relationship [default: 5]
    -r, --rates <STR>                           Per branch length unit rates of fission and fusion, comma delimited [default: 1,1] 
    -o, --outprefix <STR>                       Outprefix [default: test]
    -h, --help                                  Show this message

"""

import sys
from docopt import docopt
import pathlib
import ete3
import copy
import collections
from scipy import stats
from timeit import default_timer as timer
from source import syngraph as sg

class ParameterObj():
    def __init__(self, args):
        self.syngraph = self._get_path(args['--syngraph'])
        self.tree = self._get_tree(args['--tree'])
        self.inference = self._get_inference_method(args['--inference'])
        self.outprefix = args['--outprefix']
        self.minimum = int(args['--minimum'])
        self.rates = [float(rate) for rate in args["--rates"].split(",")]

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
        return tree

    def _get_inference_method(self, inference_string):
        if not inference_string == "parsimony" and not inference_string == "likelihood":
            sys.exit("[X] Inference method not supported, please choose from parsimony/likelihood")
        return inference_string

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

        #### TO DO LIST AND QUESTIONS
        ####
        #### Should the second traversal be recursive?
        #### How easy would including RTs be?
        #### Could adding in unassignables cause problems in later traversals?

        # copy syngraph
        traversal_0_syngraph = copy.deepcopy(syngraph)
        traversal_1_syngraph = copy.deepcopy(syngraph)

        # a function for getting the closest outgroup to a tree node
        # outgroup must be in 'available_taxa' and cannot be a descenant of the tree_node or the tree_nodes itself.
        def get_closest_outgroup(tree, tree_node, available_taxa):
            closest_taxon_so_far = "null"
            closest_distance_so_far = float("inf")
            for some_tree_node in tree.search_nodes():
                if some_tree_node.name in available_taxa:
                    if not some_tree_node.name in [descendant.name for descendant in tree_node.get_descendants()]:
                        if not some_tree_node.name == tree_node.name:
                            if tree_node.get_distance(some_tree_node.name) < closest_distance_so_far:
                                closest_distance_so_far = tree_node.get_distance(some_tree_node.name)
                                closest_taxon_so_far = some_tree_node.name
            return closest_taxon_so_far

        # a function which gives branch lengths by orientation, i.e. how much up and how much down
        # the median/target will be A, the child/outgroup/parent will be B
        def get_branch_distances(tree, tree_node_A, tree_node_B):
            mrca = tree.get_common_ancestor(tree_node_A, tree_node_B)
            up_distance = mrca.get_distance(tree_node_A)
            down_distance = mrca.get_distance(tree_node_B)
            return up_distance, down_distance

        # define which taxa are extant and so can be sampled from the start
        available_taxa = set()
        for leaf in parameterObj.tree.get_leaves():
            available_taxa.add(leaf.name)


        # first traversal forms triplets from a node's two children and an outgroup
        # syngraph is updated each iteration
        print("[+] Starting first traversal ...")
        print("[+] ========================================================================")
        for tree_node in parameterObj.tree.traverse(strategy='postorder'):
            if not tree_node.is_leaf() and not tree_node.is_root():
                child_1 = tree_node.get_children()[0].name
                child_2 = tree_node.get_children()[1].name
                outgroup = get_closest_outgroup(parameterObj.tree, tree_node, available_taxa)
                branch_lengths = collections.defaultdict(list)
                for taxon in child_1, child_2, outgroup:
                    up_distance, down_distance = get_branch_distances(parameterObj.tree, tree_node, taxon)
                    branch_lengths["up"].append(up_distance)
                    branch_lengths["down"].append(down_distance)
                print("[+] Inferring median genome for {} using data from {}, {}, and {} ...". format(tree_node.name, child_1, child_2, outgroup))
                print(branch_lengths["up"])
                print(branch_lengths["down"])
                traversal_0_syngraph = sg.median_genome(tree_node.name, traversal_0_syngraph, traversal_0_syngraph, [child_1, child_2, outgroup], branch_lengths, parameterObj.rates, parameterObj.minimum)
                available_taxa.add(tree_node.name)
        print("[=] ========================================================================")

        # second traversal forms triplets from the two children and the parent, as we now have parents (internal nodes) from the first traversal
        # if a node is a child of the root then the triplet is formed from the node's two children and an outgroup, i.e. the other child of the root
        # info is read from the first traversals syngraph but a new syngraph is written to
        print("[+] Starting second traversal ...")
        print("[+] ========================================================================")
        for tree_node in parameterObj.tree.traverse(strategy='postorder'):
            if not tree_node.is_leaf() and not tree_node.is_root():
                child_1 = tree_node.get_children()[0].name
                child_2 = tree_node.get_children()[1].name
                if tree_node.up.is_root():
                    outgroup = get_closest_outgroup(parameterObj.tree, tree_node, available_taxa)
                else:
                    outgroup = tree_node.up.name
                branch_lengths = collections.defaultdict(list)
                for taxon in child_1, child_2, outgroup:
                    up_distance, down_distance = get_branch_distances(parameterObj.tree, tree_node, taxon)
                    branch_lengths["up"].append(up_distance)
                    branch_lengths["down"].append(down_distance)
                print("[+] Inferring median genome for {} using data from {}, {}, and {} ...". format(tree_node.name, child_1, child_2, outgroup))
                print(branch_lengths["up"])
                print(branch_lengths["down"])                
                traversal_1_syngraph = sg.median_genome(tree_node.name, traversal_0_syngraph, traversal_1_syngraph, [child_1, child_2, outgroup], branch_lengths, parameterObj.rates, parameterObj.minimum)
        print("[=] ========================================================================")

        # evaluator should iterate from the children of the root to the leaves
        print("[=]\tBranch_start\tBranch_end\tBranch_length\tFusions\tFissions\tLikelihood")
        total_likelihood = 1
        for tree_node in parameterObj.tree.traverse(strategy='preorder'):
            if not tree_node.is_leaf() and not tree_node.is_root():
                child_1 = tree_node.get_children()[0].name
                child_2 = tree_node.get_children()[1].name
                for child in child_1, child_2:
                    likelihood = 1
                    parent_child_LMSs, unassignable_markers = sg.get_LMSs(traversal_1_syngraph, [tree_node.name, child], parameterObj.minimum)
                    parent_LMS_ios = sg.compact_synteny(traversal_1_syngraph, parent_child_LMSs, tree_node.name, "LMS_syngraph")
                    child_LMS_ios = sg.compact_synteny(traversal_1_syngraph, parent_child_LMSs, child, "LMS_syngraph")
                    fusions, fissions = sg.ffsd(sg.compact_synteny("null", [value for value in child_LMS_ios.values()], [value for value in parent_LMS_ios.values()], "index_index"))
                    likelihood *= stats.poisson.pmf(fissions, parameterObj.rates[0]*tree_node.get_distance(child))
                    likelihood *= stats.poisson.pmf(fusions, parameterObj.rates[1]*tree_node.get_distance(child))
                    total_likelihood *= likelihood
                    print("[=]\t{}\t{}\t{}\t{}\t{}\t{}".format(tree_node.name, child, tree_node.get_distance(child), fusions, fissions, likelihood))
        print(total_likelihood)
        print("[=] ========================================================================")

        print("[*] Total runtime: %.3fs" % (timer() - main_time))
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()