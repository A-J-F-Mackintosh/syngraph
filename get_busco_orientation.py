#!/usr/bin/env python3
# coding: utf-8

"""usage: get_busco_orientation                 -d <DIR> -o <PRE> [-h]
    
    -d, --directory <DIR>                       Path of BUSCO run directory
    -o, --outprefix <PRE						Outprefix for tsv file
    -h, --help                                  Show this message
    
"""

import os
import sys
from docopt import docopt
import pandas as pd
import collections
from timeit import default_timer as timer
import numpy as np
pd.set_option('display.max_rows', None)


def get_full_table(directory):
	for file in os.listdir(directory):
		if file.startswith("run_"):
			ft_file = os.path.join(directory, file, "full_table.tsv")
			full_table = pd.read_csv(ft_file, sep='\t', usecols=[0, 1, 2, 3, 4], 
				skiprows=3, names=['name', 'status', 'seq', 'start', 'end'], 
				dtype={'name': str, 'status': str , 'seq': str, 'start': float, 'end': float})
			full_table["orientation"] = "NaN"
			return(full_table)

def get_orientation(directory, full_table):
	for file in os.listdir(directory):
		if file.startswith("run_"):
			gff_dir = os.path.join(directory, file, "augustus_output", "gff")
			for ft_row in full_table.itertuples():
				for gff_file in os.listdir(gff_dir):
					if gff_file == ft_row.name + ".gff":
						gff_path =  os.path.join(gff_dir, gff_file)
						gff_table = pd.read_csv(gff_path, sep='\t', usecols=[0, 2, 3, 4, 6], 
							names=['seq', 'feature', 'start', 'end', 'orientation'], 
							dtype={'seq': str, 'feature': str , 'start': float, 'end': float, 'orientation': str})
						gff_table = gff_table[gff_table["feature"] == "gene"]
						gff_table = gff_table.drop_duplicates()
						for g_row in gff_table.itertuples():
							if g_row.seq == ft_row.seq:
								if g_row.start == ft_row.start:
									if g_row.end == g_row.end:
										full_table.loc[ft_row.Index,'orientation'] = g_row.orientation
	full_table.to_csv(args["--outprefix"]+'.tsv', index=False, sep="\t", header=False)

args = docopt(__doc__)
full_table = get_full_table(args["--directory"])
get_orientation(args["--directory"], full_table)