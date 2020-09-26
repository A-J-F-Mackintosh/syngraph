# Syngraph
Analysis toolkit for evolutionary linkage group analyses of genome assemblies 

# Dependencies
Best addressed via [conda](https://docs.conda.io/en/latest/miniconda.html)

```
$ conda install -c conda-forge networkx pandas docopt tqdm ete3 pygraphviz
```

# Usage
```
Usage: syngraph <module> [<args>...] [-D -V -h]

  [Modules]
    build               Build graph from orthology data (e.g. BUSCO *.full_table.tsv)
    units               Infer ancestral chromosomal units for the graph
    ffsd                Models fission and fusion over a tree, will eventually replace units
    viz                 Visualise graph/data [TBI]
    recon               Reconstruct ancestral syngraph given a tree [TBI]
    synulate            Simulate graphs [TBI]
    
  [Options]
    -h, --help          Show this screen.
    -D, --debug         Print debug information [TBI]
    -v, --version       Show version

  [Dependencies] 
    ------------------------------------------------------------------------------
    | $ conda install -c conda-forge networkx pandas docopt tqdm ete3 pygraphviz |
    ------------------------------------------------------------------------------
```

## Build a syngraph
```
syngraph build -d directory_of_tsv_files -o test
```

## Calculate the fission-fusion distance between two genomes
```
syngraph ffsd -g test.pickle -r genus_speciesA -q genus_speciesB
```

## Find the median genome for three genomes under the fission-fusion distance
```
syngraph ffsd -g test.pickle -r genus_speciesA -q genus_speciesB -Q genus_speciesC
```
