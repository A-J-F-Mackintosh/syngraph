# Syngraph
Toolkit for evolutionary analyses of linkage groups  

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
    infer               Models fission and fusion over a tree
    viz                 Visualise graph/data [TBI]
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

## Model fissions and fusions over a tree
```
syngraph infer -g test.pickle -t newick.txt -o test
```
