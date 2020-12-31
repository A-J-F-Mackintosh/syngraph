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
    ffsd                Models fission and fusion over a tree
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
syngraph ffsd -g test.pickle -t newick.txt -i likelihood -r 0.1,0.1
```
