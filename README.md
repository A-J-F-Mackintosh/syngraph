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
    build               Build graph from orthology data (BUSCO *.full_table.tsv)
    recon               Reconstruct ancestral syngraph given a tree [TBI]
    viz                 Visualise graph/data
    synulate            Simulate graphs [TBI]
    
  [Options]
    -h, --help          Show this screen.
    -D, --debug         Print debug information [TBI]
    -v, --version       Show version
```
