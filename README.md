# Requirements
In order to run this analysis, run

```
git clone git@github.com:matsengrp/historydag.git
conda config --set channel_priority false
conda env create -f environment.yml
conda activate build_usher_trees
conda install -c conda-forge google-api-python-client
ml Singularity
```

# Running
Just make sure the clades of interest are in `focus_clades.txt` (one per line)
and run:

```
./driver_script.sh
```

Even once this finishes, there may be aggregation jobs still running on the
cluster. Once `full_dag.p` exists in each subdirectory of `clades`, it's done.
