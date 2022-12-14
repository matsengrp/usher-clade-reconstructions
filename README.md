# Requirements
In order to run this analysis, run

```
git clone git@github.com:matsengrp/historydag.git
conda config --set channel_priority false
conda env create -f environment.yml
conda activate build_usher_trees
ml Singularity
```

# Running DAG Generation
Just make sure the clades of interest are in `focus_clades.txt` (one per line)
and run:

```
./driver_script.sh
```

Even once this finishes, there may be aggregation jobs still running on the
cluster. Once `full_dag.p` exists in each subdirectory of `clades`, it's done.

NOTE: The `v2` version of the driver script pushes more work onto the servers to
be run in parallel. I've found that it takes notably less time to aggregate trees
this way.

# Running Node Annotation
Run the following:
```
bash annotate_all.sh
```
This dumps a pickle of the trimmed dag and the annotated TOI to `trimmed_dag.pkl` and `annotated_toi.pk` in the appropriate
clade directory. for ecah clade 

After node annotation you can run plotting scripts from agg_mut.py with `make_annotation_plots.bash`. It
runs the following two commands:
- `python path/to/agg_mut.py plot_hists` aggregates information about all clades into a single histogram.
Examples include measuring the average ammount of certainty in each clade (`certainty_histogram.png`), and
parsimony difference between TOI and tree set (`toi_difference.png`) across all clades.
- `bash explore_annotation.bash $CLADE` generates plots about an individual clade. For example the histogram
of (uncertain) support values for each node in that clade's TOI (support_hist.png). There are several other
histograms/scatter plots that compare support to other node properties (e.g., height and clade size)

NOTE: All plots that are for a specific clade are stored in the `plots` directory of that clade (e.g., `usher-clade-reconstructions/clades/AY.36/plots/support_hist.png`). All plots that are aggregating information across all clades are stored in the `clades/plots` (e.g., `usher-clade-reconstructions/clades/plots/certainty_hisogram.png`).
