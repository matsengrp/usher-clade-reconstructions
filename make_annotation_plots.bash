#!/bin/bash
set -eu
# Init shell
eval "$(conda shell.bash hook)"
conda activate build_usher_trees

# Aggregated plots over all clades will be stored in `clades/plots`
python historydag/scripts/agg_mut.py plot-hists

# Individual plots for each clade will be stored in the plots directory OF THAT CLADE (e.g., `clades/19B/plots`)
for CLADE in $(cat focus_clades.txt); do
    bash explore_annotation.bash $CLADE
done
