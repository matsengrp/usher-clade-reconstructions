#!/bin/bash
set -eu
# Init shell
eval "$(conda shell.bash hook)"
conda activate build_usher_trees

# /home/whowards/anaconda3/envs/build_usher_trees/bin/python historydag/scripts/agg_mut.py plot-hists

for CLADE in $(cat focus_clades.txt); do
    bash explore_annotation.bash $CLADE
done
