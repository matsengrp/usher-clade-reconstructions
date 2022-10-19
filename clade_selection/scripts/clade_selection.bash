#!/bin/bash
set -eu
eval "$(conda shell.bash hook)"
conda activate build_usher_trees

~/anaconda3/envs/build_usher_trees/bin/python historydag/scripts/agg_mut.py select-clades