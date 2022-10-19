#! /bin/bash
set -eu

# NOTE: Might be unnecessary
# # Init shell
# conda list
eval "$(conda shell.bash hook)"
conda activate build_usher_trees

CLADE=$1

python historydag/scripts/agg_mut.py summarize $CLADE/full_dag.p -t $CLADE/trees -c $(basename $CLADE) -p -o $CLADE/stats.csv
