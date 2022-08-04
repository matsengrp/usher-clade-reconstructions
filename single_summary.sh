#! /bin/bash
set -eu

# NOTE: Might be unnecessary
# # Init shell
# conda list
eval "$(conda shell.bash hook)"
<<<<<<< HEAD
conda activate hdag-agg-new

CLADE=$1
=======
conda activate build_usher_trees

CLADE=$1

>>>>>>> dec82ab53879299b99681c59b1a798cc75a6c231
python historydag/scripts/agg_mut.py summarize $CLADE/full_dag.p -t $CLADE/trees -c $(basename $CLADE) -p -o $CLADE/stats.csv
