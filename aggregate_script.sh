#!/bin/bash
set -eu
# source $HOME/miniconda3/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
conda activate build_usher_trees


CLADEDIR=$1
python historydag/scripts/agg_mut.py aggregate $CLADEDIR/trees/*.pb -o $CLADEDIR/full_dag.p --refseq $CLADEDIR/reference.fasta
