#!/bin/bash
set -eu
# Init shell
eval "$(conda shell.bash hook)"
conda activate build_usher_trees

CLADE=$1
CLADEDIR=clades/$CLADE
SUBSETMAT=$CLADEDIR/subset_mat.pb
REFERENCEFASTA=$CLADEDIR/reference.fasta
UNIQUEFASTA=$CLADEDIR/unique_seqs.fasta

python historydag/scripts/agg_mut.py explore-annotation $SUBSETMAT $REFERENCEFASTA $CLADEDIR $UNIQUEFASTA
