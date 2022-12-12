#!/bin/bash
set -eu
# source $HOME/miniconda3/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
conda activate build_usher_trees


CLADEDIR=$1
USHERDOCKER=docker://quay.io/matsengrp/usher
UNIQUEFASTA=$CLADEDIR/unique_seqs.fasta

echo '4) Use reference ID and vcf to reconstruct lots of trees on these samples...'
singularity exec $USHERDOCKER ./historydag/scripts/find_trees.sh -o $CLADEDIR/trees -f $UNIQUEFASTA


echo '5) write reference sequence fasta...'
REFERENCEFASTA=$CLADEDIR/reference.fasta
REFID=$(cat $CLADEDIR/trees/refid.txt)
echo ">$REFID" > $REFERENCEFASTA
python historydag/scripts/agg_mut.py lookup-in-fasta $UNIQUEFASTA $REFID >> $REFERENCEFASTA

echo '6) submit a cluster job to aggregate trees into a DAG.'
python historydag/scripts/agg_mut.py aggregate $CLADEDIR/trees/*.pb -o $CLADEDIR/full_dag.p --refseq $CLADEDIR/reference.fasta
