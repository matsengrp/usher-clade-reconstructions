#! /bin/bash
set -eu
eval "$(conda shell.bash hook)"
conda activate build_usher_trees
ml Singularity

BIGMAT=public-latest.all.masked.pb.gz
USHERDOCKER=docker://quay.io/matsengrp/usher
CURDIR=clade_selection
SUBSETMAT=$CURDIR/data/subset_mat.pb

singularity exec $USHERDOCKER matUtils extract -i $BIGMAT -o $SUBSETMAT -O