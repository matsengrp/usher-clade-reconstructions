#!/bin/bash
set -eu
# Init shell
eval "$(conda shell.bash hook)"
conda activate build_usher_trees

# do everything for each clade in focus_clades.txt
for CLADE in $(cat focus_clades.txt); do
    CLADEDIR=clades/$CLADE

    echo 'submit a cluster job to annotate TOI nodes of ' $CLADE
    sbatch -c 1 -J $CLADE -o $CLADEDIR/annotate.log ./annotate_tree.sh $CLADE
done