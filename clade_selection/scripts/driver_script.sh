#! /bin/bash
set -eu

# Init shell
eval "$(conda shell.bash hook)"
conda activate build_usher_trees
ml Singularity

CLADEDIR=/fh/fast/matsen_e/whowards/usher-clade-reconstructions/clade_selection

sbatch -c 1 -J 2 -o $CLADEDIR/selection.log ./clade_selection/clade_selection.bash
