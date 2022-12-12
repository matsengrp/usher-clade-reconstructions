#! /bin/bash
set -eu

# Init shell
eval "$(conda shell.bash hook)"
conda activate build_usher_trees
ml Singularity

BIGMAT=public-latest.all.masked.pb.gz
REFSEQFILE=public-latest-reference.fasta
USHERDOCKER=docker://quay.io/matsengrp/usher
mkdir -p clades
# do everything for each clade in focus_clades.txt
for CLADE in $(cat focus_clades_test.txt); do
    CLADEDIR=clades/$CLADE
    echo Reconstructing $CLADE at $CLADEDIR
    mkdir -p $CLADEDIR

    echo '1) subsetting the MAT...'
    SUBSETMAT=$CLADEDIR/subset_mat.pb
    singularity exec $USHERDOCKER matUtils extract -i $BIGMAT -o $SUBSETMAT -c $CLADE > /dev/null

    echo "2) getting list of ids in subsetted MAT at ${SUBSETMAT}..."
    CLADEIDS=$CLADEDIR/clade_ids.txt
    python historydag/scripts/agg_mut.py get-leaf-ids -t $SUBSETMAT -o $CLADEIDS > /dev/null

    echo '3) building fasta file of all unique sequences in the clade...'
    UNIQUEFASTA=$CLADEDIR/unique_seqs.fasta
    python historydag/scripts/agg_mut.py find-leaf-seq $BIGMAT $REFSEQFILE -f $CLADEIDS -o $UNIQUEFASTA -u
    
    echo '4) Use reference ID and vcf to reconstruct lots of trees on these samples...'
    # A reference ID will be chosen by find_trees as the first id + seq in the given fasta file
    sbatch -c 1 -J $CLADE -o $CLADEDIR/aggregate_2.log ./v2_find_trees.sh $CLADEDIR
done
