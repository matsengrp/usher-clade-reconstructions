#! /bin/bash
set -eu

# Init shell
eval "$(conda shell.bash hook)"
conda activate build_usher_trees
ml Singularity

BIGMAT=public-2022-03-02.all.masked.pb.gz
[ -f $BIGMAT ] || wget http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2022/03/02/public-2022-03-02.all.masked.pb.gz
# get the proper version of the historydag script 'agg_mut.py'
[ -f historydag ] || {
    wget https://github.com/matsengrp/historydag/archive/d3488a63d418d1e16abd25e41fc2da5e424bded6.zip -O historydag.zip;
    unzip historydag.zip;
    mv historydag-*/ historydag;
}

REFSEQFILE=public-latest-reference.fasta
USHERDOCKER=docker://quay.io/matsengrp/usher
mkdir -p clades
# do everything for each clade in focus_clades.txt
for CLADE in $(cat focus_clades.txt); do
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
