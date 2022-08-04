#! /bin/bash
set -eu


for CLADEDIR in clades/*/; do
    sbatch -c 1 -J $(basename $CLADEDIR) -o $CLADEDIR/stats.log ./single_summary.sh $CLADEDIR
done

while :
do
    sleep 40
    NJOBS=$(squeue -l -u $USER | wc -l)
    [ $NJOBS == 2 ] && break
done

# Now all jobs are finished!
Cladepattern="clades/*/"
CLADES=( $Cladepattern )
head -1 $(echo $CLADES)/stats.csv > all_stats.csv
for CLADEDIR in clades/*/; do
    tail -n 1 $CLADEDIR/stats.csv >> all_stats.csv || true

