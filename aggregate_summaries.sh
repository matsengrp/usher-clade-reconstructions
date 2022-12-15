#!/bin/bash
set -eu

head -1 $(echo -e "$(ls -b clades/*/summary.csv | head -n1)") > aggregated_summary.csv

for CLADE in clades/*/; do
    tail -n +2 $CLADE/summary.csv >> aggregated_summary.csv
done
