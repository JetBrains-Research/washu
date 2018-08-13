#!/bin/bash

WDIR='/mnt/stripe/bio/experiments/blueprint/data/'
cd $WDIR

for full_file in $(find mincov0 -maxdepth 1 -mindepth 1 -type f); 
do
    name=$(basename $full_file)
    name=${name%????};
    out_name=sorted_data/${name}.sorted.txt
    echo $out_name
    head -1 $full_file > $out_name
    tail -n+2 $full_file | sed -e 's/ \+/\t/g' | sort -k2,2 -k3,3n  >> $out_name
done
