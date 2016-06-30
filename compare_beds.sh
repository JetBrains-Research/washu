#!/bin/bash
# Pairwise comparison of bed tracks.
# For each pair summary length of each track and intersection is reported.
# NOTE: Tracks should be sorted, otherwise multiIntersectBed will fail.
#
# Example:
#   compare_beds.sh A.bed B.bed
# Output:
#   A.bed B.bed 1000 2000 800
#
# author oleg.shpynov@jetbrains.com

for i in $@
do
  for j in $@ 
  do 
    if [[ $i > $j ]]
    then 
      LENGTH_I=`awk '{sum+=$3-$2} END {print sum}' $i`
      LENGTH_J=`awk '{sum+=$3-$2} END {print sum}' $j`
      INTERSECTION=`multiIntersectBed -i $i $j | grep '1,2' | awk '{sum+=$3-$2} END {print sum}'`
      echo -e "$i\t$j\t$LENGTH_I\t$LENGTH_J\t$INTERSECTION"
    fi 
  done
done
