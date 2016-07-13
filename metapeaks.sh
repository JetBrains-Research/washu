#!/bin/bash
# author Konstantin Zaytsev

# FILTERED data on chromosomes only, i.e. no contig
CHRFILES=()
PEAKS=()
for i in $@;
do
    tmpfile=${i}.chr_only.tmp
    grep -E "chr\d+|chrX" $i > $tmpfile
    CHRFILES+=("$tmpfile")
    peak=`cat $tmpfile | wc -l`
    PEAKS+=("$peak")
done

#module load bedtools2
echo "CHRFILES: ${CHRFILES[@]}"
echo "PEAKS: ${PEAKS[@]}"


keys=$(for i in $(seq 4 $(($# + 3)));
do
    echo "-k$i,$i"
done
)
range=`seq -s, 6 1 $(($# + 5))`

multiIntersectBed -i "${CHRFILES[@]}" |\
bedtools merge -c $range -o max |\
awk -v OFS="\t" '{for (i=4; i<= NF; i++) $i = int($i); print $0}' |\
sort $keys |\
bedtools groupby -g 4-$(($# + 3)) -c 1 -o count

rm ${CHRFILES[@]}
