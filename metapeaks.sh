#!/bin/bash
# author Konstantin Zaytsev
# author Oleg Shpynov

which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }

# FILTERED data on chromosomes only, i.e. no contig
CHRFILES=()
PEAKS=()
for i in $@;
do
    tmpfile=${i}.chr_only.tmp
    grep -E "chr[0-9]+|chrX" $i > $tmpfile
    CHRFILES+=("$tmpfile")
    peak=`cat $tmpfile | wc -l`
    PEAKS+=("$peak")
done

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
# Extract columns 4 up to the end
awk '{for (i=4; i<=NF; i++) printf("%s%s", $i, (i==NF) ? "\n" : OFS)}' |\
# Compute all the different lines and log it
awk '{ tot[$0]++ } END { for (i in tot) print i"\t"tot[i] }'

rm ${CHRFILES[@]}
