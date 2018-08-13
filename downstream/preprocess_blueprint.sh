#!/usr/bin/env bash
# Early beta of script, which used to convert blueprint *.bw files to format similar to our preprocessed
# RRBS *.txt files. Script is mainly for putting down steps which we've made in order to
# convert data and isn't ready utility script.


# Run in dir with blueprint WGBS *.bw files

# BigWig -> Bed
find . -name "*.bw" | xargs -n1 -P8 -I fname bash -c "echo fname | sed s/bw/wig/g | xargs -I nname bigWigToWig fname nname"
find . -name "*.wig" | xargs -n1 -P8 -I {} bash -c "wig2bed < '{}' > '{}.bed'"

# Join bs_call & bs_cov files for each donor:
find . -name "*.bs_call.*.bed" | xargs -I fname bash -c "echo fname | sed s/bs_call/bs_cov/g | xargs -I nname echo fname nname" > bs_call_cov.txt
cat bs_call_cov.txt | while read line; do prefix="${line%.bs_call*}"; suffix="${line##*bs_cov.}";  join -j 4 -o 1.1,1.2,1.3,1.4,1.5,2.5 $line > "$prefix.bs.$suffix"; done;

# Convert 0-based, + stranded data to 1-based stranded data
for f in $(find . -name "*.bs.*.bed"); do fname=${f##*/}; prefix=$(echo "${fname%.wig.bed}" | sed s/.CPG_methylation_calls.bs//g); cat $f | awk 'BEGIN {OFS="\t"; print "chrBase","chr","base","strand","coverage","freqC","freqT"} {print $1"."($2+1),$1,($2+1),"F",int($6), $5 * 100, (1-$5)*100}' > "methylcall.CpG.$prefix.mincov0.txt"; done;

# Filter cov 10:
for f in $(find . -name "*.mincov0.txt"); do prefix="${f%.mincov0.txt*}"; awk '{ if ($5 >= 10) print }' ${f} > "$prefix.mincov10.txt"; done;