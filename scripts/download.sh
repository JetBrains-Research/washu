#!/usr/bin/env bash
Download data from GTAC facility
# author oleg.shpynov@jetbrains.com

echo "Downloading 1901_6 and 1901_7 data"
wget -nc -r https://htcf.wustl.edu/files/G3ex3Bdw/Oltz_1901_6/
wget -nc -r https://htcf.wustl.edu/files/G3ex3Bdw/Oltz_1901_7/
find . -name *.fq.gz | xargs -i mv {} .
gunzip *.gz

# Cleanup
rm -r htcf.wustl.edu