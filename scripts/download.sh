#!/usr/bin/env bash
Download data from GTAC facility
# author oleg.shpynov@jetbrains.com

echo "Downloading 1903_1 and 1903_2 data"
wget -nc -r https://htcf.wustl.edu/files/kQePwPdD/Oltz_1903_1/
wget -nc -r https://htcf.wustl.edu/files/kQePwPdD/Oltz_1903_2/

echo "Downloading 1901_6 and 1901_7 data"
wget -nc -r https://htcf.wustl.edu/files/G3ex3Bdw/Oltz_1901_6/
wget -nc -r https://htcf.wustl.edu/files/G3ex3Bdw/Oltz_1901_7/

find . -name *.fq.gz | xargs -i cp {} .
gunzip *.gz