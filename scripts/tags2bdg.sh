#!/usr/bin/env bash
# author: oleg.shpynov@jetbrains.com

>&2 echo "tags2bdg $@"
if [ $# -lt 1 ]; then
    echo "Need 1 parameters! <TAGS>"
    exit 1
fi

TAGS=$1
cat ${TAGS} |\
    awk -v OFS='\t' 'BEGIN{C="";S=0;E=0;X=0}
    {if(C!=$1||S!=$2||E!=$3){if(X!=0){print(C,S,E,X)};C=$1;S=$2;E=$3;X=1}else{X=X+1}}
    END{if(X!=0){print(C,S,E,X)}}'
