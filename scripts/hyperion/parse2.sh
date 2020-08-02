#!/bin/bash

DIR=/home/bkoirala/repos/summer-2020-research
DATADIR=$DIR/data/outputs
JOBDATADIR=$DIR/data/jobdata

for file in $DATADIR/* ; do
    echo "FILE: ${file##*/}"
    i=0
    while read line; do
        if [[ $line =~ ^"density" && $i == 0 ]]
        then
            i=1
            echo "$line"
        elif [[ $i -ge 1 && $i -le 2 ]]
        then
            echo "$line"
            i=$((i+1))
        elif [[ $i -gt 2 ]]
        then
            break
        fi
    done < "$file"
    echo "-----------"
done
