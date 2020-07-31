#!/bin/bash

DIR=/home/bkoirala/repos/summer-2020-research
DATADIR=$DIR/data/interpolation
JOBDATADIR=$DIR/data/jobdata

# 1st line has the prime numbers
head -n 1 $JOBDATADIR/all > $DATADIR/summary

tail -n +1 $JOBDATADIR/all | while read line; do
    [[ $line =~ ':' ]] || continue
    index=${line%:*}
    dist=${line#*:}
    density=$(tail -n 1 $DATADIR/out$index)
    density=${density#*= }

    #while read dline; do
    #    echo $dline
    #    if [[ $dline =~ ^"density"  ]]; then
    #        density=${dline#*= }
    #        break
    #    fi
    #done < "$DATADIR/out$index"

    echo "$index,$dist,$density" >> $DATADIR/summary
    printf "$index\r"
done

