#!/bin/bash

DIR=/home/hfl/tdata/rnd/edu/summer-2020-number-theory/summer-2020-research
DATADIR=$DIR/data/interpolation
JOBDATADIR=$DIR/data/jobdata

# 1st line has the prime numbers
head -n 1 $JOBDATADIR/all > $DATADIR/summary

tail -n +1 $JOBDATADIR/all | while read line; do
    [[ $line =~ ':' ]] || continue
    index=${line%:*}
    dist=${line#*:}
    density=""
    printf "$index"

    while read line; do
        if [[ $line =~ ^"density"  ]]; then
            density=${line#*= }
            break
        fi
    done < "$DATADIR/out$index"
    printf ", $density\n"
    echo "$index,$dist,$density" >> $DATADIR/summary
done

