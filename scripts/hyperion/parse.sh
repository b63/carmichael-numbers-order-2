#!/bin/bash

DIR=/home/bkoirala/repos/summer-2020-research
DATADIR=$DIR/data/intrp
JOBDATADIR=$DIR/data/jobdata

# 1st line has the prime numbers
head -n 1 $JOBDATADIR/all > $DATADIR/summary

echo "$@" >> $DATADIR/summary
cat $DATADIR/out* >> $DATADIR/summary
