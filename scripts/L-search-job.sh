#!/usr/bin/env bash
#PBS -l walltime=120:00:00
#PBS -N high_density_search
#
# Email address for notifications
#PBS -M bkoirala@iwu.edu
#
# Set up a job array
#PBS -t 0-10
#
# store stderr and stdout in joboutput directory
#PBS -e jobstreams/err
#PBS -o jobstreams/out
DIR=/mnt/tdata/rnd/edu/summer-2020-number-theory/summer-2020-research
OUTPUT_DIR=$DIR/data/outputs
INPUT_DATA=$DIR/L_values
BIN=$DIR/build

PARAMS=('-l10000000')

[[ -f $OUTPUT_DIR ]] || mkdir -p $OUTPUT_DIR

i=0
while read line;
do
    IFS=' ' read -r -a vals <<< "$line"

    params=("${PARAMS[@]}" "${vals[@]}")
    echo "line $i"
    echo ${vals[@]} > $OUTPUT_DIR/out-$i
    "$BIN/calc_density" ${params[@]} >> "$OUTPUT_DIR/out-$i"

    i=$(($i+1))
done < "$INPUT_DATA"
