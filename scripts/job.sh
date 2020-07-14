#!/usr/bin/env bash
#PBS -l walltime=120:00:00
#PBS -N cofactor_interpolation
#
# Email address for notifications
#PBS -M bkoirala@iwu.edu
#
# Set up a job array
#PBS -t 0-9
#
# store stderr and stdout in joboutput directory
#PBS -e jobstreams/err
#PBS -o jobstreams/out



DIR=/home/hfl/tdata/rnd/edu/summer-2020-number-theory/summer-2020-research
OUTPUT_DIR=$DIR/data/interpolation
INPUT_DATA=$DIR/data/jobdata
PROGRESS_DATA=$DIR/data/progress
BIN=$DIR

PARAMS=('-m1' '-l10000')

[[ -f $OUTPUT_DIR ]] || mkdir -p $OUTPUT_DIR
[[ -f $PROGRESS_DATA ]] || mkdir -p $PROGRESS_DATA

if [[ ! -f "$INPUT_DATA/job$PBS_ARRAYID" ]]; then
    echo "no data found: $INPUT_DATA/job$PBS_ARRAYID" >&2
    exit 1
fi

pstring=$(head  -n 1 $INPUT_DATA/all)
IFS=', ' read -r -a primes <<< "$pstring"

while read line;
do
    IFS=':' read -r -a vals <<< "$line"
    dist_index=${vals[0]}
    powers=${vals[1]}

    params=("${PARAMS[@]}")
    IFS=', ' read -r -a array <<< "$powers"
    for index in "${!array[@]}";
    do
        params+=( "${primes[index]}^${array[index]}" )
    done

    stdbuf -e0 -o0 echo "$dist_index, ${powers[@]}" >> "$PROGRESS_DATA/job$PBS_ARRAYID"
    stdbuf -e0 -o0 $BIN/construct_P "${params[@]}" > "$OUTPUT_DIR/out$dist_index" 2>> "$PROGRESS_DATA/out$dist_index"

    i=$(($i+1))
done < "$INPUT_DATA/job$PBS_ARRAYID"
