#!/usr/bin/env bash
#PBS -l walltime=120:00:00
#PBS -N cofactor_interpolation
#
# Email address for notifications
#PBS -M bkoirala@iwu.edu
#
# Set up a job array
#PBS -t 0-47
#
# store stderr and stdout in joboutput directory
#PBS -e jobstreams/err
#PBS -o jobstreams/out



DIR=/home/bkoirala/repos/summer-2020-research
OUTPUT_DIR=$DIR/data/interpolation
INPUT_DATA=$DIR/data/jobdata
PROGRESS_DIR=$DIR/data/progress
BIN=$DIR/build

PARAMS=('-l10000000')

[[ -f $OUTPUT_DIR ]] || mkdir -p $OUTPUT_DIR
[[ -f $PROGRESS_DIR ]] || mkdir -p $PROGRESS_DIR

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

    stdbuf -e0 -o0 echo "$dist_index, ${powers[@]}" >> "$PROGRESS_DIR/job$PBS_ARRAYID"
    $BIN/calc_density "${params[@]}" > "$OUTPUT_DIR/out$dist_index" 2>/dev/null

    i=$(($i+1))
done < "$INPUT_DATA/job$PBS_ARRAYID"
