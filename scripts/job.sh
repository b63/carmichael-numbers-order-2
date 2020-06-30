#!/usr/bin/env bash
#PBS -l walltime=60:00
#PBS -N cofactor_interpolation
#
# Email address for notifications
#PBS -M bkoirala@iwu.edu
#
# Set up a job array with 100 jobs, numbered 1 to 100.
#PBS -t 0-1
 

OUTPUT_DIR=data/inter
INPUT_DATA=data/jobdata
BIN=""

PARAMS=('-m1' '-l10000')

[[ -f $OUTPUT_DIR ]] || mkdir -p $OUTPUT_DIR
[[ -z $BIN ]] || cd $BIN

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

    echo "$dist_index ${powers[@]}"
    ./construct_P "${params[@]}" > "$OUTPUT_DIR/out$dist_index"

    i=$(($i+1))
done < "$INPUT_DATA/job$PBS_ARRAYID"
