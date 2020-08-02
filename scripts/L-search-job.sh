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
<<<<<<< HEAD




DIR=/home/bkoirala/repos/summer-2020-research
OUTPUT_DIR=$DIR/data/outputs
PROGRESS_DIR=$DIR/data/progress
INPUT_DATA=$DIR/data/jobdata
BIN=$DIR/build

PARAMS=('-l10000')

[[ -f $OUTPUT_DIR ]] || mkdir -p $OUTPUT_DIR
[[ -f $PROGRESS_DIR ]] || mkdir -p $PROGRESS_DIR

if [[ ! -f "$INPUT_DATA/job$PBS_ARRAYID" ]]; then
    echo "no data found: $INPUT_DATA/job$PBS_ARRAYID" >&2
    exit 1
fi
=======
DIR=/mnt/tdata/rnd/edu/summer-2020-number-theory/summer-2020-research
OUTPUT_DIR=$DIR/data/outputs
INPUT_DATA=$DIR/L_values
BIN=$DIR/build

PARAMS=('-l10000000')

[[ -f $OUTPUT_DIR ]] || mkdir -p $OUTPUT_DIR
>>>>>>> nonrigid

i=0
while read line;
do
    IFS=' ' read -r -a vals <<< "$line"

    params=("${PARAMS[@]}" "${vals[@]}")
<<<<<<< HEAD
    stdbuf -e0 -o0 $BIN/calc_density"${params[@]}" > "$OUTPUT_DIR/out$PBS_ARRAYID-$i" 2>> $PROGRESS_DIR/job$PBS_ARRAYID

    i=$(($i+1))
done < "$INPUT_DATA/job$PBS_ARRAYID"
=======
    echo "line $i"
    echo ${vals[@]} > $OUTPUT_DIR/out-$i
    "$BIN/calc_density" ${params[@]} >> "$OUTPUT_DIR/out-$i"

    i=$(($i+1))
done < "$INPUT_DATA"
>>>>>>> nonrigid
