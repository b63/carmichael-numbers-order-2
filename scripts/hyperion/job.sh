#!/usr/bin/env bash
#PBS -l walltime=120:00:00
#PBS -N cofactor_interpolation
#
# Email address for notifications
#PBS -M bkoirala@iwu.edu
#
# Set up a job array
#PBS -t 0-2
#
# store stderr and stdout in joboutput directory
#PBS -e jobstreams/err
#PBS -o jobstreams/out



DIR=/home/bkoirala/repos/summer-2020-research
OUTPUT_DIR=$DIR/data/intrp
INPUT_DATA=$DIR/data/jobdata
PROGRESS_DIR=$DIR/data/progress
BIN=$DIR/build

[[ -f $OUTPUT_DIR ]] || mkdir -p $OUTPUT_DIR

if [[ ! -f "$INPUT_DATA/job$PBS_ARRAYID" ]]; then
    echo "no data found: $INPUT_DATA/job$PBS_ARRAYID" >&2
    exit 1
fi

PARAMS=("10000000" "$INPUT_DATA/job$PBS_ARRAYID" "$INPUT_DATA/all" "$OUTPUT_DIR/out$PBS_ARRAYID" "$PROGRESS_DIR/job$PBS_ARRAYID")
$BIN/calc_density_batch "${PARAMS[@]}"
