#!/bin/bash


DIR=/home/bkoirala/repos/summer-2020-research
OUTPUT_DIR=$DIR/data/jobdata
INPUT_FILE=""
BIN=$DIR/build
JOBS=""

OPTIND=1
while getopts ":o:C:i:" opt; do
    case $opt in
        o)
            OUTPUT_DIR=$OPTARG
            ;;
        C)
            BIN=$OPTARG
            ;;
        i)
            INPUT_FILE=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires argument" >&2
            exit 1
            ;;
    esac
done
shift "$((OPTIND-1))"

if [[ -z $INPUT_FILE ]]; then
    echo "Input file must be specified: -i" >&2
    exit 1
fi

[[ -f $OUTPUT_DIR ]] || mkdir -p $OUTPUT_DIR
[[ -z $BIN ]] || cd $BIN

# empty output directory
rm -rfv $OUTPUT_DIR/job* "$OUTPUT_DIR/all"

linei=0
while read line;
do
    echo "Writing \"${line:0:5}...\" to $OUTPUT_DIR/job$linei..."
    echo $line  > "$OUTPUT_DIR/job$linei"
    linei=$((linei+1))
done < $INPUT_FILE
echo "$linei total jobs"

