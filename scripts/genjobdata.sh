#!/bin/bash

DIR=/home/hfl/tdata/rnd/edu/summer-2020-number-theory/summer-2020-research
OUTPUT_DIR=$DIR/data/jobdata
LOW=30
HIGH=31
BIN=$DIR
JOBS=""

OPTIND=1
while getopts ":o:j:C:l:h:" opt; do
    case $opt in
        o)
            OUTPUT_DIR=$OPTARG
            ;;
        j)
            JOBS=$((OPTARG))
            ;;
        C)
            BIN=$OPTARG
            ;;
        l)
            LOW=$OPTARG
            ;;
        h)
            HIGH=$OPTARG
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

if [[ -z $JOBS ]]; then
    echo "Number of jobs must be specified: -j" >&2
    exit 1
fi

[[ -f $OUTPUT_DIR ]] || mkdir -p $OUTPUT_DIR
[[ -z $BIN ]] || cd $BIN

# empty output directory
rm -rfv $OUTPUT_DIR/job* "$OUTPUT_DIR/all"

echo "Generating distributions: $LOW<=magnitude<$HIGH, primes=" "$@" "..."
$BIN/gen_distributions "-l$LOW" "-h$HIGH" "$@" > $OUTPUT_DIR/all
total=$(wc -l $OUTPUT_DIR/all | cut -d ' ' -f 1)
each=$((total/JOBS))
each=$((each<1 ? 1 : each))

echo "Total $total"
echo "Splitting among $JOBS jobs, $each each..."

jobi=0
linei=0
tail -n +2 $OUTPUT_DIR/all | while read line;
do
    if [[ $linei -eq 0 ]]; then
        echo $line  > "$OUTPUT_DIR/job$jobi"
    else
        echo $line >> "$OUTPUT_DIR/job$jobi"
    fi
    linei=$((linei+1))

    if [[ $linei -ge $each  &&  $jobi -lt $((JOBS-1)) ]]; then
        jobi=$((jobi+1))
        linei=0
    fi
done


