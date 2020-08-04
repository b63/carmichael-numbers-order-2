#!/bin/bash

DIR=/home/bkoirala/repos/summer-2020-research
jobdata=0

OPTIND=1
while getopts ":j" opt; do
    case $opt in
        j)
            jobdata=1
            ;;
    esac
done

echo -e "Cleaning data/intrp ..."
rm -fr data/intrp
mkdir data/intrp

echo -e "Cleaning data/outputs ..."
rm -fr data/outputs
mkdir data/outputs

echo -e "\nCleaning data/progress ..."
rm -fr data/progress
mkdir  data/progress

if [[ $jobdata -eq 1 ]]; then
    echo -e "\nCleaning data/jobdata ..."
    rm -fr data/jobdata
    mkdir data/jobdata
fi

echo -e "\nCleaning jobstreams/err ..."
rm -fr jobstreams/err
[[ -d jobstreams ]] && mkdir -p jobstreams/err

echo -e "\nCleaning jobstreams/out ..."
rm -fr jobstreams/out
[[ -d jobstreams ]] && mkdir -p jobstreams/out
