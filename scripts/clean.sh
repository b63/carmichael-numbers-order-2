#!/bin/bash

DIR=/home/hfl/tdata/rnd/edu/summer-2020-number-theory/summer-2020-research
jobdata=0

OPTIND=1
while getopts ":j" opt; do
    case $opt in
        j)
            jobdata=1
            ;;
    esac
done

echo -e "Cleaning data/interpolation ..."
rm -fv data/interpolation/*
echo -e "Cleaning data/outputs ..."
rm -fv data/outputs/*
echo -e "\nCleaning data/progress ..."
rm -fv data/progress/*

if [[ $jobdata -eq 1 ]]; then
    echo -e "\nCleaning data/jobdata ..."
    rm -fv data/jobdata/*
fi

echo -e "\nCleaning jobstreams/err ..."
rm -fv jobstreams/err/*

echo -e "\nCleaning jobstreams/out ..."
rm -fv jobstreams/out/*
