#!/bin/bash                                                                                                                                                                                               

ANALYSIS=$1

source ../lpc_setup.sh $ANALYSIS met

while read p; do
    echo $p  | xargs -n 1 -P 10 python merge.py
done <signals.txt