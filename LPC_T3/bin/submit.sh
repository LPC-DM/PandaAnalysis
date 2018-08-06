#!/bin/bash                                                                                                                                                                                                 

ANALYSIS=$1
REGION=$2

source ../lpc_setup.sh $ANALYSIS $REGION
sh buildMergedInputs.sh -t -n 40
python submit.py