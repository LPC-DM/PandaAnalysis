#!/bin/bash                                                                                                                                                                                                 

ANALYSIS=$1
REGION=$2

source ../lpc_setup.sh $ANALYSIS $REGION
#sh buildMergedInputs.sh -t -n 40
sh buildMergedInputs.sh -t -n 20
eosrm ${SUBMIT_OUTDIR}/*.root 
python submit.py