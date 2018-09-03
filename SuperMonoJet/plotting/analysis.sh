#!/bin/bash                                                                                                                                                                                                 
ANALYSIS=$1
PRESEL=$2
SELECTION=$3
FROMLIMIT=$4

source ../../LPC_T3/lpc_setup.sh $ANALYSIS $PRESEL
#echo -e "python analysis.py --region $SELECTION --analysis $ANALYSIS $FROMLIMIT"
python analysis.py --region $SELECTION --analysis $ANALYSIS $FROMLIMIT
