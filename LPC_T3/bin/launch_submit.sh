#!/bin/bash                                                                                                                                                                                                                                                                                                                                                                                                          
ANALYSIS=$1

for REGION in met pho singleele singlemu diele dimu
do
    source submit.sh $ANALYSIS $REGION
done