#!/bin/bash                                                                                                                                                                                                 
ANALYSIS=$1

for REGION in met pho singleele singlemu diele dimu elemu muele
do 
    source merge.sh $ANALYSIS $REGION
done
