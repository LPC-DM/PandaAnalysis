#!/bin/bash                                                                                                                                                                                         

ANALYSIS=$1

declare -A array=( [signal]=met [wmn]=singlemu [tmn]=singlemu [tme]=singlemu [wen]=singleele [ten]=singleele [tem]=singleele [zmm]=dimu [zee]=diele [pho]=pho )

if [ $ANALYSIS != 'monojet' ]
then
    for SEL in signal wmn wen tmn ten zmm zee pho
    do
	source makeFittingForest.sh $ANALYSIS ${array[$SEL]} $SEL 
	source makeFittingForest.sh $ANALYSIS ${array[$SEL]} ${SEL}_fail 
    done
fi

if [ $ANALYSIS == 'monojet' ]
then
    for SEL in signal wmn wen tmn ten zmm zee pho tme tem
    do
#            echo -e "source makeFittingForest.sh $ANALYSIS ${array[$SEL]} ${SEL}_$TAG $FROMLIMIT"      
	source makeFittingForest.sh $ANALYSIS ${array[$SEL]} ${SEL}_0tag 
        source makeFittingForest.sh $ANALYSIS ${array[$SEL]} ${SEL}_1tag 
        source makeFittingForest.sh $ANALYSIS ${array[$SEL]} ${SEL}_2tag 
    done
fi
