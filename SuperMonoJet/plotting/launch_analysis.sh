#!/bin/bash                                                                                                                                                                                         

ANALYSIS=$1
FROMLIMIT=$4

mkdir -p $ANALYSIS

declare -A array=( [signal]=met [wmn]=singlemu [tmn]=singlemu [tme]=singlemu [wen]=singleele [ten]=singleele [tem]=singleele [zmm]=dimu [zee]=diele [pho]=pho )

if [ $ANALYSIS != 'monojet' ]
then
    for SEL in signal wmn wen tmn ten zmm zee pho
    do
	source analysis.sh $ANALYSIS ${array[$SEL]} $SEL $FROMLIMIT
	source analysis.sh $ANALYSIS ${array[$SEL]} ${SEL}_fail $FROMLIMIT
    done
fi

if [ $ANALYSIS == 'monojet' ]
then
    for SEL in signal wmn wen tmn ten zmm zee pho tme tem
    do
#            echo -e "source analysis.sh $ANALYSIS ${array[$SEL]} ${SEL}_$TAG $FROMLIMIT"      
	source analysis.sh $ANALYSIS ${array[$SEL]} ${SEL}_0tag $FROMLIMIT
        source analysis.sh $ANALYSIS ${array[$SEL]} ${SEL}_1tag $FROMLIMIT
        source analysis.sh $ANALYSIS ${array[$SEL]} ${SEL}_2tag $FROMLIMIT
    done
fi
