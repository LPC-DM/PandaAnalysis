#!/bin/bash                                                                                                                                                                                         

ANALYSIS=$1

mkdir -p $ANALYSIS

declare -A array=( [signal]=met [wmn]=singlemu [tmn]=singlemu [tme]=singlemu [wen]=singleele [ten]=singleele [tem]=singleele [zmm]=dimu [zee]=diele [pho]=pho )

BASE_COMMAND='hadd fittingForest_'
DIR=/uscms_data/d3/ahall/panda/80X-v1.5/boosted_
DIR2=/flat/limits

if [ $ANALYSIS != 'monojet' ]
then
FAIL_COMMAND="hadd ${ANALYSIS}/fittingForest_${ANALYSIS}_doublebf.root "
PASS_COMMAND="hadd ${ANALYSIS}/fittingForest_${ANALYSIS}_doublebp.root "
    for SEL in signal wmn wen tmn ten zmm zee pho
    do
	FAIL_COMMAND="${FAIL_COMMAND} ${DIR}${array[$SEL]}${DIR2}/fittingForest_${SEL}_fail.root"
	PASS_COMMAND="${PASS_COMMAND} ${DIR}${array[$SEL]}${DIR2}/fittingForest_${SEL}.root"
    done
echo "Here is the hadd commands you should run:"
echo $FAIL_COMMAND
echo $PASS_COMMAND
fi

if [ $ANALYSIS == 'monojet' ]
then
    echo "Cannot handle monojet right now."
fi

