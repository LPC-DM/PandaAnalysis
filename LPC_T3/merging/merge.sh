#!/bin/bash 

ANALYSIS=$1
REGION=$2

source ../lpc_setup.sh $ANALYSIS $REGION

if [ "$REGION" == "met" ];then
    echo MET TTbar ZtoNuNu ZJets WJets SingleTop QCD Diboson ZnunuH ZllH WmH WpH ttH ggH VBFH| xargs -n 1 -P 10 python merge.py
    while read p; do                                                                                                                                                                                     
	echo $p  | xargs -n 1 -P 10 python merge.py                                                                                                                                                    
    done <signals.txt
elif [ "$REGION" == "pho" ];then
    echo SinglePhoton GJets QCD | xargs -n 1 -P 5 python merge.py 
elif [ "$REGION" == "singleele" ];then
    echo SingleElectron TTbar ZJets WJets SingleTop QCD Diboson ZllH WmH WpH ttH | xargs -n 1 -P 10 python merge.py
elif [ "$REGION" == "singlemu" ];then
    echo MET TTbar ZJets WJets SingleTop QCD Diboson ZllH WmH WpH ttH | xargs -n 1 -P 10 python merge.py
elif [ "$REGION" == "diele" ];then
    echo SingleElectron TTbar ZJets Diboson ZllH ttH | xargs -n 1 -P 10 python merge.py
elif [ "$REGION" == "dimu" ];then
    echo MET TTbar ZJets Diboson ZllH ttH | xargs -n 1 -P 10 python merge.py
fi

