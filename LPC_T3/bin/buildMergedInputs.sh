#!/bin/bash

: ${SUBMIT_WORKDIR:?"Need to set SUBMIT_WORKDIR"}
: ${SUBMIT_OUTDIR:?"Need to set SUBMIT_OUTDIR"}
: ${SUBMIT_LOGDIR:?"Need to set SUBMIT_LOGDIR"}

PInfo -n "$0" "Cleaning up staging areas..."
WD=$SUBMIT_WORKDIR
rm -rf $WD/*
rm -rf $SUBMIT_LOGDIR/*

doTar=0
filesetSize=20
while getopts ":tn:" opt; do
  case $opt in
    t)
      doTar=1
      ;;
    n)
      filesetSize=$OPTARG
      ;;
    :)
      PError -n "$0"  "Option -n must specify number of files"
      exit 1
      ;;
  esac 
done

PInfo -n "$0" "Acquiring configuration..."
wget -O ${WD}/list.cfg $PANDA_CFG
${CMSSW_BASE}/src/PandaAnalysis/LPC_T3/bin/configBuilder.py --infile ${WD}/list.cfg --outfile ${WD}/local.cfg --nfiles $filesetSize
cp -v ${WD}/list.cfg ${WD}/list_all.cfg 
cp -v ${WD}/local.cfg ${WD}/local_all.cfg 

cd $CMSSW_BASE/
if [[ $doTar == 1 ]]; then
  PInfo -n "$0" "Tarring up CMSSW..."
  tar --exclude-vcs -chzf cmssw.tgz src python biglib bin lib objs test external # h = --dereference symlinks
  mv -v cmssw.tgz ${WD}
fi

PInfo -n "$0" "Creating executable..."
cd ${CMSSW_BASE}/src/PandaAnalysis/LPC_T3/inputs/
cp -v ${SUBMIT_TMPL} ${WD}/skim.py
chmod 775 ${WD}/skim.py

PInfo -n "$0" "Finalizing work area..."
cat $HOME/private/$USER.txt | voms-proxy-init -voms cms --valid 168:00 -pwstdin
#voms-proxy-init -voms cms --valid 168:00 
cp -v /tmp/x509up_u$UID $WD/x509up

cp -v ${CMSSW_BASE}/src/PandaAnalysis/LPC_T3/inputs/exec.sh ${WD}

#PInfo -n "$0" "Taking a snapshot of work area..."
#eosmkdir -p $SUBMIT_OUTDIR/workdir/
#xrdcp -r ${WD} root://cmseos.fnal.gov/$SUBMIT_OUTDIR/workdir/

PInfo -n "$0" "Done!"

# input files for submission: cmssw.tgz, skim.py, x509up, local.cfg. exec.sh is the executable
