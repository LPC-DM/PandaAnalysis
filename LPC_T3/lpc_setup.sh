#!/bin/bash

ANALYSIS=$1
REGION=$2
if [ -z "$ANALYSIS" ];then
    echo -e "\033[0;31m Environment is not correctly setup \033[0m"
    echo -e "Please specify Analysis: \033[0;33m boosted \033[0m ; \033[1;36m resolved \033[0m ; \033[1;35m monojet \033[0m"
    exit 0
fi

if [ "$ANALYSIS" == "boosted" ];then
    COLOR="\033[0;33m"
elif [ "$ANALYSIS" == "resolved" ];then
    COLOR="\033[1;36m"
elif [ "$ANALYSIS" == "monojet" ];then
    COLOR="\033[1;35m"
fi

if [ -z "$REGION" ];then
    echo -e "\033[0;31m Environment is not correctly setup \033[0m"
    echo -e "Please specify region: \033[0;33m met \033[0m ; \033[1;36m singleele \033[0m ; \033[1;35m singlemu \033[0m ; \033[1;35m diele \033[0m ; \033[1;35m dimu \033[0m ; \033[1;35m pho \033[0m"
    exit 0
fi
                                                                                                                            
export PATH=${PATH}:${CMSSW_BASE}/src/PandaCore/bin/

#submission number
export SUBMIT_NAME="80X-v1dot1"
#scratch space
export scratch_area="/uscms_data/d3"
export PANDA="${CMSSW_BASE}/src/PandaAnalysis"
#cfg file
export PANDA_CFG="http://sundleeb.web.cern.ch/sundleeb/panda_config/20180801_${REGION}.cfg"

#skim
export SUBMIT_TMPL="skim_${ANALYSIS}_${REGION}_tmpl.py"
#panda's 
export PANDA_FLATDIR="/uscms_data/d3/${USER}/panda/${SUBMIT_NAME}/${ANALYSIS}_${REGION}/flat/"
export SUBMIT_OUTDIR="/store/user/lpcmetx/panda/${SUBMIT_NAME}/${ANALYSIS}_${REGION}/batch/" 

#condor's
export SUBMIT_WORKDIR="${scratch_area}/lpcmetx/condor/${SUBMIT_NAME}/${ANALYSIS}_${REGION}_1/work/"
export SUBMIT_LOGDIR="${scratch_area}/lpcmetx/condor/${SUBMIT_NAME}/${ANALYSIS}_${REGION}_1/logs/"
mkdir -p $PANDA_FLATDIR $SUBMIT_WORKDIR $SUBMIT_LOGDIR
eosmkdir -p $SUBMIT_OUTDIR

export SUBMIT_CONFIG=T2  # allow running on T3 or T2. if $SUBMIT_CONFIG==T3, then only run on T3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
echo ""
echo ""
cat << "EOF"
  _____        _   _ _____                      ______ _   _          _          _____  __  __ 
 |  __ \ /\   | \ | |  __ \   /\        ____   |  ____| \ | |   /\   | |        |  __ \|  \/  |
 | |__) /  \  |  \| | |  | | /  \      / __ \  | |__  |  \| |  /  \  | |  ______| |  | | \  / |
 |  ___/ /\ \ | . ` | |  | |/ /\ \    / / _` | |  __| | . ` | / /\ \ | | |______| |  | | |\/| |
 | |  / ____ \| |\  | |__| / ____ \  | | (_| | | |    | |\  |/ ____ \| |____    | |__| | |  | |
 |_| /_/    \_\_| \_|_____/_/    \_\  \ \__,_| |_|    |_| \_/_/    \_\______|   |_____/|_|  |_|
                                       \____/                                                  
EOF
echo ""
echo "Checking ENV path"
echo "======================================================================="

for path in $PANDA_FLATDIR /eos/uscms${SUBMIT_OUTDIR} $SUBMIT_WORKDIR $SUBMIT_LOGDIR
do
if [ -e $path ];then
echo -e "Path : \033[0;32m ${path} is properly set \033[0m"
else
echo -e "Path : \033[0;31m ${path} does not exist, please fix it. \033[0m"
fi
done
echo "======================================================================"
echo "INFO"
echo "======================================================================"
echo -e "Analysis     = ${COLOR}${ANALYSIS}\033[0m"
echo "Submit Name  = ${SUBMIT_NAME}"
echo "cfg selected = ${PANDA_CFG}"
echo -e "submit tmpl  =${COLOR} ${SUBMIT_TMPL}\033[0m"
echo "======================================================================"
echo ""
echo ""
