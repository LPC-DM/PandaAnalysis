#!/bin/bash                                                                                                                                                                                                 
export PATH=${PATH}:${CMSSW_BASE}/src/PandaCore/bin/

#submission number
#export SUBMIT_NAME="v_8026_monoj_pho_v2"
export SUBMIT_NAME="v_8026_monoj_v2"
#export SUBMIT_NAME="testjob"
#scratch space
export scratch_area="/uscms_data/d3"

export PANDA="${CMSSW_BASE}/src/PandaAnalysis"

#export PANDA_CFG="http://t3serv001.mit.edu/~mcremone/histcatalog/test.cfg"
#export PANDA_CFG="http://t3serv001.mit.edu/~bmaier/stuff/ZpA0.txt"
#export PANDA_CFG="http://t3serv001.mit.edu/~bmaier/stuff/ZpBaryonic.txt"

#Monojet dataset
#export PANDA_CFG="http://t3serv001.mit.edu/~snarayan/histcatalog/20170522_004.cfg"
#Photon
export PANDA_CFG="http://t3serv001.mit.edu/~mcremone/eoscatalog/20170905.cfg"

export PANDA_FLATDIR="${scratch_area}/${USER}/panda/"${SUBMIT_NAME}"/flat/"
#export PANDA_FLATDIR="/uscms_data/d3/matteoc/panda/v_8026_0_4/flat/"
mkdir -p $PANDA_FLATDIR

#export SUBMIT_TMPL="skim_monojet_tmpl.py" ####
#export SUBMIT_TMPL="skim_vbf_tmpl.py"
export SUBMIT_TMPL="skim_monoj_tmpl.py"                                                                                                                                            

export SUBMIT_WORKDIR="${scratch_area}/${USER}/condor/"${SUBMIT_NAME}"/work/"
export SUBMIT_LOGDIR="${scratch_area}/${USER}/condor/"${SUBMIT_NAME}"/logs/"
export SKIM_CFGDIR="${scratch_area}/${USER}/skim/configs"

#EOS
export SUBMIT_OUTDIR="/store/user/${USER}/panda/"${SUBMIT_NAME}"/batch/"

export SKIM_MONOJET_FLATDIR="${scratch_area}/${USER}/skim/"${SUBMIT_NAME}"/monojet/"
export SKIM_MONOHIGGS_FLATDIR="${scratch_area}/${USER}/skim/"${SUBMIT_NAME}"/monohiggs_boosted/"
export SKIM_MONOHIGGS_RESOLVED_FLATDIR="${scratch_area}/${USER}/skim/"${SUBMIT_NAME}"/monohiggs_resolved/"

mkdir -p $SUBMIT_WORKDIR $SUBMIT_LOGDIR $SKIM_MONOJET_FLATDIR $SKIM_MONOHIGGS_FLATDIR $SKIM_MONOHIGGS_RESOLVED_FLATDIR
eosmkdir -p $SUBMIT_OUTDIR
#mkdir -p $SUBMIT_OUTDIR/locks/
#mkdir -p $SUBMIT_OUTDIR/workdir/

rm $SKIM_MONOJET_FLATDIR/*.sh
rm $SKIM_MONOHIGGS_FLATDIR/*.sh
rm $SKIM_MONOHIGGS_RESOLVED_FLATDIR/*.sh

ln -s $SKIM_CFGDIR/runSkim.sh  $SKIM_MONOJET_FLATDIR
ln -s $SKIM_CFGDIR/runSkimAll.sh  $SKIM_MONOJET_FLATDIR

ln -s $SKIM_CFGDIR/runSkim.sh  $SKIM_MONOHIGGS_FLATDIR
ln -s $SKIM_CFGDIR/runSkimAll.sh $SKIM_MONOHIGGS_FLATDIR

ln -s $SKIM_CFGDIR/runSkim.sh  $SKIM_MONOHIGGS_RESOLVED_FLATDIR
ln -s $SKIM_CFGDIR/runSkimAll.sh $SKIM_MONOHIGGS_RESOLVED_FLATDIR

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
echo "Checking validity of path"
echo "======================================================================="

for path in $PANDA_FLATDIR $SUBMIT_WORKDIR $SUBMIT_LOGDIR /eos/uscms${SUBMIT_OUTDIR} $SKIM_CFGDIR $SKIM_MONOJET_FLATDIR $SKIM_MONOHIGGS_FLATDIR $SKIM_MONOHIGGS_RESOLVED_FLATDIR
do
if [ -e $path ];then
echo "Path : ${path} is properly set"
else
echo "Path : ${path} is not exist, please fix it."
fi
done
echo "======================================================================"
echo ""
echo ""

#private production                                                                                                                                                       
#export PRIVATE_LOGDIR="${HOME}/cms/logs/monotop_private_panda/"
#export PRIVATE_PRODDIR="${HOME}/cms/hist/monotop_private_pandatree/"
#export PRIVATE_CFGDIR="${HOME}/cms/condor/monotop_private_panda/"

# fitting                                                                                                                                                                                                   
#export PANDA_FIT=/data/t3serv014/snarayan/CMSSW_7_4_7/
#export PANDA_XSECS=/home/snarayan/cms/cmssw/analysis/MonoTop_Xsec/
#export PANDA_XSECS=/eos/uscms/store/user/shoh/zprime_cross_section/
#export PANDA_PROD=/eos/uscms/store/user/shoh/miniaod/
#export PANDA_FITTING=${PANDA_FLATDIR}/fitting/
#mkdir -p $PANDA_FITTING/scans/ $PANDA_FITTING/logs/


