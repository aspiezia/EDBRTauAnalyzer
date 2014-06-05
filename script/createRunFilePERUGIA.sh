#! /bin/bash

fileToRun=$1
folder=$2
folder2=$3

echo '#!/bin/sh'
echo 'date'
echo 'hostname'
echo 'pwd'
echo 'ls'
echo "cd /data06/users/spiezia/EXO/CMSSW_5_3_13/src/Analyzer/EDBRTauAnalyzer/RISULTATI/$folder2/$folder"
echo 'source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh'
echo 'export SCRAM_ARCH=slc5_amd64_gcc462'
echo 'export VO_CMS_SW_DIR=/opt/exp_soft/cms'
echo 'source $VO_CMS_SW_DIR/cmsset_default.sh'
echo 'pwd'
echo 'cmsenv'
echo 'eval `scramv1 runtime -sh`'
echo "cmsRun $fileToRun"
