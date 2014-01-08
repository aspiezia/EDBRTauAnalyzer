#! /bin/bash

fileToRun=$1
folder=$2

echo 'MYWORKAREA=$CMSSW_BASE/src/'
echo 'cd $MYWORKAREA'
echo 'eval `scram runtime -sh`'
echo 'cd $MYWORKAREA/Analyzer/EDBRTauAnalyzer/'
echo "cd $folder"
echo ''
echo "cmsRun $fileToRun"
