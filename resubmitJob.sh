#! /bin/bash

folder=$1

cd $CMSSW_BASE/src/Analyzer/EDBRTauAnalyzer/RISULTATI/$folder

echo ""
echo ""
JOB=`grep "Connection timed out" */r* | wc -l`
echo "You are going to resubmit $JOB job"
echo ""

COMMAND1=`grep "Connection timed out" */r* | sed 's/\/run/; qsub -q cms-medium -l nodes=sl5 run/g' | sed 's/\.sh/\.sh; cd -; /g' | sed -e 's@\.e.*out@@' | sed -e "s/^/cd /" | sed  's/.run.*\.sh/&; rm -rf&\.e\*/g'`
COMMAND2=`echo $COMMAND1`
#echo "$COMMAND2"
eval "$COMMAND2"
cd $CMSSW_BASE/src/Analyzer/EDBRTauAnalyzer/
echo ""