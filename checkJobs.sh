#! /bin/bash

folder=$1

cd $CMSSW_BASE/src/Analyzer/EDBRTauAnalyzer/RISULTATI/$folder

echo ""
echo ""
JOB=`grep "T---Report end\!" */* | wc -l`
TOT=`ls -altrh */r*sh | wc -l`
RUN=`qstat | grep spiezia | grep " R " | wc -l`
PND=`qstat | grep spiezia | grep " Q " | wc -l`
RES=`grep "Connection timed out" */r* | wc -l`
echo "NUMBER OF JOBS ENDED SUCCESSFULLY = $JOB/$TOT"
echo "$RUN jobs are running"
echo "$PND jobs are pending"
if [ "$RES" -gt "0" ]; then
    echo ""
    echo "YOU NEED TO RESUBMIT $RES JOBS. YOU CAN DO IT WITH:"
    echo "source resubmitJob.sh $folder"
fi
echo ""

lines1=`ls -altrh ZH1000/r*sh 2>/dev/null | wc -l`
lines2=`grep "T---Report end\!" ZH1000/* -s | wc -l`
if [ "$lines1" -ne "0" ]; then 
    echo "ZH1000:    $lines2/$lines1 root files are ready"
fi
lines1=`ls -altrh ZH1500/r*sh 2>/dev/null  | wc -l`
lines2=`grep "T---Report end\!" ZH1500/* -s | wc -l`
if [ "$lines1" -ne "0" ]; then 
    echo "ZH1500:    $lines2/$lines1 root files are ready"
fi
lines1=`ls -altrh ZH2000/r*sh 2>/dev/null  | wc -l`
lines2=`grep "T---Report end\!" ZH2000/* -s | wc -l`
if [ "$lines1" -ne "0" ]; then 
    echo "ZH2000:    $lines2/$lines1 root files are ready"
fi
lines1=`ls -altrh ZH2500/r*sh 2>/dev/null  | wc -l`
lines2=`grep "T---Report end\!" ZH2500/* -s | wc -l`
if [ "$lines1" -ne "0" ]; then 
    echo "ZH2500:    $lines2/$lines1 root files are ready"
fi
lines1=`ls -altrh RunA/r*sh 2>/dev/null 2>/dev/null | wc -l `
lines2=`grep "T---Report end\!" RunA/* -s -s | wc -l`
if [ "$lines1" -ne "0" ]; then 
    echo "RunA:      $lines2/$lines1 root files are ready"
fi
lines1=`ls -altrh RunB/r*sh 2>/dev/null  | wc -l`
lines2=`grep "T---Report end\!" RunB/* -s | wc -l`
if [ "$lines1" -ne "0" ]; then 
    echo "RunB:      $lines2/$lines1 root files are ready"
fi
lines1=`ls -altrh RunC/r*sh 2>/dev/null  | wc -l`
lines2=`grep "T---Report end\!" RunC/* -s | wc -l`
if [ "$lines1" -ne "0" ]; then 
    echo "RunC:      $lines2/$lines1 root files are ready"
fi
lines1=`ls -altrh RunD/r*sh 2>/dev/null  | wc -l`
lines2=`grep "T---Report end\!" RunD/* -s | wc -l`
if [ "$lines1" -ne "0" ]; then 
    echo "RunD:      $lines2/$lines1 root files are ready"
fi
lines1=`ls -altrh QCD250/r*sh 2>/dev/null  | wc -l`
lines2=`grep "T---Report end\!" QCD250/* -s | wc -l`
if [ "$lines1" -ne "0" ]; then 
    echo "QCD250:    $lines2/$lines1 root files are ready"
fi
lines1=`ls -altrh QCD500/r*sh 2>/dev/null  | wc -l`
lines2=`grep "T---Report end\!" QCD500/* -s | wc -l`
if [ "$lines1" -ne "0" ]; then 
    echo "QCD500:    $lines2/$lines1 root files are ready"
fi
lines1=`ls -altrh QCD1000/r*sh 2>/dev/null  | wc -l`
lines2=`grep "T---Report end\!" QCD1000/* -s | wc -l`
if [ "$lines1" -ne "0" ]; then 
    echo "QCD1000:   $lines2/$lines1 root files are ready"
fi
lines1=`ls -altrh DY70/r*sh 2>/dev/null  | wc -l`
lines2=`grep "T---Report end\!" DY70/* -s | wc -l`
if [ "$lines1" -ne "0" ]; then 
    echo "DY70:      $lines2/$lines1 root files are ready"
fi
lines1=`ls -altrh DY100/r*sh 2>/dev/null  | wc -l`
lines2=`grep "T---Report end\!" DY100/* -s | wc -l`
if [ "$lines1" -ne "0" ]; then 
    echo "DY100:     $lines2/$lines1 root files are ready"
fi
lines1=`ls -altrh DYM50_100/r*sh 2>/dev/null  | wc -l`
lines2=`grep "T---Report end\!" DYM50_100/* -s | wc -l`
if [ "$lines1" -ne "0" ]; then 
    echo "DYM50_100: $lines2/$lines1 root files are ready"
fi
lines1=`ls -altrh DYM50_70/r*sh 2>/dev/null  | wc -l`
lines2=`grep "T---Report end\!" DYM50_70/* -s | wc -l`
if [ "$lines1" -ne "0" ]; then 
    echo "DYM50_70:  $lines2/$lines1 root files are ready"
fi
lines1=`ls -altrh WW/r*sh 2>/dev/null  | wc -l`
lines2=`grep "T---Report end\!" WW/* -s | wc -l`
if [ "$lines1" -ne "0" ]; then 
    echo "WW:        $lines2/$lines1 root files are ready"
fi
lines1=`ls -altrh ZZ/r*sh 2>/dev/null  | wc -l`
lines2=`grep "T---Report end\!" ZZ/* -s | wc -l`
if [ "$lines1" -ne "0" ]; then 
    echo "ZZ:        $lines2/$lines1 root files are ready"
fi
lines1=`ls -altrh WZ/r*sh 2>/dev/null  | wc -l`
lines2=`grep "T---Report end\!" WZ/* -s | wc -l`
if [ "$lines1" -ne "0" ]; then 
    echo "WZ:        $lines2/$lines1 root files are ready"
fi
lines1=`ls -altrh TT/r*sh 2>/dev/null  | wc -l`
lines2=`grep "T---Report end\!" TT/* -s | wc -l`
if [ "$lines1" -ne "0" ]; then 
    echo "TT:        $lines2/$lines1 root files are ready"
fi
lines1=`ls -altrh WJets180/r*sh 2>/dev/null  | wc -l`
lines2=`grep "T---Report end\!" WJets180/* -s | wc -l`
if [ "$lines1" -ne "0" ]; then 
    echo "WJets180:  $lines2/$lines1 root files are ready"
fi
lines1=`ls -altrh WJetsHT/r*sh 2>/dev/null  | wc -l`
lines2=`grep "T---Report end\!" WJetsHT/* -s | wc -l`
if [ "$lines1" -ne "0" ]; then 
    echo "WJetsHT:   $lines2/$lines1 root files are ready"
fi
echo ""

cd $CMSSW_BASE/src/Analyzer/EDBRTauAnalyzer/