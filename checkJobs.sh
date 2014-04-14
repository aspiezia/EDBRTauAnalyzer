#! /bin/bash

folder=$1

cd $CMSSW_BASE/src/Analyzer/EDBRTauAnalyzer/RISULTATI/$folder

echo ""
echo ""
JOB=`grep "T---Report end\!" */* | wc -l`
TOT=`ls -altrh */r*sh | wc -l`
echo "NUMBER OF JOBS ENDED SUCCESSFULLY = $JOB/$TOT"
echo ""
 
lines1=`ls -altrh ZH1000_SL/r*sh  | wc -l`
lines2=`ls -altrh ZH1000_SL/*root | wc -l`
echo "ZH1000:    $lines2/$lines1 root files are ready"
lines1=`ls -altrh ZH1500_SL/r*sh  | wc -l`
lines2=`ls -altrh ZH1500_SL/*root | wc -l`
echo "ZH1500:    $lines2/$lines1 root files are ready"
lines1=`ls -altrh ZH2000_SL/r*sh  | wc -l`
lines2=`ls -altrh ZH2000_SL/*root | wc -l`
echo "ZH2000:    $lines2/$lines1 root files are ready"
lines1=`ls -altrh ZH2500_SL/r*sh  | wc -l`
lines2=`ls -altrh ZH2500_SL/*root | wc -l`
echo "ZH2500:    $lines2/$lines1 root files are ready"
lines1=`ls -altrh RunA_SL/r*sh  | wc -l`
lines2=`ls -altrh RunA_SL/*root | wc -l`
echo "RunA:      $lines2/$lines1 root files are ready"
lines1=`ls -altrh RunB_SL/r*sh  | wc -l`
lines2=`ls -altrh RunB_SL/*root | wc -l`
echo "RunB:      $lines2/$lines1 root files are ready"
lines1=`ls -altrh RunC_SL/r*sh  | wc -l`
lines2=`ls -altrh RunC_SL/*root | wc -l`
echo "RunC:      $lines2/$lines1 root files are ready"
lines1=`ls -altrh RunD_SL/r*sh  | wc -l`
lines2=`ls -altrh RunD_SL/*root | wc -l`
echo "RunD:      $lines2/$lines1 root files are ready"
lines1=`ls -altrh QCD250_SL/r*sh  | wc -l`
lines2=`ls -altrh QCD250_SL/*root | wc -l`
echo "QCD250:    $lines2/$lines1 root files are ready"
lines1=`ls -altrh QCD500_SL/r*sh  | wc -l`
lines2=`ls -altrh QCD500_SL/*root | wc -l`
echo "QCD500:    $lines2/$lines1 root files are ready"
lines1=`ls -altrh QCD1000_SL/r*sh  | wc -l`
lines2=`ls -altrh QCD1000_SL/*root | wc -l`
echo "QCD1000:   $lines2/$lines1 root files are ready"
lines1=`ls -altrh DY70_SL/r*sh  | wc -l`
lines2=`ls -altrh DY70_SL/*root | wc -l`
echo "DY70:      $lines2/$lines1 root files are ready"
lines1=`ls -altrh DY100_SL/r*sh  | wc -l`
lines2=`ls -altrh DY100_SL/*root | wc -l`
echo "DY100:     $lines2/$lines1 root files are ready"
lines1=`ls -altrh DYM50_100_SL/r*sh  | wc -l`
lines2=`ls -altrh DYM50_100_SL/*root | wc -l`
echo "DYM50_100: $lines2/$lines1 root files are ready"
lines1=`ls -altrh DYM50_70_SL/r*sh  | wc -l`
lines2=`ls -altrh DYM50_70_SL/*root | wc -l`
echo "DYM50_70:  $lines2/$lines1 root files are ready"
lines1=`ls -altrh WW_SL/r*sh  | wc -l`
lines2=`ls -altrh WW_SL/*root | wc -l`
echo "WW:        $lines2/$lines1 root files are ready"
lines1=`ls -altrh ZZ_SL/r*sh  | wc -l`
lines2=`ls -altrh ZZ_SL/*root | wc -l`
echo "ZZ:        $lines2/$lines1 root files are ready"
lines1=`ls -altrh WZ_SL/r*sh  | wc -l`
lines2=`ls -altrh WZ_SL/*root | wc -l`
echo "WZ:        $lines2/$lines1 root files are ready"
lines1=`ls -altrh TT_SL/r*sh  | wc -l`
lines2=`ls -altrh TT_SL/*root | wc -l`
echo "TT:        $lines2/$lines1 root files are ready"
echo ""

cd -
