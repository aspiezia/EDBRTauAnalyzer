#! /bin/bash

index_min=$1
index_max=$2
index_OutputFile=$3
dataset=$4
index=0

echo "import FWCore.ParameterSet.Config as cms"
echo ""
echo 'process = cms.Process("Demo")'
echo ""
echo 'process.load("FWCore.MessageService.MessageLogger_cfi")'
echo "process.MessageLogger.cerr.FwkReport.reportEvery = 10000"
echo ""
echo "process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )"
echo ""
echo 'process.source = cms.Source("PoolSource",'
echo "    # replace 'myfile.root' with the source file you want to use"
echo "    fileNames = cms.untracked.vstring("


while read LINE
  do
  let index=$index+1
  if [  "$index" -lt "$index_max" ] && [ "$index" -ge "$index_min" ]
  then
  lin=" $LINE "
  echo $LINE
  fi
done < $dataset\_fileList.txt


echo "    )"
echo ")"
echo ""
echo "process.demo = cms.EDAnalyzer('FullyLeptonicAnalyzer'"
echo ")"
echo ""
echo 'process.badEventFilter = cms.EDFilter("HLTHighLevel",'
echo "                                      TriggerResultsTag ="
echo '                                      cms.InputTag("TriggerResults","","PAT"),'
echo "                                      HLTPaths ="
echo "                                      cms.vstring('primaryVertexFilterPath',"
echo "                                                  'noscrapingFilterPath',"
echo "                                                  'hcalLaserEventFilterPath',"
echo "                                                  'HBHENoiseFilterPath',"
echo "                                                  'trackingFailureFilterPath',"
echo "                                                  'CSCTightHaloFilterPath',"
echo "                                                  'eeBadScFilterPath',"
echo "                                                  'EcalDeadCellTriggerPrimitiveFilterPath'"
echo "                                                  ),"
echo "                                      eventSetupPathsKey = cms.string(''),"
echo "                                      andOr = cms.bool(False), # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true"
echo "                                      throw = cms.bool(True)   # throw exception on unknown path names"
echo "                                      )"
echo ""
echo 'process.TFileService = cms.Service("TFileService",'
echo "        fileName = cms.string('output$index_OutputFile.root')"
echo ")"
echo ""
echo "#process.p = cms.Path(process.badEventFilter + process.demo)"
echo "process.p = cms.Path(process.demo)"
echo "process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )"
