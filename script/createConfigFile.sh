#! /bin/bash

index_min=$1
index_max=$2
index_OutputFile=$3
dataset=$4
sample=$5
analysis=$6

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
if [ "$analysis" == "analyzer" ]; then
    echo "process.demo = cms.EDAnalyzer('Analyzer',"
fi
if [ "$analysis" == "fullyLeptonic" ]; then
    echo "process.demo = cms.EDAnalyzer('FullyLeptonicAnalyzer',"
fi
if [ "$analysis" == "semiLeptonic" ]; then
    echo "process.demo = cms.EDAnalyzer('SemiLeptonicAnalyzer',"
fi
if [ "$analysis" == "QCD" ]; then
    echo "process.demo = cms.EDAnalyzer('QCDAnalyzer',"
fi
if [ "$analysis" == "tauMC" ]; then
    echo "process.demo = cms.EDAnalyzer('TauMCAnalyzer',"
fi
if [ "$analysis" == "efficiency" ]; then
    echo "process.demo = cms.EDAnalyzer('Efficiency',"
fi
if [ "$analysis" == "JetIdStudy" ]; then
    echo "process.demo = cms.EDAnalyzer('JetIdStudy',"
fi
if [ "$sample" == "data" ]; then
    echo "                              isData_ = cms.untracked.bool(True),"
fi
if [ "$sample" == "MC" ]; then
    echo "                              isData_ = cms.untracked.bool(False),"
fi
echo '#                              vtxColl = cms.InputTag("primaryVertexFilter"),'
echo '                              vtxColl = cms.InputTag("goodOfflinePrimaryVertices"),'
echo '                              jetColl = cms.InputTag("selectedPatJetsCA8CHSwithQJetsForBoostedTaus"),'
echo '                              jetPrunedColl = cms.InputTag("selectedPatJetsCA8CHSprunedForBoostedTaus"),'
echo '                              electronETColl = cms.InputTag("patElectronsWithTriggerBoosted"),'
echo '                              electronColl = cms.InputTag("patElectronsWithTrigger"),'
echo '                              muonColl = cms.InputTag("patMuonsWithTriggerBoosted"),'
echo '                              tauMuTauColl = cms.InputTag("selectedPatTausMuTau"),'
echo '                              tauElTauColl = cms.InputTag("selectedPatTausEleTau"),'
echo '                              metColl = cms.InputTag("patMetShiftCorrected"),'
echo '                              uncorrmetColl = cms.InputTag("patMETs"),'
echo '                              metRawColl = cms.InputTag("patMETsRaw"),'
echo '                              ak5JetColl = cms.InputTag("patJetsWithVarCHS"),'
if [ "$dataset" == "WW" ]
then
    echo '                              NeventsTOT = cms.int32(10000431),'
    echo '                              xsec = cms.double(57.1097),'
    echo '                              lumi = cms.double(19702.),'
elif [ "$dataset" == "WZ" ]
then
    echo '                              NeventsTOT = cms.int32(10000283),'
    echo '                              xsec = cms.double(33.21),'
    echo '                              lumi = cms.double(19702.),'
elif [ "$dataset" == "ZZ" ]
then
    echo '                              NeventsTOT = cms.int32(9799908),'
    echo '                              xsec = cms.double(8.059),'
    echo '                              lumi = cms.double(19702.),'
elif [ "$dataset" == "DY100" ]
then
    echo '                              NeventsTOT = cms.int32(11764538),'
    echo '                              xsec = cms.double(32.9),'
    echo '                              lumi = cms.double(19702.),'
elif [ "$dataset" == "DY70" ]
then
    echo '                              NeventsTOT = cms.int32(12511326),'
    echo '                              xsec = cms.double(53.0),'
    echo '                              lumi = cms.double(19702.),'
elif [ "$dataset" == "DYM50_100" ]
then
    echo '                              NeventsTOT = cms.int32(1),'
    echo '                              xsec = cms.double(1),'
    echo '                              lumi = cms.double(19702.),'
elif [ "$dataset" == "DYM50_70" ]
then
    echo '                              NeventsTOT = cms.int32(1),'
    echo '                              xsec = cms.double(1),'
    echo '                              lumi = cms.double(19702.),'
elif [ "$dataset" == "QCD250" ]
then
    echo '                              NeventsTOT = cms.int32(27062078),'
    echo '                              xsec = cms.double(276000.),'
    echo '                              lumi = cms.double(19702.),'
elif [ "$dataset" == "QCD500" ]
then
    echo '                              NeventsTOT = cms.int32(30599292),'
    echo '                              xsec = cms.double(8426.),'
    echo '                              lumi = cms.double(19702.),'
elif [ "$dataset" == "QCD1000" ]
then
    echo '                              NeventsTOT = cms.int32(13843863),'
    echo '                              xsec = cms.double(204.),'
    echo '                              lumi = cms.double(19702.),'
elif [ "$dataset" == "TT" ]
then
    echo '                              NeventsTOT = cms.int32(21675970),'
    echo '                              xsec = cms.double(225.197),'
    echo '                              lumi = cms.double(19702.),'
elif [ "$dataset" == "WJets180" ]
then
    echo '                              NeventsTOT = cms.int32(9739464),'
    echo '                              xsec = cms.double(23.5),'
    echo '                              lumi = cms.double(19702.),'
elif [ "$dataset" == "WJets100" ]
then
    echo '                              NeventsTOT = cms.int32(12742382),'
    echo '                              xsec = cms.double(228.9),'
    echo '                              lumi = cms.double(19702.),'
else
    echo '                              NeventsTOT = cms.int32(1),'
    echo '                              xsec = cms.double(1.),'
    echo '                              lumi = cms.double(1.),'
fi
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
if [ "$sample" == "data" ]; then
    echo "import FWCore.PythonUtilities.LumiList as LumiList"
    echo "import FWCore.ParameterSet.Types as CfgTypes"
    echo "process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())"
    echo "JSONfile = '/data06/users/spiezia/EXO/CMSSW_5_3_13/src/Analyzer/EDBRTauAnalyzer/data/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt'"
    echo "myLumis = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')"
    echo "process.source.lumisToProcess.extend(myLumis)"
    echo ""
fi
echo "process.p = cms.Path(process.badEventFilter * process.demo)"
echo "#process.p = cms.Path(process.demo)"
echo "process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )"
