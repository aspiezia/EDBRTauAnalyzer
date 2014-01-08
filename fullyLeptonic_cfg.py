import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    #'root://osg-se.sprace.org.br//store/user/caber/BulkG_m1000_TauTau_v2/EDBR_PATtuple_edbr_zz_tautaujj_13082013/1b325ddfb984c14533be7920e22baeef/BulkG_m1000_TauTau_v2__aspiezia-BulkG_m1000_TauTau_v2-3664d28163503ca8171ba37083c39fc9__USER_3_1_zKJ.root',
    'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/aspiezia/DYJetsToLL_PtZ-100_TuneZ2star_8TeV_ext-madgraph-tarball/EDBRTauTau_DY/d697403925b5703dc02456d20f4b4874/DYJetsToLL_PtZ-100_TuneZ2star_8TeV_ext-madgraph-tarball__Summer12_DR53X-PU_S10_START53_V7C-v1__AODSIM_62_1_mO3.root',
    )
)

process.demo = cms.EDAnalyzer('FullyLeptonicAnalyzer'
)


process.badEventFilter = cms.EDFilter("HLTHighLevel",
                                      TriggerResultsTag =
                                      cms.InputTag("TriggerResults","","PAT"),
                                      HLTPaths =
                                      cms.vstring('primaryVertexFilterPath',
                                                  'noscrapingFilterPath',
                                                  'hcalLaserEventFilterPath',
                                                  'HBHENoiseFilterPath',
                                                  'trackingFailureFilterPath',
                                                  'CSCTightHaloFilterPath',
                                                  'eeBadScFilterPath',
                                                  'EcalDeadCellTriggerPrimitiveFilterPath'
                                                  ),
                                      eventSetupPathsKey = cms.string(''),
                                      andOr = cms.bool(False), # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
                                      throw = cms.bool(True)   # throw exception on unknown path names
                                      )

process.TFileService = cms.Service("TFileService",
        fileName = cms.string("output.root")
)


process.p = cms.Path(process.badEventFilter + process.demo)
#process.p = cms.Path(process.demo)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
