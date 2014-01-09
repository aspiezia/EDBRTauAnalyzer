// -*- C++ -*-
//
// Package:    SemiLeptonicAnalyzer
// Class:      SemiLeptonicAnalyzer
// 
/**\class SemiLeptonicAnalyzer SemiLeptonicAnalyzer.cc ExoDiBosonResonances/EDBRTauAnalyzer/src/SemiLeptonicAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Aniello Spiezia,21 1-007,+41227676459,
//         Created:  Mon Sep  9 13:14:05 CEST 2013
// $Id$
//
//


// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//new inclusion
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h" 
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "TauAnalysis/CandidateTools/interface/NSVfitStandaloneAlgorithm.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TMath.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/VectorUtil.h"
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"


//
// class declaration
//

class SemiLeptonicAnalyzer : public edm::EDAnalyzer {
public:
  explicit SemiLeptonicAnalyzer(const edm::ParameterSet&);
  ~SemiLeptonicAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  float MuonPFIso(pat::MuonCollection::const_iterator muon, bool highpt);
  void BtagVeto(edm::Handle<pat::JetCollection> ak5jetCands, int & nbtagsL, int & nbtagsM, int & nbtagsT, pat::JetCollection::const_iterator SelectedJet);

  TH1D* Nevents;

  TTree *TreeMuoTau;                        TTree *TreeEleTau;
  float m_jetPtMuoTau;			    float m_jetPtEleTau;
  float m_jetEtaMuoTau;			    float m_jetEtaEleTau;
  float m_jetMassMuoTau;		    float m_jetMassEleTau;
  float m_jetSubjettinessMuoTau;	    float m_jetSubjettinessEleTau;
  float m_dPhiJetMetMuoTau;		    float m_dPhiJetMetEleTau;
  float m_dRJetMetMuoTau;		    float m_dRJetMetEleTau;
  float m_dRJetLepMuoTau;		    float m_dRJetLepEleTau;
  float m_dRJetTauMuoTau;		    float m_dRJetTauEleTau;
  float m_dPhiTauMetMuoTau;		    float m_dPhiTauMetEleTau;
  float m_dRTauMetMuoTau;		    float m_dRTauMetEleTau;
  float m_dRTauLepMuoTau;		    float m_dRTauLepEleTau;
  float m_dPhiLepMetMuoTau;		    float m_dPhiLepMetEleTau;
  float m_dRLepMetMuoTau;		    float m_dRLepMetEleTau;
  float m_dRZZMuoTau;			    float m_dRZZEleTau;
  float m_tauPtMuoTau;			    float m_tauPtEleTau;
  float m_tauEtaMuoTau;			    float m_tauEtaEleTau;
  float m_lepPtMuoTau;			    float m_lepPtEleTau;
  float m_lepEtaMuoTau;			    float m_lepEtaEleTau;
  float m_lepPFIsoMuoTau;		    float m_lepPFIsoEleTau;
  float m_metMuoTau;			    float m_metEleTau;
  float m_MassSvfitTauLepMuoTau;	    float m_MassSvfitTauLepEleTau;
  float m_XMassSVFitMuoTau;		    float m_XMassSVFitEleTau;
  int m_nbtagsLMuoTau;			    int m_nbtagsLEleTau;
  int m_nbtagsMMuoTau;			    int m_nbtagsMEleTau;
  int m_nbtagsTMuoTau;			    int m_nbtagsTEleTau;
  int m_trigger650MuoTau;		    int m_trigger650EleTau;
  int m_trigger320MuoTau;		    int m_trigger320EleTau;
  int m_NVerticesMuoTau;		    int m_NVerticesEleTau;
  float m_PUWeightMuoTau;		    float m_PUWeightEleTau;
  float m_metPxMuoTau;			    float m_metPxEleTau;
  float m_metPyMuoTau;			    float m_metPyEleTau;
  float m_metEtaMuoTau;			    float m_metEtaEleTau;
  float m_metPhiMuoTau;         	    float m_metPhiEleTau;

  edm::LumiReWeighting LumiWeights_;
  bool isData; 
  bool sideband; 
  edm::InputTag vtxColl_;
  edm::InputTag jetColl_;
  edm::InputTag jetPrunedColl_;
  edm::InputTag electronColl_;
  edm::InputTag muonColl_;
  edm::InputTag tauMuTauColl_;
  edm::InputTag tauElTauColl_;
  edm::InputTag metColl_;
  edm::InputTag metRawColl_;
  edm::InputTag ak5JetColl_;

  // ----------member data ---------------------------
};

using namespace edm;
using namespace std;


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SemiLeptonicAnalyzer::SemiLeptonicAnalyzer(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
  isData = iConfig.getUntrackedParameter<bool>("isData_");
  sideband = iConfig.getUntrackedParameter<bool>("sideband_");
  vtxColl_ = iConfig.getParameter<edm::InputTag>("vtxColl"); 
  jetColl_ = iConfig.getParameter<edm::InputTag>("jetColl"); 
  jetPrunedColl_ = iConfig.getParameter<edm::InputTag>("jetPrunedColl"); 
  electronColl_ = iConfig.getParameter<edm::InputTag>("electronColl"); 
  muonColl_ = iConfig.getParameter<edm::InputTag>("muonColl"); 
  tauMuTauColl_ = iConfig.getParameter<edm::InputTag>("tauMuTauColl");
  tauElTauColl_ = iConfig.getParameter<edm::InputTag>("tauElTauColl"); 
  metColl_ = iConfig.getParameter<edm::InputTag>("metColl"); 
  metRawColl_ = iConfig.getParameter<edm::InputTag>("metRawColl"); 
  ak5JetColl_ = iConfig.getParameter<edm::InputTag>("ak5JetColl"); 


  // True number of interaction for data produced as in: https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
  TFile *da_=new TFile ("/shome/aspiezia/EXO/CMSSW_5_3_13/src/Analyzer/EDBRTauAnalyzer/data/MyDataPileupHistogram_True.root");
  TH1F *da = (TH1F*) da_->Get("pileup");
  
  // MC distribution of true number of interactions as in: https://twiki.cern.ch/twiki/bin/view/CMS/Pileup_MC_Gen_Scenarios
  Double_t dat[60] = {2.560E-06, 5.239E-06, 1.420E-05, 5.005E-05, 1.001E-04, 2.705E-04, 1.999E-03, 6.097E-03, 1.046E-02, 1.383E-02, 1.685E-02, 2.055E-02, 2.572E-02, 3.262E-02, 4.121E-02, 4.977E-02, 5.539E-02, 5.725E-02, 5.607E-02, 5.312E-02, 5.008E-02, 4.763E-02, 4.558E-02, 4.363E-02, 4.159E-02, 3.933E-02, 3.681E-02, 3.406E-02, 3.116E-02, 2.818E-02, 2.519E-02, 2.226E-02, 1.946E-02, 1.682E-02, 1.437E-02, 1.215E-02, 1.016E-02, 8.400E-03, 6.873E-03, 5.564E-03, 4.457E-03, 3.533E-03, 2.772E-03, 2.154E-03, 1.656E-03, 1.261E-03, 9.513E-04, 7.107E-04, 5.259E-04, 3.856E-04, 2.801E-04, 2.017E-04, 1.439E-04, 1.017E-04, 7.126E-05, 4.948E-05, 3.405E-05, 2.322E-05, 1.570E-05, 5.005E-06 };
  
  //PileUp weights calculation
  double d,m;
  std::vector< float > mcNum; 
  std::vector< float > dataNum;
  for (Int_t i=1; i< 50; i++){
    m=dat[i-1];
    d=da->GetBinContent(i);
    mcNum.push_back(m);
    dataNum.push_back(d); 
  }
  LumiWeights_=edm::LumiReWeighting(mcNum, dataNum);
}


SemiLeptonicAnalyzer::~SemiLeptonicAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SemiLeptonicAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 
  Nevents->Fill(1);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(vtxColl_, vertices);
  reco::Vertex primaryVertex;
  primaryVertex = vertices->at(0);

  edm::Handle<pat::JetCollection> CA8JetswithQjets;
  iEvent.getByLabel(jetColl_, CA8JetswithQjets);
  edm::Handle<pat::JetCollection> CA8JetsPruned;
  iEvent.getByLabel(jetPrunedColl_, CA8JetsPruned);

  edm::Handle<pat::MuonCollection> muoH;
  iEvent.getByLabel(muonColl_, muoH);

  edm::Handle<pat::ElectronCollection> eleH;
  iEvent.getByLabel(electronColl_, eleH);

  Handle<pat::TauCollection> tauMuTauHandle;
  iEvent.getByLabel(tauMuTauColl_,tauMuTauHandle);

  Handle<pat::TauCollection> tauElTauHandle;
  iEvent.getByLabel(tauElTauColl_,tauElTauHandle);

  edm::Handle<pat::METCollection> met;
  iEvent.getByLabel(metColl_, met);

  edm::Handle<pat::METCollection> metRaw;
  iEvent.getByLabel(metRawColl_, metRaw);

  edm::Handle<pat::JetCollection> ak5jetCands;
  iEvent.getByLabel(ak5JetColl_,ak5jetCands);

  //TRIGGER
  bool isFired_HLT_HT650 = false;
  bool isFired_HLT_PFJet320 = false;
  edm::Handle<edm::TriggerResults> trigResults;
  edm::InputTag trigResultsTag("TriggerResults","","HLT");  
  iEvent.getByLabel(trigResultsTag,trigResults);
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*trigResults);
  unsigned int TrggIndex_PFHT650_v5( trigNames.triggerIndex("HLT_PFHT650_v5"));
  unsigned int TrggIndex_PFHT650_v6( trigNames.triggerIndex("HLT_PFHT650_v6"));
  unsigned int TrggIndex_PFHT650_v7( trigNames.triggerIndex("HLT_PFHT650_v7"));
  unsigned int TrggIndex_PFHT650_v8( trigNames.triggerIndex("HLT_PFHT650_v8"));
  unsigned int TrggIndex_PFHT650_v9( trigNames.triggerIndex("HLT_PFHT650_v9"));
  unsigned int TrggIndex_PFNoPUHT650_v1( trigNames.triggerIndex("HLT_PFNoPUHT650_v1"));
  unsigned int TrggIndex_PFNoPUHT650_v3( trigNames.triggerIndex("HLT_PFNoPUHT650_v3"));
  unsigned int TrggIndex_PFNoPUHT650_v4( trigNames.triggerIndex("HLT_PFNoPUHT650_v4"));
  unsigned int TrggIndex_PFJet320_v3( trigNames.triggerIndex("HLT_PFJet320_v3"));
  unsigned int TrggIndex_PFJet320_v4( trigNames.triggerIndex("HLT_PFJet320_v4"));
  unsigned int TrggIndex_PFJet320_v5( trigNames.triggerIndex("HLT_PFJet320_v5"));
  unsigned int TrggIndex_PFJet320_v6( trigNames.triggerIndex("HLT_PFJet320_v6"));
  unsigned int TrggIndex_PFJet320_v8( trigNames.triggerIndex("HLT_PFJet320_v8"));
  unsigned int TrggIndex_PFJet320_v9( trigNames.triggerIndex("HLT_PFJet320_v9"));
  if(TrggIndex_PFHT650_v5 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFHT650_v5);
  if(TrggIndex_PFHT650_v6 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFHT650_v6);
  if(TrggIndex_PFHT650_v7 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFHT650_v7);
  if(TrggIndex_PFHT650_v8 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFHT650_v8);
  if(TrggIndex_PFHT650_v9 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFHT650_v9);
  if(TrggIndex_PFNoPUHT650_v1 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFNoPUHT650_v1);
  if(TrggIndex_PFNoPUHT650_v3 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFNoPUHT650_v3);
  if(TrggIndex_PFNoPUHT650_v4 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFNoPUHT650_v4);
  if(TrggIndex_PFJet320_v3 < trigResults->size()) isFired_HLT_PFJet320 = trigResults->accept(TrggIndex_PFJet320_v3);
  if(TrggIndex_PFJet320_v4 < trigResults->size()) isFired_HLT_PFJet320 = trigResults->accept(TrggIndex_PFJet320_v4);
  if(TrggIndex_PFJet320_v5 < trigResults->size()) isFired_HLT_PFJet320 = trigResults->accept(TrggIndex_PFJet320_v5);
  if(TrggIndex_PFJet320_v6 < trigResults->size()) isFired_HLT_PFJet320 = trigResults->accept(TrggIndex_PFJet320_v6);
  if(TrggIndex_PFJet320_v8 < trigResults->size()) isFired_HLT_PFJet320 = trigResults->accept(TrggIndex_PFJet320_v8);
  if(TrggIndex_PFJet320_v9 < trigResults->size()) isFired_HLT_PFJet320 = trigResults->accept(TrggIndex_PFJet320_v9); 

  //PILEUP WEIGHT
  double MyWeight = 1;
  if(!isData){
    //edm::EventBase* iEventB = dynamic_cast<edm::EventBase*>(&iEvent);
    Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    float Tnpv = -1;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      int BX = PVI->getBunchCrossing();
      if(BX == 0) { 
	Tnpv = PVI->getTrueNumInteractions();
	continue;
      }
    }
    MyWeight = LumiWeights_.weight( Tnpv );
  }

  
  //JET SELECTION
  pat::JetCollection::const_iterator SelectedJet;
  float massZ=-9999; float ptZ=-999; bool foundJet=false; float tau21Z=-9999;
  for(pat::JetCollection::const_iterator jet = CA8JetswithQjets->begin(); jet != CA8JetswithQjets->end(); ++jet) {
    float dRmin = 9999.; float mass = 0.;
    for(pat::JetCollection::const_iterator jetPruned = CA8JetsPruned->begin(); jetPruned != CA8JetsPruned->end(); ++jetPruned) {
      float dRtmp = ROOT::Math::VectorUtil::DeltaR(jet->p4(),jetPruned->p4());
      if(dRtmp<dRmin && dRtmp<0.8){//matching failed if greater than jet radius
	dRmin=dRtmp;
	mass=jetPruned->mass();
      }
    }
    if(jet->muonEnergyFraction()>=0.99) continue;
    if(jet->photonEnergyFraction()>=0.99) continue;
    if(jet->chargedEmEnergyFraction()>=0.99) continue;
    if(jet->neutralHadronEnergyFraction()>=0.99) continue;
    if(jet->chargedHadronEnergyFraction()<=0.00) continue;
    if(jet->pt()<400) continue;
    if(abs(jet->eta())>2.4) continue;
    if(sideband) {if(!(mass>20 && mass<70))  continue;}
    else         {if(!(mass>70 && mass<110)) continue;}
    if(jet->userFloat("tau2")/jet->userFloat("tau1")>0.75) continue;
    foundJet=true;
    if(jet->pt()>ptZ){
      massZ=mass;
      ptZ=jet->pt();
      tau21Z=jet->userFloat("tau2")/jet->userFloat("tau1")>0.75;
      SelectedJet=jet;
    }
  }
  
  //TAU SELECTION - MUTAU
  float ptTauMuTau=-99; bool foundTauMuTau=false;
  pat::TauCollection::const_iterator SelectedTauMuTau;
  for (pat::TauCollection::const_iterator patTau = tauMuTauHandle->begin(); patTau != tauMuTauHandle->end(); ++patTau ) {
    if(patTau->pt()<20) continue;
    if(abs(patTau->eta())>2.4) continue;
    if(patTau->tauID("decayModeFinding")<0.5) continue;
    if(patTau->tauID("againstMuonLoose")<0.5) continue;
    if(patTau->tauID("againstElectronLoose")<0.5) continue;
    if(patTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr")<0.5) continue;
    foundTauMuTau=true;
    if(patTau->pt()>ptTauMuTau){
      SelectedTauMuTau=patTau;
      ptTauMuTau=patTau->pt();
    }
  }
  
  //TAU SELECTION - ELTAU
  float ptTauElTau=-99; bool foundTauElTau=false;
  pat::TauCollection::const_iterator SelectedTauElTau;
  for (pat::TauCollection::const_iterator patTau = tauElTauHandle->begin(); patTau != tauElTauHandle->end(); ++patTau ) {
    if(patTau->pt()<20) continue;
    if(abs(patTau->eta())>2.4) continue;
    if(patTau->tauID("decayModeFinding")<0.5) continue;
    if(patTau->tauID("againstMuonLoose")<0.5) continue;
    if(patTau->tauID("againstElectronLoose")<0.5) continue;
    if(patTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr")<0.5) continue;
    foundTauElTau=true;
    if(patTau->pt()>ptTauElTau){
      SelectedTauElTau=patTau;
      ptTauElTau=patTau->pt();
    }
  }

  //MUON SELECTION
  float ptMuon=-99; bool foundMuon=false;
  pat::MuonCollection::const_iterator SelectedMuon;
  for(pat::MuonCollection::const_iterator muon = muoH->begin(); muon != muoH->end(); ++muon) {
    if(!(muon->isGlobalMuon())) continue;
    reco::TrackRef cktTrack = (muon::tevOptimized(*muon, 200, 17., 40., 0.25)).first;
    if(cktTrack->pt()<10) continue;
    if(abs(cktTrack->eta())>2.4) continue;
    if(abs(cktTrack->phi())>3.2) continue;
    if((cktTrack->ptError()/cktTrack->pt())>0.3) continue;
    if(muon->globalTrack()->hitPattern().numberOfValidMuonHits()<=0) continue;
    if(muon->numberOfMatches()<=1) continue;
    if(fabs(cktTrack->dxy(primaryVertex.position()))>=0.2) continue;
    if(fabs(cktTrack->dz( primaryVertex.position()))>=0.5) continue;
    if(muon->innerTrack()->hitPattern().numberOfValidPixelHits()<=0) continue;
    if(muon->innerTrack()->hitPattern().trackerLayersWithMeasurement()<=5) continue;
    foundMuon=true;
    if(cktTrack->pt()>ptMuon){
      SelectedMuon=muon;
      ptMuon=cktTrack->pt();
    }
  }

  //ELECTRON SELECTION
  float ptElectron=-99; bool foundElectron=false;
  pat::ElectronCollection::const_iterator SelectedElectron;
  for(pat::ElectronCollection::const_iterator electron = eleH->begin(); electron != eleH->end(); ++electron) {
    if(electron->pt()<10) continue;
    if(fabs(electron->superCluster()->eta())<=1.479){
      if(electron->deltaEtaSuperClusterTrackAtVtx()>=0.004) continue;
      if(electron->deltaPhiSuperClusterTrackAtVtx()>=0.030) continue;
      if(electron->sigmaIetaIeta()>=0.01) continue;
      if(electron->hadronicOverEm()>=0.12) continue;
      if(fabs(electron->gsfTrack()->dxy(primaryVertex.position()))>=0.02) continue;
      if(fabs(electron->gsfTrack()->dz(primaryVertex.position()))>=0.1) continue;
      if((fabs(1/electron->ecalEnergy() - electron->eSuperClusterOverP()/electron->ecalEnergy()))>=0.05) continue;
      if(electron->passConversionVeto()==0) continue;
      if(electron->gsfTrack()->trackerExpectedHitsInner().numberOfHits()!=0) continue;
    }
    if(fabs(electron->superCluster()->eta())>1.479 && fabs(electron->superCluster()->eta())<2.5){
      if(electron->deltaEtaSuperClusterTrackAtVtx()>=0.005) continue;
      if(electron->deltaPhiSuperClusterTrackAtVtx()>=0.020) continue;
      if(electron->sigmaIetaIeta()>=0.03) continue;
      if(electron->hadronicOverEm()>=0.10) continue;
      if(fabs(electron->gsfTrack()->dxy(primaryVertex.position()))>=0.02) continue;
      if(fabs(electron->gsfTrack()->dz(primaryVertex.position()))>=0.1) continue;
      if((fabs(1/electron->ecalEnergy() - electron->eSuperClusterOverP()/electron->ecalEnergy()))>=0.05) continue;
      if(electron->passConversionVeto()==0) continue;
      if(electron->gsfTrack()->trackerExpectedHitsInner().numberOfHits()!=0) continue;
    }
    foundMuon=true;
    if(electron->pt()>ptElectron){
      SelectedElectron=electron;
      ptElectron=electron->pt();
    }
  }
  
  //BTAG VETO
  int nbtagsL=0; int nbtagsM=0; int nbtagsT=0;
  if(foundJet) BtagVeto(ak5jetCands, nbtagsL, nbtagsM, nbtagsT, SelectedJet);

  //MU-TAU ANALYSIS
  if(foundJet && foundTauMuTau && foundMuon){
    TMatrixD covMET(2, 2); // PFMET significance matrix
    covMET[0][0] = (metRaw->front() ).getSignificanceMatrix()(0,0);
    covMET[1][0] = (metRaw->front() ).getSignificanceMatrix()(1,0);
    covMET[0][1] = (metRaw->front() ).getSignificanceMatrix()(0,1);
    covMET[1][1] = (metRaw->front() ).getSignificanceMatrix()(1,1);
    std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
    measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, SelectedTauMuTau->p4()));
    measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, SelectedMuon->p4()));
    NSVfitStandaloneAlgorithm algo(measuredTauLeptons, met->begin()->momentum(), covMET, 0);
    algo.addLogM(false);
    algo.integrateMarkovChain();
    if(algo.pt()>0){
      
      TLorentzVector SVFitTauTau; SVFitTauTau.SetPtEtaPhiM(algo.pt(), algo.eta(), algo.phi(), algo.getMass());
      math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedJet->pt(),SelectedJet->eta(),SelectedJet->phi(),massZ);
      TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E()); 
      m_jetPtMuoTau=SelectedJet->pt();
      m_jetEtaMuoTau=SelectedJet->eta();
      m_jetMassMuoTau=massZ;
      m_jetSubjettinessMuoTau=tau21Z;
      m_dPhiJetMetMuoTau=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet->p4());
      m_dRJetMetMuoTau=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedJet->p4());
      m_dRJetLepMuoTau=ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedJet->p4());
      m_dRJetTauMuoTau=ROOT::Math::VectorUtil::DeltaR(SelectedTauMuTau->p4(),SelectedJet->p4());
      m_dPhiTauMetMuoTau=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedTauMuTau->p4());
      m_dRTauMetMuoTau=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedTauMuTau->p4());
      m_dRTauLepMuoTau=ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedTauMuTau->p4());
      m_dPhiLepMetMuoTau=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuon->p4());
      m_dRLepMetMuoTau=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedMuon->p4());
      m_dRZZMuoTau=ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,SelectedJet->p4());
      m_tauPtMuoTau=SelectedTauMuTau->pt();
      m_tauEtaMuoTau=SelectedTauMuTau->eta();
      m_lepPtMuoTau=SelectedMuon->pt();
      m_lepEtaMuoTau=SelectedMuon->eta();
      m_lepPFIsoMuoTau=MuonPFIso(SelectedMuon,true);
      m_metMuoTau=met->begin()->pt();
      m_MassSvfitTauLepMuoTau=algo.getMass();
      m_XMassSVFitMuoTau=(SVFitTauTau+PrunedJet).M();
      m_nbtagsLMuoTau=nbtagsL;
      m_nbtagsMMuoTau=nbtagsM;
      m_nbtagsTMuoTau=nbtagsT;
      m_trigger320MuoTau=(int)isFired_HLT_PFJet320;
      m_trigger650MuoTau=(int)isFired_HLT_HT650;
      m_NVerticesMuoTau=vertices->size();
      m_PUWeightMuoTau=MyWeight;
      m_metPxMuoTau=met->begin()->px();
      m_metPyMuoTau=met->begin()->py();
      m_metEtaMuoTau=met->begin()->eta();
      m_metPhiMuoTau=met->begin()->phi();
      TreeMuoTau->Fill();
    }
  }


  //ELE-TAU ANALYSIS
  if(foundJet && foundTauElTau && foundElectron){
    TMatrixD covMET(2, 2); // PFMET significance matrix
    covMET[0][0] = (metRaw->front() ).getSignificanceMatrix()(0,0);
    covMET[1][0] = (metRaw->front() ).getSignificanceMatrix()(1,0);
    covMET[0][1] = (metRaw->front() ).getSignificanceMatrix()(0,1);
    covMET[1][1] = (metRaw->front() ).getSignificanceMatrix()(1,1);
    std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
    measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, SelectedTauElTau->p4()));
    measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, SelectedElectron->p4()));
    NSVfitStandaloneAlgorithm algo(measuredTauLeptons, met->begin()->momentum(), covMET, 0);
    algo.addLogM(false);
    algo.integrateMarkovChain();
    if(algo.pt()>0){
      
      TLorentzVector SVFitTauTau; SVFitTauTau.SetPtEtaPhiM(algo.pt(), algo.eta(), algo.phi(), algo.getMass());
      math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedJet->pt(),SelectedJet->eta(),SelectedJet->phi(),massZ);
      TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E()); 
      m_jetPtEleTau=SelectedJet->pt();
      m_jetEtaEleTau=SelectedJet->eta();
      m_jetMassEleTau=massZ;
      m_jetSubjettinessEleTau=tau21Z;
      m_dPhiJetMetEleTau=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet->p4());
      m_dRJetMetEleTau=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedJet->p4());
      m_dRJetLepEleTau=ROOT::Math::VectorUtil::DeltaR(SelectedElectron->p4(),SelectedJet->p4());
      m_dRJetTauEleTau=ROOT::Math::VectorUtil::DeltaR(SelectedTauElTau->p4(),SelectedJet->p4());
      m_dPhiTauMetEleTau=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedTauElTau->p4());
      m_dRTauMetEleTau=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedTauElTau->p4());
      m_dRTauLepEleTau=ROOT::Math::VectorUtil::DeltaR(SelectedElectron->p4(),SelectedTauElTau->p4());
      m_dPhiLepMetEleTau=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedElectron->p4());
      m_dRLepMetEleTau=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedElectron->p4());
      m_dRZZEleTau=ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,SelectedJet->p4());
      m_tauPtEleTau=SelectedTauElTau->pt();
      m_tauEtaEleTau=SelectedTauElTau->eta();
      m_lepPtEleTau=SelectedMuon->pt();
      m_lepEtaEleTau=SelectedMuon->eta();
      m_lepPFIsoEleTau=MuonPFIso(SelectedMuon,true);
      m_metEleTau=met->begin()->pt();
      m_MassSvfitTauLepEleTau=algo.getMass();
      m_XMassSVFitEleTau=(SVFitTauTau+PrunedJet).M();
      m_nbtagsLEleTau=nbtagsL;
      m_nbtagsMEleTau=nbtagsM;
      m_nbtagsTEleTau=nbtagsT;
      m_trigger320EleTau=(int)isFired_HLT_PFJet320;
      m_trigger650EleTau=(int)isFired_HLT_HT650;
      m_NVerticesEleTau=vertices->size();
      m_PUWeightEleTau=MyWeight;
      m_metPxEleTau=met->begin()->px();
      m_metPyEleTau=met->begin()->py();
      m_metEtaEleTau=met->begin()->eta();
      m_metPhiEleTau=met->begin()->phi();
      TreeEleTau->Fill();
    }
  }

#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
SemiLeptonicAnalyzer::beginJob()
{
  Service<TFileService> fs;
  Nevents = fs->make<TH1D>("Nevents", "Nevents", 3, -0.5, 2.5);

  TreeMuoTau = fs->make<TTree>("TreeMuoTau", "TreeMuoTau");
  TreeMuoTau->Branch("jetPtMuoTau", &m_jetPtMuoTau, "jetPtMuoTau/f");
  TreeMuoTau->Branch("jetEtaMuoTau", &m_jetEtaMuoTau, "jetEtaMuoTau/f");
  TreeMuoTau->Branch("jetMassMuoTau", &m_jetMassMuoTau, "jetMassMuoTau/f");
  TreeMuoTau->Branch("jetSubjettinessMuoTau", &m_jetSubjettinessMuoTau, "jetSubjettinessMuoTau/f");
  TreeMuoTau->Branch("dPhiJetMetMuoTau", &m_dPhiJetMetMuoTau, "dPhiJetMetMuoTau/f");
  TreeMuoTau->Branch("dRJetMetMuoTau", &m_dRJetMetMuoTau, "dRJetMetMuoTau/f");
  TreeMuoTau->Branch("dRJetLepMuoTau", &m_dRJetLepMuoTau, "dRJetLepMuoTau/f");
  TreeMuoTau->Branch("dRJetTauMuoTau", &m_dRJetTauMuoTau, "dRJetTauMuoTau/f");
  TreeMuoTau->Branch("dPhiTauMetMuoTau", &m_dPhiTauMetMuoTau, "dPhiTauMetMuoTau/f");
  TreeMuoTau->Branch("dRTauMetMuoTau", &m_dRTauMetMuoTau, "dRTauMetMuoTau/f");
  TreeMuoTau->Branch("dRTauLepMuoTau", &m_dRTauLepMuoTau, "dRTauLepMuoTau/f");
  TreeMuoTau->Branch("dPhiLepMetMuoTau", &m_dPhiLepMetMuoTau, "dPhiLepMetMuoTau/f");
  TreeMuoTau->Branch("dRLepMetMuoTau", &m_dRLepMetMuoTau, "dRLepMetMuoTau/f");
  TreeMuoTau->Branch("dRZZMuoTau", &m_dRZZMuoTau, "dRZZMuoTau/f");
  TreeMuoTau->Branch("tauPtMuoTau", &m_tauPtMuoTau, "tauPtMuoTau/f");
  TreeMuoTau->Branch("tauEtaMuoTau", &m_tauEtaMuoTau, "tauEtaMuoTau/f");
  TreeMuoTau->Branch("lepPtMuoTau", &m_lepPtMuoTau, "lepPtMuoTau/f");
  TreeMuoTau->Branch("lepEtaMuoTau", &m_lepEtaMuoTau, "lepEtaMuoTau/f");
  TreeMuoTau->Branch("lepPFIsoMuoTau", &m_lepPFIsoMuoTau, "lepPFIsoMuoTau/f");
  TreeMuoTau->Branch("metMuoTau", &m_metMuoTau, "metMuoTau/f");
  TreeMuoTau->Branch("MassSvfitTauLepMuoTau", &m_MassSvfitTauLepMuoTau, "MassSvfitTauLepMuoTau/f");
  TreeMuoTau->Branch("XMassSVFitMuoTau", &m_XMassSVFitMuoTau, "XMassSVFitMuoTau/f");
  TreeMuoTau->Branch("nbtagsLMuoTau", &m_nbtagsLMuoTau, "nbtagsLMuoTau/i");
  TreeMuoTau->Branch("nbtagsMMuoTau", &m_nbtagsMMuoTau, "nbtagsMMuoTau/i");
  TreeMuoTau->Branch("nbtagsTMuoTau", &m_nbtagsTMuoTau, "nbtagsTMuoTau/i");
  TreeMuoTau->Branch("trigger320MuoTau", &m_trigger320MuoTau, "trigger320MuoTau/i");
  TreeMuoTau->Branch("trigger650MuoTau", &m_trigger650MuoTau, "trigger650MuoTau/i");
  TreeMuoTau->Branch("NVerticesMuoTau", &m_NVerticesMuoTau, "NVerticesMuoTau/i");
  TreeMuoTau->Branch("PUWeightMuoTau", &m_PUWeightMuoTau, "PUWeightMuoTau/f");
  TreeMuoTau->Branch("metPxMuoTau", &m_metPxMuoTau, "metPxMuoTau/f");
  TreeMuoTau->Branch("metPyMuoTau", &m_metPyMuoTau, "metPyMuoTau/f");
  TreeMuoTau->Branch("metEtaMuoTau", &m_metEtaMuoTau, "metEtaMuoTau/f");
  TreeMuoTau->Branch("metPhiMuoTau", &m_metPhiMuoTau, "metPhiMuoTau/f");

  TreeEleTau = fs->make<TTree>("TreeEleTau", "TreeEleTau");
  TreeEleTau->Branch("jetPtEleTau", &m_jetPtEleTau, "jetPtEleTau/f");
  TreeEleTau->Branch("jetEtaEleTau", &m_jetEtaEleTau, "jetEtaEleTau/f");
  TreeEleTau->Branch("jetMassEleTau", &m_jetMassEleTau, "jetMassEleTau/f");
  TreeEleTau->Branch("jetSubjettinessEleTau", &m_jetSubjettinessEleTau, "jetSubjettinessEleTau/f");
  TreeEleTau->Branch("dPhiJetMetEleTau", &m_dPhiJetMetEleTau, "dPhiJetMetEleTau/f");
  TreeEleTau->Branch("dRJetMetEleTau", &m_dRJetMetEleTau, "dRJetMetEleTau/f");
  TreeEleTau->Branch("dRJetLepEleTau", &m_dRJetLepEleTau, "dRJetLepEleTau/f");
  TreeEleTau->Branch("dRJetTauEleTau", &m_dRJetTauEleTau, "dRJetTauEleTau/f");
  TreeEleTau->Branch("dPhiTauMetEleTau", &m_dPhiTauMetEleTau, "dPhiTauMetEleTau/f");
  TreeEleTau->Branch("dRTauMetEleTau", &m_dRTauMetEleTau, "dRTauMetEleTau/f");
  TreeEleTau->Branch("dRTauLepEleTau", &m_dRTauLepEleTau, "dRTauLepEleTau/f");
  TreeEleTau->Branch("dPhiLepMetEleTau", &m_dPhiLepMetEleTau, "dPhiLepMetEleTau/f");
  TreeEleTau->Branch("dRLepMetEleTau", &m_dRLepMetEleTau, "dRLepMetEleTau/f");
  TreeEleTau->Branch("dRZZEleTau", &m_dRZZEleTau, "dRZZEleTau/f");
  TreeEleTau->Branch("tauPtEleTau", &m_tauPtEleTau, "tauPtEleTau/f");
  TreeEleTau->Branch("tauEtaEleTau", &m_tauEtaEleTau, "tauEtaEleTau/f");
  TreeEleTau->Branch("lepPtEleTau", &m_lepPtEleTau, "lepPtEleTau/f");
  TreeEleTau->Branch("lepEtaEleTau", &m_lepEtaEleTau, "lepEtaEleTau/f");
  TreeEleTau->Branch("lepPFIsoEleTau", &m_lepPFIsoEleTau, "lepPFIsoEleTau/f");
  TreeEleTau->Branch("metEleTau", &m_metEleTau, "metEleTau/f");
  TreeEleTau->Branch("MassSvfitTauLepEleTau", &m_MassSvfitTauLepEleTau, "MassSvfitTauLepEleTau/f");
  TreeEleTau->Branch("XMassSVFitEleTau", &m_XMassSVFitEleTau, "XMassSVFitEleTau/f");
  TreeEleTau->Branch("nbtagsLEleTau", &m_nbtagsLEleTau, "nbtagsLEleTau/i");
  TreeEleTau->Branch("nbtagsMEleTau", &m_nbtagsMEleTau, "nbtagsMEleTau/i");
  TreeEleTau->Branch("nbtagsTEleTau", &m_nbtagsTEleTau, "nbtagsTEleTau/i");
  TreeEleTau->Branch("trigger320EleTau", &m_trigger320EleTau, "trigger320EleTau/i");
  TreeEleTau->Branch("trigger650EleTau", &m_trigger650EleTau, "trigger650EleTau/i");
  TreeEleTau->Branch("NVerticesEleTau", &m_NVerticesEleTau, "NVerticesEleTau/i");
  TreeEleTau->Branch("PUWeightEleTau", &m_PUWeightEleTau, "PUWeightEleTau/f");
  TreeEleTau->Branch("metPxEleTau", &m_metPxEleTau, "metPxEleTau/f");
  TreeEleTau->Branch("metPyEleTau", &m_metPyEleTau, "metPyEleTau/f");
  TreeEleTau->Branch("metEtaEleTau", &m_metEtaEleTau, "metEtaEleTau/f");
  TreeEleTau->Branch("metPhiEleTau", &m_metPhiEleTau, "metPhiEleTau/f");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SemiLeptonicAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
SemiLeptonicAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
SemiLeptonicAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
SemiLeptonicAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
SemiLeptonicAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SemiLeptonicAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

float SemiLeptonicAnalyzer::MuonPFIso(pat::MuonCollection::const_iterator muon, bool highpt){
  float sumChargedHadronPt = muon->pfIsolationR04().sumChargedHadronPt;
  float sumNeutralHadronEt = muon->pfIsolationR04().sumNeutralHadronEt;
  float sumPhotonEt = muon->pfIsolationR04().sumPhotonEt;
  float sumPUPt = muon->pfIsolationR04().sumPUPt;
  float iso = (sumChargedHadronPt+ max(0.,sumNeutralHadronEt+sumPhotonEt-0.5*sumPUPt))/muon->pt();
  if(highpt){
    reco::TrackRef cktTrack = (muon::tevOptimized(*muon, 200, 17., 40., 0.25)).first;
    iso = (sumChargedHadronPt+ max(0.,sumNeutralHadronEt+sumPhotonEt-0.5*sumPUPt))/cktTrack->pt();
  }
  return iso;
}

void SemiLeptonicAnalyzer::BtagVeto(edm::Handle<pat::JetCollection> ak5jetCands, int & nbtagsL, int & nbtagsM, int & nbtagsT,
				     pat::JetCollection::const_iterator SelectedJet){
  for(pat::JetCollection::const_iterator ak5 = ak5jetCands->begin(); ak5 != ak5jetCands->end(); ++ak5) {
    if(ROOT::Math::VectorUtil::DeltaR(ak5->p4(),SelectedJet->p4())<0.8) continue;
    double discCSV = ak5->bDiscriminator("combinedSecondaryVertexBJetTags");
    if(discCSV>0.244) nbtagsL++; //loose working point
    if(discCSV>0.679) nbtagsM++; //medium working point
    if(discCSV>0.898) nbtagsT++; //tight working point
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(SemiLeptonicAnalyzer);
