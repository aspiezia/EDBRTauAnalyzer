// -*- C++ -*-
//
// Package:    QCDAnalyzer
// Class:      QCDAnalyzer
// 
/**\class QCDAnalyzer QCDAnalyzer.cc ExoDiBosonResonances/EDBRTauAnalyzer/src/QCDAnalyzer.cc

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

class QCDAnalyzer : public edm::EDAnalyzer {
public:
  explicit QCDAnalyzer(const edm::ParameterSet&);
  ~QCDAnalyzer();
  
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
  float ElectronPFIso(pat::ElectronCollection::const_iterator electron, float rho);
  float MuonCorrPFIso(pat::MuonCollection::const_iterator muon, bool highpt);
  float ElectronCorrPFIso(pat::ElectronCollection::const_iterator electron, float rho);
  void BtagVeto(edm::Handle<pat::JetCollection> ak5jetCands, int & nbtagsL, int & nbtagsM, int & nbtagsT, pat::JetCollection::const_iterator SelectedJet);

  TH1D* Nevents;

  TTree *TreeSignalEff;
  float m_genEvent;

  TTree *TreeMuoTau;
  TTree *TreeEleTau;
  float m_jetPt;
  float m_jetEta;
  float m_jetMass;
  float m_jetSubjettiness;
  float m_dPhiJetMet;
  float m_dRJetMet;
  float m_dRJetLep;
  float m_dRJetTau;
  float m_dPhiTauMet;
  float m_dRTauMet;
  float m_dRTauLep;
  float m_dPhiLepMet;
  float m_dRLepMet;
  float m_dRZZVis;
  float m_dRZZEff;
  float m_dRZZSvFit;
  float m_dRZZCA;
  float m_tauPt;
  float m_tauEta;
  float m_tauCharge;
  float m_tauIso;
  int m_tauFake;
  float m_lepPt;
  float m_lepEta;
  float m_lepCharge;
  float m_lepPFIso;
  float m_lepCorrPFIso;
  float m_charge;
  float m_met;
  float m_metPhi;
  float m_uncorrmet;
  float m_uncorrmetPhi;
  float m_MassVis;
  float m_MassEff;
  float m_MassSvfit;
  float m_MassCA;
  float m_XMassVis;
  float m_XMassEff;
  float m_XMassSVFit;
  float m_XMassCA;
  int m_nbtagsL;
  int m_nbtagsM;
  int m_nbtagsT;
  int m_trigger;
  int m_trigger650;
  int m_trigger320;
  int m_sideband;
  int m_NVertices;
  float m_PUWeight;
  float m_metPx;
  float m_metPy;
  int m_NeventsTOT;
  double m_xsec;
  double m_lumi;
  float m_weight;

  edm::LumiReWeighting LumiWeights_;
  bool isData; 
  edm::InputTag vtxColl_;
  edm::InputTag jetColl_;
  edm::InputTag jetPrunedColl_;
  edm::InputTag electronColl_;
  edm::InputTag muonColl_;
  edm::InputTag tauMuTauColl_;
  edm::InputTag tauElTauColl_;
  edm::InputTag metColl_;
  edm::InputTag metRawColl_;
  edm::InputTag uncorrmetColl_;
  edm::InputTag ak5JetColl_;
  int NeventsTOT_;
  double xsec_;
  double lumi_;

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
QCDAnalyzer::QCDAnalyzer(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
  isData = iConfig.getUntrackedParameter<bool>("isData_");
  vtxColl_ = iConfig.getParameter<edm::InputTag>("vtxColl"); 
  jetColl_ = iConfig.getParameter<edm::InputTag>("jetColl"); 
  jetPrunedColl_ = iConfig.getParameter<edm::InputTag>("jetPrunedColl"); 
  electronColl_ = iConfig.getParameter<edm::InputTag>("electronColl"); 
  muonColl_ = iConfig.getParameter<edm::InputTag>("muonColl"); 
  tauMuTauColl_ = iConfig.getParameter<edm::InputTag>("tauMuTauColl");
  tauElTauColl_ = iConfig.getParameter<edm::InputTag>("tauElTauColl"); 
  metColl_ = iConfig.getParameter<edm::InputTag>("metColl"); 
  metRawColl_ = iConfig.getParameter<edm::InputTag>("metRawColl");
  uncorrmetColl_ = iConfig.getParameter<edm::InputTag>("uncorrmetColl");  
  ak5JetColl_ = iConfig.getParameter<edm::InputTag>("ak5JetColl"); 
  NeventsTOT_ = iConfig.getParameter<int>( "NeventsTOT" );
  xsec_= iConfig.getParameter<double>( "xsec" );
  lumi_= iConfig.getParameter<double>( "lumi" );


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


QCDAnalyzer::~QCDAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
QCDAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  edm::Handle<pat::METCollection> uncorrmet;
  iEvent.getByLabel(uncorrmetColl_, uncorrmet);

  edm::Handle<pat::JetCollection> ak5jetCands;
  iEvent.getByLabel(ak5JetColl_,ak5jetCands);

  edm::Handle<double> rhoHandle;
  iEvent.getByLabel("kt6PFJets", "rho", rhoHandle);
  float rho = *(rhoHandle.product());

  Handle<vector<reco::GenParticle> > genParts;
  iEvent.getByLabel("genParticles", genParts);

  //FOR SIGNAL EFFICIENCY
  float genEvent = 0;
  if(isData==false){
    int ele = 0; int muo = 0;
    vector<math::PtEtaPhiELorentzVector> genele;
    for(size_t ngenPart=0; ngenPart<genParts->size(); ngenPart++){
      const reco::GenParticle & genPart = (*genParts)[ngenPart];
      if(abs(genPart.pdgId())==15 && genPart.status()!=3){
	for(unsigned int ndaugh=0; ndaugh<genPart.numberOfDaughters(); ndaugh++){
	  const reco::Candidate * daughter = genPart.daughter(ndaugh);
	  if(abs(daughter->pdgId())==11 && daughter->status()==1) ele = ele + 1;
	  if(abs(daughter->pdgId())==13 && daughter->status()==1) muo = muo + 1;
	}
      }
    }
    if(ele==2 && muo==0)      genEvent = 1;
    else if(ele==1 && muo==1) genEvent = 2;
    else if(ele==1 && muo==0) genEvent = 3;
    else if(ele==0 && muo==2) genEvent = 4;
    else if(ele==0 && muo==1) genEvent = 5;
    else if(ele==0 && muo==0) genEvent = 6;
    else                      genEvent = 7;
  }
  m_genEvent = genEvent;
  TreeSignalEff->Fill();

  //TRIGGER
  bool isFired_HLT_HT650 = false;
  bool isFired_HLT_PFJet320 = false;
  bool isFired_HLT = false;
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
  if(isFired_HLT_PFJet320 || isFired_HLT_HT650) isFired_HLT=true;

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
    if(!(mass>20 && mass<110))  continue;
    if(jet->userFloat("tau2")/jet->userFloat("tau1")>0.75) continue;
    foundJet=true;
    if(jet->pt()>ptZ){
      massZ=mass;
      ptZ=jet->pt();
      tau21Z=jet->userFloat("tau2")/jet->userFloat("tau1");
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
    //if(patTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr")<0.5) continue;
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
    //if(patTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr")<0.5) continue;
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
    foundElectron=true;
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
    
    //FROM WHERE DOES IT COME THE RECO TAU?
    int tauFake = -999;
    if(isData==false){ 
      float dR=999; float TauPdgId=-99999;
      for(size_t ngenPart=0; ngenPart<genParts->size(); ngenPart++){
	const reco::GenParticle & genPart = (*genParts)[ngenPart];
	if(genPart.status()==3) continue;
	float dR_GenTau=ROOT::Math::VectorUtil::DeltaR(genPart.p4(),SelectedTauMuTau->p4());
	if(dR_GenTau<dR && dR_GenTau<0.2){
	  dR=dR_GenTau;
	  TauPdgId=genPart.pdgId();
	}
      }
      if(abs(TauPdgId)==15) tauFake=1;
      else if (abs(TauPdgId)==11 || abs(TauPdgId)==13) tauFake=2;
      else tauFake=3;
    }

    //SVFIT
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedJet->pt(),SelectedJet->eta(),SelectedJet->phi(),massZ);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());
    TMatrixD covMET(2, 2); // PFMET significance matrix
    covMET[0][0] = (metRaw->front() ).getSignificanceMatrix()(0,0);
    covMET[1][0] = (metRaw->front() ).getSignificanceMatrix()(1,0);
    covMET[0][1] = (metRaw->front() ).getSignificanceMatrix()(0,1);
    covMET[1][1] = (metRaw->front() ).getSignificanceMatrix()(1,1);
    std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
    measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay, SelectedTauMuTau->p4()));
    measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, SelectedMuon->p4()));
    NSVfitStandaloneAlgorithm algo(measuredTauLeptons, met->begin()->momentum(), covMET, 0);
    algo.addLogM(false);
    algo.integrateMarkovChain();
    float MassSVFit = -1;
    float XMassSVFit = -1;
    float dRJetZSVFit = -1;
    if(algo.pt()>0){
      TLorentzVector SVFitTauTau; SVFitTauTau.SetPtEtaPhiM(algo.pt(), algo.eta(), algo.phi(), algo.getMass());
      dRJetZSVFit=ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,SelectedJet->p4());
      XMassSVFit=(SVFitTauTau+PrunedJet).M();
      MassSVFit=algo.getMass();
    }

    //COLLINEAR APPROXIMATION
    TLorentzVector CATauTau; bool CA = false;
    float a = (SelectedTauMuTau->py()*met->begin()->px()-SelectedTauMuTau->px()*met->begin()->py())/
      (SelectedMuon->px()*SelectedTauMuTau->py()-SelectedMuon->py()*SelectedTauMuTau->px());
    float b = (SelectedMuon->py()*met->begin()->px()-SelectedMuon->px()*met->begin()->py())/
      (SelectedTauMuTau->px()*SelectedMuon->py()-SelectedTauMuTau->py()*SelectedMuon->px());
    if(((1+a)*(1+b))>0) {
      CATauTau.SetPxPyPzE((1+a)*SelectedMuon->px()+(1+b)*SelectedTauMuTau->px(),
			  (1+a)*SelectedMuon->py()+(1+b)*SelectedTauMuTau->py(),
			  (1+a)*SelectedMuon->pz()+(1+b)*SelectedTauMuTau->pz(),
			  (1+a)*SelectedMuon->energy()+(1+b)*SelectedTauMuTau->energy());
      CA=true;
    }
    else CA = false;
    float MassCA = -1;
    float XmassCA = -1.;
    float dRJetZCA = -1.;
    if(CA){
      MassCA = CATauTau.M();
      XmassCA = (CATauTau+PrunedJet).M();
      dRJetZCA = ROOT::Math::VectorUtil::DeltaR(CATauTau,PrunedJet);
    }

    //VISIBLE AND EFFECTIVE
    math::PtEtaPhiELorentzVector dilep; dilep = SelectedTauMuTau->p4()+SelectedMuon->p4();
    math::PtEtaPhiELorentzVector dilepmet; dilepmet = SelectedTauMuTau->p4()+SelectedMuon->p4()+met->begin()->p4();
    math::PtEtaPhiELorentzVector dilepjet; dilepjet = SelectedTauMuTau->p4()+SelectedMuon->p4()+PrunedJet_prov;
    math::PtEtaPhiELorentzVector dilepmetjet; dilepmetjet = SelectedTauMuTau->p4()+SelectedMuon->p4()+met->begin()->p4()+PrunedJet_prov;

    m_jetPt=SelectedJet->pt();
    m_jetEta=SelectedJet->eta();
    m_jetMass=massZ;
    m_jetSubjettiness=tau21Z;
    m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet->p4());
    m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedJet->p4());
    m_dRJetLep=ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedJet->p4());
    m_dRJetTau=ROOT::Math::VectorUtil::DeltaR(SelectedTauMuTau->p4(),SelectedJet->p4());
    m_dPhiTauMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedTauMuTau->p4());
    m_dRTauMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedTauMuTau->p4());
    m_dRTauLep=ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedTauMuTau->p4());
    m_dPhiLepMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuon->p4());
    m_dRLepMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedMuon->p4());
    m_dRZZVis=ROOT::Math::VectorUtil::DeltaR(dilep,PrunedJet);
    m_dRZZEff=ROOT::Math::VectorUtil::DeltaR(dilepmet,PrunedJet);
    m_dRZZSvFit=dRJetZSVFit;
    m_dRZZCA=dRJetZCA;
    m_tauPt=SelectedTauMuTau->pt();
    m_tauEta=SelectedTauMuTau->eta();
    m_tauCharge=SelectedTauMuTau->charge();
    m_tauIso=SelectedTauMuTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr");
    m_tauFake=tauFake;
    m_lepPt=SelectedMuon->pt();
    m_lepEta=SelectedMuon->eta();
    m_lepCharge=SelectedMuon->charge();
    m_lepPFIso=MuonPFIso(SelectedMuon,true);
    m_lepCorrPFIso=MuonCorrPFIso(SelectedMuon,true);
    m_charge=SelectedTauMuTau->charge()*SelectedMuon->charge();
    m_met=met->begin()->pt();
    m_metPhi=met->begin()->phi();
    m_uncorrmet=uncorrmet->begin()->pt();
    m_uncorrmetPhi=uncorrmet->begin()->phi();
    m_MassVis=dilep.mass();
    m_MassEff=dilepmet.mass();
    m_MassSvfit=MassSVFit;
    m_MassCA=MassCA;
    m_XMassVis=dilepjet.mass();
    m_XMassEff=dilepmetjet.mass();
    m_XMassSVFit=XMassSVFit;
    m_XMassCA=XmassCA;
    m_nbtagsL=nbtagsL;
    m_nbtagsM=nbtagsM;
    m_nbtagsT=nbtagsT;
    m_trigger=(int)isFired_HLT;
    m_trigger320=(int)isFired_HLT_PFJet320;
    m_trigger650=(int)isFired_HLT_HT650;
    m_sideband=(int)(massZ<70);
    m_NVertices=vertices->size();
    m_PUWeight=MyWeight;
    m_metPx=met->begin()->px();
    m_metPy=met->begin()->py();
    m_NeventsTOT=NeventsTOT_;
    m_xsec=xsec_;
    m_lumi=lumi_;
    m_weight=xsec_*lumi_/NeventsTOT_;
    m_genEvent = genEvent;
    TreeMuoTau->Fill();
  }


  //ELE-TAU ANALYSIS
  if(foundJet && foundTauElTau && foundElectron){
    
    //FROM WHERE DOES IT COME THE RECO TAU?
    int tauFake = -999;
    if(isData==false){ 
      float dR=999; float TauPdgId=-99999;
      for(size_t ngenPart=0; ngenPart<genParts->size(); ngenPart++){
	const reco::GenParticle & genPart = (*genParts)[ngenPart];
	if(genPart.status()==3) continue;
	float dR_GenTau=ROOT::Math::VectorUtil::DeltaR(genPart.p4(),SelectedTauElTau->p4());
	if(dR_GenTau<dR && dR_GenTau<0.2){
	  dR=dR_GenTau;
	  TauPdgId=genPart.pdgId();
	}
      }
      if(abs(TauPdgId)==15) tauFake=1;
      else if (abs(TauPdgId)==11 || abs(TauPdgId)==13) tauFake=2;
      else tauFake=3;
    }

    //SVFIT
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedJet->pt(),SelectedJet->eta(),SelectedJet->phi(),massZ);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());
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
    float MassSVFit = -1;
    float XMassSVFit = -1;
    float dRJetZSVFit = -1;
    if(algo.pt()>0){
      TLorentzVector SVFitTauTau; SVFitTauTau.SetPtEtaPhiM(algo.pt(), algo.eta(), algo.phi(), algo.getMass());
      dRJetZSVFit=ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,SelectedJet->p4());
      XMassSVFit=(SVFitTauTau+PrunedJet).M();
      MassSVFit=algo.getMass();
    }

    //COLLINEAR APPROXIMATION
    TLorentzVector CATauTau; bool CA = false;
    float a = (SelectedTauElTau->py()*met->begin()->px()-SelectedTauElTau->px()*met->begin()->py())/
      (SelectedElectron->px()*SelectedTauElTau->py()-SelectedElectron->py()*SelectedTauElTau->px());
    float b = (SelectedElectron->py()*met->begin()->px()-SelectedElectron->px()*met->begin()->py())/
      (SelectedTauElTau->px()*SelectedElectron->py()-SelectedTauElTau->py()*SelectedElectron->px());
    if(((1+a)*(1+b))>0) {
      CATauTau.SetPxPyPzE((1+a)*SelectedElectron->px()+(1+b)*SelectedTauElTau->px(),
			  (1+a)*SelectedElectron->py()+(1+b)*SelectedTauElTau->py(),
			  (1+a)*SelectedElectron->pz()+(1+b)*SelectedTauElTau->pz(),
			  (1+a)*SelectedElectron->energy()+(1+b)*SelectedTauElTau->energy());
      CA=true;
    }
    else CA = false;
    float MassCA = -1;
    float XmassCA = -1.;
    float dRJetZCA = -1.;
    if(CA){
      MassCA = CATauTau.M();
      XmassCA = (CATauTau+PrunedJet).M();
      dRJetZCA = ROOT::Math::VectorUtil::DeltaR(CATauTau,PrunedJet);
    }

    //VISIBLE AND EFFECTIVE
    math::PtEtaPhiELorentzVector dilep; dilep = SelectedTauElTau->p4()+SelectedElectron->p4();
    math::PtEtaPhiELorentzVector dilepmet; dilepmet = SelectedTauElTau->p4()+SelectedElectron->p4()+met->begin()->p4();
    math::PtEtaPhiELorentzVector dilepjet; dilepjet = SelectedTauElTau->p4()+SelectedElectron->p4()+PrunedJet_prov;
    math::PtEtaPhiELorentzVector dilepmetjet; dilepmetjet = SelectedTauElTau->p4()+SelectedElectron->p4()+met->begin()->p4()+PrunedJet_prov;

    m_jetPt=SelectedJet->pt();
    m_jetEta=SelectedJet->eta();
    m_jetMass=massZ;
    m_jetSubjettiness=tau21Z;
    m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet->p4());
    m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedJet->p4());
    m_dRJetLep=ROOT::Math::VectorUtil::DeltaR(SelectedElectron->p4(),SelectedJet->p4());
    m_dRJetTau=ROOT::Math::VectorUtil::DeltaR(SelectedTauElTau->p4(),SelectedJet->p4());
    m_dPhiTauMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedTauElTau->p4());
    m_dRTauMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedTauElTau->p4());
    m_dRTauLep=ROOT::Math::VectorUtil::DeltaR(SelectedElectron->p4(),SelectedTauElTau->p4());
    m_dPhiLepMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedElectron->p4());
    m_dRLepMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedElectron->p4());
    m_dRZZVis=ROOT::Math::VectorUtil::DeltaR(dilep,PrunedJet);
    m_dRZZEff=ROOT::Math::VectorUtil::DeltaR(dilepmet,PrunedJet);
    m_dRZZSvFit=dRJetZSVFit;
    m_dRZZCA=dRJetZCA;
    m_tauPt=SelectedTauElTau->pt();
    m_tauEta=SelectedTauElTau->eta();
    m_tauCharge=SelectedTauElTau->charge();
    m_tauIso=SelectedTauElTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr");
    m_tauFake=tauFake;
    m_lepPt=SelectedElectron->pt();
    m_lepEta=SelectedElectron->eta();
    m_lepCharge=SelectedElectron->charge();
    m_lepPFIso=ElectronPFIso(SelectedElectron,rho);
    m_lepCorrPFIso=ElectronCorrPFIso(SelectedElectron,rho);
    m_charge=SelectedTauElTau->charge()*SelectedElectron->charge();
    m_met=met->begin()->pt();
    m_metPhi=met->begin()->phi();
    m_uncorrmet=uncorrmet->begin()->pt();
    m_uncorrmetPhi=uncorrmet->begin()->phi();
    m_MassVis=dilep.mass();
    m_MassEff=dilepmet.mass();
    m_MassSvfit=MassSVFit;
    m_MassCA=MassCA;
    m_XMassVis=dilepjet.mass();
    m_XMassEff=dilepmetjet.mass();
    m_XMassSVFit=XMassSVFit;
    m_XMassCA=XmassCA;
    m_nbtagsL=nbtagsL;
    m_nbtagsM=nbtagsM;
    m_nbtagsT=nbtagsT;
    m_trigger=(int)isFired_HLT;
    m_trigger320=(int)isFired_HLT_PFJet320;
    m_trigger650=(int)isFired_HLT_HT650;
    m_sideband=(int)(massZ<70);
    m_NVertices=vertices->size();
    m_PUWeight=MyWeight;
    m_metPx=met->begin()->px();
    m_metPy=met->begin()->py();
    m_NeventsTOT=NeventsTOT_;
    m_xsec=xsec_;
    m_lumi=lumi_;
    m_weight=xsec_*lumi_/NeventsTOT_;
    m_genEvent = genEvent;
    TreeEleTau->Fill();
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
QCDAnalyzer::beginJob()
{
  Service<TFileService> fs;
  Nevents = fs->make<TH1D>("Nevents", "Nevents", 3, -0.5, 2.5);

  TreeSignalEff = fs->make<TTree>("TreeSignalEff", "TreeSignalEff");
  TreeSignalEff->Branch("genEvent", &m_genEvent, "genEvent/f");

  TreeMuoTau = fs->make<TTree>("TreeMuoTau", "TreeMuoTau");
  TreeMuoTau->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeMuoTau->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeMuoTau->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeMuoTau->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeMuoTau->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeMuoTau->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeMuoTau->Branch("dRJetLep", &m_dRJetLep, "dRJetLep/f");
  TreeMuoTau->Branch("dRJetTau", &m_dRJetTau, "dRJetTau/f");
  TreeMuoTau->Branch("dPhiTauMet", &m_dPhiTauMet, "dPhiTauMet/f");
  TreeMuoTau->Branch("dRTauMet", &m_dRTauMet, "dRTauMet/f");
  TreeMuoTau->Branch("dRTauLep", &m_dRTauLep, "dRTauLep/f");
  TreeMuoTau->Branch("dPhiLepMet", &m_dPhiLepMet, "dPhiLepMet/f");
  TreeMuoTau->Branch("dRLepMet", &m_dRLepMet, "dRLepMet/f");
  TreeMuoTau->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeMuoTau->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeMuoTau->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeMuoTau->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeMuoTau->Branch("tauPt", &m_tauPt, "tauPt/f");
  TreeMuoTau->Branch("tauEta", &m_tauEta, "tauEta/f");
  TreeMuoTau->Branch("tauCharge", &m_tauCharge, "tauCharge/f");
  TreeMuoTau->Branch("tauIso", &m_tauIso, "tauIso/f");
  TreeMuoTau->Branch("tauFake", &m_tauFake, "tauFake/i");
  TreeMuoTau->Branch("lepPt", &m_lepPt, "lepPt/f");
  TreeMuoTau->Branch("lepEta", &m_lepEta, "lepEta/f");
  TreeMuoTau->Branch("lepCharge", &m_lepCharge, "lepCharge/f");
  TreeMuoTau->Branch("lepPFIso", &m_lepPFIso, "lepPFIso/f");
  TreeMuoTau->Branch("lepCorrPFIso", &m_lepCorrPFIso, "lepCorrPFIso/f");
  TreeMuoTau->Branch("charge", &m_charge, "charge/f");
  TreeMuoTau->Branch("met", &m_met, "met/f");
  TreeMuoTau->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeMuoTau->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeMuoTau->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  TreeMuoTau->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeMuoTau->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeMuoTau->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeMuoTau->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeMuoTau->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeMuoTau->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeMuoTau->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeMuoTau->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeMuoTau->Branch("nbtagsL", &m_nbtagsL, "nbtagsL/i");
  TreeMuoTau->Branch("nbtagsM", &m_nbtagsM, "nbtagsM/i");
  TreeMuoTau->Branch("nbtagsT", &m_nbtagsT, "nbtagsT/i");
  TreeMuoTau->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeMuoTau->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeMuoTau->Branch("trigger", &m_trigger, "trigger/i");
  TreeMuoTau->Branch("sideband", &m_sideband, "sideband/i");
  TreeMuoTau->Branch("NVertices", &m_NVertices, "NVertices/i");
  TreeMuoTau->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeMuoTau->Branch("metPx", &m_metPx, "metPx/f");
  TreeMuoTau->Branch("metPy", &m_metPy, "metPy/f");
  TreeMuoTau->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeMuoTau->Branch("xsec", &m_xsec, "xsec/d");
  TreeMuoTau->Branch("lumi", &m_lumi, "lumi/d");
  TreeMuoTau->Branch("weight", &m_weight, "weight/f");
  TreeMuoTau->Branch("genEvent", &m_genEvent, "genEvent/f");

  TreeEleTau = fs->make<TTree>("TreeEleTau", "TreeEleTau");
  TreeEleTau->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeEleTau->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeEleTau->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeEleTau->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeEleTau->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeEleTau->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeEleTau->Branch("dRJetLep", &m_dRJetLep, "dRJetLep/f");
  TreeEleTau->Branch("dRJetTau", &m_dRJetTau, "dRJetTau/f");
  TreeEleTau->Branch("dPhiTauMet", &m_dPhiTauMet, "dPhiTauMet/f");
  TreeEleTau->Branch("dRTauMet", &m_dRTauMet, "dRTauMet/f");
  TreeEleTau->Branch("dRTauLep", &m_dRTauLep, "dRTauLep/f");
  TreeEleTau->Branch("dPhiLepMet", &m_dPhiLepMet, "dPhiLepMet/f");
  TreeEleTau->Branch("dRLepMet", &m_dRLepMet, "dRLepMet/f");
  TreeEleTau->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeEleTau->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeEleTau->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeEleTau->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeEleTau->Branch("tauPt", &m_tauPt, "tauPt/f");
  TreeEleTau->Branch("tauEta", &m_tauEta, "tauEta/f");
  TreeEleTau->Branch("tauCharge", &m_tauCharge, "tauCharge/f");
  TreeEleTau->Branch("tauIso", &m_tauIso, "tauIso/f");
  TreeEleTau->Branch("tauFake", &m_tauFake, "tauFake/i");
  TreeEleTau->Branch("lepPt", &m_lepPt, "lepPt/f");
  TreeEleTau->Branch("lepEta", &m_lepEta, "lepEta/f");
  TreeEleTau->Branch("lepCharge", &m_lepCharge, "lepCharge/f");
  TreeEleTau->Branch("lepPFIso", &m_lepPFIso, "lepPFIso/f");
  TreeEleTau->Branch("lepCorrPFIso", &m_lepCorrPFIso, "lepCorrPFIso/f");
  TreeEleTau->Branch("charge", &m_charge, "charge/f");
  TreeEleTau->Branch("met", &m_met, "met/f");
  TreeEleTau->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeEleTau->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeEleTau->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  TreeEleTau->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeEleTau->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeEleTau->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeEleTau->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeEleTau->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeEleTau->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeEleTau->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeEleTau->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeEleTau->Branch("nbtagsL", &m_nbtagsL, "nbtagsL/i");
  TreeEleTau->Branch("nbtagsM", &m_nbtagsM, "nbtagsM/i");
  TreeEleTau->Branch("nbtagsT", &m_nbtagsT, "nbtagsT/i");
  TreeEleTau->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeEleTau->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeEleTau->Branch("trigger", &m_trigger, "trigger/i");
  TreeEleTau->Branch("sideband", &m_sideband, "sideband/i");
  TreeEleTau->Branch("NVertices", &m_NVertices, "NVertices/i");
  TreeEleTau->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeEleTau->Branch("metPx", &m_metPx, "metPx/f");
  TreeEleTau->Branch("metPy", &m_metPy, "metPy/f");
  TreeEleTau->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeEleTau->Branch("xsec", &m_xsec, "xsec/d");
  TreeEleTau->Branch("lumi", &m_lumi, "lumi/d");
  TreeEleTau->Branch("weight", &m_weight, "weight/f");
  TreeEleTau->Branch("genEvent", &m_genEvent, "genEvent/f");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
QCDAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
QCDAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
QCDAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
QCDAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
QCDAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
QCDAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

float QCDAnalyzer::MuonPFIso(pat::MuonCollection::const_iterator muon, bool highpt){
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


float QCDAnalyzer::ElectronPFIso(pat::ElectronCollection::const_iterator electron, float rho){
  float chargedHadronIso = electron->chargedHadronIso();
  float neutralHadronIso = electron->neutralHadronIso();
  float photonIso = electron->photonIso();
  float thiseta = fabs(electron->superCluster()->eta());
  float Aeff=0.;
  if(thiseta<1.0) Aeff=0.13;
  if(thiseta>=1.0 && thiseta<1.479) Aeff=0.14;
  if(thiseta>=1.479 && thiseta<2.0) Aeff=0.07;
  if(thiseta>=2.0 && thiseta<2.2) Aeff=0.09;
  if(thiseta>=2.2 && thiseta<2.3) Aeff=0.11;
  if(thiseta>=2.3 && thiseta<2.4) Aeff=0.11;
  if(thiseta>=2.4) Aeff=0.14;
  float zero = 0.;
  float iso = (chargedHadronIso + max(zero, neutralHadronIso + photonIso - rho*Aeff))/electron->pt();
  return iso;
}

float QCDAnalyzer::MuonCorrPFIso(pat::MuonCollection::const_iterator muon, bool highpt){
  float sumChargedHadronPt = muon->userIsolation(pat::PfChargedHadronIso);//pfIsolationR04().sumChargedHadronPt;
  float sumNeutralHadronEt = muon->userIsolation(pat::PfNeutralHadronIso);//pfIsolationR04().sumNeutralHadronEt;
  float sumPhotonEt = muon->userIsolation(pat::PfGammaIso);//pfIsolationR04().sumPhotonEt;
  float sumPUPt = muon->userIsolation(pat::User2Iso);//pfIsolationR04().sumPUPt;
  float iso = (sumChargedHadronPt+ max(0.,sumNeutralHadronEt+sumPhotonEt-0.5*sumPUPt))/muon->pt();
  if(highpt){
    reco::TrackRef cktTrack = (muon::tevOptimized(*muon, 200, 17., 40., 0.25)).first;
    iso = (sumChargedHadronPt+ max(0.,sumNeutralHadronEt+sumPhotonEt-0.5*sumPUPt))/cktTrack->pt();
  }
  return iso;
}


float QCDAnalyzer::ElectronCorrPFIso(pat::ElectronCollection::const_iterator electron, float rho){
  float chargedHadronIso = electron->userIsolation(pat::PfChargedHadronIso);//chargedHadronIso();
  float neutralHadronIso = electron->userIsolation(pat::PfNeutralHadronIso);//neutralHadronIso();
  float photonIso = electron->userIsolation(pat::PfGammaIso);//photonIso();
  float thiseta = fabs(electron->superCluster()->eta());
  float Aeff=0.;
  if(thiseta<1.0) Aeff=0.13;
  if(thiseta>=1.0 && thiseta<1.479) Aeff=0.14;
  if(thiseta>=1.479 && thiseta<2.0) Aeff=0.07;
  if(thiseta>=2.0 && thiseta<2.2) Aeff=0.09;
  if(thiseta>=2.2 && thiseta<2.3) Aeff=0.11;
  if(thiseta>=2.3 && thiseta<2.4) Aeff=0.11;
  if(thiseta>=2.4) Aeff=0.14;
  float zero = 0.;
  float iso = (chargedHadronIso + max(zero, neutralHadronIso + photonIso - rho*Aeff))/electron->pt();
  return iso;
}

void QCDAnalyzer::BtagVeto(edm::Handle<pat::JetCollection> ak5jetCands, int & nbtagsL, int & nbtagsM, int & nbtagsT,
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
DEFINE_FWK_MODULE(QCDAnalyzer);
