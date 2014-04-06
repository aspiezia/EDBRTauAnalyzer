// -*- C++ -*-
//
// Package:    FakeRateRatio
// Class:      FakeRateRatio
// 
/**\class FakeRateRatio FakeRateRatio.cc ExoDiBosonResonances/EDBRTauAnalyzer/src/FakeRateRatio.cc

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

class FakeRateRatio : public edm::EDAnalyzer {
public:
  explicit FakeRateRatio(const edm::ParameterSet&);
  ~FakeRateRatio();
  
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

  TTree *TreeMuoTau;
  TTree *TreeEleTau;
  TTree *TreeFRLooseMuoTau;
  TTree *TreeFRTightMuoTau;
  TTree *TreeFRLooseEleTau;
  TTree *TreeFRTightEleTau;
  float m_jetPt;
  float m_jetEta;
  float m_jetMass;
  float m_jetSubjettiness;
  float m_tauPt;
  float m_tauEta;
  float m_tauIso;
  float m_decayModeFinding;
  float m_met;
  float m_metPhi;
  float m_uncorrmet;
  float m_uncorrmetPhi;
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
  float m_lepPt;
  float m_lepEta;
  float m_lepPFIso;
  float m_lepCorrPFIso;
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
FakeRateRatio::FakeRateRatio(const edm::ParameterSet& iConfig)

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


FakeRateRatio::~FakeRateRatio()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
FakeRateRatio::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
    if(!(mass>70 && mass<110))  continue;
    if(jet->userFloat("tau2")/jet->userFloat("tau1")>0.75) continue;
    foundJet=true;
    if(jet->pt()>ptZ){
      massZ=mass;
      ptZ=jet->pt();
      tau21Z=jet->userFloat("tau2")/jet->userFloat("tau1");
      SelectedJet=jet;
    }
  }

  /*
  //JET SELECTION FOR TRIGGERING NECESSARY FOR THE FR MEASUREMENT
  pat::JetCollection::const_iterator SelectedLooseJet; float ptLooseZ=-999; bool foundLooseJet=false;
  float massLooseZ=-9999; float tau21LooseZ=-9999;
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
    foundLooseJet=true;
    if(jet->pt()>ptLooseZ){
      ptLooseZ=jet->pt();
      massLooseZ=mass;
      tau21LooseZ=jet->userFloat("tau2")/jet->userFloat("tau1");
      SelectedLooseJet=jet;
    }
  }

  //TIGHT TAU SELECTION - MUTAU
  float ptTauMuTau=-99; bool foundTauMuTau=false;
  pat::TauCollection::const_iterator SelectedTauMuTau;
  for (pat::TauCollection::const_iterator patTau = tauMuTauHandle->begin(); patTau != tauMuTauHandle->end(); ++patTau ) {
    if(patTau->pt()<20) continue;
    if(abs(patTau->eta())>2.4) continue;
    if(patTau->tauID("decayModeFinding")<0.5) continue;
    if(patTau->tauID("againstMuonLoose")<0.5) continue;
    if(patTau->tauID("againstElectronLoose")<0.5) continue;
    if(patTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr")<0.5) continue;
    if(foundLooseJet){
      m_jetPt=SelectedLooseJet->pt();
      m_jetEta=SelectedLooseJet->eta();
      m_jetMass=massLooseZ;
      m_jetSubjettiness=tau21LooseZ;
      m_met=met->begin()->pt();
      m_metPhi=met->begin()->phi();
      m_uncorrmet=uncorrmet->begin()->pt();
      m_uncorrmetPhi=uncorrmet->begin()->phi();
      m_tauPt=patTau->pt();
      m_tauEta=patTau->eta();
      m_decayModeFinding=patTau->tauID("decayModeFinding");
      m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedLooseJet->p4());
      m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedLooseJet->p4());
      m_dRJetTau=ROOT::Math::VectorUtil::DeltaR(patTau->p4(),SelectedLooseJet->p4());
      m_trigger=(int)isFired_HLT;
      m_trigger320=(int)isFired_HLT_PFJet320;
      m_trigger650=(int)isFired_HLT_HT650;
      m_sideband=(int)(massZ<70);
      m_PUWeight=MyWeight;
      m_NeventsTOT=NeventsTOT_;
      m_xsec=xsec_;
      m_lumi=lumi_;
      m_weight=xsec_*lumi_/NeventsTOT_;
      TreeFRTightMuoTau->Fill();
    }
    if(foundJet) if(ROOT::Math::VectorUtil::DeltaR(patTau->p4(),SelectedJet->p4())<0.8){ continue;}
    foundTauMuTau=true;   
    if(patTau->pt()>ptTauMuTau){
      SelectedTauMuTau=patTau;
      ptTauMuTau=patTau->pt();
    }
  }

  //LOOSE TAU SELECTION - MUTAU
  float ptLooseTauMuTau=-99; bool foundLooseTauMuTau=false;
  pat::TauCollection::const_iterator SelectedLooseTauMuTau;
  for (pat::TauCollection::const_iterator patTau = tauMuTauHandle->begin(); patTau != tauMuTauHandle->end(); ++patTau ) {
    if(patTau->pt()<20) continue;
    if(abs(patTau->eta())>2.4) continue;
    if(patTau->tauID("decayModeFinding")<0.5) continue;
    if(patTau->tauID("againstMuonLoose")<0.5) continue;
    if(patTau->tauID("againstElectronLoose")<0.5) continue;
    if(patTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr")>0.5) continue;
    if(foundLooseJet){   
      m_jetPt=SelectedLooseJet->pt();
      m_jetEta=SelectedLooseJet->eta();
      m_jetMass=massLooseZ;
      m_jetSubjettiness=tau21LooseZ;
      m_met=met->begin()->pt();
      m_metPhi=met->begin()->phi();
      m_uncorrmet=uncorrmet->begin()->pt();
      m_uncorrmetPhi=uncorrmet->begin()->phi();
      m_tauPt=patTau->pt();
      m_tauEta=patTau->eta();
      m_decayModeFinding=patTau->tauID("decayModeFinding");
      m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedLooseJet->p4());
      m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedLooseJet->p4());
      m_dRJetTau=ROOT::Math::VectorUtil::DeltaR(patTau->p4(),SelectedLooseJet->p4());
      m_trigger=(int)isFired_HLT;
      m_trigger320=(int)isFired_HLT_PFJet320;
      m_trigger650=(int)isFired_HLT_HT650;
      m_sideband=(int)(massZ<70);
      m_PUWeight=MyWeight;
      m_NeventsTOT=NeventsTOT_;
      m_xsec=xsec_;
      m_lumi=lumi_;
      m_weight=xsec_*lumi_/NeventsTOT_;
      TreeFRLooseMuoTau->Fill();
    }
    if(foundJet) if(ROOT::Math::VectorUtil::DeltaR(patTau->p4(),SelectedJet->p4())<0.8){ continue;}
    foundLooseTauMuTau=true;
    if(patTau->pt()>ptLooseTauMuTau){
      SelectedLooseTauMuTau=patTau;
      ptLooseTauMuTau=patTau->pt();
    }
  }
  
  //TIGHT TAU SELECTION - ELTAU
  float ptTauElTau=-99; bool foundTauElTau=false;
  pat::TauCollection::const_iterator SelectedTauElTau;
  for (pat::TauCollection::const_iterator patTau = tauElTauHandle->begin(); patTau != tauElTauHandle->end(); ++patTau ) {
    if(patTau->pt()<20) continue;
    if(abs(patTau->eta())>2.4) continue;
    if(patTau->tauID("decayModeFinding")<0.5) continue;
    if(patTau->tauID("againstMuonLoose")<0.5) continue;
    if(patTau->tauID("againstElectronLoose")<0.5) continue;
    if(patTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr")<0.5) continue;
    if(foundLooseJet){
      m_jetPt=SelectedLooseJet->pt();
      m_jetEta=SelectedLooseJet->eta();
      m_jetMass=massLooseZ;
      m_jetSubjettiness=tau21LooseZ;
      m_met=met->begin()->pt();
      m_metPhi=met->begin()->phi();
      m_uncorrmet=uncorrmet->begin()->pt();
      m_uncorrmetPhi=uncorrmet->begin()->phi();
      m_tauPt=patTau->pt();
      m_tauEta=patTau->eta();
      m_decayModeFinding=patTau->tauID("decayModeFinding");
      m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedLooseJet->p4());
      m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedLooseJet->p4());
      m_dRJetTau=ROOT::Math::VectorUtil::DeltaR(patTau->p4(),SelectedLooseJet->p4());
      m_trigger=(int)isFired_HLT;
      m_trigger320=(int)isFired_HLT_PFJet320;
      m_trigger650=(int)isFired_HLT_HT650;
      m_sideband=(int)(massZ<70);
      m_PUWeight=MyWeight;
      m_NeventsTOT=NeventsTOT_;
      m_xsec=xsec_;
      m_lumi=lumi_;
      m_weight=xsec_*lumi_/NeventsTOT_;
      TreeFRTightEleTau->Fill();
    }
    if(foundJet) if(ROOT::Math::VectorUtil::DeltaR(patTau->p4(),SelectedJet->p4())<0.8){ continue;}
    foundTauElTau=true;
    if(patTau->pt()>ptTauElTau){
      SelectedTauElTau=patTau;
      ptTauElTau=patTau->pt();
    }
  }

  //LOOSE TAU SELECTION - ELTAU
  float ptLooseTauElTau=-99; bool foundLooseTauElTau=false;
  pat::TauCollection::const_iterator SelectedLooseTauElTau;
  for (pat::TauCollection::const_iterator patTau = tauElTauHandle->begin(); patTau != tauElTauHandle->end(); ++patTau ) {
    if(patTau->pt()<20) continue;
    if(abs(patTau->eta())>2.4) continue;
    if(patTau->tauID("decayModeFinding")<0.5) continue;
    if(patTau->tauID("againstMuonLoose")<0.5) continue;
    if(patTau->tauID("againstElectronLoose")<0.5) continue;
    if(patTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr")>0.5) continue;
    if(foundLooseJet){
      m_jetPt=SelectedLooseJet->pt();
      m_jetEta=SelectedLooseJet->eta();
      m_jetMass=massLooseZ;
      m_jetSubjettiness=tau21LooseZ;
      m_met=met->begin()->pt();
      m_metPhi=met->begin()->phi();
      m_uncorrmet=uncorrmet->begin()->pt();
      m_uncorrmetPhi=uncorrmet->begin()->phi();
      m_tauPt=patTau->pt();
      m_tauEta=patTau->eta();
      m_decayModeFinding=patTau->tauID("decayModeFinding");
      m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedLooseJet->p4());
      m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedLooseJet->p4());
      m_dRJetTau=ROOT::Math::VectorUtil::DeltaR(patTau->p4(),SelectedLooseJet->p4());
      m_trigger=(int)isFired_HLT;
      m_trigger320=(int)isFired_HLT_PFJet320;
      m_trigger650=(int)isFired_HLT_HT650;
      m_sideband=(int)(massZ<70);
      m_PUWeight=MyWeight;
      m_NeventsTOT=NeventsTOT_;
      m_xsec=xsec_;
      m_lumi=lumi_;
      m_weight=xsec_*lumi_/NeventsTOT_;
      TreeFRLooseEleTau->Fill();
    }
    if(foundJet) if(ROOT::Math::VectorUtil::DeltaR(patTau->p4(),SelectedJet->p4())<0.8){ continue;}
    foundLooseTauElTau=true;
    if(patTau->pt()>ptLooseTauElTau){
      SelectedLooseTauElTau=patTau;
      ptLooseTauElTau=patTau->pt();
    }
  }
  */

  //TAU SELECTION - MUOTAU
  float ptLooseTauMuTau=-99; bool foundLooseTauMuTau=false;
  float ptTightTauMuTau=-99; bool foundTightTauMuTau=false;
  float dRJetTauMuTau = -99;
  pat::TauCollection::const_iterator SelectedLooseTauMuTau;
  pat::TauCollection::const_iterator SelectedTightTauMuTau;
  for (pat::TauCollection::const_iterator patTau = tauMuTauHandle->begin(); patTau != tauMuTauHandle->end(); ++patTau ) {
    if(patTau->pt()<20) continue;
    if(abs(patTau->eta())>2.4) continue;
    if(patTau->tauID("decayModeFinding")<0.5) continue;
    if(patTau->tauID("againstMuonLoose")<0.5) continue;
    if(patTau->tauID("againstElectronLoose")<0.5) continue;
    if(patTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr")>0.5 && ((foundJet==true && ROOT::Math::VectorUtil::DeltaR(patTau->p4(),SelectedJet->p4())>0.8) || 
								       foundJet==false)){
      if(foundJet) dRJetTauMuTau = ROOT::Math::VectorUtil::DeltaR(patTau->p4(),SelectedJet->p4());
      m_dRJetTau = dRJetTauMuTau;
      m_tauPt=patTau->pt();
      m_tauEta=patTau->eta();
      m_met=met->begin()->pt();
      m_metPhi=met->begin()->phi();
      m_uncorrmet=uncorrmet->begin()->pt();
      m_uncorrmetPhi=uncorrmet->begin()->phi();
      m_trigger=(int)isFired_HLT;
      m_trigger320=(int)isFired_HLT_PFJet320;
      m_trigger650=(int)isFired_HLT_HT650;
      m_sideband=(int)(massZ<70);
      m_PUWeight=MyWeight;
      m_NeventsTOT=NeventsTOT_;
      m_xsec=xsec_;
      m_lumi=lumi_;
      m_weight=xsec_*lumi_/NeventsTOT_;
      TreeFRTightMuoTau->Fill();
      foundTightTauMuTau=true;
      if(patTau->pt()>ptTightTauMuTau){
	SelectedTightTauMuTau=patTau;
	ptTightTauMuTau=patTau->pt();
      }
    } else {
      if(foundJet) dRJetTauMuTau = ROOT::Math::VectorUtil::DeltaR(patTau->p4(),SelectedJet->p4());
      m_dRJetTau = dRJetTauMuTau;
      m_tauPt=patTau->pt();
      m_tauEta=patTau->eta();
      m_met=met->begin()->pt();
      m_metPhi=met->begin()->phi();
      m_uncorrmet=uncorrmet->begin()->pt();
      m_uncorrmetPhi=uncorrmet->begin()->phi();
      m_trigger=(int)isFired_HLT;
      m_trigger320=(int)isFired_HLT_PFJet320;
      m_trigger650=(int)isFired_HLT_HT650;
      m_sideband=(int)(massZ<70);
      m_PUWeight=MyWeight;
      m_NeventsTOT=NeventsTOT_;
      m_xsec=xsec_;
      m_lumi=lumi_;
      m_weight=xsec_*lumi_/NeventsTOT_;
      TreeFRLooseMuoTau->Fill();
      foundLooseTauMuTau=true;
      if(patTau->pt()>ptLooseTauMuTau){
	SelectedLooseTauMuTau=patTau;
	ptLooseTauMuTau=patTau->pt();
      }
    }
  }

  //TAU SELECTION - ELETAU
  float ptLooseTauElTau=-99; bool foundLooseTauElTau=false;
  float ptTightTauElTau=-99; bool foundTightTauElTau=false;
  float dRJetTauElTau = -99;
  pat::TauCollection::const_iterator SelectedLooseTauElTau;
  pat::TauCollection::const_iterator SelectedTightTauElTau;
  for (pat::TauCollection::const_iterator patTau = tauElTauHandle->begin(); patTau != tauElTauHandle->end(); ++patTau ) {
    if(patTau->pt()<20) continue;
    if(abs(patTau->eta())>2.4) continue;
    if(patTau->tauID("decayModeFinding")<0.5) continue;
    if(patTau->tauID("againstMuonLoose")<0.5) continue;
    if(patTau->tauID("againstElectronLoose")<0.5) continue;
    if(patTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr")>0.5 && ((foundJet==true && ROOT::Math::VectorUtil::DeltaR(patTau->p4(),SelectedJet->p4())>0.8) || 
								       foundJet==false)){
      if(foundJet) dRJetTauElTau = ROOT::Math::VectorUtil::DeltaR(patTau->p4(),SelectedJet->p4());
      m_dRJetTau = dRJetTauElTau;
      m_tauPt=patTau->pt();
      m_tauEta=patTau->eta();
      m_met=met->begin()->pt();
      m_metPhi=met->begin()->phi();
      m_uncorrmet=uncorrmet->begin()->pt();
      m_uncorrmetPhi=uncorrmet->begin()->phi();
      m_trigger=(int)isFired_HLT;
      m_trigger320=(int)isFired_HLT_PFJet320;
      m_trigger650=(int)isFired_HLT_HT650;
      m_sideband=(int)(massZ<70);
      m_PUWeight=MyWeight;
      m_NeventsTOT=NeventsTOT_;
      m_xsec=xsec_;
      m_lumi=lumi_;
      m_weight=xsec_*lumi_/NeventsTOT_;
      TreeFRTightEleTau->Fill();
      foundTightTauElTau=true;
      if(patTau->pt()>ptTightTauElTau){
	SelectedTightTauElTau=patTau;
	ptTightTauElTau=patTau->pt();
      }
    } else {
      if(foundJet) dRJetTauElTau = ROOT::Math::VectorUtil::DeltaR(patTau->p4(),SelectedJet->p4());
      m_dRJetTau = dRJetTauElTau;
      m_tauPt=patTau->pt();
      m_tauEta=patTau->eta();
      m_met=met->begin()->pt();
      m_metPhi=met->begin()->phi();
      m_uncorrmet=uncorrmet->begin()->pt();
      m_uncorrmetPhi=uncorrmet->begin()->phi();
      m_trigger=(int)isFired_HLT;
      m_trigger320=(int)isFired_HLT_PFJet320;
      m_trigger650=(int)isFired_HLT_HT650;
      m_sideband=(int)(massZ<70);
      m_PUWeight=MyWeight;
      m_NeventsTOT=NeventsTOT_;
      m_xsec=xsec_;
      m_lumi=lumi_;
      m_weight=xsec_*lumi_/NeventsTOT_;
      TreeFRLooseEleTau->Fill();
      foundLooseTauElTau=true;
      if(patTau->pt()>ptLooseTauElTau){
	SelectedLooseTauElTau=patTau;
	ptLooseTauElTau=patTau->pt();
      }
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
    if(foundJet) if(ROOT::Math::VectorUtil::DeltaR(muon->p4(),SelectedJet->p4())<0.8) continue;
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
    if(foundJet) if(ROOT::Math::VectorUtil::DeltaR(electron->p4(),SelectedJet->p4())<0.8) continue;
    foundElectron=true;
    if(electron->pt()>ptElectron){
      SelectedElectron=electron;
      ptElectron=electron->pt();
    }
  }
  
  //BTAG VETO
  int nbtagsL=0; int nbtagsM=0; int nbtagsT=0;
  if(foundJet) BtagVeto(ak5jetCands, nbtagsL, nbtagsM, nbtagsT, SelectedJet);


  if(foundJet && foundLooseTauMuTau && foundMuon){

    //SVFIT
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedJet->pt(),SelectedJet->eta(),SelectedJet->phi(),massZ);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());
    TMatrixD covMET(2, 2); // PFMET significance matrix
    covMET[0][0] = (metRaw->front() ).getSignificanceMatrix()(0,0);
    covMET[1][0] = (metRaw->front() ).getSignificanceMatrix()(1,0);
    covMET[0][1] = (metRaw->front() ).getSignificanceMatrix()(0,1);
    covMET[1][1] = (metRaw->front() ).getSignificanceMatrix()(1,1);
    std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
    measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, SelectedLooseTauMuTau->p4()));
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
    float a = (SelectedLooseTauMuTau->py()*met->begin()->px()-SelectedLooseTauMuTau->px()*met->begin()->py())/
      (SelectedMuon->px()*SelectedLooseTauMuTau->py()-SelectedMuon->py()*SelectedLooseTauMuTau->px());
    float b = (SelectedMuon->py()*met->begin()->px()-SelectedMuon->px()*met->begin()->py())/
      (SelectedLooseTauMuTau->px()*SelectedMuon->py()-SelectedLooseTauMuTau->py()*SelectedMuon->px());
    if(((1+a)*(1+b))>0 && ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedLooseTauMuTau->p4())>0.005) {
      CATauTau.SetPxPyPzE((1+a)*SelectedMuon->px()+(1+b)*SelectedLooseTauMuTau->px(),
			  (1+a)*SelectedMuon->py()+(1+b)*SelectedLooseTauMuTau->py(),
			  (1+a)*SelectedMuon->pz()+(1+b)*SelectedLooseTauMuTau->pz(),
			  (1+a)*SelectedMuon->energy()+(1+b)*SelectedLooseTauMuTau->energy());
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
    math::PtEtaPhiELorentzVector dilep; dilep = SelectedLooseTauMuTau->p4()+SelectedMuon->p4();
    math::PtEtaPhiELorentzVector dilepmet; dilepmet = SelectedLooseTauMuTau->p4()+SelectedMuon->p4()+met->begin()->p4();
    math::PtEtaPhiELorentzVector dilepjet; dilepjet = SelectedLooseTauMuTau->p4()+SelectedMuon->p4()+PrunedJet_prov;
    math::PtEtaPhiELorentzVector dilepmetjet; dilepmetjet = SelectedLooseTauMuTau->p4()+SelectedMuon->p4()+met->begin()->p4()+PrunedJet_prov;
 
    m_jetPt=SelectedJet->pt();
    m_jetEta=SelectedJet->eta();
    m_jetMass=massZ;
    m_jetSubjettiness=tau21Z;
    m_met=met->begin()->pt();
    m_metPhi=met->begin()->phi();
    m_uncorrmet=uncorrmet->begin()->pt();
    m_uncorrmetPhi=uncorrmet->begin()->phi();
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
    m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet->p4());
    m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedJet->p4());
    m_dRJetLep=ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedJet->p4());
    m_dRJetTau=ROOT::Math::VectorUtil::DeltaR(SelectedLooseTauMuTau->p4(),SelectedJet->p4());
    m_dPhiTauMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedLooseTauMuTau->p4());
    m_dRTauMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedLooseTauMuTau->p4());
    m_dRTauLep=ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedLooseTauMuTau->p4());
    m_dPhiLepMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuon->p4());
    m_dRLepMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedMuon->p4());
    m_dRZZVis=ROOT::Math::VectorUtil::DeltaR(dilep,PrunedJet);
    m_dRZZEff=ROOT::Math::VectorUtil::DeltaR(dilepmet,PrunedJet);
    m_dRZZSvFit=dRJetZSVFit;
    m_dRZZCA=dRJetZCA;
    m_tauPt=SelectedLooseTauMuTau->pt();
    m_tauEta=SelectedLooseTauMuTau->eta();
    m_tauIso=SelectedLooseTauMuTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr");
    m_decayModeFinding=SelectedLooseTauMuTau->tauID("decayModeFinding");
    m_lepPt=SelectedMuon->pt();
    m_lepEta=SelectedMuon->eta();
    m_lepPFIso=MuonPFIso(SelectedMuon,true);
    m_lepCorrPFIso=MuonCorrPFIso(SelectedMuon,true);
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
    TreeMuoTau->Fill();
  }


  if(foundJet && foundLooseTauElTau && foundElectron){ 

    //SVFIT
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedJet->pt(),SelectedJet->eta(),SelectedJet->phi(),massZ);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());
    TMatrixD covMET(2, 2); // PFMET significance matrix
    covMET[0][0] = (metRaw->front() ).getSignificanceMatrix()(0,0);
    covMET[1][0] = (metRaw->front() ).getSignificanceMatrix()(1,0);
    covMET[0][1] = (metRaw->front() ).getSignificanceMatrix()(0,1);
    covMET[1][1] = (metRaw->front() ).getSignificanceMatrix()(1,1);
    std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
    measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, SelectedLooseTauElTau->p4()));
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
    float a = (SelectedLooseTauElTau->py()*met->begin()->px()-SelectedLooseTauElTau->px()*met->begin()->py())/
      (SelectedElectron->px()*SelectedLooseTauElTau->py()-SelectedElectron->py()*SelectedLooseTauElTau->px());
    float b = (SelectedElectron->py()*met->begin()->px()-SelectedElectron->px()*met->begin()->py())/
      (SelectedLooseTauElTau->px()*SelectedElectron->py()-SelectedLooseTauElTau->py()*SelectedElectron->px());
    if(((1+a)*(1+b))>0 && ROOT::Math::VectorUtil::DeltaR(SelectedElectron->p4(),SelectedLooseTauElTau->p4())>0.005) {
      CATauTau.SetPxPyPzE((1+a)*SelectedElectron->px()+(1+b)*SelectedLooseTauElTau->px(),
			  (1+a)*SelectedElectron->py()+(1+b)*SelectedLooseTauElTau->py(),
			  (1+a)*SelectedElectron->pz()+(1+b)*SelectedLooseTauElTau->pz(),
			  (1+a)*SelectedElectron->energy()+(1+b)*SelectedLooseTauElTau->energy());
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
    math::PtEtaPhiELorentzVector dilep; dilep = SelectedLooseTauElTau->p4()+SelectedElectron->p4();
    math::PtEtaPhiELorentzVector dilepmet; dilepmet = SelectedLooseTauElTau->p4()+SelectedElectron->p4()+met->begin()->p4();
    math::PtEtaPhiELorentzVector dilepjet; dilepjet = SelectedLooseTauElTau->p4()+SelectedElectron->p4()+PrunedJet_prov;
    math::PtEtaPhiELorentzVector dilepmetjet; dilepmetjet = SelectedLooseTauElTau->p4()+SelectedElectron->p4()+met->begin()->p4()+PrunedJet_prov;

    m_jetPt=SelectedJet->pt();
    m_jetEta=SelectedJet->eta();
    m_jetMass=massZ;
    m_jetSubjettiness=tau21Z;
    m_met=met->begin()->pt();
    m_metPhi=met->begin()->phi();
    m_uncorrmet=uncorrmet->begin()->pt();
    m_uncorrmetPhi=uncorrmet->begin()->phi();
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
    m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet->p4());
    m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedJet->p4());
    m_dRJetLep=ROOT::Math::VectorUtil::DeltaR(SelectedElectron->p4(),SelectedJet->p4());
    m_dRJetTau=ROOT::Math::VectorUtil::DeltaR(SelectedLooseTauElTau->p4(),SelectedJet->p4());
    m_dPhiTauMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedLooseTauElTau->p4());
    m_dRTauMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedLooseTauElTau->p4());
    m_dRTauLep=ROOT::Math::VectorUtil::DeltaR(SelectedElectron->p4(),SelectedLooseTauElTau->p4());
    m_dPhiLepMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedElectron->p4());
    m_dRLepMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedElectron->p4());
    m_dRZZVis=ROOT::Math::VectorUtil::DeltaR(dilep,PrunedJet);
    m_dRZZEff=ROOT::Math::VectorUtil::DeltaR(dilepmet,PrunedJet);
    m_dRZZSvFit=dRJetZSVFit;
    m_dRZZCA=dRJetZCA;
    m_tauPt=SelectedLooseTauElTau->pt();
    m_tauEta=SelectedLooseTauElTau->eta();
    m_tauIso=SelectedLooseTauElTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr");
    m_decayModeFinding=SelectedLooseTauElTau->tauID("decayModeFinding");
    m_lepPt=SelectedElectron->pt();
    m_lepEta=SelectedElectron->eta();
    m_lepPFIso=ElectronPFIso(SelectedElectron,rho);
    m_lepCorrPFIso=ElectronCorrPFIso(SelectedElectron,rho);
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
FakeRateRatio::beginJob()
{
  Service<TFileService> fs;
  Nevents = fs->make<TH1D>("Nevents", "Nevents", 3, -0.5, 2.5);

  TreeMuoTau = fs->make<TTree>("TreeMuoTau", "TreeMuoTau");
  TreeMuoTau->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeMuoTau->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeMuoTau->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeMuoTau->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeMuoTau->Branch("tauPt", &m_tauPt, "tauPt/f");
  TreeMuoTau->Branch("tauEta", &m_tauEta, "tauEta/f");
  TreeMuoTau->Branch("tauIso", &m_tauIso, "tauIso/f");
  TreeMuoTau->Branch("decayModeFinding", &m_decayModeFinding, "decayModeFinding/f");
  TreeMuoTau->Branch("met", &m_met, "met/f");
  TreeMuoTau->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeMuoTau->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeMuoTau->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
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
  TreeMuoTau->Branch("lepPt", &m_lepPt, "lepPt/f");
  TreeMuoTau->Branch("lepEta", &m_lepEta, "lepEta/f");
  TreeMuoTau->Branch("lepPFIso", &m_lepPFIso, "lepPFIso/f");
  TreeMuoTau->Branch("lepCorrPFIso", &m_lepCorrPFIso, "lepCorrPFIso/f");
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

  TreeEleTau = fs->make<TTree>("TreeEleTau", "TreeEleTau");
  TreeEleTau->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeEleTau->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeEleTau->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeEleTau->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeEleTau->Branch("tauPt", &m_tauPt, "tauPt/f");
  TreeEleTau->Branch("decayModeFinding", &m_decayModeFinding, "decayModeFinding/f");
  TreeEleTau->Branch("tauEta", &m_tauEta, "tauEta/f");
  TreeEleTau->Branch("tauIso", &m_tauIso, "tauIso/f");
  TreeEleTau->Branch("met", &m_met, "met/f");
  TreeEleTau->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeEleTau->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeEleTau->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
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
  TreeEleTau->Branch("lepPt", &m_lepPt, "lepPt/f");
  TreeEleTau->Branch("lepEta", &m_lepEta, "lepEta/f");
  TreeEleTau->Branch("lepPFIso", &m_lepPFIso, "lepPFIso/f");
  TreeEleTau->Branch("lepCorrPFIso", &m_lepCorrPFIso, "lepCorrPFIso/f");
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

  TreeFRTightMuoTau = fs->make<TTree>("TreeFRTightMuoTau", "TreeFRTightMuoTau");
  //TreeFRTightMuoTau->Branch("jetPt", &m_jetPt, "jetPt/f");
  //TreeFRTightMuoTau->Branch("jetEta", &m_jetEta, "jetEta/f");
  //TreeFRTightMuoTau->Branch("jetMass", &m_jetMass, "jetMass/f");
  //TreeFRTightMuoTau->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeFRTightMuoTau->Branch("tauPt", &m_tauPt, "tauPt/f");
  TreeFRTightMuoTau->Branch("tauEta", &m_tauEta, "tauEta/f");
  //TreeFRTightMuoTau->Branch("decayModeFinding", &m_decayModeFinding, "decayModeFinding/f");
  TreeFRTightMuoTau->Branch("met", &m_met, "met/f");
  TreeFRTightMuoTau->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeFRTightMuoTau->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeFRTightMuoTau->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  //TreeFRTightMuoTau->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  //TreeFRTightMuoTau->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeFRTightMuoTau->Branch("dRJetTau", &m_dRJetTau, "dRJetTau/f");
  TreeFRTightMuoTau->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeFRTightMuoTau->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeFRTightMuoTau->Branch("trigger", &m_trigger, "trigger/i");
  TreeFRTightMuoTau->Branch("sideband", &m_sideband, "sideband/i");
  TreeFRTightMuoTau->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeFRTightMuoTau->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeFRTightMuoTau->Branch("xsec", &m_xsec, "xsec/d");
  TreeFRTightMuoTau->Branch("lumi", &m_lumi, "lumi/d");
  TreeFRTightMuoTau->Branch("weight", &m_weight, "weight/f");

  TreeFRLooseMuoTau = fs->make<TTree>("TreeFRLooseMuoTau", "TreeFRLooseMuoTau");
  //TreeFRLooseMuoTau->Branch("jetPt", &m_jetPt, "jetPt/f");
  //TreeFRLooseMuoTau->Branch("jetEta", &m_jetEta, "jetEta/f");
  //TreeFRLooseMuoTau->Branch("jetMass", &m_jetMass, "jetMass/f");
  //TreeFRLooseMuoTau->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeFRLooseMuoTau->Branch("tauPt", &m_tauPt, "tauPt/f");
  TreeFRLooseMuoTau->Branch("tauEta", &m_tauEta, "tauEta/f");
  //TreeFRLooseMuoTau->Branch("decayModeFinding", &m_decayModeFinding, "decayModeFinding/f");
  TreeFRLooseMuoTau->Branch("met", &m_met, "met/f");
  TreeFRLooseMuoTau->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeFRLooseMuoTau->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeFRLooseMuoTau->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  //TreeFRLooseMuoTau->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  //TreeFRLooseMuoTau->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeFRLooseMuoTau->Branch("dRJetTau", &m_dRJetTau, "dRJetTau/f");
  TreeFRLooseMuoTau->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeFRLooseMuoTau->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeFRLooseMuoTau->Branch("trigger", &m_trigger, "trigger/i");
  TreeFRLooseMuoTau->Branch("sideband", &m_sideband, "sideband/i");
  TreeFRLooseMuoTau->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeFRLooseMuoTau->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeFRLooseMuoTau->Branch("xsec", &m_xsec, "xsec/d");
  TreeFRLooseMuoTau->Branch("lumi", &m_lumi, "lumi/d");
  TreeFRLooseMuoTau->Branch("weight", &m_weight, "weight/f");

  TreeFRTightEleTau = fs->make<TTree>("TreeFRTightEleTau", "TreeFRTightEleTau");
  //TreeFRTightEleTau->Branch("jetPt", &m_jetPt, "jetPt/f");
  //TreeFRTightEleTau->Branch("jetEta", &m_jetEta, "jetEta/f");
  //TreeFRTightEleTau->Branch("jetMass", &m_jetMass, "jetMass/f");
  //TreeFRTightEleTau->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeFRTightEleTau->Branch("tauPt", &m_tauPt, "tauPt/f");
  TreeFRTightEleTau->Branch("tauEta", &m_tauEta, "tauEta/f");
  //TreeFRTightEleTau->Branch("decayModeFinding", &m_decayModeFinding, "decayModeFinding/f");
  TreeFRTightEleTau->Branch("met", &m_met, "met/f");
  TreeFRTightEleTau->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeFRTightEleTau->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeFRTightEleTau->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  TreeFRTightEleTau->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  //TreeFRTightEleTau->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  //TreeFRTightEleTau->Branch("dRJetTau", &m_dRJetTau, "dRJetTau/f");
  TreeFRTightEleTau->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeFRTightEleTau->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeFRTightEleTau->Branch("trigger", &m_trigger, "trigger/i");
  TreeFRTightEleTau->Branch("sideband", &m_sideband, "sideband/i");
  TreeFRTightEleTau->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeFRTightEleTau->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeFRTightEleTau->Branch("xsec", &m_xsec, "xsec/d");
  TreeFRTightEleTau->Branch("lumi", &m_lumi, "lumi/d");
  TreeFRTightEleTau->Branch("weight", &m_weight, "weight/f");

  TreeFRLooseEleTau = fs->make<TTree>("TreeFRLooseEleTau", "TreeFRLooseEleTau");
  //TreeFRLooseEleTau->Branch("jetPt", &m_jetPt, "jetPt/f");
  //TreeFRLooseEleTau->Branch("jetEta", &m_jetEta, "jetEta/f");
  //TreeFRLooseEleTau->Branch("jetMass", &m_jetMass, "jetMass/f");
  //TreeFRLooseEleTau->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeFRLooseEleTau->Branch("tauPt", &m_tauPt, "tauPt/f");
  TreeFRLooseEleTau->Branch("tauEta", &m_tauEta, "tauEta/f");
  //TreeFRLooseEleTau->Branch("decayModeFinding", &m_decayModeFinding, "decayModeFinding/f");
  TreeFRLooseEleTau->Branch("met", &m_met, "met/f");
  TreeFRLooseEleTau->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeFRLooseEleTau->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeFRLooseEleTau->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  //TreeFRLooseEleTau->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  //TreeFRLooseEleTau->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeFRLooseEleTau->Branch("dRJetTau", &m_dRJetTau, "dRJetTau/f");
  TreeFRLooseEleTau->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeFRLooseEleTau->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeFRLooseEleTau->Branch("trigger", &m_trigger, "trigger/i");
  TreeFRLooseEleTau->Branch("sideband", &m_sideband, "sideband/i");
  TreeFRLooseEleTau->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeFRLooseEleTau->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeFRLooseEleTau->Branch("xsec", &m_xsec, "xsec/d");
  TreeFRLooseEleTau->Branch("lumi", &m_lumi, "lumi/d");
  TreeFRLooseEleTau->Branch("weight", &m_weight, "weight/f");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FakeRateRatio::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
FakeRateRatio::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
FakeRateRatio::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
FakeRateRatio::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
FakeRateRatio::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FakeRateRatio::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

float FakeRateRatio::MuonPFIso(pat::MuonCollection::const_iterator muon, bool highpt){
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


float FakeRateRatio::ElectronPFIso(pat::ElectronCollection::const_iterator electron, float rho){
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

float FakeRateRatio::MuonCorrPFIso(pat::MuonCollection::const_iterator muon, bool highpt){
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


float FakeRateRatio::ElectronCorrPFIso(pat::ElectronCollection::const_iterator electron, float rho){
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

void FakeRateRatio::BtagVeto(edm::Handle<pat::JetCollection> ak5jetCands, int & nbtagsL, int & nbtagsM, int & nbtagsT,
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
DEFINE_FWK_MODULE(FakeRateRatio);
