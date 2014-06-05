// -*- C++ -*-
//
// Package:    Trigger
// Class:      Trigger
// 
/**\class Trigger Trigger.cc ExoDiBosonResonances/EDBRTauTrigger/src/Trigger.cc

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

class Trigger : public edm::EDAnalyzer {
public:
  explicit Trigger(const edm::ParameterSet&);
  ~Trigger();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  TTree *Tree;
  int m_run;
  int m_lumi;
  int m_event;
  float m_met;
  int m_trigger650;
  int m_trigger320;
  int m_trigger40;
  int m_trigger80;
  int m_trigger140;
  int m_trigger200;
  int m_trigger260;

  TTree *TreeTaggedJet;
  int m_runTagged;
  int m_lumiTagged;
  int m_eventTagged;
  float m_jetMassTagged;
  float m_jetPtTagged;
  float m_jetTau21Tagged;
  float m_metTagged;
  int m_trigger650Tagged;
  int m_trigger320Tagged;
  int m_trigger40Tagged;
  int m_trigger80Tagged;
  int m_trigger140Tagged;
  int m_trigger200Tagged;
  int m_trigger260Tagged;

  TTree *TreeNOTaggedJet;
  int m_runNO;
  int m_lumiNO;
  int m_eventNO;
  float m_jetMassNO;
  float m_jetPtNO;
  float m_jetTau21NO;
  float m_metNO;
  int m_trigger650NO;
  int m_trigger320NO;
  int m_trigger40NO;
  int m_trigger80NO;
  int m_trigger140NO;
  int m_trigger200NO;
  int m_trigger260NO;

  bool isData; 
  edm::InputTag vtxColl_;
  edm::InputTag jetColl_;
  edm::InputTag jetPrunedColl_;
  edm::InputTag electronETColl_;
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
Trigger::Trigger(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
  isData = iConfig.getUntrackedParameter<bool>("isData_");
  vtxColl_ = iConfig.getParameter<edm::InputTag>("vtxColl"); 
  jetColl_ = iConfig.getParameter<edm::InputTag>("jetColl"); 
  jetPrunedColl_ = iConfig.getParameter<edm::InputTag>("jetPrunedColl"); 
  electronETColl_ = iConfig.getParameter<edm::InputTag>("electronETColl"); 
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
}


Trigger::~Trigger()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Trigger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<pat::JetCollection> CA8JetswithQjets;
  iEvent.getByLabel(jetColl_, CA8JetswithQjets);
  edm::Handle<pat::JetCollection> CA8JetsPruned;
  iEvent.getByLabel(jetPrunedColl_, CA8JetsPruned);

  edm::Handle<pat::METCollection> met;
  iEvent.getByLabel(metColl_, met);

  edm::Handle<pat::METCollection> metRaw;
  iEvent.getByLabel(metRawColl_, metRaw);

  edm::Handle<pat::METCollection> uncorrmet;
  iEvent.getByLabel(uncorrmetColl_, uncorrmet);

  //TRIGGER
  bool isFired_HLT_HT650 = false;
  bool isFired_HLT_PFJet320 = false;
  bool isFired_HLT_PFJet40 = false;
  bool isFired_HLT_PFJet80 = false;
  bool isFired_HLT_PFJet140 = false;
  bool isFired_HLT_PFJet200 = false;
  bool isFired_HLT_PFJet260 = false;
  edm::Handle<edm::TriggerResults> trigResults;
  edm::InputTag trigResultsTag("TriggerResults","","HLT");  
  iEvent.getByLabel(trigResultsTag,trigResults);
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*trigResults);
  //####---####---####---####---####---####---####---####---####---####---####---####---####---####---####---//
  unsigned int TrggIndex_PFHT650_v5( trigNames.triggerIndex("HLT_PFHT650_v5"));
  unsigned int TrggIndex_PFHT650_v6( trigNames.triggerIndex("HLT_PFHT650_v6"));
  unsigned int TrggIndex_PFHT650_v7( trigNames.triggerIndex("HLT_PFHT650_v7"));
  unsigned int TrggIndex_PFHT650_v8( trigNames.triggerIndex("HLT_PFHT650_v8"));
  unsigned int TrggIndex_PFHT650_v9( trigNames.triggerIndex("HLT_PFHT650_v9"));
  unsigned int TrggIndex_PFNoPUHT650_v1( trigNames.triggerIndex("HLT_PFNoPUHT650_v1"));
  unsigned int TrggIndex_PFNoPUHT650_v3( trigNames.triggerIndex("HLT_PFNoPUHT650_v3"));
  unsigned int TrggIndex_PFNoPUHT650_v4( trigNames.triggerIndex("HLT_PFNoPUHT650_v4"));
  //####---####---####---####---####---####---####---####---####---####---####---####---####---####---####---//
  unsigned int TrggIndex_PFJet320_v3( trigNames.triggerIndex("HLT_PFJet320_v3"));
  unsigned int TrggIndex_PFJet320_v4( trigNames.triggerIndex("HLT_PFJet320_v4"));
  unsigned int TrggIndex_PFJet320_v5( trigNames.triggerIndex("HLT_PFJet320_v5"));
  unsigned int TrggIndex_PFJet320_v6( trigNames.triggerIndex("HLT_PFJet320_v6"));
  unsigned int TrggIndex_PFJet320_v8( trigNames.triggerIndex("HLT_PFJet320_v8"));
  unsigned int TrggIndex_PFJet320_v9( trigNames.triggerIndex("HLT_PFJet320_v9"));
  //####---####---####---####---####---####---####---####---####---####---####---####---####---####---####---//
  unsigned int TrggIndex_PFJet40_v3( trigNames.triggerIndex("HLT_PFJet40_v3"));
  unsigned int TrggIndex_PFJet40_v4( trigNames.triggerIndex("HLT_PFJet40_v4"));
  unsigned int TrggIndex_PFJet40_v5( trigNames.triggerIndex("HLT_PFJet40_v5"));
  unsigned int TrggIndex_PFJet40_v6( trigNames.triggerIndex("HLT_PFJet40_v6"));
  unsigned int TrggIndex_PFJet40_v7( trigNames.triggerIndex("HLT_PFJet40_v7"));
  unsigned int TrggIndex_PFJet40_v8( trigNames.triggerIndex("HLT_PFJet40_v8"));
  //####---####---####---####---####---####---####---####---####---####---####---####---####---####---####---//
  unsigned int TrggIndex_PFJet80_v3( trigNames.triggerIndex("HLT_PFJet80_v3"));
  unsigned int TrggIndex_PFJet80_v4( trigNames.triggerIndex("HLT_PFJet80_v4"));
  unsigned int TrggIndex_PFJet80_v5( trigNames.triggerIndex("HLT_PFJet80_v5"));
  unsigned int TrggIndex_PFJet80_v6( trigNames.triggerIndex("HLT_PFJet80_v6"));
  unsigned int TrggIndex_PFJet80_v8( trigNames.triggerIndex("HLT_PFJet80_v8"));
  unsigned int TrggIndex_PFJet80_v9( trigNames.triggerIndex("HLT_PFJet80_v9"));
  //####---####---####---####---####---####---####---####---####---####---####---####---####---####---####---//
  unsigned int TrggIndex_PFJet140_v3( trigNames.triggerIndex("HLT_PFJet140_v3"));
  unsigned int TrggIndex_PFJet140_v4( trigNames.triggerIndex("HLT_PFJet140_v4"));
  unsigned int TrggIndex_PFJet140_v5( trigNames.triggerIndex("HLT_PFJet140_v5"));
  unsigned int TrggIndex_PFJet140_v6( trigNames.triggerIndex("HLT_PFJet140_v6"));
  unsigned int TrggIndex_PFJet140_v8( trigNames.triggerIndex("HLT_PFJet140_v8"));
  unsigned int TrggIndex_PFJet140_v9( trigNames.triggerIndex("HLT_PFJet140_v9"));
  //####---####---####---####---####---####---####---####---####---####---####---####---####---####---####---//
  unsigned int TrggIndex_PFJet200_v3( trigNames.triggerIndex("HLT_PFJet200_v3"));
  unsigned int TrggIndex_PFJet200_v4( trigNames.triggerIndex("HLT_PFJet200_v4"));
  unsigned int TrggIndex_PFJet200_v5( trigNames.triggerIndex("HLT_PFJet200_v5"));
  unsigned int TrggIndex_PFJet200_v6( trigNames.triggerIndex("HLT_PFJet200_v6"));
  unsigned int TrggIndex_PFJet200_v8( trigNames.triggerIndex("HLT_PFJet200_v8"));
  unsigned int TrggIndex_PFJet200_v9( trigNames.triggerIndex("HLT_PFJet200_v9"));
  //####---####---####---####---####---####---####---####---####---####---####---####---####---####---####---//
  unsigned int TrggIndex_PFJet260_v3( trigNames.triggerIndex("HLT_PFJet260_v3"));
  unsigned int TrggIndex_PFJet260_v4( trigNames.triggerIndex("HLT_PFJet260_v4"));
  unsigned int TrggIndex_PFJet260_v5( trigNames.triggerIndex("HLT_PFJet260_v5"));
  unsigned int TrggIndex_PFJet260_v6( trigNames.triggerIndex("HLT_PFJet260_v6"));
  unsigned int TrggIndex_PFJet260_v8( trigNames.triggerIndex("HLT_PFJet260_v8"));
  unsigned int TrggIndex_PFJet260_v9( trigNames.triggerIndex("HLT_PFJet260_v9"));
  //####---####---####---####---####---####---####---####---####---####---####---####---####---####---####---//
  if(TrggIndex_PFHT650_v5 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFHT650_v5);
  if(TrggIndex_PFHT650_v6 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFHT650_v6);
  if(TrggIndex_PFHT650_v7 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFHT650_v7);
  if(TrggIndex_PFHT650_v8 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFHT650_v8);
  if(TrggIndex_PFHT650_v9 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFHT650_v9);
  if(TrggIndex_PFNoPUHT650_v1 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFNoPUHT650_v1);
  if(TrggIndex_PFNoPUHT650_v3 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFNoPUHT650_v3);
  if(TrggIndex_PFNoPUHT650_v4 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFNoPUHT650_v4);
  //####---####---####---####---####---####---####---####---####---####---####---####---####---####---####---//
  if(TrggIndex_PFJet320_v3 < trigResults->size()) isFired_HLT_PFJet320 = trigResults->accept(TrggIndex_PFJet320_v3);
  if(TrggIndex_PFJet320_v4 < trigResults->size()) isFired_HLT_PFJet320 = trigResults->accept(TrggIndex_PFJet320_v4);
  if(TrggIndex_PFJet320_v5 < trigResults->size()) isFired_HLT_PFJet320 = trigResults->accept(TrggIndex_PFJet320_v5);
  if(TrggIndex_PFJet320_v6 < trigResults->size()) isFired_HLT_PFJet320 = trigResults->accept(TrggIndex_PFJet320_v6);
  if(TrggIndex_PFJet320_v8 < trigResults->size()) isFired_HLT_PFJet320 = trigResults->accept(TrggIndex_PFJet320_v8);
  if(TrggIndex_PFJet320_v9 < trigResults->size()) isFired_HLT_PFJet320 = trigResults->accept(TrggIndex_PFJet320_v9);
  //####---####---####---####---####---####---####---####---####---####---####---####---####---####---####---// 
  if(TrggIndex_PFJet40_v3 < trigResults->size()) isFired_HLT_PFJet40 = trigResults->accept(TrggIndex_PFJet40_v3);
  if(TrggIndex_PFJet40_v4 < trigResults->size()) isFired_HLT_PFJet40 = trigResults->accept(TrggIndex_PFJet40_v4);
  if(TrggIndex_PFJet40_v5 < trigResults->size()) isFired_HLT_PFJet40 = trigResults->accept(TrggIndex_PFJet40_v5);
  if(TrggIndex_PFJet40_v6 < trigResults->size()) isFired_HLT_PFJet40 = trigResults->accept(TrggIndex_PFJet40_v6);
  if(TrggIndex_PFJet40_v7 < trigResults->size()) isFired_HLT_PFJet40 = trigResults->accept(TrggIndex_PFJet40_v7);
  if(TrggIndex_PFJet40_v8 < trigResults->size()) isFired_HLT_PFJet40 = trigResults->accept(TrggIndex_PFJet40_v8); 
  //####---####---####---####---####---####---####---####---####---####---####---####---####---####---####---//
  if(TrggIndex_PFJet80_v3 < trigResults->size()) isFired_HLT_PFJet80 = trigResults->accept(TrggIndex_PFJet80_v3);
  if(TrggIndex_PFJet80_v4 < trigResults->size()) isFired_HLT_PFJet80 = trigResults->accept(TrggIndex_PFJet80_v4);
  if(TrggIndex_PFJet80_v5 < trigResults->size()) isFired_HLT_PFJet80 = trigResults->accept(TrggIndex_PFJet80_v5);
  if(TrggIndex_PFJet80_v6 < trigResults->size()) isFired_HLT_PFJet80 = trigResults->accept(TrggIndex_PFJet80_v6);
  if(TrggIndex_PFJet80_v8 < trigResults->size()) isFired_HLT_PFJet80 = trigResults->accept(TrggIndex_PFJet80_v8);
  if(TrggIndex_PFJet80_v9 < trigResults->size()) isFired_HLT_PFJet80 = trigResults->accept(TrggIndex_PFJet80_v9); 
  //####---####---####---####---####---####---####---####---####---####---####---####---####---####---####---//
  if(TrggIndex_PFJet140_v3 < trigResults->size()) isFired_HLT_PFJet140 = trigResults->accept(TrggIndex_PFJet140_v3);
  if(TrggIndex_PFJet140_v4 < trigResults->size()) isFired_HLT_PFJet140 = trigResults->accept(TrggIndex_PFJet140_v4);
  if(TrggIndex_PFJet140_v5 < trigResults->size()) isFired_HLT_PFJet140 = trigResults->accept(TrggIndex_PFJet140_v5);
  if(TrggIndex_PFJet140_v6 < trigResults->size()) isFired_HLT_PFJet140 = trigResults->accept(TrggIndex_PFJet140_v6);
  if(TrggIndex_PFJet140_v8 < trigResults->size()) isFired_HLT_PFJet140 = trigResults->accept(TrggIndex_PFJet140_v8);
  if(TrggIndex_PFJet140_v9 < trigResults->size()) isFired_HLT_PFJet140 = trigResults->accept(TrggIndex_PFJet140_v9); 
  //####---####---####---####---####---####---####---####---####---####---####---####---####---####---####---//
  if(TrggIndex_PFJet200_v3 < trigResults->size()) isFired_HLT_PFJet200 = trigResults->accept(TrggIndex_PFJet200_v3);
  if(TrggIndex_PFJet200_v4 < trigResults->size()) isFired_HLT_PFJet200 = trigResults->accept(TrggIndex_PFJet200_v4);
  if(TrggIndex_PFJet200_v5 < trigResults->size()) isFired_HLT_PFJet200 = trigResults->accept(TrggIndex_PFJet200_v5);
  if(TrggIndex_PFJet200_v6 < trigResults->size()) isFired_HLT_PFJet200 = trigResults->accept(TrggIndex_PFJet200_v6);
  if(TrggIndex_PFJet200_v8 < trigResults->size()) isFired_HLT_PFJet200 = trigResults->accept(TrggIndex_PFJet200_v8);
  if(TrggIndex_PFJet200_v9 < trigResults->size()) isFired_HLT_PFJet200 = trigResults->accept(TrggIndex_PFJet200_v9); 
  //####---####---####---####---####---####---####---####---####---####---####---####---####---####---####---//
  if(TrggIndex_PFJet260_v3 < trigResults->size()) isFired_HLT_PFJet260 = trigResults->accept(TrggIndex_PFJet260_v3);
  if(TrggIndex_PFJet260_v4 < trigResults->size()) isFired_HLT_PFJet260 = trigResults->accept(TrggIndex_PFJet260_v4);
  if(TrggIndex_PFJet260_v5 < trigResults->size()) isFired_HLT_PFJet260 = trigResults->accept(TrggIndex_PFJet260_v5);
  if(TrggIndex_PFJet260_v6 < trigResults->size()) isFired_HLT_PFJet260 = trigResults->accept(TrggIndex_PFJet260_v6);
  if(TrggIndex_PFJet260_v8 < trigResults->size()) isFired_HLT_PFJet260 = trigResults->accept(TrggIndex_PFJet260_v8);
  if(TrggIndex_PFJet260_v9 < trigResults->size()) isFired_HLT_PFJet260 = trigResults->accept(TrggIndex_PFJet260_v9);
  //####---####---####---####---####---####---####---####---####---####---####---####---####---####---####---// 
  
  bool foundJet=false;
  pat::JetCollection::const_iterator SelectedJet;
  float prunedMass=-99.; 
  float tau21Z=-99.; 
  float ptZ=-99.;
  for(pat::JetCollection::const_iterator jet = CA8JetswithQjets->begin(); jet != CA8JetswithQjets->end(); ++jet) {
    float dRmin = 9999.; float mass = 0.;
    for(pat::JetCollection::const_iterator jetPruned = CA8JetsPruned->begin(); jetPruned != CA8JetsPruned->end(); ++jetPruned) {
      float dRtmp = ROOT::Math::VectorUtil::DeltaR(jet->p4(),jetPruned->p4());
      if(dRtmp<dRmin && dRtmp<0.8 ){//matching failed if greater than jet radius
        dRmin=dRtmp;
        mass=jetPruned->mass();
      }
    }
    if(jet->muonEnergyFraction()>=0.99) continue;
    if(jet->photonEnergyFraction()>=0.99) continue;
    if(jet->chargedEmEnergyFraction()>=0.99) continue;
    if(jet->neutralHadronEnergyFraction()>=0.99) continue;
    if(jet->chargedHadronEnergyFraction()<=0.00) continue;
    //if(jet->pt()<400) continue;
    if(abs(jet->eta())>2.4) continue;
    if(!(mass>70 && mass<110))  continue;
    if(jet->userFloat("tau2")/jet->userFloat("tau1")>0.75) continue;
    foundJet=true;
    if(jet->pt()>ptZ){
      prunedMass=mass;
      ptZ=jet->pt();
      tau21Z=jet->userFloat("tau2")/jet->userFloat("tau1");
      SelectedJet=jet;
    }
  }
  
  bool foundJetNOTagged=false;
  pat::JetCollection::const_iterator SelectedJetNOTagged;
  float prunedMassNOTagged=-99.; 
  float tau21ZNOTagged=-99.; 
  float ptZNOTagged=-99.;
  for(pat::JetCollection::const_iterator jet = CA8JetswithQjets->begin(); jet != CA8JetswithQjets->end(); ++jet) {
    float dRmin = 9999.; float mass = 0.;
    for(pat::JetCollection::const_iterator jetPruned = CA8JetsPruned->begin(); jetPruned != CA8JetsPruned->end(); ++jetPruned) {
      float dRtmp = ROOT::Math::VectorUtil::DeltaR(jet->p4(),jetPruned->p4());
      if(dRtmp<dRmin && dRtmp<0.8 ){//matching failed if greater than jet radius
        dRmin=dRtmp;
        mass=jetPruned->mass();
      }
    }
    if(jet->muonEnergyFraction()>=0.99) continue;
    if(jet->photonEnergyFraction()>=0.99) continue;
    if(jet->chargedEmEnergyFraction()>=0.99) continue;
    if(jet->neutralHadronEnergyFraction()>=0.99) continue;
    if(jet->chargedHadronEnergyFraction()<=0.00) continue;
    //if(jet->pt()<400) continue;
    if(abs(jet->eta())>2.4) continue;
    //if(!(mass>massMin && mass<massMax))  continue;
    //if(jet->userFloat("tau2")/jet->userFloat("tau1")>0.75) continue;
    foundJetNOTagged=true;
    if(jet->pt()>ptZNOTagged){
      prunedMassNOTagged=mass;
      ptZNOTagged=jet->pt();
      tau21ZNOTagged=jet->userFloat("tau2")/jet->userFloat("tau1");
      SelectedJetNOTagged=jet;
    }
  }

  
  m_run=iEvent.id().run();
  m_lumi=iEvent.id().luminosityBlock();
  m_event=iEvent.id().event();
  m_met=met->begin()->pt();
  m_trigger40=(int)isFired_HLT_PFJet40;
  m_trigger80=(int)isFired_HLT_PFJet80;
  m_trigger140=(int)isFired_HLT_PFJet140;
  m_trigger200=(int)isFired_HLT_PFJet200;
  m_trigger260=(int)isFired_HLT_PFJet260;
  m_trigger320=(int)isFired_HLT_PFJet320;
  m_trigger650=(int)isFired_HLT_HT650;
  Tree->Fill();

  if(foundJet){
    m_runTagged=iEvent.id().run();
    m_lumiTagged=iEvent.id().luminosityBlock();
    m_eventTagged=iEvent.id().event();
    m_jetMassTagged=prunedMass;
    m_jetPtTagged=SelectedJet->pt();
    m_jetTau21Tagged=tau21Z;
    m_metTagged=met->begin()->pt();
    m_trigger40Tagged=(int)isFired_HLT_PFJet40;
    m_trigger80Tagged=(int)isFired_HLT_PFJet80;
    m_trigger140Tagged=(int)isFired_HLT_PFJet140;
    m_trigger200Tagged=(int)isFired_HLT_PFJet200;
    m_trigger260Tagged=(int)isFired_HLT_PFJet260;
    m_trigger320Tagged=(int)isFired_HLT_PFJet320;
    m_trigger650Tagged=(int)isFired_HLT_HT650;
    TreeTaggedJet->Fill();
  }
  
  if(foundJetNOTagged){
    m_runNO=iEvent.id().run();
    m_lumiNO=iEvent.id().luminosityBlock();
    m_eventNO=iEvent.id().event();
    m_jetMassNO=prunedMassNOTagged;
    m_jetPtNO=SelectedJetNOTagged->pt();
    m_jetTau21NO=tau21ZNOTagged;
    m_metNO=met->begin()->pt();
    m_trigger40NO=(int)isFired_HLT_PFJet40;
    m_trigger80NO=(int)isFired_HLT_PFJet80;
    m_trigger140NO=(int)isFired_HLT_PFJet140;
    m_trigger200NO=(int)isFired_HLT_PFJet200;
    m_trigger260NO=(int)isFired_HLT_PFJet260;
    m_trigger320NO=(int)isFired_HLT_PFJet320;
    m_trigger650NO=(int)isFired_HLT_HT650;
    TreeNOTaggedJet->Fill();
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
Trigger::beginJob()
{
  Service<TFileService> fs;
  Tree = fs->make<TTree>("Tree", "Tree");
  Tree->Branch("run", &m_run, "run/i");
  Tree->Branch("event", &m_event, "event/i");
  Tree->Branch("lumi", &m_lumi, "lumi/i");
  Tree->Branch("met", &m_met, "met/f");
  Tree->Branch("trigger40", &m_trigger40, "trigger40/i");
  Tree->Branch("trigger80", &m_trigger80, "trigger80/i");
  Tree->Branch("trigger140", &m_trigger140, "trigger140/i");
  Tree->Branch("trigger200", &m_trigger200, "trigger200/i");
  Tree->Branch("trigger260", &m_trigger260, "trigger260/i");
  Tree->Branch("trigger320", &m_trigger320, "trigger320/i");
  Tree->Branch("trigger650", &m_trigger650, "trigger650/i");

  TreeTaggedJet = fs->make<TTree>("TreeTaggedJet", "TreeTaggedJet");
  TreeTaggedJet->Branch("run", &m_runTagged, "run/i");
  TreeTaggedJet->Branch("event", &m_eventTagged, "event/i");
  TreeTaggedJet->Branch("lumi", &m_lumiTagged, "lumi/i");
  TreeTaggedJet->Branch("jetMass", &m_jetMassTagged, "jetMass/f");
  TreeTaggedJet->Branch("jetPt", &m_jetPtTagged, "jetPt/f");
  TreeTaggedJet->Branch("jetTau21", &m_jetTau21Tagged, "jetTau21/f");
  TreeTaggedJet->Branch("met", &m_metTagged, "met/f");
  TreeTaggedJet->Branch("trigger40", &m_trigger40Tagged, "trigger40/i");
  TreeTaggedJet->Branch("trigger80", &m_trigger80Tagged, "trigger80/i");
  TreeTaggedJet->Branch("trigger140", &m_trigger140Tagged, "trigger140/i");
  TreeTaggedJet->Branch("trigger200", &m_trigger200Tagged, "trigger200/i");
  TreeTaggedJet->Branch("trigger260", &m_trigger260Tagged, "trigger260/i");
  TreeTaggedJet->Branch("trigger320", &m_trigger320Tagged, "trigger320/i");
  TreeTaggedJet->Branch("trigger650", &m_trigger650Tagged, "trigger650/i");

  TreeNOTaggedJet = fs->make<TTree>("TreeNOTaggedJet", "TreeNOTaggedJet");
  TreeNOTaggedJet->Branch("run", &m_runNO, "run/i");
  TreeNOTaggedJet->Branch("event", &m_eventNO, "event/i");
  TreeNOTaggedJet->Branch("lumi", &m_lumiNO, "lumi/i");
  TreeNOTaggedJet->Branch("jetMass", &m_jetMassNO, "jetMass/f");
  TreeNOTaggedJet->Branch("jetPt", &m_jetPtNO, "jetPt/f");
  TreeNOTaggedJet->Branch("jetTau21", &m_jetTau21NO, "jetTau21/f");
  TreeNOTaggedJet->Branch("met", &m_metNO, "met/f");
  TreeNOTaggedJet->Branch("trigger40", &m_trigger40NO, "trigger40/i");
  TreeNOTaggedJet->Branch("trigger80", &m_trigger80NO, "trigger80/i");
  TreeNOTaggedJet->Branch("trigger140", &m_trigger140NO, "trigger140/i");
  TreeNOTaggedJet->Branch("trigger200", &m_trigger200NO, "trigger200/i");
  TreeNOTaggedJet->Branch("trigger260", &m_trigger260NO, "trigger260/i");
  TreeNOTaggedJet->Branch("trigger320", &m_trigger320NO, "trigger320/i");
  TreeNOTaggedJet->Branch("trigger650", &m_trigger650NO, "trigger650/i");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
Trigger::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
Trigger::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
Trigger::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Trigger::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Trigger::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Trigger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(Trigger);
