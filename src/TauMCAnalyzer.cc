// -*- C++ -*-
//
// Package:    TauMCAnalyzer
// Class:      TauMCAnalyzer
// 
/**\class TauMCAnalyzer TauMCAnalyzer.cc BulkG_TauTau/TauMCAnalyzer/src/TauMCAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Aniello Spiezia,21 1-007,+41227676459,
//         Created:  Mon Jul 15 17:08:06 CEST 2013
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "TauAnalysis/CandidateTools/interface/NSVfitStandaloneAlgorithm.h"

#include "TMath.h"
#include "TTree.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/VectorUtil.h"
#include "Math/GenVector/LorentzVector.h"
#include <vector>
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "TLorentzVector.h"



//
// class declaration
//

class TauMCAnalyzer : public edm::EDAnalyzer {
public:
  explicit TauMCAnalyzer(const edm::ParameterSet&);
  ~TauMCAnalyzer();
  
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
  float m_deltaRTauTau;

  TTree *TreeUsualMuoTau;
  TTree *TreeUsualEleTau;
  TTree *TreeMuoTau;
  TTree *TreeEleTau;
  TTree *TreeQCDUsual;
  TTree *TreeQCDMuoTau;
  TTree *TreeQCDEleTau;
  float m_MassSVFit;
  float m_ptSVFit;
  float m_MEtSigmaParlZ;
  float m_MEtSigmaPerpZ;
  float m_MEtPullParlZ;
  float m_MEtPullPerpZ;
  float m_metErrParl;
  float m_metErrPerp;
  float m_metGen;
  float m_met;
  float m_tauPtGen;
  float m_tauEtaGen;
  float m_tauMassGen;
  float m_tauMass;
  float m_deltaRReco;
  float m_deltaRGen;
  int m_trigger;
  int m_trigger650;
  int m_trigger320;
  float m_deltaR;
  float m_tauPt;
  float m_tauEta;
  float m_tauVertexZ;
  float m_tauDecay01;
  float m_tauDecay02;
  float m_tauDecay03;
  float m_tauIso01;
  float m_tauIso02;
  float m_tauIso03;
  float m_tauIso04;
  float m_tauIso05;
  float m_tauIso06;
  float m_tauIso07;
  float m_tauIso08;
  float m_tauIso09;
  float m_tauIso10;
  float m_tauIso11;
  float m_tauIso12;
  float m_tauIso13;
  float m_tauIso14;
  float m_tauIso15;
  float m_tauIso16;
  float m_tauIso17;
  float m_tauIso18;
  float m_tauIso19;
  float m_tauIso20;
  float m_tauIso21;
  float m_tauIso22;
  float m_tauIso23;
  float m_tauIso24;
  float m_tauIso25;
  float m_tauIso26;
  float m_tauIso27;
  float m_tauIso28;
  float m_tauIso29;
  float m_tauIso30;
  float m_tauIso31;
  float m_tauIso32;
  float m_tauIso33;
  float m_tauIso34;
  float m_tauIso35;
  float m_tauIso36;
  float m_tauIso37;
  float m_tauIso38;
  float m_tauIso39;
  float m_tauIso40;
  float m_tauAgainstElectron01;
  float m_tauAgainstElectron02;
  float m_tauAgainstElectron03;
  float m_tauAgainstElectron04;
  float m_tauAgainstElectron05;
  float m_tauAgainstElectron06;
  float m_tauAgainstElectron07;
  float m_tauAgainstElectron08;
  float m_tauAgainstElectron09;
  float m_tauAgainstElectron10;
  float m_tauAgainstElectron11;
  float m_tauAgainstMuon01;
  float m_tauAgainstMuon02;
  float m_tauAgainstMuon03;
  float m_tauAgainstMuon04;
  float m_tauAgainstMuon05;
  float m_tauAgainstMuon06;
  float m_tauAgainstMuon07;
  float m_tauAgainstMuon08;
  float m_tauAgainstMuon09;
  float m_tauAgainstMuon10;
  float m_tauAgainstMuon11;
  float m_tauAgainstMuon12;
  float m_tauPuCorrPtSum;

  TTree * TreeEffMuoTau;
  TTree * TreeEffEleTau;
  TTree * TreeEffUsualMuoTau;
  TTree * TreeEffUsualEleTau;
  int m_uno;  int m_unoU;
  int m_due;  int m_dueU;
  int m_tre;  int m_treU;
  int m_qua;  int m_quaU;
  int m_cin;  int m_cinU;
  int m_sei;  int m_seiU;
  int m_set;  int m_setU;
  
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
TauMCAnalyzer::TauMCAnalyzer(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
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

}


TauMCAnalyzer::~TauMCAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TauMCAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  Handle<vector<reco::GenParticle> > genParts;
  iEvent.getByLabel("genParticles", genParts);

  edm::Handle<pat::ElectronCollection> eleH;
  iEvent.getByLabel(electronColl_, eleH);

  edm::Handle<pat::MuonCollection> muoH;
  iEvent.getByLabel(muonColl_, muoH);

  Handle<pat::TauCollection> tauMuTauHandle;
  iEvent.getByLabel(tauMuTauColl_,tauMuTauHandle);

  Handle<pat::TauCollection> tauElTauHandle;
  iEvent.getByLabel(tauElTauColl_,tauElTauHandle);

  Handle<pat::TauCollection> tauHandle;
  iEvent.getByLabel("selectedPatTaus",tauHandle);
  
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel("offlinePrimaryVertices", vertices);
  
  edm::Handle<vector<reco::GenJet> > genjets;
  iEvent.getByLabel("ak5GenJetsNoNu", genjets);

  //Handle<reco::GenMETCollection> genMets;
  //iEvent.getByLabel("genMetTrue",genMets);

  edm::Handle<pat::METCollection> met;
  iEvent.getByLabel(metColl_, met);

  edm::Handle<pat::METCollection> metRaw;
  iEvent.getByLabel(metRawColl_, metRaw);

  edm::Handle<pat::JetCollection> ak5jetCands;
  iEvent.getByLabel(ak5JetColl_,ak5jetCands);

  bool tauele1 = false; bool taumuo1 = false;
  bool tauele2 = false; bool taumuo2 = false;
  vector<math::PtEtaPhiELorentzVector> gentau;
  vector<math::PtEtaPhiELorentzVector> gentauHad;
  math::PtEtaPhiELorentzVector genMuo;   bool foundMuo   = false;
  math::PtEtaPhiELorentzVector genEle;   bool foundEle   = false;
  math::PtEtaPhiELorentzVector genHiggs; bool foundHiggs = false;
  for(size_t ngenPart=0; ngenPart<genParts->size(); ngenPart++){
    const reco::GenParticle & genPart = (*genParts)[ngenPart];
    if(abs(genPart.pdgId())==2212) continue;
    //if((!genPart.pt()>20)) continue;
    const reco::Candidate * mom = genPart.mother();
    int tauCharge=0;
    if(abs(genPart.pdgId())==15 && genPart.status()!=3 && abs(mom->pdgId())==15 && genPart.pt()>20){
      math::PtEtaPhiELorentzVector gentau_prov;
      math::PtEtaPhiELorentzVector gentauHad_prov;
      bool genTauHadBool = false;
      gentau_prov=genPart.p4();
      for(unsigned int ndaugh=0; ndaugh<genPart.numberOfDaughters(); ndaugh++){
	const reco::Candidate * daughter = genPart.daughter(ndaugh);
	if(genPart.pdgId()== 15 && abs(daughter->pdgId())==11) tauele1 = true;
	if(genPart.pdgId()== 15 && abs(daughter->pdgId())==13) taumuo1 = true;
	if(genPart.pdgId()==-15 && abs(daughter->pdgId())==11) tauele2 = true;
	if(genPart.pdgId()==-15 && abs(daughter->pdgId())==13) taumuo2 = true;
	if(abs(daughter->pdgId())!=11 && abs(daughter->pdgId())!=12 && abs(daughter->pdgId())!=13 && abs(daughter->pdgId())!=14 
	   && abs(daughter->pdgId())!=15 && abs(daughter->pdgId())!=16 && abs(daughter->pdgId())!=22){
	  gentauHad_prov=genPart.p4();
	  tauCharge=genPart.pdgId();
	  genTauHadBool = true;
	}
      }
      for(unsigned int ndaugh=0; ndaugh<genPart.numberOfDaughters(); ndaugh++){
	const reco::Candidate * daughter = genPart.daughter(ndaugh);
	if((abs(daughter->pdgId())==12 || abs(daughter->pdgId())==14 ||abs(daughter->pdgId())==16) && daughter->status()==1){
	  gentau_prov    = gentau_prov    - daughter->p4();
	  if(genTauHadBool && genPart.pdgId()==tauCharge) gentauHad_prov = gentauHad_prov - daughter->p4();
	}
      }
      gentau.push_back(gentau_prov);
      if(genTauHadBool) gentauHad.push_back(gentauHad_prov);
    }

    //LOOK FOR THE ELECTRON
    if(abs(genPart.pdgId())==11 && genPart.status()==1 && abs(mom->pdgId())==15 && genPart.pt()>10){
      genEle = genPart.p4();
      foundEle = true;
    }

    //LOOK FOR THE MUON
    if(abs(genPart.pdgId())==13 && genPart.status()==1 && abs(mom->pdgId())==15 && genPart.pt()>10){
      genMuo = genPart.p4();
      foundMuo = true;
    }

    //LOOK FOR THE HIGGS
    if(abs(genPart.pdgId())==25 && genPart.status()!=3){
      genHiggs = genPart.p4();
      foundHiggs = true;
    }
  }

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

  int uno = 0;  int unoU = 0;
  int due = 0;	int dueU = 0;
  int tre = 0;	int treU = 0;
  int qua = 0;	int quaU = 0;
  int cin = 0;	int cinU = 0;
  int sei = 0;	int seiU = 0;
  int set = 0;	int setU = 0;

  //GENERAL PLOTS
  if(!isData && gentau.size()==2){
    m_deltaRTauTau = ROOT::Math::VectorUtil::DeltaR(gentau[0],gentau[1]);
    Tree->Fill();
  }

  //MUO-TAU CHANNEL
  if(!isData && gentau.size()==2 && ((taumuo1 && !taumuo2 && !tauele2) || (!taumuo1 && !tauele1 && taumuo2))){

    float deltaRCleaned=99.; bool matchedCleaned = false;
    pat::TauCollection::const_iterator SelectedTauMuTau;
    for (pat::TauCollection::const_iterator patTau = tauMuTauHandle->begin(); patTau != tauMuTauHandle->end(); ++patTau ) {
      //if(!(patTau->pt()>20)) continue;
      if(ROOT::Math::VectorUtil::DeltaR(gentauHad[0],patTau->p4())<0.5 && ROOT::Math::VectorUtil::DeltaR(gentauHad[0],patTau->p4())<deltaRCleaned){
	deltaRCleaned = ROOT::Math::VectorUtil::DeltaR(gentauHad[0],patTau->p4());
	SelectedTauMuTau=patTau;
	matchedCleaned = true;
      }
    }

    float deltaRUsual=99.; bool matchedUsual = false;
    pat::TauCollection::const_iterator SelectedTau;
    for (pat::TauCollection::const_iterator patTau = tauHandle->begin(); patTau != tauHandle->end(); ++patTau ) {
      //if(!(patTau->pt()>20)) continue;
      if(ROOT::Math::VectorUtil::DeltaR(gentauHad[0],patTau->p4())<0.5 && ROOT::Math::VectorUtil::DeltaR(gentauHad[0],patTau->p4())<deltaRUsual){
	deltaRUsual = ROOT::Math::VectorUtil::DeltaR(gentauHad[0],patTau->p4());
	SelectedTau=patTau;
	matchedUsual = true;
      }
    }

    float deltaRMuo=99.; bool matchedMuo = false;
    pat::MuonCollection::const_iterator recoMuo;
    for(pat::MuonCollection::const_iterator muon = muoH->begin(); muon != muoH->end(); ++muon) {
      if(!(muon->pt()>10)) continue;
      if(ROOT::Math::VectorUtil::DeltaR(genMuo,muon->p4())<0.3 && ROOT::Math::VectorUtil::DeltaR(genMuo,muon->p4())<deltaRMuo){
	deltaRMuo = ROOT::Math::VectorUtil::DeltaR(genMuo,muon->p4());
	recoMuo=muon;
	matchedMuo = true;
      }
    }

    uno = 1;
    if(matchedCleaned) due=1;
    if(matchedCleaned && SelectedTauMuTau->pt()>20) tre=1;
    if(matchedCleaned && SelectedTauMuTau->pt()>20 && SelectedTauMuTau->tauID("decayModeFindingNewDMs")>0.5) qua=1;
    if(matchedCleaned && SelectedTauMuTau->pt()>20 && SelectedTauMuTau->tauID("decayModeFindingNewDMs")>0.5 && SelectedTauMuTau->tauID("againstMuonLoose")>0.5) cin=1;
    if(matchedCleaned && SelectedTauMuTau->pt()>20 && SelectedTauMuTau->tauID("decayModeFindingNewDMs")>0.5 && SelectedTauMuTau->tauID("againstMuonLoose")>0.5 && 
       SelectedTauMuTau->tauID("againstElectronLoose")>0.5) sei=1;
    if(matchedCleaned && SelectedTauMuTau->pt()>20 && SelectedTauMuTau->tauID("decayModeFindingNewDMs")>0.5 && SelectedTauMuTau->tauID("againstMuonLoose")>0.5 && 
       SelectedTauMuTau->tauID("againstElectronLoose")>0.5 && SelectedTauMuTau->tauID("byVLooseIsolationMVA3newDMwoLT")>0.5) set=1;
    m_uno = uno;
    m_due = due;
    m_tre = tre;
    m_qua = qua;
    m_cin = cin;
    m_sei = sei;
    m_set = set;
    TreeEffMuoTau->Fill();
    unoU = 1;
    if(matchedUsual) dueU=1;
    if(matchedUsual && SelectedTau->pt()>20) treU=1;
    if(matchedUsual && SelectedTau->pt()>20 && SelectedTau->tauID("decayModeFindingNewDMs")>0.5) quaU=1;
    if(matchedUsual && SelectedTau->pt()>20 && SelectedTau->tauID("decayModeFindingNewDMs")>0.5 && SelectedTau->tauID("againstMuonLoose")>0.5) cinU=1;
    if(matchedUsual && SelectedTau->pt()>20 && SelectedTau->tauID("decayModeFindingNewDMs")>0.5 && SelectedTau->tauID("againstMuonLoose")>0.5 && 
       SelectedTau->tauID("againstElectronLoose")>0.5) seiU=1;
    if(matchedUsual && SelectedTau->pt()>20 && SelectedTau->tauID("decayModeFindingNewDMs")>0.5 && SelectedTau->tauID("againstMuonLoose")>0.5 && 
       SelectedTau->tauID("againstElectronLoose")>0.5 && SelectedTau->tauID("byVLooseIsolationMVA3newDMwoLT")>0.5) setU=1;
    m_unoU = unoU;
    m_dueU = dueU;
    m_treU = treU;
    m_quaU = quaU;
    m_cinU = cinU;
    m_seiU = seiU;
    m_setU = setU;
    TreeEffUsualMuoTau->Fill();

    if(matchedCleaned){
      float dRmin=999.; pat::JetCollection::const_iterator TauJet; bool TauMatched = false;
      for(pat::JetCollection::const_iterator ak5 = ak5jetCands->begin(); ak5 != ak5jetCands->end(); ++ak5) {
	if(!(ak5->pt()>15)) continue;
	if(ROOT::Math::VectorUtil::DeltaR(SelectedTauMuTau->p4Jet(),ak5->p4())<0.3 && ROOT::Math::VectorUtil::DeltaR(SelectedTauMuTau->p4Jet(),ak5->p4())<dRmin){
	  dRmin=ROOT::Math::VectorUtil::DeltaR(SelectedTauMuTau->p4Jet(),ak5->p4());
	  TauJet=ak5;
	  TauMatched=true;
	}
      }

      float deltaRReco=99.;
      if(matchedMuo) deltaRReco=ROOT::Math::VectorUtil::DeltaR(SelectedTauMuTau->p4(),recoMuo->p4());

      /*
      //corrPFMET = rawPFMET + sum(rawJetPt-corrJetPt) [*]
      float PX = 0.;
      float PY = 0.;
      //pat::MET *CorrectedMet;
      //TLorentzVector CorrectedMet_prov; 
      if(TauMatched){
	PX = met->begin()->px() + TauJet->px()*(1 - TauJet->jecFactor(0));
	PY = met->begin()->py() + TauJet->py()*(1 - TauJet->jecFactor(0));
      } else {
	PX = met->begin()->px();
	PY = met->begin()->py();
      }
      reco::Candidate::Vector CorrectedMet(PX,PY,0.);
      reco::Candidate::LorentzVector CorrectedMetP4(PX,PY,0.,sqrt(PX*PX+PY*PY));

      reco::Candidate::LorentzVector dMET = CorrectedMetP4 - genMets->begin()->p4();
      double cosHiggsPhi = TMath::Cos(genHiggs.phi());
      double sinHiggsPhi = TMath::Sin(genHiggs.phi());
      double px = dMET.px();
      double py = dMET.py();
      double pxInZetaFrame =  cosHiggsPhi*px + sinHiggsPhi*py;
      double pyInZetaFrame = -sinHiggsPhi*px + cosHiggsPhi*py;
      reco::Candidate::LorentzVector metErr_rotated(pxInZetaFrame, pyInZetaFrame, dMET.pz(), dMET.E());
      double metErrParl = metErr_rotated.px();
      double metErrPerp = metErr_rotated.py();

      TMatrixD metCov(2, 2); // PFMET significance matrix
      metCov[0][0] = (metRaw->front() ).getSignificanceMatrix()(0,0);
      metCov[1][0] = (metRaw->front() ).getSignificanceMatrix()(1,0);
      metCov[0][1] = (metRaw->front() ).getSignificanceMatrix()(0,1);
      metCov[1][1] = (metRaw->front() ).getSignificanceMatrix()(1,1);
      double cov00 = metCov(0,0);
      double cov01 = metCov(0,1);
      double cov10 = metCov(1,0);
      double cov11 = metCov(1,1);
      TMatrixD metCov_rotated(2,2);
      metCov_rotated(0,0) =  cosHiggsPhi*cosHiggsPhi*cov00 + cosHiggsPhi*sinHiggsPhi*cov01 + sinHiggsPhi*cosHiggsPhi*cov10 + sinHiggsPhi*sinHiggsPhi*cov11;
      metCov_rotated(0,1) = -sinHiggsPhi*cosHiggsPhi*cov00 + cosHiggsPhi*cosHiggsPhi*cov01 - sinHiggsPhi*sinHiggsPhi*cov10 + cosHiggsPhi*sinHiggsPhi*cov11;
      metCov_rotated(1,0) = -cosHiggsPhi*sinHiggsPhi*cov00 - sinHiggsPhi*sinHiggsPhi*cov01 + cosHiggsPhi*cosHiggsPhi*cov10 + sinHiggsPhi*cosHiggsPhi*cov11;
      metCov_rotated(1,1) =  sinHiggsPhi*sinHiggsPhi*cov00 - cosHiggsPhi*sinHiggsPhi*cov01 - sinHiggsPhi*cosHiggsPhi*cov10 + cosHiggsPhi*cosHiggsPhi*cov11;
      double metSigmaParl = TMath::Sqrt(TMath::Abs(metCov_rotated(0, 0)));
      double metSigmaPerp = TMath::Sqrt(TMath::Abs(metCov_rotated(1, 1)));
      float MEtPullParlZ=-9999.;
      float MEtPullPerpZ=-9999.;
      if(metSigmaParl > 0.) MEtPullParlZ=metErrParl/metSigmaParl;  
      if(metSigmaPerp > 0.) MEtPullPerpZ=metErrPerp/metSigmaPerp;

      //SVFIT
      float MassSVFit = -1;
      float ptSVFit = -1;
      if(matchedMuo){
	std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
	measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay, SelectedTauMuTau->p4()));
	measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, recoMuo->p4()));
	NSVfitStandaloneAlgorithm algo(measuredTauLeptons, CorrectedMet, metCov, 0);
	algo.addLogM(false);
	algo.integrateMarkovChain();
	if(algo.pt()>0){
	  MassSVFit=algo.getMass();
	  ptSVFit=algo.pt();
	}
      }

      m_MassSVFit=MassSVFit;
      m_ptSVFit=ptSVFit;
      m_MEtSigmaParlZ = metSigmaParl;
      m_MEtSigmaPerpZ = metSigmaPerp;
      m_MEtPullParlZ  = MEtPullParlZ;
      m_MEtPullPerpZ  = MEtPullPerpZ;
      m_metErrParl    = metErrParl;
      m_metErrPerp    = metErrPerp;
      m_metGen        = genMets->begin()->pt();
      m_met           = sqrt(PX*PX+PY*PY);
      */

      m_tauPtGen      = gentauHad[0].pt();
      m_tauEtaGen     = gentauHad[0].eta();
      m_tauMassGen    = gentauHad[0].mass();
      m_tauMass       = SelectedTauMuTau->mass();
      m_deltaRReco    = deltaRReco;
      m_deltaRGen     = ROOT::Math::VectorUtil::DeltaR(gentauHad[0],genMuo);
      m_tauPt=SelectedTauMuTau->pt();
      m_tauEta=SelectedTauMuTau->eta();
      m_deltaR=deltaRCleaned;
      m_trigger=(int)isFired_HLT;
      m_trigger320=(int)isFired_HLT_PFJet320;
      m_trigger650=(int)isFired_HLT_HT650;
      m_tauDecay01 = SelectedTauMuTau->tauID("decayModeFinding");
      m_tauDecay02 = SelectedTauMuTau->tauID("decayModeFindingNewDMs");
      m_tauDecay03 = SelectedTauMuTau->tauID("decayModeFindingOldDMs");
      m_tauIso01   = SelectedTauMuTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr");
      m_tauIso02   = SelectedTauMuTau->tauID("byLooseIsolation");
      m_tauIso03   = SelectedTauMuTau->tauID("byLooseCombinedIsolationDeltaBetaCorr");
      m_tauIso04   = SelectedTauMuTau->tauID("byMediumCombinedIsolationDeltaBetaCorr");
      m_tauIso05   = SelectedTauMuTau->tauID("byTightCombinedIsolationDeltaBetaCorr");
      m_tauIso06   = SelectedTauMuTau->tauID("byCombinedIsolationDeltaBetaCorrRaw");
      m_tauIso07   = SelectedTauMuTau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
      m_tauIso08   = SelectedTauMuTau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
      m_tauIso09   = SelectedTauMuTau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
      m_tauIso10   = SelectedTauMuTau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
      m_tauIso11   = SelectedTauMuTau->tauID("chargedIsoPtSum");
      m_tauIso12   = SelectedTauMuTau->tauID("neutralIsoPtSum");
      m_tauIso13   = SelectedTauMuTau->tauID("byIsolationMVA3oldDMwoLTraw");
      m_tauIso14   = SelectedTauMuTau->tauID("byVLooseIsolationMVA3oldDMwoLT");
      m_tauIso15   = SelectedTauMuTau->tauID("byLooseIsolationMVA3oldDMwoLT");
      m_tauIso16   = SelectedTauMuTau->tauID("byMediumIsolationMVA3oldDMwoLT");
      m_tauIso17   = SelectedTauMuTau->tauID("byTightIsolationMVA3oldDMwoLT");
      m_tauIso18   = SelectedTauMuTau->tauID("byVTightIsolationMVA3oldDMwoLT");
      m_tauIso19   = SelectedTauMuTau->tauID("byVVTightIsolationMVA3oldDMwoLT");
      m_tauIso20   = SelectedTauMuTau->tauID("byIsolationMVA3oldDMwLTraw");
      m_tauIso21   = SelectedTauMuTau->tauID("byVLooseIsolationMVA3oldDMwLT");
      m_tauIso22   = SelectedTauMuTau->tauID("byLooseIsolationMVA3oldDMwLT");
      m_tauIso23   = SelectedTauMuTau->tauID("byMediumIsolationMVA3oldDMwLT");
      m_tauIso24   = SelectedTauMuTau->tauID("byTightIsolationMVA3oldDMwLT");
      m_tauIso25   = SelectedTauMuTau->tauID("byVTightIsolationMVA3oldDMwLT");
      m_tauIso26   = SelectedTauMuTau->tauID("byVVTightIsolationMVA3oldDMwLT");
      m_tauIso27   = SelectedTauMuTau->tauID("byIsolationMVA3newDMwoLTraw");
      m_tauIso28   = SelectedTauMuTau->tauID("byVLooseIsolationMVA3newDMwoLT");
      m_tauIso29   = SelectedTauMuTau->tauID("byLooseIsolationMVA3newDMwoLT");
      m_tauIso30   = SelectedTauMuTau->tauID("byMediumIsolationMVA3newDMwoLT");
      m_tauIso31   = SelectedTauMuTau->tauID("byTightIsolationMVA3newDMwoLT");
      m_tauIso32   = SelectedTauMuTau->tauID("byVTightIsolationMVA3newDMwoLT");
      m_tauIso33   = SelectedTauMuTau->tauID("byVVTightIsolationMVA3newDMwoLT");
      m_tauIso34   = SelectedTauMuTau->tauID("byIsolationMVA3newDMwLTraw");
      m_tauIso35   = SelectedTauMuTau->tauID("byVLooseIsolationMVA3newDMwLT");
      m_tauIso36   = SelectedTauMuTau->tauID("byLooseIsolationMVA3newDMwLT");
      m_tauIso37   = SelectedTauMuTau->tauID("byMediumIsolationMVA3newDMwLT");
      m_tauIso38   = SelectedTauMuTau->tauID("byTightIsolationMVA3newDMwLT");
      m_tauIso39   = SelectedTauMuTau->tauID("byVTightIsolationMVA3newDMwLT");
      m_tauIso40   = SelectedTauMuTau->tauID("byVVTightIsolationMVA3newDMwLT");
      m_tauAgainstElectron01 = SelectedTauMuTau->tauID("againstElectronLoose");
      m_tauAgainstElectron02 = SelectedTauMuTau->tauID("againstElectronMedium");
      m_tauAgainstElectron03 = SelectedTauMuTau->tauID("againstElectronTight");
      m_tauAgainstElectron04 = SelectedTauMuTau->tauID("againstElectronMVA5raw");
      m_tauAgainstElectron05 = SelectedTauMuTau->tauID("againstElectronMVA5category");
      m_tauAgainstElectron06 = SelectedTauMuTau->tauID("againstElectronVLooseMVA5");
      m_tauAgainstElectron07 = SelectedTauMuTau->tauID("againstElectronLooseMVA5");
      m_tauAgainstElectron08 = SelectedTauMuTau->tauID("againstElectronMediumMVA5");
      m_tauAgainstElectron09 = SelectedTauMuTau->tauID("againstElectronTightMVA5");
      m_tauAgainstElectron10 = SelectedTauMuTau->tauID("againstElectronVTightMVA5");
      m_tauAgainstElectron11 = SelectedTauMuTau->tauID("againstElectronDeadECAL");
      m_tauAgainstMuon01 = SelectedTauMuTau->tauID("againstMuonLoose");
      m_tauAgainstMuon02 = SelectedTauMuTau->tauID("againstMuonMedium");
      m_tauAgainstMuon03 = SelectedTauMuTau->tauID("againstMuonTight");
      m_tauAgainstMuon04 = SelectedTauMuTau->tauID("againstMuonLoose2");
      m_tauAgainstMuon05 = SelectedTauMuTau->tauID("againstMuonMedium2");
      m_tauAgainstMuon06 = SelectedTauMuTau->tauID("againstMuonTight2");
      m_tauAgainstMuon07 = SelectedTauMuTau->tauID("againstMuonLoose3");
      m_tauAgainstMuon08 = SelectedTauMuTau->tauID("againstMuonTight3");
      m_tauAgainstMuon09 = SelectedTauMuTau->tauID("againstMuonMVAraw");
      m_tauAgainstMuon10 = SelectedTauMuTau->tauID("againstMuonLooseMVA");
      m_tauAgainstMuon11 = SelectedTauMuTau->tauID("againstMuonMediumMVA");
      m_tauAgainstMuon12 = SelectedTauMuTau->tauID("againstMuonTightMVA");
      m_tauPuCorrPtSum   = SelectedTauMuTau->tauID("puCorrPtSum");
      TreeMuoTau->Fill();
    }

    if(matchedUsual){

      // "met->begin()->"   ---->   "CorrectedMet->"
      float dRmin=999.; pat::JetCollection::const_iterator TauJet; bool TauMatched = false;
      for(pat::JetCollection::const_iterator ak5 = ak5jetCands->begin(); ak5 != ak5jetCands->end(); ++ak5) {
	//if(!(ak5->pt()>20)) continue;
	if(ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4Jet(),ak5->p4())<0.3 && ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4Jet(),ak5->p4())<dRmin){
	  dRmin=ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4Jet(),ak5->p4());
	  TauJet=ak5;
	  TauMatched=true;
	}
      }

      float deltaRReco=99.;
      if(matchedMuo) deltaRReco=ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(),recoMuo->p4());

      /*
      //corrPFMET = rawPFMET + sum(rawJetPt-corrJetPt) [*]
      float PX = 0.;
      float PY = 0.;
      //pat::MET *CorrectedMet;
      //TLorentzVector CorrectedMet_prov; 
      if(TauMatched){
	PX = met->begin()->px() + TauJet->px()*(1 - TauJet->jecFactor(0));
	PY = met->begin()->py() + TauJet->py()*(1 - TauJet->jecFactor(0));
      } else {
	PX = met->begin()->px();
	PY = met->begin()->py();
      }
      reco::Candidate::Vector CorrectedMet(PX,PY,0.);
      reco::Candidate::LorentzVector CorrectedMetP4(PX,PY,0.,sqrt(PX*PX+PY*PY));

      reco::Candidate::LorentzVector dMET = CorrectedMetP4 - genMets->begin()->p4();
      double cosHiggsPhi = TMath::Cos(genHiggs.phi());
      double sinHiggsPhi = TMath::Sin(genHiggs.phi());
      double px = dMET.px();
      double py = dMET.py();
      double pxInZetaFrame =  cosHiggsPhi*px + sinHiggsPhi*py;
      double pyInZetaFrame = -sinHiggsPhi*px + cosHiggsPhi*py;
      reco::Candidate::LorentzVector metErr_rotated(pxInZetaFrame, pyInZetaFrame, dMET.pz(), dMET.E());
      double metErrParl = metErr_rotated.px();
      double metErrPerp = metErr_rotated.py();

      TMatrixD metCov(2, 2); // PFMET significance matrix
      metCov[0][0] = (metRaw->front() ).getSignificanceMatrix()(0,0);
      metCov[1][0] = (metRaw->front() ).getSignificanceMatrix()(1,0);
      metCov[0][1] = (metRaw->front() ).getSignificanceMatrix()(0,1);
      metCov[1][1] = (metRaw->front() ).getSignificanceMatrix()(1,1);
      double cov00 = metCov(0,0);
      double cov01 = metCov(0,1);
      double cov10 = metCov(1,0);
      double cov11 = metCov(1,1);
      TMatrixD metCov_rotated(2,2);
      metCov_rotated(0,0) =  cosHiggsPhi*cosHiggsPhi*cov00 + cosHiggsPhi*sinHiggsPhi*cov01 + sinHiggsPhi*cosHiggsPhi*cov10 + sinHiggsPhi*sinHiggsPhi*cov11;
      metCov_rotated(0,1) = -sinHiggsPhi*cosHiggsPhi*cov00 + cosHiggsPhi*cosHiggsPhi*cov01 - sinHiggsPhi*sinHiggsPhi*cov10 + cosHiggsPhi*sinHiggsPhi*cov11;
      metCov_rotated(1,0) = -cosHiggsPhi*sinHiggsPhi*cov00 - sinHiggsPhi*sinHiggsPhi*cov01 + cosHiggsPhi*cosHiggsPhi*cov10 + sinHiggsPhi*cosHiggsPhi*cov11;
      metCov_rotated(1,1) =  sinHiggsPhi*sinHiggsPhi*cov00 - cosHiggsPhi*sinHiggsPhi*cov01 - sinHiggsPhi*cosHiggsPhi*cov10 + cosHiggsPhi*cosHiggsPhi*cov11;
      double metSigmaParl = TMath::Sqrt(TMath::Abs(metCov_rotated(0, 0)));
      double metSigmaPerp = TMath::Sqrt(TMath::Abs(metCov_rotated(1, 1)));
      float MEtPullParlZ=-9999.;
      float MEtPullPerpZ=-9999.;
      if(metSigmaParl > 0.) MEtPullParlZ=metErrParl/metSigmaParl;  
      if(metSigmaPerp > 0.) MEtPullPerpZ=metErrPerp/metSigmaPerp;

      //SVFIT
      float MassSVFit = -1;
      float ptSVFit = -1;
      if(matchedMuo){
	std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
	measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay, SelectedTau->p4()));
	measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, recoMuo->p4()));
	NSVfitStandaloneAlgorithm algo(measuredTauLeptons, CorrectedMet, metCov, 0);
	algo.addLogM(false);
	algo.integrateMarkovChain();
	if(algo.pt()>0){
	  MassSVFit=algo.getMass();
	  ptSVFit=algo.pt();
	}
      }
      
      m_MassSVFit=MassSVFit;
      m_ptSVFit=ptSVFit;
      m_MEtSigmaParlZ = metSigmaParl;
      m_MEtSigmaPerpZ = metSigmaPerp;
      m_MEtPullParlZ  = MEtPullParlZ;
      m_MEtPullPerpZ  = MEtPullPerpZ;
      m_metErrParl    = metErrParl;
      m_metErrPerp    = metErrPerp;
      m_metGen        = genMets->begin()->pt();
      m_met           = sqrt(PX*PX+PY*PY);
      */

      m_tauPtGen      = gentauHad[0].pt();
      m_tauEtaGen     = gentauHad[0].eta();
      m_tauMassGen    = gentauHad[0].mass();
      m_tauMass       = SelectedTau->mass();
      m_deltaRReco    = deltaRReco;
      m_deltaRGen     = ROOT::Math::VectorUtil::DeltaR(gentauHad[0],genMuo);
      m_tauPt=SelectedTau->pt();
      m_tauEta=SelectedTau->eta();
      m_deltaR=deltaRUsual;
      m_trigger=(int)isFired_HLT;
      m_trigger320=(int)isFired_HLT_PFJet320;
      m_trigger650=(int)isFired_HLT_HT650;
      m_tauDecay01 = SelectedTau->tauID("decayModeFinding");
      m_tauDecay02 = SelectedTau->tauID("decayModeFindingNewDMs");
      m_tauDecay03 = SelectedTau->tauID("decayModeFindingOldDMs");
      m_tauIso01   = SelectedTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr");
      m_tauIso02   = SelectedTau->tauID("byLooseIsolation");
      m_tauIso03   = SelectedTau->tauID("byLooseCombinedIsolationDeltaBetaCorr");
      m_tauIso04   = SelectedTau->tauID("byMediumCombinedIsolationDeltaBetaCorr");
      m_tauIso05   = SelectedTau->tauID("byTightCombinedIsolationDeltaBetaCorr");
      m_tauIso06   = SelectedTau->tauID("byCombinedIsolationDeltaBetaCorrRaw");
      m_tauIso07   = SelectedTau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
      m_tauIso08   = SelectedTau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
      m_tauIso09   = SelectedTau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
      m_tauIso10   = SelectedTau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
      m_tauIso11   = SelectedTau->tauID("chargedIsoPtSum");
      m_tauIso12   = SelectedTau->tauID("neutralIsoPtSum");
      m_tauIso13   = SelectedTau->tauID("byIsolationMVA3oldDMwoLTraw");
      m_tauIso14   = SelectedTau->tauID("byVLooseIsolationMVA3oldDMwoLT");
      m_tauIso15   = SelectedTau->tauID("byLooseIsolationMVA3oldDMwoLT");
      m_tauIso16   = SelectedTau->tauID("byMediumIsolationMVA3oldDMwoLT");
      m_tauIso17   = SelectedTau->tauID("byTightIsolationMVA3oldDMwoLT");
      m_tauIso18   = SelectedTau->tauID("byVTightIsolationMVA3oldDMwoLT");
      m_tauIso19   = SelectedTau->tauID("byVVTightIsolationMVA3oldDMwoLT");
      m_tauIso20   = SelectedTau->tauID("byIsolationMVA3oldDMwLTraw");
      m_tauIso21   = SelectedTau->tauID("byVLooseIsolationMVA3oldDMwLT");
      m_tauIso22   = SelectedTau->tauID("byLooseIsolationMVA3oldDMwLT");
      m_tauIso23   = SelectedTau->tauID("byMediumIsolationMVA3oldDMwLT");
      m_tauIso24   = SelectedTau->tauID("byTightIsolationMVA3oldDMwLT");
      m_tauIso25   = SelectedTau->tauID("byVTightIsolationMVA3oldDMwLT");
      m_tauIso26   = SelectedTau->tauID("byVVTightIsolationMVA3oldDMwLT");
      m_tauIso27   = SelectedTau->tauID("byIsolationMVA3newDMwoLTraw");
      m_tauIso28   = SelectedTau->tauID("byVLooseIsolationMVA3newDMwoLT");
      m_tauIso29   = SelectedTau->tauID("byLooseIsolationMVA3newDMwoLT");
      m_tauIso30   = SelectedTau->tauID("byMediumIsolationMVA3newDMwoLT");
      m_tauIso31   = SelectedTau->tauID("byTightIsolationMVA3newDMwoLT");
      m_tauIso32   = SelectedTau->tauID("byVTightIsolationMVA3newDMwoLT");
      m_tauIso33   = SelectedTau->tauID("byVVTightIsolationMVA3newDMwoLT");
      m_tauIso34   = SelectedTau->tauID("byIsolationMVA3newDMwLTraw");
      m_tauIso35   = SelectedTau->tauID("byVLooseIsolationMVA3newDMwLT");
      m_tauIso36   = SelectedTau->tauID("byLooseIsolationMVA3newDMwLT");
      m_tauIso37   = SelectedTau->tauID("byMediumIsolationMVA3newDMwLT");
      m_tauIso38   = SelectedTau->tauID("byTightIsolationMVA3newDMwLT");
      m_tauIso39   = SelectedTau->tauID("byVTightIsolationMVA3newDMwLT");
      m_tauIso40   = SelectedTau->tauID("byVVTightIsolationMVA3newDMwLT");
      m_tauAgainstElectron01 = SelectedTau->tauID("againstElectronLoose");
      m_tauAgainstElectron02 = SelectedTau->tauID("againstElectronMedium");
      m_tauAgainstElectron03 = SelectedTau->tauID("againstElectronTight");
      m_tauAgainstElectron04 = SelectedTau->tauID("againstElectronMVA5raw");
      m_tauAgainstElectron05 = SelectedTau->tauID("againstElectronMVA5category");
      m_tauAgainstElectron06 = SelectedTau->tauID("againstElectronVLooseMVA5");
      m_tauAgainstElectron07 = SelectedTau->tauID("againstElectronLooseMVA5");
      m_tauAgainstElectron08 = SelectedTau->tauID("againstElectronMediumMVA5");
      m_tauAgainstElectron09 = SelectedTau->tauID("againstElectronTightMVA5");
      m_tauAgainstElectron10 = SelectedTau->tauID("againstElectronVTightMVA5");
      m_tauAgainstElectron11 = SelectedTau->tauID("againstElectronDeadECAL");
      m_tauAgainstMuon01 = SelectedTau->tauID("againstMuonLoose");
      m_tauAgainstMuon02 = SelectedTau->tauID("againstMuonMedium");
      m_tauAgainstMuon03 = SelectedTau->tauID("againstMuonTight");
      m_tauAgainstMuon04 = SelectedTau->tauID("againstMuonLoose2");
      m_tauAgainstMuon05 = SelectedTau->tauID("againstMuonMedium2");
      m_tauAgainstMuon06 = SelectedTau->tauID("againstMuonTight2");
      m_tauAgainstMuon07 = SelectedTau->tauID("againstMuonLoose3");
      m_tauAgainstMuon08 = SelectedTau->tauID("againstMuonTight3");
      m_tauAgainstMuon09 = SelectedTau->tauID("againstMuonMVAraw");
      m_tauAgainstMuon10 = SelectedTau->tauID("againstMuonLooseMVA");
      m_tauAgainstMuon11 = SelectedTau->tauID("againstMuonMediumMVA");
      m_tauAgainstMuon12 = SelectedTau->tauID("againstMuonTightMVA");
      m_tauPuCorrPtSum   = SelectedTau->tauID("puCorrPtSum");
      TreeUsualMuoTau->Fill();
    }
  }
	
  
  //ELE-TAU CHANNEL
  if(!isData && gentau.size()==2 && ((tauele1 && !taumuo2 && !tauele2) || (!taumuo1 && !tauele1 && tauele2))){

    float deltaRCleaned=99.; bool matchedCleaned = false;
    pat::TauCollection::const_iterator SelectedTauElTau;
    for (pat::TauCollection::const_iterator patTau = tauElTauHandle->begin(); patTau != tauElTauHandle->end(); ++patTau ) {
      //if(!(patTau->pt()>20)) continue;
      if(ROOT::Math::VectorUtil::DeltaR(gentauHad[0],patTau->p4())<0.5 && ROOT::Math::VectorUtil::DeltaR(gentauHad[0],patTau->p4())<deltaRCleaned){
	deltaRCleaned = ROOT::Math::VectorUtil::DeltaR(gentauHad[0],patTau->p4());
	SelectedTauElTau=patTau;
	matchedCleaned = true;
      }
    }
    float deltaRUsual=99.; bool matchedUsual = false;
    pat::TauCollection::const_iterator SelectedTau;
    for (pat::TauCollection::const_iterator patTau = tauHandle->begin(); patTau != tauHandle->end(); ++patTau ) {
      //if(!(patTau->pt()>20)) continue;
      if(ROOT::Math::VectorUtil::DeltaR(gentauHad[0],patTau->p4())<0.5 && ROOT::Math::VectorUtil::DeltaR(gentauHad[0],patTau->p4())<deltaRUsual){
	deltaRUsual = ROOT::Math::VectorUtil::DeltaR(gentauHad[0],patTau->p4());
	SelectedTau=patTau;
	matchedUsual = true;
      }
    } 
    
    float deltaREle=99.; bool matchedEle = false;
    pat::ElectronCollection::const_iterator recoEle;
    for(pat::ElectronCollection::const_iterator electron = eleH->begin(); electron != eleH->end(); ++electron) {
      if(!(electron->pt()>10)) continue;
      if(ROOT::Math::VectorUtil::DeltaR(genEle,electron->p4())<0.3 && ROOT::Math::VectorUtil::DeltaR(genEle,electron->p4())<deltaREle){
	deltaREle = ROOT::Math::VectorUtil::DeltaR(genEle,electron->p4());
	recoEle=electron;
	matchedEle = true;
      }
    }

    uno = 1;
    if(matchedCleaned) due=1;
    if(matchedCleaned && SelectedTauElTau->pt()>20) tre=1;
    if(matchedCleaned && SelectedTauElTau->pt()>20 && SelectedTauElTau->tauID("decayModeFindingNewDMs")>0.5) qua=1;
    if(matchedCleaned && SelectedTauElTau->pt()>20 && SelectedTauElTau->tauID("decayModeFindingNewDMs")>0.5 && SelectedTauElTau->tauID("againstMuonLoose")>0.5) cin=1;
    if(matchedCleaned && SelectedTauElTau->pt()>20 && SelectedTauElTau->tauID("decayModeFindingNewDMs")>0.5 && SelectedTauElTau->tauID("againstMuonLoose")>0.5 && 
       SelectedTauElTau->tauID("againstElectronLoose")>0.5) sei=1;
    if(matchedCleaned && SelectedTauElTau->pt()>20 && SelectedTauElTau->tauID("decayModeFindingNewDMs")>0.5 && SelectedTauElTau->tauID("againstMuonLoose")>0.5 && 
       SelectedTauElTau->tauID("againstElectronLoose")>0.5 && SelectedTauElTau->tauID("byVLooseIsolationMVA3newDMwoLT")>0.5) set=1;
    m_uno = uno;
    m_due = due;
    m_tre = tre;
    m_qua = qua;
    m_cin = cin;
    m_sei = sei;
    m_set = set;
    TreeEffEleTau->Fill();
    unoU = 1;
    if(matchedUsual) dueU=1;
    if(matchedUsual && SelectedTau->pt()>20) treU=1;
    if(matchedUsual && SelectedTau->pt()>20 && SelectedTau->tauID("decayModeFindingNewDMs")>0.5) quaU=1;
    if(matchedUsual && SelectedTau->pt()>20 && SelectedTau->tauID("decayModeFindingNewDMs")>0.5 && SelectedTau->tauID("againstMuonLoose")>0.5) cinU=1;
    if(matchedUsual && SelectedTau->pt()>20 && SelectedTau->tauID("decayModeFindingNewDMs")>0.5 && SelectedTau->tauID("againstMuonLoose")>0.5 && 
       SelectedTau->tauID("againstElectronLoose")>0.5) seiU=1;
    if(matchedUsual && SelectedTau->pt()>20 && SelectedTau->tauID("decayModeFindingNewDMs")>0.5 && SelectedTau->tauID("againstMuonLoose")>0.5 && 
       SelectedTau->tauID("againstElectronLoose")>0.5 && SelectedTau->tauID("byVLooseIsolationMVA3newDMwoLT")>0.5) setU=1;
    m_unoU = unoU;
    m_dueU = dueU;
    m_treU = treU;
    m_quaU = quaU;
    m_cinU = cinU;
    m_seiU = seiU;
    m_setU = setU;
    TreeEffUsualEleTau->Fill();

    if(matchedCleaned){

      // "met->begin()->"   ---->   "CorrectedMet->"
      float dRmin=999.; pat::JetCollection::const_iterator TauJet; bool TauMatched = false;
      for(pat::JetCollection::const_iterator ak5 = ak5jetCands->begin(); ak5 != ak5jetCands->end(); ++ak5) {
	//if(!(ak5->pt()>20)) continue;
	if(ROOT::Math::VectorUtil::DeltaR(SelectedTauElTau->p4Jet(),ak5->p4())<0.3 && ROOT::Math::VectorUtil::DeltaR(SelectedTauElTau->p4Jet(),ak5->p4())<dRmin){
	  dRmin=ROOT::Math::VectorUtil::DeltaR(SelectedTauElTau->p4Jet(),ak5->p4());
	  TauJet=ak5;
	  TauMatched=true;
	}
      }

      float deltaRReco=99.;
      if(matchedEle) deltaRReco=ROOT::Math::VectorUtil::DeltaR(SelectedTauElTau->p4(),recoEle->p4());

      /*
      //corrPFMET = rawPFMET + sum(rawJetPt-corrJetPt) [*]
      float PX = 0.;
      float PY = 0.;
      //pat::MET *CorrectedMet;
      //TLorentzVector CorrectedMet_prov; 
      if(TauMatched){
	PX = met->begin()->px() + TauJet->px()*(1 - TauJet->jecFactor(0));
	PY = met->begin()->py() + TauJet->py()*(1 - TauJet->jecFactor(0));
      } else {
	PX = met->begin()->px();
	PY = met->begin()->py();
      }
      reco::Candidate::Vector CorrectedMet(PX,PY,0.);
      reco::Candidate::LorentzVector CorrectedMetP4(PX,PY,0.,sqrt(PX*PX+PY*PY));

      reco::Candidate::LorentzVector dMET = CorrectedMetP4 - genMets->begin()->p4();
      double cosHiggsPhi = TMath::Cos(genHiggs.phi());
      double sinHiggsPhi = TMath::Sin(genHiggs.phi());
      double px = dMET.px();
      double py = dMET.py();
      double pxInZetaFrame =  cosHiggsPhi*px + sinHiggsPhi*py;
      double pyInZetaFrame = -sinHiggsPhi*px + cosHiggsPhi*py;
      reco::Candidate::LorentzVector metErr_rotated(pxInZetaFrame, pyInZetaFrame, dMET.pz(), dMET.E());
      double metErrParl = metErr_rotated.px();
      double metErrPerp = metErr_rotated.py();

      TMatrixD metCov(2, 2); // PFMET significance matrix
      metCov[0][0] = (metRaw->front() ).getSignificanceMatrix()(0,0);
      metCov[1][0] = (metRaw->front() ).getSignificanceMatrix()(1,0);
      metCov[0][1] = (metRaw->front() ).getSignificanceMatrix()(0,1);
      metCov[1][1] = (metRaw->front() ).getSignificanceMatrix()(1,1);
      double cov00 = metCov(0,0);
      double cov01 = metCov(0,1);
      double cov10 = metCov(1,0);
      double cov11 = metCov(1,1);
      TMatrixD metCov_rotated(2,2);
      metCov_rotated(0,0) =  cosHiggsPhi*cosHiggsPhi*cov00 + cosHiggsPhi*sinHiggsPhi*cov01 + sinHiggsPhi*cosHiggsPhi*cov10 + sinHiggsPhi*sinHiggsPhi*cov11;
      metCov_rotated(0,1) = -sinHiggsPhi*cosHiggsPhi*cov00 + cosHiggsPhi*cosHiggsPhi*cov01 - sinHiggsPhi*sinHiggsPhi*cov10 + cosHiggsPhi*sinHiggsPhi*cov11;
      metCov_rotated(1,0) = -cosHiggsPhi*sinHiggsPhi*cov00 - sinHiggsPhi*sinHiggsPhi*cov01 + cosHiggsPhi*cosHiggsPhi*cov10 + sinHiggsPhi*cosHiggsPhi*cov11;
      metCov_rotated(1,1) =  sinHiggsPhi*sinHiggsPhi*cov00 - cosHiggsPhi*sinHiggsPhi*cov01 - sinHiggsPhi*cosHiggsPhi*cov10 + cosHiggsPhi*cosHiggsPhi*cov11;
      double metSigmaParl = TMath::Sqrt(TMath::Abs(metCov_rotated(0, 0)));
      double metSigmaPerp = TMath::Sqrt(TMath::Abs(metCov_rotated(1, 1)));
      float MEtPullParlZ=-9999.;
      float MEtPullPerpZ=-9999.;
      if(metSigmaParl > 0.) MEtPullParlZ=metErrParl/metSigmaParl;  
      if(metSigmaPerp > 0.) MEtPullPerpZ=metErrPerp/metSigmaPerp;

      //SVFIT
      float MassSVFit = -1;
      float ptSVFit = -1;
      if(matchedEle){
	std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
	measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay, SelectedTauElTau->p4()));
	measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, recoEle->p4()));
	NSVfitStandaloneAlgorithm algo(measuredTauLeptons, CorrectedMet, metCov, 0);
	algo.addLogM(false);
	algo.integrateMarkovChain();
	if(algo.pt()>0){
	  MassSVFit=algo.getMass();
	  ptSVFit=algo.pt();
	}
      }
      
      m_MassSVFit=MassSVFit;
      m_ptSVFit=ptSVFit;
      m_MEtSigmaParlZ = metSigmaParl;
      m_MEtSigmaPerpZ = metSigmaPerp;
      m_MEtPullParlZ  = MEtPullParlZ;
      m_MEtPullPerpZ  = MEtPullPerpZ;
      m_metErrParl    = metErrParl;
      m_metErrPerp    = metErrPerp;
      m_metGen        = genMets->begin()->pt();
      m_met           = sqrt(PX*PX+PY*PY);
      */

      m_tauPtGen      = gentauHad[0].pt();
      m_tauEtaGen     = gentauHad[0].eta();
      m_tauMassGen    = gentauHad[0].mass();
      m_tauMass       = SelectedTauElTau->mass();
      m_deltaRReco    = deltaRReco;
      m_deltaRGen     = ROOT::Math::VectorUtil::DeltaR(gentauHad[0],genEle);
      m_tauPt=SelectedTauElTau->pt();
      m_tauEta=SelectedTauElTau->eta();
      m_deltaR=deltaRCleaned;
      m_trigger=(int)isFired_HLT;
      m_trigger320=(int)isFired_HLT_PFJet320;
      m_trigger650=(int)isFired_HLT_HT650;
      m_tauDecay01 = SelectedTauElTau->tauID("decayModeFinding");
      m_tauDecay02 = SelectedTauElTau->tauID("decayModeFindingNewDMs");
      m_tauDecay03 = SelectedTauElTau->tauID("decayModeFindingOldDMs");
      m_tauIso01   = SelectedTauElTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr");
      m_tauIso02   = SelectedTauElTau->tauID("byLooseIsolation");
      m_tauIso03   = SelectedTauElTau->tauID("byLooseCombinedIsolationDeltaBetaCorr");
      m_tauIso04   = SelectedTauElTau->tauID("byMediumCombinedIsolationDeltaBetaCorr");
      m_tauIso05   = SelectedTauElTau->tauID("byTightCombinedIsolationDeltaBetaCorr");
      m_tauIso06   = SelectedTauElTau->tauID("byCombinedIsolationDeltaBetaCorrRaw");
      m_tauIso07   = SelectedTauElTau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
      m_tauIso08   = SelectedTauElTau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
      m_tauIso09   = SelectedTauElTau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
      m_tauIso10   = SelectedTauElTau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
      m_tauIso11   = SelectedTauElTau->tauID("chargedIsoPtSum");
      m_tauIso12   = SelectedTauElTau->tauID("neutralIsoPtSum");
      m_tauIso13   = SelectedTauElTau->tauID("byIsolationMVA3oldDMwoLTraw");
      m_tauIso14   = SelectedTauElTau->tauID("byVLooseIsolationMVA3oldDMwoLT");
      m_tauIso15   = SelectedTauElTau->tauID("byLooseIsolationMVA3oldDMwoLT");
      m_tauIso16   = SelectedTauElTau->tauID("byMediumIsolationMVA3oldDMwoLT");
      m_tauIso17   = SelectedTauElTau->tauID("byTightIsolationMVA3oldDMwoLT");
      m_tauIso18   = SelectedTauElTau->tauID("byVTightIsolationMVA3oldDMwoLT");
      m_tauIso19   = SelectedTauElTau->tauID("byVVTightIsolationMVA3oldDMwoLT");
      m_tauIso20   = SelectedTauElTau->tauID("byIsolationMVA3oldDMwLTraw");
      m_tauIso21   = SelectedTauElTau->tauID("byVLooseIsolationMVA3oldDMwLT");
      m_tauIso22   = SelectedTauElTau->tauID("byLooseIsolationMVA3oldDMwLT");
      m_tauIso23   = SelectedTauElTau->tauID("byMediumIsolationMVA3oldDMwLT");
      m_tauIso24   = SelectedTauElTau->tauID("byTightIsolationMVA3oldDMwLT");
      m_tauIso25   = SelectedTauElTau->tauID("byVTightIsolationMVA3oldDMwLT");
      m_tauIso26   = SelectedTauElTau->tauID("byVVTightIsolationMVA3oldDMwLT");
      m_tauIso27   = SelectedTauElTau->tauID("byIsolationMVA3newDMwoLTraw");
      m_tauIso28   = SelectedTauElTau->tauID("byVLooseIsolationMVA3newDMwoLT");
      m_tauIso29   = SelectedTauElTau->tauID("byLooseIsolationMVA3newDMwoLT");
      m_tauIso30   = SelectedTauElTau->tauID("byMediumIsolationMVA3newDMwoLT");
      m_tauIso31   = SelectedTauElTau->tauID("byTightIsolationMVA3newDMwoLT");
      m_tauIso32   = SelectedTauElTau->tauID("byVTightIsolationMVA3newDMwoLT");
      m_tauIso33   = SelectedTauElTau->tauID("byVVTightIsolationMVA3newDMwoLT");
      m_tauIso34   = SelectedTauElTau->tauID("byIsolationMVA3newDMwLTraw");
      m_tauIso35   = SelectedTauElTau->tauID("byVLooseIsolationMVA3newDMwLT");
      m_tauIso36   = SelectedTauElTau->tauID("byLooseIsolationMVA3newDMwLT");
      m_tauIso37   = SelectedTauElTau->tauID("byMediumIsolationMVA3newDMwLT");
      m_tauIso38   = SelectedTauElTau->tauID("byTightIsolationMVA3newDMwLT");
      m_tauIso39   = SelectedTauElTau->tauID("byVTightIsolationMVA3newDMwLT");
      m_tauIso40   = SelectedTauElTau->tauID("byVVTightIsolationMVA3newDMwLT");
      m_tauAgainstElectron01 = SelectedTauElTau->tauID("againstElectronLoose");
      m_tauAgainstElectron02 = SelectedTauElTau->tauID("againstElectronMedium");
      m_tauAgainstElectron03 = SelectedTauElTau->tauID("againstElectronTight");
      m_tauAgainstElectron04 = SelectedTauElTau->tauID("againstElectronMVA5raw");
      m_tauAgainstElectron05 = SelectedTauElTau->tauID("againstElectronMVA5category");
      m_tauAgainstElectron06 = SelectedTauElTau->tauID("againstElectronVLooseMVA5");
      m_tauAgainstElectron07 = SelectedTauElTau->tauID("againstElectronLooseMVA5");
      m_tauAgainstElectron08 = SelectedTauElTau->tauID("againstElectronMediumMVA5");
      m_tauAgainstElectron09 = SelectedTauElTau->tauID("againstElectronTightMVA5");
      m_tauAgainstElectron10 = SelectedTauElTau->tauID("againstElectronVTightMVA5");
      m_tauAgainstElectron11 = SelectedTauElTau->tauID("againstElectronDeadECAL");
      m_tauAgainstMuon01 = SelectedTauElTau->tauID("againstMuonLoose");
      m_tauAgainstMuon02 = SelectedTauElTau->tauID("againstMuonMedium");
      m_tauAgainstMuon03 = SelectedTauElTau->tauID("againstMuonTight");
      m_tauAgainstMuon04 = SelectedTauElTau->tauID("againstMuonLoose2");
      m_tauAgainstMuon05 = SelectedTauElTau->tauID("againstMuonMedium2");
      m_tauAgainstMuon06 = SelectedTauElTau->tauID("againstMuonTight2");
      m_tauAgainstMuon07 = SelectedTauElTau->tauID("againstMuonLoose3");
      m_tauAgainstMuon08 = SelectedTauElTau->tauID("againstMuonTight3");
      m_tauAgainstMuon09 = SelectedTauElTau->tauID("againstMuonMVAraw");
      m_tauAgainstMuon10 = SelectedTauElTau->tauID("againstMuonLooseMVA");
      m_tauAgainstMuon11 = SelectedTauElTau->tauID("againstMuonMediumMVA");
      m_tauAgainstMuon12 = SelectedTauElTau->tauID("againstMuonTightMVA");
      m_tauPuCorrPtSum   = SelectedTauElTau->tauID("puCorrPtSum");
      TreeEleTau->Fill();
    }

    if(matchedUsual){

      // "met->begin()->"   ---->   "CorrectedMet->"
      float dRmin=999.; pat::JetCollection::const_iterator TauJet; bool TauMatched = false;
      for(pat::JetCollection::const_iterator ak5 = ak5jetCands->begin(); ak5 != ak5jetCands->end(); ++ak5) {
	//if(!(ak5->pt()>20)) continue;
	if(ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4Jet(),ak5->p4())<0.3 && ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4Jet(),ak5->p4())<dRmin){
	  dRmin=ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4Jet(),ak5->p4());
	  TauJet=ak5;
	  TauMatched=true;
	}
      }

      float deltaRReco=99.;
      if(matchedEle) deltaRReco=ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(),recoEle->p4());

      /*
      //corrPFMET = rawPFMET + sum(rawJetPt-corrJetPt) [*]
      float PX = 0.;
      float PY = 0.;
      //pat::MET *CorrectedMet;
      //TLorentzVector CorrectedMet_prov; 
      if(TauMatched){
	PX = met->begin()->px() + TauJet->px()*(1 - TauJet->jecFactor(0));
	PY = met->begin()->py() + TauJet->py()*(1 - TauJet->jecFactor(0));
      } else {
	PX = met->begin()->px();
	PY = met->begin()->py();
      }
      reco::Candidate::Vector CorrectedMet(PX,PY,0.);
      reco::Candidate::LorentzVector CorrectedMetP4(PX,PY,0.,sqrt(PX*PX+PY*PY));

      reco::Candidate::LorentzVector dMET = CorrectedMetP4 - genMets->begin()->p4();
      double cosHiggsPhi = TMath::Cos(genHiggs.phi());
      double sinHiggsPhi = TMath::Sin(genHiggs.phi());
      double px = dMET.px();
      double py = dMET.py();
      double pxInZetaFrame =  cosHiggsPhi*px + sinHiggsPhi*py;
      double pyInZetaFrame = -sinHiggsPhi*px + cosHiggsPhi*py;
      reco::Candidate::LorentzVector metErr_rotated(pxInZetaFrame, pyInZetaFrame, dMET.pz(), dMET.E());
      double metErrParl = metErr_rotated.px();
      double metErrPerp = metErr_rotated.py();

      TMatrixD metCov(2, 2); // PFMET significance matrix
      metCov[0][0] = (metRaw->front() ).getSignificanceMatrix()(0,0);
      metCov[1][0] = (metRaw->front() ).getSignificanceMatrix()(1,0);
      metCov[0][1] = (metRaw->front() ).getSignificanceMatrix()(0,1);
      metCov[1][1] = (metRaw->front() ).getSignificanceMatrix()(1,1);
      double cov00 = metCov(0,0);
      double cov01 = metCov(0,1);
      double cov10 = metCov(1,0);
      double cov11 = metCov(1,1);
      TMatrixD metCov_rotated(2,2);
      metCov_rotated(0,0) =  cosHiggsPhi*cosHiggsPhi*cov00 + cosHiggsPhi*sinHiggsPhi*cov01 + sinHiggsPhi*cosHiggsPhi*cov10 + sinHiggsPhi*sinHiggsPhi*cov11;
      metCov_rotated(0,1) = -sinHiggsPhi*cosHiggsPhi*cov00 + cosHiggsPhi*cosHiggsPhi*cov01 - sinHiggsPhi*sinHiggsPhi*cov10 + cosHiggsPhi*sinHiggsPhi*cov11;
      metCov_rotated(1,0) = -cosHiggsPhi*sinHiggsPhi*cov00 - sinHiggsPhi*sinHiggsPhi*cov01 + cosHiggsPhi*cosHiggsPhi*cov10 + sinHiggsPhi*cosHiggsPhi*cov11;
      metCov_rotated(1,1) =  sinHiggsPhi*sinHiggsPhi*cov00 - cosHiggsPhi*sinHiggsPhi*cov01 - sinHiggsPhi*cosHiggsPhi*cov10 + cosHiggsPhi*cosHiggsPhi*cov11;
      double metSigmaParl = TMath::Sqrt(TMath::Abs(metCov_rotated(0, 0)));
      double metSigmaPerp = TMath::Sqrt(TMath::Abs(metCov_rotated(1, 1)));
      float MEtPullParlZ=-9999.;
      float MEtPullPerpZ=-9999.;
      if(metSigmaParl > 0.) MEtPullParlZ=metErrParl/metSigmaParl;  
      if(metSigmaPerp > 0.) MEtPullPerpZ=metErrPerp/metSigmaPerp;

      //SVFIT
      float MassSVFit = -1;
      float ptSVFit = -1;
      if(matchedEle){
	std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
	measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay, SelectedTau->p4()));
	measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, recoEle->p4()));
	NSVfitStandaloneAlgorithm algo(measuredTauLeptons, CorrectedMet, metCov, 0);
	algo.addLogM(false);
	algo.integrateMarkovChain();
	if(algo.pt()>0){
	  MassSVFit=algo.getMass();
	  ptSVFit=algo.pt();
	}
      }
      
      m_MassSVFit=MassSVFit;
      m_ptSVFit=ptSVFit;      
      m_MEtSigmaParlZ = metSigmaParl;
      m_MEtSigmaPerpZ = metSigmaPerp;
      m_MEtPullParlZ  = MEtPullParlZ;
      m_MEtPullPerpZ  = MEtPullPerpZ;
      m_metErrParl    = metErrParl;
      m_metErrPerp    = metErrPerp;
      m_metGen        = genMets->begin()->pt();
      m_met           = sqrt(PX*PX+PY*PY);
      */

      m_tauPtGen      = gentauHad[0].pt();
      m_tauEtaGen     = gentauHad[0].eta();
      m_tauMassGen    = gentauHad[0].mass();
      m_tauMass       = SelectedTau->mass();
      m_deltaRReco    = deltaRReco;
      m_deltaRGen     = ROOT::Math::VectorUtil::DeltaR(gentauHad[0],genEle);
      m_tauPt=SelectedTau->pt();
      m_tauEta=SelectedTau->eta();
      m_deltaR=deltaRUsual;
      m_trigger=(int)isFired_HLT;
      m_trigger320=(int)isFired_HLT_PFJet320;
      m_trigger650=(int)isFired_HLT_HT650;
      m_tauDecay01 = SelectedTau->tauID("decayModeFinding");
      m_tauDecay02 = SelectedTau->tauID("decayModeFindingNewDMs");
      m_tauDecay03 = SelectedTau->tauID("decayModeFindingOldDMs");
      m_tauIso01   = SelectedTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr");
      m_tauIso02   = SelectedTau->tauID("byLooseIsolation");
      m_tauIso03   = SelectedTau->tauID("byLooseCombinedIsolationDeltaBetaCorr");
      m_tauIso04   = SelectedTau->tauID("byMediumCombinedIsolationDeltaBetaCorr");
      m_tauIso05   = SelectedTau->tauID("byTightCombinedIsolationDeltaBetaCorr");
      m_tauIso06   = SelectedTau->tauID("byCombinedIsolationDeltaBetaCorrRaw");
      m_tauIso07   = SelectedTau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
      m_tauIso08   = SelectedTau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
      m_tauIso09   = SelectedTau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
      m_tauIso10   = SelectedTau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
      m_tauIso11   = SelectedTau->tauID("chargedIsoPtSum");
      m_tauIso12   = SelectedTau->tauID("neutralIsoPtSum");
      m_tauIso13   = SelectedTau->tauID("byIsolationMVA3oldDMwoLTraw");
      m_tauIso14   = SelectedTau->tauID("byVLooseIsolationMVA3oldDMwoLT");
      m_tauIso15   = SelectedTau->tauID("byLooseIsolationMVA3oldDMwoLT");
      m_tauIso16   = SelectedTau->tauID("byMediumIsolationMVA3oldDMwoLT");
      m_tauIso17   = SelectedTau->tauID("byTightIsolationMVA3oldDMwoLT");
      m_tauIso18   = SelectedTau->tauID("byVTightIsolationMVA3oldDMwoLT");
      m_tauIso19   = SelectedTau->tauID("byVVTightIsolationMVA3oldDMwoLT");
      m_tauIso20   = SelectedTau->tauID("byIsolationMVA3oldDMwLTraw");
      m_tauIso21   = SelectedTau->tauID("byVLooseIsolationMVA3oldDMwLT");
      m_tauIso22   = SelectedTau->tauID("byLooseIsolationMVA3oldDMwLT");
      m_tauIso23   = SelectedTau->tauID("byMediumIsolationMVA3oldDMwLT");
      m_tauIso24   = SelectedTau->tauID("byTightIsolationMVA3oldDMwLT");
      m_tauIso25   = SelectedTau->tauID("byVTightIsolationMVA3oldDMwLT");
      m_tauIso26   = SelectedTau->tauID("byVVTightIsolationMVA3oldDMwLT");
      m_tauIso27   = SelectedTau->tauID("byIsolationMVA3newDMwoLTraw");
      m_tauIso28   = SelectedTau->tauID("byVLooseIsolationMVA3newDMwoLT");
      m_tauIso29   = SelectedTau->tauID("byLooseIsolationMVA3newDMwoLT");
      m_tauIso30   = SelectedTau->tauID("byMediumIsolationMVA3newDMwoLT");
      m_tauIso31   = SelectedTau->tauID("byTightIsolationMVA3newDMwoLT");
      m_tauIso32   = SelectedTau->tauID("byVTightIsolationMVA3newDMwoLT");
      m_tauIso33   = SelectedTau->tauID("byVVTightIsolationMVA3newDMwoLT");
      m_tauIso34   = SelectedTau->tauID("byIsolationMVA3newDMwLTraw");
      m_tauIso35   = SelectedTau->tauID("byVLooseIsolationMVA3newDMwLT");
      m_tauIso36   = SelectedTau->tauID("byLooseIsolationMVA3newDMwLT");
      m_tauIso37   = SelectedTau->tauID("byMediumIsolationMVA3newDMwLT");
      m_tauIso38   = SelectedTau->tauID("byTightIsolationMVA3newDMwLT");
      m_tauIso39   = SelectedTau->tauID("byVTightIsolationMVA3newDMwLT");
      m_tauIso40   = SelectedTau->tauID("byVVTightIsolationMVA3newDMwLT");
      m_tauAgainstElectron01 = SelectedTau->tauID("againstElectronLoose");
      m_tauAgainstElectron02 = SelectedTau->tauID("againstElectronMedium");
      m_tauAgainstElectron03 = SelectedTau->tauID("againstElectronTight");
      m_tauAgainstElectron04 = SelectedTau->tauID("againstElectronMVA5raw");
      m_tauAgainstElectron05 = SelectedTau->tauID("againstElectronMVA5category");
      m_tauAgainstElectron06 = SelectedTau->tauID("againstElectronVLooseMVA5");
      m_tauAgainstElectron07 = SelectedTau->tauID("againstElectronLooseMVA5");
      m_tauAgainstElectron08 = SelectedTau->tauID("againstElectronMediumMVA5");
      m_tauAgainstElectron09 = SelectedTau->tauID("againstElectronTightMVA5");
      m_tauAgainstElectron10 = SelectedTau->tauID("againstElectronVTightMVA5");
      m_tauAgainstElectron11 = SelectedTau->tauID("againstElectronDeadECAL");
      m_tauAgainstMuon01 = SelectedTau->tauID("againstMuonLoose");
      m_tauAgainstMuon02 = SelectedTau->tauID("againstMuonMedium");
      m_tauAgainstMuon03 = SelectedTau->tauID("againstMuonTight");
      m_tauAgainstMuon04 = SelectedTau->tauID("againstMuonLoose2");
      m_tauAgainstMuon05 = SelectedTau->tauID("againstMuonMedium2");
      m_tauAgainstMuon06 = SelectedTau->tauID("againstMuonTight2");
      m_tauAgainstMuon07 = SelectedTau->tauID("againstMuonLoose3");
      m_tauAgainstMuon08 = SelectedTau->tauID("againstMuonTight3");
      m_tauAgainstMuon09 = SelectedTau->tauID("againstMuonMVAraw");
      m_tauAgainstMuon10 = SelectedTau->tauID("againstMuonLooseMVA");
      m_tauAgainstMuon11 = SelectedTau->tauID("againstMuonMediumMVA");
      m_tauAgainstMuon12 = SelectedTau->tauID("againstMuonTightMVA");
      m_tauPuCorrPtSum   = SelectedTau->tauID("puCorrPtSum");
      TreeUsualEleTau->Fill();
    }
  }	
  
  if(!isData){
    vector<math::PtEtaPhiELorentzVector> genele_QCD;
    vector<math::PtEtaPhiELorentzVector> genmuo_QCD;
    math::PtEtaPhiELorentzVector gen_prov;
    bool realTau = false;
    for(size_t ngenPart=0; ngenPart<genParts->size(); ngenPart++){
      const reco::GenParticle & genPart = (*genParts)[ngenPart];
      gen_prov=genPart.p4();
      if(abs(genPart.pdgId())==15) realTau=true;
      if(abs(genPart.pdgId())==11 && genPart.status()!=3) genele_QCD.push_back(gen_prov);
      if(abs(genPart.pdgId())==13 && genPart.status()!=3) genmuo_QCD.push_back(gen_prov);
    }
    
    //QCD - USUAL COLLECTION
    for (pat::TauCollection::const_iterator patTau = tauHandle->begin(); patTau != tauHandle->end(); ++patTau ) {
      float deltaRGen=99.; bool genJetMatch = false;
      for(vector<reco::GenJet>::const_iterator genjet=genjets->begin(); genjet!=genjets->end(); genjet++){
	if(ROOT::Math::VectorUtil::DeltaR(genjet->p4(),patTau->p4Jet())<0.3 && ROOT::Math::VectorUtil::DeltaR(genjet->p4(),patTau->p4Jet())<deltaRGen){
	  genJetMatch = true;
	  deltaRGen = ROOT::Math::VectorUtil::DeltaR(genjet->p4(),patTau->p4Jet());
	}
      }
      
      float taudzvertex = -9999;
      if(vertices.isValid() && vertices->size())taudzvertex = TMath::Abs(patTau->vertex().z()-(*vertices)[0].position().z());

      //if(!(patTau->pt()>20))    continue;
      if(!(realTau==false))     continue;
      if(!(genJetMatch==true))  continue;
      m_deltaR=deltaRGen;
      m_trigger=(int)isFired_HLT;
      m_trigger320=(int)isFired_HLT_PFJet320;
      m_trigger650=(int)isFired_HLT_HT650;
      m_tauPt=patTau->p4Jet().pt();
      m_tauEta=patTau->p4Jet().eta();
      m_tauVertexZ=taudzvertex;
      m_tauDecay01 = patTau->tauID("decayModeFinding");
      m_tauDecay02 = patTau->tauID("decayModeFindingNewDMs");
      m_tauDecay03 = patTau->tauID("decayModeFindingOldDMs");
      m_tauIso01   = patTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr");
      m_tauIso02   = patTau->tauID("byLooseIsolation");
      m_tauIso03   = patTau->tauID("byLooseCombinedIsolationDeltaBetaCorr");
      m_tauIso04   = patTau->tauID("byMediumCombinedIsolationDeltaBetaCorr");
      m_tauIso05   = patTau->tauID("byTightCombinedIsolationDeltaBetaCorr");
      m_tauIso06   = patTau->tauID("byCombinedIsolationDeltaBetaCorrRaw");
      m_tauIso07   = patTau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
      m_tauIso08   = patTau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
      m_tauIso09   = patTau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
      m_tauIso10   = patTau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
      m_tauIso11   = patTau->tauID("chargedIsoPtSum");
      m_tauIso12   = patTau->tauID("neutralIsoPtSum");
      m_tauIso13   = patTau->tauID("byIsolationMVA3oldDMwoLTraw");
      m_tauIso14   = patTau->tauID("byVLooseIsolationMVA3oldDMwoLT");
      m_tauIso15   = patTau->tauID("byLooseIsolationMVA3oldDMwoLT");
      m_tauIso16   = patTau->tauID("byMediumIsolationMVA3oldDMwoLT");
      m_tauIso17   = patTau->tauID("byTightIsolationMVA3oldDMwoLT");
      m_tauIso18   = patTau->tauID("byVTightIsolationMVA3oldDMwoLT");
      m_tauIso19   = patTau->tauID("byVVTightIsolationMVA3oldDMwoLT");
      m_tauIso20   = patTau->tauID("byIsolationMVA3oldDMwLTraw");
      m_tauIso21   = patTau->tauID("byVLooseIsolationMVA3oldDMwLT");
      m_tauIso22   = patTau->tauID("byLooseIsolationMVA3oldDMwLT");
      m_tauIso23   = patTau->tauID("byMediumIsolationMVA3oldDMwLT");
      m_tauIso24   = patTau->tauID("byTightIsolationMVA3oldDMwLT");
      m_tauIso25   = patTau->tauID("byVTightIsolationMVA3oldDMwLT");
      m_tauIso26   = patTau->tauID("byVVTightIsolationMVA3oldDMwLT");
      m_tauIso27   = patTau->tauID("byIsolationMVA3newDMwoLTraw");
      m_tauIso28   = patTau->tauID("byVLooseIsolationMVA3newDMwoLT");
      m_tauIso29   = patTau->tauID("byLooseIsolationMVA3newDMwoLT");
      m_tauIso30   = patTau->tauID("byMediumIsolationMVA3newDMwoLT");
      m_tauIso31   = patTau->tauID("byTightIsolationMVA3newDMwoLT");
      m_tauIso32   = patTau->tauID("byVTightIsolationMVA3newDMwoLT");
      m_tauIso33   = patTau->tauID("byVVTightIsolationMVA3newDMwoLT");
      m_tauIso34   = patTau->tauID("byIsolationMVA3newDMwLTraw");
      m_tauIso35   = patTau->tauID("byVLooseIsolationMVA3newDMwLT");
      m_tauIso36   = patTau->tauID("byLooseIsolationMVA3newDMwLT");
      m_tauIso37   = patTau->tauID("byMediumIsolationMVA3newDMwLT");
      m_tauIso38   = patTau->tauID("byTightIsolationMVA3newDMwLT");
      m_tauIso39   = patTau->tauID("byVTightIsolationMVA3newDMwLT");
      m_tauIso40   = patTau->tauID("byVVTightIsolationMVA3newDMwLT");
      m_tauAgainstElectron01 = patTau->tauID("againstElectronLoose");
      m_tauAgainstElectron02 = patTau->tauID("againstElectronMedium");
      m_tauAgainstElectron03 = patTau->tauID("againstElectronTight");
      m_tauAgainstElectron04 = patTau->tauID("againstElectronMVA5raw");
      m_tauAgainstElectron05 = patTau->tauID("againstElectronMVA5category");
      m_tauAgainstElectron06 = patTau->tauID("againstElectronVLooseMVA5");
      m_tauAgainstElectron07 = patTau->tauID("againstElectronLooseMVA5");
      m_tauAgainstElectron08 = patTau->tauID("againstElectronMediumMVA5");
      m_tauAgainstElectron09 = patTau->tauID("againstElectronTightMVA5");
      m_tauAgainstElectron10 = patTau->tauID("againstElectronVTightMVA5");
      m_tauAgainstElectron11 = patTau->tauID("againstElectronDeadECAL");
      m_tauAgainstMuon01 = patTau->tauID("againstMuonLoose");
      m_tauAgainstMuon02 = patTau->tauID("againstMuonMedium");
      m_tauAgainstMuon03 = patTau->tauID("againstMuonTight");
      m_tauAgainstMuon04 = patTau->tauID("againstMuonLoose2");
      m_tauAgainstMuon05 = patTau->tauID("againstMuonMedium2");
      m_tauAgainstMuon06 = patTau->tauID("againstMuonTight2");
      m_tauAgainstMuon07 = patTau->tauID("againstMuonLoose3");
      m_tauAgainstMuon08 = patTau->tauID("againstMuonTight3");
      m_tauAgainstMuon09 = patTau->tauID("againstMuonMVAraw");
      m_tauAgainstMuon10 = patTau->tauID("againstMuonLooseMVA");
      m_tauAgainstMuon11 = patTau->tauID("againstMuonMediumMVA");
      m_tauAgainstMuon12 = patTau->tauID("againstMuonTightMVA");
      m_tauPuCorrPtSum   = patTau->tauID("puCorrPtSum");
      TreeQCDUsual->Fill();
    }

    //QCD - MUOTAU COLLECTION
    for (pat::TauCollection::const_iterator patTau = tauMuTauHandle->begin(); patTau != tauMuTauHandle->end(); ++patTau ) {
      float deltaRGen=99.; bool genJetMatch = false;
      for(vector<reco::GenJet>::const_iterator genjet=genjets->begin(); genjet!=genjets->end(); genjet++){
	if(ROOT::Math::VectorUtil::DeltaR(genjet->p4(),patTau->p4Jet())<0.3 && ROOT::Math::VectorUtil::DeltaR(genjet->p4(),patTau->p4Jet())<deltaRGen){
	  genJetMatch = true;
	  deltaRGen = ROOT::Math::VectorUtil::DeltaR(genjet->p4(),patTau->p4Jet());
	}
      }
      
      float taudzvertex = -9999;
      if(vertices.isValid() && vertices->size())taudzvertex = TMath::Abs(patTau->vertex().z()-(*vertices)[0].position().z());

      //if(!(patTau->pt()>20))    continue;
      if(!(realTau==false))     continue;
      if(!(genJetMatch==true))  continue;
      m_deltaR=deltaRGen;
      m_trigger=(int)isFired_HLT;
      m_trigger320=(int)isFired_HLT_PFJet320;
      m_trigger650=(int)isFired_HLT_HT650;
      m_tauPt=patTau->p4Jet().pt();
      m_tauEta=patTau->p4Jet().eta();
      m_tauVertexZ=taudzvertex;
      m_tauDecay01 = patTau->tauID("decayModeFinding");
      m_tauDecay02 = patTau->tauID("decayModeFindingNewDMs");
      m_tauDecay03 = patTau->tauID("decayModeFindingOldDMs");
      m_tauIso01   = patTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr");
      m_tauIso02   = patTau->tauID("byLooseIsolation");
      m_tauIso03   = patTau->tauID("byLooseCombinedIsolationDeltaBetaCorr");
      m_tauIso04   = patTau->tauID("byMediumCombinedIsolationDeltaBetaCorr");
      m_tauIso05   = patTau->tauID("byTightCombinedIsolationDeltaBetaCorr");
      m_tauIso06   = patTau->tauID("byCombinedIsolationDeltaBetaCorrRaw");
      m_tauIso07   = patTau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
      m_tauIso08   = patTau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
      m_tauIso09   = patTau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
      m_tauIso10   = patTau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
      m_tauIso11   = patTau->tauID("chargedIsoPtSum");
      m_tauIso12   = patTau->tauID("neutralIsoPtSum");
      m_tauIso13   = patTau->tauID("byIsolationMVA3oldDMwoLTraw");
      m_tauIso14   = patTau->tauID("byVLooseIsolationMVA3oldDMwoLT");
      m_tauIso15   = patTau->tauID("byLooseIsolationMVA3oldDMwoLT");
      m_tauIso16   = patTau->tauID("byMediumIsolationMVA3oldDMwoLT");
      m_tauIso17   = patTau->tauID("byTightIsolationMVA3oldDMwoLT");
      m_tauIso18   = patTau->tauID("byVTightIsolationMVA3oldDMwoLT");
      m_tauIso19   = patTau->tauID("byVVTightIsolationMVA3oldDMwoLT");
      m_tauIso20   = patTau->tauID("byIsolationMVA3oldDMwLTraw");
      m_tauIso21   = patTau->tauID("byVLooseIsolationMVA3oldDMwLT");
      m_tauIso22   = patTau->tauID("byLooseIsolationMVA3oldDMwLT");
      m_tauIso23   = patTau->tauID("byMediumIsolationMVA3oldDMwLT");
      m_tauIso24   = patTau->tauID("byTightIsolationMVA3oldDMwLT");
      m_tauIso25   = patTau->tauID("byVTightIsolationMVA3oldDMwLT");
      m_tauIso26   = patTau->tauID("byVVTightIsolationMVA3oldDMwLT");
      m_tauIso27   = patTau->tauID("byIsolationMVA3newDMwoLTraw");
      m_tauIso28   = patTau->tauID("byVLooseIsolationMVA3newDMwoLT");
      m_tauIso29   = patTau->tauID("byLooseIsolationMVA3newDMwoLT");
      m_tauIso30   = patTau->tauID("byMediumIsolationMVA3newDMwoLT");
      m_tauIso31   = patTau->tauID("byTightIsolationMVA3newDMwoLT");
      m_tauIso32   = patTau->tauID("byVTightIsolationMVA3newDMwoLT");
      m_tauIso33   = patTau->tauID("byVVTightIsolationMVA3newDMwoLT");
      m_tauIso34   = patTau->tauID("byIsolationMVA3newDMwLTraw");
      m_tauIso35   = patTau->tauID("byVLooseIsolationMVA3newDMwLT");
      m_tauIso36   = patTau->tauID("byLooseIsolationMVA3newDMwLT");
      m_tauIso37   = patTau->tauID("byMediumIsolationMVA3newDMwLT");
      m_tauIso38   = patTau->tauID("byTightIsolationMVA3newDMwLT");
      m_tauIso39   = patTau->tauID("byVTightIsolationMVA3newDMwLT");
      m_tauIso40   = patTau->tauID("byVVTightIsolationMVA3newDMwLT");
      m_tauAgainstElectron01 = patTau->tauID("againstElectronLoose");
      m_tauAgainstElectron02 = patTau->tauID("againstElectronMedium");
      m_tauAgainstElectron03 = patTau->tauID("againstElectronTight");
      m_tauAgainstElectron04 = patTau->tauID("againstElectronMVA5raw");
      m_tauAgainstElectron05 = patTau->tauID("againstElectronMVA5category");
      m_tauAgainstElectron06 = patTau->tauID("againstElectronVLooseMVA5");
      m_tauAgainstElectron07 = patTau->tauID("againstElectronLooseMVA5");
      m_tauAgainstElectron08 = patTau->tauID("againstElectronMediumMVA5");
      m_tauAgainstElectron09 = patTau->tauID("againstElectronTightMVA5");
      m_tauAgainstElectron10 = patTau->tauID("againstElectronVTightMVA5");
      m_tauAgainstElectron11 = patTau->tauID("againstElectronDeadECAL");
      m_tauAgainstMuon01 = patTau->tauID("againstMuonLoose");
      m_tauAgainstMuon02 = patTau->tauID("againstMuonMedium");
      m_tauAgainstMuon03 = patTau->tauID("againstMuonTight");
      m_tauAgainstMuon04 = patTau->tauID("againstMuonLoose2");
      m_tauAgainstMuon05 = patTau->tauID("againstMuonMedium2");
      m_tauAgainstMuon06 = patTau->tauID("againstMuonTight2");
      m_tauAgainstMuon07 = patTau->tauID("againstMuonLoose3");
      m_tauAgainstMuon08 = patTau->tauID("againstMuonTight3");
      m_tauAgainstMuon09 = patTau->tauID("againstMuonMVAraw");
      m_tauAgainstMuon10 = patTau->tauID("againstMuonLooseMVA");
      m_tauAgainstMuon11 = patTau->tauID("againstMuonMediumMVA");
      m_tauAgainstMuon12 = patTau->tauID("againstMuonTightMVA");
      m_tauPuCorrPtSum   = patTau->tauID("puCorrPtSum");
      TreeQCDMuoTau->Fill();
    }

    //QCD ELETAU COLLECTION
    for (pat::TauCollection::const_iterator patTau = tauElTauHandle->begin(); patTau != tauElTauHandle->end(); ++patTau ) {
      float deltaRGen=99.; bool genJetMatch = false;
      for(vector<reco::GenJet>::const_iterator genjet=genjets->begin(); genjet!=genjets->end(); genjet++){
	if(ROOT::Math::VectorUtil::DeltaR(genjet->p4(),patTau->p4Jet())<0.3 && ROOT::Math::VectorUtil::DeltaR(genjet->p4(),patTau->p4Jet())<deltaRGen){
	  genJetMatch = true;
	  deltaRGen = ROOT::Math::VectorUtil::DeltaR(genjet->p4(),patTau->p4Jet());
	}
      }
      
      float taudzvertex = -9999;
      if(vertices.isValid() && vertices->size())taudzvertex = TMath::Abs(patTau->vertex().z()-(*vertices)[0].position().z());

      //if(!(patTau->pt()>20))    continue;
      if(!(realTau==false))     continue;
      if(!(genJetMatch==true))  continue;
      m_deltaR=deltaRGen;
      m_trigger=(int)isFired_HLT;
      m_trigger320=(int)isFired_HLT_PFJet320;
      m_trigger650=(int)isFired_HLT_HT650;
      m_tauPt=patTau->p4Jet().pt();
      m_tauEta=patTau->p4Jet().eta();
      m_tauVertexZ=taudzvertex;
      m_tauDecay01 = patTau->tauID("decayModeFinding");
      m_tauDecay02 = patTau->tauID("decayModeFindingNewDMs");
      m_tauDecay03 = patTau->tauID("decayModeFindingOldDMs");
      m_tauIso01   = patTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr");
      m_tauIso02   = patTau->tauID("byLooseIsolation");
      m_tauIso03   = patTau->tauID("byLooseCombinedIsolationDeltaBetaCorr");
      m_tauIso04   = patTau->tauID("byMediumCombinedIsolationDeltaBetaCorr");
      m_tauIso05   = patTau->tauID("byTightCombinedIsolationDeltaBetaCorr");
      m_tauIso06   = patTau->tauID("byCombinedIsolationDeltaBetaCorrRaw");
      m_tauIso07   = patTau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
      m_tauIso08   = patTau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
      m_tauIso09   = patTau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
      m_tauIso10   = patTau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
      m_tauIso11   = patTau->tauID("chargedIsoPtSum");
      m_tauIso12   = patTau->tauID("neutralIsoPtSum");
      m_tauIso13   = patTau->tauID("byIsolationMVA3oldDMwoLTraw");
      m_tauIso14   = patTau->tauID("byVLooseIsolationMVA3oldDMwoLT");
      m_tauIso15   = patTau->tauID("byLooseIsolationMVA3oldDMwoLT");
      m_tauIso16   = patTau->tauID("byMediumIsolationMVA3oldDMwoLT");
      m_tauIso17   = patTau->tauID("byTightIsolationMVA3oldDMwoLT");
      m_tauIso18   = patTau->tauID("byVTightIsolationMVA3oldDMwoLT");
      m_tauIso19   = patTau->tauID("byVVTightIsolationMVA3oldDMwoLT");
      m_tauIso20   = patTau->tauID("byIsolationMVA3oldDMwLTraw");
      m_tauIso21   = patTau->tauID("byVLooseIsolationMVA3oldDMwLT");
      m_tauIso22   = patTau->tauID("byLooseIsolationMVA3oldDMwLT");
      m_tauIso23   = patTau->tauID("byMediumIsolationMVA3oldDMwLT");
      m_tauIso24   = patTau->tauID("byTightIsolationMVA3oldDMwLT");
      m_tauIso25   = patTau->tauID("byVTightIsolationMVA3oldDMwLT");
      m_tauIso26   = patTau->tauID("byVVTightIsolationMVA3oldDMwLT");
      m_tauIso27   = patTau->tauID("byIsolationMVA3newDMwoLTraw");
      m_tauIso28   = patTau->tauID("byVLooseIsolationMVA3newDMwoLT");
      m_tauIso29   = patTau->tauID("byLooseIsolationMVA3newDMwoLT");
      m_tauIso30   = patTau->tauID("byMediumIsolationMVA3newDMwoLT");
      m_tauIso31   = patTau->tauID("byTightIsolationMVA3newDMwoLT");
      m_tauIso32   = patTau->tauID("byVTightIsolationMVA3newDMwoLT");
      m_tauIso33   = patTau->tauID("byVVTightIsolationMVA3newDMwoLT");
      m_tauIso34   = patTau->tauID("byIsolationMVA3newDMwLTraw");
      m_tauIso35   = patTau->tauID("byVLooseIsolationMVA3newDMwLT");
      m_tauIso36   = patTau->tauID("byLooseIsolationMVA3newDMwLT");
      m_tauIso37   = patTau->tauID("byMediumIsolationMVA3newDMwLT");
      m_tauIso38   = patTau->tauID("byTightIsolationMVA3newDMwLT");
      m_tauIso39   = patTau->tauID("byVTightIsolationMVA3newDMwLT");
      m_tauIso40   = patTau->tauID("byVVTightIsolationMVA3newDMwLT");
      m_tauAgainstElectron01 = patTau->tauID("againstElectronLoose");
      m_tauAgainstElectron02 = patTau->tauID("againstElectronMedium");
      m_tauAgainstElectron03 = patTau->tauID("againstElectronTight");
      m_tauAgainstElectron04 = patTau->tauID("againstElectronMVA5raw");
      m_tauAgainstElectron05 = patTau->tauID("againstElectronMVA5category");
      m_tauAgainstElectron06 = patTau->tauID("againstElectronVLooseMVA5");
      m_tauAgainstElectron07 = patTau->tauID("againstElectronLooseMVA5");
      m_tauAgainstElectron08 = patTau->tauID("againstElectronMediumMVA5");
      m_tauAgainstElectron09 = patTau->tauID("againstElectronTightMVA5");
      m_tauAgainstElectron10 = patTau->tauID("againstElectronVTightMVA5");
      m_tauAgainstElectron11 = patTau->tauID("againstElectronDeadECAL");
      m_tauAgainstMuon01 = patTau->tauID("againstMuonLoose");
      m_tauAgainstMuon02 = patTau->tauID("againstMuonMedium");
      m_tauAgainstMuon03 = patTau->tauID("againstMuonTight");
      m_tauAgainstMuon04 = patTau->tauID("againstMuonLoose2");
      m_tauAgainstMuon05 = patTau->tauID("againstMuonMedium2");
      m_tauAgainstMuon06 = patTau->tauID("againstMuonTight2");
      m_tauAgainstMuon07 = patTau->tauID("againstMuonLoose3");
      m_tauAgainstMuon08 = patTau->tauID("againstMuonTight3");
      m_tauAgainstMuon09 = patTau->tauID("againstMuonMVAraw");
      m_tauAgainstMuon10 = patTau->tauID("againstMuonLooseMVA");
      m_tauAgainstMuon11 = patTau->tauID("againstMuonMediumMVA");
      m_tauAgainstMuon12 = patTau->tauID("againstMuonTightMVA");
      m_tauPuCorrPtSum   = patTau->tauID("puCorrPtSum");
      TreeQCDEleTau->Fill();
    }
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
TauMCAnalyzer::beginJob()
{
  
  Service<TFileService> fs;

  TreeEffMuoTau = fs->make<TTree>("TreeEffMuoTau", "TreeEffMuoTau");
  TreeEffMuoTau->Branch("uno", &m_uno, "uno/i"); 
  TreeEffMuoTau->Branch("due", &m_due, "due/i"); 
  TreeEffMuoTau->Branch("tre", &m_tre, "tre/i"); 
  TreeEffMuoTau->Branch("qua", &m_qua, "qua/i"); 
  TreeEffMuoTau->Branch("cin", &m_cin, "cin/i"); 
  TreeEffMuoTau->Branch("sei", &m_sei, "sei/i"); 
  TreeEffMuoTau->Branch("set", &m_set, "set/i"); 
  TreeEffUsualMuoTau = fs->make<TTree>("TreeEffUsualMuoTau", "TreeEffUsualMuoTau");
  TreeEffUsualMuoTau->Branch("uno", &m_unoU, "uno/i"); 
  TreeEffUsualMuoTau->Branch("due", &m_dueU, "due/i"); 
  TreeEffUsualMuoTau->Branch("tre", &m_treU, "tre/i"); 
  TreeEffUsualMuoTau->Branch("qua", &m_quaU, "qua/i"); 
  TreeEffUsualMuoTau->Branch("cin", &m_cinU, "cin/i"); 
  TreeEffUsualMuoTau->Branch("sei", &m_seiU, "sei/i"); 
  TreeEffUsualMuoTau->Branch("set", &m_setU, "set/i"); 
  TreeEffEleTau = fs->make<TTree>("TreeEffEleTau", "TreeEffEleTau");
  TreeEffEleTau->Branch("uno", &m_uno, "uno/i"); 
  TreeEffEleTau->Branch("due", &m_due, "due/i"); 
  TreeEffEleTau->Branch("tre", &m_tre, "tre/i"); 
  TreeEffEleTau->Branch("qua", &m_qua, "qua/i"); 
  TreeEffEleTau->Branch("cin", &m_cin, "cin/i"); 
  TreeEffEleTau->Branch("sei", &m_sei, "sei/i"); 
  TreeEffEleTau->Branch("set", &m_set, "set/i"); 
  TreeEffUsualEleTau = fs->make<TTree>("TreeEffUsualEleTau", "TreeEffUsualEleTau");
  TreeEffUsualEleTau->Branch("uno", &m_unoU, "uno/i"); 
  TreeEffUsualEleTau->Branch("due", &m_dueU, "due/i"); 
  TreeEffUsualEleTau->Branch("tre", &m_treU, "tre/i"); 
  TreeEffUsualEleTau->Branch("qua", &m_quaU, "qua/i"); 
  TreeEffUsualEleTau->Branch("cin", &m_cinU, "cin/i"); 
  TreeEffUsualEleTau->Branch("sei", &m_seiU, "sei/i"); 
  TreeEffUsualEleTau->Branch("set", &m_setU, "set/i"); 

  Tree = fs->make<TTree>("Tree", "Tree");
  Tree->Branch("deltaRTauTau", &m_deltaRTauTau, "deltaRTauTau/f");
 
  TreeMuoTau = fs->make<TTree>("TreeMuoTau", "TreeMuoTau");
  TreeMuoTau->Branch("MassSVFit", &m_MassSVFit, "MassSVFit/f");
  TreeMuoTau->Branch("ptSVFit", &m_ptSVFit, "ptSVFit/f");
  TreeMuoTau->Branch("MEtSigmaParlZ", &m_MEtSigmaParlZ, "MEtSigmaParlZ/f");
  TreeMuoTau->Branch("MEtSigmaPerpZ", &m_MEtSigmaPerpZ, "MEtSigmaPerpZ/f");
  TreeMuoTau->Branch("MEtPullParlZ", &m_MEtPullParlZ, "MEtPullParlZ/f");
  TreeMuoTau->Branch("MEtPullPerpZ", &m_MEtPullPerpZ, "MEtPullPerpZ/f");
  TreeMuoTau->Branch("metErrParl", &m_metErrParl, "metErrParl/f");
  TreeMuoTau->Branch("metErrPerp", &m_metErrPerp, "metErrPerp/f");
  TreeMuoTau->Branch("metGen", &m_metGen, "metGen/f");
  TreeMuoTau->Branch("met", &m_met, "met/f");
  TreeMuoTau->Branch("tauPtGen", &m_tauPtGen, "tauPtGen/f");
  TreeMuoTau->Branch("tauEtaGen", &m_tauEtaGen, "tauEtaGen/f");
  TreeMuoTau->Branch("tauMassGen", &m_tauMassGen, "tauMassGen/f");
  TreeMuoTau->Branch("tauMass", &m_tauMass, "tauMass/f");
  TreeMuoTau->Branch("deltaRReco", &m_deltaRReco, "deltaRReco/f");
  TreeMuoTau->Branch("deltaRGen", &m_deltaRGen, "deltaRGen/f");
  TreeMuoTau->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeMuoTau->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeMuoTau->Branch("trigger", &m_trigger, "trigger/i");
  TreeMuoTau->Branch("deltaR", &m_deltaR, "deltaR/f");
  TreeMuoTau->Branch("tauPt", &m_tauPt, "tauPt/f");
  TreeMuoTau->Branch("tauEta", &m_tauEta, "tauEta/f");
  TreeMuoTau->Branch("tauDecay01", &m_tauDecay01, "tauDecay01/f");
  TreeMuoTau->Branch("tauDecay02", &m_tauDecay02, "tauDecay02/f");
  TreeMuoTau->Branch("tauDecay03", &m_tauDecay03, "tauDecay03/f");
  TreeMuoTau->Branch("tauIso01", &m_tauIso01, "tauIso01/f");
  TreeMuoTau->Branch("tauIso02", &m_tauIso02, "tauIso02/f");
  TreeMuoTau->Branch("tauIso03", &m_tauIso03, "tauIso03/f");
  TreeMuoTau->Branch("tauIso04", &m_tauIso04, "tauIso04/f");
  TreeMuoTau->Branch("tauIso05", &m_tauIso05, "tauIso05/f");
  TreeMuoTau->Branch("tauIso06", &m_tauIso06, "tauIso06/f");
  TreeMuoTau->Branch("tauIso07", &m_tauIso07, "tauIso07/f");
  TreeMuoTau->Branch("tauIso08", &m_tauIso08, "tauIso08/f");
  TreeMuoTau->Branch("tauIso09", &m_tauIso09, "tauIso09/f");
  TreeMuoTau->Branch("tauIso10", &m_tauIso10, "tauIso10/f");
  TreeMuoTau->Branch("tauIso11", &m_tauIso11, "tauIso11/f");
  TreeMuoTau->Branch("tauIso12", &m_tauIso12, "tauIso12/f");
  TreeMuoTau->Branch("tauIso13", &m_tauIso13, "tauIso13/f");
  TreeMuoTau->Branch("tauIso14", &m_tauIso14, "tauIso14/f");
  TreeMuoTau->Branch("tauIso15", &m_tauIso15, "tauIso15/f");
  TreeMuoTau->Branch("tauIso16", &m_tauIso16, "tauIso16/f");
  TreeMuoTau->Branch("tauIso17", &m_tauIso17, "tauIso17/f");
  TreeMuoTau->Branch("tauIso18", &m_tauIso18, "tauIso18/f");
  TreeMuoTau->Branch("tauIso19", &m_tauIso19, "tauIso19/f");
  TreeMuoTau->Branch("tauIso20", &m_tauIso20, "tauIso20/f");
  TreeMuoTau->Branch("tauIso21", &m_tauIso21, "tauIso21/f");
  TreeMuoTau->Branch("tauIso22", &m_tauIso22, "tauIso22/f");
  TreeMuoTau->Branch("tauIso23", &m_tauIso23, "tauIso23/f");
  TreeMuoTau->Branch("tauIso24", &m_tauIso24, "tauIso24/f");
  TreeMuoTau->Branch("tauIso25", &m_tauIso25, "tauIso25/f");
  TreeMuoTau->Branch("tauIso26", &m_tauIso26, "tauIso26/f");
  TreeMuoTau->Branch("tauIso27", &m_tauIso27, "tauIso27/f");
  TreeMuoTau->Branch("tauIso28", &m_tauIso28, "tauIso28/f");
  TreeMuoTau->Branch("tauIso29", &m_tauIso29, "tauIso29/f");
  TreeMuoTau->Branch("tauIso30", &m_tauIso30, "tauIso30/f");
  TreeMuoTau->Branch("tauIso31", &m_tauIso31, "tauIso31/f");
  TreeMuoTau->Branch("tauIso32", &m_tauIso32, "tauIso32/f");
  TreeMuoTau->Branch("tauIso33", &m_tauIso33, "tauIso33/f");
  TreeMuoTau->Branch("tauIso34", &m_tauIso34, "tauIso34/f");
  TreeMuoTau->Branch("tauIso35", &m_tauIso35, "tauIso35/f");
  TreeMuoTau->Branch("tauIso36", &m_tauIso36, "tauIso36/f");
  TreeMuoTau->Branch("tauIso37", &m_tauIso37, "tauIso37/f");
  TreeMuoTau->Branch("tauIso38", &m_tauIso38, "tauIso38/f");
  TreeMuoTau->Branch("tauIso39", &m_tauIso39, "tauIso39/f");
  TreeMuoTau->Branch("tauIso40", &m_tauIso40, "tauIso40/f");
  TreeMuoTau->Branch("tauAgainstElectron01", &m_tauAgainstElectron01, "tauAgainstElectron01/f");
  TreeMuoTau->Branch("tauAgainstElectron02", &m_tauAgainstElectron02, "tauAgainstElectron02/f");
  TreeMuoTau->Branch("tauAgainstElectron03", &m_tauAgainstElectron03, "tauAgainstElectron03/f");
  TreeMuoTau->Branch("tauAgainstElectron04", &m_tauAgainstElectron04, "tauAgainstElectron04/f");
  TreeMuoTau->Branch("tauAgainstElectron05", &m_tauAgainstElectron05, "tauAgainstElectron05/f");
  TreeMuoTau->Branch("tauAgainstElectron06", &m_tauAgainstElectron06, "tauAgainstElectron06/f");
  TreeMuoTau->Branch("tauAgainstElectron07", &m_tauAgainstElectron07, "tauAgainstElectron07/f");
  TreeMuoTau->Branch("tauAgainstElectron08", &m_tauAgainstElectron08, "tauAgainstElectron08/f");
  TreeMuoTau->Branch("tauAgainstElectron09", &m_tauAgainstElectron09, "tauAgainstElectron09/f");
  TreeMuoTau->Branch("tauAgainstElectron10", &m_tauAgainstElectron10, "tauAgainstElectron10/f");
  TreeMuoTau->Branch("tauAgainstElectron11", &m_tauAgainstElectron11, "tauAgainstElectron11/f");
  TreeMuoTau->Branch("tauAgainstMuon01", &m_tauAgainstMuon01, "tauAgainstMuon01/f");
  TreeMuoTau->Branch("tauAgainstMuon02", &m_tauAgainstMuon02, "tauAgainstMuon02/f");
  TreeMuoTau->Branch("tauAgainstMuon03", &m_tauAgainstMuon03, "tauAgainstMuon03/f");
  TreeMuoTau->Branch("tauAgainstMuon04", &m_tauAgainstMuon04, "tauAgainstMuon04/f");
  TreeMuoTau->Branch("tauAgainstMuon05", &m_tauAgainstMuon05, "tauAgainstMuon05/f");
  TreeMuoTau->Branch("tauAgainstMuon06", &m_tauAgainstMuon06, "tauAgainstMuon06/f");
  TreeMuoTau->Branch("tauAgainstMuon07", &m_tauAgainstMuon07, "tauAgainstMuon07/f");
  TreeMuoTau->Branch("tauAgainstMuon08", &m_tauAgainstMuon08, "tauAgainstMuon08/f");
  TreeMuoTau->Branch("tauAgainstMuon09", &m_tauAgainstMuon09, "tauAgainstMuon09/f");
  TreeMuoTau->Branch("tauAgainstMuon10", &m_tauAgainstMuon10, "tauAgainstMuon10/f");
  TreeMuoTau->Branch("tauAgainstMuon11", &m_tauAgainstMuon11, "tauAgainstMuon11/f");
  TreeMuoTau->Branch("tauAgainstMuon12", &m_tauAgainstMuon12, "tauAgainstMuon12/f");
  TreeMuoTau->Branch("tauPuCorrPtSum", &m_tauPuCorrPtSum, "tauPuCorrPtSum/f");
 

  TreeEleTau = fs->make<TTree>("TreeEleTau", "TreeEleTau");
  TreeEleTau->Branch("MassSVFit", &m_MassSVFit, "MassSVFit/f");
  TreeEleTau->Branch("ptSVFit", &m_ptSVFit, "ptSVFit/f");
  TreeEleTau->Branch("MEtSigmaParlZ", &m_MEtSigmaParlZ, "MEtSigmaParlZ/f");
  TreeEleTau->Branch("MEtSigmaPerpZ", &m_MEtSigmaPerpZ, "MEtSigmaPerpZ/f");
  TreeEleTau->Branch("MEtPullParlZ", &m_MEtPullParlZ, "MEtPullParlZ/f");
  TreeEleTau->Branch("MEtPullPerpZ", &m_MEtPullPerpZ, "MEtPullPerpZ/f");
  TreeEleTau->Branch("metErrParl", &m_metErrParl, "metErrParl/f");
  TreeEleTau->Branch("metErrPerp", &m_metErrPerp, "metErrPerp/f");
  TreeEleTau->Branch("metGen", &m_metGen, "metGen/f");
  TreeEleTau->Branch("met", &m_met, "met/f");
  TreeEleTau->Branch("tauPtGen", &m_tauPtGen, "tauPtGen/f");
  TreeEleTau->Branch("tauEtaGen", &m_tauEtaGen, "tauEtaGen/f");
  TreeEleTau->Branch("tauMassGen", &m_tauMassGen, "tauMassGen/f");
  TreeEleTau->Branch("tauMass", &m_tauMass, "tauMass/f");
  TreeEleTau->Branch("deltaRReco", &m_deltaRReco, "deltaRReco/f");
  TreeEleTau->Branch("deltaRGen", &m_deltaRGen, "deltaRGen/f");
  TreeEleTau->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeEleTau->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeEleTau->Branch("trigger", &m_trigger, "trigger/i");
  TreeEleTau->Branch("deltaR", &m_deltaR, "deltaR/f");
  TreeEleTau->Branch("tauPt", &m_tauPt, "tauPt/f");
  TreeEleTau->Branch("tauEta", &m_tauEta, "tauEta/f");
  TreeEleTau->Branch("tauDecay01", &m_tauDecay01, "tauDecay01/f");
  TreeEleTau->Branch("tauDecay02", &m_tauDecay02, "tauDecay02/f");
  TreeEleTau->Branch("tauDecay03", &m_tauDecay03, "tauDecay03/f");
  TreeEleTau->Branch("tauIso01", &m_tauIso01, "tauIso01/f");
  TreeEleTau->Branch("tauIso02", &m_tauIso02, "tauIso02/f");
  TreeEleTau->Branch("tauIso03", &m_tauIso03, "tauIso03/f");
  TreeEleTau->Branch("tauIso04", &m_tauIso04, "tauIso04/f");
  TreeEleTau->Branch("tauIso05", &m_tauIso05, "tauIso05/f");
  TreeEleTau->Branch("tauIso06", &m_tauIso06, "tauIso06/f");
  TreeEleTau->Branch("tauIso07", &m_tauIso07, "tauIso07/f");
  TreeEleTau->Branch("tauIso08", &m_tauIso08, "tauIso08/f");
  TreeEleTau->Branch("tauIso09", &m_tauIso09, "tauIso09/f");
  TreeEleTau->Branch("tauIso10", &m_tauIso10, "tauIso10/f");
  TreeEleTau->Branch("tauIso11", &m_tauIso11, "tauIso11/f");
  TreeEleTau->Branch("tauIso12", &m_tauIso12, "tauIso12/f");
  TreeEleTau->Branch("tauIso13", &m_tauIso13, "tauIso13/f");
  TreeEleTau->Branch("tauIso14", &m_tauIso14, "tauIso14/f");
  TreeEleTau->Branch("tauIso15", &m_tauIso15, "tauIso15/f");
  TreeEleTau->Branch("tauIso16", &m_tauIso16, "tauIso16/f");
  TreeEleTau->Branch("tauIso17", &m_tauIso17, "tauIso17/f");
  TreeEleTau->Branch("tauIso18", &m_tauIso18, "tauIso18/f");
  TreeEleTau->Branch("tauIso19", &m_tauIso19, "tauIso19/f");
  TreeEleTau->Branch("tauIso20", &m_tauIso20, "tauIso20/f");
  TreeEleTau->Branch("tauIso21", &m_tauIso21, "tauIso21/f");
  TreeEleTau->Branch("tauIso22", &m_tauIso22, "tauIso22/f");
  TreeEleTau->Branch("tauIso23", &m_tauIso23, "tauIso23/f");
  TreeEleTau->Branch("tauIso24", &m_tauIso24, "tauIso24/f");
  TreeEleTau->Branch("tauIso25", &m_tauIso25, "tauIso25/f");
  TreeEleTau->Branch("tauIso26", &m_tauIso26, "tauIso26/f");
  TreeEleTau->Branch("tauIso27", &m_tauIso27, "tauIso27/f");
  TreeEleTau->Branch("tauIso28", &m_tauIso28, "tauIso28/f");
  TreeEleTau->Branch("tauIso29", &m_tauIso29, "tauIso29/f");
  TreeEleTau->Branch("tauIso30", &m_tauIso30, "tauIso30/f");
  TreeEleTau->Branch("tauIso31", &m_tauIso31, "tauIso31/f");
  TreeEleTau->Branch("tauIso32", &m_tauIso32, "tauIso32/f");
  TreeEleTau->Branch("tauIso33", &m_tauIso33, "tauIso33/f");
  TreeEleTau->Branch("tauIso34", &m_tauIso34, "tauIso34/f");
  TreeEleTau->Branch("tauIso35", &m_tauIso35, "tauIso35/f");
  TreeEleTau->Branch("tauIso36", &m_tauIso36, "tauIso36/f");
  TreeEleTau->Branch("tauIso37", &m_tauIso37, "tauIso37/f");
  TreeEleTau->Branch("tauIso38", &m_tauIso38, "tauIso38/f");
  TreeEleTau->Branch("tauIso39", &m_tauIso39, "tauIso39/f");
  TreeEleTau->Branch("tauIso40", &m_tauIso40, "tauIso40/f");
  TreeEleTau->Branch("tauAgainstElectron01", &m_tauAgainstElectron01, "tauAgainstElectron01/f");
  TreeEleTau->Branch("tauAgainstElectron02", &m_tauAgainstElectron02, "tauAgainstElectron02/f");
  TreeEleTau->Branch("tauAgainstElectron03", &m_tauAgainstElectron03, "tauAgainstElectron03/f");
  TreeEleTau->Branch("tauAgainstElectron04", &m_tauAgainstElectron04, "tauAgainstElectron04/f");
  TreeEleTau->Branch("tauAgainstElectron05", &m_tauAgainstElectron05, "tauAgainstElectron05/f");
  TreeEleTau->Branch("tauAgainstElectron06", &m_tauAgainstElectron06, "tauAgainstElectron06/f");
  TreeEleTau->Branch("tauAgainstElectron07", &m_tauAgainstElectron07, "tauAgainstElectron07/f");
  TreeEleTau->Branch("tauAgainstElectron08", &m_tauAgainstElectron08, "tauAgainstElectron08/f");
  TreeEleTau->Branch("tauAgainstElectron09", &m_tauAgainstElectron09, "tauAgainstElectron09/f");
  TreeEleTau->Branch("tauAgainstElectron10", &m_tauAgainstElectron10, "tauAgainstElectron10/f");
  TreeEleTau->Branch("tauAgainstElectron11", &m_tauAgainstElectron11, "tauAgainstElectron11/f");
  TreeEleTau->Branch("tauAgainstMuon01", &m_tauAgainstMuon01, "tauAgainstMuon01/f");
  TreeEleTau->Branch("tauAgainstMuon02", &m_tauAgainstMuon02, "tauAgainstMuon02/f");
  TreeEleTau->Branch("tauAgainstMuon03", &m_tauAgainstMuon03, "tauAgainstMuon03/f");
  TreeEleTau->Branch("tauAgainstMuon04", &m_tauAgainstMuon04, "tauAgainstMuon04/f");
  TreeEleTau->Branch("tauAgainstMuon05", &m_tauAgainstMuon05, "tauAgainstMuon05/f");
  TreeEleTau->Branch("tauAgainstMuon06", &m_tauAgainstMuon06, "tauAgainstMuon06/f");
  TreeEleTau->Branch("tauAgainstMuon07", &m_tauAgainstMuon07, "tauAgainstMuon07/f");
  TreeEleTau->Branch("tauAgainstMuon08", &m_tauAgainstMuon08, "tauAgainstMuon08/f");
  TreeEleTau->Branch("tauAgainstMuon09", &m_tauAgainstMuon09, "tauAgainstMuon09/f");
  TreeEleTau->Branch("tauAgainstMuon10", &m_tauAgainstMuon10, "tauAgainstMuon10/f");
  TreeEleTau->Branch("tauAgainstMuon11", &m_tauAgainstMuon11, "tauAgainstMuon11/f");
  TreeEleTau->Branch("tauAgainstMuon12", &m_tauAgainstMuon12, "tauAgainstMuon12/f");
  TreeEleTau->Branch("tauPuCorrPtSum", &m_tauPuCorrPtSum, "tauPuCorrPtSum/f");
 

  TreeQCDMuoTau = fs->make<TTree>("TreeQCDMuoTau", "TreeQCDMuoTau");
  TreeQCDMuoTau->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeQCDMuoTau->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeQCDMuoTau->Branch("trigger", &m_trigger, "trigger/i");
  TreeQCDMuoTau->Branch("deltaR", &m_deltaR, "deltaR/f");
  TreeQCDMuoTau->Branch("tauPt", &m_tauPt, "tauPt/f");
  TreeQCDMuoTau->Branch("tauEta", &m_tauEta, "tauEta/f");
  TreeQCDMuoTau->Branch("tauVertexZ", &m_tauVertexZ, "tauVertexZ/f");
  TreeQCDMuoTau->Branch("tauDecay01", &m_tauDecay01, "tauDecay01/f");
  TreeQCDMuoTau->Branch("tauDecay02", &m_tauDecay02, "tauDecay02/f");
  TreeQCDMuoTau->Branch("tauDecay03", &m_tauDecay03, "tauDecay03/f");
  TreeQCDMuoTau->Branch("tauIso01", &m_tauIso01, "tauIso01/f");
  TreeQCDMuoTau->Branch("tauIso02", &m_tauIso02, "tauIso02/f");
  TreeQCDMuoTau->Branch("tauIso03", &m_tauIso03, "tauIso03/f");
  TreeQCDMuoTau->Branch("tauIso04", &m_tauIso04, "tauIso04/f");
  TreeQCDMuoTau->Branch("tauIso05", &m_tauIso05, "tauIso05/f");
  TreeQCDMuoTau->Branch("tauIso06", &m_tauIso06, "tauIso06/f");
  TreeQCDMuoTau->Branch("tauIso07", &m_tauIso07, "tauIso07/f");
  TreeQCDMuoTau->Branch("tauIso08", &m_tauIso08, "tauIso08/f");
  TreeQCDMuoTau->Branch("tauIso09", &m_tauIso09, "tauIso09/f");
  TreeQCDMuoTau->Branch("tauIso10", &m_tauIso10, "tauIso10/f");
  TreeQCDMuoTau->Branch("tauIso11", &m_tauIso11, "tauIso11/f");
  TreeQCDMuoTau->Branch("tauIso12", &m_tauIso12, "tauIso12/f");
  TreeQCDMuoTau->Branch("tauIso13", &m_tauIso13, "tauIso13/f");
  TreeQCDMuoTau->Branch("tauIso14", &m_tauIso14, "tauIso14/f");
  TreeQCDMuoTau->Branch("tauIso15", &m_tauIso15, "tauIso15/f");
  TreeQCDMuoTau->Branch("tauIso16", &m_tauIso16, "tauIso16/f");
  TreeQCDMuoTau->Branch("tauIso17", &m_tauIso17, "tauIso17/f");
  TreeQCDMuoTau->Branch("tauIso18", &m_tauIso18, "tauIso18/f");
  TreeQCDMuoTau->Branch("tauIso19", &m_tauIso19, "tauIso19/f");
  TreeQCDMuoTau->Branch("tauIso20", &m_tauIso20, "tauIso20/f");
  TreeQCDMuoTau->Branch("tauIso21", &m_tauIso21, "tauIso21/f");
  TreeQCDMuoTau->Branch("tauIso22", &m_tauIso22, "tauIso22/f");
  TreeQCDMuoTau->Branch("tauIso23", &m_tauIso23, "tauIso23/f");
  TreeQCDMuoTau->Branch("tauIso24", &m_tauIso24, "tauIso24/f");
  TreeQCDMuoTau->Branch("tauIso25", &m_tauIso25, "tauIso25/f");
  TreeQCDMuoTau->Branch("tauIso26", &m_tauIso26, "tauIso26/f");
  TreeQCDMuoTau->Branch("tauIso27", &m_tauIso27, "tauIso27/f");
  TreeQCDMuoTau->Branch("tauIso28", &m_tauIso28, "tauIso28/f");
  TreeQCDMuoTau->Branch("tauIso29", &m_tauIso29, "tauIso29/f");
  TreeQCDMuoTau->Branch("tauIso30", &m_tauIso30, "tauIso30/f");
  TreeQCDMuoTau->Branch("tauIso31", &m_tauIso31, "tauIso31/f");
  TreeQCDMuoTau->Branch("tauIso32", &m_tauIso32, "tauIso32/f");
  TreeQCDMuoTau->Branch("tauIso33", &m_tauIso33, "tauIso33/f");
  TreeQCDMuoTau->Branch("tauIso34", &m_tauIso34, "tauIso34/f");
  TreeQCDMuoTau->Branch("tauIso35", &m_tauIso35, "tauIso35/f");
  TreeQCDMuoTau->Branch("tauIso36", &m_tauIso36, "tauIso36/f");
  TreeQCDMuoTau->Branch("tauIso37", &m_tauIso37, "tauIso37/f");
  TreeQCDMuoTau->Branch("tauIso38", &m_tauIso38, "tauIso38/f");
  TreeQCDMuoTau->Branch("tauIso39", &m_tauIso39, "tauIso39/f");
  TreeQCDMuoTau->Branch("tauIso40", &m_tauIso40, "tauIso40/f");
  TreeQCDMuoTau->Branch("tauAgainstElectron01", &m_tauAgainstElectron01, "tauAgainstElectron01/f");
  TreeQCDMuoTau->Branch("tauAgainstElectron02", &m_tauAgainstElectron02, "tauAgainstElectron02/f");
  TreeQCDMuoTau->Branch("tauAgainstElectron03", &m_tauAgainstElectron03, "tauAgainstElectron03/f");
  TreeQCDMuoTau->Branch("tauAgainstElectron04", &m_tauAgainstElectron04, "tauAgainstElectron04/f");
  TreeQCDMuoTau->Branch("tauAgainstElectron05", &m_tauAgainstElectron05, "tauAgainstElectron05/f");
  TreeQCDMuoTau->Branch("tauAgainstElectron06", &m_tauAgainstElectron06, "tauAgainstElectron06/f");
  TreeQCDMuoTau->Branch("tauAgainstElectron07", &m_tauAgainstElectron07, "tauAgainstElectron07/f");
  TreeQCDMuoTau->Branch("tauAgainstElectron08", &m_tauAgainstElectron08, "tauAgainstElectron08/f");
  TreeQCDMuoTau->Branch("tauAgainstElectron09", &m_tauAgainstElectron09, "tauAgainstElectron09/f");
  TreeQCDMuoTau->Branch("tauAgainstElectron10", &m_tauAgainstElectron10, "tauAgainstElectron10/f");
  TreeQCDMuoTau->Branch("tauAgainstElectron11", &m_tauAgainstElectron11, "tauAgainstElectron11/f");
  TreeQCDMuoTau->Branch("tauAgainstMuon01", &m_tauAgainstMuon01, "tauAgainstMuon01/f");
  TreeQCDMuoTau->Branch("tauAgainstMuon02", &m_tauAgainstMuon02, "tauAgainstMuon02/f");
  TreeQCDMuoTau->Branch("tauAgainstMuon03", &m_tauAgainstMuon03, "tauAgainstMuon03/f");
  TreeQCDMuoTau->Branch("tauAgainstMuon04", &m_tauAgainstMuon04, "tauAgainstMuon04/f");
  TreeQCDMuoTau->Branch("tauAgainstMuon05", &m_tauAgainstMuon05, "tauAgainstMuon05/f");
  TreeQCDMuoTau->Branch("tauAgainstMuon06", &m_tauAgainstMuon06, "tauAgainstMuon06/f");
  TreeQCDMuoTau->Branch("tauAgainstMuon07", &m_tauAgainstMuon07, "tauAgainstMuon07/f");
  TreeQCDMuoTau->Branch("tauAgainstMuon08", &m_tauAgainstMuon08, "tauAgainstMuon08/f");
  TreeQCDMuoTau->Branch("tauAgainstMuon09", &m_tauAgainstMuon09, "tauAgainstMuon09/f");
  TreeQCDMuoTau->Branch("tauAgainstMuon10", &m_tauAgainstMuon10, "tauAgainstMuon10/f");
  TreeQCDMuoTau->Branch("tauAgainstMuon11", &m_tauAgainstMuon11, "tauAgainstMuon11/f");
  TreeQCDMuoTau->Branch("tauAgainstMuon12", &m_tauAgainstMuon12, "tauAgainstMuon12/f");
  TreeQCDMuoTau->Branch("tauPuCorrPtSum", &m_tauPuCorrPtSum, "tauPuCorrPtSum/f");
 

  TreeQCDEleTau = fs->make<TTree>("TreeQCDEleTau", "TreeQCDEleTau");
  TreeQCDEleTau->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeQCDEleTau->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeQCDEleTau->Branch("trigger", &m_trigger, "trigger/i");
  TreeQCDEleTau->Branch("deltaR", &m_deltaR, "deltaR/f");
  TreeQCDEleTau->Branch("tauPt", &m_tauPt, "tauPt/f");
  TreeQCDEleTau->Branch("tauEta", &m_tauEta, "tauEta/f");
  TreeQCDEleTau->Branch("tauVertexZ", &m_tauVertexZ, "tauVertexZ/f");
  TreeQCDEleTau->Branch("tauDecay01", &m_tauDecay01, "tauDecay01/f");
  TreeQCDEleTau->Branch("tauDecay02", &m_tauDecay02, "tauDecay02/f");
  TreeQCDEleTau->Branch("tauDecay03", &m_tauDecay03, "tauDecay03/f");
  TreeQCDEleTau->Branch("tauIso01", &m_tauIso01, "tauIso01/f");
  TreeQCDEleTau->Branch("tauIso02", &m_tauIso02, "tauIso02/f");
  TreeQCDEleTau->Branch("tauIso03", &m_tauIso03, "tauIso03/f");
  TreeQCDEleTau->Branch("tauIso04", &m_tauIso04, "tauIso04/f");
  TreeQCDEleTau->Branch("tauIso05", &m_tauIso05, "tauIso05/f");
  TreeQCDEleTau->Branch("tauIso06", &m_tauIso06, "tauIso06/f");
  TreeQCDEleTau->Branch("tauIso07", &m_tauIso07, "tauIso07/f");
  TreeQCDEleTau->Branch("tauIso08", &m_tauIso08, "tauIso08/f");
  TreeQCDEleTau->Branch("tauIso09", &m_tauIso09, "tauIso09/f");
  TreeQCDEleTau->Branch("tauIso10", &m_tauIso10, "tauIso10/f");
  TreeQCDEleTau->Branch("tauIso11", &m_tauIso11, "tauIso11/f");
  TreeQCDEleTau->Branch("tauIso12", &m_tauIso12, "tauIso12/f");
  TreeQCDEleTau->Branch("tauIso13", &m_tauIso13, "tauIso13/f");
  TreeQCDEleTau->Branch("tauIso14", &m_tauIso14, "tauIso14/f");
  TreeQCDEleTau->Branch("tauIso15", &m_tauIso15, "tauIso15/f");
  TreeQCDEleTau->Branch("tauIso16", &m_tauIso16, "tauIso16/f");
  TreeQCDEleTau->Branch("tauIso17", &m_tauIso17, "tauIso17/f");
  TreeQCDEleTau->Branch("tauIso18", &m_tauIso18, "tauIso18/f");
  TreeQCDEleTau->Branch("tauIso19", &m_tauIso19, "tauIso19/f");
  TreeQCDEleTau->Branch("tauIso20", &m_tauIso20, "tauIso20/f");
  TreeQCDEleTau->Branch("tauIso21", &m_tauIso21, "tauIso21/f");
  TreeQCDEleTau->Branch("tauIso22", &m_tauIso22, "tauIso22/f");
  TreeQCDEleTau->Branch("tauIso23", &m_tauIso23, "tauIso23/f");
  TreeQCDEleTau->Branch("tauIso24", &m_tauIso24, "tauIso24/f");
  TreeQCDEleTau->Branch("tauIso25", &m_tauIso25, "tauIso25/f");
  TreeQCDEleTau->Branch("tauIso26", &m_tauIso26, "tauIso26/f");
  TreeQCDEleTau->Branch("tauIso27", &m_tauIso27, "tauIso27/f");
  TreeQCDEleTau->Branch("tauIso28", &m_tauIso28, "tauIso28/f");
  TreeQCDEleTau->Branch("tauIso29", &m_tauIso29, "tauIso29/f");
  TreeQCDEleTau->Branch("tauIso30", &m_tauIso30, "tauIso30/f");
  TreeQCDEleTau->Branch("tauIso31", &m_tauIso31, "tauIso31/f");
  TreeQCDEleTau->Branch("tauIso32", &m_tauIso32, "tauIso32/f");
  TreeQCDEleTau->Branch("tauIso33", &m_tauIso33, "tauIso33/f");
  TreeQCDEleTau->Branch("tauIso34", &m_tauIso34, "tauIso34/f");
  TreeQCDEleTau->Branch("tauIso35", &m_tauIso35, "tauIso35/f");
  TreeQCDEleTau->Branch("tauIso36", &m_tauIso36, "tauIso36/f");
  TreeQCDEleTau->Branch("tauIso37", &m_tauIso37, "tauIso37/f");
  TreeQCDEleTau->Branch("tauIso38", &m_tauIso38, "tauIso38/f");
  TreeQCDEleTau->Branch("tauIso39", &m_tauIso39, "tauIso39/f");
  TreeQCDEleTau->Branch("tauIso40", &m_tauIso40, "tauIso40/f");
  TreeQCDEleTau->Branch("tauAgainstElectron01", &m_tauAgainstElectron01, "tauAgainstElectron01/f");
  TreeQCDEleTau->Branch("tauAgainstElectron02", &m_tauAgainstElectron02, "tauAgainstElectron02/f");
  TreeQCDEleTau->Branch("tauAgainstElectron03", &m_tauAgainstElectron03, "tauAgainstElectron03/f");
  TreeQCDEleTau->Branch("tauAgainstElectron04", &m_tauAgainstElectron04, "tauAgainstElectron04/f");
  TreeQCDEleTau->Branch("tauAgainstElectron05", &m_tauAgainstElectron05, "tauAgainstElectron05/f");
  TreeQCDEleTau->Branch("tauAgainstElectron06", &m_tauAgainstElectron06, "tauAgainstElectron06/f");
  TreeQCDEleTau->Branch("tauAgainstElectron07", &m_tauAgainstElectron07, "tauAgainstElectron07/f");
  TreeQCDEleTau->Branch("tauAgainstElectron08", &m_tauAgainstElectron08, "tauAgainstElectron08/f");
  TreeQCDEleTau->Branch("tauAgainstElectron09", &m_tauAgainstElectron09, "tauAgainstElectron09/f");
  TreeQCDEleTau->Branch("tauAgainstElectron10", &m_tauAgainstElectron10, "tauAgainstElectron10/f");
  TreeQCDEleTau->Branch("tauAgainstElectron11", &m_tauAgainstElectron11, "tauAgainstElectron11/f");
  TreeQCDEleTau->Branch("tauAgainstMuon01", &m_tauAgainstMuon01, "tauAgainstMuon01/f");
  TreeQCDEleTau->Branch("tauAgainstMuon02", &m_tauAgainstMuon02, "tauAgainstMuon02/f");
  TreeQCDEleTau->Branch("tauAgainstMuon03", &m_tauAgainstMuon03, "tauAgainstMuon03/f");
  TreeQCDEleTau->Branch("tauAgainstMuon04", &m_tauAgainstMuon04, "tauAgainstMuon04/f");
  TreeQCDEleTau->Branch("tauAgainstMuon05", &m_tauAgainstMuon05, "tauAgainstMuon05/f");
  TreeQCDEleTau->Branch("tauAgainstMuon06", &m_tauAgainstMuon06, "tauAgainstMuon06/f");
  TreeQCDEleTau->Branch("tauAgainstMuon07", &m_tauAgainstMuon07, "tauAgainstMuon07/f");
  TreeQCDEleTau->Branch("tauAgainstMuon08", &m_tauAgainstMuon08, "tauAgainstMuon08/f");
  TreeQCDEleTau->Branch("tauAgainstMuon09", &m_tauAgainstMuon09, "tauAgainstMuon09/f");
  TreeQCDEleTau->Branch("tauAgainstMuon10", &m_tauAgainstMuon10, "tauAgainstMuon10/f");
  TreeQCDEleTau->Branch("tauAgainstMuon11", &m_tauAgainstMuon11, "tauAgainstMuon11/f");
  TreeQCDEleTau->Branch("tauAgainstMuon12", &m_tauAgainstMuon12, "tauAgainstMuon12/f");
  TreeQCDEleTau->Branch("tauPuCorrPtSum", &m_tauPuCorrPtSum, "tauPuCorrPtSum/f");
 

  TreeUsualMuoTau = fs->make<TTree>("TreeUsualMuoTau", "TreeUsualMuoTau");
  TreeUsualMuoTau->Branch("MassSVFit", &m_MassSVFit, "MassSVFit/f");
  TreeUsualMuoTau->Branch("ptSVFit", &m_ptSVFit, "ptSVFit/f");
  TreeUsualMuoTau->Branch("MEtSigmaParlZ", &m_MEtSigmaParlZ, "MEtSigmaParlZ/f");
  TreeUsualMuoTau->Branch("MEtSigmaPerpZ", &m_MEtSigmaPerpZ, "MEtSigmaPerpZ/f");
  TreeUsualMuoTau->Branch("MEtPullParlZ", &m_MEtPullParlZ, "MEtPullParlZ/f");
  TreeUsualMuoTau->Branch("MEtPullPerpZ", &m_MEtPullPerpZ, "MEtPullPerpZ/f");
  TreeUsualMuoTau->Branch("metErrParl", &m_metErrParl, "metErrParl/f");
  TreeUsualMuoTau->Branch("metErrPerp", &m_metErrPerp, "metErrPerp/f");
  TreeUsualMuoTau->Branch("metGen", &m_metGen, "metGen/f");
  TreeUsualMuoTau->Branch("met", &m_met, "met/f");
  TreeUsualMuoTau->Branch("tauPtGen", &m_tauPtGen, "tauPtGen/f");
  TreeUsualMuoTau->Branch("tauEtaGen", &m_tauEtaGen, "tauEtaGen/f");
  TreeUsualMuoTau->Branch("tauMassGen", &m_tauMassGen, "tauMassGen/f");
  TreeUsualMuoTau->Branch("tauMass", &m_tauMass, "tauMass/f");
  TreeUsualMuoTau->Branch("deltaRReco", &m_deltaRReco, "deltaRReco/f");
  TreeUsualMuoTau->Branch("deltaRGen", &m_deltaRGen, "deltaRGen/f");
  TreeUsualMuoTau->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeUsualMuoTau->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeUsualMuoTau->Branch("trigger", &m_trigger, "trigger/i");
  TreeUsualMuoTau->Branch("deltaR", &m_deltaR, "deltaR/f");
  TreeUsualMuoTau->Branch("tauPt", &m_tauPt, "tauPt/f");
  TreeUsualMuoTau->Branch("tauEta", &m_tauEta, "tauEta/f");
  TreeUsualMuoTau->Branch("tauDecay01", &m_tauDecay01, "tauDecay01/f");
  TreeUsualMuoTau->Branch("tauDecay02", &m_tauDecay02, "tauDecay02/f");
  TreeUsualMuoTau->Branch("tauDecay03", &m_tauDecay03, "tauDecay03/f");
  TreeUsualMuoTau->Branch("tauIso01", &m_tauIso01, "tauIso01/f");
  TreeUsualMuoTau->Branch("tauIso02", &m_tauIso02, "tauIso02/f");
  TreeUsualMuoTau->Branch("tauIso03", &m_tauIso03, "tauIso03/f");
  TreeUsualMuoTau->Branch("tauIso04", &m_tauIso04, "tauIso04/f");
  TreeUsualMuoTau->Branch("tauIso05", &m_tauIso05, "tauIso05/f");
  TreeUsualMuoTau->Branch("tauIso06", &m_tauIso06, "tauIso06/f");
  TreeUsualMuoTau->Branch("tauIso07", &m_tauIso07, "tauIso07/f");
  TreeUsualMuoTau->Branch("tauIso08", &m_tauIso08, "tauIso08/f");
  TreeUsualMuoTau->Branch("tauIso09", &m_tauIso09, "tauIso09/f");
  TreeUsualMuoTau->Branch("tauIso10", &m_tauIso10, "tauIso10/f");
  TreeUsualMuoTau->Branch("tauIso11", &m_tauIso11, "tauIso11/f");
  TreeUsualMuoTau->Branch("tauIso12", &m_tauIso12, "tauIso12/f");
  TreeUsualMuoTau->Branch("tauIso13", &m_tauIso13, "tauIso13/f");
  TreeUsualMuoTau->Branch("tauIso14", &m_tauIso14, "tauIso14/f");
  TreeUsualMuoTau->Branch("tauIso15", &m_tauIso15, "tauIso15/f");
  TreeUsualMuoTau->Branch("tauIso16", &m_tauIso16, "tauIso16/f");
  TreeUsualMuoTau->Branch("tauIso17", &m_tauIso17, "tauIso17/f");
  TreeUsualMuoTau->Branch("tauIso18", &m_tauIso18, "tauIso18/f");
  TreeUsualMuoTau->Branch("tauIso19", &m_tauIso19, "tauIso19/f");
  TreeUsualMuoTau->Branch("tauIso20", &m_tauIso20, "tauIso20/f");
  TreeUsualMuoTau->Branch("tauIso21", &m_tauIso21, "tauIso21/f");
  TreeUsualMuoTau->Branch("tauIso22", &m_tauIso22, "tauIso22/f");
  TreeUsualMuoTau->Branch("tauIso23", &m_tauIso23, "tauIso23/f");
  TreeUsualMuoTau->Branch("tauIso24", &m_tauIso24, "tauIso24/f");
  TreeUsualMuoTau->Branch("tauIso25", &m_tauIso25, "tauIso25/f");
  TreeUsualMuoTau->Branch("tauIso26", &m_tauIso26, "tauIso26/f");
  TreeUsualMuoTau->Branch("tauIso27", &m_tauIso27, "tauIso27/f");
  TreeUsualMuoTau->Branch("tauIso28", &m_tauIso28, "tauIso28/f");
  TreeUsualMuoTau->Branch("tauIso29", &m_tauIso29, "tauIso29/f");
  TreeUsualMuoTau->Branch("tauIso30", &m_tauIso30, "tauIso30/f");
  TreeUsualMuoTau->Branch("tauIso31", &m_tauIso31, "tauIso31/f");
  TreeUsualMuoTau->Branch("tauIso32", &m_tauIso32, "tauIso32/f");
  TreeUsualMuoTau->Branch("tauIso33", &m_tauIso33, "tauIso33/f");
  TreeUsualMuoTau->Branch("tauIso34", &m_tauIso34, "tauIso34/f");
  TreeUsualMuoTau->Branch("tauIso35", &m_tauIso35, "tauIso35/f");
  TreeUsualMuoTau->Branch("tauIso36", &m_tauIso36, "tauIso36/f");
  TreeUsualMuoTau->Branch("tauIso37", &m_tauIso37, "tauIso37/f");
  TreeUsualMuoTau->Branch("tauIso38", &m_tauIso38, "tauIso38/f");
  TreeUsualMuoTau->Branch("tauIso39", &m_tauIso39, "tauIso39/f");
  TreeUsualMuoTau->Branch("tauIso40", &m_tauIso40, "tauIso40/f");
  TreeUsualMuoTau->Branch("tauAgainstElectron01", &m_tauAgainstElectron01, "tauAgainstElectron01/f");
  TreeUsualMuoTau->Branch("tauAgainstElectron02", &m_tauAgainstElectron02, "tauAgainstElectron02/f");
  TreeUsualMuoTau->Branch("tauAgainstElectron03", &m_tauAgainstElectron03, "tauAgainstElectron03/f");
  TreeUsualMuoTau->Branch("tauAgainstElectron04", &m_tauAgainstElectron04, "tauAgainstElectron04/f");
  TreeUsualMuoTau->Branch("tauAgainstElectron05", &m_tauAgainstElectron05, "tauAgainstElectron05/f");
  TreeUsualMuoTau->Branch("tauAgainstElectron06", &m_tauAgainstElectron06, "tauAgainstElectron06/f");
  TreeUsualMuoTau->Branch("tauAgainstElectron07", &m_tauAgainstElectron07, "tauAgainstElectron07/f");
  TreeUsualMuoTau->Branch("tauAgainstElectron08", &m_tauAgainstElectron08, "tauAgainstElectron08/f");
  TreeUsualMuoTau->Branch("tauAgainstElectron09", &m_tauAgainstElectron09, "tauAgainstElectron09/f");
  TreeUsualMuoTau->Branch("tauAgainstElectron10", &m_tauAgainstElectron10, "tauAgainstElectron10/f");
  TreeUsualMuoTau->Branch("tauAgainstElectron11", &m_tauAgainstElectron11, "tauAgainstElectron11/f");
  TreeUsualMuoTau->Branch("tauAgainstMuon01", &m_tauAgainstMuon01, "tauAgainstMuon01/f");
  TreeUsualMuoTau->Branch("tauAgainstMuon02", &m_tauAgainstMuon02, "tauAgainstMuon02/f");
  TreeUsualMuoTau->Branch("tauAgainstMuon03", &m_tauAgainstMuon03, "tauAgainstMuon03/f");
  TreeUsualMuoTau->Branch("tauAgainstMuon04", &m_tauAgainstMuon04, "tauAgainstMuon04/f");
  TreeUsualMuoTau->Branch("tauAgainstMuon05", &m_tauAgainstMuon05, "tauAgainstMuon05/f");
  TreeUsualMuoTau->Branch("tauAgainstMuon06", &m_tauAgainstMuon06, "tauAgainstMuon06/f");
  TreeUsualMuoTau->Branch("tauAgainstMuon07", &m_tauAgainstMuon07, "tauAgainstMuon07/f");
  TreeUsualMuoTau->Branch("tauAgainstMuon08", &m_tauAgainstMuon08, "tauAgainstMuon08/f");
  TreeUsualMuoTau->Branch("tauAgainstMuon09", &m_tauAgainstMuon09, "tauAgainstMuon09/f");
  TreeUsualMuoTau->Branch("tauAgainstMuon10", &m_tauAgainstMuon10, "tauAgainstMuon10/f");
  TreeUsualMuoTau->Branch("tauAgainstMuon11", &m_tauAgainstMuon11, "tauAgainstMuon11/f");
  TreeUsualMuoTau->Branch("tauAgainstMuon12", &m_tauAgainstMuon12, "tauAgainstMuon12/f");
  TreeUsualMuoTau->Branch("tauPuCorrPtSum", &m_tauPuCorrPtSum, "tauPuCorrPtSum/f");
 

  TreeUsualEleTau = fs->make<TTree>("TreeUsualEleTau", "TreeUsualEleTau");
  TreeUsualEleTau->Branch("MassSVFit", &m_MassSVFit, "MassSVFit/f");
  TreeUsualEleTau->Branch("ptSVFit", &m_ptSVFit, "ptSVFit/f");
  TreeUsualEleTau->Branch("MEtSigmaParlZ", &m_MEtSigmaParlZ, "MEtSigmaParlZ/f");
  TreeUsualEleTau->Branch("MEtSigmaPerpZ", &m_MEtSigmaPerpZ, "MEtSigmaPerpZ/f");
  TreeUsualEleTau->Branch("MEtPullParlZ", &m_MEtPullParlZ, "MEtPullParlZ/f");
  TreeUsualEleTau->Branch("MEtPullPerpZ", &m_MEtPullPerpZ, "MEtPullPerpZ/f");
  TreeUsualEleTau->Branch("metErrParl", &m_metErrParl, "metErrParl/f");
  TreeUsualEleTau->Branch("metErrPerp", &m_metErrPerp, "metErrPerp/f");
  TreeUsualEleTau->Branch("metGen", &m_metGen, "metGen/f");
  TreeUsualEleTau->Branch("met", &m_met, "met/f");
  TreeUsualEleTau->Branch("tauPtGen", &m_tauPtGen, "tauPtGen/f");
  TreeUsualEleTau->Branch("tauEtaGen", &m_tauEtaGen, "tauEtaGen/f");
  TreeUsualEleTau->Branch("tauMassGen", &m_tauMassGen, "tauMassGen/f");
  TreeUsualEleTau->Branch("tauMass", &m_tauMass, "tauMass/f");
  TreeUsualEleTau->Branch("deltaRReco", &m_deltaRReco, "deltaRReco/f");
  TreeUsualEleTau->Branch("deltaRGen", &m_deltaRGen, "deltaRGen/f");
  TreeUsualEleTau->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeUsualEleTau->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeUsualEleTau->Branch("trigger", &m_trigger, "trigger/i");
  TreeUsualEleTau->Branch("deltaR", &m_deltaR, "deltaR/f");
  TreeUsualEleTau->Branch("tauPt", &m_tauPt, "tauPt/f");
  TreeUsualEleTau->Branch("tauEta", &m_tauEta, "tauEta/f");
  TreeUsualEleTau->Branch("tauDecay01", &m_tauDecay01, "tauDecay01/f");
  TreeUsualEleTau->Branch("tauDecay02", &m_tauDecay02, "tauDecay02/f");
  TreeUsualEleTau->Branch("tauDecay03", &m_tauDecay03, "tauDecay03/f");
  TreeUsualEleTau->Branch("tauIso01", &m_tauIso01, "tauIso01/f");
  TreeUsualEleTau->Branch("tauIso02", &m_tauIso02, "tauIso02/f");
  TreeUsualEleTau->Branch("tauIso03", &m_tauIso03, "tauIso03/f");
  TreeUsualEleTau->Branch("tauIso04", &m_tauIso04, "tauIso04/f");
  TreeUsualEleTau->Branch("tauIso05", &m_tauIso05, "tauIso05/f");
  TreeUsualEleTau->Branch("tauIso06", &m_tauIso06, "tauIso06/f");
  TreeUsualEleTau->Branch("tauIso07", &m_tauIso07, "tauIso07/f");
  TreeUsualEleTau->Branch("tauIso08", &m_tauIso08, "tauIso08/f");
  TreeUsualEleTau->Branch("tauIso09", &m_tauIso09, "tauIso09/f");
  TreeUsualEleTau->Branch("tauIso10", &m_tauIso10, "tauIso10/f");
  TreeUsualEleTau->Branch("tauIso11", &m_tauIso11, "tauIso11/f");
  TreeUsualEleTau->Branch("tauIso12", &m_tauIso12, "tauIso12/f");
  TreeUsualEleTau->Branch("tauIso13", &m_tauIso13, "tauIso13/f");
  TreeUsualEleTau->Branch("tauIso14", &m_tauIso14, "tauIso14/f");
  TreeUsualEleTau->Branch("tauIso15", &m_tauIso15, "tauIso15/f");
  TreeUsualEleTau->Branch("tauIso16", &m_tauIso16, "tauIso16/f");
  TreeUsualEleTau->Branch("tauIso17", &m_tauIso17, "tauIso17/f");
  TreeUsualEleTau->Branch("tauIso18", &m_tauIso18, "tauIso18/f");
  TreeUsualEleTau->Branch("tauIso19", &m_tauIso19, "tauIso19/f");
  TreeUsualEleTau->Branch("tauIso20", &m_tauIso20, "tauIso20/f");
  TreeUsualEleTau->Branch("tauIso21", &m_tauIso21, "tauIso21/f");
  TreeUsualEleTau->Branch("tauIso22", &m_tauIso22, "tauIso22/f");
  TreeUsualEleTau->Branch("tauIso23", &m_tauIso23, "tauIso23/f");
  TreeUsualEleTau->Branch("tauIso24", &m_tauIso24, "tauIso24/f");
  TreeUsualEleTau->Branch("tauIso25", &m_tauIso25, "tauIso25/f");
  TreeUsualEleTau->Branch("tauIso26", &m_tauIso26, "tauIso26/f");
  TreeUsualEleTau->Branch("tauIso27", &m_tauIso27, "tauIso27/f");
  TreeUsualEleTau->Branch("tauIso28", &m_tauIso28, "tauIso28/f");
  TreeUsualEleTau->Branch("tauIso29", &m_tauIso29, "tauIso29/f");
  TreeUsualEleTau->Branch("tauIso30", &m_tauIso30, "tauIso30/f");
  TreeUsualEleTau->Branch("tauIso31", &m_tauIso31, "tauIso31/f");
  TreeUsualEleTau->Branch("tauIso32", &m_tauIso32, "tauIso32/f");
  TreeUsualEleTau->Branch("tauIso33", &m_tauIso33, "tauIso33/f");
  TreeUsualEleTau->Branch("tauIso34", &m_tauIso34, "tauIso34/f");
  TreeUsualEleTau->Branch("tauIso35", &m_tauIso35, "tauIso35/f");
  TreeUsualEleTau->Branch("tauIso36", &m_tauIso36, "tauIso36/f");
  TreeUsualEleTau->Branch("tauIso37", &m_tauIso37, "tauIso37/f");
  TreeUsualEleTau->Branch("tauIso38", &m_tauIso38, "tauIso38/f");
  TreeUsualEleTau->Branch("tauIso39", &m_tauIso39, "tauIso39/f");
  TreeUsualEleTau->Branch("tauIso40", &m_tauIso40, "tauIso40/f");
  TreeUsualEleTau->Branch("tauAgainstElectron01", &m_tauAgainstElectron01, "tauAgainstElectron01/f");
  TreeUsualEleTau->Branch("tauAgainstElectron02", &m_tauAgainstElectron02, "tauAgainstElectron02/f");
  TreeUsualEleTau->Branch("tauAgainstElectron03", &m_tauAgainstElectron03, "tauAgainstElectron03/f");
  TreeUsualEleTau->Branch("tauAgainstElectron04", &m_tauAgainstElectron04, "tauAgainstElectron04/f");
  TreeUsualEleTau->Branch("tauAgainstElectron05", &m_tauAgainstElectron05, "tauAgainstElectron05/f");
  TreeUsualEleTau->Branch("tauAgainstElectron06", &m_tauAgainstElectron06, "tauAgainstElectron06/f");
  TreeUsualEleTau->Branch("tauAgainstElectron07", &m_tauAgainstElectron07, "tauAgainstElectron07/f");
  TreeUsualEleTau->Branch("tauAgainstElectron08", &m_tauAgainstElectron08, "tauAgainstElectron08/f");
  TreeUsualEleTau->Branch("tauAgainstElectron09", &m_tauAgainstElectron09, "tauAgainstElectron09/f");
  TreeUsualEleTau->Branch("tauAgainstElectron10", &m_tauAgainstElectron10, "tauAgainstElectron10/f");
  TreeUsualEleTau->Branch("tauAgainstElectron11", &m_tauAgainstElectron11, "tauAgainstElectron11/f");
  TreeUsualEleTau->Branch("tauAgainstMuon01", &m_tauAgainstMuon01, "tauAgainstMuon01/f");
  TreeUsualEleTau->Branch("tauAgainstMuon02", &m_tauAgainstMuon02, "tauAgainstMuon02/f");
  TreeUsualEleTau->Branch("tauAgainstMuon03", &m_tauAgainstMuon03, "tauAgainstMuon03/f");
  TreeUsualEleTau->Branch("tauAgainstMuon04", &m_tauAgainstMuon04, "tauAgainstMuon04/f");
  TreeUsualEleTau->Branch("tauAgainstMuon05", &m_tauAgainstMuon05, "tauAgainstMuon05/f");
  TreeUsualEleTau->Branch("tauAgainstMuon06", &m_tauAgainstMuon06, "tauAgainstMuon06/f");
  TreeUsualEleTau->Branch("tauAgainstMuon07", &m_tauAgainstMuon07, "tauAgainstMuon07/f");
  TreeUsualEleTau->Branch("tauAgainstMuon08", &m_tauAgainstMuon08, "tauAgainstMuon08/f");
  TreeUsualEleTau->Branch("tauAgainstMuon09", &m_tauAgainstMuon09, "tauAgainstMuon09/f");
  TreeUsualEleTau->Branch("tauAgainstMuon10", &m_tauAgainstMuon10, "tauAgainstMuon10/f");
  TreeUsualEleTau->Branch("tauAgainstMuon11", &m_tauAgainstMuon11, "tauAgainstMuon11/f");
  TreeUsualEleTau->Branch("tauAgainstMuon12", &m_tauAgainstMuon12, "tauAgainstMuon12/f");
  TreeUsualEleTau->Branch("tauPuCorrPtSum", &m_tauPuCorrPtSum, "tauPuCorrPtSum/f");
 

  TreeQCDUsual = fs->make<TTree>("TreeQCDUsual", "TreeQCDUsual");
  TreeQCDUsual->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeQCDUsual->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeQCDUsual->Branch("trigger", &m_trigger, "trigger/i");
  TreeQCDUsual->Branch("deltaR", &m_deltaR, "deltaR/f");
  TreeQCDUsual->Branch("tauPt", &m_tauPt, "tauPt/f");
  TreeQCDUsual->Branch("tauEta", &m_tauEta, "tauEta/f");
  TreeQCDUsual->Branch("tauVertexZ", &m_tauVertexZ, "tauVertexZ/f");
  TreeQCDUsual->Branch("tauDecay01", &m_tauDecay01, "tauDecay01/f");
  TreeQCDUsual->Branch("tauDecay02", &m_tauDecay02, "tauDecay02/f");
  TreeQCDUsual->Branch("tauDecay03", &m_tauDecay03, "tauDecay03/f");
  TreeQCDUsual->Branch("tauIso01", &m_tauIso01, "tauIso01/f");
  TreeQCDUsual->Branch("tauIso02", &m_tauIso02, "tauIso02/f");
  TreeQCDUsual->Branch("tauIso03", &m_tauIso03, "tauIso03/f");
  TreeQCDUsual->Branch("tauIso04", &m_tauIso04, "tauIso04/f");
  TreeQCDUsual->Branch("tauIso05", &m_tauIso05, "tauIso05/f");
  TreeQCDUsual->Branch("tauIso06", &m_tauIso06, "tauIso06/f");
  TreeQCDUsual->Branch("tauIso07", &m_tauIso07, "tauIso07/f");
  TreeQCDUsual->Branch("tauIso08", &m_tauIso08, "tauIso08/f");
  TreeQCDUsual->Branch("tauIso09", &m_tauIso09, "tauIso09/f");
  TreeQCDUsual->Branch("tauIso10", &m_tauIso10, "tauIso10/f");
  TreeQCDUsual->Branch("tauIso11", &m_tauIso11, "tauIso11/f");
  TreeQCDUsual->Branch("tauIso12", &m_tauIso12, "tauIso12/f");
  TreeQCDUsual->Branch("tauIso13", &m_tauIso13, "tauIso13/f");
  TreeQCDUsual->Branch("tauIso14", &m_tauIso14, "tauIso14/f");
  TreeQCDUsual->Branch("tauIso15", &m_tauIso15, "tauIso15/f");
  TreeQCDUsual->Branch("tauIso16", &m_tauIso16, "tauIso16/f");
  TreeQCDUsual->Branch("tauIso17", &m_tauIso17, "tauIso17/f");
  TreeQCDUsual->Branch("tauIso18", &m_tauIso18, "tauIso18/f");
  TreeQCDUsual->Branch("tauIso19", &m_tauIso19, "tauIso19/f");
  TreeQCDUsual->Branch("tauIso20", &m_tauIso20, "tauIso20/f");
  TreeQCDUsual->Branch("tauIso21", &m_tauIso21, "tauIso21/f");
  TreeQCDUsual->Branch("tauIso22", &m_tauIso22, "tauIso22/f");
  TreeQCDUsual->Branch("tauIso23", &m_tauIso23, "tauIso23/f");
  TreeQCDUsual->Branch("tauIso24", &m_tauIso24, "tauIso24/f");
  TreeQCDUsual->Branch("tauIso25", &m_tauIso25, "tauIso25/f");
  TreeQCDUsual->Branch("tauIso26", &m_tauIso26, "tauIso26/f");
  TreeQCDUsual->Branch("tauIso27", &m_tauIso27, "tauIso27/f");
  TreeQCDUsual->Branch("tauIso28", &m_tauIso28, "tauIso28/f");
  TreeQCDUsual->Branch("tauIso29", &m_tauIso29, "tauIso29/f");
  TreeQCDUsual->Branch("tauIso30", &m_tauIso30, "tauIso30/f");
  TreeQCDUsual->Branch("tauIso31", &m_tauIso31, "tauIso31/f");
  TreeQCDUsual->Branch("tauIso32", &m_tauIso32, "tauIso32/f");
  TreeQCDUsual->Branch("tauIso33", &m_tauIso33, "tauIso33/f");
  TreeQCDUsual->Branch("tauIso34", &m_tauIso34, "tauIso34/f");
  TreeQCDUsual->Branch("tauIso35", &m_tauIso35, "tauIso35/f");
  TreeQCDUsual->Branch("tauIso36", &m_tauIso36, "tauIso36/f");
  TreeQCDUsual->Branch("tauIso37", &m_tauIso37, "tauIso37/f");
  TreeQCDUsual->Branch("tauIso38", &m_tauIso38, "tauIso38/f");
  TreeQCDUsual->Branch("tauIso39", &m_tauIso39, "tauIso39/f");
  TreeQCDUsual->Branch("tauIso40", &m_tauIso40, "tauIso40/f");
  TreeQCDUsual->Branch("tauAgainstElectron01", &m_tauAgainstElectron01, "tauAgainstElectron01/f");
  TreeQCDUsual->Branch("tauAgainstElectron02", &m_tauAgainstElectron02, "tauAgainstElectron02/f");
  TreeQCDUsual->Branch("tauAgainstElectron03", &m_tauAgainstElectron03, "tauAgainstElectron03/f");
  TreeQCDUsual->Branch("tauAgainstElectron04", &m_tauAgainstElectron04, "tauAgainstElectron04/f");
  TreeQCDUsual->Branch("tauAgainstElectron05", &m_tauAgainstElectron05, "tauAgainstElectron05/f");
  TreeQCDUsual->Branch("tauAgainstElectron06", &m_tauAgainstElectron06, "tauAgainstElectron06/f");
  TreeQCDUsual->Branch("tauAgainstElectron07", &m_tauAgainstElectron07, "tauAgainstElectron07/f");
  TreeQCDUsual->Branch("tauAgainstElectron08", &m_tauAgainstElectron08, "tauAgainstElectron08/f");
  TreeQCDUsual->Branch("tauAgainstElectron09", &m_tauAgainstElectron09, "tauAgainstElectron09/f");
  TreeQCDUsual->Branch("tauAgainstElectron10", &m_tauAgainstElectron10, "tauAgainstElectron10/f");
  TreeQCDUsual->Branch("tauAgainstElectron11", &m_tauAgainstElectron11, "tauAgainstElectron11/f");
  TreeQCDUsual->Branch("tauAgainstMuon01", &m_tauAgainstMuon01, "tauAgainstMuon01/f");
  TreeQCDUsual->Branch("tauAgainstMuon02", &m_tauAgainstMuon02, "tauAgainstMuon02/f");
  TreeQCDUsual->Branch("tauAgainstMuon03", &m_tauAgainstMuon03, "tauAgainstMuon03/f");
  TreeQCDUsual->Branch("tauAgainstMuon04", &m_tauAgainstMuon04, "tauAgainstMuon04/f");
  TreeQCDUsual->Branch("tauAgainstMuon05", &m_tauAgainstMuon05, "tauAgainstMuon05/f");
  TreeQCDUsual->Branch("tauAgainstMuon06", &m_tauAgainstMuon06, "tauAgainstMuon06/f");
  TreeQCDUsual->Branch("tauAgainstMuon07", &m_tauAgainstMuon07, "tauAgainstMuon07/f");
  TreeQCDUsual->Branch("tauAgainstMuon08", &m_tauAgainstMuon08, "tauAgainstMuon08/f");
  TreeQCDUsual->Branch("tauAgainstMuon09", &m_tauAgainstMuon09, "tauAgainstMuon09/f");
  TreeQCDUsual->Branch("tauAgainstMuon10", &m_tauAgainstMuon10, "tauAgainstMuon10/f");
  TreeQCDUsual->Branch("tauAgainstMuon11", &m_tauAgainstMuon11, "tauAgainstMuon11/f");
  TreeQCDUsual->Branch("tauAgainstMuon12", &m_tauAgainstMuon12, "tauAgainstMuon12/f");
  TreeQCDUsual->Branch("tauPuCorrPtSum", &m_tauPuCorrPtSum, "tauPuCorrPtSum/f");

  return;
 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TauMCAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
TauMCAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TauMCAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TauMCAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TauMCAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TauMCAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TauMCAnalyzer);
