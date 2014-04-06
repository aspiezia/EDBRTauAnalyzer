// -*- C++ -*-
//
// Package: Efficiency
// Class: Efficiency
//
/**\class Efficiency Efficiency.cc ExoDiBosonResonances/EDBRTauAnalyzer/src/Efficiency.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author: Aniello Spiezia,21 1-007,+41227676459,
// Created: Mon Sep 9 13:14:05 CEST 2013
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

//new inclusion
#include "TH1.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TauAnalysis/CandidateTools/interface/NSVfitStandaloneAlgorithm.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "TLorentzVector.h"
//
// class declaration
//

class Efficiency : public edm::EDAnalyzer {
public:
  explicit Efficiency(const edm::ParameterSet&);
  ~Efficiency();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  //FUNCTION
  void SelectJet(edm::Handle<pat::JetCollection> CA8JetswithQjets, edm::Handle<pat::JetCollection> CA8JetsPruned, bool & foundSelectedJet,
                 pat::JetCollection::const_iterator & SelectedJet, float & massZ, float & tau21Z, float & ptZ, float massMin, float massMax,
		 bool foundJet, std::vector<reco::GenJet>::const_iterator GenJeT, bool & matched);
  void SelectTau(edm::Handle<pat::TauCollection> tauHandle, pat::JetCollection::const_iterator SelectedJet, bool & foundTau, 
		 pat::TauCollection::const_iterator & SelectedTau, float & ptTau, bool foundJet,
		 bool foundPart, math::PtEtaPhiELorentzVector genPart, bool & matched);
  void SelectHighptMuon(edm::Handle<pat::MuonCollection> muoH, pat::JetCollection::const_iterator SelectedJet, bool & foundMuon,
			pat::MuonCollection::const_iterator & SelectedMuon, float & ptMuon, bool foundJet, reco::Vertex primaryVertex,
			bool foundMuo, math::PtEtaPhiELorentzVector genMuo_mt, bool & matched);
  void SelectElectron(edm::Handle<pat::ElectronCollection> eleH, pat::JetCollection::const_iterator SelectedJet, bool & foundElectron,
		      pat::ElectronCollection::const_iterator & SelectedElectron, float & ptElectron, bool foundJet, reco::Vertex primaryVertex,
		      bool foundEle_et, math::PtEtaPhiELorentzVector genEle_et, bool & matched);
  void SelectTrackerMuon(edm::Handle<pat::MuonCollection> muoH, bool & foundMuon, pat::MuonCollection::const_iterator & SelectedMuon, float & ptElectron, 
			 reco::Vertex primaryVertex, pat::JetCollection::const_iterator SelectedJet, 
			 bool foundJet, bool foundMuo1_mm, math::PtEtaPhiELorentzVector genMuo1_mm, bool & matched);
  void SelectTightMuon(edm::Handle<pat::MuonCollection> muoH, pat::JetCollection::const_iterator SelectedJet, bool & foundMuon,
		       pat::MuonCollection::const_iterator & SelectedMuon, float & ptMuon, bool foundJet, reco::Vertex primaryVertex,
		       bool foundMuo, math::PtEtaPhiELorentzVector genMuo_mt, bool & matched);
  void SelectTrackerGlobalID(pat::MuonCollection::const_iterator SelectedMuo1, pat::MuonCollection::const_iterator SelectedMuo2, 
			     bool matched1, bool matched2, reco::Vertex primaryVertex, bool & hasAtLeastOneHighPtMuo);
  
  float MuonPFIso(pat::MuonCollection::const_iterator muon, bool highpt);
  float ElectronPFIso(pat::ElectronCollection::const_iterator electron, float rho);
  float MuonCorrPFIso(pat::MuonCollection::const_iterator muon, bool highpt);
  float ElectronCorrPFIso(pat::ElectronCollection::const_iterator electron, float rho);
  float MuonDETIso(pat::MuonCollection::const_iterator SelectedMuo1, pat::MuonCollection::const_iterator SelectedMuo2, bool secondMuon);


  TTree *Tree; 
  TTree *TreeQCDMuon; 
  TTree *TreeQCDElectron; 
 
  float m_JetMuon_Pt;
  float m_JetMuon_Eta;
  float m_JetMuon_PFIso;
  float m_JetMuon_CorrPFIso;
  float m_JetElectron_Pt;
  float m_JetElectron_Eta;
  float m_JetElectron_PFIso;
  float m_JetElectron_CorrPFIso;

  int   m_Nele;
  int   m_Nmuo;
  int   m_genJet;
  int   m_genTau_et;
  int   m_genTau_mt;
  int   m_genEle_et;
  int   m_genEle1_ee;
  int   m_genEle2_ee;
  int   m_genEle_em;
  int   m_genMuo_mt;
  int   m_genMuo1_mm;
  int   m_genMuo2_mm;
  int   m_genMuo_em;
  int   m_matchedJet;
  int   m_matchedTau_et;
  int   m_matchedTau_mt;
  int   m_matchedEle_et;
  int   m_matchedEle1_ee;
  int   m_matchedEle2_ee;
  int   m_matchedEle_em;
  int   m_matchedTightMuo_mt;
  int   m_matchedTightMuo1_mm;
  int   m_matchedTightMuo2_mm;
  int   m_matchedTightMuo_em;
  int   m_matchedHighptMuo_mt;
  int   m_matchedHighptMuo1_mm;
  int   m_matchedHighptMuo2_mm;
  int   m_matchedHighptMuo_em;
  int   m_matchedTrackerMuo1_mm;
  int   m_matchedTrackerMuo2_mm;
  int   m_recoJet;
  int   m_recoTau_mt;
  int   m_recoTau_et;
  int   m_recoEle_et;
  int   m_recoEle1_ee;
  int   m_recoEle2_ee;
  int   m_recoEle_em;
  int   m_recoTightMuo_mt;
  int   m_recoTightMuo1_mm;
  int   m_recoTightMuo2_mm;
  int   m_recoTightMuo_em;
  int   m_recoHighptMuo_mt;
  int   m_recoHighptMuo1_mm;
  int   m_recoHighptMuo2_mm;
  int   m_recoHighptMuo_em;
  int   m_recoTrackerMuo1_mm;
  int   m_recoTrackerMuo2_mm;
  int   m_hasAtLeastOneHighPtMuo;
  float m_genJet_Pt;
  float m_genJet_Eta;
  float m_genTau_et_Pt;
  float m_genTau_et_Eta;
  float m_genTau_mt_Pt;
  float m_genTau_mt_eta;
  float m_genEle_et_Pt;
  float m_genEle_et_Eta;
  float m_genEle1_ee_Pt;
  float m_genEle1_ee_Eta;
  float m_genEle2_ee_Pt;
  float m_genEle2_ee_Eta;
  float m_genEle_em_Pt;
  float m_genEle_em_Eta;
  float m_genMuo_mt_Pt;
  float m_genMuo_mt_Eta;
  float m_genMuo1_mm_Pt;
  float m_genMuo1_mm_Eta;
  float m_genMuo2_mm_Pt;
  float m_genMuo2_mm_Eta;
  float m_genMuo_em_Pt;
  float m_genMuo_em_Eta;
  float m_recoJet_Pt;
  float m_recoJet_Eta;
  float m_recoTau_et_Pt;
  float m_recoTau_et_Eta;
  float m_recoTau_mt_Pt;
  float m_recoTau_mt_eta;
  float m_recoEle_et_Pt;
  float m_recoEle_et_Eta;
  float m_recoEle1_ee_Pt;
  float m_recoEle1_ee_Eta;
  float m_recoEle2_ee_Pt;
  float m_recoEle2_ee_Eta;
  float m_recoEle_em_Pt;
  float m_recoEle_em_Eta;
  float m_recoTightMuo_mt_Pt;
  float m_recoTightMuo_mt_Eta;
  float m_recoTightMuo1_mm_Pt;
  float m_recoTightMuo1_mm_Eta;
  float m_recoTightMuo2_mm_Pt;
  float m_recoTightMuo2_mm_Eta;
  float m_recoTightMuo_em_Pt;
  float m_recoTightMuo_em_Eta;
  float m_recoHighptMuo_mt_Pt;
  float m_recoHighptMuo_mt_Eta;
  float m_recoHighptMuo1_mm_Pt;
  float m_recoHighptMuo1_mm_Eta;
  float m_recoHighptMuo2_mm_Pt;
  float m_recoHighptMuo2_mm_Eta;
  float m_recoHighptMuo_em_Pt;
  float m_recoHighptMuo_em_Eta;
  float m_recoTrackerMuo1_mm_Pt;
  float m_recoTrackerMuo1_mm_Eta;
  float m_recoTrackerMuo2_mm_Pt;
  float m_recoTrackerMuo2_mm_Eta;
  float m_recoEle_et_PFIso;
  float m_recoEle_et_CorrPFIso;
  float m_recoEle1_ee_PFIso;
  float m_recoEle1_ee_CorrPFIso;
  float m_recoEle2_ee_PFIso;
  float m_recoEle2_ee_CorrPFIso;
  float m_recoEle_em_PFIso;
  float m_recoEle_em_CorrPFIso;
  float m_recoTightMuo_mt_PFIso;
  float m_recoTightMuo_mt_CorrPFIso;
  float m_recoTightMuo1_mm_PFIso;
  float m_recoTightMuo1_mm_CorrPFIso;
  float m_recoTightMuo1_mm_DetIso;
  float m_recoTightMuo2_mm_PFIso;
  float m_recoTightMuo2_mm_CorrPFIso;
  float m_recoTightMuo2_mm_DetIso;
  float m_recoTightMuo_em_PFIso;
  float m_recoTightMuo_em_CorrPFIso;
  float m_recoHighptMuo_mt_PFIso;
  float m_recoHighptMuo_mt_CorrPFIso;
  float m_recoHighptMuo1_mm_PFIso;
  float m_recoHighptMuo1_mm_CorrPFIso;
  float m_recoHighptMuo1_mm_DetIso;
  float m_recoHighptMuo2_mm_PFIso;
  float m_recoHighptMuo2_mm_CorrPFIso;
  float m_recoHighptMuo2_mm_DetIso;
  float m_recoHighptMuo_em_PFIso;
  float m_recoHighptMuo_em_CorrPFIso;
  float m_recoTrackerMuo1_mm_PFIso;
  float m_recoTrackerMuo1_mm_CorrPFIso;
  float m_recoTrackerMuo1_mm_DetIso;
  float m_recoTrackerMuo2_mm_PFIso;
  float m_recoTrackerMuo2_mm_CorrPFIso;
  float m_recoTrackerMuo2_mm_DetIso;
  float m_reco_em_deltaR;
  float m_reco_et_deltaR;
  float m_reco_ee_deltaR;
  float m_reco_mt_deltaR;
  float m_reco_mm_deltaR;
  float m_dRJet;
  float m_dRTau_et;
  float m_dRTau_mt;
  float m_dREle_et;
  float m_dREle1_ee;
  float m_dREle2_ee;
  float m_dREle_em;
  float m_dRMuo_mt;
  float m_dRMuo1_mm;
  float m_dRMuo2_mm;
  float m_dRMuo_em;

  edm::LumiReWeighting LumiWeights_;
  bool isData; 
  edm::InputTag vtxColl_;
  edm::InputTag jetColl_;
  edm::InputTag jetPrunedColl_;
  edm::InputTag metColl_;
  edm::InputTag electronColl_;
  edm::InputTag muonColl_;
  edm::InputTag tauMuTauColl_;
  edm::InputTag tauElTauColl_;
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
using NSVfitStandalone::Vector;
using NSVfitStandalone::LorentzVector;
using NSVfitStandalone::MeasuredTauLepton;


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Efficiency::Efficiency(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
  isData = iConfig.getUntrackedParameter<bool>("isData_");
  vtxColl_ = iConfig.getParameter<edm::InputTag>("vtxColl"); 
  jetColl_ = iConfig.getParameter<edm::InputTag>("jetColl"); 
  jetPrunedColl_ = iConfig.getParameter<edm::InputTag>("jetPrunedColl"); 
  metColl_ = iConfig.getParameter<edm::InputTag>("metColl"); 
  electronColl_ = iConfig.getParameter<edm::InputTag>("electronColl"); 
  muonColl_ = iConfig.getParameter<edm::InputTag>("muonColl"); 
  tauMuTauColl_ = iConfig.getParameter<edm::InputTag>("tauMuTauColl");
  tauElTauColl_ = iConfig.getParameter<edm::InputTag>("tauElTauColl"); 
  metRawColl_ = iConfig.getParameter<edm::InputTag>("metRawColl"); 
  uncorrmetColl_ = iConfig.getParameter<edm::InputTag>("uncorrmetColl"); 
  ak5JetColl_ = iConfig.getParameter<edm::InputTag>("ak5JetColl");
  NeventsTOT_ = iConfig.getParameter<int>( "NeventsTOT" );
  xsec_= iConfig.getParameter<double>( "xsec" );
  lumi_= iConfig.getParameter<double>( "lumi" );

}


Efficiency::~Efficiency()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event ------------
void
Efficiency::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(vtxColl_, vertices);
  reco::Vertex primaryVertex;
  primaryVertex = vertices->at(0);

  edm::Handle<pat::JetCollection> CA8JetswithQjets;
  iEvent.getByLabel(jetColl_, CA8JetswithQjets);
  edm::Handle<pat::JetCollection> CA8JetsPruned;
  iEvent.getByLabel(jetPrunedColl_, CA8JetsPruned);

  edm::Handle<pat::ElectronCollection> eleH;
  iEvent.getByLabel(electronColl_, eleH);

  edm::Handle<pat::MuonCollection> muoH;
  iEvent.getByLabel(muonColl_, muoH);

  edm::Handle<pat::METCollection> met;
  iEvent.getByLabel(metColl_, met);

  edm::Handle<pat::METCollection> metRaw;
  iEvent.getByLabel(metRawColl_, metRaw);

  edm::Handle<pat::METCollection> uncorrmet;
  iEvent.getByLabel(uncorrmetColl_, uncorrmet);

  edm::Handle<pat::JetCollection> ak5jetCands;
  iEvent.getByLabel(ak5JetColl_,ak5jetCands);

  Handle<pat::TauCollection> tauMuTauHandle;
  iEvent.getByLabel(tauMuTauColl_,tauMuTauHandle);

  Handle<pat::TauCollection> tauElTauHandle;
  iEvent.getByLabel(tauElTauColl_,tauElTauHandle);

  edm::Handle<double> rhoHandle;
  iEvent.getByLabel("kt6PFJets", "rho", rhoHandle);
  float rho = *(rhoHandle.product());

  Handle<vector<reco::GenParticle> > genParts;
  iEvent.getByLabel("genParticles", genParts);

  edm::Handle<vector<reco::GenJet> > genjets;
  iEvent.getByLabel("ak5GenJetsNoNu", genjets);

  
  //DEFINING THE GEN EVENT
  int ele = 0; int muo = 0;
  for(size_t ngenPart=0; ngenPart<genParts->size(); ngenPart++){
    const reco::GenParticle & genPart = (*genParts)[ngenPart];
    if(abs(genPart.pdgId())==2212) continue;
    const reco::Candidate * mom = genPart.mother();
    if(abs(genPart.pdgId())==15 && genPart.status()!=3 && (abs(mom->pdgId())==25 || abs(mom->pdgId())==15)){
      for(unsigned int ndaugh=0; ndaugh<genPart.numberOfDaughters(); ndaugh++){
	const reco::Candidate * daughter = genPart.daughter(ndaugh);
	if(abs(daughter->pdgId())==11 && daughter->status()==1) ele = ele + 1;
	if(abs(daughter->pdgId())==13 && daughter->status()==1) muo = muo + 1;
      }
    }
  }
  
  math::PtEtaPhiELorentzVector genZ;       bool foundZ    = false;
  math::PtEtaPhiELorentzVector genTau_et;  bool foundTau_et  = false;
  math::PtEtaPhiELorentzVector genTau_mt;  bool foundTau_mt  = false;
  math::PtEtaPhiELorentzVector genMuo_mt;  bool foundMuo_mt  = false;
  math::PtEtaPhiELorentzVector genMuo1_mm; bool foundMuo1_mm = false;
  math::PtEtaPhiELorentzVector genMuo2_mm; bool foundMuo2_mm = false;
  math::PtEtaPhiELorentzVector genMuo_em;  bool foundMuo_em = false;
  math::PtEtaPhiELorentzVector genEle_et;  bool foundEle_et  = false;
  math::PtEtaPhiELorentzVector genEle1_ee; bool foundEle1_ee = false;
  math::PtEtaPhiELorentzVector genEle2_ee; bool foundEle2_ee = false;
  math::PtEtaPhiELorentzVector genEle_em;  bool foundEle_em = false;
  for(size_t ngenPart=0; ngenPart<genParts->size(); ngenPart++){
    const reco::GenParticle & genPart = (*genParts)[ngenPart];
    if(abs(genPart.pdgId())==2212) continue;
    const reco::Candidate * mom = genPart.mother();

  
    //LOOK FOR THE Z
    if(abs(genPart.pdgId())==23 && genPart.status()!=3) {
      genZ = genPart.p4();
      foundZ = true;
    }
  
    //LOOK FOR THE HADRONIC TAU - MUOTAU
    if((ele==0 && muo==1) && abs(genPart.pdgId())==15 && genPart.status()!=3 && (abs(mom->pdgId())==25 || abs(mom->pdgId())==15) && genPart.pt()>20){
      for(unsigned int ndaugh=0; ndaugh<genPart.numberOfDaughters(); ndaugh++){
	const reco::Candidate * daughter = genPart.daughter(ndaugh);
	if(abs(daughter->pdgId())!=11 && abs(daughter->pdgId())!=12 && abs(daughter->pdgId())!=13 && abs(daughter->pdgId())!=14 
	   && abs(daughter->pdgId())!=15 && abs(daughter->pdgId())!=16){
	  genTau_mt=genPart.p4();
	  foundTau_mt = true;
	}
      }
      for(unsigned int ndaugh=0; ndaugh<genPart.numberOfDaughters(); ndaugh++){
	const reco::Candidate * daughter = genPart.daughter(ndaugh);
	if((abs(daughter->pdgId())==12 || abs(daughter->pdgId())==14 ||abs(daughter->pdgId())==16) && daughter->status()==1){
	  if(foundTau_mt) genTau_mt = genTau_mt - daughter->p4();
	}
      }
    }
  
    //LOOK FOR THE HADRONIC TAU - ELETAU
    if(((ele==1 && muo==0) || (ele==0 && muo==1)) && abs(genPart.pdgId())==15 && genPart.status()!=3 && (abs(mom->pdgId())==25 || abs(mom->pdgId())==15) && genPart.pt()>20){
      for(unsigned int ndaugh=0; ndaugh<genPart.numberOfDaughters(); ndaugh++){
	const reco::Candidate * daughter = genPart.daughter(ndaugh);
	if(abs(daughter->pdgId())!=11 && abs(daughter->pdgId())!=12 && abs(daughter->pdgId())!=13 && abs(daughter->pdgId())!=14 
	   && abs(daughter->pdgId())!=15 && abs(daughter->pdgId())!=16){
	  genTau_et=genPart.p4();
	  foundTau_et = true;
	}
      }
      for(unsigned int ndaugh=0; ndaugh<genPart.numberOfDaughters(); ndaugh++){
	const reco::Candidate * daughter = genPart.daughter(ndaugh);
	if((abs(daughter->pdgId())==12 || abs(daughter->pdgId())==14 ||abs(daughter->pdgId())==16) && daughter->status()==1){
	  if(foundTau_et) genTau_et = genTau_et - daughter->p4();
	}
      }
    }

    //LOOK FOR THE ELECTRON
    if((ele==1 && muo==0) && abs(genPart.pdgId())==11 && genPart.status()==1 && (abs(mom->pdgId())==15 || mom->status()!=3) && genPart.pt()>10){
      genEle_et = genPart.p4();
      foundEle_et = true;
    }

    //LOOK FOR THE MUON
    if((ele==0 && muo==1) && abs(genPart.pdgId())==13 && genPart.status()==1 && (abs(mom->pdgId())==15 || mom->status()!=3) && genPart.pt()>10){
      genMuo_mt = genPart.p4();
      foundMuo_mt = true;
    }
  
    //LOOK FOR THE ELECTRONS
    if((ele==2 && muo==0) && genPart.pdgId()==11  && genPart.status()==1 && (abs(mom->pdgId())==15 || mom->status()!=3) && genPart.pt()>10){
      genEle1_ee = genPart.p4();
      foundEle1_ee = true;
    }
    if((ele==2 && muo==0) && genPart.pdgId()==-11 && genPart.status()==1 && (abs(mom->pdgId())==15 || mom->status()!=3) && genPart.pt()>10){
      genEle2_ee = genPart.p4();
      foundEle2_ee = true;
    }

    //LOOK FOR THE MUONS
    if((ele==0 && muo==2) && genPart.pdgId()==13  && genPart.status()==1 && (abs(mom->pdgId())==15 || mom->status()!=3) && genPart.pt()>10){
      genMuo1_mm = genPart.p4();
      foundMuo1_mm = true;
    }
    if((ele==0 && muo==2) && genPart.pdgId()==-13 && genPart.status()==1 && (abs(mom->pdgId())==15 || mom->status()!=3) && genPart.pt()>10){
      genMuo2_mm = genPart.p4();
      foundMuo2_mm = true;
    }

    //LOOK FOR THE ELE-MUO
    if((ele==1 && muo==1) && abs(genPart.pdgId())==13  && genPart.status()==1 && (abs(mom->pdgId())==15 || mom->status()!=3) && genPart.pt()>10){
      genMuo_em = genPart.p4();
      foundMuo_em = true;
    }
    if((ele==1 && muo==1) && abs(genPart.pdgId())==11 && genPart.status()==1 && (abs(mom->pdgId())==15 || mom->status()!=3) && genPart.pt()>10){
      genEle_em = genPart.p4();
      foundEle_em = true;
    }
  }

  //MATCH THE GEN Z TO THE GEN JET
  vector<reco::GenJet>::const_iterator GenJeT; bool foundJet=false;
  for(vector<reco::GenJet>::const_iterator genjet=genjets->begin(); genjet!=genjets->end(); genjet++){
    if(foundZ) {if(ROOT::Math::VectorUtil::DeltaR(genZ,genjet->p4())<0.3) {GenJeT=genjet; foundJet=true;}}
  }

  
  //JET SELECTION - SR
  pat::JetCollection::const_iterator SelectedJet; bool matchedJet = false;
  float massZ=-9999; float ptZ=-999; bool foundSelectedJet=false; float tau21Z=-9999;
  SelectJet(CA8JetswithQjets, CA8JetsPruned, foundSelectedJet, SelectedJet, massZ, tau21Z, ptZ, 70, 110, foundJet, GenJeT, matchedJet);
  
  //TAU SELECTION - MT - SR
  float ptTau_mt=-99; bool foundSelectedTau_mt=false;
  pat::TauCollection::const_iterator SelectedTau_mt; bool matchedTau_mt = false;
  SelectTau(tauMuTauHandle, SelectedJet, foundSelectedTau_mt, SelectedTau_mt, ptTau_mt, foundSelectedJet, foundTau_mt, genTau_mt, matchedTau_mt);
  
  //TAU SELECTION - ET - SR
  float ptTau_et=-99; bool foundSelectedTau_et=false;
  pat::TauCollection::const_iterator SelectedTau_et; bool matchedTau_et = false;
  SelectTau(tauElTauHandle, SelectedJet, foundSelectedTau_et, SelectedTau_et, ptTau_et, foundSelectedJet, foundTau_et, genTau_et, matchedTau_et); 
  
  //ELECTRON SELECTION - ET - SR
  float ptElectron_et=-99; bool foundSelectedEle_et=false;
  pat::ElectronCollection::const_iterator SelectedElectron_et; bool matchedEle_et = false;
  SelectElectron(eleH, SelectedJet, foundSelectedEle_et, SelectedElectron_et, ptElectron_et, foundSelectedJet, primaryVertex, foundEle_et, genEle_et, matchedEle_et);
  
  //ELECTRON SELECTION - EE1 - SR
  float ptElectron1_ee=-99; bool foundSelectedEle1_ee=false;
  pat::ElectronCollection::const_iterator SelectedElectron1_ee; bool matchedEle1_ee = false;
  SelectElectron(eleH, SelectedJet, foundSelectedEle1_ee, SelectedElectron1_ee, ptElectron1_ee, foundSelectedJet, primaryVertex, foundEle1_ee, genEle1_ee, matchedEle1_ee);
  
  //ELECTRON SELECTION - EE2 - SR
  float ptElectron2_ee=-99; bool foundSelectedEle2_ee=false;
  pat::ElectronCollection::const_iterator SelectedElectron2_ee; bool matchedEle2_ee = false;
  SelectElectron(eleH, SelectedJet, foundSelectedEle2_ee, SelectedElectron2_ee, ptElectron2_ee, foundSelectedJet, primaryVertex, foundEle2_ee, genEle2_ee, matchedEle2_ee);
  
  //ELECTRON SELECTION - EM - SR
  float ptElectron_em=-99; bool foundSelectedEle_em=false;
  pat::ElectronCollection::const_iterator SelectedElectron_em; bool matchedEle_em = false;
  SelectElectron(eleH, SelectedJet, foundSelectedEle_em, SelectedElectron_em, ptElectron_em, foundSelectedJet, primaryVertex, foundEle_em, genEle_em, matchedEle_em);

  //MUON SELECTION - MT - SR
  float ptHighptMuon_mt=-99; bool foundSelectedHighptMuo_mt=false;
  pat::MuonCollection::const_iterator SelectedHighptMuon_mt; bool matchedHighptMuo_mt = false;
  SelectHighptMuon(muoH,SelectedJet,foundSelectedHighptMuo_mt,SelectedHighptMuon_mt,ptHighptMuon_mt,foundSelectedJet,primaryVertex,foundMuo_mt,genMuo_mt,matchedHighptMuo_mt);
  float ptTightMuon_mt=-99; bool foundSelectedTightMuo_mt=false;
  pat::MuonCollection::const_iterator SelectedTightMuon_mt; bool matchedTightMuo_mt = false;
  SelectTightMuon(muoH,SelectedJet,foundSelectedTightMuo_mt,SelectedTightMuon_mt,ptTightMuon_mt,foundSelectedJet,primaryVertex,foundMuo_mt,genMuo_mt,matchedTightMuo_mt);
  
  //MUON SELECTION - MM1 - SR
  float ptTrackerMuon1_mm=-99; bool foundSelectedTrackerMuo1_mm=false;
  pat::MuonCollection::const_iterator SelectedTrackerMuon1_mm; bool matchedTrackerMuo1_mm = false;
  SelectTrackerMuon(muoH,foundSelectedTrackerMuo1_mm,SelectedTrackerMuon1_mm,ptTrackerMuon1_mm,primaryVertex,SelectedJet,foundSelectedJet,foundMuo1_mm,genMuo1_mm,matchedTrackerMuo1_mm);
  float ptHighptMuon1_mm=-99; bool foundSelectedHighptMuo1_mm=false;
  pat::MuonCollection::const_iterator SelectedHighptMuon1_mm; bool matchedHighptMuo1_mm = false;
  SelectHighptMuon(muoH,SelectedJet,foundSelectedHighptMuo1_mm,SelectedHighptMuon1_mm,ptHighptMuon1_mm,foundSelectedJet,primaryVertex,foundMuo1_mm,genMuo1_mm,matchedHighptMuo1_mm);
  float ptTightMuon1_mm=-99; bool foundSelectedTightMuo1_mm=false;
  pat::MuonCollection::const_iterator SelectedTightMuon1_mm; bool matchedTightMuo1_mm = false;
  SelectTightMuon(muoH,SelectedJet,foundSelectedTightMuo1_mm,SelectedTightMuon1_mm,ptTightMuon1_mm,foundSelectedJet,primaryVertex,foundMuo1_mm,genMuo1_mm,matchedTightMuo1_mm);

  //MUON SELECTION - MM2 - SR
  float ptTrackerMuon2_mm=-99; bool foundSelectedTrackerMuo2_mm=false;
  pat::MuonCollection::const_iterator SelectedTrackerMuon2_mm; bool matchedTrackerMuo2_mm = false;
  SelectTrackerMuon(muoH,foundSelectedTrackerMuo2_mm,SelectedTrackerMuon2_mm,ptTrackerMuon2_mm,primaryVertex,SelectedJet,foundSelectedJet,foundMuo2_mm,genMuo2_mm,matchedTrackerMuo2_mm);
  float ptHighptMuon2_mm=-99; bool foundSelectedHighptMuo2_mm=false;
  pat::MuonCollection::const_iterator SelectedHighptMuon2_mm; bool matchedHighptMuo2_mm = false;
  SelectHighptMuon(muoH,SelectedJet,foundSelectedHighptMuo2_mm,SelectedHighptMuon2_mm,ptHighptMuon2_mm,foundSelectedJet,primaryVertex,foundMuo2_mm,genMuo2_mm,matchedHighptMuo2_mm);
  float ptTightMuon2_mm=-99; bool foundSelectedTightMuo2_mm=false;
  pat::MuonCollection::const_iterator SelectedTightMuon2_mm; bool matchedTightMuo2_mm = false;
  SelectTightMuon(muoH,SelectedJet,foundSelectedTightMuo2_mm,SelectedTightMuon2_mm,ptTightMuon2_mm,foundSelectedJet,primaryVertex,foundMuo2_mm,genMuo2_mm,matchedTightMuo2_mm);

  //MUON-MUON SELECTION
  bool hasAtLeastOneHighPtMuo=false;
  SelectTrackerGlobalID(SelectedTrackerMuon1_mm, SelectedTrackerMuon2_mm, foundSelectedTrackerMuo1_mm, foundSelectedTrackerMuo2_mm, primaryVertex, hasAtLeastOneHighPtMuo);

  //MUON SELECTION - EM - SR
  float ptHighptMuon_em=-99; bool foundSelectedHighptMuo_em=false;
  pat::MuonCollection::const_iterator SelectedHighptMuon_em; bool matchedHighptMuo_em = false;
  SelectHighptMuon(muoH,SelectedJet,foundSelectedHighptMuo_em,SelectedHighptMuon_em,ptHighptMuon_em,foundSelectedJet,primaryVertex,foundMuo_em,genMuo_em,matchedHighptMuo_em);
  float ptTightMuon_em=-99; bool foundSelectedTightMuo_em=false;
  pat::MuonCollection::const_iterator SelectedTightMuon_em; bool matchedTightMuo_em = false;
  SelectTightMuon(muoH,SelectedJet,foundSelectedTightMuo_em,SelectedTightMuon_em,ptTightMuon_em,foundSelectedJet,primaryVertex,foundMuo_em,genMuo_em,matchedTightMuo_em);

  float genJet_Pt   = -99.;
  float genJet_Eta  = -99.;
  float genTau_mt_Pt   = -99.;
  float genTau_mt_eta  = -99.;
  float genTau_et_Pt   = -99.;
  float genTau_et_Eta  = -99.;
  float genEle_et_Pt   = -99.;
  float genEle_et_Eta  = -99.;
  float genEle1_ee_Pt  = -99.;
  float genEle1_ee_Eta = -99.;
  float genEle2_ee_Pt  = -99.;
  float genEle2_ee_Eta = -99.;
  float genEle_em_Pt  = -99.;
  float genEle_em_Eta = -99.;
  float genMuo_mt_Pt   = -99.;
  float genMuo_mt_Eta  = -99.;
  float genMuo1_mm_Pt  = -99.;
  float genMuo1_mm_Eta = -99.;
  float genMuo2_mm_Pt  = -99.;
  float genMuo2_mm_Eta = -99.;
  float genMuo_em_Pt  = -99.;
  float genMuo_em_Eta = -99.;

  float recoJet_Pt   = -99.;
  float recoJet_Eta  = -99.;
  float recoTau_et_Pt   = -99.;
  float recoTau_et_Eta  = -99.;
  float recoTau_mt_Pt   = -99.;
  float recoTau_mt_eta  = -99.;
  float recoEle_et_Pt   = -99.;
  float recoEle_et_Eta  = -99.;
  float recoEle1_ee_Pt  = -99.;
  float recoEle1_ee_Eta = -99.;
  float recoEle2_ee_Pt  = -99.;
  float recoEle2_ee_Eta = -99.;
  float recoEle_em_Pt   = -99.;
  float recoEle_em_Eta  = -99.;
  float recoTightMuo_mt_Pt   = -99.;
  float recoTightMuo_mt_Eta  = -99.;
  float recoTightMuo1_mm_Pt  = -99.;
  float recoTightMuo1_mm_Eta = -99.;
  float recoTightMuo2_mm_Pt  = -99.;
  float recoTightMuo2_mm_Eta = -99.;
  float recoTightMuo_em_Pt   = -99.;
  float recoTightMuo_em_Eta  = -99.;
  float recoHighptMuo_mt_Pt   = -99.;
  float recoHighptMuo_mt_Eta  = -99.;
  float recoHighptMuo1_mm_Pt  = -99.;
  float recoHighptMuo1_mm_Eta = -99.;
  float recoHighptMuo2_mm_Pt  = -99.;
  float recoHighptMuo2_mm_Eta = -99.;
  float recoHighptMuo_em_Pt   = -99.;
  float recoHighptMuo_em_Eta  = -99.;
  float recoTrackerMuo1_mm_Pt  = -99.;
  float recoTrackerMuo1_mm_Eta = -99.;
  float recoTrackerMuo2_mm_Pt  = -99.;
  float recoTrackerMuo2_mm_Eta = -99.;
  float recoEle_et_PFIso      = -99.;
  float recoEle_et_CorrPFIso  = -99.;
  float recoEle1_ee_PFIso     = -99.;
  float recoEle1_ee_CorrPFIso = -99.;
  float recoEle2_ee_PFIso     = -99.;
  float recoEle2_ee_CorrPFIso = -99.;
  float recoEle_em_PFIso      = -99.;
  float recoEle_em_CorrPFIso  = -99.;
  float recoTightMuo_mt_PFIso      = -99.;
  float recoTightMuo_mt_CorrPFIso  = -99.;
  float recoTightMuo1_mm_PFIso     = -99.;
  float recoTightMuo1_mm_CorrPFIso = -99.;
  float recoTightMuo1_mm_DetIso    = -99.;
  float recoTightMuo2_mm_PFIso     = -99.;
  float recoTightMuo2_mm_CorrPFIso = -99.;
  float recoTightMuo2_mm_DetIso    = -99.;
  float recoTightMuo_em_PFIso      = -99.;
  float recoTightMuo_em_CorrPFIso  = -99.;
  float recoHighptMuo_mt_PFIso      = -99.;
  float recoHighptMuo_mt_CorrPFIso  = -99.;
  float recoHighptMuo1_mm_PFIso     = -99.;
  float recoHighptMuo1_mm_CorrPFIso = -99.;
  float recoHighptMuo1_mm_DetIso    = -99.;
  float recoHighptMuo2_mm_PFIso     = -99.;
  float recoHighptMuo2_mm_CorrPFIso = -99.;
  float recoHighptMuo2_mm_DetIso    = -99.;
  float recoHighptMuo_em_PFIso      = -99.;
  float recoHighptMuo_em_CorrPFIso  = -99.;
  float recoTrackerMuo1_mm_PFIso     = -99.;
  float recoTrackerMuo1_mm_CorrPFIso = -99.;
  float recoTrackerMuo1_mm_DetIso    = -99.;
  float recoTrackerMuo2_mm_PFIso     = -99.;
  float recoTrackerMuo2_mm_CorrPFIso = -99.;
  float recoTrackerMuo2_mm_DetIso    = -99.;
  float reco_em_deltaR = -99;
  float reco_et_deltaR = -99;
  float reco_ee_deltaR = -99;
  float reco_mt_deltaR = -99;
  float reco_mm_deltaR = -99;

  float dRJet   = 99.;
  float dRTau_et   = 99.;
  float dRTau_mt   = 99.;
  float dREle_et   = 99.;
  float dREle1_ee  = 99.;
  float dREle2_ee  = 99.;
  float dREle_em  = 99.;
  float dRMuo_mt   = 99.;
  float dRMuo1_mm  = 99.;
  float dRMuo2_mm  = 99.;
  float dRMuo_em  = 99.;

  if(foundJet){
    genJet_Pt = GenJeT->pt();
    genJet_Eta = GenJeT->eta();
  }
  if(foundTau_et){
    genTau_et_Pt = genTau_et.pt();
    genTau_et_Eta = genTau_et.eta();
  }
  if(foundTau_mt){
    genTau_mt_Pt = genTau_mt.pt();
    genTau_mt_eta = genTau_mt.eta();
  }
  if(foundEle_et){
    genEle_et_Pt = genEle_et.pt();
    genEle_et_Eta = genEle_et.eta();
  }
  if(foundEle1_ee){
    genEle1_ee_Pt = genEle1_ee.pt();
    genEle1_ee_Eta = genEle1_ee.eta();
  }
  if(foundEle2_ee){
    genEle2_ee_Pt = genEle2_ee.pt();
    genEle2_ee_Eta = genEle2_ee.eta();
  }
  if(foundEle_em){
    genEle_em_Pt = genEle_em.pt();
    genEle_em_Eta = genEle_em.eta();
  }
  if(foundMuo_mt){
    genMuo_mt_Pt = genMuo_mt.pt();
    genMuo_mt_Eta = genMuo_mt.eta();
  }
  if(foundMuo1_mm){
    genMuo1_mm_Pt = genMuo1_mm.pt();
    genMuo1_mm_Eta = genMuo1_mm.eta();
  }
  if(foundMuo2_mm){
    genMuo2_mm_Pt = genMuo2_mm.pt();
    genMuo2_mm_Eta = genMuo2_mm.eta();
  }
  if(foundMuo_em){
    genMuo_em_Pt = genMuo_em.pt();
    genMuo_em_Eta = genMuo_em.eta();
  }

  if(foundSelectedJet){
    recoJet_Pt = SelectedJet->pt();
    recoJet_Eta = SelectedJet->eta();
  }
  if(foundSelectedTau_mt){
    recoTau_mt_Pt = SelectedTau_mt->pt();
    recoTau_mt_eta = SelectedTau_mt->eta();
  }
  if(foundSelectedTau_et){
    recoTau_et_Pt = SelectedTau_et->pt();
    recoTau_et_Eta = SelectedTau_et->eta();
  }
  if(foundSelectedEle_et){
    recoEle_et_Pt = SelectedElectron_et->pt();
    recoEle_et_Eta = SelectedElectron_et->eta();
    recoEle_et_PFIso=ElectronPFIso(SelectedElectron_et, rho);
    recoEle_et_CorrPFIso=ElectronCorrPFIso(SelectedElectron_et, rho);
  }
  if(foundSelectedEle1_ee){
    recoEle1_ee_Pt = SelectedElectron1_ee->pt();
    recoEle1_ee_Eta = SelectedElectron1_ee->eta();
    recoEle1_ee_PFIso=ElectronPFIso(SelectedElectron1_ee, rho);
    recoEle1_ee_CorrPFIso=ElectronCorrPFIso(SelectedElectron1_ee, rho);
  }
  if(foundSelectedEle2_ee){
    recoEle2_ee_Pt = SelectedElectron2_ee->pt();
    recoEle2_ee_Eta = SelectedElectron2_ee->eta();
    recoEle2_ee_PFIso=ElectronPFIso(SelectedElectron2_ee, rho);
    recoEle2_ee_CorrPFIso=ElectronCorrPFIso(SelectedElectron2_ee, rho);
  }
  if(foundSelectedEle_em){
    recoEle_em_Pt = SelectedElectron_em->pt();
    recoEle_em_Eta = SelectedElectron_em->eta();
    recoEle_em_PFIso=ElectronPFIso(SelectedElectron_em, rho);
    recoEle_em_CorrPFIso=ElectronCorrPFIso(SelectedElectron_em, rho);
  }

  if(foundSelectedTightMuo_mt){
    recoTightMuo_mt_Pt = SelectedTightMuon_mt->pt();
    recoTightMuo_mt_Eta =SelectedTightMuon_mt->eta();
    recoTightMuo_mt_PFIso=MuonPFIso(SelectedTightMuon_mt,true);
    recoTightMuo_mt_CorrPFIso=MuonCorrPFIso(SelectedTightMuon_mt,true);
  }
  if(foundSelectedTightMuo1_mm){
    recoTightMuo1_mm_Pt = SelectedTightMuon1_mm->pt();
    recoTightMuo1_mm_Eta =SelectedTightMuon1_mm->eta();
    recoTightMuo1_mm_PFIso=MuonPFIso(SelectedTightMuon1_mm,true);
    recoTightMuo1_mm_CorrPFIso=MuonCorrPFIso(SelectedTightMuon1_mm,true);
    recoTightMuo1_mm_DetIso=MuonDETIso(SelectedTightMuon1_mm,SelectedTightMuon2_mm, foundSelectedTightMuo2_mm);
  }
  if(foundSelectedTightMuo2_mm){
    recoTightMuo2_mm_Pt = SelectedTightMuon2_mm->pt();
    recoTightMuo2_mm_Eta =SelectedTightMuon2_mm->eta();
    recoTightMuo2_mm_PFIso=MuonPFIso(SelectedTightMuon2_mm,true);
    recoTightMuo2_mm_CorrPFIso=MuonCorrPFIso(SelectedTightMuon2_mm,true);
    recoTightMuo2_mm_DetIso=MuonDETIso(SelectedTightMuon2_mm,SelectedTightMuon1_mm, foundSelectedTightMuo1_mm);
  }
  if(foundSelectedTightMuo_em){
    recoTightMuo_em_Pt = SelectedTightMuon_em->pt();
    recoTightMuo_em_Eta = SelectedTightMuon_em->eta();
    recoTightMuo_em_PFIso=MuonPFIso(SelectedTightMuon_em,true);
    recoTightMuo_em_CorrPFIso=MuonCorrPFIso(SelectedTightMuon_em,true);
  }

  if(foundSelectedHighptMuo_mt){
    recoHighptMuo_mt_Pt = SelectedHighptMuon_mt->pt();
    recoHighptMuo_mt_Eta =SelectedHighptMuon_mt->eta();
    recoHighptMuo_mt_PFIso=MuonPFIso(SelectedHighptMuon_mt,true);
    recoHighptMuo_mt_CorrPFIso=MuonCorrPFIso(SelectedHighptMuon_mt,true);
  }
  if(foundSelectedHighptMuo1_mm){
    recoHighptMuo1_mm_Pt = SelectedHighptMuon1_mm->pt();
    recoHighptMuo1_mm_Eta =SelectedHighptMuon1_mm->eta();
    recoHighptMuo1_mm_PFIso=MuonPFIso(SelectedHighptMuon1_mm,true);
    recoHighptMuo1_mm_CorrPFIso=MuonCorrPFIso(SelectedHighptMuon1_mm,true);
    recoHighptMuo1_mm_DetIso=MuonDETIso(SelectedHighptMuon1_mm,SelectedHighptMuon2_mm, foundSelectedHighptMuo2_mm);
  }
  if(foundSelectedHighptMuo2_mm){
    recoHighptMuo2_mm_Pt = SelectedHighptMuon2_mm->pt();
    recoHighptMuo2_mm_Eta =SelectedHighptMuon2_mm->eta();
    recoHighptMuo2_mm_PFIso=MuonPFIso(SelectedHighptMuon2_mm,true);
    recoHighptMuo2_mm_CorrPFIso=MuonCorrPFIso(SelectedHighptMuon2_mm,true);
    recoHighptMuo2_mm_DetIso=MuonDETIso(SelectedHighptMuon2_mm,SelectedHighptMuon1_mm, foundSelectedHighptMuo1_mm);
  }
  if(foundSelectedHighptMuo_em){
    recoHighptMuo_em_Pt = SelectedHighptMuon_em->pt();
    recoHighptMuo_em_Eta = SelectedHighptMuon_em->eta();
    recoHighptMuo_em_PFIso=MuonPFIso(SelectedHighptMuon_em,true);
    recoHighptMuo_em_CorrPFIso=MuonCorrPFIso(SelectedHighptMuon_em,true);
  }

  if(foundSelectedTrackerMuo1_mm){
    recoTrackerMuo1_mm_Pt = SelectedTrackerMuon1_mm->pt();
    recoTrackerMuo1_mm_Eta =SelectedTrackerMuon1_mm->eta();
    recoTrackerMuo1_mm_PFIso=MuonPFIso(SelectedTrackerMuon1_mm,true);
    recoTrackerMuo1_mm_CorrPFIso=MuonCorrPFIso(SelectedTrackerMuon1_mm,true);
    recoTrackerMuo1_mm_DetIso=MuonDETIso(SelectedTrackerMuon1_mm,SelectedTrackerMuon2_mm, foundSelectedTrackerMuo2_mm);
  }
  if(foundSelectedTrackerMuo2_mm){
    recoTrackerMuo2_mm_Pt = SelectedTrackerMuon2_mm->pt();
    recoTrackerMuo2_mm_Eta =SelectedTrackerMuon2_mm->eta();
    recoTrackerMuo2_mm_PFIso=MuonPFIso(SelectedTrackerMuon2_mm,true);
    recoTrackerMuo2_mm_CorrPFIso=MuonCorrPFIso(SelectedTrackerMuon2_mm,true);
    recoTrackerMuo2_mm_DetIso=MuonDETIso(SelectedTrackerMuon2_mm,SelectedTrackerMuon1_mm, foundSelectedTrackerMuo1_mm);
  }

  if(foundSelectedEle_em && foundSelectedHighptMuo_em){
    reco_em_deltaR=ROOT::Math::VectorUtil::DeltaR(SelectedElectron_em->p4(),SelectedHighptMuon_em->p4());
  }
  if(foundSelectedEle_et && foundSelectedTau_et){
    reco_et_deltaR=ROOT::Math::VectorUtil::DeltaR(SelectedElectron_et->p4(),SelectedTau_et->p4());
  }
  if(foundSelectedHighptMuo_mt && foundSelectedTau_mt){
    reco_mt_deltaR=ROOT::Math::VectorUtil::DeltaR(SelectedHighptMuon_mt->p4(),SelectedTau_mt->p4());
  }
  if(foundSelectedTrackerMuo1_mm && foundSelectedTrackerMuo2_mm){
    reco_mm_deltaR=ROOT::Math::VectorUtil::DeltaR(SelectedTrackerMuon1_mm->p4(),SelectedTrackerMuon2_mm->p4());
  }
  if(foundSelectedEle1_ee && foundSelectedEle2_ee){
    reco_ee_deltaR=ROOT::Math::VectorUtil::DeltaR(SelectedElectron1_ee->p4(),SelectedElectron2_ee->p4());
  }

  if(foundSelectedJet     && foundJet)     dRJet     = ROOT::Math::VectorUtil::DeltaR(GenJeT->p4(),SelectedJet->p4());
  if(foundSelectedTau_mt  && foundTau_mt)  dRTau_mt  = ROOT::Math::VectorUtil::DeltaR(genTau_mt,   SelectedTau_mt->p4());
  if(foundSelectedTau_et  && foundTau_et)  dRTau_et  = ROOT::Math::VectorUtil::DeltaR(genTau_et,   SelectedTau_et->p4());
  if(foundSelectedEle_et  && foundEle_et)  dREle_et  = ROOT::Math::VectorUtil::DeltaR(genEle_et,   SelectedElectron_et->p4());
  if(foundSelectedEle1_ee && foundEle1_ee) dREle1_ee = ROOT::Math::VectorUtil::DeltaR(genEle1_ee,  SelectedElectron1_ee->p4());
  if(foundSelectedEle2_ee && foundEle2_ee) dREle2_ee = ROOT::Math::VectorUtil::DeltaR(genEle2_ee,  SelectedElectron2_ee->p4());
  if(foundSelectedEle_em  && foundEle_em)  dREle_em  = ROOT::Math::VectorUtil::DeltaR(genEle_em,   SelectedElectron_em->p4());
  if(foundSelectedHighptMuo_mt  && foundMuo_mt)   dRMuo_mt  = ROOT::Math::VectorUtil::DeltaR(genMuo_mt,   SelectedHighptMuon_mt->p4());
  if(foundSelectedTrackerMuo1_mm && foundMuo1_mm) dRMuo1_mm = ROOT::Math::VectorUtil::DeltaR(genMuo1_mm,  SelectedTrackerMuon1_mm->p4());
  if(foundSelectedTrackerMuo2_mm && foundMuo2_mm) dRMuo2_mm = ROOT::Math::VectorUtil::DeltaR(genMuo2_mm,  SelectedTrackerMuon2_mm->p4());
  if(foundSelectedHighptMuo_em  && foundMuo_em)   dRMuo_em  = ROOT::Math::VectorUtil::DeltaR(genMuo_em,   SelectedHighptMuon_em->p4());
 
  m_Nele = (int)ele;
  m_Nmuo = (int)muo;
  m_genJet   = (int)foundJet;
  m_genTau_et   = (int)foundTau_et;
  m_genTau_mt   = (int)foundTau_mt;
  m_genEle_et   = (int)foundEle_et;
  m_genEle1_ee  = (int)foundEle1_ee;
  m_genEle2_ee  = (int)foundEle2_ee;
  m_genEle_em  = (int)foundEle_em;
  m_genMuo_mt   = (int)foundMuo_mt;
  m_genMuo1_mm  = (int)foundMuo1_mm;
  m_genMuo2_mm  = (int)foundMuo2_mm;
  m_genMuo_em  = (int)foundMuo_em;
  m_matchedJet   = (int)matchedJet;
  m_matchedTau_et   = (int)matchedTau_et;
  m_matchedTau_mt   = (int)matchedTau_mt;
  m_matchedEle_et   = (int)matchedEle_et;
  m_matchedEle1_ee  = (int)matchedEle1_ee;
  m_matchedEle2_ee  = (int)matchedEle2_ee;
  m_matchedEle_em  = (int)matchedEle_em;
  m_matchedTightMuo_mt    = (int)matchedTightMuo_mt;
  m_matchedTightMuo1_mm   = (int)matchedTightMuo1_mm;
  m_matchedTightMuo2_mm   = (int)matchedTightMuo2_mm;
  m_matchedTightMuo_em    = (int)matchedTightMuo_em;
  m_matchedHighptMuo_mt   = (int)matchedHighptMuo_mt;
  m_matchedHighptMuo1_mm  = (int)matchedHighptMuo1_mm;
  m_matchedHighptMuo2_mm  = (int)matchedHighptMuo2_mm;
  m_matchedHighptMuo_em   = (int)matchedHighptMuo_em;
  m_matchedTrackerMuo1_mm = (int)matchedTrackerMuo1_mm;
  m_matchedTrackerMuo2_mm = (int)matchedTrackerMuo2_mm;
  m_recoJet  = (int)foundSelectedJet;
  m_recoTau_mt  = (int)foundSelectedTau_mt;
  m_recoTau_et  = (int)foundSelectedTau_et;
  m_recoEle_et  = (int)foundSelectedEle_et;
  m_recoEle1_ee = (int)foundSelectedEle1_ee;
  m_recoEle2_ee = (int)foundSelectedEle2_ee;
  m_recoEle_em  = (int)foundSelectedEle_em;
  m_recoTightMuo_mt    = (int)foundSelectedTightMuo_mt;
  m_recoTightMuo1_mm   = (int)foundSelectedTightMuo1_mm;
  m_recoTightMuo2_mm   = (int)foundSelectedTightMuo2_mm;
  m_recoTightMuo_em    = (int)foundSelectedTightMuo_em;
  m_recoHighptMuo_mt   = (int)foundSelectedHighptMuo_mt;
  m_recoHighptMuo1_mm  = (int)foundSelectedHighptMuo1_mm;
  m_recoHighptMuo2_mm  = (int)foundSelectedHighptMuo2_mm;
  m_recoHighptMuo_em   = (int)foundSelectedHighptMuo_em;
  m_recoTrackerMuo1_mm = (int)foundSelectedTrackerMuo1_mm;
  m_recoTrackerMuo2_mm = (int)foundSelectedTrackerMuo2_mm;
  m_hasAtLeastOneHighPtMuo = (int)hasAtLeastOneHighPtMuo;
  m_genJet_Pt   = genJet_Pt;
  m_genJet_Eta  = genJet_Eta;
  m_genTau_et_Pt   = genTau_et_Pt;
  m_genTau_et_Eta  = genTau_et_Eta;
  m_genTau_mt_Pt   = genTau_mt_Pt;
  m_genTau_mt_eta  = genTau_mt_eta;
  m_genEle_et_Pt   = genEle_et_Pt;
  m_genEle_et_Eta  = genEle_et_Eta;
  m_genEle1_ee_Pt  = genEle1_ee_Pt;
  m_genEle1_ee_Eta = genEle1_ee_Eta;
  m_genEle2_ee_Pt  = genEle2_ee_Pt;
  m_genEle2_ee_Eta = genEle2_ee_Eta;
  m_genEle_em_Pt  = genEle_em_Pt;
  m_genEle_em_Eta = genEle_em_Eta;
  m_genMuo_mt_Pt   = genMuo_mt_Pt;
  m_genMuo_mt_Eta  = genMuo_mt_Eta;
  m_genMuo1_mm_Pt  = genMuo1_mm_Pt;
  m_genMuo1_mm_Eta = genMuo1_mm_Eta;
  m_genMuo2_mm_Pt  = genMuo2_mm_Pt;
  m_genMuo2_mm_Eta = genMuo2_mm_Eta;
  m_genMuo_em_Pt  = genMuo_em_Pt;
  m_genMuo_em_Eta = genMuo_em_Eta;
  m_recoJet_Pt   = recoJet_Pt;
  m_recoJet_Eta  = recoJet_Eta;
  m_recoTau_et_Pt   = recoTau_et_Pt;
  m_recoTau_et_Eta  = recoTau_et_Eta;
  m_recoTau_mt_Pt   = recoTau_mt_Pt;
  m_recoTau_mt_eta  = recoTau_mt_eta;
  m_recoEle_et_Pt   = recoEle_et_Pt;
  m_recoEle_et_Eta  = recoEle_et_Eta;
  m_recoEle1_ee_Pt  = recoEle1_ee_Pt;
  m_recoEle1_ee_Eta = recoEle1_ee_Eta;
  m_recoEle2_ee_Pt  = recoEle2_ee_Pt;
  m_recoEle2_ee_Eta = recoEle2_ee_Eta;
  m_recoEle_em_Pt  = recoEle_em_Pt;
  m_recoEle_em_Eta = recoEle_em_Eta;
  m_recoTightMuo_mt_Pt   = recoTightMuo_mt_Pt;
  m_recoTightMuo_mt_Eta  = recoTightMuo_mt_Eta;
  m_recoTightMuo1_mm_Pt  = recoTightMuo1_mm_Pt;
  m_recoTightMuo1_mm_Eta = recoTightMuo1_mm_Eta;
  m_recoTightMuo2_mm_Pt  = recoTightMuo2_mm_Pt;
  m_recoTightMuo2_mm_Eta = recoTightMuo2_mm_Eta;
  m_recoTightMuo_em_Pt   = recoTightMuo_em_Pt;
  m_recoTightMuo_em_Eta  = recoTightMuo_em_Eta;
  m_recoHighptMuo_mt_Pt   = recoHighptMuo_mt_Pt;
  m_recoHighptMuo_mt_Eta  = recoHighptMuo_mt_Eta;
  m_recoHighptMuo1_mm_Pt  = recoHighptMuo1_mm_Pt;
  m_recoHighptMuo1_mm_Eta = recoHighptMuo1_mm_Eta;
  m_recoHighptMuo2_mm_Pt  = recoHighptMuo2_mm_Pt;
  m_recoHighptMuo2_mm_Eta = recoHighptMuo2_mm_Eta;
  m_recoHighptMuo_em_Pt   = recoHighptMuo_em_Pt;
  m_recoHighptMuo_em_Eta  = recoHighptMuo_em_Eta;
  m_recoTrackerMuo1_mm_Pt  = recoTrackerMuo1_mm_Pt;
  m_recoTrackerMuo1_mm_Eta = recoTrackerMuo1_mm_Eta;
  m_recoTrackerMuo2_mm_Pt  = recoTrackerMuo2_mm_Pt;
  m_recoTrackerMuo2_mm_Eta = recoTrackerMuo2_mm_Eta;
  m_recoEle_et_PFIso   = recoEle_et_PFIso;
  m_recoEle_et_CorrPFIso  = recoEle_et_CorrPFIso;
  m_recoEle1_ee_PFIso  = recoEle1_ee_PFIso;
  m_recoEle1_ee_CorrPFIso = recoEle1_ee_CorrPFIso;
  m_recoEle2_ee_PFIso  = recoEle2_ee_PFIso;
  m_recoEle2_ee_CorrPFIso = recoEle2_ee_CorrPFIso;
  m_recoEle_em_PFIso  = recoEle_em_PFIso;
  m_recoEle_em_CorrPFIso = recoEle_em_CorrPFIso;
  m_recoTightMuo_mt_PFIso   = recoTightMuo_mt_PFIso;
  m_recoTightMuo_mt_CorrPFIso  = recoTightMuo_mt_CorrPFIso;
  m_recoTightMuo1_mm_PFIso  = recoTightMuo1_mm_PFIso;
  m_recoTightMuo1_mm_CorrPFIso = recoTightMuo1_mm_CorrPFIso;
  m_recoTightMuo2_mm_PFIso  = recoTightMuo2_mm_PFIso;
  m_recoTightMuo2_mm_CorrPFIso = recoTightMuo2_mm_CorrPFIso;
  m_recoTightMuo_em_PFIso  = recoTightMuo_em_PFIso;
  m_recoTightMuo_em_CorrPFIso = recoTightMuo_em_CorrPFIso;
  m_recoHighptMuo_mt_PFIso   = recoHighptMuo_mt_PFIso;
  m_recoHighptMuo_mt_CorrPFIso  = recoHighptMuo_mt_CorrPFIso;
  m_recoHighptMuo1_mm_PFIso  = recoHighptMuo1_mm_PFIso;
  m_recoHighptMuo1_mm_CorrPFIso = recoHighptMuo1_mm_CorrPFIso;
  m_recoHighptMuo2_mm_PFIso  = recoHighptMuo2_mm_PFIso;
  m_recoHighptMuo2_mm_CorrPFIso = recoHighptMuo2_mm_CorrPFIso;
  m_recoHighptMuo_em_PFIso  = recoHighptMuo_em_PFIso;
  m_recoHighptMuo_em_CorrPFIso = recoHighptMuo_em_CorrPFIso;
  m_recoTrackerMuo1_mm_PFIso  = recoTrackerMuo1_mm_PFIso;
  m_recoTrackerMuo1_mm_CorrPFIso = recoTrackerMuo1_mm_CorrPFIso;
  m_recoTrackerMuo2_mm_PFIso  = recoTrackerMuo2_mm_PFIso;
  m_recoTrackerMuo2_mm_CorrPFIso = recoTrackerMuo2_mm_CorrPFIso;
  m_recoTightMuo1_mm_DetIso   = recoTightMuo1_mm_DetIso;
  m_recoTightMuo2_mm_DetIso   = recoTightMuo2_mm_DetIso;
  m_recoHighptMuo1_mm_DetIso  = recoHighptMuo1_mm_DetIso;
  m_recoHighptMuo2_mm_DetIso  = recoHighptMuo2_mm_DetIso;
  m_recoTrackerMuo1_mm_DetIso = recoTrackerMuo1_mm_DetIso;
  m_recoTrackerMuo2_mm_DetIso = recoTrackerMuo2_mm_DetIso;
  m_reco_em_deltaR = reco_em_deltaR; 
  m_reco_et_deltaR = reco_et_deltaR;
  m_reco_ee_deltaR = reco_ee_deltaR;
  m_reco_mt_deltaR = reco_mt_deltaR;
  m_reco_mm_deltaR = reco_mm_deltaR;
  m_dRJet       = dRJet;
  m_dRTau_et  = dRTau_et;
  m_dRTau_mt  = dRTau_mt;
  m_dREle_et   = dREle_et;
  m_dREle1_ee  = dREle1_ee;
  m_dREle2_ee  = dREle2_ee;
  m_dREle_em  = dREle_em;
  m_dRMuo_mt   = dRMuo_mt;
  m_dRMuo1_mm  = dRMuo1_mm;
  m_dRMuo2_mm  = dRMuo2_mm;
  m_dRMuo_em  = dRMuo_em;
  Tree->Fill();

  //MUON ISOLATION FOR QCD JETS
  for(pat::MuonCollection::const_iterator muon = muoH->begin(); muon != muoH->end(); ++muon) {
    m_JetMuon_Pt        = muon->pt();
    m_JetMuon_Eta       = muon->eta();
    m_JetMuon_PFIso     = MuonPFIso(muon,true);
    m_JetMuon_CorrPFIso = MuonCorrPFIso(muon,true);
    TreeQCDMuon->Fill();
  }

  //ELECTRON ISOLATION FOR QCD JETS
  for(pat::ElectronCollection::const_iterator electron = eleH->begin(); electron != eleH->end(); ++electron) {
    m_JetElectron_Pt        = electron->pt();
    m_JetElectron_Eta       = electron->eta();
    m_JetElectron_PFIso     = ElectronPFIso(electron, rho);
    m_JetElectron_CorrPFIso = ElectronCorrPFIso(electron, rho);
    TreeQCDElectron->Fill();
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


// ------------ method called once each job just before starting event loop ------------
void
Efficiency::beginJob()
{
  Service<TFileService> fs;

  TreeQCDMuon = fs->make<TTree>("TreeQCDMuon", "TreeQCDMuon");
  TreeQCDMuon->Branch("JetMuon_Pt",        &m_JetMuon_Pt,        "JetMuon_Pt/f");
  TreeQCDMuon->Branch("JetMuon_Eta",       &m_JetMuon_Eta,       "JetMuon_Eta/f");
  TreeQCDMuon->Branch("JetMuon_PFIso",     &m_JetMuon_PFIso,     "JetMuon_PFIso/f");
  TreeQCDMuon->Branch("JetMuon_CorrPFIso", &m_JetMuon_CorrPFIso, "JetMuon_CorrPFIso/f");

  TreeQCDElectron = fs->make<TTree>("TreeQCDElectron", "TreeQCDElectron");
  TreeQCDElectron->Branch("JetElectron_Pt",        &m_JetElectron_Pt,        "JetElectron_Pt/f");
  TreeQCDElectron->Branch("JetElectron_Eta",       &m_JetElectron_Eta,       "JetElectron_Eta/f");
  TreeQCDElectron->Branch("JetElectron_PFIso",     &m_JetElectron_PFIso,     "JetElectron_PFIso/f");
  TreeQCDElectron->Branch("JetElectron_CorrPFIso", &m_JetElectron_CorrPFIso, "JetElectron_CorrPFIso/f");

  Tree = fs->make<TTree>("Tree", "Tree");
  Tree->Branch("Nele", &m_Nele, "Nele/i");
  Tree->Branch("Nmuo", &m_Nmuo, "Nmuo/i");
  Tree->Branch("genJet", &m_genJet, "genJet/i");
  Tree->Branch("genTau_et", &m_genTau_et, "genTau_et/i");
  Tree->Branch("genTau_mt", &m_genTau_mt, "genTau_mt/i");
  Tree->Branch("genEle_et", &m_genEle_et, "genEle_et/i");
  Tree->Branch("genEle1_ee", &m_genEle1_ee, "genEle1_ee/i");
  Tree->Branch("genEle2_ee", &m_genEle2_ee, "genEle2_ee/i");
  Tree->Branch("genEle_em", &m_genEle_em, "genEle_em/i");
  Tree->Branch("genMuo_mt", &m_genMuo_mt, "genMuo_mt/i");
  Tree->Branch("genMuo1_mm", &m_genMuo1_mm, "genMuo1_mm/i");
  Tree->Branch("genMuo2_mm", &m_genMuo2_mm, "genMuo2_mm/i");
  Tree->Branch("genMuo_em", &m_genMuo_em, "genMuo_em/i");
  Tree->Branch("matchedJet", &m_matchedJet, "matchedJet/i");
  Tree->Branch("matchedTau_et", &m_matchedTau_et, "matchedTau_et/i");
  Tree->Branch("matchedTau_mt", &m_matchedTau_mt, "matchedTau_mt/i");
  Tree->Branch("matchedEle_et", &m_matchedEle_et, "matchedEle_et/i");
  Tree->Branch("matchedEle1_ee", &m_matchedEle1_ee, "matchedEle1_ee/i");
  Tree->Branch("matchedEle2_ee", &m_matchedEle2_ee, "matchedEle2_ee/i");
  Tree->Branch("matchedEle_em", &m_matchedEle_em, "matchedEle_em/i");
  Tree->Branch("matchedTightMuo_mt", &m_matchedTightMuo_mt, "matchedTightMuo_mt/i");
  Tree->Branch("matchedTightMuo1_mm", &m_matchedTightMuo1_mm, "matchedTightMuo1_mm/i");
  Tree->Branch("matchedTightMuo2_mm", &m_matchedTightMuo2_mm, "matchedTightMuo2_mm/i");
  Tree->Branch("matchedTightMuo_em", &m_matchedTightMuo_em, "matchedTightMuo_em/i");
  Tree->Branch("matchedHighptMuo_mt", &m_matchedHighptMuo_mt, "matchedHighptMuo_mt/i");
  Tree->Branch("matchedHighptMuo1_mm", &m_matchedHighptMuo1_mm, "matchedHighptMuo1_mm/i");
  Tree->Branch("matchedHighptMuo2_mm", &m_matchedHighptMuo2_mm, "matchedHighptMuo2_mm/i");
  Tree->Branch("matchedHighptMuo_em", &m_matchedHighptMuo_em, "matchedHighptMuo_em/i");
  Tree->Branch("matchedTrackerMuo1_mm", &m_matchedTrackerMuo1_mm, "matchedTrackerMuo1_mm/i");
  Tree->Branch("matchedTrackerMuo2_mm", &m_matchedTrackerMuo2_mm, "matchedTrackerMuo2_mm/i");
  Tree->Branch("recoJet", &m_recoJet, "recoJet/i");
  Tree->Branch("recoTau_mt", &m_recoTau_mt, "recoTau_mt/i");
  Tree->Branch("recoTau_et", &m_recoTau_et, "recoTau_et/i");
  Tree->Branch("recoEle_et", &m_recoEle_et, "recoEle_et/i");
  Tree->Branch("recoEle1_ee", &m_recoEle1_ee, "recoEle1_ee/i");
  Tree->Branch("recoEle2_ee", &m_recoEle2_ee, "recoEle2_ee/i");
  Tree->Branch("recoEle_em", &m_recoEle_em, "recoEle_em/i");
  Tree->Branch("recoTightMuo_mt",  &m_recoTightMuo_mt,  "recoTightMuo_mt/i");
  Tree->Branch("recoTightMuo1_mm", &m_recoTightMuo1_mm, "recoTightMuo1_mm/i");
  Tree->Branch("recoTightMuo2_mm", &m_recoTightMuo2_mm, "recoTightMuo2_mm/i");
  Tree->Branch("recoTightMuo_em",  &m_recoTightMuo_em,  "recoTightMuo_em/i");
  Tree->Branch("recoHighptMuo_mt",  &m_recoHighptMuo_mt,  "recoHighptMuo_mt/i");
  Tree->Branch("recoHighptMuo1_mm", &m_recoHighptMuo1_mm, "recoHighptMuo1_mm/i");
  Tree->Branch("recoHighptMuo2_mm", &m_recoHighptMuo2_mm, "recoHighptMuo2_mm/i");
  Tree->Branch("recoHighptMuo_em",  &m_recoHighptMuo_em,  "recoHighptMuo_em/i");
  Tree->Branch("recoTrackerMuo1_mm", &m_recoTrackerMuo1_mm, "recoTrackerMuo1_mm/i");
  Tree->Branch("recoTrackerMuo2_mm", &m_recoTrackerMuo2_mm, "recoTrackerMuo2_mm/i");
  Tree->Branch("hasAtLeastOneHighPtMuo", &m_hasAtLeastOneHighPtMuo, "hasAtLeastOneHighPtMuo/i");
  Tree->Branch("genJet_Pt", &m_genJet_Pt, "genJet_Pt/f");
  Tree->Branch("genJet_Eta", &m_genJet_Eta, "genJet_Eta/f");
  Tree->Branch("genTau_et_Pt", &m_genTau_et_Pt, "genTau_et_Pt/f");
  Tree->Branch("genTau_et_Eta", &m_genTau_et_Eta, "genTau_et_Eta/f");
  Tree->Branch("genTau_mt_Pt", &m_genTau_mt_Pt, "genTau_mt_Pt/f");
  Tree->Branch("genTau_mt_eta", &m_genTau_mt_eta, "genTau_mt_eta/f");
  Tree->Branch("genEle_et_Pt", &m_genEle_et_Pt, "genEle_et_Pt/f");
  Tree->Branch("genEle_et_Eta", &m_genEle_et_Eta, "genEle_et_Eta/f");
  Tree->Branch("genEle1_ee_Pt", &m_genEle1_ee_Pt, "genEle1_ee_Pt/f");
  Tree->Branch("genEle1_ee_Eta", &m_genEle1_ee_Eta, "genEle1_ee_Eta/f");
  Tree->Branch("genEle2_ee_Pt", &m_genEle2_ee_Pt, "genEle2_ee_Pt/f");
  Tree->Branch("genEle2_ee_Eta", &m_genEle2_ee_Eta, "genEle2_ee_Eta/f");
  Tree->Branch("genEle_em_Pt", &m_genEle_em_Pt, "genEle_em_Pt/f");
  Tree->Branch("genEle_em_Eta", &m_genEle_em_Eta, "genEle_em_Eta/f");
  Tree->Branch("genMuo_mt_Pt", &m_genMuo_mt_Pt, "genMuo_mt_Pt/f");
  Tree->Branch("genMuo_mt_Eta", &m_genMuo_mt_Eta, "genMuo_mt_Eta/f");
  Tree->Branch("genMuo1_mm_Pt", &m_genMuo1_mm_Pt, "genMuo1_mm_Pt/f");
  Tree->Branch("genMuo1_mm_Eta", &m_genMuo1_mm_Eta, "genMuo1_mm_Eta/f");
  Tree->Branch("genMuo2_mm_Pt", &m_genMuo2_mm_Pt, "genMuo2_mm_Pt/f");
  Tree->Branch("genMuo2_mm_Eta", &m_genMuo2_mm_Eta, "genMuo2_mm_Eta/f");
  Tree->Branch("genMuo_em_Pt", &m_genMuo_em_Pt, "genMuo_em_Pt/f");
  Tree->Branch("genMuo_em_Eta", &m_genMuo_em_Eta, "genMuo_em_Eta/f");
  Tree->Branch("recoJet_Pt", &m_recoJet_Pt, "recoJet_Pt/f");
  Tree->Branch("recoJet_Eta", &m_recoJet_Eta, "recoJet_Eta/f");
  Tree->Branch("recoTau_et_Pt", &m_recoTau_et_Pt, "recoTau_et_Pt/f");
  Tree->Branch("recoTau_et_Eta",&m_recoTau_et_Eta,"recoTau_et_Eta/f");
  Tree->Branch("recoTau_mt_Pt", &m_recoTau_mt_Pt, "recoTau_mt_Pt/f");
  Tree->Branch("recoTau_mt_eta",&m_recoTau_mt_eta,"recoTau_mt_eta/f");
  Tree->Branch("recoEle_et_Pt", &m_recoEle_et_Pt, "recoEle_et_Pt/f");
  Tree->Branch("recoEle_et_Eta", &m_recoEle_et_Eta, "recoEle_et_Eta/f");
  Tree->Branch("recoEle1_ee_Pt", &m_recoEle1_ee_Pt, "recoEle1_ee_Pt/f");
  Tree->Branch("recoEle1_ee_Eta", &m_recoEle1_ee_Eta, "recoEle1_ee_Eta/f");
  Tree->Branch("recoEle2_ee_Pt", &m_recoEle2_ee_Pt, "recoEle2_ee_Pt/f");
  Tree->Branch("recoEle2_ee_Eta", &m_recoEle2_ee_Eta, "recoEle2_ee_Eta/f");
  Tree->Branch("recoEle_em_Pt", &m_recoEle_em_Pt, "recoEle_em_Pt/f");
  Tree->Branch("recoEle_em_Eta", &m_recoEle_em_Eta, "recoEle_em_Eta/f");
  Tree->Branch("recoTightMuo_mt_Pt", &m_recoTightMuo_mt_Pt, "recoTightMuo_mt_Pt/f");
  Tree->Branch("recoTightMuo_mt_Eta", &m_recoTightMuo_mt_Eta, "recoTightMuo_mt_Eta/f");
  Tree->Branch("recoTightMuo1_mm_Pt", &m_recoTightMuo1_mm_Pt, "recoTightMuo1_mm_Pt/f");
  Tree->Branch("recoTightMuo1_mm_Eta", &m_recoTightMuo1_mm_Eta, "recoTightMuo1_mm_Eta/f");
  Tree->Branch("recoTightMuo2_mm_Pt", &m_recoTightMuo2_mm_Pt, "recoTightMuo2_mm_Pt/f");
  Tree->Branch("recoTightMuo2_mm_Eta", &m_recoTightMuo2_mm_Eta, "recoTightMuo2_mm_Eta/f");
  Tree->Branch("recoTightMuo_em_Pt", &m_recoTightMuo_em_Pt, "recoTightMuo_em_Pt/f");
  Tree->Branch("recoTightMuo_em_Eta", &m_recoTightMuo_em_Eta, "recoTightMuo_em_Eta/f");
  Tree->Branch("recoHighptMuo_mt_Pt", &m_recoHighptMuo_mt_Pt, "recoHighptMuo_mt_Pt/f");
  Tree->Branch("recoHighptMuo_mt_Eta", &m_recoHighptMuo_mt_Eta, "recoHighptMuo_mt_Eta/f");
  Tree->Branch("recoHighptMuo1_mm_Pt", &m_recoHighptMuo1_mm_Pt, "recoHighptMuo1_mm_Pt/f");
  Tree->Branch("recoHighptMuo1_mm_Eta", &m_recoHighptMuo1_mm_Eta, "recoHighptMuo1_mm_Eta/f");
  Tree->Branch("recoHighptMuo2_mm_Pt", &m_recoHighptMuo2_mm_Pt, "recoHighptMuo2_mm_Pt/f");
  Tree->Branch("recoHighptMuo2_mm_Eta", &m_recoHighptMuo2_mm_Eta, "recoHighptMuo2_mm_Eta/f");
  Tree->Branch("recoHighptMuo_em_Pt", &m_recoHighptMuo_em_Pt, "recoHighptMuo_em_Pt/f");
  Tree->Branch("recoHighptMuo_em_Eta", &m_recoHighptMuo_em_Eta, "recoHighptMuo_em_Eta/f");
  Tree->Branch("recoTrackerMuo1_mm_Pt", &m_recoTrackerMuo1_mm_Pt, "recoTrackerMuo1_mm_Pt/f");
  Tree->Branch("recoTrackerMuo1_mm_Eta", &m_recoTrackerMuo1_mm_Eta, "recoTrackerMuo1_mm_Eta/f");
  Tree->Branch("recoTrackerMuo2_mm_Pt", &m_recoTrackerMuo2_mm_Pt, "recoTrackerMuo2_mm_Pt/f");
  Tree->Branch("recoTrackerMuo2_mm_Eta", &m_recoTrackerMuo2_mm_Eta, "recoTrackerMuo2_mm_Eta/f");
  Tree->Branch("recoEle_et_PFIso", &m_recoEle_et_PFIso, "recoEle_et_PFIso/f");
  Tree->Branch("recoEle_et_CorrPFIso", &m_recoEle_et_CorrPFIso, "recoEle_et_CorrPFIso/f");
  Tree->Branch("recoEle1_ee_PFIso", &m_recoEle1_ee_PFIso, "recoEle1_ee_PFIso/f");
  Tree->Branch("recoEle1_ee_CorrPFIso", &m_recoEle1_ee_CorrPFIso, "recoEle1_ee_CorrPFIso/f");
  Tree->Branch("recoEle2_ee_PFIso", &m_recoEle2_ee_PFIso, "recoEle2_ee_PFIso/f");
  Tree->Branch("recoEle2_ee_CorrPFIso", &m_recoEle2_ee_CorrPFIso, "recoEle2_ee_CorrPFIso/f");
  Tree->Branch("recoEle_em_PFIso", &m_recoEle_em_PFIso, "recoEle_em_PFIso/f");
  Tree->Branch("recoEle_em_CorrPFIso", &m_recoEle_em_CorrPFIso, "recoEle_em_CorrPFIso/f");
  Tree->Branch("recoTightMuo_mt_PFIso", &m_recoTightMuo_mt_PFIso, "recoTightMuo_mt_PFIso/f");
  Tree->Branch("recoTightMuo_mt_CorrPFIso", &m_recoTightMuo_mt_CorrPFIso, "recoTightMuo_mt_CorrPFIso/f");
  Tree->Branch("recoTightMuo1_mm_PFIso", &m_recoTightMuo1_mm_PFIso, "recoTightMuo1_mm_PFIso/f");
  Tree->Branch("recoTightMuo1_mm_CorrPFIso", &m_recoTightMuo1_mm_CorrPFIso, "recoTightMuo1_mm_CorrPFIso/f");
  Tree->Branch("recoTightMuo2_mm_PFIso", &m_recoTightMuo2_mm_PFIso, "recoTightMuo2_mm_PFIso/f");
  Tree->Branch("recoTightMuo2_mm_CorrPFIso", &m_recoTightMuo2_mm_CorrPFIso, "recoTightMuo2_mm_CorrPFIso/f");
  Tree->Branch("recoTightMuo_em_PFIso", &m_recoTightMuo_em_PFIso, "recoTightMuo_em_PFIso/f");
  Tree->Branch("recoTightMuo_em_CorrPFIso", &m_recoTightMuo_em_CorrPFIso, "recoTightMuo_em_CorrPFIso/f");
  Tree->Branch("recoHighptMuo_mt_PFIso", &m_recoHighptMuo_mt_PFIso, "recoHighptMuo_mt_PFIso/f");
  Tree->Branch("recoHighptMuo_mt_CorrPFIso", &m_recoHighptMuo_mt_CorrPFIso, "recoHighptMuo_mt_CorrPFIso/f");
  Tree->Branch("recoHighptMuo1_mm_PFIso", &m_recoHighptMuo1_mm_PFIso, "recoHighptMuo1_mm_PFIso/f");
  Tree->Branch("recoHighptMuo1_mm_CorrPFIso", &m_recoHighptMuo1_mm_CorrPFIso, "recoHighptMuo1_mm_CorrPFIso/f");
  Tree->Branch("recoHighptMuo2_mm_PFIso", &m_recoHighptMuo2_mm_PFIso, "recoHighptMuo2_mm_PFIso/f");
  Tree->Branch("recoHighptMuo2_mm_CorrPFIso", &m_recoHighptMuo2_mm_CorrPFIso, "recoHighptMuo2_mm_CorrPFIso/f");
  Tree->Branch("recoHighptMuo_em_PFIso", &m_recoHighptMuo_em_PFIso, "recoHighptMuo_em_PFIso/f");
  Tree->Branch("recoHighptMuo_em_CorrPFIso", &m_recoHighptMuo_em_CorrPFIso, "recoHighptMuo_em_CorrPFIso/f");
  Tree->Branch("recoTrackerMuo1_mm_PFIso", &m_recoTrackerMuo1_mm_PFIso, "recoTrackerMuo1_mm_PFIso/f");
  Tree->Branch("recoTrackerMuo1_mm_CorrPFIso", &m_recoTrackerMuo1_mm_CorrPFIso, "recoTrackerMuo1_mm_CorrPFIso/f");
  Tree->Branch("recoTrackerMuo2_mm_PFIso", &m_recoTrackerMuo2_mm_PFIso, "recoTrackerMuo2_mm_PFIso/f");
  Tree->Branch("recoTrackerMuo2_mm_CorrPFIso", &m_recoTrackerMuo2_mm_CorrPFIso, "recoTrackerMuo2_mm_CorrPFIso/f");
  Tree->Branch("recoTightMuo1_mm_DetIso",   &m_recoTightMuo1_mm_DetIso,   "recoTightMuo1_mm_DetIso/f");
  Tree->Branch("recoTightMuo2_mm_DetIso",   &m_recoTightMuo2_mm_DetIso,   "recoTightMuo2_mm_DetIso/f");
  Tree->Branch("recoHighptMuo1_mm_DetIso",  &m_recoHighptMuo1_mm_DetIso,  "recoHighptMuo1_mm_DetIso/f");
  Tree->Branch("recoHighptMuo2_mm_DetIso",  &m_recoHighptMuo2_mm_DetIso,  "recoHighptMuo2_mm_DetIso/f");
  Tree->Branch("recoTrackerMuo1_mm_DetIso", &m_recoTrackerMuo1_mm_DetIso, "recoTrackerMuo1_mm_DetIso/f");
  Tree->Branch("recoTrackerMuo2_mm_DetIso", &m_recoTrackerMuo2_mm_DetIso, "recoTrackerMuo2_mm_DetIso/f");
  Tree->Branch("reco_em_deltaR", &m_reco_em_deltaR, "reco_em_deltaR/f");
  Tree->Branch("reco_et_deltaR", &m_reco_et_deltaR, "reco_et_deltaR/f");
  Tree->Branch("reco_ee_deltaR", &m_reco_ee_deltaR, "reco_ee_deltaR/f");
  Tree->Branch("reco_mt_deltaR", &m_reco_mt_deltaR, "reco_mt_deltaR/f");
  Tree->Branch("reco_mm_deltaR", &m_reco_mm_deltaR, "reco_mm_deltaR/f");
  Tree->Branch("dRJet", &m_dRJet, "dRJet/f");
  Tree->Branch("dRTau_et", &m_dRTau_et, "dRTau_et/f");
  Tree->Branch("dRTau_mt", &m_dRTau_mt, "dRTau_mt/f");
  Tree->Branch("dREle_et", &m_dREle_et, "dREle_et/f");
  Tree->Branch("dREle1_ee", &m_dREle1_ee, "dREle1_ee/f");
  Tree->Branch("dREle2_ee", &m_dREle2_ee, "dREle2_ee/f");
  Tree->Branch("dREle_em", &m_dREle_em, "dREle_em/f");
  Tree->Branch("dRMuo_mt", &m_dRMuo_mt, "dRMuo_mt/f");
  Tree->Branch("dRMuo1_mm", &m_dRMuo1_mm, "dRMuo1_mm/f");
  Tree->Branch("dRMuo2_mm", &m_dRMuo2_mm, "dRMuo2_mm/f");
  Tree->Branch("dRMuo_em", &m_dRMuo_em, "dRMuo_em/f");

}

// ------------ method called once each job just after ending the event loop ------------
void
Efficiency::endJob()
{
}

// ------------ method called when starting to processes a run ------------
void
Efficiency::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run ------------
void
Efficiency::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block ------------
void
Efficiency::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block ------------
void
Efficiency::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module ------------
void
Efficiency::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void Efficiency::SelectJet(edm::Handle<pat::JetCollection> CA8JetswithQjets,
			   edm::Handle<pat::JetCollection> CA8JetsPruned,
			   bool & foundSelectedJet,
			   pat::JetCollection::const_iterator & SelectedJet,
			   float & massZ, float & tau21Z, float & ptZ,
			   float massMin, float massMax,
			   bool foundJet, vector<reco::GenJet>::const_iterator GenJeT, bool & matched){

  if(foundJet){
    float dRGenReco = 9999.;
    pat::JetCollection::const_iterator SelectedJet_prov;
    for(pat::JetCollection::const_iterator jet = CA8JetswithQjets->begin(); jet != CA8JetswithQjets->end(); ++jet) {
      if(ROOT::Math::VectorUtil::DeltaR(jet->p4(),GenJeT->p4())<0.3 && ROOT::Math::VectorUtil::DeltaR(jet->p4(),GenJeT->p4())<dRGenReco){
	dRGenReco = ROOT::Math::VectorUtil::DeltaR(jet->p4(),GenJeT->p4());
	SelectedJet_prov = jet;
	matched = true;
      }
    }

    if(matched){
      float dRmin = 9999.; float mass = 0.;
      for(pat::JetCollection::const_iterator jetPruned = CA8JetsPruned->begin(); jetPruned != CA8JetsPruned->end(); ++jetPruned) {
	float dRtmp = ROOT::Math::VectorUtil::DeltaR(SelectedJet_prov->p4(),jetPruned->p4());
	if(dRtmp<dRmin && dRtmp<0.8 ){//matching failed if greater than jet radius
	  dRmin=dRtmp;
	  mass=jetPruned->mass();
	}
      }
      if(SelectedJet_prov->muonEnergyFraction()<0.99 && 
	 SelectedJet_prov->photonEnergyFraction()<0.99 &&
	 SelectedJet_prov->chargedEmEnergyFraction()<0.99 && 
	 SelectedJet_prov->neutralHadronEnergyFraction()<0.99 &&
	 SelectedJet_prov->chargedHadronEnergyFraction()>0.00 && 
	 SelectedJet_prov->pt()>400 && 
	 abs(SelectedJet_prov->eta())<2.4 && 
	 (mass>massMin && mass<massMax) && 
	 (SelectedJet_prov->userFloat("tau2")/SelectedJet_prov->userFloat("tau1"))<0.75){
	foundSelectedJet=true;
	if(SelectedJet_prov->pt()>ptZ){
	  massZ=mass;
	  ptZ=SelectedJet_prov->pt();
	  tau21Z=SelectedJet_prov->userFloat("tau2")/SelectedJet_prov->userFloat("tau1");
	  SelectedJet=SelectedJet_prov;
	}
      }
    }
  }
}

void Efficiency::SelectTau(edm::Handle<pat::TauCollection> tauHandle,
			   pat::JetCollection::const_iterator SelectedJet,
			   bool & foundTau,
			   pat::TauCollection::const_iterator & SelectedTau,
			   float & ptTau, bool foundJet,
			   bool foundPart, math::PtEtaPhiELorentzVector genPart, bool & matched){

  if(foundPart){
    float dRGenReco = 9999.;
    pat::TauCollection::const_iterator SelectedTau_prov;
    for (pat::TauCollection::const_iterator patTau = tauHandle->begin(); patTau != tauHandle->end(); ++patTau ) {
      if(ROOT::Math::VectorUtil::DeltaR(patTau->p4(),genPart)<0.3 && ROOT::Math::VectorUtil::DeltaR(patTau->p4(),genPart)<dRGenReco){
	dRGenReco = ROOT::Math::VectorUtil::DeltaR(patTau->p4(),genPart);
	SelectedTau_prov = patTau;
	matched = true;
      }
    }
    

    if(matched){
      if(SelectedTau_prov->pt()>20 && 
	 abs(SelectedTau_prov->eta())<2.4 &&
	 SelectedTau_prov->tauID("decayModeFinding")>0.5 && 
	 SelectedTau_prov->tauID("againstMuonLoose")>0.5 && 
	 SelectedTau_prov->tauID("againstElectronLoose")>0.5 &&
	 SelectedTau_prov->tauID("byVLooseCombinedIsolationDeltaBetaCorr")>0.5){
	if(foundJet){
	  if(ROOT::Math::VectorUtil::DeltaR(SelectedTau_prov->p4(),SelectedJet->p4())>0.8){
	    foundTau=true;
	    if(SelectedTau_prov->pt()>ptTau){
	      SelectedTau=SelectedTau_prov;
	      ptTau=SelectedTau_prov->pt();
	    }
	  }
	}
      }
    }
  }
}

void Efficiency::SelectHighptMuon(edm::Handle<pat::MuonCollection> muoH,
			    pat::JetCollection::const_iterator SelectedJet, bool & foundMuon,
			    pat::MuonCollection::const_iterator & SelectedMuon,
			    float & ptMuon, bool foundJet, reco::Vertex primaryVertex,
			    bool foundMuo, math::PtEtaPhiELorentzVector genMuo_mt, bool & matched){
  if(foundMuo){
    float dRGenReco = 9999.;
    pat::MuonCollection::const_iterator SelectedMuon_prov;
    for(pat::MuonCollection::const_iterator muon = muoH->begin(); muon != muoH->end(); ++muon) {
      if(ROOT::Math::VectorUtil::DeltaR(muon->p4(),genMuo_mt)<0.3 && ROOT::Math::VectorUtil::DeltaR(muon->p4(),genMuo_mt)<dRGenReco){
	dRGenReco = ROOT::Math::VectorUtil::DeltaR(muon->p4(),genMuo_mt);
	SelectedMuon_prov = muon;
	matched = true;
      }
    }


    if(matched){
      if(SelectedMuon_prov->isGlobalMuon()){
	reco::TrackRef cktTrack = (muon::tevOptimized(*SelectedMuon_prov, 200, 17., 40., 0.25)).first;
	if(cktTrack->pt()>10 &&
	   abs(cktTrack->eta())<2.4 && 
	   abs(cktTrack->phi())<3.2 && 
	   (cktTrack->ptError()/cktTrack->pt())<0.3 &&
	   SelectedMuon_prov->globalTrack()->hitPattern().numberOfValidMuonHits()>0 && 
	   SelectedMuon_prov->numberOfMatches()>1 && 
	   fabs(cktTrack->dxy(primaryVertex.position()))<0.2 && 
	   fabs(cktTrack->dz( primaryVertex.position()))<0.5 &&
	   SelectedMuon_prov->innerTrack()->hitPattern().numberOfValidPixelHits()>0 &&
	   SelectedMuon_prov->innerTrack()->hitPattern().trackerLayersWithMeasurement()>5){
	  if(foundJet) {
	    if(ROOT::Math::VectorUtil::DeltaR(SelectedMuon_prov->p4(),SelectedJet->p4())>0.8){
	      foundMuon=true;
	      if(cktTrack->pt()>ptMuon){
		SelectedMuon=SelectedMuon_prov;
		ptMuon=cktTrack->pt();
	      }
	    }
	  }
	}
      }
    }
  }
}
void Efficiency::SelectElectron(edm::Handle<pat::ElectronCollection> eleH,
				pat::JetCollection::const_iterator SelectedJet, bool & foundElectron,
				pat::ElectronCollection::const_iterator & SelectedElectron,
				float & ptElectron, bool foundJet, reco::Vertex primaryVertex,
				bool foundEle_et, math::PtEtaPhiELorentzVector genEle_et, bool & matched){
  if(foundEle_et){
    float dRGenReco = 9999.;
    pat::ElectronCollection::const_iterator SelectedElectron_prov;
    for(pat::ElectronCollection::const_iterator electron = eleH->begin(); electron != eleH->end(); ++electron) {
      if(ROOT::Math::VectorUtil::DeltaR(electron->p4(),genEle_et)<0.3 && ROOT::Math::VectorUtil::DeltaR(electron->p4(),genEle_et)<dRGenReco){
	dRGenReco = ROOT::Math::VectorUtil::DeltaR(electron->p4(),genEle_et);
	SelectedElectron_prov = electron;
	matched = true;
      }
    }
    

    if(matched){
      if(SelectedElectron_prov->pt()>10){
	if(fabs(SelectedElectron_prov->superCluster()->eta())<=1.479){
	  if(SelectedElectron_prov->deltaEtaSuperClusterTrackAtVtx()<0.004 &&
	     SelectedElectron_prov->deltaPhiSuperClusterTrackAtVtx()<0.030 &&
	     SelectedElectron_prov->sigmaIetaIeta()<0.01 &&
	     SelectedElectron_prov->hadronicOverEm()<0.12 &&
	     fabs(SelectedElectron_prov->gsfTrack()->dxy(primaryVertex.position()))<0.02 &&
	     fabs(SelectedElectron_prov->gsfTrack()->dz(primaryVertex.position()))<0.1 &&
	     (fabs(1/SelectedElectron_prov->ecalEnergy() - SelectedElectron_prov->eSuperClusterOverP()/SelectedElectron_prov->ecalEnergy()))<0.05 &&
	     SelectedElectron_prov->passConversionVeto()!=0 &&
	     SelectedElectron_prov->gsfTrack()->trackerExpectedHitsInner().numberOfHits()==0){
	    if(foundJet){
	      if(ROOT::Math::VectorUtil::DeltaR(SelectedElectron_prov->p4(),SelectedJet->p4())>0.8){
		foundElectron=true;
		if(SelectedElectron_prov->pt()>ptElectron){
		  SelectedElectron=SelectedElectron_prov;
		  ptElectron=SelectedElectron_prov->pt();
		}
	      }
	    }
	  }
	}
	if(fabs(SelectedElectron_prov->superCluster()->eta())>1.479 && fabs(SelectedElectron_prov->superCluster()->eta())<2.5){
	  if(SelectedElectron_prov->deltaEtaSuperClusterTrackAtVtx()<0.005  &&
	     SelectedElectron_prov->deltaPhiSuperClusterTrackAtVtx()<0.020  &&
	     SelectedElectron_prov->sigmaIetaIeta()<0.03  &&
	     SelectedElectron_prov->hadronicOverEm()<0.10  &&
	     fabs(SelectedElectron_prov->gsfTrack()->dxy(primaryVertex.position()))<0.02  &&
	     fabs(SelectedElectron_prov->gsfTrack()->dz(primaryVertex.position()))<0.1  &&
	     (fabs(1/SelectedElectron_prov->ecalEnergy() - SelectedElectron_prov->eSuperClusterOverP()/SelectedElectron_prov->ecalEnergy()))<0.05  &&
	     SelectedElectron_prov->passConversionVeto()!=0  &&
	     SelectedElectron_prov->gsfTrack()->trackerExpectedHitsInner().numberOfHits()==0){
	    if(foundJet){
	      if(ROOT::Math::VectorUtil::DeltaR(SelectedElectron_prov->p4(),SelectedJet->p4())>0.8){
		foundElectron=true;
		if(SelectedElectron_prov->pt()>ptElectron){
		  SelectedElectron=SelectedElectron_prov;
		  ptElectron=SelectedElectron_prov->pt();
		}
	      }
	    }
	  }
	}
      }
    }
  }
}


void Efficiency::SelectTrackerMuon(edm::Handle<pat::MuonCollection> muoH, bool & foundMuon, pat::MuonCollection::const_iterator & SelectedMuon, 
				   float & ptMuon, reco::Vertex primaryVertex, pat::JetCollection::const_iterator SelectedJet, bool foundJet, 
				   bool foundMuo, math::PtEtaPhiELorentzVector genMuo, bool & matched){
  if(foundMuo){
    float dRGenReco = 9999.;
    pat::MuonCollection::const_iterator SelectedMuon_prov;
    for(pat::MuonCollection::const_iterator muon = muoH->begin(); muon != muoH->end(); ++muon) {
      if(ROOT::Math::VectorUtil::DeltaR(muon->p4(),genMuo)<0.3 && ROOT::Math::VectorUtil::DeltaR(muon->p4(),genMuo)<dRGenReco){
	dRGenReco = ROOT::Math::VectorUtil::DeltaR(muon->p4(),genMuo);
	SelectedMuon_prov = muon;
	matched = true;
      }
    }
    

    if(matched){
      if(SelectedMuon_prov->pt()>10 && 
	 abs(SelectedMuon_prov->eta())<2.4 &&
	 abs(SelectedMuon_prov->phi())<3.2 &&
	 (SelectedMuon_prov->isTrackerMuon()) &&
	 SelectedMuon_prov->numberOfMatches()>1 &&
	 fabs(SelectedMuon_prov->muonBestTrack()->dz(primaryVertex.position()))<0.5 &&
	 fabs(SelectedMuon_prov->dB())<0.2  &&
	 SelectedMuon_prov->innerTrack()->hitPattern().numberOfValidPixelHits()>0 &&
	 SelectedMuon_prov->innerTrack()->hitPattern().trackerLayersWithMeasurement()>5 &&
	 (SelectedMuon_prov->muonBestTrack()->ptError()/SelectedMuon_prov->muonBestTrack()->pt())<=0.3){
	if(foundJet){
	  if(ROOT::Math::VectorUtil::DeltaR(SelectedMuon_prov->p4(),SelectedJet->p4())>0.8){
	    foundMuon=true;
	    if(SelectedMuon_prov->pt()>ptMuon){
	      SelectedMuon=SelectedMuon_prov;
	      ptMuon=SelectedMuon_prov->pt();
	    }
	  }
	}
      }
    }
  }
}


void Efficiency::SelectTightMuon(edm::Handle<pat::MuonCollection> muoH,
				 pat::JetCollection::const_iterator SelectedJet, bool & foundMuon,
				 pat::MuonCollection::const_iterator & SelectedMuon,
				 float & ptMuon, bool foundJet, reco::Vertex primaryVertex,
				 bool foundMuo, math::PtEtaPhiELorentzVector genMuo_mt, bool & matched){
  if(foundMuo){
    float dRGenReco = 9999.;
    pat::MuonCollection::const_iterator SelectedMuon_prov;
    for(pat::MuonCollection::const_iterator muon = muoH->begin(); muon != muoH->end(); ++muon) {
      if(ROOT::Math::VectorUtil::DeltaR(muon->p4(),genMuo_mt)<0.3 && ROOT::Math::VectorUtil::DeltaR(muon->p4(),genMuo_mt)<dRGenReco){
	dRGenReco = ROOT::Math::VectorUtil::DeltaR(muon->p4(),genMuo_mt);
	SelectedMuon_prov = muon;
	matched = true;
      }
    }
    if(matched){
      if(SelectedMuon_prov->pt()>10 && 
	 abs(SelectedMuon_prov->eta())<2.4 && 
	 abs(SelectedMuon_prov->phi())<3.2 && 
	 SelectedMuon_prov->isGlobalMuon() && 
	 SelectedMuon_prov->isPFMuon() && 
	 SelectedMuon_prov->globalTrack()->normalizedChi2()<10 && 
	 SelectedMuon_prov->globalTrack()->hitPattern().numberOfValidMuonHits()>0 && 
	 SelectedMuon_prov->numberOfMatches()>1 && 
	 fabs(SelectedMuon_prov->dB())<0.2  && 
	 fabs(SelectedMuon_prov->muonBestTrack()->dz(primaryVertex.position()))<0.5 && 
	 SelectedMuon_prov->innerTrack()->hitPattern().numberOfValidPixelHits()>0 && 
	 SelectedMuon_prov->innerTrack()->hitPattern().trackerLayersWithMeasurement()>5){
	if(foundJet) {
	  if(ROOT::Math::VectorUtil::DeltaR(SelectedMuon_prov->p4(),SelectedJet->p4())>0.8){
	    foundMuon=true;
	    if(SelectedMuon_prov->pt()>ptMuon){
	      SelectedMuon=SelectedMuon_prov;
	      ptMuon=SelectedMuon_prov->pt();
	    }
	  }
	}
      }
    }
  }
}

void Efficiency::SelectTrackerGlobalID(pat::MuonCollection::const_iterator SelectedMuo1, pat::MuonCollection::const_iterator SelectedMuo2, 
				       bool matched1, bool matched2, reco::Vertex primaryVertex, bool & hasAtLeastOneHighPtMuo){
  if(matched1 && matched2){
    vector<pat::MuonCollection::const_iterator> SelectedMuo;
    SelectedMuo.push_back(SelectedMuo1);
    SelectedMuo.push_back(SelectedMuo2);
    for(unsigned int i=0; i<SelectedMuo.size(); i++){
      if(!(SelectedMuo[i]->isGlobalMuon())) continue;
      reco::TrackRef cktTrack = (muon::tevOptimized(*(SelectedMuo[i]), 200, 17., 40., 0.25)).first;
      if(cktTrack->pt()<10) continue;
      if(abs(cktTrack->eta())>2.4) continue;
      if(abs(cktTrack->phi())>3.2) continue;
      if((cktTrack->ptError()/cktTrack->pt())>0.3) continue;
      if(SelectedMuo[i]->globalTrack()->hitPattern().numberOfValidMuonHits()<=0) continue;
      if(SelectedMuo[i]->numberOfMatches()<=1) continue;
      if(fabs(cktTrack->dxy(primaryVertex.position()))>=0.2) continue;
      if(fabs(cktTrack->dz( primaryVertex.position()))>=0.5) continue;
      if(SelectedMuo[i]->innerTrack()->hitPattern().numberOfValidPixelHits()<=0) continue;
      if(SelectedMuo[i]->innerTrack()->hitPattern().trackerLayersWithMeasurement()<=5) continue;
      hasAtLeastOneHighPtMuo = true;
    }
  }
}

float Efficiency::MuonPFIso(pat::MuonCollection::const_iterator muon, bool highpt){
  float sumChargedHadronPt = muon->pfIsolationR04().sumChargedHadronPt;
  float sumNeutralHadronEt = muon->pfIsolationR04().sumNeutralHadronEt;
  float sumPhotonEt = muon->pfIsolationR04().sumPhotonEt;
  float sumPUPt = muon->pfIsolationR04().sumPUPt;
  float iso = (sumChargedHadronPt+ max(0.,sumNeutralHadronEt+sumPhotonEt-0.5*sumPUPt))/muon->pt();
  if(highpt && muon->isGlobalMuon()){
    reco::TrackRef cktTrack = (muon::tevOptimized(*muon, 200, 17., 40., 0.25)).first;
    iso = (sumChargedHadronPt+ max(0.,sumNeutralHadronEt+sumPhotonEt-0.5*sumPUPt))/cktTrack->pt();
  }
  return iso;
}


float Efficiency::ElectronPFIso(pat::ElectronCollection::const_iterator electron, float rho){
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

float Efficiency::MuonCorrPFIso(pat::MuonCollection::const_iterator muon, bool highpt){
  float sumChargedHadronPt = muon->userIsolation(pat::PfChargedHadronIso);
  float sumNeutralHadronEt = muon->userIsolation(pat::PfNeutralHadronIso);
  float sumPhotonEt = muon->userIsolation(pat::PfGammaIso);
  float sumPUPt = muon->userIsolation(pat::User2Iso);
  float iso = (sumChargedHadronPt+ max(0.,sumNeutralHadronEt+sumPhotonEt-0.5*sumPUPt))/muon->pt();
  if(highpt && muon->isGlobalMuon()){
    reco::TrackRef cktTrack = (muon::tevOptimized(*muon, 200, 17., 40., 0.25)).first;
    iso = (sumChargedHadronPt+ max(0.,sumNeutralHadronEt+sumPhotonEt-0.5*sumPUPt))/cktTrack->pt();
  }
  return iso;
}


float Efficiency::ElectronCorrPFIso(pat::ElectronCollection::const_iterator electron, float rho){
  float chargedHadronIso = electron->userIsolation(pat::PfChargedHadronIso);
  float neutralHadronIso = electron->userIsolation(pat::PfNeutralHadronIso);
  float photonIso = electron->userIsolation(pat::PfGammaIso);
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

float Efficiency::MuonDETIso(pat::MuonCollection::const_iterator SelectedMuo1, pat::MuonCollection::const_iterator SelectedMuo2, bool secondMuon){
  float isovar = SelectedMuo1->trackIso()/SelectedMuo1->pt();
  if(secondMuon){
    double dR = ROOT::Math::VectorUtil::DeltaR(SelectedMuo1->p4(), SelectedMuo2->p4());
    if(dR < 0.3) isovar = isovar - SelectedMuo2->track()->pt()/SelectedMuo1->pt();
  }
  return isovar;
}

//define this as a plug-in
DEFINE_FWK_MODULE(Efficiency);
