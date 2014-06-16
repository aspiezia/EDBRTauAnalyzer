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
  void SelectTau(edm::Handle<pat::TauCollection> tauHandle, pat::JetCollection::const_iterator SelectedJet,
		 pat::TauCollection::const_iterator & SelectedTau, bool foundJet,
		 bool foundPart, math::PtEtaPhiELorentzVector genPart, bool & matched);
  void SelectElectron(edm::Handle<pat::ElectronCollection> eleH, pat::JetCollection::const_iterator SelectedJet,
		      pat::ElectronCollection::const_iterator & SelectedEle,
		      bool foundEle_et, math::PtEtaPhiELorentzVector genEle_et, bool & matched, bool & isoHEEP, float rho);
  void SelectMuon(edm::Handle<pat::MuonCollection> muoH, pat::JetCollection::const_iterator SelectedJet,
		  pat::MuonCollection::const_iterator & SelectedMuon,bool foundMuo, math::PtEtaPhiELorentzVector genMuo_mt, bool & matched);
  
  float MuonPFIso(pat::MuonCollection::const_iterator muon, bool highpt);
  float ElectronPFIso(pat::ElectronCollection::const_iterator electron, float rho);
  float MuonCorrPFIso(pat::MuonCollection::const_iterator muon, bool highpt);
  float ElectronCorrPFIso(pat::ElectronCollection::const_iterator electron, float rho);
  float MuonDETIso(pat::MuonCollection::const_iterator SelectedMuo1, pat::MuonCollection::const_iterator SelectedMuo2, bool secondMuon);
  bool ElectronDETIso(pat::ElectronCollection::const_iterator electron, float rho);

  void FillTree(int category, TTree *Tree, pat::JetCollection::const_iterator SelectedJet, pat::TauCollection::const_iterator SelectedTau,
		pat::MuonCollection::const_iterator SelectedMuon, pat::MuonCollection::const_iterator SelectedMuon1, pat::MuonCollection::const_iterator SelectedMuon2,
		pat::ElectronCollection::const_iterator SelectedElectron, pat::ElectronCollection::const_iterator SelectedElectron1, 
		pat::ElectronCollection::const_iterator SelectedElectron2, reco::Vertex primaryVertex, math::PtEtaPhiELorentzVector genLep1, 
		math::PtEtaPhiELorentzVector genLep2, std::vector<reco::GenJet>::const_iterator genJet, float rho, bool isoHEEP1, bool isoHEEP2,
		bool foundGenJet, bool foundSelectedJet, bool foundGenLep1, bool foundSelectedLep1,bool foundGenLep2, bool foundSelectedLep2);
  

  TTree *TreeEleMuo; 
  TTree *TreeMuoMuo; 
  TTree *TreeEleEle; 
  TTree *TreeMuoTau; 
  TTree *TreeEleTau; 
  TTree *TreeQCDMuon; 
  TTree *TreeQCDElectron; 

  float m_genJetPt;
  float m_genJetEta;
  float m_genJetPhi;
  float m_genLep1Pt;
  float m_genLep2Pt;
  float m_genLep1Eta;
  float m_genLep2Eta;
  float m_genLep1Phi;
  float m_genLep2Phi;
  float m_recoJetPt;
  float m_recoJetEta;
  float m_recoJetPhi;
  float m_recoDeltaR;
  float m_dRGenRecoLep1;	
  float m_dRGenRecoLep2;
  float m_genDeltaR;
  float m_dRGenRecoJet;
  float m_recoLep1dRJetF;
  float m_recoLep2dRJetF;
  float m_recoLep1PtF;
  float m_recoLep1EtaF;
  float m_recoLep1PhiF;
  float m_recoLep2PtF;
  float m_recoLep2EtaF;
  float m_recoLep2PhiF;
  int m_recoLep1Pt;
  int m_recoLep1Eta;
  int m_recoLep1Phi;
  int m_recoLep1HEEP;
  int m_recoLep1dRJet;
  int m_recoLep1PFIso;
  int m_recoLep1CorrPFIso;
  int m_recoLep1DETIso;
  int m_recoLep1HEEPIso;
  int m_recoLep1dEtaIn;
  int m_recoLep1dPhiIn;
  int m_recoLep1SigmaIEIE;
  int m_recoLep1HoverE;
  int m_recoLep1dxy;
  int m_recoLep1dz;
  int m_recoLep1EandP;
  int m_recoLep1Conv1;
  int m_recoLep1Conv2;
  int m_recoLep1EcalDriven;
  int m_recoLep1dPhiInHEEP;
  int m_recoLep1HoverEHEEP;
  int m_recoLep1TkSumPt;
  int m_recoLep1LostHits;
  int m_recoLep1dEtaInHEEP;
  int m_recoLep1dxyHEEP;
  int m_recoLep1e5xe5;
  int m_recoLep1SigmaIEIEHEEP;
  int m_recoLep1IsPFTight;
  int m_recoLep1Chi2Tight;
  int m_recoLep1MuonHitsTight;
  int m_recoLep1MatchesTight;
  int m_recoLep1dxyTight;
  int m_recoLep1dzTight;
  int m_recoLep1PixelHitsTight;
  int m_recoLep1TrackerLTight;
  int m_recoLep1IsGlobal;
  int m_recoLep1ptErr;
  int m_recoLep1MuonHits;
  int m_recoLep1Matches;
  int m_recoLep1PixelHits;
  int m_recoLep1TrackerL;
  int m_recoLep1IsPFMuon;
  int m_recoLep1Chi2;
  int m_recoLep1Discr1;
  int m_recoLep1Discr2;
  int m_recoLep1Discr3;
  int m_recoLep1Discr4;
  int m_recoLep2Pt;
  int m_recoLep2Eta;
  int m_recoLep2Phi;
  int m_recoLep2HEEP;
  int m_recoLep2dRJet;
  int m_recoLep2PFIso;
  int m_recoLep2CorrPFIso;
  int m_recoLep2DETIso;
  int m_recoLep2HEEPIso;
  int m_recoLep2dEtaIn;
  int m_recoLep2dPhiIn;
  int m_recoLep2SigmaIEIE;
  int m_recoLep2HoverE;
  int m_recoLep2dxy;
  int m_recoLep2dz;
  int m_recoLep2EandP;
  int m_recoLep2Conv1;
  int m_recoLep2Conv2;
  int m_recoLep2EcalDriven;
  int m_recoLep2dPhiInHEEP;
  int m_recoLep2HoverEHEEP;
  int m_recoLep2TkSumPt;
  int m_recoLep2LostHits;
  int m_recoLep2dEtaInHEEP;
  int m_recoLep2dxyHEEP;
  int m_recoLep2e5xe5;
  int m_recoLep2SigmaIEIEHEEP;
  int m_recoLep2IsPFTight;
  int m_recoLep2Chi2Tight;
  int m_recoLep2MuonHitsTight;
  int m_recoLep2MatchesTight;
  int m_recoLep2dxyTight;
  int m_recoLep2dzTight;
  int m_recoLep2PixelHitsTight;
  int m_recoLep2TrackerLTight;
  int m_recoLep2IsGlobal;
  int m_recoLep2ptErr;
  int m_recoLep2MuonHits;
  int m_recoLep2Matches;
  int m_recoLep2PixelHits;
  int m_recoLep2TrackerL;
  int m_recoLep2IsPFMuon;
  int m_recoLep2Chi2;
 
  float m_JetMuon_Pt;
  float m_JetMuon_Eta;
  float m_JetMuon_PFIso;
  float m_JetMuon_CorrPFIso;
  float m_JetElectron_Pt;
  float m_JetElectron_Eta;
  float m_JetElectron_PFIso;
  float m_JetElectron_CorrPFIso;

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

  math::PtEtaPhiELorentzVector genZ;       bool foundGenZ    = false;
  math::PtEtaPhiELorentzVector genTau_et;  bool foundGenTau_et  = false;
  math::PtEtaPhiELorentzVector genTau_mt;  bool foundGenTau_mt  = false;
  math::PtEtaPhiELorentzVector genMuo_mt;  bool foundGenMuo_mt  = false;
  math::PtEtaPhiELorentzVector genMuo1_mm; bool foundGenMuo1_mm = false;
  math::PtEtaPhiELorentzVector genMuo2_mm; bool foundGenMuo2_mm = false;
  math::PtEtaPhiELorentzVector genMuo_em;  bool foundGenMuo_em = false;
  math::PtEtaPhiELorentzVector genEle_et;  bool foundGenEle_et  = false;
  math::PtEtaPhiELorentzVector genEle1_ee; bool foundGenEle1_ee = false;
  math::PtEtaPhiELorentzVector genEle2_ee; bool foundGenEle2_ee = false;
  math::PtEtaPhiELorentzVector genEle_em;  bool foundGenEle_em = false;
  math::PtEtaPhiELorentzVector genEle;     bool foundGenEle = false;
  math::PtEtaPhiELorentzVector genMuo;     bool foundGenMuo = false;
  for(size_t ngenPart=0; ngenPart<genParts->size(); ngenPart++){
    const reco::GenParticle & genPart = (*genParts)[ngenPart];
    if(abs(genPart.pdgId())==2212) continue;
    const reco::Candidate * mom = genPart.mother();
    //cout<<genPart.pdgId()<<" "<<genPart.status()<<" "<<mom->pdgId()<<endl;
  
    //LOOK FOR THE ELECTRON - DRELL-YAN
    if(abs(genPart.pdgId())==11 && genPart.status()!=3 && (abs(mom->pdgId())==23 || abs(mom->pdgId())==11)) {
      genEle = genPart.p4();
      foundGenEle = true;
    }
  
    //LOOK FOR THE MUON - DRELL-YAN
    if(abs(genPart.pdgId())==13 && genPart.status()!=3 && (abs(mom->pdgId())==23 || abs(mom->pdgId())==13)) {
      genMuo = genPart.p4();
      foundGenMuo = true;
    }
  
    //LOOK FOR THE Z
    if(abs(genPart.pdgId())==23 && genPart.status()!=3) {
      genZ = genPart.p4();
      foundGenZ = true;
    }
  
    //LOOK FOR THE HADRONIC TAU - MUOTAU
    int tauChargeMT=0;
    if((ele==0 && muo==1) && abs(genPart.pdgId())==15 && genPart.status()!=3 && abs(mom->pdgId())==15 && genPart.pt()>20 && fabs(genPart.eta())<2.3){
	
      for(unsigned int ndaugh=0; ndaugh<genPart.numberOfDaughters(); ndaugh++){
	const reco::Candidate * daughter = genPart.daughter(ndaugh);
	if(abs(daughter->pdgId())!=11 && abs(daughter->pdgId())!=12 && abs(daughter->pdgId())!=13 && abs(daughter->pdgId())!=14 
	   && abs(daughter->pdgId())!=15 && abs(daughter->pdgId())!=16 && abs(daughter->pdgId())!=22){
	  genTau_mt=genPart.p4();
	  tauChargeMT=genPart.pdgId();
	  foundGenTau_mt = true;
	}
      }
      for(unsigned int ndaugh=0; ndaugh<genPart.numberOfDaughters(); ndaugh++){
	const reco::Candidate * daughter = genPart.daughter(ndaugh);
	if((abs(daughter->pdgId())==12 || abs(daughter->pdgId())==14 ||abs(daughter->pdgId())==16) && daughter->status()==1){
	  if(foundGenTau_mt && genPart.pdgId()==tauChargeMT) genTau_mt = genTau_mt - daughter->p4();
	}
      }
    }
  
    //LOOK FOR THE HADRONIC TAU - ELETAU
    int tauChargeET=0;
    if((ele==1 && muo==0) && abs(genPart.pdgId())==15 && genPart.status()!=3 && abs(mom->pdgId())==15 && genPart.pt()>20 && fabs(genPart.eta())<2.3){
      for(unsigned int ndaugh=0; ndaugh<genPart.numberOfDaughters(); ndaugh++){
	const reco::Candidate * daughter = genPart.daughter(ndaugh);
	if(abs(daughter->pdgId())!=11 && abs(daughter->pdgId())!=12 && abs(daughter->pdgId())!=13 && abs(daughter->pdgId())!=14 
	   && abs(daughter->pdgId())!=15 && abs(daughter->pdgId())!=16 && abs(daughter->pdgId())!=22){
	  genTau_et=genPart.p4();
	  tauChargeET=genPart.pdgId();
	  foundGenTau_et = true;
	}
      }
      for(unsigned int ndaugh=0; ndaugh<genPart.numberOfDaughters(); ndaugh++){
	const reco::Candidate * daughter = genPart.daughter(ndaugh);
	if((abs(daughter->pdgId())==12 || abs(daughter->pdgId())==14 ||abs(daughter->pdgId())==16) && daughter->status()==1){
	  if(foundGenTau_et && genPart.pdgId()==tauChargeET) genTau_et = genTau_et - daughter->p4();
	}
      }
    }

    //LOOK FOR THE ELECTRON - ELETAU
    if((ele==1 && muo==0) && abs(genPart.pdgId())==11 && genPart.status()==1 && abs(mom->pdgId())==15 && genPart.pt()>10 && fabs(genPart.eta())<2.5){
      genEle_et = genPart.p4();
      foundGenEle_et = true;
    }

    //LOOK FOR THE MUON - MUOTAU
    if((ele==0 && muo==1) && abs(genPart.pdgId())==13 && genPart.status()==1 && abs(mom->pdgId())==15 && genPart.pt()>10 && fabs(genPart.eta())<2.4){
      genMuo_mt = genPart.p4();
      foundGenMuo_mt = true;
    }
  
    //LOOK FOR THE ELECTRONS - ELEELE
    if((ele==2 && muo==0) && genPart.pdgId()==11  && genPart.status()==1 && abs(mom->pdgId())==15 && genPart.pt()>10 && fabs(genPart.eta())<2.5){
      genEle1_ee = genPart.p4();
      foundGenEle1_ee = true;
    }
    if((ele==2 && muo==0) && genPart.pdgId()==-11 && genPart.status()==1 && abs(mom->pdgId())==15 && genPart.pt()>10 && fabs(genPart.eta())<2.5){
      genEle2_ee = genPart.p4();
      foundGenEle2_ee = true;
    }

    //LOOK FOR THE MUONS - MUOMUO
    if((ele==0 && muo==2) && genPart.pdgId()==13  && genPart.status()==1 && abs(mom->pdgId())==15 && genPart.pt()>10 && fabs(genPart.eta())<2.4){
      genMuo1_mm = genPart.p4();
      foundGenMuo1_mm = true;
    }
    if((ele==0 && muo==2) && genPart.pdgId()==-13 && genPart.status()==1 && abs(mom->pdgId())==15 && genPart.pt()>10 && fabs(genPart.eta())<2.4){
      genMuo2_mm = genPart.p4();
      foundGenMuo2_mm = true;
    }

    //LOOK FOR THE MUON - ELEMUO
    if((ele==1 && muo==1) && abs(genPart.pdgId())==13  && genPart.status()==1 && abs(mom->pdgId())==15 && genPart.pt()>10 && fabs(genPart.eta())<2.4){
      genMuo_em = genPart.p4();
      foundGenMuo_em = true;
    }

    //LOOK FOR THE ELECTRON - ELEMUO
    if((ele==1 && muo==1) && abs(genPart.pdgId())==11 && genPart.status()==1 && abs(mom->pdgId())==15 && genPart.pt()>10 && fabs(genPart.eta())<2.5){
      genEle_em = genPart.p4();
      foundGenEle_em = true;
    }
  }

  //MATCH THE GEN Z TO THE GEN JET
  vector<reco::GenJet>::const_iterator genJet; bool foundGenJet=false;
  for(vector<reco::GenJet>::const_iterator genjet=genjets->begin(); genjet!=genjets->end(); genjet++){
    if(foundGenZ) {if(ROOT::Math::VectorUtil::DeltaR(genZ,genjet->p4())<0.3) {genJet=genjet; foundGenJet=true;}}
  }

  //ORDER PT OF LEPTONS
  math::PtEtaPhiELorentzVector genEle_prov;
  if(foundGenEle1_ee && foundGenEle2_ee){
    if(genEle1_ee.pt()<genEle2_ee.pt()){
      genEle_prov=genEle1_ee;
      genEle1_ee=genEle2_ee;
      genEle2_ee=genEle_prov;
    }
  }
  math::PtEtaPhiELorentzVector genMuo_prov;
  if(foundGenMuo1_mm && foundGenMuo2_mm){
    if(genMuo1_mm.pt()<genMuo2_mm.pt()){
      genMuo_prov=genMuo1_mm;
      genMuo1_mm=genMuo2_mm;
      genMuo2_mm=genMuo_prov;
    }
  }

  
  //JET SELECTION
  pat::JetCollection::const_iterator SelectedJet; bool matchedJet = false;
  float massZ=-9999; float ptZ=-999; bool foundSelectedJet=false; float tau21Z=-9999;
  SelectJet(CA8JetswithQjets, CA8JetsPruned, foundSelectedJet, SelectedJet, massZ, tau21Z, ptZ, 70, 110, foundGenJet, genJet, matchedJet);

  //TAU SELECTION - MT
  pat::TauCollection::const_iterator SelectedTau_mt; bool matchedTau_mt = false;
  SelectTau(tauMuTauHandle, SelectedJet, SelectedTau_mt, foundSelectedJet, foundGenTau_mt, genTau_mt, matchedTau_mt);

  //TAU SELECTION - ET
  pat::TauCollection::const_iterator SelectedTau_et; bool matchedTau_et = false;
  SelectTau(tauElTauHandle, SelectedJet, SelectedTau_et, foundSelectedJet, foundGenTau_et, genTau_et, matchedTau_et); 

  //ELECTRON SELECTION - ET
  pat::ElectronCollection::const_iterator SelectedElectron_et; bool matchedEle_et = false; bool isoHEEP_et=false;
  SelectElectron(eleH,SelectedJet,SelectedElectron_et,foundGenEle_et,genEle_et,matchedEle_et,isoHEEP_et,rho);

  //ELECTRON SELECTION - EE1
  pat::ElectronCollection::const_iterator SelectedEle1_ee; bool matchedEle1_ee = false; bool isoHEEP_ee1=false;
  SelectElectron(eleH,SelectedJet,SelectedEle1_ee,foundGenEle1_ee,genEle1_ee,matchedEle1_ee,isoHEEP_ee1,rho);

  //ELECTRON SELECTION - EE2
  pat::ElectronCollection::const_iterator SelectedEle2_ee; bool matchedEle2_ee = false; bool isoHEEP_ee2=false;
  SelectElectron(eleH,SelectedJet,SelectedEle2_ee,foundGenEle2_ee,genEle2_ee,matchedEle2_ee,isoHEEP_ee2,rho);

  //ELECTRON SELECTION - EM
  pat::ElectronCollection::const_iterator SelectedElectron_em; bool matchedEle_em = false; bool isoHEEP_em=false;
  SelectElectron(eleH,SelectedJet,SelectedElectron_em,foundGenEle_em,genEle_em,matchedEle_em,isoHEEP_em,rho);

  //MUON SELECTION - MT
  pat::MuonCollection::const_iterator SelectedMuon_mt; bool matchedMuo_mt = false;
  SelectMuon(muoH,SelectedJet,SelectedMuon_mt,foundGenMuo_mt,genMuo_mt,matchedMuo_mt);

  //MUON SELECTION - MM1
  pat::MuonCollection::const_iterator SelectedMuon1_mm; bool matchedMuo1_mm = false;
  SelectMuon(muoH,SelectedJet,SelectedMuon1_mm,foundGenMuo1_mm,genMuo1_mm,matchedMuo1_mm);

  //MUON SELECTION - MM2
  pat::MuonCollection::const_iterator SelectedMuon2_mm; bool matchedMuo2_mm = false;
  SelectMuon(muoH,SelectedJet,SelectedMuon2_mm,foundGenMuo2_mm,genMuo2_mm,matchedMuo2_mm);

  //MUON SELECTION - EM
  pat::MuonCollection::const_iterator SelectedMuon_em; bool matchedMuo_em = false;
  SelectMuon(muoH,SelectedJet,SelectedMuon_em,foundGenMuo_em,genMuo_em,matchedMuo_em);

  pat::TauCollection::const_iterator SelectedTauFake;
  pat::MuonCollection::const_iterator SelectedMuonFake;
  pat::MuonCollection::const_iterator SelectedMuon1Fake;
  pat::MuonCollection::const_iterator SelectedMuon2Fake;
  pat::ElectronCollection::const_iterator SelectedElectronFake;
  pat::ElectronCollection::const_iterator SelectedElectron1Fake;
  pat::ElectronCollection::const_iterator SelectedElectron2Fake;
  bool isoHEEPFake=false;
  if(ele==1 && muo==1) FillTree(0,TreeEleMuo,SelectedJet,SelectedTauFake,SelectedMuon_em,SelectedMuon1Fake,SelectedMuon2Fake,SelectedElectron_em,SelectedElectron1Fake,
				SelectedElectron2Fake,primaryVertex,genEle_em,genMuo_em,genJet,rho,isoHEEP_em,isoHEEPFake,foundGenJet,foundSelectedJet,foundGenEle_em, 
				matchedEle_em,foundGenMuo_em, matchedMuo_em);
  if(ele==0 && muo==2) FillTree(1,TreeMuoMuo,SelectedJet,SelectedTauFake,SelectedMuonFake,SelectedMuon1_mm,SelectedMuon2_mm,SelectedElectronFake,SelectedElectron1Fake,
				SelectedElectron2Fake,primaryVertex,genMuo1_mm,genMuo2_mm,genJet,rho,isoHEEPFake,isoHEEPFake,foundGenJet,foundSelectedJet,foundGenMuo1_mm, 
				matchedMuo1_mm,foundGenMuo2_mm, matchedMuo2_mm);
  if(ele==2 && muo==0) FillTree(2,TreeEleEle,SelectedJet,SelectedTauFake,SelectedMuonFake,SelectedMuon1Fake,SelectedMuon2Fake,SelectedElectronFake,SelectedEle1_ee,
				SelectedEle2_ee,primaryVertex,genEle1_ee,genEle2_ee,genJet,rho,isoHEEP_ee1,isoHEEP_ee2,foundGenJet,foundSelectedJet,foundGenEle1_ee, 
				matchedEle1_ee,foundGenEle2_ee, matchedEle2_ee);
  if(ele==0 && muo==1) FillTree(3,TreeMuoTau,SelectedJet,SelectedTau_mt,SelectedMuon_mt,SelectedMuon1Fake,SelectedMuon2Fake,SelectedElectronFake,SelectedElectron1Fake,
				SelectedElectron2Fake,primaryVertex,genTau_mt,genMuo_mt,genJet,rho,isoHEEPFake,isoHEEPFake,foundGenJet,foundSelectedJet,foundGenTau_mt, 
				matchedTau_mt,foundGenMuo_mt, matchedMuo_mt);
  if(ele==1 && muo==0) FillTree(4,TreeEleTau,SelectedJet,SelectedTau_et,SelectedMuonFake,SelectedMuon1Fake,SelectedMuon2Fake,SelectedElectron_et,SelectedElectron1Fake,
				SelectedElectron2Fake,primaryVertex,genTau_et,genEle_et,genJet,rho,isoHEEP_et,isoHEEP_et,foundGenJet,foundSelectedJet,foundGenTau_et, 
				matchedTau_et,foundGenEle_et, matchedEle_et);

  //ELECTRON SELECTION - DRELL-YAN
  pat::ElectronCollection::const_iterator SelectedEle; bool matchedEle = false; bool isoHEEP=false;
  SelectElectron(eleH,SelectedJet,SelectedEle,foundGenEle,genEle,matchedEle,isoHEEP,rho);
  //MUON SELECTION - DRELL-YAN
  pat::MuonCollection::const_iterator SelectedMuon; bool matchedMuo = false;
  SelectMuon(muoH,SelectedJet,SelectedMuon,foundGenMuo,genMuo,matchedMuo);

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

  TreeEleMuo = fs->make<TTree>("TreeEleMuo", "TreeEleMuo");
  TreeMuoMuo = fs->make<TTree>("TreeMuoMuo", "TreeMuoMuo");
  TreeEleEle = fs->make<TTree>("TreeEleEle", "TreeEleEle");
  TreeMuoTau = fs->make<TTree>("TreeMuoTau", "TreeMuoTau");
  TreeEleTau = fs->make<TTree>("TreeEleTau", "TreeEleTau");

  TreeEleMuo->Branch("genJetPt", &m_genJetPt,"genJetPt/f");
  TreeEleMuo->Branch("genJetEta", &m_genJetEta,"genJetEta/f");
  TreeEleMuo->Branch("genJetPhi", &m_genJetPhi,"genJetPhi/f");
  TreeEleMuo->Branch("genLep1Pt", &m_genLep1Pt,"genLep1Pt/f");
  TreeEleMuo->Branch("genLep2Pt", &m_genLep2Pt,"genLep2Pt/f");
  TreeEleMuo->Branch("genLep1Eta", &m_genLep1Eta,"genLep1Eta/f");
  TreeEleMuo->Branch("genLep2Eta", &m_genLep2Eta,"genLep2Eta/f");
  TreeEleMuo->Branch("genLep1Phi", &m_genLep1Phi,"genLep1Phi/f");
  TreeEleMuo->Branch("genLep2Phi", &m_genLep2Phi,"genLep2Phi/f");
  TreeEleMuo->Branch("recoJetPt", &m_recoJetPt,"recoJetPt/f");
  TreeEleMuo->Branch("recoJetEta", &m_recoJetEta,"recoJetEta/f");
  TreeEleMuo->Branch("recoJetPhi", &m_recoJetPhi,"recoJetPhi/f");
  TreeEleMuo->Branch("recoDeltaR", &m_recoDeltaR,"recoDeltaR/f");
  TreeEleMuo->Branch("dRGenRecoLep1", &m_dRGenRecoLep1,"dRGenRecoLep1/f");
  TreeEleMuo->Branch("dRGenRecoLep2", &m_dRGenRecoLep2,"dRGenRecoLep2/f");
  TreeEleMuo->Branch("genDeltaR", &m_genDeltaR,"genDeltaR/f");
  TreeEleMuo->Branch("dRGenRecoJet", &m_dRGenRecoJet,"dRGenRecoJet/f");
  TreeEleMuo->Branch("recoLep1dRJetF", &m_recoLep1dRJetF,"recoLep1dRJetF/f");
  TreeEleMuo->Branch("recoLep2dRJetF", &m_recoLep2dRJetF,"recoLep2dRJetF/f");
  TreeEleMuo->Branch("recoLep1PtF", &m_recoLep1PtF,"recoLep1PtF/f");
  TreeEleMuo->Branch("recoLep1EtaF", &m_recoLep1EtaF,"recoLep1EtaF/f");
  TreeEleMuo->Branch("recoLep1PhiF", &m_recoLep1PhiF,"recoLep1PhiF/f");
  TreeEleMuo->Branch("recoLep2PtF", &m_recoLep2PtF,"recoLep2PtF/f");
  TreeEleMuo->Branch("recoLep2EtaF", &m_recoLep2EtaF,"recoLep2EtaF/f");
  TreeEleMuo->Branch("recoLep2PhiF", &m_recoLep2PhiF,"recoLep2PhiF/f");
  TreeEleMuo->Branch("recoLep1Pt", &m_recoLep1Pt,"recoLep1Pt/i");
  TreeEleMuo->Branch("recoLep1Eta", &m_recoLep1Eta,"recoLep1Eta/i");
  TreeEleMuo->Branch("recoLep1Phi", &m_recoLep1Phi,"recoLep1Phi/i");
  TreeEleMuo->Branch("recoLep1HEEP", &m_recoLep1HEEP,"recoLep1HEEP/i");
  TreeEleMuo->Branch("recoLep1dRJet", &m_recoLep1dRJet,"recoLep1dRJet/i");
  TreeEleMuo->Branch("recoLep1PFIso", &m_recoLep1PFIso,"recoLep1PFIso/i");
  TreeEleMuo->Branch("recoLep1CorrPFIso", &m_recoLep1CorrPFIso,"recoLep1CorrPFIso/i");
  TreeEleMuo->Branch("recoLep1DETIso", &m_recoLep1DETIso,"recoLep1DETIso/i");
  TreeEleMuo->Branch("recoLep1HEEPIso", &m_recoLep1HEEPIso,"recoLep1HEEPIso/i");
  TreeEleMuo->Branch("recoLep1dEtaIn", &m_recoLep1dEtaIn,"recoLep1dEtaIn/i");
  TreeEleMuo->Branch("recoLep1dPhiIn", &m_recoLep1dPhiIn,"recoLep1dPhiIn/i");
  TreeEleMuo->Branch("recoLep1SigmaIEIE", &m_recoLep1SigmaIEIE,"recoLep1SigmaIEIE/i");
  TreeEleMuo->Branch("recoLep1HoverE", &m_recoLep1HoverE,"recoLep1HoverE/i");
  TreeEleMuo->Branch("recoLep1dxy", &m_recoLep1dxy,"recoLep1dxy/i");
  TreeEleMuo->Branch("recoLep1dz", &m_recoLep1dz,"recoLep1dz/i");
  TreeEleMuo->Branch("recoLep1EandP", &m_recoLep1EandP,"recoLep1EandP/i");
  TreeEleMuo->Branch("recoLep1Conv1", &m_recoLep1Conv1,"recoLep1Conv1/i");
  TreeEleMuo->Branch("recoLep1Conv2", &m_recoLep1Conv2,"recoLep1Conv2/i");
  TreeEleMuo->Branch("recoLep1EcalDriven", &m_recoLep1EcalDriven,"recoLep1EcalDriven/i");
  TreeEleMuo->Branch("recoLep1dPhiInHEEP", &m_recoLep1dPhiInHEEP,"recoLep1dPhiInHEEP/i");
  TreeEleMuo->Branch("recoLep1HoverEHEEP", &m_recoLep1HoverEHEEP,"recoLep1HoverEHEEP/i");
  TreeEleMuo->Branch("recoLep1TkSumPt", &m_recoLep1TkSumPt,"recoLep1TkSumPt/i");
  TreeEleMuo->Branch("recoLep1LostHits", &m_recoLep1LostHits,"recoLep1LostHits/i");
  TreeEleMuo->Branch("recoLep1dEtaInHEEP", &m_recoLep1dEtaInHEEP,"recoLep1dEtaInHEEP/i");
  TreeEleMuo->Branch("recoLep1dxyHEEP", &m_recoLep1dxyHEEP,"recoLep1dxyHEEP/i");
  TreeEleMuo->Branch("recoLep1e5xe5", &m_recoLep1e5xe5,"recoLep1e5xe5/i");
  TreeEleMuo->Branch("recoLep1SigmaIEIEHEEP", &m_recoLep1SigmaIEIEHEEP,"recoLep1SigmaIEIEHEEP/i");
  TreeEleMuo->Branch("recoLep1IsGlobal", &m_recoLep1IsGlobal,"recoLep1IsGlobal/i");
  TreeEleMuo->Branch("recoLep1IsPFTight", &m_recoLep1IsPFTight,"recoLep1IsPFTight/i");
  TreeEleMuo->Branch("recoLep1Chi2Tight", &m_recoLep1Chi2Tight,"recoLep1Chi2Tight/i");
  TreeEleMuo->Branch("recoLep1MuonHitsTight", &m_recoLep1MuonHitsTight,"recoLep1MuonHitsTight/i");
  TreeEleMuo->Branch("recoLep1MatchesTight", &m_recoLep1MatchesTight,"recoLep1MatchesTight/i");
  TreeEleMuo->Branch("recoLep1dxyTight", &m_recoLep1dxyTight,"recoLep1dxyTight/i");
  TreeEleMuo->Branch("recoLep1dzTight", &m_recoLep1dzTight,"recoLep1dzTight/i");
  TreeEleMuo->Branch("recoLep1PixelHitsTight", &m_recoLep1PixelHitsTight,"recoLep1PixelHitsTight/i");
  TreeEleMuo->Branch("recoLep1TrackerLTight", &m_recoLep1TrackerLTight,"recoLep1TrackerLTight/i");
  TreeEleMuo->Branch("recoLep1ptErr", &m_recoLep1ptErr,"recoLep1ptErr/i");
  TreeEleMuo->Branch("recoLep1MuonHits", &m_recoLep1MuonHits,"recoLep1MuonHits/i");
  TreeEleMuo->Branch("recoLep1Matches", &m_recoLep1Matches,"recoLep1Matches/i");
  TreeEleMuo->Branch("recoLep1PixelHits", &m_recoLep1PixelHits,"recoLep1PixelHits/i");
  TreeEleMuo->Branch("recoLep1TrackerL", &m_recoLep1TrackerL,"recoLep1TrackerL/i");
  TreeEleMuo->Branch("recoLep1IsPFMuon", &m_recoLep1IsPFMuon,"recoLep1IsPFMuon/i");
  TreeEleMuo->Branch("recoLep1Chi2", &m_recoLep1Chi2,"recoLep1Chi2/i");
  TreeEleMuo->Branch("recoLep1Discr1", &m_recoLep1Discr1,"recoLep1Discr1/i");
  TreeEleMuo->Branch("recoLep1Discr2", &m_recoLep1Discr2,"recoLep1Discr2/i");
  TreeEleMuo->Branch("recoLep1Discr3", &m_recoLep1Discr3,"recoLep1Discr3/i");
  TreeEleMuo->Branch("recoLep1Discr4", &m_recoLep1Discr4,"recoLep1Discr4/i");
  TreeEleMuo->Branch("recoLep2Pt", &m_recoLep2Pt,"recoLep2Pt/i");
  TreeEleMuo->Branch("recoLep2Eta", &m_recoLep2Eta,"recoLep2Eta/i");
  TreeEleMuo->Branch("recoLep2Phi", &m_recoLep2Phi,"recoLep2Phi/i");
  TreeEleMuo->Branch("recoLep2HEEP", &m_recoLep2HEEP,"recoLep2HEEP/i");
  TreeEleMuo->Branch("recoLep2dRJet", &m_recoLep2dRJet,"recoLep2dRJet/i");
  TreeEleMuo->Branch("recoLep2PFIso", &m_recoLep2PFIso,"recoLep2PFIso/i");
  TreeEleMuo->Branch("recoLep2CorrPFIso", &m_recoLep2CorrPFIso,"recoLep2CorrPFIso/i");
  TreeEleMuo->Branch("recoLep2DETIso", &m_recoLep2DETIso,"recoLep2DETIso/i");
  TreeEleMuo->Branch("recoLep2HEEPIso", &m_recoLep2HEEPIso,"recoLep2HEEPIso/i");
  TreeEleMuo->Branch("recoLep2dEtaIn", &m_recoLep2dEtaIn,"recoLep2dEtaIn/i");
  TreeEleMuo->Branch("recoLep2dPhiIn", &m_recoLep2dPhiIn,"recoLep2dPhiIn/i");
  TreeEleMuo->Branch("recoLep2SigmaIEIE", &m_recoLep2SigmaIEIE,"recoLep2SigmaIEIE/i");
  TreeEleMuo->Branch("recoLep2HoverE", &m_recoLep2HoverE,"recoLep2HoverE/i");
  TreeEleMuo->Branch("recoLep2dxy", &m_recoLep2dxy,"recoLep2dxy/i");
  TreeEleMuo->Branch("recoLep2dz", &m_recoLep2dz,"recoLep2dz/i");
  TreeEleMuo->Branch("recoLep2EandP", &m_recoLep2EandP,"recoLep2EandP/i");
  TreeEleMuo->Branch("recoLep2Conv1", &m_recoLep2Conv1,"recoLep2Conv1/i");
  TreeEleMuo->Branch("recoLep2Conv2", &m_recoLep2Conv2,"recoLep2Conv2/i");
  TreeEleMuo->Branch("recoLep2EcalDriven", &m_recoLep2EcalDriven,"recoLep2EcalDriven/i");
  TreeEleMuo->Branch("recoLep2dPhiInHEEP", &m_recoLep2dPhiInHEEP,"recoLep2dPhiInHEEP/i");
  TreeEleMuo->Branch("recoLep2HoverEHEEP", &m_recoLep2HoverEHEEP,"recoLep2HoverEHEEP/i");
  TreeEleMuo->Branch("recoLep2TkSumPt", &m_recoLep2TkSumPt,"recoLep2TkSumPt/i");
  TreeEleMuo->Branch("recoLep2LostHits", &m_recoLep2LostHits,"recoLep2LostHits/i");
  TreeEleMuo->Branch("recoLep2dEtaInHEEP", &m_recoLep2dEtaInHEEP,"recoLep2dEtaInHEEP/i");
  TreeEleMuo->Branch("recoLep2dxyHEEP", &m_recoLep2dxyHEEP,"recoLep2dxyHEEP/i");
  TreeEleMuo->Branch("recoLep2e5xe5", &m_recoLep2e5xe5,"recoLep2e5xe5/i");
  TreeEleMuo->Branch("recoLep2SigmaIEIEHEEP", &m_recoLep2SigmaIEIEHEEP,"recoLep2SigmaIEIEHEEP/i");
  TreeEleMuo->Branch("recoLep2IsGlobal", &m_recoLep2IsGlobal,"recoLep2IsGlobal/i");
  TreeEleMuo->Branch("recoLep2IsPFTight", &m_recoLep2IsPFTight,"recoLep2IsPFTight/i");
  TreeEleMuo->Branch("recoLep2Chi2Tight", &m_recoLep2Chi2Tight,"recoLep2Chi2Tight/i");
  TreeEleMuo->Branch("recoLep2MuonHitsTight", &m_recoLep2MuonHitsTight,"recoLep2MuonHitsTight/i");
  TreeEleMuo->Branch("recoLep2MatchesTight", &m_recoLep2MatchesTight,"recoLep2MatchesTight/i");
  TreeEleMuo->Branch("recoLep2dxyTight", &m_recoLep2dxyTight,"recoLep2dxyTight/i");
  TreeEleMuo->Branch("recoLep2dzTight", &m_recoLep2dzTight,"recoLep2dzTight/i");
  TreeEleMuo->Branch("recoLep2PixelHitsTight", &m_recoLep2PixelHitsTight,"recoLep2PixelHitsTight/i");
  TreeEleMuo->Branch("recoLep2TrackerLTight", &m_recoLep2TrackerLTight,"recoLep2TrackerLTight/i");
  TreeEleMuo->Branch("recoLep2ptErr", &m_recoLep2ptErr,"recoLep2ptErr/i");
  TreeEleMuo->Branch("recoLep2MuonHits", &m_recoLep2MuonHits,"recoLep2MuonHits/i");
  TreeEleMuo->Branch("recoLep2Matches", &m_recoLep2Matches,"recoLep2Matches/i");
  TreeEleMuo->Branch("recoLep2PixelHits", &m_recoLep2PixelHits,"recoLep2PixelHits/i");
  TreeEleMuo->Branch("recoLep2TrackerL", &m_recoLep2TrackerL,"recoLep2TrackerL/i");
  TreeEleMuo->Branch("recoLep2IsPFMuon", &m_recoLep2IsPFMuon,"recoLep2IsPFMuon/i");
  TreeEleMuo->Branch("recoLep2Chi2", &m_recoLep2Chi2,"recoLep2Chi2/i");

  TreeMuoMuo->Branch("genJetPt", &m_genJetPt,"genJetPt/f");
  TreeMuoMuo->Branch("genJetEta", &m_genJetEta,"genJetEta/f");
  TreeMuoMuo->Branch("genJetPhi", &m_genJetPhi,"genJetPhi/f");
  TreeMuoMuo->Branch("genLep1Pt", &m_genLep1Pt,"genLep1Pt/f");
  TreeMuoMuo->Branch("genLep2Pt", &m_genLep2Pt,"genLep2Pt/f");
  TreeMuoMuo->Branch("genLep1Eta", &m_genLep1Eta,"genLep1Eta/f");
  TreeMuoMuo->Branch("genLep2Eta", &m_genLep2Eta,"genLep2Eta/f");
  TreeMuoMuo->Branch("genLep1Phi", &m_genLep1Phi,"genLep1Phi/f");
  TreeMuoMuo->Branch("genLep2Phi", &m_genLep2Phi,"genLep2Phi/f");
  TreeMuoMuo->Branch("recoJetPt", &m_recoJetPt,"recoJetPt/f");
  TreeMuoMuo->Branch("recoJetEta", &m_recoJetEta,"recoJetEta/f");
  TreeMuoMuo->Branch("recoJetPhi", &m_recoJetPhi,"recoJetPhi/f");
  TreeMuoMuo->Branch("recoDeltaR", &m_recoDeltaR,"recoDeltaR/f");
  TreeMuoMuo->Branch("dRGenRecoLep1", &m_dRGenRecoLep1,"dRGenRecoLep1/f");
  TreeMuoMuo->Branch("dRGenRecoLep2", &m_dRGenRecoLep2,"dRGenRecoLep2/f");
  TreeMuoMuo->Branch("genDeltaR", &m_genDeltaR,"genDeltaR/f");
  TreeMuoMuo->Branch("dRGenRecoJet", &m_dRGenRecoJet,"dRGenRecoJet/f");
  TreeMuoMuo->Branch("recoLep1dRJetF", &m_recoLep1dRJetF,"recoLep1dRJetF/f");
  TreeMuoMuo->Branch("recoLep2dRJetF", &m_recoLep2dRJetF,"recoLep2dRJetF/f");
  TreeMuoMuo->Branch("recoLep1PtF", &m_recoLep1PtF,"recoLep1PtF/f");
  TreeMuoMuo->Branch("recoLep1EtaF", &m_recoLep1EtaF,"recoLep1EtaF/f");
  TreeMuoMuo->Branch("recoLep1PhiF", &m_recoLep1PhiF,"recoLep1PhiF/f");
  TreeMuoMuo->Branch("recoLep2PtF", &m_recoLep2PtF,"recoLep2PtF/f");
  TreeMuoMuo->Branch("recoLep2EtaF", &m_recoLep2EtaF,"recoLep2EtaF/f");
  TreeMuoMuo->Branch("recoLep2PhiF", &m_recoLep2PhiF,"recoLep2PhiF/f");
  TreeMuoMuo->Branch("recoLep1Pt", &m_recoLep1Pt,"recoLep1Pt/i");
  TreeMuoMuo->Branch("recoLep1Eta", &m_recoLep1Eta,"recoLep1Eta/i");
  TreeMuoMuo->Branch("recoLep1Phi", &m_recoLep1Phi,"recoLep1Phi/i");
  TreeMuoMuo->Branch("recoLep1HEEP", &m_recoLep1HEEP,"recoLep1HEEP/i");
  TreeMuoMuo->Branch("recoLep1dRJet", &m_recoLep1dRJet,"recoLep1dRJet/i");
  TreeMuoMuo->Branch("recoLep1PFIso", &m_recoLep1PFIso,"recoLep1PFIso/i");
  TreeMuoMuo->Branch("recoLep1CorrPFIso", &m_recoLep1CorrPFIso,"recoLep1CorrPFIso/i");
  TreeMuoMuo->Branch("recoLep1DETIso", &m_recoLep1DETIso,"recoLep1DETIso/i");
  TreeMuoMuo->Branch("recoLep1HEEPIso", &m_recoLep1HEEPIso,"recoLep1HEEPIso/i");
  TreeMuoMuo->Branch("recoLep1dEtaIn", &m_recoLep1dEtaIn,"recoLep1dEtaIn/i");
  TreeMuoMuo->Branch("recoLep1dPhiIn", &m_recoLep1dPhiIn,"recoLep1dPhiIn/i");
  TreeMuoMuo->Branch("recoLep1SigmaIEIE", &m_recoLep1SigmaIEIE,"recoLep1SigmaIEIE/i");
  TreeMuoMuo->Branch("recoLep1HoverE", &m_recoLep1HoverE,"recoLep1HoverE/i");
  TreeMuoMuo->Branch("recoLep1dxy", &m_recoLep1dxy,"recoLep1dxy/i");
  TreeMuoMuo->Branch("recoLep1dz", &m_recoLep1dz,"recoLep1dz/i");
  TreeMuoMuo->Branch("recoLep1EandP", &m_recoLep1EandP,"recoLep1EandP/i");
  TreeMuoMuo->Branch("recoLep1Conv1", &m_recoLep1Conv1,"recoLep1Conv1/i");
  TreeMuoMuo->Branch("recoLep1Conv2", &m_recoLep1Conv2,"recoLep1Conv2/i");
  TreeMuoMuo->Branch("recoLep1EcalDriven", &m_recoLep1EcalDriven,"recoLep1EcalDriven/i");
  TreeMuoMuo->Branch("recoLep1dPhiInHEEP", &m_recoLep1dPhiInHEEP,"recoLep1dPhiInHEEP/i");
  TreeMuoMuo->Branch("recoLep1HoverEHEEP", &m_recoLep1HoverEHEEP,"recoLep1HoverEHEEP/i");
  TreeMuoMuo->Branch("recoLep1TkSumPt", &m_recoLep1TkSumPt,"recoLep1TkSumPt/i");
  TreeMuoMuo->Branch("recoLep1LostHits", &m_recoLep1LostHits,"recoLep1LostHits/i");
  TreeMuoMuo->Branch("recoLep1dEtaInHEEP", &m_recoLep1dEtaInHEEP,"recoLep1dEtaInHEEP/i");
  TreeMuoMuo->Branch("recoLep1dxyHEEP", &m_recoLep1dxyHEEP,"recoLep1dxyHEEP/i");
  TreeMuoMuo->Branch("recoLep1e5xe5", &m_recoLep1e5xe5,"recoLep1e5xe5/i");
  TreeMuoMuo->Branch("recoLep1SigmaIEIEHEEP", &m_recoLep1SigmaIEIEHEEP,"recoLep1SigmaIEIEHEEP/i");
  TreeMuoMuo->Branch("recoLep1IsGlobal", &m_recoLep1IsGlobal,"recoLep1IsGlobal/i");
  TreeMuoMuo->Branch("recoLep1IsPFTight", &m_recoLep1IsPFTight,"recoLep1IsPFTight/i");
  TreeMuoMuo->Branch("recoLep1Chi2Tight", &m_recoLep1Chi2Tight,"recoLep1Chi2Tight/i");
  TreeMuoMuo->Branch("recoLep1MuonHitsTight", &m_recoLep1MuonHitsTight,"recoLep1MuonHitsTight/i");
  TreeMuoMuo->Branch("recoLep1MatchesTight", &m_recoLep1MatchesTight,"recoLep1MatchesTight/i");
  TreeMuoMuo->Branch("recoLep1dxyTight", &m_recoLep1dxyTight,"recoLep1dxyTight/i");
  TreeMuoMuo->Branch("recoLep1dzTight", &m_recoLep1dzTight,"recoLep1dzTight/i");
  TreeMuoMuo->Branch("recoLep1PixelHitsTight", &m_recoLep1PixelHitsTight,"recoLep1PixelHitsTight/i");
  TreeMuoMuo->Branch("recoLep1TrackerLTight", &m_recoLep1TrackerLTight,"recoLep1TrackerLTight/i");
  TreeMuoMuo->Branch("recoLep1ptErr", &m_recoLep1ptErr,"recoLep1ptErr/i");
  TreeMuoMuo->Branch("recoLep1MuonHits", &m_recoLep1MuonHits,"recoLep1MuonHits/i");
  TreeMuoMuo->Branch("recoLep1Matches", &m_recoLep1Matches,"recoLep1Matches/i");
  TreeMuoMuo->Branch("recoLep1PixelHits", &m_recoLep1PixelHits,"recoLep1PixelHits/i");
  TreeMuoMuo->Branch("recoLep1TrackerL", &m_recoLep1TrackerL,"recoLep1TrackerL/i");
  TreeMuoMuo->Branch("recoLep1IsPFMuon", &m_recoLep1IsPFMuon,"recoLep1IsPFMuon/i");
  TreeMuoMuo->Branch("recoLep1Chi2", &m_recoLep1Chi2,"recoLep1Chi2/i");
  TreeMuoMuo->Branch("recoLep1Discr1", &m_recoLep1Discr1,"recoLep1Discr1/i");
  TreeMuoMuo->Branch("recoLep1Discr2", &m_recoLep1Discr2,"recoLep1Discr2/i");
  TreeMuoMuo->Branch("recoLep1Discr3", &m_recoLep1Discr3,"recoLep1Discr3/i");
  TreeMuoMuo->Branch("recoLep1Discr4", &m_recoLep1Discr4,"recoLep1Discr4/i");
  TreeMuoMuo->Branch("recoLep2Pt", &m_recoLep2Pt,"recoLep2Pt/i");
  TreeMuoMuo->Branch("recoLep2Eta", &m_recoLep2Eta,"recoLep2Eta/i");
  TreeMuoMuo->Branch("recoLep2Phi", &m_recoLep2Phi,"recoLep2Phi/i");
  TreeMuoMuo->Branch("recoLep2HEEP", &m_recoLep2HEEP,"recoLep2HEEP/i");
  TreeMuoMuo->Branch("recoLep2dRJet", &m_recoLep2dRJet,"recoLep2dRJet/i");
  TreeMuoMuo->Branch("recoLep2PFIso", &m_recoLep2PFIso,"recoLep2PFIso/i");
  TreeMuoMuo->Branch("recoLep2CorrPFIso", &m_recoLep2CorrPFIso,"recoLep2CorrPFIso/i");
  TreeMuoMuo->Branch("recoLep2DETIso", &m_recoLep2DETIso,"recoLep2DETIso/i");
  TreeMuoMuo->Branch("recoLep2HEEPIso", &m_recoLep2HEEPIso,"recoLep2HEEPIso/i");
  TreeMuoMuo->Branch("recoLep2dEtaIn", &m_recoLep2dEtaIn,"recoLep2dEtaIn/i");
  TreeMuoMuo->Branch("recoLep2dPhiIn", &m_recoLep2dPhiIn,"recoLep2dPhiIn/i");
  TreeMuoMuo->Branch("recoLep2SigmaIEIE", &m_recoLep2SigmaIEIE,"recoLep2SigmaIEIE/i");
  TreeMuoMuo->Branch("recoLep2HoverE", &m_recoLep2HoverE,"recoLep2HoverE/i");
  TreeMuoMuo->Branch("recoLep2dxy", &m_recoLep2dxy,"recoLep2dxy/i");
  TreeMuoMuo->Branch("recoLep2dz", &m_recoLep2dz,"recoLep2dz/i");
  TreeMuoMuo->Branch("recoLep2EandP", &m_recoLep2EandP,"recoLep2EandP/i");
  TreeMuoMuo->Branch("recoLep2Conv1", &m_recoLep2Conv1,"recoLep2Conv1/i");
  TreeMuoMuo->Branch("recoLep2Conv2", &m_recoLep2Conv2,"recoLep2Conv2/i");
  TreeMuoMuo->Branch("recoLep2EcalDriven", &m_recoLep2EcalDriven,"recoLep2EcalDriven/i");
  TreeMuoMuo->Branch("recoLep2dPhiInHEEP", &m_recoLep2dPhiInHEEP,"recoLep2dPhiInHEEP/i");
  TreeMuoMuo->Branch("recoLep2HoverEHEEP", &m_recoLep2HoverEHEEP,"recoLep2HoverEHEEP/i");
  TreeMuoMuo->Branch("recoLep2TkSumPt", &m_recoLep2TkSumPt,"recoLep2TkSumPt/i");
  TreeMuoMuo->Branch("recoLep2LostHits", &m_recoLep2LostHits,"recoLep2LostHits/i");
  TreeMuoMuo->Branch("recoLep2dEtaInHEEP", &m_recoLep2dEtaInHEEP,"recoLep2dEtaInHEEP/i");
  TreeMuoMuo->Branch("recoLep2dxyHEEP", &m_recoLep2dxyHEEP,"recoLep2dxyHEEP/i");
  TreeMuoMuo->Branch("recoLep2e5xe5", &m_recoLep2e5xe5,"recoLep2e5xe5/i");
  TreeMuoMuo->Branch("recoLep2SigmaIEIEHEEP", &m_recoLep2SigmaIEIEHEEP,"recoLep2SigmaIEIEHEEP/i");
  TreeMuoMuo->Branch("recoLep2IsGlobal", &m_recoLep2IsGlobal,"recoLep2IsGlobal/i");
  TreeMuoMuo->Branch("recoLep2IsPFTight", &m_recoLep2IsPFTight,"recoLep2IsPFTight/i");
  TreeMuoMuo->Branch("recoLep2Chi2Tight", &m_recoLep2Chi2Tight,"recoLep2Chi2Tight/i");
  TreeMuoMuo->Branch("recoLep2MuonHitsTight", &m_recoLep2MuonHitsTight,"recoLep2MuonHitsTight/i");
  TreeMuoMuo->Branch("recoLep2MatchesTight", &m_recoLep2MatchesTight,"recoLep2MatchesTight/i");
  TreeMuoMuo->Branch("recoLep2dxyTight", &m_recoLep2dxyTight,"recoLep2dxyTight/i");
  TreeMuoMuo->Branch("recoLep2dzTight", &m_recoLep2dzTight,"recoLep2dzTight/i");
  TreeMuoMuo->Branch("recoLep2PixelHitsTight", &m_recoLep2PixelHitsTight,"recoLep2PixelHitsTight/i");
  TreeMuoMuo->Branch("recoLep2TrackerLTight", &m_recoLep2TrackerLTight,"recoLep2TrackerLTight/i");
  TreeMuoMuo->Branch("recoLep2ptErr", &m_recoLep2ptErr,"recoLep2ptErr/i");
  TreeMuoMuo->Branch("recoLep2MuonHits", &m_recoLep2MuonHits,"recoLep2MuonHits/i");
  TreeMuoMuo->Branch("recoLep2Matches", &m_recoLep2Matches,"recoLep2Matches/i");
  TreeMuoMuo->Branch("recoLep2PixelHits", &m_recoLep2PixelHits,"recoLep2PixelHits/i");
  TreeMuoMuo->Branch("recoLep2TrackerL", &m_recoLep2TrackerL,"recoLep2TrackerL/i");
  TreeMuoMuo->Branch("recoLep2IsPFMuon", &m_recoLep2IsPFMuon,"recoLep2IsPFMuon/i");
  TreeMuoMuo->Branch("recoLep2Chi2", &m_recoLep2Chi2,"recoLep2Chi2/i");

  TreeEleEle->Branch("genJetPt", &m_genJetPt,"genJetPt/f");
  TreeEleEle->Branch("genJetEta", &m_genJetEta,"genJetEta/f");
  TreeEleEle->Branch("genJetPhi", &m_genJetPhi,"genJetPhi/f");
  TreeEleEle->Branch("genLep1Pt", &m_genLep1Pt,"genLep1Pt/f");
  TreeEleEle->Branch("genLep2Pt", &m_genLep2Pt,"genLep2Pt/f");
  TreeEleEle->Branch("genLep1Eta", &m_genLep1Eta,"genLep1Eta/f");
  TreeEleEle->Branch("genLep2Eta", &m_genLep2Eta,"genLep2Eta/f");
  TreeEleEle->Branch("genLep1Phi", &m_genLep1Phi,"genLep1Phi/f");
  TreeEleEle->Branch("genLep2Phi", &m_genLep2Phi,"genLep2Phi/f");
  TreeEleEle->Branch("recoJetPt", &m_recoJetPt,"recoJetPt/f");
  TreeEleEle->Branch("recoJetEta", &m_recoJetEta,"recoJetEta/f");
  TreeEleEle->Branch("recoJetPhi", &m_recoJetPhi,"recoJetPhi/f");
  TreeEleEle->Branch("recoDeltaR", &m_recoDeltaR,"recoDeltaR/f");
  TreeEleEle->Branch("dRGenRecoLep1", &m_dRGenRecoLep1,"dRGenRecoLep1/f");
  TreeEleEle->Branch("dRGenRecoLep2", &m_dRGenRecoLep2,"dRGenRecoLep2/f");
  TreeEleEle->Branch("genDeltaR", &m_genDeltaR,"genDeltaR/f");
  TreeEleEle->Branch("dRGenRecoJet", &m_dRGenRecoJet,"dRGenRecoJet/f");
  TreeEleEle->Branch("recoLep1dRJetF", &m_recoLep1dRJetF,"recoLep1dRJetF/f");
  TreeEleEle->Branch("recoLep2dRJetF", &m_recoLep2dRJetF,"recoLep2dRJetF/f");
  TreeEleEle->Branch("recoLep1PtF", &m_recoLep1PtF,"recoLep1PtF/f");
  TreeEleEle->Branch("recoLep1EtaF", &m_recoLep1EtaF,"recoLep1EtaF/f");
  TreeEleEle->Branch("recoLep1PhiF", &m_recoLep1PhiF,"recoLep1PhiF/f");
  TreeEleEle->Branch("recoLep2PtF", &m_recoLep2PtF,"recoLep2PtF/f");
  TreeEleEle->Branch("recoLep2EtaF", &m_recoLep2EtaF,"recoLep2EtaF/f");
  TreeEleEle->Branch("recoLep2PhiF", &m_recoLep2PhiF,"recoLep2PhiF/f");
  TreeEleEle->Branch("recoLep1Pt", &m_recoLep1Pt,"recoLep1Pt/i");
  TreeEleEle->Branch("recoLep1Eta", &m_recoLep1Eta,"recoLep1Eta/i");
  TreeEleEle->Branch("recoLep1Phi", &m_recoLep1Phi,"recoLep1Phi/i");
  TreeEleEle->Branch("recoLep1HEEP", &m_recoLep1HEEP,"recoLep1HEEP/i");
  TreeEleEle->Branch("recoLep1dRJet", &m_recoLep1dRJet,"recoLep1dRJet/i");
  TreeEleEle->Branch("recoLep1PFIso", &m_recoLep1PFIso,"recoLep1PFIso/i");
  TreeEleEle->Branch("recoLep1CorrPFIso", &m_recoLep1CorrPFIso,"recoLep1CorrPFIso/i");
  TreeEleEle->Branch("recoLep1DETIso", &m_recoLep1DETIso,"recoLep1DETIso/i");
  TreeEleEle->Branch("recoLep1HEEPIso", &m_recoLep1HEEPIso,"recoLep1HEEPIso/i");
  TreeEleEle->Branch("recoLep1dEtaIn", &m_recoLep1dEtaIn,"recoLep1dEtaIn/i");
  TreeEleEle->Branch("recoLep1dPhiIn", &m_recoLep1dPhiIn,"recoLep1dPhiIn/i");
  TreeEleEle->Branch("recoLep1SigmaIEIE", &m_recoLep1SigmaIEIE,"recoLep1SigmaIEIE/i");
  TreeEleEle->Branch("recoLep1HoverE", &m_recoLep1HoverE,"recoLep1HoverE/i");
  TreeEleEle->Branch("recoLep1dxy", &m_recoLep1dxy,"recoLep1dxy/i");
  TreeEleEle->Branch("recoLep1dz", &m_recoLep1dz,"recoLep1dz/i");
  TreeEleEle->Branch("recoLep1EandP", &m_recoLep1EandP,"recoLep1EandP/i");
  TreeEleEle->Branch("recoLep1Conv1", &m_recoLep1Conv1,"recoLep1Conv1/i");
  TreeEleEle->Branch("recoLep1Conv2", &m_recoLep1Conv2,"recoLep1Conv2/i");
  TreeEleEle->Branch("recoLep1EcalDriven", &m_recoLep1EcalDriven,"recoLep1EcalDriven/i");
  TreeEleEle->Branch("recoLep1dPhiInHEEP", &m_recoLep1dPhiInHEEP,"recoLep1dPhiInHEEP/i");
  TreeEleEle->Branch("recoLep1HoverEHEEP", &m_recoLep1HoverEHEEP,"recoLep1HoverEHEEP/i");
  TreeEleEle->Branch("recoLep1TkSumPt", &m_recoLep1TkSumPt,"recoLep1TkSumPt/i");
  TreeEleEle->Branch("recoLep1LostHits", &m_recoLep1LostHits,"recoLep1LostHits/i");
  TreeEleEle->Branch("recoLep1dEtaInHEEP", &m_recoLep1dEtaInHEEP,"recoLep1dEtaInHEEP/i");
  TreeEleEle->Branch("recoLep1dxyHEEP", &m_recoLep1dxyHEEP,"recoLep1dxyHEEP/i");
  TreeEleEle->Branch("recoLep1e5xe5", &m_recoLep1e5xe5,"recoLep1e5xe5/i");
  TreeEleEle->Branch("recoLep1SigmaIEIEHEEP", &m_recoLep1SigmaIEIEHEEP,"recoLep1SigmaIEIEHEEP/i");
  TreeEleEle->Branch("recoLep1IsGlobal", &m_recoLep1IsGlobal,"recoLep1IsGlobal/i");
  TreeEleEle->Branch("recoLep1IsPFTight", &m_recoLep1IsPFTight,"recoLep1IsPFTight/i");
  TreeEleEle->Branch("recoLep1Chi2Tight", &m_recoLep1Chi2Tight,"recoLep1Chi2Tight/i");
  TreeEleEle->Branch("recoLep1MuonHitsTight", &m_recoLep1MuonHitsTight,"recoLep1MuonHitsTight/i");
  TreeEleEle->Branch("recoLep1MatchesTight", &m_recoLep1MatchesTight,"recoLep1MatchesTight/i");
  TreeEleEle->Branch("recoLep1dxyTight", &m_recoLep1dxyTight,"recoLep1dxyTight/i");
  TreeEleEle->Branch("recoLep1dzTight", &m_recoLep1dzTight,"recoLep1dzTight/i");
  TreeEleEle->Branch("recoLep1PixelHitsTight", &m_recoLep1PixelHitsTight,"recoLep1PixelHitsTight/i");
  TreeEleEle->Branch("recoLep1TrackerLTight", &m_recoLep1TrackerLTight,"recoLep1TrackerLTight/i");
  TreeEleEle->Branch("recoLep1ptErr", &m_recoLep1ptErr,"recoLep1ptErr/i");
  TreeEleEle->Branch("recoLep1MuonHits", &m_recoLep1MuonHits,"recoLep1MuonHits/i");
  TreeEleEle->Branch("recoLep1Matches", &m_recoLep1Matches,"recoLep1Matches/i");
  TreeEleEle->Branch("recoLep1PixelHits", &m_recoLep1PixelHits,"recoLep1PixelHits/i");
  TreeEleEle->Branch("recoLep1TrackerL", &m_recoLep1TrackerL,"recoLep1TrackerL/i");
  TreeEleEle->Branch("recoLep1IsPFMuon", &m_recoLep1IsPFMuon,"recoLep1IsPFMuon/i");
  TreeEleEle->Branch("recoLep1Chi2", &m_recoLep1Chi2,"recoLep1Chi2/i");
  TreeEleEle->Branch("recoLep1Discr1", &m_recoLep1Discr1,"recoLep1Discr1/i");
  TreeEleEle->Branch("recoLep1Discr2", &m_recoLep1Discr2,"recoLep1Discr2/i");
  TreeEleEle->Branch("recoLep1Discr3", &m_recoLep1Discr3,"recoLep1Discr3/i");
  TreeEleEle->Branch("recoLep1Discr4", &m_recoLep1Discr4,"recoLep1Discr4/i");
  TreeEleEle->Branch("recoLep2Pt", &m_recoLep2Pt,"recoLep2Pt/i");
  TreeEleEle->Branch("recoLep2Eta", &m_recoLep2Eta,"recoLep2Eta/i");
  TreeEleEle->Branch("recoLep2Phi", &m_recoLep2Phi,"recoLep2Phi/i");
  TreeEleEle->Branch("recoLep2HEEP", &m_recoLep2HEEP,"recoLep2HEEP/i");
  TreeEleEle->Branch("recoLep2dRJet", &m_recoLep2dRJet,"recoLep2dRJet/i");
  TreeEleEle->Branch("recoLep2PFIso", &m_recoLep2PFIso,"recoLep2PFIso/i");
  TreeEleEle->Branch("recoLep2CorrPFIso", &m_recoLep2CorrPFIso,"recoLep2CorrPFIso/i");
  TreeEleEle->Branch("recoLep2DETIso", &m_recoLep2DETIso,"recoLep2DETIso/i");
  TreeEleEle->Branch("recoLep2HEEPIso", &m_recoLep2HEEPIso,"recoLep2HEEPIso/i");
  TreeEleEle->Branch("recoLep2dEtaIn", &m_recoLep2dEtaIn,"recoLep2dEtaIn/i");
  TreeEleEle->Branch("recoLep2dPhiIn", &m_recoLep2dPhiIn,"recoLep2dPhiIn/i");
  TreeEleEle->Branch("recoLep2SigmaIEIE", &m_recoLep2SigmaIEIE,"recoLep2SigmaIEIE/i");
  TreeEleEle->Branch("recoLep2HoverE", &m_recoLep2HoverE,"recoLep2HoverE/i");
  TreeEleEle->Branch("recoLep2dxy", &m_recoLep2dxy,"recoLep2dxy/i");
  TreeEleEle->Branch("recoLep2dz", &m_recoLep2dz,"recoLep2dz/i");
  TreeEleEle->Branch("recoLep2EandP", &m_recoLep2EandP,"recoLep2EandP/i");
  TreeEleEle->Branch("recoLep2Conv1", &m_recoLep2Conv1,"recoLep2Conv1/i");
  TreeEleEle->Branch("recoLep2Conv2", &m_recoLep2Conv2,"recoLep2Conv2/i");
  TreeEleEle->Branch("recoLep2EcalDriven", &m_recoLep2EcalDriven,"recoLep2EcalDriven/i");
  TreeEleEle->Branch("recoLep2dPhiInHEEP", &m_recoLep2dPhiInHEEP,"recoLep2dPhiInHEEP/i");
  TreeEleEle->Branch("recoLep2HoverEHEEP", &m_recoLep2HoverEHEEP,"recoLep2HoverEHEEP/i");
  TreeEleEle->Branch("recoLep2TkSumPt", &m_recoLep2TkSumPt,"recoLep2TkSumPt/i");
  TreeEleEle->Branch("recoLep2LostHits", &m_recoLep2LostHits,"recoLep2LostHits/i");
  TreeEleEle->Branch("recoLep2dEtaInHEEP", &m_recoLep2dEtaInHEEP,"recoLep2dEtaInHEEP/i");
  TreeEleEle->Branch("recoLep2dxyHEEP", &m_recoLep2dxyHEEP,"recoLep2dxyHEEP/i");
  TreeEleEle->Branch("recoLep2e5xe5", &m_recoLep2e5xe5,"recoLep2e5xe5/i");
  TreeEleEle->Branch("recoLep2SigmaIEIEHEEP", &m_recoLep2SigmaIEIEHEEP,"recoLep2SigmaIEIEHEEP/i");
  TreeEleEle->Branch("recoLep2IsGlobal", &m_recoLep2IsGlobal,"recoLep2IsGlobal/i");
  TreeEleEle->Branch("recoLep2IsPFTight", &m_recoLep2IsPFTight,"recoLep2IsPFTight/i");
  TreeEleEle->Branch("recoLep2Chi2Tight", &m_recoLep2Chi2Tight,"recoLep2Chi2Tight/i");
  TreeEleEle->Branch("recoLep2MuonHitsTight", &m_recoLep2MuonHitsTight,"recoLep2MuonHitsTight/i");
  TreeEleEle->Branch("recoLep2MatchesTight", &m_recoLep2MatchesTight,"recoLep2MatchesTight/i");
  TreeEleEle->Branch("recoLep2dxyTight", &m_recoLep2dxyTight,"recoLep2dxyTight/i");
  TreeEleEle->Branch("recoLep2dzTight", &m_recoLep2dzTight,"recoLep2dzTight/i");
  TreeEleEle->Branch("recoLep2PixelHitsTight", &m_recoLep2PixelHitsTight,"recoLep2PixelHitsTight/i");
  TreeEleEle->Branch("recoLep2TrackerLTight", &m_recoLep2TrackerLTight,"recoLep2TrackerLTight/i");
  TreeEleEle->Branch("recoLep2ptErr", &m_recoLep2ptErr,"recoLep2ptErr/i");
  TreeEleEle->Branch("recoLep2MuonHits", &m_recoLep2MuonHits,"recoLep2MuonHits/i");
  TreeEleEle->Branch("recoLep2Matches", &m_recoLep2Matches,"recoLep2Matches/i");
  TreeEleEle->Branch("recoLep2PixelHits", &m_recoLep2PixelHits,"recoLep2PixelHits/i");
  TreeEleEle->Branch("recoLep2TrackerL", &m_recoLep2TrackerL,"recoLep2TrackerL/i");
  TreeEleEle->Branch("recoLep2IsPFMuon", &m_recoLep2IsPFMuon,"recoLep2IsPFMuon/i");
  TreeEleEle->Branch("recoLep2Chi2", &m_recoLep2Chi2,"recoLep2Chi2/i");

  TreeMuoTau->Branch("genJetPt", &m_genJetPt,"genJetPt/f");
  TreeMuoTau->Branch("genJetEta", &m_genJetEta,"genJetEta/f");
  TreeMuoTau->Branch("genJetPhi", &m_genJetPhi,"genJetPhi/f");
  TreeMuoTau->Branch("genLep1Pt", &m_genLep1Pt,"genLep1Pt/f");
  TreeMuoTau->Branch("genLep2Pt", &m_genLep2Pt,"genLep2Pt/f");
  TreeMuoTau->Branch("genLep1Eta", &m_genLep1Eta,"genLep1Eta/f");
  TreeMuoTau->Branch("genLep2Eta", &m_genLep2Eta,"genLep2Eta/f");
  TreeMuoTau->Branch("genLep1Phi", &m_genLep1Phi,"genLep1Phi/f");
  TreeMuoTau->Branch("genLep2Phi", &m_genLep2Phi,"genLep2Phi/f");
  TreeMuoTau->Branch("recoJetPt", &m_recoJetPt,"recoJetPt/f");
  TreeMuoTau->Branch("recoJetEta", &m_recoJetEta,"recoJetEta/f");
  TreeMuoTau->Branch("recoJetPhi", &m_recoJetPhi,"recoJetPhi/f");
  TreeMuoTau->Branch("recoDeltaR", &m_recoDeltaR,"recoDeltaR/f");
  TreeMuoTau->Branch("dRGenRecoLep1", &m_dRGenRecoLep1,"dRGenRecoLep1/f");
  TreeMuoTau->Branch("dRGenRecoLep2", &m_dRGenRecoLep2,"dRGenRecoLep2/f");
  TreeMuoTau->Branch("genDeltaR", &m_genDeltaR,"genDeltaR/f");
  TreeMuoTau->Branch("dRGenRecoJet", &m_dRGenRecoJet,"dRGenRecoJet/f");
  TreeMuoTau->Branch("recoLep1dRJetF", &m_recoLep1dRJetF,"recoLep1dRJetF/f");
  TreeMuoTau->Branch("recoLep2dRJetF", &m_recoLep2dRJetF,"recoLep2dRJetF/f");
  TreeMuoTau->Branch("recoLep1PtF", &m_recoLep1PtF,"recoLep1PtF/f");
  TreeMuoTau->Branch("recoLep1EtaF", &m_recoLep1EtaF,"recoLep1EtaF/f");
  TreeMuoTau->Branch("recoLep1PhiF", &m_recoLep1PhiF,"recoLep1PhiF/f");
  TreeMuoTau->Branch("recoLep2PtF", &m_recoLep2PtF,"recoLep2PtF/f");
  TreeMuoTau->Branch("recoLep2EtaF", &m_recoLep2EtaF,"recoLep2EtaF/f");
  TreeMuoTau->Branch("recoLep2PhiF", &m_recoLep2PhiF,"recoLep2PhiF/f");
  TreeMuoTau->Branch("recoLep1Pt", &m_recoLep1Pt,"recoLep1Pt/i");
  TreeMuoTau->Branch("recoLep1Eta", &m_recoLep1Eta,"recoLep1Eta/i");
  TreeMuoTau->Branch("recoLep1Phi", &m_recoLep1Phi,"recoLep1Phi/i");
  TreeMuoTau->Branch("recoLep1HEEP", &m_recoLep1HEEP,"recoLep1HEEP/i");
  TreeMuoTau->Branch("recoLep1dRJet", &m_recoLep1dRJet,"recoLep1dRJet/i");
  TreeMuoTau->Branch("recoLep1PFIso", &m_recoLep1PFIso,"recoLep1PFIso/i");
  TreeMuoTau->Branch("recoLep1CorrPFIso", &m_recoLep1CorrPFIso,"recoLep1CorrPFIso/i");
  TreeMuoTau->Branch("recoLep1DETIso", &m_recoLep1DETIso,"recoLep1DETIso/i");
  TreeMuoTau->Branch("recoLep1HEEPIso", &m_recoLep1HEEPIso,"recoLep1HEEPIso/i");
  TreeMuoTau->Branch("recoLep1dEtaIn", &m_recoLep1dEtaIn,"recoLep1dEtaIn/i");
  TreeMuoTau->Branch("recoLep1dPhiIn", &m_recoLep1dPhiIn,"recoLep1dPhiIn/i");
  TreeMuoTau->Branch("recoLep1SigmaIEIE", &m_recoLep1SigmaIEIE,"recoLep1SigmaIEIE/i");
  TreeMuoTau->Branch("recoLep1HoverE", &m_recoLep1HoverE,"recoLep1HoverE/i");
  TreeMuoTau->Branch("recoLep1dxy", &m_recoLep1dxy,"recoLep1dxy/i");
  TreeMuoTau->Branch("recoLep1dz", &m_recoLep1dz,"recoLep1dz/i");
  TreeMuoTau->Branch("recoLep1EandP", &m_recoLep1EandP,"recoLep1EandP/i");
  TreeMuoTau->Branch("recoLep1Conv1", &m_recoLep1Conv1,"recoLep1Conv1/i");
  TreeMuoTau->Branch("recoLep1Conv2", &m_recoLep1Conv2,"recoLep1Conv2/i");
  TreeMuoTau->Branch("recoLep1EcalDriven", &m_recoLep1EcalDriven,"recoLep1EcalDriven/i");
  TreeMuoTau->Branch("recoLep1dPhiInHEEP", &m_recoLep1dPhiInHEEP,"recoLep1dPhiInHEEP/i");
  TreeMuoTau->Branch("recoLep1HoverEHEEP", &m_recoLep1HoverEHEEP,"recoLep1HoverEHEEP/i");
  TreeMuoTau->Branch("recoLep1TkSumPt", &m_recoLep1TkSumPt,"recoLep1TkSumPt/i");
  TreeMuoTau->Branch("recoLep1LostHits", &m_recoLep1LostHits,"recoLep1LostHits/i");
  TreeMuoTau->Branch("recoLep1dEtaInHEEP", &m_recoLep1dEtaInHEEP,"recoLep1dEtaInHEEP/i");
  TreeMuoTau->Branch("recoLep1dxyHEEP", &m_recoLep1dxyHEEP,"recoLep1dxyHEEP/i");
  TreeMuoTau->Branch("recoLep1e5xe5", &m_recoLep1e5xe5,"recoLep1e5xe5/i");
  TreeMuoTau->Branch("recoLep1SigmaIEIEHEEP", &m_recoLep1SigmaIEIEHEEP,"recoLep1SigmaIEIEHEEP/i");
  TreeMuoTau->Branch("recoLep1IsGlobal", &m_recoLep1IsGlobal,"recoLep1IsGlobal/i");
  TreeMuoTau->Branch("recoLep1IsPFTight", &m_recoLep1IsPFTight,"recoLep1IsPFTight/i");
  TreeMuoTau->Branch("recoLep1Chi2Tight", &m_recoLep1Chi2Tight,"recoLep1Chi2Tight/i");
  TreeMuoTau->Branch("recoLep1MuonHitsTight", &m_recoLep1MuonHitsTight,"recoLep1MuonHitsTight/i");
  TreeMuoTau->Branch("recoLep1MatchesTight", &m_recoLep1MatchesTight,"recoLep1MatchesTight/i");
  TreeMuoTau->Branch("recoLep1dxyTight", &m_recoLep1dxyTight,"recoLep1dxyTight/i");
  TreeMuoTau->Branch("recoLep1dzTight", &m_recoLep1dzTight,"recoLep1dzTight/i");
  TreeMuoTau->Branch("recoLep1PixelHitsTight", &m_recoLep1PixelHitsTight,"recoLep1PixelHitsTight/i");
  TreeMuoTau->Branch("recoLep1TrackerLTight", &m_recoLep1TrackerLTight,"recoLep1TrackerLTight/i");
  TreeMuoTau->Branch("recoLep1ptErr", &m_recoLep1ptErr,"recoLep1ptErr/i");
  TreeMuoTau->Branch("recoLep1MuonHits", &m_recoLep1MuonHits,"recoLep1MuonHits/i");
  TreeMuoTau->Branch("recoLep1Matches", &m_recoLep1Matches,"recoLep1Matches/i");
  TreeMuoTau->Branch("recoLep1PixelHits", &m_recoLep1PixelHits,"recoLep1PixelHits/i");
  TreeMuoTau->Branch("recoLep1TrackerL", &m_recoLep1TrackerL,"recoLep1TrackerL/i");
  TreeMuoTau->Branch("recoLep1IsPFMuon", &m_recoLep1IsPFMuon,"recoLep1IsPFMuon/i");
  TreeMuoTau->Branch("recoLep1Chi2", &m_recoLep1Chi2,"recoLep1Chi2/i");
  TreeMuoTau->Branch("recoLep1Discr1", &m_recoLep1Discr1,"recoLep1Discr1/i");
  TreeMuoTau->Branch("recoLep1Discr2", &m_recoLep1Discr2,"recoLep1Discr2/i");
  TreeMuoTau->Branch("recoLep1Discr3", &m_recoLep1Discr3,"recoLep1Discr3/i");
  TreeMuoTau->Branch("recoLep1Discr4", &m_recoLep1Discr4,"recoLep1Discr4/i");
  TreeMuoTau->Branch("recoLep2Pt", &m_recoLep2Pt,"recoLep2Pt/i");
  TreeMuoTau->Branch("recoLep2Eta", &m_recoLep2Eta,"recoLep2Eta/i");
  TreeMuoTau->Branch("recoLep2Phi", &m_recoLep2Phi,"recoLep2Phi/i");
  TreeMuoTau->Branch("recoLep2HEEP", &m_recoLep2HEEP,"recoLep2HEEP/i");
  TreeMuoTau->Branch("recoLep2dRJet", &m_recoLep2dRJet,"recoLep2dRJet/i");
  TreeMuoTau->Branch("recoLep2PFIso", &m_recoLep2PFIso,"recoLep2PFIso/i");
  TreeMuoTau->Branch("recoLep2CorrPFIso", &m_recoLep2CorrPFIso,"recoLep2CorrPFIso/i");
  TreeMuoTau->Branch("recoLep2DETIso", &m_recoLep2DETIso,"recoLep2DETIso/i");
  TreeMuoTau->Branch("recoLep2HEEPIso", &m_recoLep2HEEPIso,"recoLep2HEEPIso/i");
  TreeMuoTau->Branch("recoLep2dEtaIn", &m_recoLep2dEtaIn,"recoLep2dEtaIn/i");
  TreeMuoTau->Branch("recoLep2dPhiIn", &m_recoLep2dPhiIn,"recoLep2dPhiIn/i");
  TreeMuoTau->Branch("recoLep2SigmaIEIE", &m_recoLep2SigmaIEIE,"recoLep2SigmaIEIE/i");
  TreeMuoTau->Branch("recoLep2HoverE", &m_recoLep2HoverE,"recoLep2HoverE/i");
  TreeMuoTau->Branch("recoLep2dxy", &m_recoLep2dxy,"recoLep2dxy/i");
  TreeMuoTau->Branch("recoLep2dz", &m_recoLep2dz,"recoLep2dz/i");
  TreeMuoTau->Branch("recoLep2EandP", &m_recoLep2EandP,"recoLep2EandP/i");
  TreeMuoTau->Branch("recoLep2Conv1", &m_recoLep2Conv1,"recoLep2Conv1/i");
  TreeMuoTau->Branch("recoLep2Conv2", &m_recoLep2Conv2,"recoLep2Conv2/i");
  TreeMuoTau->Branch("recoLep2EcalDriven", &m_recoLep2EcalDriven,"recoLep2EcalDriven/i");
  TreeMuoTau->Branch("recoLep2dPhiInHEEP", &m_recoLep2dPhiInHEEP,"recoLep2dPhiInHEEP/i");
  TreeMuoTau->Branch("recoLep2HoverEHEEP", &m_recoLep2HoverEHEEP,"recoLep2HoverEHEEP/i");
  TreeMuoTau->Branch("recoLep2TkSumPt", &m_recoLep2TkSumPt,"recoLep2TkSumPt/i");
  TreeMuoTau->Branch("recoLep2LostHits", &m_recoLep2LostHits,"recoLep2LostHits/i");
  TreeMuoTau->Branch("recoLep2dEtaInHEEP", &m_recoLep2dEtaInHEEP,"recoLep2dEtaInHEEP/i");
  TreeMuoTau->Branch("recoLep2dxyHEEP", &m_recoLep2dxyHEEP,"recoLep2dxyHEEP/i");
  TreeMuoTau->Branch("recoLep2e5xe5", &m_recoLep2e5xe5,"recoLep2e5xe5/i");
  TreeMuoTau->Branch("recoLep2SigmaIEIEHEEP", &m_recoLep2SigmaIEIEHEEP,"recoLep2SigmaIEIEHEEP/i");
  TreeMuoTau->Branch("recoLep2IsGlobal", &m_recoLep2IsGlobal,"recoLep2IsGlobal/i");
  TreeMuoTau->Branch("recoLep2IsPFTight", &m_recoLep2IsPFTight,"recoLep2IsPFTight/i");
  TreeMuoTau->Branch("recoLep2Chi2Tight", &m_recoLep2Chi2Tight,"recoLep2Chi2Tight/i");
  TreeMuoTau->Branch("recoLep2MuonHitsTight", &m_recoLep2MuonHitsTight,"recoLep2MuonHitsTight/i");
  TreeMuoTau->Branch("recoLep2MatchesTight", &m_recoLep2MatchesTight,"recoLep2MatchesTight/i");
  TreeMuoTau->Branch("recoLep2dxyTight", &m_recoLep2dxyTight,"recoLep2dxyTight/i");
  TreeMuoTau->Branch("recoLep2dzTight", &m_recoLep2dzTight,"recoLep2dzTight/i");
  TreeMuoTau->Branch("recoLep2PixelHitsTight", &m_recoLep2PixelHitsTight,"recoLep2PixelHitsTight/i");
  TreeMuoTau->Branch("recoLep2TrackerLTight", &m_recoLep2TrackerLTight,"recoLep2TrackerLTight/i");
  TreeMuoTau->Branch("recoLep2ptErr", &m_recoLep2ptErr,"recoLep2ptErr/i");
  TreeMuoTau->Branch("recoLep2MuonHits", &m_recoLep2MuonHits,"recoLep2MuonHits/i");
  TreeMuoTau->Branch("recoLep2Matches", &m_recoLep2Matches,"recoLep2Matches/i");
  TreeMuoTau->Branch("recoLep2PixelHits", &m_recoLep2PixelHits,"recoLep2PixelHits/i");
  TreeMuoTau->Branch("recoLep2TrackerL", &m_recoLep2TrackerL,"recoLep2TrackerL/i");
  TreeMuoTau->Branch("recoLep2IsPFMuon", &m_recoLep2IsPFMuon,"recoLep2IsPFMuon/i");
  TreeMuoTau->Branch("recoLep2Chi2", &m_recoLep2Chi2,"recoLep2Chi2/i");

  TreeEleTau->Branch("genJetPt", &m_genJetPt,"genJetPt/f");
  TreeEleTau->Branch("genJetEta", &m_genJetEta,"genJetEta/f");
  TreeEleTau->Branch("genJetPhi", &m_genJetPhi,"genJetPhi/f");
  TreeEleTau->Branch("genLep2Pt", &m_genLep1Pt,"genLep1Pt/f");
  TreeEleTau->Branch("genLep2Pt", &m_genLep2Pt,"genLep2Pt/f");
  TreeEleTau->Branch("genLep1Eta", &m_genLep1Eta,"genLep1Eta/f");
  TreeEleTau->Branch("genLep2Eta", &m_genLep2Eta,"genLep2Eta/f");
  TreeEleTau->Branch("genLep1Phi", &m_genLep1Phi,"genLep1Phi/f");
  TreeEleTau->Branch("genLep2Phi", &m_genLep2Phi,"genLep2Phi/f");
  TreeEleTau->Branch("recoJetPt", &m_recoJetPt,"recoJetPt/f");
  TreeEleTau->Branch("recoJetEta", &m_recoJetEta,"recoJetEta/f");
  TreeEleTau->Branch("recoJetPhi", &m_recoJetPhi,"recoJetPhi/f");
  TreeEleTau->Branch("recoDeltaR", &m_recoDeltaR,"recoDeltaR/f");
  TreeEleTau->Branch("dRGenRecoLep1", &m_dRGenRecoLep1,"dRGenRecoLep1/f");
  TreeEleTau->Branch("dRGenRecoLep2", &m_dRGenRecoLep2,"dRGenRecoLep2/f");
  TreeEleTau->Branch("genDeltaR", &m_genDeltaR,"genDeltaR/f");
  TreeEleTau->Branch("dRGenRecoJet", &m_dRGenRecoJet,"dRGenRecoJet/f");
  TreeEleTau->Branch("recoLep1dRJetF", &m_recoLep1dRJetF,"recoLep1dRJetF/f");
  TreeEleTau->Branch("recoLep2dRJetF", &m_recoLep2dRJetF,"recoLep2dRJetF/f");
  TreeEleTau->Branch("recoLep1PtF", &m_recoLep1PtF,"recoLep1PtF/f");
  TreeEleTau->Branch("recoLep1EtaF", &m_recoLep1EtaF,"recoLep1EtaF/f");
  TreeEleTau->Branch("recoLep1PhiF", &m_recoLep1PhiF,"recoLep1PhiF/f");
  TreeEleTau->Branch("recoLep2PtF", &m_recoLep2PtF,"recoLep2PtF/f");
  TreeEleTau->Branch("recoLep2EtaF", &m_recoLep2EtaF,"recoLep2EtaF/f");
  TreeEleTau->Branch("recoLep2PhiF", &m_recoLep2PhiF,"recoLep2PhiF/f");
  TreeEleTau->Branch("recoLep1Pt", &m_recoLep1Pt,"recoLep1Pt/i");
  TreeEleTau->Branch("recoLep1Eta", &m_recoLep1Eta,"recoLep1Eta/i");
  TreeEleTau->Branch("recoLep1Phi", &m_recoLep1Phi,"recoLep1Phi/i");
  TreeEleTau->Branch("recoLep1HEEP", &m_recoLep1HEEP,"recoLep1HEEP/i");
  TreeEleTau->Branch("recoLep1dRJet", &m_recoLep1dRJet,"recoLep1dRJet/i");
  TreeEleTau->Branch("recoLep1PFIso", &m_recoLep1PFIso,"recoLep1PFIso/i");
  TreeEleTau->Branch("recoLep1CorrPFIso", &m_recoLep1CorrPFIso,"recoLep1CorrPFIso/i");
  TreeEleTau->Branch("recoLep1DETIso", &m_recoLep1DETIso,"recoLep1DETIso/i");
  TreeEleTau->Branch("recoLep1HEEPIso", &m_recoLep1HEEPIso,"recoLep1HEEPIso/i");
  TreeEleTau->Branch("recoLep1dEtaIn", &m_recoLep1dEtaIn,"recoLep1dEtaIn/i");
  TreeEleTau->Branch("recoLep1dPhiIn", &m_recoLep1dPhiIn,"recoLep1dPhiIn/i");
  TreeEleTau->Branch("recoLep1SigmaIEIE", &m_recoLep1SigmaIEIE,"recoLep1SigmaIEIE/i");
  TreeEleTau->Branch("recoLep1HoverE", &m_recoLep1HoverE,"recoLep1HoverE/i");
  TreeEleTau->Branch("recoLep1dxy", &m_recoLep1dxy,"recoLep1dxy/i");
  TreeEleTau->Branch("recoLep1dz", &m_recoLep1dz,"recoLep1dz/i");
  TreeEleTau->Branch("recoLep1EandP", &m_recoLep1EandP,"recoLep1EandP/i");
  TreeEleTau->Branch("recoLep1Conv1", &m_recoLep1Conv1,"recoLep1Conv1/i");
  TreeEleTau->Branch("recoLep1Conv2", &m_recoLep1Conv2,"recoLep1Conv2/i");
  TreeEleTau->Branch("recoLep1EcalDriven", &m_recoLep1EcalDriven,"recoLep1EcalDriven/i");
  TreeEleTau->Branch("recoLep1dPhiInHEEP", &m_recoLep1dPhiInHEEP,"recoLep1dPhiInHEEP/i");
  TreeEleTau->Branch("recoLep1HoverEHEEP", &m_recoLep1HoverEHEEP,"recoLep1HoverEHEEP/i");
  TreeEleTau->Branch("recoLep1TkSumPt", &m_recoLep1TkSumPt,"recoLep1TkSumPt/i");
  TreeEleTau->Branch("recoLep1LostHits", &m_recoLep1LostHits,"recoLep1LostHits/i");
  TreeEleTau->Branch("recoLep1dEtaInHEEP", &m_recoLep1dEtaInHEEP,"recoLep1dEtaInHEEP/i");
  TreeEleTau->Branch("recoLep1dxyHEEP", &m_recoLep1dxyHEEP,"recoLep1dxyHEEP/i");
  TreeEleTau->Branch("recoLep1e5xe5", &m_recoLep1e5xe5,"recoLep1e5xe5/i");
  TreeEleTau->Branch("recoLep1SigmaIEIEHEEP", &m_recoLep1SigmaIEIEHEEP,"recoLep1SigmaIEIEHEEP/i");
  TreeEleTau->Branch("recoLep1IsGlobal", &m_recoLep1IsGlobal,"recoLep1IsGlobal/i");
  TreeEleTau->Branch("recoLep1IsPFTight", &m_recoLep1IsPFTight,"recoLep1IsPFTight/i");
  TreeEleTau->Branch("recoLep1Chi2Tight", &m_recoLep1Chi2Tight,"recoLep1Chi2Tight/i");
  TreeEleTau->Branch("recoLep1MuonHitsTight", &m_recoLep1MuonHitsTight,"recoLep1MuonHitsTight/i");
  TreeEleTau->Branch("recoLep1MatchesTight", &m_recoLep1MatchesTight,"recoLep1MatchesTight/i");
  TreeEleTau->Branch("recoLep1dxyTight", &m_recoLep1dxyTight,"recoLep1dxyTight/i");
  TreeEleTau->Branch("recoLep1dzTight", &m_recoLep1dzTight,"recoLep1dzTight/i");
  TreeEleTau->Branch("recoLep1PixelHitsTight", &m_recoLep1PixelHitsTight,"recoLep1PixelHitsTight/i");
  TreeEleTau->Branch("recoLep1TrackerLTight", &m_recoLep1TrackerLTight,"recoLep1TrackerLTight/i");
  TreeEleTau->Branch("recoLep1ptErr", &m_recoLep1ptErr,"recoLep1ptErr/i");
  TreeEleTau->Branch("recoLep1MuonHits", &m_recoLep1MuonHits,"recoLep1MuonHits/i");
  TreeEleTau->Branch("recoLep1Matches", &m_recoLep1Matches,"recoLep1Matches/i");
  TreeEleTau->Branch("recoLep1PixelHits", &m_recoLep1PixelHits,"recoLep1PixelHits/i");
  TreeEleTau->Branch("recoLep1TrackerL", &m_recoLep1TrackerL,"recoLep1TrackerL/i");
  TreeEleTau->Branch("recoLep1IsPFMuon", &m_recoLep1IsPFMuon,"recoLep1IsPFMuon/i");
  TreeEleTau->Branch("recoLep1Chi2", &m_recoLep1Chi2,"recoLep1Chi2/i");
  TreeEleTau->Branch("recoLep1Discr1", &m_recoLep1Discr1,"recoLep1Discr1/i");
  TreeEleTau->Branch("recoLep1Discr2", &m_recoLep1Discr2,"recoLep1Discr2/i");
  TreeEleTau->Branch("recoLep1Discr3", &m_recoLep1Discr3,"recoLep1Discr3/i");
  TreeEleTau->Branch("recoLep1Discr4", &m_recoLep1Discr4,"recoLep1Discr4/i");
  TreeEleTau->Branch("recoLep2Pt", &m_recoLep2Pt,"recoLep2Pt/i");
  TreeEleTau->Branch("recoLep2Eta", &m_recoLep2Eta,"recoLep2Eta/i");
  TreeEleTau->Branch("recoLep2Phi", &m_recoLep2Phi,"recoLep2Phi/i");
  TreeEleTau->Branch("recoLep2HEEP", &m_recoLep2HEEP,"recoLep2HEEP/i");
  TreeEleTau->Branch("recoLep2dRJet", &m_recoLep2dRJet,"recoLep2dRJet/i");
  TreeEleTau->Branch("recoLep2PFIso", &m_recoLep2PFIso,"recoLep2PFIso/i");
  TreeEleTau->Branch("recoLep2CorrPFIso", &m_recoLep2CorrPFIso,"recoLep2CorrPFIso/i");
  TreeEleTau->Branch("recoLep2DETIso", &m_recoLep2DETIso,"recoLep2DETIso/i");
  TreeEleTau->Branch("recoLep2HEEPIso", &m_recoLep2HEEPIso,"recoLep2HEEPIso/i");
  TreeEleTau->Branch("recoLep2dEtaIn", &m_recoLep2dEtaIn,"recoLep2dEtaIn/i");
  TreeEleTau->Branch("recoLep2dPhiIn", &m_recoLep2dPhiIn,"recoLep2dPhiIn/i");
  TreeEleTau->Branch("recoLep2SigmaIEIE", &m_recoLep2SigmaIEIE,"recoLep2SigmaIEIE/i");
  TreeEleTau->Branch("recoLep2HoverE", &m_recoLep2HoverE,"recoLep2HoverE/i");
  TreeEleTau->Branch("recoLep2dxy", &m_recoLep2dxy,"recoLep2dxy/i");
  TreeEleTau->Branch("recoLep2dz", &m_recoLep2dz,"recoLep2dz/i");
  TreeEleTau->Branch("recoLep2EandP", &m_recoLep2EandP,"recoLep2EandP/i");
  TreeEleTau->Branch("recoLep2Conv1", &m_recoLep2Conv1,"recoLep2Conv1/i");
  TreeEleTau->Branch("recoLep2Conv2", &m_recoLep2Conv2,"recoLep2Conv2/i");
  TreeEleTau->Branch("recoLep2EcalDriven", &m_recoLep2EcalDriven,"recoLep2EcalDriven/i");
  TreeEleTau->Branch("recoLep2dPhiInHEEP", &m_recoLep2dPhiInHEEP,"recoLep2dPhiInHEEP/i");
  TreeEleTau->Branch("recoLep2HoverEHEEP", &m_recoLep2HoverEHEEP,"recoLep2HoverEHEEP/i");
  TreeEleTau->Branch("recoLep2TkSumPt", &m_recoLep2TkSumPt,"recoLep2TkSumPt/i");
  TreeEleTau->Branch("recoLep2LostHits", &m_recoLep2LostHits,"recoLep2LostHits/i");
  TreeEleTau->Branch("recoLep2dEtaInHEEP", &m_recoLep2dEtaInHEEP,"recoLep2dEtaInHEEP/i");
  TreeEleTau->Branch("recoLep2dxyHEEP", &m_recoLep2dxyHEEP,"recoLep2dxyHEEP/i");
  TreeEleTau->Branch("recoLep2e5xe5", &m_recoLep2e5xe5,"recoLep2e5xe5/i");
  TreeEleTau->Branch("recoLep2SigmaIEIEHEEP", &m_recoLep2SigmaIEIEHEEP,"recoLep2SigmaIEIEHEEP/i");
  TreeEleTau->Branch("recoLep2IsGlobal", &m_recoLep2IsGlobal,"recoLep2IsGlobal/i");
  TreeEleTau->Branch("recoLep2IsPFTight", &m_recoLep2IsPFTight,"recoLep2IsPFTight/i");
  TreeEleTau->Branch("recoLep2Chi2Tight", &m_recoLep2Chi2Tight,"recoLep2Chi2Tight/i");
  TreeEleTau->Branch("recoLep2MuonHitsTight", &m_recoLep2MuonHitsTight,"recoLep2MuonHitsTight/i");
  TreeEleTau->Branch("recoLep2MatchesTight", &m_recoLep2MatchesTight,"recoLep2MatchesTight/i");
  TreeEleTau->Branch("recoLep2dxyTight", &m_recoLep2dxyTight,"recoLep2dxyTight/i");
  TreeEleTau->Branch("recoLep2dzTight", &m_recoLep2dzTight,"recoLep2dzTight/i");
  TreeEleTau->Branch("recoLep2PixelHitsTight", &m_recoLep2PixelHitsTight,"recoLep2PixelHitsTight/i");
  TreeEleTau->Branch("recoLep2TrackerLTight", &m_recoLep2TrackerLTight,"recoLep2TrackerLTight/i");
  TreeEleTau->Branch("recoLep2ptErr", &m_recoLep2ptErr,"recoLep2ptErr/i");
  TreeEleTau->Branch("recoLep2MuonHits", &m_recoLep2MuonHits,"recoLep2MuonHits/i");
  TreeEleTau->Branch("recoLep2Matches", &m_recoLep2Matches,"recoLep2Matches/i");
  TreeEleTau->Branch("recoLep2PixelHits", &m_recoLep2PixelHits,"recoLep2PixelHits/i");
  TreeEleTau->Branch("recoLep2TrackerL", &m_recoLep2TrackerL,"recoLep2TrackerL/i");
  TreeEleTau->Branch("recoLep2IsPFMuon", &m_recoLep2IsPFMuon,"recoLep2IsPFMuon/i");
  TreeEleTau->Branch("recoLep2Chi2", &m_recoLep2Chi2,"recoLep2Chi2/i");

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
			   bool foundGenJet, vector<reco::GenJet>::const_iterator genJet, bool & matched){

  if(foundGenJet){
    float dRGenReco = 9999.;
    pat::JetCollection::const_iterator SelectedJet_prov;
    for(pat::JetCollection::const_iterator jet = CA8JetswithQjets->begin(); jet != CA8JetswithQjets->end(); ++jet) {
      if(ROOT::Math::VectorUtil::DeltaR(jet->p4(),genJet->p4())<0.3 && ROOT::Math::VectorUtil::DeltaR(jet->p4(),genJet->p4())<dRGenReco){
	dRGenReco = ROOT::Math::VectorUtil::DeltaR(jet->p4(),genJet->p4());
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
	 fabs(SelectedJet_prov->eta())<2.4 && 
	 (mass>massMin && mass<massMax) && 
	 (SelectedJet_prov->userFloat("tau2")/SelectedJet_prov->userFloat("tau1"))<0.75){
	foundSelectedJet=true;
	massZ=mass;
	ptZ=SelectedJet_prov->pt();
	tau21Z=SelectedJet_prov->userFloat("tau2")/SelectedJet_prov->userFloat("tau1");
	SelectedJet=SelectedJet_prov;
      }
    }
  }
}

void Efficiency::SelectTau(edm::Handle<pat::TauCollection> tauHandle,
			   pat::JetCollection::const_iterator SelectedJet,
			   pat::TauCollection::const_iterator & SelectedTau, bool foundJet,
			   bool foundGenPart, math::PtEtaPhiELorentzVector genPart, bool & matched){
  if(foundGenPart){
    float dRGenReco = 9999.;
    for (pat::TauCollection::const_iterator patTau = tauHandle->begin(); patTau != tauHandle->end(); ++patTau ) {
      if(ROOT::Math::VectorUtil::DeltaR(patTau->p4(),genPart)<0.3 && ROOT::Math::VectorUtil::DeltaR(patTau->p4(),genPart)<dRGenReco){
	dRGenReco = ROOT::Math::VectorUtil::DeltaR(patTau->p4(),genPart);
	matched=true;
	SelectedTau=patTau;
      }
    }
  }
}

void Efficiency::SelectMuon(edm::Handle<pat::MuonCollection> muoH,
			    pat::JetCollection::const_iterator SelectedJet,
			    pat::MuonCollection::const_iterator & SelectedMuon, bool foundGenMuo, 
			    math::PtEtaPhiELorentzVector genMuo_mt, bool & matched){
  if(foundGenMuo){
    float dRGenReco = 9999.;
    for(pat::MuonCollection::const_iterator muon = muoH->begin(); muon != muoH->end(); ++muon) {
      if(ROOT::Math::VectorUtil::DeltaR(muon->p4(),genMuo_mt)<0.3 && ROOT::Math::VectorUtil::DeltaR(muon->p4(),genMuo_mt)<dRGenReco){
	dRGenReco = ROOT::Math::VectorUtil::DeltaR(muon->p4(),genMuo_mt);
	matched = true;
	SelectedMuon=muon;
      }
    }
  }
}

void Efficiency::SelectElectron(edm::Handle<pat::ElectronCollection> eleH,
				pat::JetCollection::const_iterator SelectedJet,
				pat::ElectronCollection::const_iterator & SelectedEle,
				bool foundGenEle_et, math::PtEtaPhiELorentzVector genEle_et, bool & matched, bool & isoHEEP, float rho){
  if(foundGenEle_et){
    float dRGenReco = 9999.;
    pat::ElectronCollection::const_iterator SelectedEle_prov;
    for(pat::ElectronCollection::const_iterator electron = eleH->begin(); electron != eleH->end(); ++electron) {
      if(ROOT::Math::VectorUtil::DeltaR(electron->p4(),genEle_et)<0.3 && ROOT::Math::VectorUtil::DeltaR(electron->p4(),genEle_et)<dRGenReco){
	dRGenReco = ROOT::Math::VectorUtil::DeltaR(electron->p4(),genEle_et);
	SelectedEle_prov = electron;
	matched = true;
      }
    }
    if(matched){
      SelectedEle=SelectedEle_prov;
      float Et=SelectedEle_prov->caloEnergy()*sin(SelectedEle_prov->p4().theta());
      float eta  = fabs(SelectedEle_prov->caloPosition().eta());
      float iso = 999.;
      float isoCut = 0;
      if(eta<1.442){
	iso = SelectedEle_prov->dr03EcalRecHitSumEt() + SelectedEle_prov->dr03HcalDepth1TowerSumEt();
      	isoCut = 2 + 0.03*Et + 0.28*rho;
	if(iso<isoCut) isoHEEP=true;
      } else if(eta>1.56 && eta<2.5){
	iso = SelectedEle_prov->dr03EcalRecHitSumEt() + SelectedEle_prov->dr03HcalDepth1TowerSumEt();
      	if(Et<=50) isoCut = 2.5 + 0.28*rho;
      	else       isoCut = 2.5 + 0.03*(Et-50.) + 0.28*rho;
	if(iso<isoCut) isoHEEP=true;
      }
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

bool Efficiency::ElectronDETIso(pat::ElectronCollection::const_iterator electron, float rho){
  bool iso = false;
  if(electron->userIso(0) < 5.0){
    bool inBarrel = electron->isEE();
    double ECALIsol = electron->userIso(1);
    double HCALIsol = electron->userIso(2);
    double sumCaloEt = ECALIsol + HCALIsol;
    double sumCaloEtLimit = -1;
    double et = electron->et();
    if(inBarrel) sumCaloEtLimit = 2.0 + 0.03*et + 0.28*rho;
    if(!inBarrel){
      if(et<50) sumCaloEtLimit = 2.5 + 0.28*rho;
      else sumCaloEtLimit = 2.5 + 0.03*(et-50) + 0.28*rho;
    }
    if(sumCaloEt<sumCaloEtLimit) iso = true;
  }
  return iso;
}




void Efficiency::FillTree(int category, TTree *Tree, pat::JetCollection::const_iterator SelectedJet, pat::TauCollection::const_iterator SelectedTau,
			  pat::MuonCollection::const_iterator SelectedMuon,pat::MuonCollection::const_iterator SelectedMuon1,pat::MuonCollection::const_iterator SelectedMuon2,
			  pat::ElectronCollection::const_iterator SelectedElectron, pat::ElectronCollection::const_iterator SelectedElectron1, 
			  pat::ElectronCollection::const_iterator SelectedElectron2, reco::Vertex primaryVertex, math::PtEtaPhiELorentzVector genLep1, 
			  math::PtEtaPhiELorentzVector genLep2, vector<reco::GenJet>::const_iterator genJet, float rho, bool isoHEEP1, bool isoHEEP2,
			  bool foundGenJet, bool foundSelectedJet, bool foundGenLep1, bool foundSelectedLep1,bool foundGenLep2, bool foundSelectedLep2){
  float genJetPt         = -99.;
  float genJetEta        = -99.;
  float genJetPhi        = -99.;
  float genLep1Pt        = -99.;
  float genLep2Pt        = -99.;
  float genLep1Eta       = -99.;
  float genLep2Eta       = -99.;
  float genLep1Phi       = -99.;
  float genLep2Phi       = -99.;
  float recoJetPt        = -99.;
  float recoJetEta       = -99.;
  float recoJetPhi       = -99.;
  float recoDeltaR       = -99.;
  float dRGenRecoLep1	 = -99.;
  float dRGenRecoLep2    = -99.;
  float genDeltaR        = -99.;
  float dRGenRecoJet     = -99.;
  float recoLep1dRJetF   = -99.;
  float recoLep2dRJetF   = -99.;
  float recoLep1PtF      = -99.;
  float recoLep1EtaF     = -99.;
  float recoLep1PhiF     = -99.;
  float recoLep2PtF      = -99.;
  float recoLep2EtaF     = -99.;
  float recoLep2PhiF     = -99.;
  int recoLep1Pt         = 0;
  int recoLep1Eta        = 0;
  int recoLep1Phi        = 0;
  int recoLep1HEEP       = 0;
  int recoLep1dRJet      = 0;
  int recoLep1PFIso      = 0;
  int recoLep1CorrPFIso  = 0;
  int recoLep1DETIso     = 0;
  int recoLep1HEEPIso    = 0;
  int recoLep1dEtaIn     = 0;
  int recoLep1dPhiIn     = 0;
  int recoLep1SigmaIEIE  = 0;
  int recoLep1HoverE     = 0;
  int recoLep1dxy        = 0;
  int recoLep1dz         = 0;
  int recoLep1EandP      = 0;
  int recoLep1Conv1      = 0;
  int recoLep1Conv2      = 0;
  int recoLep1EcalDriven = 0;
  int recoLep1dPhiInHEEP = 0;
  int recoLep1HoverEHEEP = 0;
  int recoLep1TkSumPt    = 0;
  int recoLep1LostHits   = 0;
  int recoLep1dEtaInHEEP = 0;
  int recoLep1dxyHEEP    = 0;
  int recoLep1e5xe5      = 0;
  int recoLep1SigmaIEIEHEEP  = 0;
  int recoLep1IsGlobal       = 0;
  int recoLep1IsPFTight      = 0;
  int recoLep1Chi2Tight      = 0;
  int recoLep1MuonHitsTight  = 0;
  int recoLep1MatchesTight   = 0;
  int recoLep1dxyTight       = 0;
  int recoLep1dzTight        = 0;
  int recoLep1PixelHitsTight = 0;
  int recoLep1TrackerLTight  = 0;
  int recoLep1ptErr      = 0;
  int recoLep1MuonHits   = 0;
  int recoLep1Matches    = 0;
  int recoLep1PixelHits  = 0;
  int recoLep1TrackerL   = 0;
  int recoLep1IsPFMuon   = 0;
  int recoLep1Chi2       = 0;
  int recoLep1Discr1     = 0;
  int recoLep1Discr2     = 0;
  int recoLep1Discr3     = 0;
  int recoLep1Discr4     = 0;
  int recoLep2Pt         = 0;
  int recoLep2Eta        = 0;
  int recoLep2Phi        = 0;
  int recoLep2HEEP       = 0;
  int recoLep2dRJet      = 0;
  int recoLep2PFIso      = 0;
  int recoLep2CorrPFIso  = 0;
  int recoLep2DETIso     = 0;
  int recoLep2HEEPIso    = 0;
  int recoLep2dEtaIn     = 0;
  int recoLep2dPhiIn     = 0;
  int recoLep2SigmaIEIE  = 0;
  int recoLep2HoverE     = 0;
  int recoLep2dxy        = 0;
  int recoLep2dz         = 0;
  int recoLep2EandP      = 0;
  int recoLep2Conv1      = 0;
  int recoLep2Conv2      = 0;
  int recoLep2EcalDriven = 0;
  int recoLep2dPhiInHEEP = 0;
  int recoLep2HoverEHEEP = 0;
  int recoLep2TkSumPt    = 0;
  int recoLep2LostHits   = 0;
  int recoLep2dEtaInHEEP = 0;
  int recoLep2dxyHEEP    = 0;
  int recoLep2e5xe5      = 0;
  int recoLep2SigmaIEIEHEEP  = 0;
  int recoLep2IsGlobal       = 0;
  int recoLep2IsPFTight      = 0;
  int recoLep2Chi2Tight      = 0;
  int recoLep2MuonHitsTight  = 0;
  int recoLep2MatchesTight   = 0;
  int recoLep2dxyTight       = 0;
  int recoLep2dzTight        = 0;
  int recoLep2PixelHitsTight = 0;
  int recoLep2TrackerLTight  = 0;
  int recoLep2ptErr      = 0;
  int recoLep2MuonHits   = 0;
  int recoLep2Matches    = 0;
  int recoLep2PixelHits  = 0;
  int recoLep2TrackerL   = 0;
  int recoLep2IsPFMuon   = 0;
  int recoLep2Chi2       = 0;

  if(foundGenJet){
    genJetPt  = genJet->pt();
    genJetEta = genJet->eta();
    genJetPhi = genJet->phi();
  }
  if(foundGenLep1){
    genLep1Pt     = genLep1.pt();
    genLep1Eta    = genLep1.eta();
    genLep1Phi    = genLep1.phi();
  }
  if(foundGenLep2){
    genLep2Pt     = genLep2.pt();
    genLep2Eta    = genLep2.eta();
    genLep2Phi    = genLep2.phi();
  }
  if(foundGenLep1 && foundGenLep2)    genDeltaR     = ROOT::Math::VectorUtil::DeltaR(genLep1,genLep2);
  if(foundGenJet && foundSelectedJet) dRGenRecoJet  = ROOT::Math::VectorUtil::DeltaR(genJet->p4(),SelectedJet->p4());
  if(foundSelectedJet){
    recoJetPt     = SelectedJet->pt();
    recoJetEta    = SelectedJet->eta();
    recoJetPhi    = SelectedJet->phi();
  }
  
  if(category==0){//ELECTRON+MUON
    //ELECTRON
    if(foundSelectedLep1){
      recoLep1PtF  = SelectedElectron->pt();
      recoLep1EtaF = SelectedElectron->eta();
      recoLep1PhiF = SelectedElectron->phi();
      if(SelectedElectron->pt()>10) recoLep1Pt = 1;
      if(fabs(SelectedElectron->superCluster()->eta())<2.5) recoLep1Eta = 1;
      if(fabs(SelectedElectron->phi())<3.2) recoLep1Phi = 1;
      if(SelectedElectron->userInt("HEEPId")==0) recoLep1HEEP = 1;
      if(ElectronPFIso(SelectedElectron, rho)<0.1) recoLep1PFIso = 1;
      if(ElectronCorrPFIso(SelectedElectron, rho)<0.1) recoLep1CorrPFIso = 1;
      if(ElectronDETIso(SelectedElectron, rho)) recoLep1DETIso = 1;
      if(isoHEEP1) recoLep1HEEPIso = 1;
      if(fabs(SelectedElectron->gsfTrack()->dxy(primaryVertex.position()))<0.02) recoLep1dxy = 1;
      if(fabs(SelectedElectron->gsfTrack()->dz(primaryVertex.position()))<0.1) recoLep1dz = 1;
      if((fabs(1/SelectedElectron->ecalEnergy() - SelectedElectron->eSuperClusterOverP()/SelectedElectron->ecalEnergy()))<0.05) recoLep1EandP = 1;
      if(SelectedElectron->passConversionVeto()!=0) recoLep1Conv1 = 1;
      if(SelectedElectron->gsfTrack()->trackerExpectedHitsInner().numberOfHits()==0) recoLep1Conv2 = 1;
      if(SelectedElectron->ecalDriven()==1) recoLep1EcalDriven=1;////////-----1
      if(fabs(SelectedElectron->deltaPhiSuperClusterTrackAtVtx())<0.060) recoLep1dPhiInHEEP = 1;////////-----3
      if(SelectedElectron->hadronicOverEm()<0.05)                        recoLep1HoverEHEEP = 1;////////-----4
      if(SelectedElectron->dr03TkSumPt()<5) recoLep1TkSumPt=1;////////-----6
      if(SelectedElectron->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits()<2) recoLep1LostHits=1;////////-----7
      if(fabs(SelectedElectron->superCluster()->eta())<=1.479){
	if(fabs(SelectedElectron->deltaEtaSuperClusterTrackAtVtx())<0.004) recoLep1dEtaIn = 1;
	if(fabs(SelectedElectron->deltaPhiSuperClusterTrackAtVtx())<0.030) recoLep1dPhiIn = 1;
	if(SelectedElectron->sigmaIetaIeta()<0.01)                         recoLep1SigmaIEIE = 1;
	if(SelectedElectron->hadronicOverEm()<0.12)                        recoLep1HoverE = 1;
	if(fabs(SelectedElectron->deltaEtaSuperClusterTrackAtVtx())<0.005) recoLep1dEtaInHEEP = 1;////////-----2
	if(fabs(SelectedElectron->gsfTrack()->dxy(primaryVertex.position()))<0.02) recoLep1dxyHEEP = 1;////////-----8
	if((SelectedElectron->e2x5Max()/SelectedElectron->e5x5()>0.94 || SelectedElectron->e1x5()/SelectedElectron->e5x5()>0.83)) recoLep1e5xe5 = 1;////////-----5
	recoLep1SigmaIEIEHEEP = 1;////////-----9
      }
      if(fabs(SelectedElectron->superCluster()->eta())>1.479 && fabs(SelectedElectron->superCluster()->eta())<2.5){
	if(fabs(SelectedElectron->deltaEtaSuperClusterTrackAtVtx())<0.005) recoLep1dEtaIn = 1;
	if(fabs(SelectedElectron->deltaPhiSuperClusterTrackAtVtx())<0.020) recoLep1dPhiIn = 1;
	if(SelectedElectron->sigmaIetaIeta()<0.03)                       recoLep1SigmaIEIE = 1;
	if(SelectedElectron->hadronicOverEm()<0.10)                      recoLep1HoverE = 1;
	if(fabs(SelectedElectron->deltaEtaSuperClusterTrackAtVtx())<0.007) recoLep1dEtaInHEEP = 1;////////-----2
	if(fabs(SelectedElectron->gsfTrack()->dxy(primaryVertex.position()))<0.05) recoLep1dxyHEEP = 1;////////-----8
	recoLep1e5xe5 = 1;////////-----5
	if(SelectedElectron->sigmaIetaIeta()<0.03)                       recoLep1SigmaIEIEHEEP = 1;////////-----9
      }
      if(foundSelectedJet) { 
	recoLep1dRJetF = ROOT::Math::VectorUtil::DeltaR(SelectedElectron->p4(),SelectedJet->p4());
	if(ROOT::Math::VectorUtil::DeltaR(SelectedElectron->p4(),SelectedJet->p4())>0.8) recoLep1dRJet = 1;
      }
      if(foundGenLep1)     dRGenRecoLep1 = ROOT::Math::VectorUtil::DeltaR(genLep1, SelectedElectron->p4());
    }
    if(foundSelectedLep2){
      //MUON
      recoLep2PtF  = SelectedMuon->pt();
      recoLep2EtaF = SelectedMuon->eta();
      recoLep2PhiF = SelectedMuon->phi();
      if(SelectedMuon->pt()>10) recoLep2Pt = 1;
      if(fabs(SelectedMuon->eta())<2.4) recoLep2Eta = 1;
      if(fabs(SelectedMuon->phi())<3.2) recoLep2Phi = 1;
      if(MuonPFIso(SelectedMuon, true)<0.2) recoLep2PFIso = 1;
      if(MuonCorrPFIso(SelectedMuon, true)<0.2) recoLep2CorrPFIso = 1;
      if(SelectedMuon->isGlobalMuon()) recoLep2IsGlobal = 1;
      if(SelectedMuon->isGlobalMuon()){
	if(SelectedMuon->isPFMuon()) recoLep2IsPFTight = 1;
	if(SelectedMuon->globalTrack()->normalizedChi2() < 10.) recoLep2Chi2Tight = 1;
	if(SelectedMuon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0) recoLep2MuonHitsTight = 1;
	if(SelectedMuon->numberOfMatches() > 1) recoLep2MatchesTight = 1;
	if(fabs(SelectedMuon->muonBestTrack()->dxy(primaryVertex.position())) < 0.2 ) recoLep2dxyTight = 1;
	if(fabs(SelectedMuon->muonBestTrack()->dz(primaryVertex.position())) < 0.5) recoLep2dzTight = 1;
	if(SelectedMuon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0) recoLep2PixelHitsTight = 1;
	if(SelectedMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5) recoLep2TrackerLTight = 1;
	reco::TrackRef cktTrack = (muon::tevOptimized(*SelectedMuon, 200, 17., 40., 0.25)).first;
	if((cktTrack->ptError()/cktTrack->pt())<0.3) recoLep2ptErr = 1;
	if(SelectedMuon->globalTrack()->hitPattern().numberOfValidMuonHits()>0) recoLep2MuonHits = 1;
	if(SelectedMuon->numberOfMatches()>1) recoLep2Matches = 1;
	if(fabs(cktTrack->dxy(primaryVertex.position()))<0.2) recoLep2dxy = 1;
	if(fabs(cktTrack->dz(primaryVertex.position()))<0.5) recoLep2dz = 1;
	if(SelectedMuon->innerTrack()->hitPattern().numberOfValidPixelHits()>0) recoLep2PixelHits = 1;
	if(SelectedMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement()>5) recoLep2TrackerL = 1;
	if(SelectedMuon->isPFMuon()) recoLep2IsPFMuon = 1;
	if(SelectedMuon->globalTrack()->normalizedChi2()<10) recoLep2Chi2 = 1;
      }
      if(foundSelectedJet) { 
	recoLep2dRJetF = ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedJet->p4());
	if(ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedJet->p4())>0.8) recoLep2dRJet = 1;
      }
      if(foundGenLep2)     dRGenRecoLep2 = ROOT::Math::VectorUtil::DeltaR(genLep2, SelectedMuon->p4());
    }
    if(foundSelectedLep1 && foundSelectedLep2) recoDeltaR = ROOT::Math::VectorUtil::DeltaR(SelectedElectron->p4(),SelectedMuon->p4());

  } else if(category==1){//MUON+MUON
    if(foundSelectedLep1){
      //MUON1
      recoLep1PtF  = SelectedMuon1->pt();
      recoLep1EtaF = SelectedMuon1->eta();
      recoLep1PhiF = SelectedMuon1->phi();
      if(SelectedMuon1->pt()>10) recoLep1Pt = 1;
      if(fabs(SelectedMuon1->eta())<2.4) recoLep1Eta = 1;
      if(fabs(SelectedMuon1->phi())<3.2) recoLep1Phi = 1;
      if(MuonPFIso(SelectedMuon1, true)<0.2) recoLep1PFIso = 1;
      if(MuonCorrPFIso(SelectedMuon1, true)<0.2) recoLep1CorrPFIso = 1;
      if(foundSelectedLep2)  { if(MuonDETIso(SelectedMuon1,SelectedMuon2,true)<0.1) recoLep1DETIso = 1;}
      if(!foundSelectedLep2) recoLep1DETIso = 1;
      if(SelectedMuon1->isGlobalMuon()) recoLep1IsGlobal = 1;
      if(SelectedMuon1->isGlobalMuon()){
	if(SelectedMuon1->isPFMuon()) recoLep1IsPFTight = 1;
	if(SelectedMuon1->globalTrack()->normalizedChi2() < 10.) recoLep1Chi2Tight = 1;
	if(SelectedMuon1->globalTrack()->hitPattern().numberOfValidMuonHits() > 0) recoLep1MuonHitsTight = 1;
	if(SelectedMuon1->numberOfMatches() > 1) recoLep1MatchesTight = 1;
	if(fabs(SelectedMuon1->muonBestTrack()->dxy(primaryVertex.position())) < 0.2 ) recoLep1dxyTight = 1;
	if(fabs(SelectedMuon1->muonBestTrack()->dz(primaryVertex.position())) < 0.5) recoLep1dzTight = 1;
	if(SelectedMuon1->innerTrack()->hitPattern().numberOfValidPixelHits() > 0) recoLep1PixelHitsTight = 1;
	if(SelectedMuon1->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5) recoLep1TrackerLTight = 1;
	reco::TrackRef cktTrack1 = (muon::tevOptimized(*SelectedMuon1, 200, 17., 40., 0.25)).first;
	if((cktTrack1->ptError()/cktTrack1->pt())<0.3) recoLep1ptErr = 1;
	if(SelectedMuon1->globalTrack()->hitPattern().numberOfValidMuonHits()>0) recoLep1MuonHits = 1;
	if(SelectedMuon1->numberOfMatches()>1) recoLep1Matches = 1;
	if(fabs(cktTrack1->dxy(primaryVertex.position()))<0.2) recoLep1dxy = 1;
	if(fabs(cktTrack1->dz(primaryVertex.position()))<0.5) recoLep1dz = 1;
	if(SelectedMuon1->innerTrack()->hitPattern().numberOfValidPixelHits()>0) recoLep1PixelHits = 1;
	if(SelectedMuon1->innerTrack()->hitPattern().trackerLayersWithMeasurement()>5) recoLep1TrackerL = 1;
	if(SelectedMuon1->isPFMuon()) recoLep1IsPFMuon = 1;
	if(SelectedMuon1->globalTrack()->normalizedChi2()<10) recoLep1Chi2 = 1;
      }
      if(foundSelectedJet) { 
	recoLep1dRJetF = ROOT::Math::VectorUtil::DeltaR(SelectedMuon1->p4(),SelectedJet->p4());
	if(ROOT::Math::VectorUtil::DeltaR(SelectedMuon1->p4(),SelectedJet->p4())>0.8) recoLep1dRJet = 1;
      }
      if(foundGenLep1)     dRGenRecoLep1 = ROOT::Math::VectorUtil::DeltaR(genLep1, SelectedMuon1->p4());
    }
    if(foundSelectedLep2){
      //MUON2
      recoLep2PtF  = SelectedMuon2->pt();
      recoLep2EtaF = SelectedMuon2->eta();
      recoLep2PhiF = SelectedMuon2->phi();
      if(SelectedMuon2->pt()>10) recoLep2Pt = 1;
      if(fabs(SelectedMuon2->eta())<2.4) recoLep2Eta = 1;
      if(fabs(SelectedMuon2->phi())<3.2) recoLep2Phi = 1;
      if(MuonPFIso(SelectedMuon2, true)<0.2) recoLep2PFIso = 1;
      if(MuonCorrPFIso(SelectedMuon2, true)<0.2) recoLep2CorrPFIso = 1;
      if(foundSelectedLep1){if(MuonDETIso(SelectedMuon2,SelectedMuon1,true)<0.1) recoLep2DETIso = 1;}
      if(!foundSelectedLep1) recoLep2DETIso=1;
      if(SelectedMuon2->isGlobalMuon()) recoLep2IsGlobal = 1;
      if(SelectedMuon2->isGlobalMuon()){
	if(SelectedMuon2->isPFMuon()) recoLep2IsPFTight = 1;
	if(SelectedMuon2->globalTrack()->normalizedChi2() < 10.) recoLep2Chi2Tight = 1;
	if(SelectedMuon2->globalTrack()->hitPattern().numberOfValidMuonHits() > 0) recoLep2MuonHitsTight = 1;
	if(SelectedMuon2->numberOfMatches() > 1) recoLep2MatchesTight = 1;
	if(fabs(SelectedMuon2->muonBestTrack()->dxy(primaryVertex.position())) < 0.2 ) recoLep2dxyTight = 1;
	if(fabs(SelectedMuon2->muonBestTrack()->dz(primaryVertex.position())) < 0.5) recoLep2dzTight = 1;
	if(SelectedMuon2->innerTrack()->hitPattern().numberOfValidPixelHits() > 0) recoLep2PixelHitsTight = 1;
	if(SelectedMuon2->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5) recoLep2TrackerLTight = 1;
	reco::TrackRef cktTrack2 = (muon::tevOptimized(*SelectedMuon2, 200, 17., 40., 0.25)).first;
	if((cktTrack2->ptError()/cktTrack2->pt())<0.3) recoLep2ptErr = 1;
	if(SelectedMuon2->globalTrack()->hitPattern().numberOfValidMuonHits()>0) recoLep2MuonHits = 1;
	if(SelectedMuon2->numberOfMatches()>1) recoLep2Matches = 1;
	if(fabs(cktTrack2->dxy(primaryVertex.position()))<0.2) recoLep2dxy = 1;
	if(fabs(cktTrack2->dz(primaryVertex.position()))<0.5) recoLep2dz = 1;
	if(SelectedMuon2->innerTrack()->hitPattern().numberOfValidPixelHits()>0) recoLep2PixelHits = 1;
	if(SelectedMuon2->innerTrack()->hitPattern().trackerLayersWithMeasurement()>5) recoLep2TrackerL = 1;
	if(SelectedMuon2->isPFMuon()) recoLep2IsPFMuon = 1;
	if(SelectedMuon2->globalTrack()->normalizedChi2()<10) recoLep2Chi2 = 1;
      }
      if(foundSelectedJet) { 
	recoLep2dRJetF = ROOT::Math::VectorUtil::DeltaR(SelectedMuon2->p4(),SelectedJet->p4());
	if(ROOT::Math::VectorUtil::DeltaR(SelectedMuon2->p4(),SelectedJet->p4())>0.8) recoLep2dRJet = 1;
      }
      if(foundGenLep2)     dRGenRecoLep2 = ROOT::Math::VectorUtil::DeltaR(genLep2, SelectedMuon2->p4());
    }
    if(foundSelectedLep1 && foundSelectedLep2) recoDeltaR=ROOT::Math::VectorUtil::DeltaR(SelectedMuon1->p4(),SelectedMuon2->p4());

    
  } else if(category==2){//ELECTRON+ELECTRON
    if(foundSelectedLep1){
      //ELECTRON
      recoLep1PtF  = SelectedElectron1->pt();
      recoLep1EtaF = SelectedElectron1->eta();
      recoLep1PhiF = SelectedElectron1->phi();
      if(SelectedElectron1->pt()>10) recoLep1Pt = 1;
      if(fabs(SelectedElectron1->superCluster()->eta())<2.5) recoLep1Eta = 1;
      if(fabs(SelectedElectron1->phi())<3.2) recoLep1Phi = 1;
      if(SelectedElectron1->userInt("HEEPId")==0) recoLep1HEEP = 1;
      if(ElectronPFIso(SelectedElectron1, rho)<0.1) recoLep1PFIso = 1;
      if(ElectronCorrPFIso(SelectedElectron1, rho)<0.1) recoLep1CorrPFIso = 1;
      if(ElectronDETIso(SelectedElectron1, rho)) recoLep1DETIso = 1;
      if(isoHEEP1) recoLep1HEEPIso = 1;
      if(fabs(SelectedElectron1->gsfTrack()->dxy(primaryVertex.position()))<0.02) recoLep1dxy = 1;
      if(fabs(SelectedElectron1->gsfTrack()->dz(primaryVertex.position()))<0.1) recoLep1dz = 1;
      if((fabs(1/SelectedElectron1->ecalEnergy() - SelectedElectron1->eSuperClusterOverP()/SelectedElectron1->ecalEnergy()))<0.05) recoLep1EandP = 1;
      if(SelectedElectron1->passConversionVeto()!=0) recoLep1Conv1 = 1;
      if(SelectedElectron1->gsfTrack()->trackerExpectedHitsInner().numberOfHits()==0) recoLep1Conv2 = 1;
      if(SelectedElectron1->ecalDriven()==1) recoLep1EcalDriven=1;////////-----1
      if(fabs(SelectedElectron1->deltaPhiSuperClusterTrackAtVtx())<0.060) recoLep1dPhiInHEEP = 1;////////-----3
      if(SelectedElectron1->hadronicOverEm()<0.05)                        recoLep1HoverEHEEP = 1;////////-----4
      if(SelectedElectron1->dr03TkSumPt()<5) recoLep1TkSumPt=1;////////-----6
      if(SelectedElectron1->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits()<2) recoLep1LostHits=1;////////-----7
      if(fabs(SelectedElectron1->superCluster()->eta())<=1.479){
	if(fabs(SelectedElectron1->deltaEtaSuperClusterTrackAtVtx())<0.004) recoLep1dEtaIn = 1;
	if(fabs(SelectedElectron1->deltaPhiSuperClusterTrackAtVtx())<0.030) recoLep1dPhiIn = 1;
	if(SelectedElectron1->sigmaIetaIeta()<0.01) recoLep1SigmaIEIE = 1;
	if(SelectedElectron1->hadronicOverEm()<0.12) recoLep1HoverE = 1;
	if(fabs(SelectedElectron1->deltaEtaSuperClusterTrackAtVtx())<0.005) recoLep1dEtaInHEEP = 1;////////-----2
	if(fabs(SelectedElectron1->gsfTrack()->dxy(primaryVertex.position()))<0.02) recoLep1dxyHEEP = 1;////////-----8
	if((SelectedElectron1->e2x5Max()/SelectedElectron1->e5x5()>0.94 || SelectedElectron1->e1x5()/SelectedElectron1->e5x5()>0.83)) recoLep1e5xe5 = 1;////////-----5
	recoLep1SigmaIEIEHEEP = 1;////////-----9
      }
      if(fabs(SelectedElectron1->superCluster()->eta())>1.479 && fabs(SelectedElectron1->superCluster()->eta())<2.5){
	if(fabs(SelectedElectron1->deltaEtaSuperClusterTrackAtVtx())<0.005) recoLep1dEtaIn = 1;
	if(fabs(SelectedElectron1->deltaPhiSuperClusterTrackAtVtx())<0.020) recoLep1dPhiIn = 1;
	if(SelectedElectron1->sigmaIetaIeta()<0.03) recoLep1SigmaIEIE = 1;
	if(SelectedElectron1->hadronicOverEm()<0.10) recoLep1HoverE = 1;
	if(fabs(SelectedElectron1->deltaEtaSuperClusterTrackAtVtx())<0.007) recoLep1dEtaInHEEP = 1;////////-----2
	if(fabs(SelectedElectron1->gsfTrack()->dxy(primaryVertex.position()))<0.05) recoLep1dxyHEEP = 1;////////-----8
	recoLep1e5xe5 = 1;////////-----5
	if(SelectedElectron1->sigmaIetaIeta()<0.03)                       recoLep1SigmaIEIEHEEP = 1;////////-----9
      }
      if(foundSelectedJet) {
	recoLep1dRJetF = ROOT::Math::VectorUtil::DeltaR(SelectedElectron1->p4(),SelectedJet->p4());
	if(ROOT::Math::VectorUtil::DeltaR(SelectedElectron1->p4(),SelectedJet->p4())>0.8) recoLep1dRJet = 1;
      }
      if(foundGenLep1)     dRGenRecoLep1 = ROOT::Math::VectorUtil::DeltaR(genLep1, SelectedElectron1->p4());
    }
    if(foundSelectedLep2){
      //ELECTRON2
      recoLep2PtF  = SelectedElectron2->pt();
      recoLep2EtaF = SelectedElectron2->eta();
      recoLep2PhiF = SelectedElectron2->phi();
      if(SelectedElectron2->pt()>10) recoLep2Pt = 1;
      if(fabs(SelectedElectron2->superCluster()->eta())<2.5) recoLep2Eta = 1;
      if(fabs(SelectedElectron2->phi())<3.2) recoLep2Phi = 1;
      if(SelectedElectron2->userInt("HEEPId")==0) recoLep2HEEP = 1;
      if(ElectronPFIso(SelectedElectron2, rho)<0.1) recoLep2PFIso = 1;
      if(ElectronCorrPFIso(SelectedElectron2, rho)<0.1) recoLep2CorrPFIso = 1;
      if(ElectronDETIso(SelectedElectron2, rho)) recoLep2DETIso = 1;
      if(isoHEEP2) recoLep2HEEPIso = 1;
      if(fabs(SelectedElectron2->gsfTrack()->dxy(primaryVertex.position()))<0.02) recoLep2dxy = 1;
      if(fabs(SelectedElectron2->gsfTrack()->dz(primaryVertex.position()))<0.1) recoLep2dz = 1;
      if((fabs(1/SelectedElectron2->ecalEnergy() - SelectedElectron2->eSuperClusterOverP()/SelectedElectron2->ecalEnergy()))<0.05) recoLep2EandP = 1;
      if(SelectedElectron2->passConversionVeto()!=0) recoLep2Conv1 = 1;
      if(SelectedElectron2->gsfTrack()->trackerExpectedHitsInner().numberOfHits()==0) recoLep2Conv2 = 1;
      if(SelectedElectron2->ecalDriven()==1) recoLep2EcalDriven=1;////////-----1
      if(fabs(SelectedElectron2->deltaPhiSuperClusterTrackAtVtx())<0.060) recoLep2dPhiInHEEP = 1;////////-----3
      if(SelectedElectron2->hadronicOverEm()<0.05)                        recoLep2HoverEHEEP = 1;////////-----4
      if(SelectedElectron2->dr03TkSumPt()<5) recoLep2TkSumPt=1;////////-----6
      if(SelectedElectron2->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits()<2) recoLep2LostHits=1;////////-----7
      if(fabs(SelectedElectron2->superCluster()->eta())<=1.479){
	if(fabs(SelectedElectron2->deltaEtaSuperClusterTrackAtVtx())<0.004) recoLep2dEtaIn = 1;
	if(fabs(SelectedElectron2->deltaPhiSuperClusterTrackAtVtx())<0.030) recoLep2dPhiIn = 1;
	if(SelectedElectron2->sigmaIetaIeta()<0.01) recoLep2SigmaIEIE = 1;
	if(SelectedElectron2->hadronicOverEm()<0.12) recoLep2HoverE = 1;
	if(fabs(SelectedElectron2->deltaEtaSuperClusterTrackAtVtx())<0.005) recoLep2dEtaInHEEP = 1;////////-----2
	if(fabs(SelectedElectron2->gsfTrack()->dxy(primaryVertex.position()))<0.02) recoLep2dxyHEEP = 1;////////-----8
	if((SelectedElectron2->e2x5Max()/SelectedElectron2->e5x5()>0.94 || SelectedElectron2->e1x5()/SelectedElectron2->e5x5()>0.83)) recoLep2e5xe5 = 1;////////-----5
	recoLep2SigmaIEIEHEEP = 1;////////-----9
      }
      if(fabs(SelectedElectron2->superCluster()->eta())>1.479 && fabs(SelectedElectron2->superCluster()->eta())<2.5){
	if(fabs(SelectedElectron2->deltaEtaSuperClusterTrackAtVtx())<0.005) recoLep2dEtaIn = 1;
	if(fabs(SelectedElectron2->deltaPhiSuperClusterTrackAtVtx())<0.020) recoLep2dPhiIn = 1;
	if(SelectedElectron2->sigmaIetaIeta()<0.03) recoLep2SigmaIEIE = 1;
	if(SelectedElectron2->hadronicOverEm()<0.10) recoLep2HoverE = 1;
	if(fabs(SelectedElectron2->deltaEtaSuperClusterTrackAtVtx())<0.007) recoLep2dEtaInHEEP = 1;////////-----2
	if(fabs(SelectedElectron2->gsfTrack()->dxy(primaryVertex.position()))<0.05) recoLep2dxyHEEP = 1;////////-----8
	recoLep2e5xe5 = 1;////////-----5
	if(SelectedElectron2->sigmaIetaIeta()<0.03)                       recoLep2SigmaIEIEHEEP = 1;////////-----9
      }
      if(foundSelectedJet) {
	recoLep2dRJetF = ROOT::Math::VectorUtil::DeltaR(SelectedElectron2->p4(),SelectedJet->p4());
	if(ROOT::Math::VectorUtil::DeltaR(SelectedElectron2->p4(),SelectedJet->p4())>0.8) recoLep2dRJet = 1;
      }
      if(foundGenLep2)     dRGenRecoLep2 = ROOT::Math::VectorUtil::DeltaR(genLep2, SelectedElectron2->p4());
    }
    if(foundSelectedLep1 && foundSelectedLep2) recoDeltaR=ROOT::Math::VectorUtil::DeltaR(SelectedElectron1->p4(),SelectedElectron2->p4());

    
  } else if(category==3){//MUON+TAU
    if(foundSelectedLep1){
      //TAU
      recoLep1PtF  = SelectedTau->pt();
      recoLep1EtaF = SelectedTau->eta();
      recoLep1PhiF = SelectedTau->phi();
      if(SelectedTau->pt()>20) recoLep1Pt = 1;
      if(fabs(SelectedTau->eta())<2.3) recoLep1Eta = 1;
      if(fabs(SelectedTau->phi())<3.2) recoLep1Phi = 1;
      if(SelectedTau->tauID("decayModeFindingNewDMs")>0.5) recoLep1Discr1 = 1;
      if (SelectedTau->tauID("againstMuonLoose")>0.5)recoLep1Discr2 = 1;
      if(SelectedTau->tauID("againstElectronLoose")>0.5) recoLep1Discr3 = 1;
      if(SelectedTau->tauID("byVLooseIsolationMVA3newDMwoLT")>0.5) recoLep1Discr4 = 1;
      if(foundSelectedJet){ 
	recoLep1dRJetF = ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(),SelectedJet->p4());
	if(ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(), SelectedJet->p4())>0.8) recoLep1dRJet = 1;
      }
      if(foundGenLep1)     dRGenRecoLep1 = ROOT::Math::VectorUtil::DeltaR(genLep1, SelectedTau->p4());
    }
    if(foundSelectedLep2){
      //MUON
      recoLep2PtF  = SelectedMuon->pt();
      recoLep2EtaF = SelectedMuon->eta();
      recoLep2PhiF = SelectedMuon->phi();
      if(SelectedMuon->pt()>10) recoLep2Pt = 1;
      if(fabs(SelectedMuon->eta())<2.4) recoLep2Eta = 1;
      if(fabs(SelectedMuon->phi())<3.2) recoLep2Phi = 1;
      if(MuonPFIso(SelectedMuon, true)<0.2) recoLep2PFIso = 1;
      if(MuonCorrPFIso(SelectedMuon, true)<0.2) recoLep2CorrPFIso = 1;
      if(SelectedMuon->isGlobalMuon()) recoLep2IsGlobal = 1;
      if(SelectedMuon->isGlobalMuon()){
	if(SelectedMuon->isPFMuon()) recoLep2IsPFTight = 1;
	if(SelectedMuon->globalTrack()->normalizedChi2() < 10.) recoLep2Chi2Tight = 1;
	if(SelectedMuon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0) recoLep2MuonHitsTight = 1;
	if(SelectedMuon->numberOfMatches() > 1) recoLep2MatchesTight = 1;
	if(fabs(SelectedMuon->muonBestTrack()->dxy(primaryVertex.position())) < 0.2 ) recoLep2dxyTight = 1;
	if(fabs(SelectedMuon->muonBestTrack()->dz(primaryVertex.position())) < 0.5) recoLep2dzTight = 1;
	if(SelectedMuon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0) recoLep2PixelHitsTight = 1;
	if(SelectedMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5) recoLep2TrackerLTight = 1;
	reco::TrackRef cktTrack = (muon::tevOptimized(*SelectedMuon, 200, 17., 40., 0.25)).first;
	if((cktTrack->ptError()/cktTrack->pt())<0.3) recoLep2ptErr = 1;
	if(SelectedMuon->globalTrack()->hitPattern().numberOfValidMuonHits()>0) recoLep2MuonHits = 1;
	if(SelectedMuon->numberOfMatches()>1) recoLep2Matches = 1;
	if(fabs(cktTrack->dxy(primaryVertex.position()))<0.2) recoLep2dxy = 1;
	if(fabs(cktTrack->dz(primaryVertex.position()))<0.5) recoLep2dz = 1;
	if(SelectedMuon->innerTrack()->hitPattern().numberOfValidPixelHits()>0) recoLep2PixelHits = 1;
	if(SelectedMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement()>5) recoLep2TrackerL = 1;
	if(SelectedMuon->isPFMuon()) recoLep2IsPFMuon = 1;
	if(SelectedMuon->globalTrack()->normalizedChi2()<10) recoLep2Chi2 = 1;
      }
      if(foundSelectedJet){ 
	recoLep2dRJetF = ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedJet->p4());
	if(ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedJet->p4())>0.8) recoLep2dRJet = 1;
      }
      if(foundGenLep2)     dRGenRecoLep2 = ROOT::Math::VectorUtil::DeltaR(genLep2, SelectedMuon->p4());
    }
    if(foundSelectedLep1 && foundSelectedLep2) recoDeltaR=ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(),SelectedMuon->p4());

    
  } else if(category==4){//ELECTRON+TAU
    if(foundSelectedLep1){
      //TAU
      recoLep1PtF  = SelectedTau->pt();
      recoLep1EtaF = SelectedTau->eta();
      recoLep1PhiF = SelectedTau->phi();
      if(SelectedTau->pt()>20) recoLep1Pt = 1;
      if(fabs(SelectedTau->eta())<2.3) recoLep1Eta = 1;
      if(fabs(SelectedTau->phi())<3.2) recoLep1Phi = 1;
      if(SelectedTau->tauID("decayModeFindingNewDMs")>0.5) recoLep1Discr1 = 1;
      if(SelectedTau->tauID("againstMuonLoose")>0.5) recoLep1Discr2 = 1;
      if(SelectedTau->tauID("againstElectronLoose")>0.5) recoLep1Discr3 = 1;
      if(SelectedTau->tauID("byVLooseIsolationMVA3newDMwoLT")>0.5) recoLep1Discr4 = 1;
      if(foundSelectedJet){ 
	recoLep1dRJetF = ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(),SelectedJet->p4());
	if(ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(),SelectedJet->p4())>0.8) recoLep1dRJet = 1;
      }
      if(foundGenLep1)     dRGenRecoLep1 = ROOT::Math::VectorUtil::DeltaR(genLep1, SelectedTau->p4());
    }
    if(foundSelectedLep2){
      //ELECTRON
      recoLep2PtF  = SelectedElectron->pt();
      recoLep2EtaF = SelectedElectron->eta();
      recoLep2PhiF = SelectedElectron->phi();
      if(SelectedElectron->pt()>10) recoLep2Pt = 1;
      if(fabs(SelectedElectron->superCluster()->eta())<2.5) recoLep2Eta = 1;
      if(fabs(SelectedElectron->phi())<3.2) recoLep2Phi = 1;
      if(SelectedElectron->userInt("HEEPId")==0) recoLep2HEEP = 1;
      if(ElectronPFIso(SelectedElectron, rho)<0.1) recoLep2PFIso = 1;
      if(ElectronCorrPFIso(SelectedElectron, rho)<0.1) recoLep2CorrPFIso = 1;
      if(ElectronDETIso(SelectedElectron, rho)) recoLep2DETIso = 1;
      if(isoHEEP2) recoLep2HEEPIso = 1;
      if(fabs(SelectedElectron->gsfTrack()->dxy(primaryVertex.position()))<0.02) recoLep2dxy = 1;
      if(fabs(SelectedElectron->gsfTrack()->dz(primaryVertex.position()))<0.1) recoLep2dz = 1;
      if((fabs(1/SelectedElectron->ecalEnergy() - SelectedElectron->eSuperClusterOverP()/SelectedElectron->ecalEnergy()))<0.05) recoLep2EandP = 1;
      if(SelectedElectron->passConversionVeto()!=0) recoLep2Conv1 = 1;
      if(SelectedElectron->gsfTrack()->trackerExpectedHitsInner().numberOfHits()==0) recoLep2Conv2 = 1;
      if(SelectedElectron->ecalDriven()==1) recoLep2EcalDriven=1;////////-----1
      if(fabs(SelectedElectron->deltaPhiSuperClusterTrackAtVtx())<0.060) recoLep2dPhiInHEEP = 1;////////-----3
      if(SelectedElectron->hadronicOverEm()<0.05)                        recoLep2HoverEHEEP = 1;////////-----4
      if(SelectedElectron->dr03TkSumPt()<5) recoLep2TkSumPt=1;////////-----6
      if(SelectedElectron->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits()<2) recoLep2LostHits=1;////////-----7
      if(fabs(SelectedElectron->superCluster()->eta())<=1.479){
	if(fabs(SelectedElectron->deltaEtaSuperClusterTrackAtVtx())<0.004) recoLep2dEtaIn = 1;
	if(fabs(SelectedElectron->deltaPhiSuperClusterTrackAtVtx())<0.030) recoLep2dPhiIn = 1;
	if(SelectedElectron->sigmaIetaIeta()<0.01) recoLep2SigmaIEIE = 1;
	if(SelectedElectron->hadronicOverEm()<0.12) recoLep2HoverE = 1;
	if(fabs(SelectedElectron->deltaEtaSuperClusterTrackAtVtx())<0.005) recoLep2dEtaInHEEP = 1;////////-----2
	if(fabs(SelectedElectron->gsfTrack()->dxy(primaryVertex.position()))<0.02) recoLep2dxyHEEP = 1;////////-----8
	if((SelectedElectron->e2x5Max()/SelectedElectron->e5x5()>0.94 || SelectedElectron->e1x5()/SelectedElectron->e5x5()>0.83)) recoLep2e5xe5 = 1;////////-----5
	recoLep2SigmaIEIEHEEP = 1;////////-----9
      }
      if(fabs(SelectedElectron->superCluster()->eta())>1.479 && fabs(SelectedElectron->superCluster()->eta())<2.5){
	if(fabs(SelectedElectron->deltaEtaSuperClusterTrackAtVtx())<0.005) recoLep2dEtaIn = 1;
	if(fabs(SelectedElectron->deltaPhiSuperClusterTrackAtVtx())<0.020) recoLep2dPhiIn = 1;
	if(SelectedElectron->sigmaIetaIeta()<0.03) recoLep2SigmaIEIE = 1;
	if(SelectedElectron->hadronicOverEm()<0.10) recoLep2HoverE = 1;
	if(fabs(SelectedElectron->deltaEtaSuperClusterTrackAtVtx())<0.007) recoLep2dEtaInHEEP = 1;////////-----2
	if(fabs(SelectedElectron->gsfTrack()->dxy(primaryVertex.position()))<0.05) recoLep2dxyHEEP = 1;////////-----8
	recoLep2e5xe5 = 1;////////-----5
	if(SelectedElectron->sigmaIetaIeta()<0.03)                       recoLep2SigmaIEIEHEEP = 1;////////-----9
      }
      if(foundSelectedJet){ 
	recoLep2dRJetF = ROOT::Math::VectorUtil::DeltaR(SelectedElectron->p4(),SelectedJet->p4());
	if(ROOT::Math::VectorUtil::DeltaR(SelectedElectron->p4(),SelectedJet->p4())>0.8) recoLep2dRJet = 1;
      }
      if(foundGenLep2)     dRGenRecoLep2 = ROOT::Math::VectorUtil::DeltaR(genLep2, SelectedElectron->p4());
    }
    if(foundSelectedLep1 && foundSelectedLep2) recoDeltaR = ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(),SelectedElectron->p4());
  }

  m_genJetPt          = (float)genJetPt;
  m_genJetEta         = (float)genJetEta;
  m_genJetPhi         = (float)genJetPhi;
  m_genLep1Pt         = (float)genLep1Pt;
  m_genLep2Pt         = (float)genLep2Pt;
  m_genLep1Eta        = (float)genLep1Eta;
  m_genLep2Eta        = (float)genLep2Eta;
  m_genLep1Phi        = (float)genLep1Phi;
  m_genLep2Phi        = (float)genLep2Phi;
  m_recoJetPt         = (float)recoJetPt;
  m_recoJetEta        = (float)recoJetEta;
  m_recoJetPhi        = (float)recoJetPhi;
  m_recoDeltaR        = (float)recoDeltaR;
  m_dRGenRecoLep1     = (float)dRGenRecoLep1;	
  m_dRGenRecoLep2     = (float)dRGenRecoLep2;
  m_genDeltaR         = (float)genDeltaR;
  m_dRGenRecoJet      = (float)dRGenRecoJet;
  m_recoLep1dRJetF    = (float)recoLep1dRJetF;
  m_recoLep2dRJetF    = (float)recoLep2dRJetF;
  m_recoLep1PtF       = (float)recoLep1PtF;
  m_recoLep1EtaF      = (float)recoLep1EtaF;
  m_recoLep1PhiF      = (float)recoLep1PhiF;
  m_recoLep2PtF       = (float)recoLep2PtF;
  m_recoLep2EtaF      = (float)recoLep2EtaF;
  m_recoLep2PhiF      = (float)recoLep2PhiF;
  m_recoLep1Pt        = (int)recoLep1Pt;
  m_recoLep1Eta       = (int)recoLep1Eta;
  m_recoLep1Phi       = (int)recoLep1Phi;
  m_recoLep1HEEP      = (int)recoLep1HEEP;
  m_recoLep1dRJet     = (int)recoLep1dRJet;
  m_recoLep1PFIso     = (int)recoLep1PFIso;
  m_recoLep1CorrPFIso = (int)recoLep1CorrPFIso;
  m_recoLep1DETIso    = (int)recoLep1DETIso;
  m_recoLep1HEEPIso   = (int)recoLep1HEEPIso;
  m_recoLep1dEtaIn    = (int)recoLep1dEtaIn;
  m_recoLep1dPhiIn    = (int)recoLep1dPhiIn;
  m_recoLep1SigmaIEIE = (int)recoLep1SigmaIEIE;
  m_recoLep1HoverE    = (int)recoLep1HoverE;
  m_recoLep1dxy       = (int)recoLep1dxy;
  m_recoLep1dz        = (int)recoLep1dz;
  m_recoLep1EandP     = (int)recoLep1EandP;
  m_recoLep1Conv1     = (int)recoLep1Conv1;
  m_recoLep1Conv2     = (int)recoLep1Conv2;
  m_recoLep1EcalDriven= (int)recoLep1EcalDriven;
  m_recoLep1dPhiInHEEP= (int)recoLep1dPhiInHEEP;
  m_recoLep1HoverEHEEP= (int)recoLep1HoverEHEEP;
  m_recoLep1TkSumPt   = (int)recoLep1TkSumPt;
  m_recoLep1LostHits  = (int)recoLep1LostHits;
  m_recoLep1dEtaInHEEP= (int)recoLep1dEtaInHEEP;
  m_recoLep1dxyHEEP   = (int)recoLep1dxyHEEP;
  m_recoLep1e5xe5     = (int)recoLep1e5xe5;
  m_recoLep1SigmaIEIEHEEP= (int)recoLep1SigmaIEIEHEEP;
  m_recoLep1IsGlobal  = (int)recoLep1IsGlobal;
  m_recoLep1IsPFTight      = (int)recoLep1IsPFTight     ;
  m_recoLep1Chi2Tight      = (int)recoLep1Chi2Tight     ;
  m_recoLep1MuonHitsTight  = (int)recoLep1MuonHitsTight ;
  m_recoLep1MatchesTight   = (int)recoLep1MatchesTight  ;
  m_recoLep1dxyTight       = (int)recoLep1dxyTight      ;
  m_recoLep1dzTight        = (int)recoLep1dzTight       ;
  m_recoLep1PixelHitsTight = (int)recoLep1PixelHitsTight;
  m_recoLep1TrackerLTight  = (int)recoLep1TrackerLTight ;
  m_recoLep1ptErr     = (int)recoLep1ptErr;
  m_recoLep1MuonHits  = (int)recoLep1MuonHits;
  m_recoLep1Matches   = (int)recoLep1Matches;
  m_recoLep1PixelHits = (int)recoLep1PixelHits;
  m_recoLep1TrackerL  = (int)recoLep1TrackerL;
  m_recoLep1IsPFMuon  = (int)recoLep1IsPFMuon;
  m_recoLep1Chi2      = (int)recoLep1Chi2;
  m_recoLep1Discr1    = (int)recoLep1Discr1;
  m_recoLep1Discr2    = (int)recoLep1Discr2;
  m_recoLep1Discr3    = (int)recoLep1Discr3;
  m_recoLep1Discr4    = (int)recoLep1Discr4;
  m_recoLep2Pt        = (int)recoLep2Pt;
  m_recoLep2Eta       = (int)recoLep2Eta;
  m_recoLep2Phi       = (int)recoLep2Phi;
  m_recoLep2HEEP      = (int)recoLep2HEEP;
  m_recoLep2dRJet     = (int)recoLep2dRJet;
  m_recoLep2PFIso     = (int)recoLep2PFIso;
  m_recoLep2CorrPFIso = (int)recoLep2CorrPFIso;
  m_recoLep2DETIso    = (int)recoLep2DETIso;
  m_recoLep2HEEPIso   = (int)recoLep2HEEPIso;
  m_recoLep2dEtaIn    = (int)recoLep2dEtaIn;
  m_recoLep2dPhiIn    = (int)recoLep2dPhiIn;
  m_recoLep2SigmaIEIE = (int)recoLep2SigmaIEIE;
  m_recoLep2HoverE    = (int)recoLep2HoverE;
  m_recoLep2dxy       = (int)recoLep2dxy;
  m_recoLep2dz        = (int)recoLep2dz;
  m_recoLep2EandP     = (int)recoLep2EandP;
  m_recoLep2Conv1     = (int)recoLep2Conv1;
  m_recoLep2Conv2     = (int)recoLep2Conv2;
  m_recoLep2EcalDriven= (int)recoLep2EcalDriven;
  m_recoLep2dPhiInHEEP= (int)recoLep2dPhiInHEEP;
  m_recoLep2HoverEHEEP= (int)recoLep2HoverEHEEP;
  m_recoLep2TkSumPt   = (int)recoLep2TkSumPt;
  m_recoLep2LostHits  = (int)recoLep2LostHits;
  m_recoLep2dEtaInHEEP= (int)recoLep2dEtaInHEEP;
  m_recoLep2dxyHEEP   = (int)recoLep2dxyHEEP;
  m_recoLep2e5xe5     = (int)recoLep2e5xe5;
  m_recoLep2SigmaIEIEHEEP= (int)recoLep2SigmaIEIEHEEP;
  m_recoLep2IsGlobal  = (int)recoLep2IsGlobal;
  m_recoLep2IsPFTight      = (int)recoLep2IsPFTight     ;
  m_recoLep2Chi2Tight      = (int)recoLep2Chi2Tight     ;
  m_recoLep2MuonHitsTight  = (int)recoLep2MuonHitsTight ;
  m_recoLep2MatchesTight   = (int)recoLep2MatchesTight  ;
  m_recoLep2dxyTight       = (int)recoLep2dxyTight      ;
  m_recoLep2dzTight        = (int)recoLep2dzTight       ;
  m_recoLep2PixelHitsTight = (int)recoLep2PixelHitsTight;
  m_recoLep2TrackerLTight  = (int)recoLep2TrackerLTight ;
  m_recoLep2ptErr     = (int)recoLep2ptErr;
  m_recoLep2MuonHits  = (int)recoLep2MuonHits;
  m_recoLep2Matches   = (int)recoLep2Matches;
  m_recoLep2PixelHits = (int)recoLep2PixelHits;
  m_recoLep2TrackerL  = (int)recoLep2TrackerL;
  m_recoLep2IsPFMuon  = (int)recoLep2IsPFMuon;
  m_recoLep2Chi2      = (int)recoLep2Chi2;
  Tree->Fill();
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(Efficiency);
