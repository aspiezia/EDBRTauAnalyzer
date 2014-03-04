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

  //FUNCTION
  void SelectJet(edm::Handle<pat::JetCollection> CA8JetswithQjets, edm::Handle<pat::JetCollection> CA8JetsPruned, bool & foundJet,
                 pat::JetCollection::const_iterator & SelectedJet, float & massZ, float & tau21Z, float & ptZ, float massMin, float massMax);
  void SelectTau(edm::Handle<pat::TauCollection> tauHandle, pat::JetCollection::const_iterator SelectedJet, bool & foundTau, 
		 pat::TauCollection::const_iterator & SelectedTau, float & ptTau, bool foundJet);
  void SelectMuon(edm::Handle<pat::MuonCollection> muoH, pat::JetCollection::const_iterator SelectedJet, bool & foundMuon,
		  pat::MuonCollection::const_iterator & SelectedMuon, float & ptMuon, bool foundJet, reco::Vertex primaryVertex);
  void SelectElectron(edm::Handle<pat::ElectronCollection> eleH, pat::JetCollection::const_iterator SelectedJet, bool & foundElectron,
		      pat::ElectronCollection::const_iterator & SelectedElectron, float & ptElectron, bool foundJet, reco::Vertex primaryVertex);
  void Efficiency(float & genEvent, bool isData, edm::Handle<std::vector<reco::GenParticle> > genParts);
  void svfit(edm::Handle<pat::METCollection> metRaw, edm::Handle<pat::METCollection> met, LorentzVector SelectedTau, LorentzVector SelectedMuon, TLorentzVector PrunedJet,
	     float & MassSVFit, float & XMassSVFit, float & dRJetZSVFit, float & ptSVFit, pat::JetCollection::const_iterator SelectedJet);
  
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
  TTree *TreeSB1MuoTau;
  TTree *TreeSB1EleTau;
  TTree *TreeSB2MuoTau;
  TTree *TreeSB2EleTau;
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
  float m_PtSvfit;
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
SemiLeptonicAnalyzer::SemiLeptonicAnalyzer(const edm::ParameterSet& iConfig)

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
  Efficiency(genEvent, isData, genParts);
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

  
  //JET SELECTION - SR
  pat::JetCollection::const_iterator SelectedJet;
  float massZ=-9999; float ptZ=-999; bool foundJet=false; float tau21Z=-9999;
  SelectJet(CA8JetswithQjets, CA8JetsPruned, foundJet, SelectedJet, massZ, tau21Z, ptZ, 70, 110);
  
  //JET SELECTION - SB1
  pat::JetCollection::const_iterator SelectedSB1Jet;
  float massSB1Z=-9999; float ptSB1Z=-999; bool foundSB1Jet=false; float tau21SB1Z=-9999;
  SelectJet(CA8JetswithQjets, CA8JetsPruned, foundSB1Jet, SelectedSB1Jet, massSB1Z, tau21SB1Z, ptSB1Z, 20, 70);
  
  //JET SELECTION - SB2
  pat::JetCollection::const_iterator SelectedSB2Jet;
  float massSB2Z=-9999; float ptSB2Z=-999; bool foundSB2Jet=false; float tau21SB2Z=-9999;
  SelectJet(CA8JetswithQjets, CA8JetsPruned, foundSB2Jet, SelectedSB2Jet, massSB2Z, tau21SB2Z, ptSB2Z, 110, 99999);
  
  //TAU SELECTION - MUTAU - SR
  float ptTauMuTau=-99; bool foundTauMuTau=false;
  pat::TauCollection::const_iterator SelectedTauMuTau;
  SelectTau(tauMuTauHandle, SelectedJet, foundTauMuTau, SelectedTauMuTau, ptTauMuTau, foundJet);
  
  //TAU SELECTION - MUTAU - SB1
  float ptSB1TauMuTau=-99; bool foundSB1TauMuTau=false;
  pat::TauCollection::const_iterator SelectedSB1TauMuTau;
  SelectTau(tauMuTauHandle, SelectedSB1Jet, foundSB1TauMuTau, SelectedSB1TauMuTau, ptSB1TauMuTau, foundSB1Jet);
  
  //TAU SELECTION - MUTAU - SB2
  float ptSB2TauMuTau=-99; bool foundSB2TauMuTau=false;
  pat::TauCollection::const_iterator SelectedSB2TauMuTau;
  SelectTau(tauMuTauHandle, SelectedSB2Jet, foundSB2TauMuTau, SelectedSB2TauMuTau, ptSB2TauMuTau, foundSB2Jet); 
  
  //TAU SELECTION - ELTAU - SR
  float ptTauElTau=-99; bool foundTauElTau=false;
  pat::TauCollection::const_iterator SelectedTauElTau;
  SelectTau(tauElTauHandle, SelectedJet, foundTauElTau, SelectedTauElTau, ptTauElTau, foundJet);
  
  //TAU SELECTION - ELTAU - SB1
  float ptSB1TauElTau=-99; bool foundSB1TauElTau=false;
  pat::TauCollection::const_iterator SelectedSB1TauElTau;
  SelectTau(tauElTauHandle, SelectedSB1Jet, foundSB1TauElTau, SelectedSB1TauElTau, ptSB1TauElTau, foundSB1Jet);
  
  //TAU SELECTION - ELTAU - SB2
  float ptSB2TauElTau=-99; bool foundSB2TauElTau=false;
  pat::TauCollection::const_iterator SelectedSB2TauElTau;
  SelectTau(tauElTauHandle, SelectedSB2Jet, foundSB2TauElTau, SelectedSB2TauElTau, ptSB2TauElTau, foundSB2Jet); 

  //MUON SELECTION - SR
  float ptMuon=-99; bool foundMuon=false;
  pat::MuonCollection::const_iterator SelectedMuon;
  SelectMuon(muoH, SelectedJet, foundMuon, SelectedMuon, ptMuon, foundJet, primaryVertex);

  //MUON SELECTION - SB1
  float ptSB1Muon=-99; bool foundSB1Muon=false;
  pat::MuonCollection::const_iterator SelectedSB1Muon;
  SelectMuon(muoH, SelectedSB1Jet, foundSB1Muon, SelectedSB1Muon, ptSB1Muon, foundSB1Jet, primaryVertex);

  //MUON SELECTION - SB2
  float ptSB2Muon=-99; bool foundSB2Muon=false;
  pat::MuonCollection::const_iterator SelectedSB2Muon;
  SelectMuon(muoH, SelectedSB2Jet, foundSB2Muon, SelectedSB2Muon, ptSB2Muon, foundSB2Jet, primaryVertex);

  //ELECTRON SELECTION - SR
  float ptElectron=-99; bool foundElectron=false;
  pat::ElectronCollection::const_iterator SelectedElectron;
  SelectElectron(eleH, SelectedJet, foundElectron, SelectedElectron, ptElectron, foundJet, primaryVertex);

  //ELECTRON SELECTION - SB1
  float ptSB1Electron=-99; bool foundSB1Electron=false;
  pat::ElectronCollection::const_iterator SelectedSB1Electron;
  SelectElectron(eleH, SelectedSB1Jet, foundSB1Electron, SelectedSB1Electron, ptSB1Electron, foundSB1Jet, primaryVertex);

  //ELECTRON SELECTION - SB2
  float ptSB2Electron=-99; bool foundSB2Electron=false;
  pat::ElectronCollection::const_iterator SelectedSB2Electron;
  SelectElectron(eleH, SelectedSB2Jet, foundSB2Electron, SelectedSB2Electron, ptSB2Electron, foundSB2Jet, primaryVertex);
  
  //BTAG VETO - SR
  int nbtagsL=0; int nbtagsM=0; int nbtagsT=0;
  if(foundJet) BtagVeto(ak5jetCands, nbtagsL, nbtagsM, nbtagsT, SelectedJet);

  //BTAG VETO - SB1
  int nSB1btagsL=0; int nSB1btagsM=0; int nSB1btagsT=0;
  if(foundSB1Jet) BtagVeto(ak5jetCands, nSB1btagsL, nSB1btagsM, nSB1btagsT, SelectedSB1Jet);

  //BTAG VETO - SB2
  int nSB2btagsL=0; int nSB2btagsM=0; int nSB2btagsT=0;
  if(foundSB2Jet) BtagVeto(ak5jetCands, nSB2btagsL, nSB2btagsM, nSB2btagsT, SelectedSB2Jet);

  //MU-TAU ANALYSIS
  if(foundJet && foundTauMuTau && foundMuon){ 
    //SVFIT
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedJet->pt(),SelectedJet->eta(),SelectedJet->phi(),massZ);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());
    float MassSVFit = -1;
    float XMassSVFit = -1;
    float dRJetZSVFit = -1;
    float ptSVFit = -1;
    svfit(metRaw, met, SelectedTauMuTau->p4(), SelectedMuon->p4(), PrunedJet, MassSVFit, XMassSVFit, dRJetZSVFit, ptSVFit, SelectedJet);
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
    m_PtSvfit=ptSVFit;
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
    //SVFIT
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedJet->pt(),SelectedJet->eta(),SelectedJet->phi(),massZ);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());
    float MassSVFit = -1;
    float XMassSVFit = -1;
    float dRJetZSVFit = -1;
    float ptSVFit = -1;
    svfit(metRaw, met, SelectedTauElTau->p4(), SelectedElectron->p4(), PrunedJet, MassSVFit, XMassSVFit, dRJetZSVFit, ptSVFit, SelectedJet);
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
    m_PtSvfit=ptSVFit;
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

  //MU-TAU ANALYSIS - SB1
  if(foundSB1Jet && foundSB1TauMuTau && foundSB1Muon){ 
    //SVFIT
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedSB1Jet->pt(),SelectedSB1Jet->eta(),SelectedSB1Jet->phi(),massSB1Z);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());
    float MassSVFit = -1;
    float XMassSVFit = -1;
    float dRJetZSVFit = -1;
    float ptSVFit = -1;
    svfit(metRaw, met, SelectedSB1TauMuTau->p4(), SelectedSB1Muon->p4(), PrunedJet, MassSVFit, XMassSVFit, dRJetZSVFit, ptSVFit, SelectedSB1Jet);
    //COLLINEAR APPROXIMATION
    TLorentzVector CATauTau; bool CA = false;
    float a = (SelectedSB1TauMuTau->py()*met->begin()->px()-SelectedSB1TauMuTau->px()*met->begin()->py())/
      (SelectedSB1Muon->px()*SelectedSB1TauMuTau->py()-SelectedSB1Muon->py()*SelectedSB1TauMuTau->px());
    float b = (SelectedSB1Muon->py()*met->begin()->px()-SelectedSB1Muon->px()*met->begin()->py())/
      (SelectedSB1TauMuTau->px()*SelectedSB1Muon->py()-SelectedSB1TauMuTau->py()*SelectedSB1Muon->px());
    if(((1+a)*(1+b))>0) {
      CATauTau.SetPxPyPzE((1+a)*SelectedSB1Muon->px()+(1+b)*SelectedSB1TauMuTau->px(),
			  (1+a)*SelectedSB1Muon->py()+(1+b)*SelectedSB1TauMuTau->py(),
			  (1+a)*SelectedSB1Muon->pz()+(1+b)*SelectedSB1TauMuTau->pz(),
			  (1+a)*SelectedSB1Muon->energy()+(1+b)*SelectedSB1TauMuTau->energy());
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
    math::PtEtaPhiELorentzVector dilep; dilep = SelectedSB1TauMuTau->p4()+SelectedSB1Muon->p4();
    math::PtEtaPhiELorentzVector dilepmet; dilepmet = SelectedSB1TauMuTau->p4()+SelectedSB1Muon->p4()+met->begin()->p4();
    math::PtEtaPhiELorentzVector dilepjet; dilepjet = SelectedSB1TauMuTau->p4()+SelectedSB1Muon->p4()+PrunedJet_prov;
    math::PtEtaPhiELorentzVector dilepmetjet; dilepmetjet = SelectedSB1TauMuTau->p4()+SelectedSB1Muon->p4()+met->begin()->p4()+PrunedJet_prov;
    m_jetPt=SelectedSB1Jet->pt();
    m_jetEta=SelectedSB1Jet->eta();
    m_jetMass=massSB1Z;
    m_jetSubjettiness=tau21SB1Z;
    m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB1Jet->p4());
    m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB1Jet->p4());
    m_dRJetLep=ROOT::Math::VectorUtil::DeltaR(SelectedSB1Muon->p4(),SelectedSB1Jet->p4());
    m_dRJetTau=ROOT::Math::VectorUtil::DeltaR(SelectedSB1TauMuTau->p4(),SelectedSB1Jet->p4());
    m_dPhiTauMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB1TauMuTau->p4());
    m_dRTauMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB1TauMuTau->p4());
    m_dRTauLep=ROOT::Math::VectorUtil::DeltaR(SelectedSB1Muon->p4(),SelectedSB1TauMuTau->p4());
    m_dPhiLepMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB1Muon->p4());
    m_dRLepMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB1Muon->p4());
    m_dRZZVis=ROOT::Math::VectorUtil::DeltaR(dilep,PrunedJet);
    m_dRZZEff=ROOT::Math::VectorUtil::DeltaR(dilepmet,PrunedJet);
    m_dRZZSvFit=dRJetZSVFit;
    m_dRZZCA=dRJetZCA;
    m_tauPt=SelectedSB1TauMuTau->pt();
    m_tauEta=SelectedSB1TauMuTau->eta();
    m_tauCharge=SelectedSB1TauMuTau->charge();
    m_lepPt=SelectedSB1Muon->pt();
    m_lepEta=SelectedSB1Muon->eta();
    m_lepCharge=SelectedSB1Muon->charge();
    m_lepPFIso=MuonPFIso(SelectedSB1Muon,true);
    m_lepCorrPFIso=MuonCorrPFIso(SelectedSB1Muon,true);
    m_charge=SelectedSB1TauMuTau->charge()*SelectedSB1Muon->charge();
    m_met=met->begin()->pt();
    m_metPhi=met->begin()->phi();
    m_uncorrmet=uncorrmet->begin()->pt();
    m_uncorrmetPhi=uncorrmet->begin()->phi();
    m_PtSvfit=ptSVFit;
    m_MassVis=dilep.mass();
    m_MassEff=dilepmet.mass();
    m_MassSvfit=MassSVFit;
    m_MassCA=MassCA;
    m_XMassVis=dilepjet.mass();
    m_XMassEff=dilepmetjet.mass();
    m_XMassSVFit=XMassSVFit;
    m_XMassCA=XmassCA;
    m_nbtagsL=nSB1btagsL;
    m_nbtagsM=nSB1btagsM;
    m_nbtagsT=nSB1btagsT;
    m_trigger=(int)isFired_HLT;
    m_trigger320=(int)isFired_HLT_PFJet320;
    m_trigger650=(int)isFired_HLT_HT650;
    m_sideband=(int)(massSB1Z<70);
    m_NVertices=vertices->size();
    m_PUWeight=MyWeight;
    m_metPx=met->begin()->px();
    m_metPy=met->begin()->py();
    m_NeventsTOT=NeventsTOT_;
    m_xsec=xsec_;
    m_lumi=lumi_;
    m_weight=xsec_*lumi_/NeventsTOT_;
    m_genEvent = genEvent;
    TreeSB1MuoTau->Fill();
  }

  //ELE-TAU ANALYSIS - SB1
  if(foundSB1Jet && foundSB1TauElTau && foundSB1Electron){
    //SVFIT
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedSB1Jet->pt(),SelectedSB1Jet->eta(),SelectedSB1Jet->phi(),massSB1Z);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());
    float MassSVFit = -1;
    float XMassSVFit = -1;
    float dRJetZSVFit = -1;
    float ptSVFit = -1;
    svfit(metRaw, met, SelectedSB1TauElTau->p4(), SelectedSB1Electron->p4(), PrunedJet, MassSVFit, XMassSVFit, dRJetZSVFit, ptSVFit, SelectedSB1Jet);
    //COLLINEAR APPROXIMATION
    TLorentzVector CATauTau; bool CA = false;
    float a = (SelectedSB1TauElTau->py()*met->begin()->px()-SelectedSB1TauElTau->px()*met->begin()->py())/
      (SelectedSB1Electron->px()*SelectedSB1TauElTau->py()-SelectedSB1Electron->py()*SelectedSB1TauElTau->px());
    float b = (SelectedSB1Electron->py()*met->begin()->px()-SelectedSB1Electron->px()*met->begin()->py())/
      (SelectedSB1TauElTau->px()*SelectedSB1Electron->py()-SelectedSB1TauElTau->py()*SelectedSB1Electron->px());
    if(((1+a)*(1+b))>0) {
      CATauTau.SetPxPyPzE((1+a)*SelectedSB1Electron->px()+(1+b)*SelectedSB1TauElTau->px(),
			  (1+a)*SelectedSB1Electron->py()+(1+b)*SelectedSB1TauElTau->py(),
			  (1+a)*SelectedSB1Electron->pz()+(1+b)*SelectedSB1TauElTau->pz(),
			  (1+a)*SelectedSB1Electron->energy()+(1+b)*SelectedSB1TauElTau->energy());
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
    math::PtEtaPhiELorentzVector dilep; dilep = SelectedSB1TauElTau->p4()+SelectedSB1Electron->p4();
    math::PtEtaPhiELorentzVector dilepmet; dilepmet = SelectedSB1TauElTau->p4()+SelectedSB1Electron->p4()+met->begin()->p4();
    math::PtEtaPhiELorentzVector dilepjet; dilepjet = SelectedSB1TauElTau->p4()+SelectedSB1Electron->p4()+PrunedJet_prov;
    math::PtEtaPhiELorentzVector dilepmetjet; dilepmetjet = SelectedSB1TauElTau->p4()+SelectedSB1Electron->p4()+met->begin()->p4()+PrunedJet_prov;
    m_jetPt=SelectedSB1Jet->pt();
    m_jetEta=SelectedSB1Jet->eta();
    m_jetMass=massSB1Z;
    m_jetSubjettiness=tau21SB1Z;
    m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB1Jet->p4());
    m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB1Jet->p4());
    m_dRJetLep=ROOT::Math::VectorUtil::DeltaR(SelectedSB1Electron->p4(),SelectedSB1Jet->p4());
    m_dRJetTau=ROOT::Math::VectorUtil::DeltaR(SelectedSB1TauElTau->p4(),SelectedSB1Jet->p4());
    m_dPhiTauMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB1TauElTau->p4());
    m_dRTauMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB1TauElTau->p4());
    m_dRTauLep=ROOT::Math::VectorUtil::DeltaR(SelectedSB1Electron->p4(),SelectedSB1TauElTau->p4());
    m_dPhiLepMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB1Electron->p4());
    m_dRLepMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB1Electron->p4());
    m_dRZZVis=ROOT::Math::VectorUtil::DeltaR(dilep,PrunedJet);
    m_dRZZEff=ROOT::Math::VectorUtil::DeltaR(dilepmet,PrunedJet);
    m_dRZZSvFit=dRJetZSVFit;
    m_dRZZCA=dRJetZCA;
    m_tauPt=SelectedSB1TauElTau->pt();
    m_tauEta=SelectedSB1TauElTau->eta();
    m_tauCharge=SelectedSB1TauElTau->charge();
    m_lepPt=SelectedSB1Electron->pt();
    m_lepEta=SelectedSB1Electron->eta();
    m_lepCharge=SelectedSB1Electron->charge();
    m_lepPFIso=ElectronPFIso(SelectedSB1Electron,rho);
    m_lepCorrPFIso=ElectronCorrPFIso(SelectedSB1Electron,rho);
    m_charge=SelectedSB1TauElTau->charge()*SelectedSB1Electron->charge();
    m_met=met->begin()->pt();
    m_metPhi=met->begin()->phi();
    m_uncorrmet=uncorrmet->begin()->pt();
    m_uncorrmetPhi=uncorrmet->begin()->phi();
    m_PtSvfit=ptSVFit;
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
    m_sideband=(int)(massSB1Z<70);
    m_NVertices=vertices->size();
    m_PUWeight=MyWeight;
    m_metPx=met->begin()->px();
    m_metPy=met->begin()->py();
    m_NeventsTOT=NeventsTOT_;
    m_xsec=xsec_;
    m_lumi=lumi_;
    m_weight=xsec_*lumi_/NeventsTOT_;
    m_genEvent = genEvent;
    TreeSB1EleTau->Fill();
  }

  //MU-TAU ANALYSIS - SB2
  if(foundSB2Jet && foundSB2TauMuTau && foundSB2Muon){ 
    //SVFIT
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedSB2Jet->pt(),SelectedSB2Jet->eta(),SelectedSB2Jet->phi(),massSB2Z);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());
    float MassSVFit = -1;
    float XMassSVFit = -1;
    float dRJetZSVFit = -1;
    float ptSVFit = -1;
    svfit(metRaw, met, SelectedSB2TauMuTau->p4(), SelectedSB2Muon->p4(), PrunedJet, MassSVFit, XMassSVFit, dRJetZSVFit, ptSVFit, SelectedSB2Jet);
    //COLLINEAR APPROXIMATION
    TLorentzVector CATauTau; bool CA = false;
    float a = (SelectedSB2TauMuTau->py()*met->begin()->px()-SelectedSB2TauMuTau->px()*met->begin()->py())/
      (SelectedSB2Muon->px()*SelectedSB2TauMuTau->py()-SelectedSB2Muon->py()*SelectedSB2TauMuTau->px());
    float b = (SelectedSB2Muon->py()*met->begin()->px()-SelectedSB2Muon->px()*met->begin()->py())/
      (SelectedSB2TauMuTau->px()*SelectedSB2Muon->py()-SelectedSB2TauMuTau->py()*SelectedSB2Muon->px());
    if(((1+a)*(1+b))>0) {
      CATauTau.SetPxPyPzE((1+a)*SelectedSB2Muon->px()+(1+b)*SelectedSB2TauMuTau->px(),
			  (1+a)*SelectedSB2Muon->py()+(1+b)*SelectedSB2TauMuTau->py(),
			  (1+a)*SelectedSB2Muon->pz()+(1+b)*SelectedSB2TauMuTau->pz(),
			  (1+a)*SelectedSB2Muon->energy()+(1+b)*SelectedSB2TauMuTau->energy());
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
    math::PtEtaPhiELorentzVector dilep; dilep = SelectedSB2TauMuTau->p4()+SelectedSB2Muon->p4();
    math::PtEtaPhiELorentzVector dilepmet; dilepmet = SelectedSB2TauMuTau->p4()+SelectedSB2Muon->p4()+met->begin()->p4();
    math::PtEtaPhiELorentzVector dilepjet; dilepjet = SelectedSB2TauMuTau->p4()+SelectedSB2Muon->p4()+PrunedJet_prov;
    math::PtEtaPhiELorentzVector dilepmetjet; dilepmetjet = SelectedSB2TauMuTau->p4()+SelectedSB2Muon->p4()+met->begin()->p4()+PrunedJet_prov;
    m_jetPt=SelectedSB2Jet->pt();
    m_jetEta=SelectedSB2Jet->eta();
    m_jetMass=massSB2Z;
    m_jetSubjettiness=tau21SB2Z;
    m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB2Jet->p4());
    m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB2Jet->p4());
    m_dRJetLep=ROOT::Math::VectorUtil::DeltaR(SelectedSB2Muon->p4(),SelectedSB2Jet->p4());
    m_dRJetTau=ROOT::Math::VectorUtil::DeltaR(SelectedSB2TauMuTau->p4(),SelectedSB2Jet->p4());
    m_dPhiTauMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB2TauMuTau->p4());
    m_dRTauMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB2TauMuTau->p4());
    m_dRTauLep=ROOT::Math::VectorUtil::DeltaR(SelectedSB2Muon->p4(),SelectedSB2TauMuTau->p4());
    m_dPhiLepMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB2Muon->p4());
    m_dRLepMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB2Muon->p4());
    m_dRZZVis=ROOT::Math::VectorUtil::DeltaR(dilep,PrunedJet);
    m_dRZZEff=ROOT::Math::VectorUtil::DeltaR(dilepmet,PrunedJet);
    m_dRZZSvFit=dRJetZSVFit;
    m_dRZZCA=dRJetZCA;
    m_tauPt=SelectedSB2TauMuTau->pt();
    m_tauEta=SelectedSB2TauMuTau->eta();
    m_tauCharge=SelectedSB2TauMuTau->charge();
    m_lepPt=SelectedSB2Muon->pt();
    m_lepEta=SelectedSB2Muon->eta();
    m_lepCharge=SelectedSB2Muon->charge();
    m_lepPFIso=MuonPFIso(SelectedSB2Muon,true);
    m_lepCorrPFIso=MuonCorrPFIso(SelectedSB2Muon,true);
    m_charge=SelectedSB2TauMuTau->charge()*SelectedSB2Muon->charge();
    m_met=met->begin()->pt();
    m_metPhi=met->begin()->phi();
    m_uncorrmet=uncorrmet->begin()->pt();
    m_uncorrmetPhi=uncorrmet->begin()->phi();
    m_PtSvfit=ptSVFit;
    m_MassVis=dilep.mass();
    m_MassEff=dilepmet.mass();
    m_MassSvfit=MassSVFit;
    m_MassCA=MassCA;
    m_XMassVis=dilepjet.mass();
    m_XMassEff=dilepmetjet.mass();
    m_XMassSVFit=XMassSVFit;
    m_XMassCA=XmassCA;
    m_nbtagsL=nSB2btagsL;
    m_nbtagsM=nSB2btagsM;
    m_nbtagsT=nSB2btagsT;
    m_trigger=(int)isFired_HLT;
    m_trigger320=(int)isFired_HLT_PFJet320;
    m_trigger650=(int)isFired_HLT_HT650;
    m_sideband=(int)(massSB2Z<70);
    m_NVertices=vertices->size();
    m_PUWeight=MyWeight;
    m_metPx=met->begin()->px();
    m_metPy=met->begin()->py();
    m_NeventsTOT=NeventsTOT_;
    m_xsec=xsec_;
    m_lumi=lumi_;
    m_weight=xsec_*lumi_/NeventsTOT_;
    m_genEvent = genEvent;
    TreeSB2MuoTau->Fill();
  }

  //ELE-TAU ANALYSIS - SB2
  if(foundSB2Jet && foundSB2TauElTau && foundSB2Electron){
    //SVFIT
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedSB2Jet->pt(),SelectedSB2Jet->eta(),SelectedSB2Jet->phi(),massSB2Z);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());
    float MassSVFit = -1;
    float XMassSVFit = -1;
    float dRJetZSVFit = -1;
    float ptSVFit = -1;
    svfit(metRaw, met, SelectedSB2TauElTau->p4(), SelectedSB2Electron->p4(), PrunedJet, MassSVFit, XMassSVFit, dRJetZSVFit, ptSVFit, SelectedSB2Jet);
    //COLLINEAR APPROXIMATION
    TLorentzVector CATauTau; bool CA = false;
    float a = (SelectedSB2TauElTau->py()*met->begin()->px()-SelectedSB2TauElTau->px()*met->begin()->py())/
      (SelectedSB2Electron->px()*SelectedSB2TauElTau->py()-SelectedSB2Electron->py()*SelectedSB2TauElTau->px());
    float b = (SelectedSB2Electron->py()*met->begin()->px()-SelectedSB2Electron->px()*met->begin()->py())/
      (SelectedSB2TauElTau->px()*SelectedSB2Electron->py()-SelectedSB2TauElTau->py()*SelectedSB2Electron->px());
    if(((1+a)*(1+b))>0) {
      CATauTau.SetPxPyPzE((1+a)*SelectedSB2Electron->px()+(1+b)*SelectedSB2TauElTau->px(),
			  (1+a)*SelectedSB2Electron->py()+(1+b)*SelectedSB2TauElTau->py(),
			  (1+a)*SelectedSB2Electron->pz()+(1+b)*SelectedSB2TauElTau->pz(),
			  (1+a)*SelectedSB2Electron->energy()+(1+b)*SelectedSB2TauElTau->energy());
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
    math::PtEtaPhiELorentzVector dilep; dilep = SelectedSB2TauElTau->p4()+SelectedSB2Electron->p4();
    math::PtEtaPhiELorentzVector dilepmet; dilepmet = SelectedSB2TauElTau->p4()+SelectedSB2Electron->p4()+met->begin()->p4();
    math::PtEtaPhiELorentzVector dilepjet; dilepjet = SelectedSB2TauElTau->p4()+SelectedSB2Electron->p4()+PrunedJet_prov;
    math::PtEtaPhiELorentzVector dilepmetjet; dilepmetjet = SelectedSB2TauElTau->p4()+SelectedSB2Electron->p4()+met->begin()->p4()+PrunedJet_prov;
    m_jetPt=SelectedSB2Jet->pt();
    m_jetEta=SelectedSB2Jet->eta();
    m_jetMass=massSB2Z;
    m_jetSubjettiness=tau21SB2Z;
    m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB2Jet->p4());
    m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB2Jet->p4());
    m_dRJetLep=ROOT::Math::VectorUtil::DeltaR(SelectedSB2Electron->p4(),SelectedSB2Jet->p4());
    m_dRJetTau=ROOT::Math::VectorUtil::DeltaR(SelectedSB2TauElTau->p4(),SelectedSB2Jet->p4());
    m_dPhiTauMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB2TauElTau->p4());
    m_dRTauMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB2TauElTau->p4());
    m_dRTauLep=ROOT::Math::VectorUtil::DeltaR(SelectedSB2Electron->p4(),SelectedSB2TauElTau->p4());
    m_dPhiLepMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB2Electron->p4());
    m_dRLepMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB2Electron->p4());
    m_dRZZVis=ROOT::Math::VectorUtil::DeltaR(dilep,PrunedJet);
    m_dRZZEff=ROOT::Math::VectorUtil::DeltaR(dilepmet,PrunedJet);
    m_dRZZSvFit=dRJetZSVFit;
    m_dRZZCA=dRJetZCA;
    m_tauPt=SelectedSB2TauElTau->pt();
    m_tauEta=SelectedSB2TauElTau->eta();
    m_tauCharge=SelectedSB2TauElTau->charge();
    m_lepPt=SelectedSB2Electron->pt();
    m_lepEta=SelectedSB2Electron->eta();
    m_lepCharge=SelectedSB2Electron->charge();
    m_lepPFIso=ElectronPFIso(SelectedSB2Electron,rho);
    m_lepCorrPFIso=ElectronCorrPFIso(SelectedSB2Electron,rho);
    m_charge=SelectedSB2TauElTau->charge()*SelectedSB2Electron->charge();
    m_met=met->begin()->pt();
    m_metPhi=met->begin()->phi();
    m_uncorrmet=uncorrmet->begin()->pt();
    m_uncorrmetPhi=uncorrmet->begin()->phi();
    m_PtSvfit=ptSVFit;
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
    m_sideband=(int)(massSB2Z<70);
    m_NVertices=vertices->size();
    m_PUWeight=MyWeight;
    m_metPx=met->begin()->px();
    m_metPy=met->begin()->py();
    m_NeventsTOT=NeventsTOT_;
    m_xsec=xsec_;
    m_lumi=lumi_;
    m_weight=xsec_*lumi_/NeventsTOT_;
    m_genEvent = genEvent;
    TreeSB2EleTau->Fill();
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
  TreeMuoTau->Branch("PtSvfit", &m_PtSvfit, "PtSvfit/f");
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
  TreeEleTau->Branch("PtSvfit", &m_PtSvfit, "PtSvfit/f");
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

  TreeSB1MuoTau = fs->make<TTree>("TreeSB1MuoTau", "TreeSB1MuoTau");
  TreeSB1MuoTau->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB1MuoTau->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB1MuoTau->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB1MuoTau->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB1MuoTau->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB1MuoTau->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB1MuoTau->Branch("dRJetLep", &m_dRJetLep, "dRJetLep/f");
  TreeSB1MuoTau->Branch("dRJetTau", &m_dRJetTau, "dRJetTau/f");
  TreeSB1MuoTau->Branch("dPhiTauMet", &m_dPhiTauMet, "dPhiTauMet/f");
  TreeSB1MuoTau->Branch("dRTauMet", &m_dRTauMet, "dRTauMet/f");
  TreeSB1MuoTau->Branch("dRTauLep", &m_dRTauLep, "dRTauLep/f");
  TreeSB1MuoTau->Branch("dPhiLepMet", &m_dPhiLepMet, "dPhiLepMet/f");
  TreeSB1MuoTau->Branch("dRLepMet", &m_dRLepMet, "dRLepMet/f");
  TreeSB1MuoTau->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB1MuoTau->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB1MuoTau->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB1MuoTau->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB1MuoTau->Branch("tauPt", &m_tauPt, "tauPt/f");
  TreeSB1MuoTau->Branch("tauEta", &m_tauEta, "tauEta/f");
  TreeSB1MuoTau->Branch("tauCharge", &m_tauCharge, "tauCharge/f");
  TreeSB1MuoTau->Branch("lepPt", &m_lepPt, "lepPt/f");
  TreeSB1MuoTau->Branch("lepEta", &m_lepEta, "lepEta/f");
  TreeSB1MuoTau->Branch("lepCharge", &m_lepCharge, "lepCharge/f");
  TreeSB1MuoTau->Branch("lepPFIso", &m_lepPFIso, "lepPFIso/f");
  TreeSB1MuoTau->Branch("lepCorrPFIso", &m_lepCorrPFIso, "lepCorrPFIso/f");
  TreeSB1MuoTau->Branch("charge", &m_charge, "charge/f");
  TreeSB1MuoTau->Branch("met", &m_met, "met/f");
  TreeSB1MuoTau->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeSB1MuoTau->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeSB1MuoTau->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  TreeSB1MuoTau->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeSB1MuoTau->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeSB1MuoTau->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeSB1MuoTau->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeSB1MuoTau->Branch("PtSvfit", &m_PtSvfit, "PtSvfit/f");
  TreeSB1MuoTau->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeSB1MuoTau->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeSB1MuoTau->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeSB1MuoTau->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeSB1MuoTau->Branch("nbtagsL", &m_nbtagsL, "nbtagsL/i");
  TreeSB1MuoTau->Branch("nbtagsM", &m_nbtagsM, "nbtagsM/i");
  TreeSB1MuoTau->Branch("nbtagsT", &m_nbtagsT, "nbtagsT/i");
  TreeSB1MuoTau->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeSB1MuoTau->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeSB1MuoTau->Branch("trigger", &m_trigger, "trigger/i");
  TreeSB1MuoTau->Branch("sideband", &m_sideband, "sideband/i");
  TreeSB1MuoTau->Branch("NVertices", &m_NVertices, "NVertices/i");
  TreeSB1MuoTau->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeSB1MuoTau->Branch("metPx", &m_metPx, "metPx/f");
  TreeSB1MuoTau->Branch("metPy", &m_metPy, "metPy/f");
  TreeSB1MuoTau->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeSB1MuoTau->Branch("xsec", &m_xsec, "xsec/d");
  TreeSB1MuoTau->Branch("lumi", &m_lumi, "lumi/d");
  TreeSB1MuoTau->Branch("weight", &m_weight, "weight/f");
  TreeSB1MuoTau->Branch("genEvent", &m_genEvent, "genEvent/f");

  TreeSB1EleTau = fs->make<TTree>("TreeSB1EleTau", "TreeSB1EleTau");
  TreeSB1EleTau->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB1EleTau->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB1EleTau->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB1EleTau->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB1EleTau->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB1EleTau->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB1EleTau->Branch("dRJetLep", &m_dRJetLep, "dRJetLep/f");
  TreeSB1EleTau->Branch("dRJetTau", &m_dRJetTau, "dRJetTau/f");
  TreeSB1EleTau->Branch("dPhiTauMet", &m_dPhiTauMet, "dPhiTauMet/f");
  TreeSB1EleTau->Branch("dRTauMet", &m_dRTauMet, "dRTauMet/f");
  TreeSB1EleTau->Branch("dRTauLep", &m_dRTauLep, "dRTauLep/f");
  TreeSB1EleTau->Branch("dPhiLepMet", &m_dPhiLepMet, "dPhiLepMet/f");
  TreeSB1EleTau->Branch("dRLepMet", &m_dRLepMet, "dRLepMet/f");
  TreeSB1EleTau->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB1EleTau->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB1EleTau->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB1EleTau->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB1EleTau->Branch("tauPt", &m_tauPt, "tauPt/f");
  TreeSB1EleTau->Branch("tauEta", &m_tauEta, "tauEta/f");
  TreeSB1EleTau->Branch("tauCharge", &m_tauCharge, "tauCharge/f");
  TreeSB1EleTau->Branch("lepPt", &m_lepPt, "lepPt/f");
  TreeSB1EleTau->Branch("lepEta", &m_lepEta, "lepEta/f");
  TreeSB1EleTau->Branch("lepCharge", &m_lepCharge, "lepCharge/f");
  TreeSB1EleTau->Branch("lepPFIso", &m_lepPFIso, "lepPFIso/f");
  TreeSB1EleTau->Branch("lepCorrPFIso", &m_lepCorrPFIso, "lepCorrPFIso/f");
  TreeSB1EleTau->Branch("charge", &m_charge, "charge/f");
  TreeSB1EleTau->Branch("met", &m_met, "met/f");
  TreeSB1EleTau->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeSB1EleTau->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeSB1EleTau->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  TreeSB1EleTau->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeSB1EleTau->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeSB1EleTau->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeSB1EleTau->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeSB1EleTau->Branch("PtSvfit", &m_PtSvfit, "PtSvfit/f");
  TreeSB1EleTau->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeSB1EleTau->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeSB1EleTau->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeSB1EleTau->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeSB1EleTau->Branch("nbtagsL", &m_nbtagsL, "nbtagsL/i");
  TreeSB1EleTau->Branch("nbtagsM", &m_nbtagsM, "nbtagsM/i");
  TreeSB1EleTau->Branch("nbtagsT", &m_nbtagsT, "nbtagsT/i");
  TreeSB1EleTau->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeSB1EleTau->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeSB1EleTau->Branch("trigger", &m_trigger, "trigger/i");
  TreeSB1EleTau->Branch("sideband", &m_sideband, "sideband/i");
  TreeSB1EleTau->Branch("NVertices", &m_NVertices, "NVertices/i");
  TreeSB1EleTau->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeSB1EleTau->Branch("metPx", &m_metPx, "metPx/f");
  TreeSB1EleTau->Branch("metPy", &m_metPy, "metPy/f");
  TreeSB1EleTau->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeSB1EleTau->Branch("xsec", &m_xsec, "xsec/d");
  TreeSB1EleTau->Branch("lumi", &m_lumi, "lumi/d");
  TreeSB1EleTau->Branch("weight", &m_weight, "weight/f");
  TreeSB1EleTau->Branch("genEvent", &m_genEvent, "genEvent/f");

  TreeSB2MuoTau = fs->make<TTree>("TreeSB2MuoTau", "TreeSB2MuoTau");
  TreeSB2MuoTau->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB2MuoTau->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB2MuoTau->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB2MuoTau->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB2MuoTau->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB2MuoTau->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB2MuoTau->Branch("dRJetLep", &m_dRJetLep, "dRJetLep/f");
  TreeSB2MuoTau->Branch("dRJetTau", &m_dRJetTau, "dRJetTau/f");
  TreeSB2MuoTau->Branch("dPhiTauMet", &m_dPhiTauMet, "dPhiTauMet/f");
  TreeSB2MuoTau->Branch("dRTauMet", &m_dRTauMet, "dRTauMet/f");
  TreeSB2MuoTau->Branch("dRTauLep", &m_dRTauLep, "dRTauLep/f");
  TreeSB2MuoTau->Branch("dPhiLepMet", &m_dPhiLepMet, "dPhiLepMet/f");
  TreeSB2MuoTau->Branch("dRLepMet", &m_dRLepMet, "dRLepMet/f");
  TreeSB2MuoTau->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB2MuoTau->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB2MuoTau->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB2MuoTau->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB2MuoTau->Branch("tauPt", &m_tauPt, "tauPt/f");
  TreeSB2MuoTau->Branch("tauEta", &m_tauEta, "tauEta/f");
  TreeSB2MuoTau->Branch("tauCharge", &m_tauCharge, "tauCharge/f");
  TreeSB2MuoTau->Branch("lepPt", &m_lepPt, "lepPt/f");
  TreeSB2MuoTau->Branch("lepEta", &m_lepEta, "lepEta/f");
  TreeSB2MuoTau->Branch("lepCharge", &m_lepCharge, "lepCharge/f");
  TreeSB2MuoTau->Branch("lepPFIso", &m_lepPFIso, "lepPFIso/f");
  TreeSB2MuoTau->Branch("lepCorrPFIso", &m_lepCorrPFIso, "lepCorrPFIso/f");
  TreeSB2MuoTau->Branch("charge", &m_charge, "charge/f");
  TreeSB2MuoTau->Branch("met", &m_met, "met/f");
  TreeSB2MuoTau->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeSB2MuoTau->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeSB2MuoTau->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  TreeSB2MuoTau->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeSB2MuoTau->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeSB2MuoTau->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeSB2MuoTau->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeSB2MuoTau->Branch("PtSvfit", &m_PtSvfit, "PtSvfit/f");
  TreeSB2MuoTau->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeSB2MuoTau->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeSB2MuoTau->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeSB2MuoTau->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeSB2MuoTau->Branch("nbtagsL", &m_nbtagsL, "nbtagsL/i");
  TreeSB2MuoTau->Branch("nbtagsM", &m_nbtagsM, "nbtagsM/i");
  TreeSB2MuoTau->Branch("nbtagsT", &m_nbtagsT, "nbtagsT/i");
  TreeSB2MuoTau->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeSB2MuoTau->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeSB2MuoTau->Branch("trigger", &m_trigger, "trigger/i");
  TreeSB2MuoTau->Branch("sideband", &m_sideband, "sideband/i");
  TreeSB2MuoTau->Branch("NVertices", &m_NVertices, "NVertices/i");
  TreeSB2MuoTau->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeSB2MuoTau->Branch("metPx", &m_metPx, "metPx/f");
  TreeSB2MuoTau->Branch("metPy", &m_metPy, "metPy/f");
  TreeSB2MuoTau->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeSB2MuoTau->Branch("xsec", &m_xsec, "xsec/d");
  TreeSB2MuoTau->Branch("lumi", &m_lumi, "lumi/d");
  TreeSB2MuoTau->Branch("weight", &m_weight, "weight/f");
  TreeSB2MuoTau->Branch("genEvent", &m_genEvent, "genEvent/f");

  TreeSB2EleTau = fs->make<TTree>("TreeSB2EleTau", "TreeSB2EleTau");
  TreeSB2EleTau->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB2EleTau->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB2EleTau->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB2EleTau->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB2EleTau->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB2EleTau->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB2EleTau->Branch("dRJetLep", &m_dRJetLep, "dRJetLep/f");
  TreeSB2EleTau->Branch("dRJetTau", &m_dRJetTau, "dRJetTau/f");
  TreeSB2EleTau->Branch("dPhiTauMet", &m_dPhiTauMet, "dPhiTauMet/f");
  TreeSB2EleTau->Branch("dRTauMet", &m_dRTauMet, "dRTauMet/f");
  TreeSB2EleTau->Branch("dRTauLep", &m_dRTauLep, "dRTauLep/f");
  TreeSB2EleTau->Branch("dPhiLepMet", &m_dPhiLepMet, "dPhiLepMet/f");
  TreeSB2EleTau->Branch("dRLepMet", &m_dRLepMet, "dRLepMet/f");
  TreeSB2EleTau->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB2EleTau->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB2EleTau->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB2EleTau->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB2EleTau->Branch("tauPt", &m_tauPt, "tauPt/f");
  TreeSB2EleTau->Branch("tauEta", &m_tauEta, "tauEta/f");
  TreeSB2EleTau->Branch("tauCharge", &m_tauCharge, "tauCharge/f");
  TreeSB2EleTau->Branch("lepPt", &m_lepPt, "lepPt/f");
  TreeSB2EleTau->Branch("lepEta", &m_lepEta, "lepEta/f");
  TreeSB2EleTau->Branch("lepCharge", &m_lepCharge, "lepCharge/f");
  TreeSB2EleTau->Branch("lepPFIso", &m_lepPFIso, "lepPFIso/f");
  TreeSB2EleTau->Branch("lepCorrPFIso", &m_lepCorrPFIso, "lepCorrPFIso/f");
  TreeSB2EleTau->Branch("charge", &m_charge, "charge/f");
  TreeSB2EleTau->Branch("met", &m_met, "met/f");
  TreeSB2EleTau->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeSB2EleTau->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeSB2EleTau->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  TreeSB2EleTau->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeSB2EleTau->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeSB2EleTau->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeSB2EleTau->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeSB2EleTau->Branch("PtSvfit", &m_PtSvfit, "PtSvfit/f");
  TreeSB2EleTau->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeSB2EleTau->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeSB2EleTau->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeSB2EleTau->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeSB2EleTau->Branch("nbtagsL", &m_nbtagsL, "nbtagsL/i");
  TreeSB2EleTau->Branch("nbtagsM", &m_nbtagsM, "nbtagsM/i");
  TreeSB2EleTau->Branch("nbtagsT", &m_nbtagsT, "nbtagsT/i");
  TreeSB2EleTau->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeSB2EleTau->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeSB2EleTau->Branch("trigger", &m_trigger, "trigger/i");
  TreeSB2EleTau->Branch("sideband", &m_sideband, "sideband/i");
  TreeSB2EleTau->Branch("NVertices", &m_NVertices, "NVertices/i");
  TreeSB2EleTau->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeSB2EleTau->Branch("metPx", &m_metPx, "metPx/f");
  TreeSB2EleTau->Branch("metPy", &m_metPy, "metPy/f");
  TreeSB2EleTau->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeSB2EleTau->Branch("xsec", &m_xsec, "xsec/d");
  TreeSB2EleTau->Branch("lumi", &m_lumi, "lumi/d");
  TreeSB2EleTau->Branch("weight", &m_weight, "weight/f");
  TreeSB2EleTau->Branch("genEvent", &m_genEvent, "genEvent/f");
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

void SemiLeptonicAnalyzer::SelectJet(edm::Handle<pat::JetCollection> CA8JetswithQjets,
				     edm::Handle<pat::JetCollection> CA8JetsPruned,
				     bool & foundJet,
				     pat::JetCollection::const_iterator & SelectedJet,
				     float & massZ, float & tau21Z, float & ptZ,
				     float massMin, float massMax){

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
    if(jet->pt()<400) continue;
    if(abs(jet->eta())>2.4) continue;
    if(!(mass>massMin && mass<massMax))  continue;
    if(jet->userFloat("tau2")/jet->userFloat("tau1")>0.75) continue;
    foundJet=true;
    if(jet->pt()>ptZ){
      massZ=mass;
      ptZ=jet->pt();
      tau21Z=jet->userFloat("tau2")/jet->userFloat("tau1");
      SelectedJet=jet;
    }
  }
}

void SemiLeptonicAnalyzer::SelectTau(edm::Handle<pat::TauCollection> tauHandle,
					  pat::JetCollection::const_iterator SelectedJet,
					  bool & foundTau,
					  pat::TauCollection::const_iterator & SelectedTau,
					  float & ptTau, bool foundJet){
  
  for (pat::TauCollection::const_iterator patTau = tauHandle->begin(); patTau != tauHandle->end(); ++patTau ) {
    if(patTau->pt()<20) continue;
    if(abs(patTau->eta())>2.4) continue;
    if(patTau->tauID("decayModeFinding")<0.5) continue;
    if(patTau->tauID("againstMuonLoose")<0.5) continue;
    if(patTau->tauID("againstElectronLoose")<0.5) continue;
    if(patTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr")<0.5) continue;
    if(foundJet){
      if(ROOT::Math::VectorUtil::DeltaR(patTau->p4(),SelectedJet->p4())>0.8){
	foundTau=true;
	if(patTau->pt()>ptTau){
	  SelectedTau=patTau;
	  ptTau=patTau->pt();
	}
      }
    }
  }
}

void SemiLeptonicAnalyzer::SelectMuon(edm::Handle<pat::MuonCollection> muoH,
				      pat::JetCollection::const_iterator SelectedJet, bool & foundMuon,
				      pat::MuonCollection::const_iterator & SelectedMuon,
				      float & ptMuon, bool foundJet, reco::Vertex primaryVertex){
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
    if(foundJet) {
      if(ROOT::Math::VectorUtil::DeltaR(muon->p4(),SelectedJet->p4())>0.8){
	foundMuon=true;
	if(cktTrack->pt()>ptMuon){
	  SelectedMuon=muon;
	  ptMuon=cktTrack->pt();
	}
      }
    }
  }
}

void SemiLeptonicAnalyzer::SelectElectron(edm::Handle<pat::ElectronCollection> eleH,
					  pat::JetCollection::const_iterator SelectedJet, bool & foundElectron,
					  pat::ElectronCollection::const_iterator & SelectedElectron,
					  float & ptElectron, bool foundJet, reco::Vertex primaryVertex){
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
    if(foundJet){
      if(ROOT::Math::VectorUtil::DeltaR(electron->p4(),SelectedJet->p4())>0.8){
	foundElectron=true;
	if(electron->pt()>ptElectron){
	  SelectedElectron=electron;
	  ptElectron=electron->pt();
	}
      }
    }
  }
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


float SemiLeptonicAnalyzer::ElectronPFIso(pat::ElectronCollection::const_iterator electron, float rho){
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

float SemiLeptonicAnalyzer::MuonCorrPFIso(pat::MuonCollection::const_iterator muon, bool highpt){
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


float SemiLeptonicAnalyzer::ElectronCorrPFIso(pat::ElectronCollection::const_iterator electron, float rho){
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

void SemiLeptonicAnalyzer::Efficiency(float & genEvent, bool isData, Handle<vector<reco::GenParticle> > genParts){
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
}

void SemiLeptonicAnalyzer::svfit(edm::Handle<pat::METCollection> metRaw, edm::Handle<pat::METCollection> met, LorentzVector SelectedTau, 
				 LorentzVector SelectedMuon, TLorentzVector PrunedJet, float & MassSVFit, float & XMassSVFit, float & dRJetZSVFit, float & ptSVFit, 
				 pat::JetCollection::const_iterator SelectedJet){
  TMatrixD covMET(2, 2); // PFMET significance matrix
  covMET[0][0] = (metRaw->front() ).getSignificanceMatrix()(0,0);
  covMET[1][0] = (metRaw->front() ).getSignificanceMatrix()(1,0);
  covMET[0][1] = (metRaw->front() ).getSignificanceMatrix()(0,1);
  covMET[1][1] = (metRaw->front() ).getSignificanceMatrix()(1,1);
  std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay, SelectedTau));
  measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, SelectedMuon));
  NSVfitStandaloneAlgorithm algo(measuredTauLeptons, met->begin()->momentum(), covMET, 0);
  algo.addLogM(false);
  algo.integrateMarkovChain();
  if(algo.pt()>0){
    TLorentzVector SVFitTauTau; SVFitTauTau.SetPtEtaPhiM(algo.pt(), algo.eta(), algo.phi(), algo.getMass());
    dRJetZSVFit=ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,SelectedJet->p4());
    XMassSVFit=(SVFitTauTau+PrunedJet).M();
    MassSVFit=algo.getMass();
    ptSVFit=algo.pt();
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(SemiLeptonicAnalyzer);
