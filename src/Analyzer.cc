// -*- C++ -*-
//
// Package:    Analyzer
// Class:      Analyzer
// 
/**\class Analyzer Analyzer.cc ExoDiBosonResonances/EDBRTauAnalyzer/src/Analyzer.cc

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

class Analyzer : public edm::EDAnalyzer {
public:
  explicit Analyzer(const edm::ParameterSet&);
  ~Analyzer();
  
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
                 pat::JetCollection::const_iterator & SelectedJet, float & prunedMass, float & tau21Z, float & ptZ, float massMin, float massMax);
  void SelectTau(edm::Handle<pat::TauCollection> tauHandle, pat::JetCollection::const_iterator SelectedJet, bool & foundTau, 
		 pat::TauCollection::const_iterator & SelectedTau, float & ptTau, bool foundJet);
  void SelectMuon(edm::Handle<pat::MuonCollection> muoH, pat::JetCollection::const_iterator SelectedJet, bool & foundMuon,
		  pat::MuonCollection::const_iterator & SelectedMuon, float & ptMuon, bool foundJet, reco::Vertex primaryVertex, bool fully);
  void SelectElectron(edm::Handle<pat::ElectronCollection> eleH, pat::JetCollection::const_iterator SelectedJet, bool & foundElectron,
		      pat::ElectronCollection::const_iterator & SelectedElectron, float & ptElectron, bool foundJet, reco::Vertex primaryVertex, 
		      float lep2Pt, float rho, bool fully);
  void SelectMuonMM(pat::MuonCollection::const_iterator & SelectedMuon1, pat::MuonCollection::const_iterator & SelectedMuon2, bool & foundMuonMM, 
		    std::vector<pat::MuonCollection::const_iterator> SelectedTrackerMuo, pat::JetCollection::const_iterator SelectedJet, bool foundJet, 
		    reco::Vertex primaryVertex);
  void SelectTrackerMuon(edm::Handle<pat::MuonCollection> muoH, std::vector<pat::MuonCollection::const_iterator> & SelectedMuo, reco::Vertex primaryVertex);
  void Efficiency(float & genEvent, bool isData, edm::Handle<std::vector<reco::GenParticle> > genParts);
  void svfit(edm::Handle<pat::METCollection> metRaw, edm::Handle<pat::METCollection> met, LorentzVector SelectedTau, LorentzVector SelectedMuon, TLorentzVector PrunedJet,
	     float & MassSVFit, float & XMassSVFit, float & dRJetZSVFit, float & ptSVFit, pat::JetCollection::const_iterator SelectedJet, int category);
  
  float MuonPFIso(pat::MuonCollection::const_iterator muon, bool highpt);
  float ElectronPFIso(pat::ElectronCollection::const_iterator electron, float rho);
  float MuonCorrPFIso(pat::MuonCollection::const_iterator muon, bool highpt);
  float ElectronCorrPFIso(pat::ElectronCollection::const_iterator electron, float rho);
  float MuonDETIso(pat::MuonCollection::const_iterator SelectedMuo, std::vector<pat::MuonCollection::const_iterator> SelectedHighptMuo);
  void BtagVeto(edm::Handle<pat::JetCollection> ak5jetCands, int & nbtagsL, int & nbtagsM, int & nbtagsT, pat::JetCollection::const_iterator SelectedJet);
  void FillTree(int category, TTree *Tree, pat::JetCollection::const_iterator SelectedJet, pat::TauCollection::const_iterator SelectedTau,
		pat::MuonCollection::const_iterator SelectedMuo, pat::MuonCollection::const_iterator SelectedMuo1, pat::MuonCollection::const_iterator SelectedMuo2,
		pat::ElectronCollection::const_iterator SelectedEle, pat::ElectronCollection::const_iterator SelectedEle1, 
		pat::ElectronCollection::const_iterator SelectedEle2, edm::Handle<reco::VertexCollection> vertices,
		edm::Handle<pat::METCollection> metRaw, edm::Handle<pat::METCollection> met, edm::Handle<pat::METCollection> uncorrmet,
		float prunedMass, int nbtagsL, int nbtagsM, int nbtagsT, bool isFired_HLT, bool isFired_HLT_PFJet320, bool isFired_HLT_HT650,
		double MyWeight, float genEvent, float rho, int EleMuo, int MuoMuo, int EleEle, int MuoTau, int EleTau);
  
  TH1D* Nevents;

  TTree *TreeSignalEff;
  float m_genEvent;

  TTree *TreeEleMuo;
  TTree *TreeMuoMuo;
  TTree *TreeEleEle;
  TTree *TreeMuoTau;
  TTree *TreeEleTau;
  TTree *TreeSB1EleMuo;
  TTree *TreeSB1MuoMuo;
  TTree *TreeSB1EleEle;
  TTree *TreeSB1MuoTau;
  TTree *TreeSB1EleTau;
  TTree *TreeSB2EleMuo;
  TTree *TreeSB2MuoMuo;
  TTree *TreeSB2EleEle;
  TTree *TreeSB2MuoTau;
  TTree *TreeSB2EleTau;
  float m_jetPt;
  float m_jetEta;
  float m_jetMass;
  float m_jetSubjettiness;
  float m_dPhiJetMet;
  float m_dRJetMet;
  float m_dRJetLep2;
  float m_dRJetLep1;
  float m_dPhiLep1Met;
  float m_dRLep1Met;
  float m_dRLep1Lep2;
  float m_dPhiLep2Met;
  float m_dRLep2Met;
  float m_dRZZVis;
  float m_dRZZEff;
  float m_dRZZSvFit;
  float m_dRZZCA;
  float m_lep1Pt;
  float m_lep1Eta;
  float m_lep1Charge;
  float m_lep1PFIso;
  float m_lep1CorrPFIso;
  float m_lep2Pt;
  float m_lep2Eta;
  float m_lep2Charge;
  float m_lep2PFIso;
  float m_lep2CorrPFIso;
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
  int m_EleMuo;
  int m_MuoMuo;
  int m_EleEle;
  int m_MuoTau;
  int m_EleTau;

  edm::LumiReWeighting LumiWeights_;
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
Analyzer::Analyzer(const edm::ParameterSet& iConfig)

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


Analyzer::~Analyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  edm::Handle<pat::ElectronCollection> eleETH;
  iEvent.getByLabel(electronETColl_, eleETH);

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
  float prunedMass=-9999; float ptZ=-999; bool foundJet=false; float tau21Z=-9999;
  SelectJet(CA8JetswithQjets, CA8JetsPruned, foundJet, SelectedJet, prunedMass, tau21Z, ptZ, 70, 110);

  //JET SELECTION - SB1
  pat::JetCollection::const_iterator SelectedSB1Jet;
  float prunedMassSB1=-9999; float ptSB1Z=-999; bool foundSB1Jet=false; float tau21SB1Z=-9999;
  SelectJet(CA8JetswithQjets, CA8JetsPruned, foundSB1Jet, SelectedSB1Jet, prunedMassSB1, tau21SB1Z, ptSB1Z, 20, 70);

  //JET SELECTION - SB2
  pat::JetCollection::const_iterator SelectedSB2Jet;
  float prunedMassSB2=-9999; float ptSB2Z=-999; bool foundSB2Jet=false; float tau21SB2Z=-9999;
  SelectJet(CA8JetswithQjets, CA8JetsPruned, foundSB2Jet, SelectedSB2Jet, prunedMassSB2, tau21SB2Z, ptSB2Z, 110, 99999);

  //TAU SELECTION - MUOTAU - SR
  float ptTauMT=-99; bool foundTauMT=false;
  pat::TauCollection::const_iterator SelectedTauMT;
  SelectTau(tauMuTauHandle, SelectedJet, foundTauMT, SelectedTauMT, ptTauMT, foundJet);

  //TAU SELECTION - MUOTAU - SB1
  float ptSB1TauMT=-99; bool foundSB1TauMT=false;
  pat::TauCollection::const_iterator SelectedSB1TauMT;
  SelectTau(tauMuTauHandle, SelectedSB1Jet, foundSB1TauMT, SelectedSB1TauMT, ptSB1TauMT, foundSB1Jet);

  //TAU SELECTION - MUOTAU - SB2
  float ptSB2TauMT=-99; bool foundSB2TauMT=false;
  pat::TauCollection::const_iterator SelectedSB2TauMT;
  SelectTau(tauMuTauHandle, SelectedSB2Jet, foundSB2TauMT, SelectedSB2TauMT, ptSB2TauMT, foundSB2Jet); 

  //TAU SELECTION - ELETAU - SR
  float ptTauET=-99; bool foundTauET=false;
  pat::TauCollection::const_iterator SelectedTauET;
  SelectTau(tauElTauHandle, SelectedJet, foundTauET, SelectedTauET, ptTauET, foundJet);

  //TAU SELECTION - ELETAU - SB1
  float ptSB1TauET=-99; bool foundSB1TauET=false;
  pat::TauCollection::const_iterator SelectedSB1TauET;
  SelectTau(tauElTauHandle, SelectedSB1Jet, foundSB1TauET, SelectedSB1TauET, ptSB1TauET, foundSB1Jet);

  //TAU SELECTION - ELETAU - SB2
  float ptSB2TauET=-99; bool foundSB2TauET=false;
  pat::TauCollection::const_iterator SelectedSB2TauET;
  SelectTau(tauElTauHandle, SelectedSB2Jet, foundSB2TauET, SelectedSB2TauET, ptSB2TauET, foundSB2Jet); 

  //MUON SELECTION - MUOTAU - SR
  float ptMuonMT=-99; bool foundMuonMT=false;
  pat::MuonCollection::const_iterator SelectedMuonMT;
  SelectMuon(muoH, SelectedJet, foundMuonMT, SelectedMuonMT, ptMuonMT, foundJet, primaryVertex, false);

  //MUON SELECTION - MUOTAU - SB1
  float ptSB1MuonMT=-99; bool foundSB1MuonMT=false;
  pat::MuonCollection::const_iterator SelectedSB1MuonMT;
  SelectMuon(muoH, SelectedSB1Jet, foundSB1MuonMT, SelectedSB1MuonMT, ptSB1MuonMT, foundSB1Jet, primaryVertex, false);

  //MUON SELECTION - MUOTAU - SB2
  float ptSB2MuonMT=-99; bool foundSB2MuonMT=false;
  pat::MuonCollection::const_iterator SelectedSB2MuonMT;
  SelectMuon(muoH, SelectedSB2Jet, foundSB2MuonMT, SelectedSB2MuonMT, ptSB2MuonMT, foundSB2Jet, primaryVertex, false);

  //MUON SELECTION - ELEMUO - SR
  float ptMuonEM=-99; bool foundMuonEM=false;
  pat::MuonCollection::const_iterator SelectedMuonEM;
  SelectMuon(muoH, SelectedJet, foundMuonEM, SelectedMuonEM, ptMuonEM, foundJet, primaryVertex, false);

  //MUON SELECTION - ELEMUO - SB1
  float ptSB1MuonEM=-99; bool foundSB1MuonEM=false;
  pat::MuonCollection::const_iterator SelectedSB1MuonEM;
  SelectMuon(muoH, SelectedSB1Jet, foundSB1MuonEM, SelectedSB1MuonEM, ptSB1MuonEM, foundSB1Jet, primaryVertex, false);

  //MUON SELECTION - ELEMUO - SB2
  float ptSB2MuonEM=-99; bool foundSB2MuonEM=false;
  pat::MuonCollection::const_iterator SelectedSB2MuonEM;
  SelectMuon(muoH, SelectedSB2Jet, foundSB2MuonEM, SelectedSB2MuonEM, ptSB2MuonEM, foundSB2Jet, primaryVertex, false);

  //ELECTRON SELECTION - ELETAU - SR
  float ptElectronET=-99; bool foundElectronET=false;
  pat::ElectronCollection::const_iterator SelectedElectronET;
  SelectElectron(eleETH, SelectedJet, foundElectronET, SelectedElectronET, ptElectronET, foundJet, primaryVertex, -1, rho, false);

  //ELECTRON SELECTION - ELETAU - SB1
  float ptSB1ElectronET=-99; bool foundSB1ElectronET=false;
  pat::ElectronCollection::const_iterator SelectedSB1ElectronET;
  SelectElectron(eleETH, SelectedSB1Jet, foundSB1ElectronET, SelectedSB1ElectronET, ptSB1ElectronET, foundSB1Jet, primaryVertex, -1, rho, false);

  //ELECTRON SELECTION - ELETAU - SB2
  float ptSB2ElectronET=-99; bool foundSB2ElectronET=false;
  pat::ElectronCollection::const_iterator SelectedSB2ElectronET;
  SelectElectron(eleETH, SelectedSB2Jet, foundSB2ElectronET, SelectedSB2ElectronET, ptSB2ElectronET, foundSB2Jet, primaryVertex, -1, rho, false);

  //ELECTRON SELECTION - ELEELE1 - SR
  float ptElectron1EE=-99; bool foundElectron1EE=false;
  pat::ElectronCollection::const_iterator SelectedElectron1EE;
  SelectElectron(eleH, SelectedJet, foundElectron1EE, SelectedElectron1EE, ptElectron1EE, foundJet, primaryVertex, -1, rho, true);

  //ELECTRON SELECTION - ELEELE1 - SB1
  float ptSB1Electron1EE=-99; bool foundSB1Electron1EE=false;
  pat::ElectronCollection::const_iterator SelectedSB1Electron1EE;
  SelectElectron(eleH, SelectedSB1Jet, foundSB1Electron1EE, SelectedSB1Electron1EE, ptSB1Electron1EE, foundSB1Jet, primaryVertex, -1, rho, true);

  //ELECTRON SELECTION - ELEELE1 - SB2
  float ptSB2Electron1EE=-99; bool foundSB2Electron1EE=false;
  pat::ElectronCollection::const_iterator SelectedSB2Electron1EE;
  SelectElectron(eleH, SelectedSB2Jet, foundSB2Electron1EE, SelectedSB2Electron1EE, ptSB2Electron1EE, foundSB2Jet, primaryVertex, -1, rho, true);

  //ELECTRON SELECTION - ELEELE2 - SR
  float ptElectron2EE=-99; bool foundElectron2EE=false;
  pat::ElectronCollection::const_iterator SelectedElectron2EE;
  SelectElectron(eleH, SelectedJet, foundElectron2EE, SelectedElectron2EE, ptElectron2EE, foundJet, primaryVertex, ptElectron1EE, rho, true);

  //ELECTRON SELECTION - ELEELE2 - SB1
  float ptSB1Electron2EE=-99; bool foundSB1Electron2EE=false;
  pat::ElectronCollection::const_iterator SelectedSB1Electron2EE;
  SelectElectron(eleH, SelectedSB1Jet, foundSB1Electron2EE, SelectedSB1Electron2EE, ptSB1Electron2EE, foundSB1Jet, primaryVertex, ptSB1Electron1EE, rho, true);

  //ELECTRON SELECTION - ELEELE2 - SB2
  float ptSB2Electron2EE=-99; bool foundSB2Electron2EE=false;
  pat::ElectronCollection::const_iterator SelectedSB2Electron2EE;
  SelectElectron(eleH, SelectedSB2Jet, foundSB2Electron2EE, SelectedSB2Electron2EE, ptSB2Electron2EE, foundSB2Jet, primaryVertex, ptSB2Electron1EE, rho, true);

  //ELECTRON SELECTION - ELEMUO - SR
  float ptElectronEM=-99; bool foundElectronEM=false;
  pat::ElectronCollection::const_iterator SelectedElectronEM;
  SelectElectron(eleH, SelectedJet, foundElectronEM, SelectedElectronEM, ptElectronEM, foundJet, primaryVertex, -1, rho, true);

  //ELECTRON SELECTION - ELEMUO - SB1
  float ptSB1ElectronEM=-99; bool foundSB1ElectronEM=false;
  pat::ElectronCollection::const_iterator SelectedSB1ElectronEM;
  SelectElectron(eleH, SelectedSB1Jet, foundSB1ElectronEM, SelectedSB1ElectronEM, ptSB1ElectronEM, foundSB1Jet, primaryVertex, -1, rho, true);

  //ELECTRON SELECTION - ELEMUO - SB2
  float ptSB2ElectronEM=-99; bool foundSB2ElectronEM=false;
  pat::ElectronCollection::const_iterator SelectedSB2ElectronEM;
  SelectElectron(eleH, SelectedSB2Jet, foundSB2ElectronEM, SelectedSB2ElectronEM, ptSB2ElectronEM, foundSB2Jet, primaryVertex, -1, rho, true);

  //MUON SELECTION - MUOMUO - SR
  vector<pat::MuonCollection::const_iterator> SelectedTrackerMuo;
  SelectTrackerMuon(muoH, SelectedTrackerMuo, primaryVertex);
  bool foundMuonMM=false; 
  pat::MuonCollection::const_iterator SelectedMuon1MM;
  pat::MuonCollection::const_iterator SelectedMuon2MM;
  SelectMuonMM(SelectedMuon1MM, SelectedMuon2MM, foundMuonMM, SelectedTrackerMuo, SelectedJet, foundJet, primaryVertex);

  //MUON SELECTION - MUOMUO - SB1
  bool foundSB1MuonMM=false;
  pat::MuonCollection::const_iterator SelectedSB1Muon1MM;
  pat::MuonCollection::const_iterator SelectedSB1Muon2MM;
  SelectMuonMM(SelectedSB1Muon1MM, SelectedSB1Muon2MM, foundSB1MuonMM, SelectedTrackerMuo, SelectedSB1Jet, foundSB1Jet, primaryVertex);

  //MUON SELECTION - MUOMUO - SB2
  bool foundSB2MuonMM=false;
  pat::MuonCollection::const_iterator SelectedSB2Muon1MM;
  pat::MuonCollection::const_iterator SelectedSB2Muon2MM;
  SelectMuonMM(SelectedSB2Muon1MM, SelectedSB2Muon2MM, foundSB2MuonMM, SelectedTrackerMuo, SelectedSB2Jet, foundSB2Jet, primaryVertex);

  //BTAG VETO - SR
  int nbtagsL=0; int nbtagsM=0; int nbtagsT=0;
  if(foundJet) BtagVeto(ak5jetCands, nbtagsL, nbtagsM, nbtagsT, SelectedJet);

  //BTAG VETO - SB1
  int nSB1btagsL=0; int nSB1btagsM=0; int nSB1btagsT=0;
  if(foundSB1Jet) BtagVeto(ak5jetCands, nSB1btagsL, nSB1btagsM, nSB1btagsT, SelectedSB1Jet);

  //BTAG VETO - SB2
  int nSB2btagsL=0; int nSB2btagsM=0; int nSB2btagsT=0;
  if(foundSB2Jet) BtagVeto(ak5jetCands, nSB2btagsL, nSB2btagsM, nSB2btagsT, SelectedSB2Jet);

  pat::TauCollection::const_iterator      SelectedTauFake;
  pat::MuonCollection::const_iterator     SelectedMuonFake;
  pat::MuonCollection::const_iterator     SelectedMuon1Fake;
  pat::MuonCollection::const_iterator     SelectedMuon2Fake;
  pat::ElectronCollection::const_iterator SelectedElectronFake;
  pat::ElectronCollection::const_iterator SelectedElectron1Fake;
  pat::ElectronCollection::const_iterator SelectedElectron2Fake;

  int EleMuo = 0; int MuoMuo = 0; int EleEle = 0; int MuoTau = 0; int EleTau = 0;
  if(foundJet && foundMuonEM && foundElectronEM) EleMuo=1;
  if(foundJet && foundMuonMM) MuoMuo=1;
  if(foundJet && foundElectron1EE && foundElectron2EE) EleEle=1;
  if(foundJet && foundMuonMT && foundTauMT) MuoTau=1;
  if(foundJet && foundElectronET && foundTauET) EleTau=1;

  //ELE-MUO ANALYSIS - SR
  if(foundJet && foundMuonEM && foundElectronEM){
    cout<<"EleMuo "<<iEvent.id().event()<<"; muon pt "<<SelectedMuonEM->pt()<<"; ele pt "<<SelectedElectronEM->pt()<<
      "; jet pt "<<SelectedJet->pt()<<"; jet mass "<<prunedMass<<"; met "<<met->begin()->pt()<<endl;
    FillTree(0,TreeEleMuo,SelectedJet,SelectedTauFake,SelectedMuonEM,SelectedMuon1Fake,SelectedMuon2Fake,SelectedElectronEM,SelectedElectron1Fake,
	     SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMass, nbtagsL, nbtagsM, nbtagsT, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuo, MuoMuo, EleEle, MuoTau, EleTau);
  }
  //MUO-MUO ANALYSIS - SR
  if(foundJet && foundMuonMM){
    cout<<"MuoMuo "<<iEvent.id().event()<<"; muon1 pt "<<SelectedMuon1MM->pt()<<"; muon2 pt "<<SelectedMuon2MM->pt()<<
      "; jet pt "<<SelectedJet->pt()<<"; jet mass "<<prunedMass<<"; met "<<met->begin()->pt()<<endl;
    FillTree(1,TreeMuoMuo,SelectedJet,SelectedTauFake,SelectedMuonFake,SelectedMuon1MM,SelectedMuon2MM,SelectedElectronFake,SelectedElectron1Fake,
	     SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMass, nbtagsL, nbtagsM, nbtagsT, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuo, MuoMuo, EleEle, MuoTau, EleTau);
  }
  //ELE-ELE ANALYSIS - SR
  if(foundJet && foundElectron1EE && foundElectron2EE){
    cout<<"EleEle "<<iEvent.id().event()<<"; ele1 pt "<<SelectedElectron1EE->pt()<<"; ele2 pt "<<SelectedElectron2EE->pt()<<
      "; jet pt "<<SelectedJet->pt()<<"; jet mass "<<prunedMass<<"; met "<<met->begin()->pt()<<endl;
    FillTree(2,TreeEleEle,SelectedJet,SelectedTauFake,SelectedMuonFake,SelectedMuon1Fake,SelectedMuon2Fake,SelectedElectronFake,SelectedElectron1EE,
	     SelectedElectron2EE, vertices, metRaw, met, uncorrmet, prunedMass, nbtagsL, nbtagsM, nbtagsT, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuo, MuoMuo, EleEle, MuoTau, EleTau);
  }
  //MUO-TAU ANALYSIS - SR
  if(foundJet && foundMuonMT && foundTauMT){
    cout<<"MuoTau "<<iEvent.id().event()<<"; tau pt "<<SelectedTauMT->pt()<<"; muon pt "<<SelectedMuonMT->pt()<<
      "; jet pt "<<SelectedJet->pt()<<"; jet mass "<<prunedMass<<"; met "<<met->begin()->pt()<<endl;
    FillTree(3,TreeMuoTau,SelectedJet,SelectedTauMT,SelectedMuonMT,SelectedMuon1Fake,SelectedMuon2Fake,SelectedElectronFake,SelectedElectron1Fake,
	     SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMass, nbtagsL, nbtagsM, nbtagsT, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuo, MuoMuo, EleEle, MuoTau, EleTau);
  }
  //ELE-TAU ANALYSIS - SR
  if(foundJet && foundElectronET && foundTauET){
    cout<<"EleTau "<<iEvent.id().event()<<"; tau pt "<<SelectedTauET->pt()<<"; electron pt "<<SelectedElectronET->pt()<<
      "; jet pt "<<SelectedJet->pt()<<"; jet mass "<<prunedMass<<"; met "<<met->begin()->pt()<<endl;
    FillTree(4,TreeEleTau,SelectedJet,SelectedTauET,SelectedMuonFake,SelectedMuon1Fake,SelectedMuon2Fake,SelectedElectronET,SelectedElectron1Fake,
	     SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMass, nbtagsL, nbtagsM, nbtagsT, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuo, MuoMuo, EleEle, MuoTau, EleTau);
  }


  int EleMuoSB1 = 0; int MuoMuoSB1 = 0; int EleEleSB1 = 0; int MuoTauSB1 = 0; int EleTauSB1 = 0;
  if(foundSB1Jet && foundSB1MuonEM && foundSB1ElectronEM) EleMuo=1;
  if(foundSB1Jet && foundSB1MuonMM) MuoMuoSB1=1;
  if(foundSB1Jet && foundSB1Electron1EE && foundSB1Electron2EE) EleEleSB1=1;
  if(foundSB1Jet && foundSB1MuonMT && foundSB1TauMT) MuoTauSB1=1;
  if(foundSB1Jet && foundSB1ElectronET && foundSB1TauET) EleTauSB1=1;

  //ELE-MUO ANALYSIS - SB1
  if(foundSB1Jet && foundSB1MuonEM && foundSB1ElectronEM){
    FillTree(0,TreeSB1EleMuo,SelectedSB1Jet,SelectedTauFake,SelectedSB1MuonEM,SelectedMuon1Fake,SelectedMuon2Fake,SelectedSB1ElectronEM,SelectedElectron1Fake,
	     SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMassSB1, nbtagsL, nbtagsM, nbtagsT, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB1, MuoMuoSB1, EleEleSB1, MuoTauSB1, EleTauSB1);
  }
  //MUO-MUO ANALYSIS - SB1
  if(foundSB1Jet && foundSB1MuonMM){
    FillTree(1,TreeSB1MuoMuo,SelectedSB1Jet,SelectedTauFake,SelectedMuonFake,SelectedSB1Muon1MM,SelectedSB1Muon2MM,SelectedElectronFake,SelectedElectron1Fake,
	     SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMassSB1, nbtagsL, nbtagsM, nbtagsT, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB1, MuoMuoSB1, EleEleSB1, MuoTauSB1, EleTauSB1);
  }
  //ELE-ELE ANALYSIS - SB1
  if(foundSB1Jet && foundSB1Electron1EE && foundSB1Electron2EE){
    FillTree(2,TreeSB1EleEle,SelectedSB1Jet,SelectedTauFake,SelectedMuonFake,SelectedMuon1Fake,SelectedMuon2Fake,SelectedElectronFake,SelectedSB1Electron1EE,
	     SelectedSB1Electron2EE, vertices, metRaw, met, uncorrmet, prunedMassSB1, nbtagsL, nbtagsM, nbtagsT, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB1, MuoMuoSB1, EleEleSB1, MuoTauSB1, EleTauSB1);
  }
  //MUO-TAU ANALYSIS - SB1
  if(foundSB1Jet && foundSB1MuonMT && foundSB1TauMT){
    FillTree(3,TreeSB1MuoTau,SelectedSB1Jet,SelectedSB1TauMT,SelectedSB1MuonMT,SelectedMuon1Fake,SelectedMuon2Fake,SelectedElectronFake,SelectedElectron1Fake,
	     SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMassSB1, nbtagsL, nbtagsM, nbtagsT, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB1, MuoMuoSB1, EleEleSB1, MuoTauSB1, EleTauSB1);
  }
  //ELE-TAU ANALYSIS - SB1
  if(foundSB1Jet && foundSB1ElectronET && foundSB1TauET){
    FillTree(4,TreeSB1EleTau,SelectedSB1Jet,SelectedSB1TauET,SelectedMuonFake,SelectedMuon1Fake,SelectedMuon2Fake,SelectedSB1ElectronET,SelectedElectron1Fake,
	     SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMassSB1, nbtagsL, nbtagsM, nbtagsT, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB1, MuoMuoSB1, EleEleSB1, MuoTauSB1, EleTauSB1);
  }


  int EleMuoSB2 = 0; int MuoMuoSB2 = 0; int EleEleSB2 = 0; int MuoTauSB2 = 0; int EleTauSB2 = 0;
  if(foundSB2Jet && foundSB2MuonEM && foundSB2ElectronEM) EleMuo=1;
  if(foundSB2Jet && foundSB2MuonMM) MuoMuoSB2=1;
  if(foundSB2Jet && foundSB2Electron1EE && foundSB2Electron2EE) EleEleSB2=1;
  if(foundSB2Jet && foundSB2MuonMT && foundSB2TauMT) MuoTauSB2=1;
  if(foundSB2Jet && foundSB2ElectronET && foundSB2TauET) EleTauSB2=1;

  //ELE-MUO ANALYSIS - SB2
  if(foundSB2Jet && foundSB2MuonEM && foundSB2ElectronEM){
    FillTree(0,TreeSB2EleMuo,SelectedSB2Jet,SelectedTauFake,SelectedSB2MuonEM,SelectedMuon1Fake,SelectedMuon2Fake,SelectedSB2ElectronEM,SelectedElectron1Fake,
	     SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMassSB2, nbtagsL, nbtagsM, nbtagsT, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB2, MuoMuoSB2, EleEleSB2, MuoTauSB2, EleTauSB2);
  }
  //MUO-MUO ANALYSIS - SB2
  if(foundSB2Jet && foundSB2MuonMM){
    FillTree(1,TreeSB2MuoMuo,SelectedSB2Jet,SelectedTauFake,SelectedMuonFake,SelectedSB2Muon1MM,SelectedSB2Muon2MM,SelectedElectronFake,SelectedElectron1Fake,
	     SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMassSB2, nbtagsL, nbtagsM, nbtagsT, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB2, MuoMuoSB2, EleEleSB2, MuoTauSB2, EleTauSB2);
  }
  //ELE-ELE ANALYSIS - SB2
  if(foundSB2Jet && foundSB2Electron1EE && foundSB2Electron2EE){
    FillTree(2,TreeSB2EleEle,SelectedSB2Jet,SelectedTauFake,SelectedMuonFake,SelectedMuon1Fake,SelectedMuon2Fake,SelectedElectronFake,SelectedSB2Electron1EE,
	     SelectedSB2Electron2EE, vertices, metRaw, met, uncorrmet, prunedMassSB2, nbtagsL, nbtagsM, nbtagsT, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB2, MuoMuoSB2, EleEleSB2, MuoTauSB2, EleTauSB2);
  }
  //MUO-TAU ANALYSIS - SB2
  if(foundSB2Jet && foundSB2MuonMT && foundSB2TauMT){
    FillTree(3,TreeSB2MuoTau,SelectedSB2Jet,SelectedSB2TauMT,SelectedSB2MuonMT,SelectedMuon1Fake,SelectedMuon2Fake,SelectedElectronFake,SelectedElectron1Fake,
	     SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMassSB2, nbtagsL, nbtagsM, nbtagsT, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB2, MuoMuoSB2, EleEleSB2, MuoTauSB2, EleTauSB2);
  }
  //ELE-TAU ANALYSIS - SB2
  if(foundSB2Jet && foundSB2ElectronET && foundSB2TauET){
    FillTree(4,TreeSB2EleTau,SelectedSB2Jet,SelectedSB2TauET,SelectedMuonFake,SelectedMuon1Fake,SelectedMuon2Fake,SelectedSB2ElectronET,SelectedElectron1Fake,
	     SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMassSB2, nbtagsL, nbtagsM, nbtagsT, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB2, MuoMuoSB2, EleEleSB2, MuoTauSB2, EleTauSB2);
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
Analyzer::beginJob()
{
  Service<TFileService> fs;
  Nevents = fs->make<TH1D>("Nevents", "Nevents", 3, -0.5, 2.5);

  TreeSignalEff = fs->make<TTree>("TreeSignalEff", "TreeSignalEff");
  TreeSignalEff->Branch("genEvent", &m_genEvent, "genEvent/f");

  TreeEleMuo = fs->make<TTree>("TreeEleMuo", "TreeEleMuo");
  TreeEleMuo->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeEleMuo->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeEleMuo->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeEleMuo->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeEleMuo->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeEleMuo->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeEleMuo->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeEleMuo->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeEleMuo->Branch("dPhiLep1Met", &m_dPhiLep1Met, "dPhiLep1Met/f");
  TreeEleMuo->Branch("dRLep1Met", &m_dRLep1Met, "dRLep1Met/f");
  TreeEleMuo->Branch("dRLep1Lep2", &m_dRLep1Lep2, "dRLep1Lep2/f");
  TreeEleMuo->Branch("dPhiLep2Met", &m_dPhiLep2Met, "dPhiLep2Met/f");
  TreeEleMuo->Branch("dRLep2Met", &m_dRLep2Met, "dRLep2Met/f");
  TreeEleMuo->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeEleMuo->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeEleMuo->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeEleMuo->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeEleMuo->Branch("lep1Pt", &m_lep1Pt, "lep1Pt/f");
  TreeEleMuo->Branch("lep1Eta", &m_lep1Eta, "lep1Eta/f");
  TreeEleMuo->Branch("lep1Charge", &m_lep1Charge, "lep1Charge/f");
  TreeEleMuo->Branch("lep1PFIso", &m_lep1PFIso, "lep1PFIso/f");
  TreeEleMuo->Branch("lep1CorrPFIso", &m_lep1CorrPFIso, "lep1CorrPFIso/f");
  TreeEleMuo->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeEleMuo->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeEleMuo->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeEleMuo->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeEleMuo->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeEleMuo->Branch("charge", &m_charge, "charge/f");
  TreeEleMuo->Branch("met", &m_met, "met/f");
  TreeEleMuo->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeEleMuo->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeEleMuo->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  TreeEleMuo->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeEleMuo->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeEleMuo->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeEleMuo->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeEleMuo->Branch("PtSvfit", &m_PtSvfit, "PtSvfit/f");
  TreeEleMuo->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeEleMuo->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeEleMuo->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeEleMuo->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeEleMuo->Branch("nbtagsL", &m_nbtagsL, "nbtagsL/i");
  TreeEleMuo->Branch("nbtagsM", &m_nbtagsM, "nbtagsM/i");
  TreeEleMuo->Branch("nbtagsT", &m_nbtagsT, "nbtagsT/i");
  TreeEleMuo->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeEleMuo->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeEleMuo->Branch("trigger", &m_trigger, "trigger/i");
  TreeEleMuo->Branch("sideband", &m_sideband, "sideband/i");
  TreeEleMuo->Branch("NVertices", &m_NVertices, "NVertices/i");
  TreeEleMuo->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeEleMuo->Branch("metPx", &m_metPx, "metPx/f");
  TreeEleMuo->Branch("metPy", &m_metPy, "metPy/f");
  TreeEleMuo->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeEleMuo->Branch("xsec", &m_xsec, "xsec/d");
  TreeEleMuo->Branch("lumi", &m_lumi, "lumi/d");
  TreeEleMuo->Branch("weight", &m_weight, "weight/f");
  TreeEleMuo->Branch("genEvent", &m_genEvent, "genEvent/f");
  TreeEleMuo->Branch("EleMuo", &m_EleMuo, "EleMuo/i");
  TreeEleMuo->Branch("MuoMuo", &m_MuoMuo, "MuoMuo/i");
  TreeEleMuo->Branch("EleEle", &m_EleEle, "EleEle/i");
  TreeEleMuo->Branch("MuoTau", &m_MuoTau, "MuoTau/i");
  TreeEleMuo->Branch("EleTau", &m_EleTau, "EleTau/i");

  TreeMuoMuo = fs->make<TTree>("TreeMuoMuo", "TreeMuoMuo");
  TreeMuoMuo->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeMuoMuo->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeMuoMuo->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeMuoMuo->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeMuoMuo->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeMuoMuo->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeMuoMuo->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeMuoMuo->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeMuoMuo->Branch("dPhiLep1Met", &m_dPhiLep1Met, "dPhiLep1Met/f");
  TreeMuoMuo->Branch("dRLep1Met", &m_dRLep1Met, "dRLep1Met/f");
  TreeMuoMuo->Branch("dRLep1Lep2", &m_dRLep1Lep2, "dRLep1Lep2/f");
  TreeMuoMuo->Branch("dPhiLep2Met", &m_dPhiLep2Met, "dPhiLep2Met/f");
  TreeMuoMuo->Branch("dRLep2Met", &m_dRLep2Met, "dRLep2Met/f");
  TreeMuoMuo->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeMuoMuo->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeMuoMuo->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeMuoMuo->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeMuoMuo->Branch("lep1Pt", &m_lep1Pt, "lep1Pt/f");
  TreeMuoMuo->Branch("lep1Eta", &m_lep1Eta, "lep1Eta/f");
  TreeMuoMuo->Branch("lep1Charge", &m_lep1Charge, "lep1Charge/f");
  TreeMuoMuo->Branch("lep1PFIso", &m_lep1PFIso, "lep1PFIso/f");
  TreeMuoMuo->Branch("lep1CorrPFIso", &m_lep1CorrPFIso, "lep1CorrPFIso/f");
  TreeMuoMuo->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeMuoMuo->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeMuoMuo->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeMuoMuo->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeMuoMuo->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeMuoMuo->Branch("charge", &m_charge, "charge/f");
  TreeMuoMuo->Branch("met", &m_met, "met/f");
  TreeMuoMuo->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeMuoMuo->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeMuoMuo->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  TreeMuoMuo->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeMuoMuo->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeMuoMuo->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeMuoMuo->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeMuoMuo->Branch("PtSvfit", &m_PtSvfit, "PtSvfit/f");
  TreeMuoMuo->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeMuoMuo->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeMuoMuo->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeMuoMuo->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeMuoMuo->Branch("nbtagsL", &m_nbtagsL, "nbtagsL/i");
  TreeMuoMuo->Branch("nbtagsM", &m_nbtagsM, "nbtagsM/i");
  TreeMuoMuo->Branch("nbtagsT", &m_nbtagsT, "nbtagsT/i");
  TreeMuoMuo->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeMuoMuo->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeMuoMuo->Branch("trigger", &m_trigger, "trigger/i");
  TreeMuoMuo->Branch("sideband", &m_sideband, "sideband/i");
  TreeMuoMuo->Branch("NVertices", &m_NVertices, "NVertices/i");
  TreeMuoMuo->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeMuoMuo->Branch("metPx", &m_metPx, "metPx/f");
  TreeMuoMuo->Branch("metPy", &m_metPy, "metPy/f");
  TreeMuoMuo->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeMuoMuo->Branch("xsec", &m_xsec, "xsec/d");
  TreeMuoMuo->Branch("lumi", &m_lumi, "lumi/d");
  TreeMuoMuo->Branch("weight", &m_weight, "weight/f");
  TreeMuoMuo->Branch("genEvent", &m_genEvent, "genEvent/f");
  TreeMuoMuo->Branch("EleMuo", &m_EleMuo, "EleMuo/i");
  TreeMuoMuo->Branch("MuoMuo", &m_MuoMuo, "MuoMuo/i");
  TreeMuoMuo->Branch("EleEle", &m_EleEle, "EleEle/i");
  TreeMuoMuo->Branch("MuoTau", &m_MuoTau, "MuoTau/i");
  TreeMuoMuo->Branch("EleTau", &m_EleTau, "EleTau/i");

  TreeEleEle = fs->make<TTree>("TreeEleEle", "TreeEleEle");
  TreeEleEle->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeEleEle->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeEleEle->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeEleEle->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeEleEle->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeEleEle->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeEleEle->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeEleEle->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeEleEle->Branch("dPhiLep1Met", &m_dPhiLep1Met, "dPhiLep1Met/f");
  TreeEleEle->Branch("dRLep1Met", &m_dRLep1Met, "dRLep1Met/f");
  TreeEleEle->Branch("dRLep1Lep2", &m_dRLep1Lep2, "dRLep1Lep2/f");
  TreeEleEle->Branch("dPhiLep2Met", &m_dPhiLep2Met, "dPhiLep2Met/f");
  TreeEleEle->Branch("dRLep2Met", &m_dRLep2Met, "dRLep2Met/f");
  TreeEleEle->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeEleEle->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeEleEle->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeEleEle->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeEleEle->Branch("lep1Pt", &m_lep1Pt, "lep1Pt/f");
  TreeEleEle->Branch("lep1Eta", &m_lep1Eta, "lep1Eta/f");
  TreeEleEle->Branch("lep1Charge", &m_lep1Charge, "lep1Charge/f");
  TreeEleEle->Branch("lep1PFIso", &m_lep1PFIso, "lep1PFIso/f");
  TreeEleEle->Branch("lep1CorrPFIso", &m_lep1CorrPFIso, "lep1CorrPFIso/f");
  TreeEleEle->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeEleEle->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeEleEle->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeEleEle->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeEleEle->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeEleEle->Branch("charge", &m_charge, "charge/f");
  TreeEleEle->Branch("met", &m_met, "met/f");
  TreeEleEle->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeEleEle->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeEleEle->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  TreeEleEle->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeEleEle->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeEleEle->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeEleEle->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeEleEle->Branch("PtSvfit", &m_PtSvfit, "PtSvfit/f");
  TreeEleEle->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeEleEle->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeEleEle->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeEleEle->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeEleEle->Branch("nbtagsL", &m_nbtagsL, "nbtagsL/i");
  TreeEleEle->Branch("nbtagsM", &m_nbtagsM, "nbtagsM/i");
  TreeEleEle->Branch("nbtagsT", &m_nbtagsT, "nbtagsT/i");
  TreeEleEle->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeEleEle->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeEleEle->Branch("trigger", &m_trigger, "trigger/i");
  TreeEleEle->Branch("sideband", &m_sideband, "sideband/i");
  TreeEleEle->Branch("NVertices", &m_NVertices, "NVertices/i");
  TreeEleEle->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeEleEle->Branch("metPx", &m_metPx, "metPx/f");
  TreeEleEle->Branch("metPy", &m_metPy, "metPy/f");
  TreeEleEle->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeEleEle->Branch("xsec", &m_xsec, "xsec/d");
  TreeEleEle->Branch("lumi", &m_lumi, "lumi/d");
  TreeEleEle->Branch("weight", &m_weight, "weight/f");
  TreeEleEle->Branch("genEvent", &m_genEvent, "genEvent/f");
  TreeEleEle->Branch("EleMuo", &m_EleMuo, "EleMuo/i");
  TreeEleEle->Branch("MuoMuo", &m_MuoMuo, "MuoMuo/i");
  TreeEleEle->Branch("EleEle", &m_EleEle, "EleEle/i");
  TreeEleEle->Branch("MuoTau", &m_MuoTau, "MuoTau/i");
  TreeEleEle->Branch("EleTau", &m_EleTau, "EleTau/i");

  TreeMuoTau = fs->make<TTree>("TreeMuoTau", "TreeMuoTau");
  TreeMuoTau->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeMuoTau->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeMuoTau->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeMuoTau->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeMuoTau->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeMuoTau->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeMuoTau->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeMuoTau->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeMuoTau->Branch("dPhiLep1Met", &m_dPhiLep1Met, "dPhiLep1Met/f");
  TreeMuoTau->Branch("dRLep1Met", &m_dRLep1Met, "dRLep1Met/f");
  TreeMuoTau->Branch("dRLep1Lep2", &m_dRLep1Lep2, "dRLep1Lep2/f");
  TreeMuoTau->Branch("dPhiLep2Met", &m_dPhiLep2Met, "dPhiLep2Met/f");
  TreeMuoTau->Branch("dRLep2Met", &m_dRLep2Met, "dRLep2Met/f");
  TreeMuoTau->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeMuoTau->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeMuoTau->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeMuoTau->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeMuoTau->Branch("lep1Pt", &m_lep1Pt, "lep1Pt/f");
  TreeMuoTau->Branch("lep1Eta", &m_lep1Eta, "lep1Eta/f");
  TreeMuoTau->Branch("lep1Charge", &m_lep1Charge, "lep1Charge/f");
  TreeMuoTau->Branch("lep1PFIso", &m_lep1PFIso, "lep1PFIso/f");
  TreeMuoTau->Branch("lep1CorrPFIso", &m_lep1CorrPFIso, "lep1CorrPFIso/f");
  TreeMuoTau->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeMuoTau->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeMuoTau->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeMuoTau->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeMuoTau->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
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
  TreeMuoTau->Branch("EleMuo", &m_EleMuo, "EleMuo/i");
  TreeMuoTau->Branch("MuoMuo", &m_MuoMuo, "MuoMuo/i");
  TreeMuoTau->Branch("EleEle", &m_EleEle, "EleEle/i");
  TreeMuoTau->Branch("MuoTau", &m_MuoTau, "MuoTau/i");
  TreeMuoTau->Branch("EleTau", &m_EleTau, "EleTau/i");

  TreeEleTau = fs->make<TTree>("TreeEleTau", "TreeEleTau");
  TreeEleTau->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeEleTau->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeEleTau->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeEleTau->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeEleTau->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeEleTau->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeEleTau->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeEleTau->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeEleTau->Branch("dPhiLep1Met", &m_dPhiLep1Met, "dPhiLep1Met/f");
  TreeEleTau->Branch("dRLep1Met", &m_dRLep1Met, "dRLep1Met/f");
  TreeEleTau->Branch("dRLep1Lep2", &m_dRLep1Lep2, "dRLep1Lep2/f");
  TreeEleTau->Branch("dPhiLep2Met", &m_dPhiLep2Met, "dPhiLep2Met/f");
  TreeEleTau->Branch("dRLep2Met", &m_dRLep2Met, "dRLep2Met/f");
  TreeEleTau->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeEleTau->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeEleTau->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeEleTau->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeEleTau->Branch("lep1Pt", &m_lep1Pt, "lep1Pt/f");
  TreeEleTau->Branch("lep1Eta", &m_lep1Eta, "lep1Eta/f");
  TreeEleTau->Branch("lep1Charge", &m_lep1Charge, "lep1Charge/f");
  TreeEleTau->Branch("lep1PFIso", &m_lep1PFIso, "lep1PFIso/f");
  TreeEleTau->Branch("lep1CorrPFIso", &m_lep1CorrPFIso, "lep1CorrPFIso/f");
  TreeEleTau->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeEleTau->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeEleTau->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeEleTau->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeEleTau->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
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
  TreeEleTau->Branch("EleMuo", &m_EleMuo, "EleMuo/i");
  TreeEleTau->Branch("MuoMuo", &m_MuoMuo, "MuoMuo/i");
  TreeEleTau->Branch("EleEle", &m_EleEle, "EleEle/i");
  TreeEleTau->Branch("MuoTau", &m_MuoTau, "MuoTau/i");
  TreeEleTau->Branch("EleTau", &m_EleTau, "EleTau/i");

  TreeSB1EleMuo = fs->make<TTree>("TreeSB1EleMuo", "TreeSB1EleMuo");
  TreeSB1EleMuo->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB1EleMuo->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB1EleMuo->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB1EleMuo->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB1EleMuo->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB1EleMuo->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB1EleMuo->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeSB1EleMuo->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeSB1EleMuo->Branch("dPhiLep1Met", &m_dPhiLep1Met, "dPhiLep1Met/f");
  TreeSB1EleMuo->Branch("dRLep1Met", &m_dRLep1Met, "dRLep1Met/f");
  TreeSB1EleMuo->Branch("dRLep1Lep2", &m_dRLep1Lep2, "dRLep1Lep2/f");
  TreeSB1EleMuo->Branch("dPhiLep2Met", &m_dPhiLep2Met, "dPhiLep2Met/f");
  TreeSB1EleMuo->Branch("dRLep2Met", &m_dRLep2Met, "dRLep2Met/f");
  TreeSB1EleMuo->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB1EleMuo->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB1EleMuo->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB1EleMuo->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB1EleMuo->Branch("lep1Pt", &m_lep1Pt, "lep1Pt/f");
  TreeSB1EleMuo->Branch("lep1Eta", &m_lep1Eta, "lep1Eta/f");
  TreeSB1EleMuo->Branch("lep1Charge", &m_lep1Charge, "lep1Charge/f");
  TreeSB1EleMuo->Branch("lep1PFIso", &m_lep1PFIso, "lep1PFIso/f");
  TreeSB1EleMuo->Branch("lep1CorrPFIso", &m_lep1CorrPFIso, "lep1CorrPFIso/f");
  TreeSB1EleMuo->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB1EleMuo->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB1EleMuo->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB1EleMuo->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB1EleMuo->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeSB1EleMuo->Branch("charge", &m_charge, "charge/f");
  TreeSB1EleMuo->Branch("met", &m_met, "met/f");
  TreeSB1EleMuo->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeSB1EleMuo->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeSB1EleMuo->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  TreeSB1EleMuo->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeSB1EleMuo->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeSB1EleMuo->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeSB1EleMuo->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeSB1EleMuo->Branch("PtSvfit", &m_PtSvfit, "PtSvfit/f");
  TreeSB1EleMuo->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeSB1EleMuo->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeSB1EleMuo->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeSB1EleMuo->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeSB1EleMuo->Branch("nbtagsL", &m_nbtagsL, "nbtagsL/i");
  TreeSB1EleMuo->Branch("nbtagsM", &m_nbtagsM, "nbtagsM/i");
  TreeSB1EleMuo->Branch("nbtagsT", &m_nbtagsT, "nbtagsT/i");
  TreeSB1EleMuo->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeSB1EleMuo->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeSB1EleMuo->Branch("trigger", &m_trigger, "trigger/i");
  TreeSB1EleMuo->Branch("sideband", &m_sideband, "sideband/i");
  TreeSB1EleMuo->Branch("NVertices", &m_NVertices, "NVertices/i");
  TreeSB1EleMuo->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeSB1EleMuo->Branch("metPx", &m_metPx, "metPx/f");
  TreeSB1EleMuo->Branch("metPy", &m_metPy, "metPy/f");
  TreeSB1EleMuo->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeSB1EleMuo->Branch("xsec", &m_xsec, "xsec/d");
  TreeSB1EleMuo->Branch("lumi", &m_lumi, "lumi/d");
  TreeSB1EleMuo->Branch("weight", &m_weight, "weight/f");
  TreeSB1EleMuo->Branch("genEvent", &m_genEvent, "genEvent/f");
  TreeSB1EleMuo->Branch("EleMuo", &m_EleMuo, "EleMuo/i");
  TreeSB1EleMuo->Branch("MuoMuo", &m_MuoMuo, "MuoMuo/i");
  TreeSB1EleMuo->Branch("EleEle", &m_EleEle, "EleEle/i");
  TreeSB1EleMuo->Branch("MuoTau", &m_MuoTau, "MuoTau/i");
  TreeSB1EleMuo->Branch("EleTau", &m_EleTau, "EleTau/i");

  TreeSB1MuoMuo = fs->make<TTree>("TreeSB1MuoMuo", "TreeSB1MuoMuo");
  TreeSB1MuoMuo->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB1MuoMuo->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB1MuoMuo->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB1MuoMuo->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB1MuoMuo->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB1MuoMuo->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB1MuoMuo->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeSB1MuoMuo->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeSB1MuoMuo->Branch("dPhiLep1Met", &m_dPhiLep1Met, "dPhiLep1Met/f");
  TreeSB1MuoMuo->Branch("dRLep1Met", &m_dRLep1Met, "dRLep1Met/f");
  TreeSB1MuoMuo->Branch("dRLep1Lep2", &m_dRLep1Lep2, "dRLep1Lep2/f");
  TreeSB1MuoMuo->Branch("dPhiLep2Met", &m_dPhiLep2Met, "dPhiLep2Met/f");
  TreeSB1MuoMuo->Branch("dRLep2Met", &m_dRLep2Met, "dRLep2Met/f");
  TreeSB1MuoMuo->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB1MuoMuo->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB1MuoMuo->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB1MuoMuo->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB1MuoMuo->Branch("lep1Pt", &m_lep1Pt, "lep1Pt/f");
  TreeSB1MuoMuo->Branch("lep1Eta", &m_lep1Eta, "lep1Eta/f");
  TreeSB1MuoMuo->Branch("lep1Charge", &m_lep1Charge, "lep1Charge/f");
  TreeSB1MuoMuo->Branch("lep1PFIso", &m_lep1PFIso, "lep1PFIso/f");
  TreeSB1MuoMuo->Branch("lep1CorrPFIso", &m_lep1CorrPFIso, "lep1CorrPFIso/f");
  TreeSB1MuoMuo->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB1MuoMuo->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB1MuoMuo->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB1MuoMuo->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB1MuoMuo->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeSB1MuoMuo->Branch("charge", &m_charge, "charge/f");
  TreeSB1MuoMuo->Branch("met", &m_met, "met/f");
  TreeSB1MuoMuo->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeSB1MuoMuo->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeSB1MuoMuo->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  TreeSB1MuoMuo->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeSB1MuoMuo->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeSB1MuoMuo->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeSB1MuoMuo->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeSB1MuoMuo->Branch("PtSvfit", &m_PtSvfit, "PtSvfit/f");
  TreeSB1MuoMuo->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeSB1MuoMuo->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeSB1MuoMuo->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeSB1MuoMuo->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeSB1MuoMuo->Branch("nbtagsL", &m_nbtagsL, "nbtagsL/i");
  TreeSB1MuoMuo->Branch("nbtagsM", &m_nbtagsM, "nbtagsM/i");
  TreeSB1MuoMuo->Branch("nbtagsT", &m_nbtagsT, "nbtagsT/i");
  TreeSB1MuoMuo->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeSB1MuoMuo->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeSB1MuoMuo->Branch("trigger", &m_trigger, "trigger/i");
  TreeSB1MuoMuo->Branch("sideband", &m_sideband, "sideband/i");
  TreeSB1MuoMuo->Branch("NVertices", &m_NVertices, "NVertices/i");
  TreeSB1MuoMuo->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeSB1MuoMuo->Branch("metPx", &m_metPx, "metPx/f");
  TreeSB1MuoMuo->Branch("metPy", &m_metPy, "metPy/f");
  TreeSB1MuoMuo->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeSB1MuoMuo->Branch("xsec", &m_xsec, "xsec/d");
  TreeSB1MuoMuo->Branch("lumi", &m_lumi, "lumi/d");
  TreeSB1MuoMuo->Branch("weight", &m_weight, "weight/f");
  TreeSB1MuoMuo->Branch("genEvent", &m_genEvent, "genEvent/f");
  TreeSB1MuoMuo->Branch("EleMuo", &m_EleMuo, "EleMuo/i");
  TreeSB1MuoMuo->Branch("MuoMuo", &m_MuoMuo, "MuoMuo/i");
  TreeSB1MuoMuo->Branch("EleEle", &m_EleEle, "EleEle/i");
  TreeSB1MuoMuo->Branch("MuoTau", &m_MuoTau, "MuoTau/i");
  TreeSB1MuoMuo->Branch("EleTau", &m_EleTau, "EleTau/i");

  TreeSB1EleEle = fs->make<TTree>("TreeSB1EleEle", "TreeSB1EleEle");
  TreeSB1EleEle->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB1EleEle->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB1EleEle->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB1EleEle->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB1EleEle->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB1EleEle->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB1EleEle->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeSB1EleEle->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeSB1EleEle->Branch("dPhiLep1Met", &m_dPhiLep1Met, "dPhiLep1Met/f");
  TreeSB1EleEle->Branch("dRLep1Met", &m_dRLep1Met, "dRLep1Met/f");
  TreeSB1EleEle->Branch("dRLep1Lep2", &m_dRLep1Lep2, "dRLep1Lep2/f");
  TreeSB1EleEle->Branch("dPhiLep2Met", &m_dPhiLep2Met, "dPhiLep2Met/f");
  TreeSB1EleEle->Branch("dRLep2Met", &m_dRLep2Met, "dRLep2Met/f");
  TreeSB1EleEle->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB1EleEle->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB1EleEle->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB1EleEle->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB1EleEle->Branch("lep1Pt", &m_lep1Pt, "lep1Pt/f");
  TreeSB1EleEle->Branch("lep1Eta", &m_lep1Eta, "lep1Eta/f");
  TreeSB1EleEle->Branch("lep1Charge", &m_lep1Charge, "lep1Charge/f");
  TreeSB1EleEle->Branch("lep1PFIso", &m_lep1PFIso, "lep1PFIso/f");
  TreeSB1EleEle->Branch("lep1CorrPFIso", &m_lep1CorrPFIso, "lep1CorrPFIso/f");
  TreeSB1EleEle->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB1EleEle->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB1EleEle->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB1EleEle->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB1EleEle->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeSB1EleEle->Branch("charge", &m_charge, "charge/f");
  TreeSB1EleEle->Branch("met", &m_met, "met/f");
  TreeSB1EleEle->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeSB1EleEle->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeSB1EleEle->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  TreeSB1EleEle->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeSB1EleEle->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeSB1EleEle->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeSB1EleEle->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeSB1EleEle->Branch("PtSvfit", &m_PtSvfit, "PtSvfit/f");
  TreeSB1EleEle->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeSB1EleEle->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeSB1EleEle->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeSB1EleEle->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeSB1EleEle->Branch("nbtagsL", &m_nbtagsL, "nbtagsL/i");
  TreeSB1EleEle->Branch("nbtagsM", &m_nbtagsM, "nbtagsM/i");
  TreeSB1EleEle->Branch("nbtagsT", &m_nbtagsT, "nbtagsT/i");
  TreeSB1EleEle->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeSB1EleEle->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeSB1EleEle->Branch("trigger", &m_trigger, "trigger/i");
  TreeSB1EleEle->Branch("sideband", &m_sideband, "sideband/i");
  TreeSB1EleEle->Branch("NVertices", &m_NVertices, "NVertices/i");
  TreeSB1EleEle->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeSB1EleEle->Branch("metPx", &m_metPx, "metPx/f");
  TreeSB1EleEle->Branch("metPy", &m_metPy, "metPy/f");
  TreeSB1EleEle->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeSB1EleEle->Branch("xsec", &m_xsec, "xsec/d");
  TreeSB1EleEle->Branch("lumi", &m_lumi, "lumi/d");
  TreeSB1EleEle->Branch("weight", &m_weight, "weight/f");
  TreeSB1EleEle->Branch("genEvent", &m_genEvent, "genEvent/f");
  TreeSB1EleEle->Branch("EleMuo", &m_EleMuo, "EleMuo/i");
  TreeSB1EleEle->Branch("MuoMuo", &m_MuoMuo, "MuoMuo/i");
  TreeSB1EleEle->Branch("EleEle", &m_EleEle, "EleEle/i");
  TreeSB1EleEle->Branch("MuoTau", &m_MuoTau, "MuoTau/i");
  TreeSB1EleEle->Branch("EleTau", &m_EleTau, "EleTau/i");

  TreeSB1MuoTau = fs->make<TTree>("TreeSB1MuoTau", "TreeSB1MuoTau");
  TreeSB1MuoTau->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB1MuoTau->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB1MuoTau->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB1MuoTau->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB1MuoTau->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB1MuoTau->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB1MuoTau->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeSB1MuoTau->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeSB1MuoTau->Branch("dPhiLep1Met", &m_dPhiLep1Met, "dPhiLep1Met/f");
  TreeSB1MuoTau->Branch("dRLep1Met", &m_dRLep1Met, "dRLep1Met/f");
  TreeSB1MuoTau->Branch("dRLep1Lep2", &m_dRLep1Lep2, "dRLep1Lep2/f");
  TreeSB1MuoTau->Branch("dPhiLep2Met", &m_dPhiLep2Met, "dPhiLep2Met/f");
  TreeSB1MuoTau->Branch("dRLep2Met", &m_dRLep2Met, "dRLep2Met/f");
  TreeSB1MuoTau->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB1MuoTau->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB1MuoTau->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB1MuoTau->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB1MuoTau->Branch("lep1Pt", &m_lep1Pt, "lep1Pt/f");
  TreeSB1MuoTau->Branch("lep1Eta", &m_lep1Eta, "lep1Eta/f");
  TreeSB1MuoTau->Branch("lep1Charge", &m_lep1Charge, "lep1Charge/f");
  TreeSB1MuoTau->Branch("lep1PFIso", &m_lep1PFIso, "lep1PFIso/f");
  TreeSB1MuoTau->Branch("lep1CorrPFIso", &m_lep1CorrPFIso, "lep1CorrPFIso/f");
  TreeSB1MuoTau->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB1MuoTau->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB1MuoTau->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB1MuoTau->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB1MuoTau->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
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
  TreeSB1MuoTau->Branch("EleMuo", &m_EleMuo, "EleMuo/i");
  TreeSB1MuoTau->Branch("MuoMuo", &m_MuoMuo, "MuoMuo/i");
  TreeSB1MuoTau->Branch("EleEle", &m_EleEle, "EleEle/i");
  TreeSB1MuoTau->Branch("MuoTau", &m_MuoTau, "MuoTau/i");
  TreeSB1MuoTau->Branch("EleTau", &m_EleTau, "EleTau/i");

  TreeSB1EleTau = fs->make<TTree>("TreeSB1EleTau", "TreeSB1EleTau");
  TreeSB1EleTau->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB1EleTau->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB1EleTau->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB1EleTau->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB1EleTau->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB1EleTau->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB1EleTau->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeSB1EleTau->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeSB1EleTau->Branch("dPhiLep1Met", &m_dPhiLep1Met, "dPhiLep1Met/f");
  TreeSB1EleTau->Branch("dRLep1Met", &m_dRLep1Met, "dRLep1Met/f");
  TreeSB1EleTau->Branch("dRLep1Lep2", &m_dRLep1Lep2, "dRLep1Lep2/f");
  TreeSB1EleTau->Branch("dPhiLep2Met", &m_dPhiLep2Met, "dPhiLep2Met/f");
  TreeSB1EleTau->Branch("dRLep2Met", &m_dRLep2Met, "dRLep2Met/f");
  TreeSB1EleTau->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB1EleTau->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB1EleTau->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB1EleTau->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB1EleTau->Branch("lep1Pt", &m_lep1Pt, "lep1Pt/f");
  TreeSB1EleTau->Branch("lep1Eta", &m_lep1Eta, "lep1Eta/f");
  TreeSB1EleTau->Branch("lep1Charge", &m_lep1Charge, "lep1Charge/f");
  TreeSB1EleTau->Branch("lep1PFIso", &m_lep1PFIso, "lep1PFIso/f");
  TreeSB1EleTau->Branch("lep1CorrPFIso", &m_lep1CorrPFIso, "lep1CorrPFIso/f");
  TreeSB1EleTau->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB1EleTau->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB1EleTau->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB1EleTau->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB1EleTau->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
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
  TreeSB1EleTau->Branch("EleMuo", &m_EleMuo, "EleMuo/i");
  TreeSB1EleTau->Branch("MuoMuo", &m_MuoMuo, "MuoMuo/i");
  TreeSB1EleTau->Branch("EleEle", &m_EleEle, "EleEle/i");
  TreeSB1EleTau->Branch("MuoTau", &m_MuoTau, "MuoTau/i");
  TreeSB1EleTau->Branch("EleTau", &m_EleTau, "EleTau/i");

  TreeSB2EleMuo = fs->make<TTree>("TreeSB2EleMuo", "TreeSB2EleMuo");
  TreeSB2EleMuo->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB2EleMuo->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB2EleMuo->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB2EleMuo->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB2EleMuo->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB2EleMuo->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB2EleMuo->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeSB2EleMuo->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeSB2EleMuo->Branch("dPhiLep1Met", &m_dPhiLep1Met, "dPhiLep1Met/f");
  TreeSB2EleMuo->Branch("dRLep1Met", &m_dRLep1Met, "dRLep1Met/f");
  TreeSB2EleMuo->Branch("dRLep1Lep2", &m_dRLep1Lep2, "dRLep1Lep2/f");
  TreeSB2EleMuo->Branch("dPhiLep2Met", &m_dPhiLep2Met, "dPhiLep2Met/f");
  TreeSB2EleMuo->Branch("dRLep2Met", &m_dRLep2Met, "dRLep2Met/f");
  TreeSB2EleMuo->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB2EleMuo->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB2EleMuo->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB2EleMuo->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB2EleMuo->Branch("lep1Pt", &m_lep1Pt, "lep1Pt/f");
  TreeSB2EleMuo->Branch("lep1Eta", &m_lep1Eta, "lep1Eta/f");
  TreeSB2EleMuo->Branch("lep1Charge", &m_lep1Charge, "lep1Charge/f");
  TreeSB2EleMuo->Branch("lep1PFIso", &m_lep1PFIso, "lep1PFIso/f");
  TreeSB2EleMuo->Branch("lep1CorrPFIso", &m_lep1CorrPFIso, "lep1CorrPFIso/f");
  TreeSB2EleMuo->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB2EleMuo->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB2EleMuo->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB2EleMuo->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB2EleMuo->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeSB2EleMuo->Branch("charge", &m_charge, "charge/f");
  TreeSB2EleMuo->Branch("met", &m_met, "met/f");
  TreeSB2EleMuo->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeSB2EleMuo->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeSB2EleMuo->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  TreeSB2EleMuo->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeSB2EleMuo->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeSB2EleMuo->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeSB2EleMuo->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeSB2EleMuo->Branch("PtSvfit", &m_PtSvfit, "PtSvfit/f");
  TreeSB2EleMuo->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeSB2EleMuo->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeSB2EleMuo->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeSB2EleMuo->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeSB2EleMuo->Branch("nbtagsL", &m_nbtagsL, "nbtagsL/i");
  TreeSB2EleMuo->Branch("nbtagsM", &m_nbtagsM, "nbtagsM/i");
  TreeSB2EleMuo->Branch("nbtagsT", &m_nbtagsT, "nbtagsT/i");
  TreeSB2EleMuo->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeSB2EleMuo->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeSB2EleMuo->Branch("trigger", &m_trigger, "trigger/i");
  TreeSB2EleMuo->Branch("sideband", &m_sideband, "sideband/i");
  TreeSB2EleMuo->Branch("NVertices", &m_NVertices, "NVertices/i");
  TreeSB2EleMuo->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeSB2EleMuo->Branch("metPx", &m_metPx, "metPx/f");
  TreeSB2EleMuo->Branch("metPy", &m_metPy, "metPy/f");
  TreeSB2EleMuo->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeSB2EleMuo->Branch("xsec", &m_xsec, "xsec/d");
  TreeSB2EleMuo->Branch("lumi", &m_lumi, "lumi/d");
  TreeSB2EleMuo->Branch("weight", &m_weight, "weight/f");
  TreeSB2EleMuo->Branch("genEvent", &m_genEvent, "genEvent/f");
  TreeSB2EleMuo->Branch("EleMuo", &m_EleMuo, "EleMuo/i");
  TreeSB2EleMuo->Branch("MuoMuo", &m_MuoMuo, "MuoMuo/i");
  TreeSB2EleMuo->Branch("EleEle", &m_EleEle, "EleEle/i");
  TreeSB2EleMuo->Branch("MuoTau", &m_MuoTau, "MuoTau/i");
  TreeSB2EleMuo->Branch("EleTau", &m_EleTau, "EleTau/i");

  TreeSB2MuoMuo = fs->make<TTree>("TreeSB2MuoMuo", "TreeSB2MuoMuo");
  TreeSB2MuoMuo->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB2MuoMuo->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB2MuoMuo->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB2MuoMuo->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB2MuoMuo->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB2MuoMuo->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB2MuoMuo->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeSB2MuoMuo->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeSB2MuoMuo->Branch("dPhiLep1Met", &m_dPhiLep1Met, "dPhiLep1Met/f");
  TreeSB2MuoMuo->Branch("dRLep1Met", &m_dRLep1Met, "dRLep1Met/f");
  TreeSB2MuoMuo->Branch("dRLep1Lep2", &m_dRLep1Lep2, "dRLep1Lep2/f");
  TreeSB2MuoMuo->Branch("dPhiLep2Met", &m_dPhiLep2Met, "dPhiLep2Met/f");
  TreeSB2MuoMuo->Branch("dRLep2Met", &m_dRLep2Met, "dRLep2Met/f");
  TreeSB2MuoMuo->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB2MuoMuo->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB2MuoMuo->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB2MuoMuo->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB2MuoMuo->Branch("lep1Pt", &m_lep1Pt, "lep1Pt/f");
  TreeSB2MuoMuo->Branch("lep1Eta", &m_lep1Eta, "lep1Eta/f");
  TreeSB2MuoMuo->Branch("lep1Charge", &m_lep1Charge, "lep1Charge/f");
  TreeSB2MuoMuo->Branch("lep1PFIso", &m_lep1PFIso, "lep1PFIso/f");
  TreeSB2MuoMuo->Branch("lep1CorrPFIso", &m_lep1CorrPFIso, "lep1CorrPFIso/f");
  TreeSB2MuoMuo->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB2MuoMuo->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB2MuoMuo->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB2MuoMuo->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB2MuoMuo->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeSB2MuoMuo->Branch("charge", &m_charge, "charge/f");
  TreeSB2MuoMuo->Branch("met", &m_met, "met/f");
  TreeSB2MuoMuo->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeSB2MuoMuo->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeSB2MuoMuo->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  TreeSB2MuoMuo->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeSB2MuoMuo->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeSB2MuoMuo->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeSB2MuoMuo->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeSB2MuoMuo->Branch("PtSvfit", &m_PtSvfit, "PtSvfit/f");
  TreeSB2MuoMuo->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeSB2MuoMuo->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeSB2MuoMuo->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeSB2MuoMuo->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeSB2MuoMuo->Branch("nbtagsL", &m_nbtagsL, "nbtagsL/i");
  TreeSB2MuoMuo->Branch("nbtagsM", &m_nbtagsM, "nbtagsM/i");
  TreeSB2MuoMuo->Branch("nbtagsT", &m_nbtagsT, "nbtagsT/i");
  TreeSB2MuoMuo->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeSB2MuoMuo->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeSB2MuoMuo->Branch("trigger", &m_trigger, "trigger/i");
  TreeSB2MuoMuo->Branch("sideband", &m_sideband, "sideband/i");
  TreeSB2MuoMuo->Branch("NVertices", &m_NVertices, "NVertices/i");
  TreeSB2MuoMuo->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeSB2MuoMuo->Branch("metPx", &m_metPx, "metPx/f");
  TreeSB2MuoMuo->Branch("metPy", &m_metPy, "metPy/f");
  TreeSB2MuoMuo->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeSB2MuoMuo->Branch("xsec", &m_xsec, "xsec/d");
  TreeSB2MuoMuo->Branch("lumi", &m_lumi, "lumi/d");
  TreeSB2MuoMuo->Branch("weight", &m_weight, "weight/f");
  TreeSB2MuoMuo->Branch("genEvent", &m_genEvent, "genEvent/f");
  TreeSB2MuoMuo->Branch("EleMuo", &m_EleMuo, "EleMuo/i");
  TreeSB2MuoMuo->Branch("MuoMuo", &m_MuoMuo, "MuoMuo/i");
  TreeSB2MuoMuo->Branch("EleEle", &m_EleEle, "EleEle/i");
  TreeSB2MuoMuo->Branch("MuoTau", &m_MuoTau, "MuoTau/i");
  TreeSB2MuoMuo->Branch("EleTau", &m_EleTau, "EleTau/i");

  TreeSB2EleEle = fs->make<TTree>("TreeSB2EleEle", "TreeSB2EleEle");
  TreeSB2EleEle->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB2EleEle->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB2EleEle->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB2EleEle->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB2EleEle->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB2EleEle->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB2EleEle->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeSB2EleEle->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeSB2EleEle->Branch("dPhiLep1Met", &m_dPhiLep1Met, "dPhiLep1Met/f");
  TreeSB2EleEle->Branch("dRLep1Met", &m_dRLep1Met, "dRLep1Met/f");
  TreeSB2EleEle->Branch("dRLep1Lep2", &m_dRLep1Lep2, "dRLep1Lep2/f");
  TreeSB2EleEle->Branch("dPhiLep2Met", &m_dPhiLep2Met, "dPhiLep2Met/f");
  TreeSB2EleEle->Branch("dRLep2Met", &m_dRLep2Met, "dRLep2Met/f");
  TreeSB2EleEle->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB2EleEle->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB2EleEle->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB2EleEle->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB2EleEle->Branch("lep1Pt", &m_lep1Pt, "lep1Pt/f");
  TreeSB2EleEle->Branch("lep1Eta", &m_lep1Eta, "lep1Eta/f");
  TreeSB2EleEle->Branch("lep1Charge", &m_lep1Charge, "lep1Charge/f");
  TreeSB2EleEle->Branch("lep1PFIso", &m_lep1PFIso, "lep1PFIso/f");
  TreeSB2EleEle->Branch("lep1CorrPFIso", &m_lep1CorrPFIso, "lep1CorrPFIso/f");
  TreeSB2EleEle->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB2EleEle->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB2EleEle->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB2EleEle->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB2EleEle->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeSB2EleEle->Branch("charge", &m_charge, "charge/f");
  TreeSB2EleEle->Branch("met", &m_met, "met/f");
  TreeSB2EleEle->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeSB2EleEle->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeSB2EleEle->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  TreeSB2EleEle->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeSB2EleEle->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeSB2EleEle->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeSB2EleEle->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeSB2EleEle->Branch("PtSvfit", &m_PtSvfit, "PtSvfit/f");
  TreeSB2EleEle->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeSB2EleEle->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeSB2EleEle->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeSB2EleEle->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeSB2EleEle->Branch("nbtagsL", &m_nbtagsL, "nbtagsL/i");
  TreeSB2EleEle->Branch("nbtagsM", &m_nbtagsM, "nbtagsM/i");
  TreeSB2EleEle->Branch("nbtagsT", &m_nbtagsT, "nbtagsT/i");
  TreeSB2EleEle->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeSB2EleEle->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeSB2EleEle->Branch("trigger", &m_trigger, "trigger/i");
  TreeSB2EleEle->Branch("sideband", &m_sideband, "sideband/i");
  TreeSB2EleEle->Branch("NVertices", &m_NVertices, "NVertices/i");
  TreeSB2EleEle->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeSB2EleEle->Branch("metPx", &m_metPx, "metPx/f");
  TreeSB2EleEle->Branch("metPy", &m_metPy, "metPy/f");
  TreeSB2EleEle->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeSB2EleEle->Branch("xsec", &m_xsec, "xsec/d");
  TreeSB2EleEle->Branch("lumi", &m_lumi, "lumi/d");
  TreeSB2EleEle->Branch("weight", &m_weight, "weight/f");
  TreeSB2EleEle->Branch("genEvent", &m_genEvent, "genEvent/f");
  TreeSB2EleEle->Branch("EleMuo", &m_EleMuo, "EleMuo/i");
  TreeSB2EleEle->Branch("MuoMuo", &m_MuoMuo, "MuoMuo/i");
  TreeSB2EleEle->Branch("EleEle", &m_EleEle, "EleEle/i");
  TreeSB2EleEle->Branch("MuoTau", &m_MuoTau, "MuoTau/i");
  TreeSB2EleEle->Branch("EleTau", &m_EleTau, "EleTau/i");

  TreeSB2MuoTau = fs->make<TTree>("TreeSB2MuoTau", "TreeSB2MuoTau");
  TreeSB2MuoTau->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB2MuoTau->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB2MuoTau->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB2MuoTau->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB2MuoTau->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB2MuoTau->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB2MuoTau->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeSB2MuoTau->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeSB2MuoTau->Branch("dPhiLep1Met", &m_dPhiLep1Met, "dPhiLep1Met/f");
  TreeSB2MuoTau->Branch("dRLep1Met", &m_dRLep1Met, "dRLep1Met/f");
  TreeSB2MuoTau->Branch("dRLep1Lep2", &m_dRLep1Lep2, "dRLep1Lep2/f");
  TreeSB2MuoTau->Branch("dPhiLep2Met", &m_dPhiLep2Met, "dPhiLep2Met/f");
  TreeSB2MuoTau->Branch("dRLep2Met", &m_dRLep2Met, "dRLep2Met/f");
  TreeSB2MuoTau->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB2MuoTau->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB2MuoTau->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB2MuoTau->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB2MuoTau->Branch("lep1Pt", &m_lep1Pt, "lep1Pt/f");
  TreeSB2MuoTau->Branch("lep1Eta", &m_lep1Eta, "lep1Eta/f");
  TreeSB2MuoTau->Branch("lep1Charge", &m_lep1Charge, "lep1Charge/f");
  TreeSB2MuoTau->Branch("lep1PFIso", &m_lep1PFIso, "lep1PFIso/f");
  TreeSB2MuoTau->Branch("lep1CorrPFIso", &m_lep1CorrPFIso, "lep1CorrPFIso/f");
  TreeSB2MuoTau->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB2MuoTau->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB2MuoTau->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB2MuoTau->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB2MuoTau->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
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
  TreeSB2MuoTau->Branch("EleMuo", &m_EleMuo, "EleMuo/i");
  TreeSB2MuoTau->Branch("MuoMuo", &m_MuoMuo, "MuoMuo/i");
  TreeSB2MuoTau->Branch("EleEle", &m_EleEle, "EleEle/i");
  TreeSB2MuoTau->Branch("MuoTau", &m_MuoTau, "MuoTau/i");
  TreeSB2MuoTau->Branch("EleTau", &m_EleTau, "EleTau/i");

  TreeSB2EleTau = fs->make<TTree>("TreeSB2EleTau", "TreeSB2EleTau");
  TreeSB2EleTau->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB2EleTau->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB2EleTau->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB2EleTau->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB2EleTau->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB2EleTau->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB2EleTau->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeSB2EleTau->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeSB2EleTau->Branch("dPhiLep1Met", &m_dPhiLep1Met, "dPhiLep1Met/f");
  TreeSB2EleTau->Branch("dRLep1Met", &m_dRLep1Met, "dRLep1Met/f");
  TreeSB2EleTau->Branch("dRLep1Lep2", &m_dRLep1Lep2, "dRLep1Lep2/f");
  TreeSB2EleTau->Branch("dPhiLep2Met", &m_dPhiLep2Met, "dPhiLep2Met/f");
  TreeSB2EleTau->Branch("dRLep2Met", &m_dRLep2Met, "dRLep2Met/f");
  TreeSB2EleTau->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB2EleTau->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB2EleTau->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB2EleTau->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB2EleTau->Branch("lep1Pt", &m_lep1Pt, "lep1Pt/f");
  TreeSB2EleTau->Branch("lep1Eta", &m_lep1Eta, "lep1Eta/f");
  TreeSB2EleTau->Branch("lep1Charge", &m_lep1Charge, "lep1Charge/f");
  TreeSB2EleTau->Branch("lep1PFIso", &m_lep1PFIso, "lep1PFIso/f");
  TreeSB2EleTau->Branch("lep1CorrPFIso", &m_lep1CorrPFIso, "lep1CorrPFIso/f");
  TreeSB2EleTau->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB2EleTau->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB2EleTau->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB2EleTau->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB2EleTau->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
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
  TreeSB2EleTau->Branch("EleMuo", &m_EleMuo, "EleMuo/i");
  TreeSB2EleTau->Branch("MuoMuo", &m_MuoMuo, "MuoMuo/i");
  TreeSB2EleTau->Branch("EleEle", &m_EleEle, "EleEle/i");
  TreeSB2EleTau->Branch("MuoTau", &m_MuoTau, "MuoTau/i");
  TreeSB2EleTau->Branch("EleTau", &m_EleTau, "EleTau/i");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Analyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
Analyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
Analyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Analyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Analyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void Analyzer::SelectJet(edm::Handle<pat::JetCollection> CA8JetswithQjets,
			 edm::Handle<pat::JetCollection> CA8JetsPruned,
			 bool & foundJet,
			 pat::JetCollection::const_iterator & SelectedJet,
			 float & prunedMass, float & tau21Z, float & ptZ,
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
      prunedMass=mass;
      ptZ=jet->pt();
      tau21Z=jet->userFloat("tau2")/jet->userFloat("tau1");
      SelectedJet=jet;
    }
  }
}

void Analyzer::SelectTau(edm::Handle<pat::TauCollection> tauHandle,
			 pat::JetCollection::const_iterator SelectedJet,
			 bool & foundTau,
			 pat::TauCollection::const_iterator & SelectedTau,
			 float & ptTau, bool foundJet){
  
  for (pat::TauCollection::const_iterator patTau = tauHandle->begin(); patTau != tauHandle->end(); ++patTau ) {
    if(patTau->pt()<20) continue;
    if(abs(patTau->eta())>2.4) continue;
    if(patTau->tauID("decayModeFindingNewDMs")<0.5) continue;
    if(patTau->tauID("againstMuonLoose")<0.5) continue;
    if(patTau->tauID("againstElectronLoose")<0.5) continue;
    if(patTau->tauID("byVLooseIsolationMVA3newDMwoLT")<0.5) continue;
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

void Analyzer::SelectMuon(edm::Handle<pat::MuonCollection> muoH,
			  pat::JetCollection::const_iterator SelectedJet, bool & foundMuon,
			  pat::MuonCollection::const_iterator & SelectedMuon,
			  float & ptMuon, bool foundJet, reco::Vertex primaryVertex, bool fully){
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
    if(fully){if(MuonPFIso(muon, true)>0.2) continue;}
    //else     {if(MuonCorrPFIso(muon, true)>0.2) continue;}
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

void Analyzer::SelectElectron(edm::Handle<pat::ElectronCollection> eleH,
			      pat::JetCollection::const_iterator SelectedJet, bool & foundElectron,
			      pat::ElectronCollection::const_iterator & SelectedElectron,
			      float & ptElectron, bool foundJet, reco::Vertex primaryVertex, 
			      float lep2Pt, float rho, bool fully){
  bool passEle=false;
  for(pat::ElectronCollection::const_iterator electron = eleH->begin(); electron != eleH->end(); ++electron) {
    if(electron->pt()<10) continue;
    if(electron->pt()==lep2Pt) continue;
    if(fully){if(ElectronPFIso(electron,rho)>0.1) continue;}
    //else     {if(ElectronCorrPFIso(electron,rho)>0.1) continue;}
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
      passEle=true;
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
      passEle=true;
    }
    if(foundJet && passEle){
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

void Analyzer::SelectMuonMM(pat::MuonCollection::const_iterator & SelectedMuon1, 
			    pat::MuonCollection::const_iterator & SelectedMuon2, 
			    bool & foundMuonMM, std::vector<pat::MuonCollection::const_iterator> SelectedTrackerMuo, 
			    pat::JetCollection::const_iterator SelectedJet, bool foundJet, reco::Vertex primaryVertex){
  bool hasAtLeastOneHighPtMuo=false; float pt1=-99; float pt2=-99;
  vector<pat::MuonCollection::const_iterator> SelectedMuo;
  pat::MuonCollection::const_iterator Muon1;
  pat::MuonCollection::const_iterator Muon2;
  for(unsigned int i=0; i<SelectedTrackerMuo.size(); i++){
    if(!foundJet) continue;
    if((MuonDETIso(SelectedTrackerMuo[i],SelectedTrackerMuo))>0.2) continue;
    if(ROOT::Math::VectorUtil::DeltaR(SelectedTrackerMuo[i]->p4(),SelectedJet->p4())<0.8) continue;
    SelectedMuo.push_back(SelectedTrackerMuo[i]);
  }
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
    if(SelectedMuo[i]->pt()>pt1){
      Muon1=SelectedMuo[i];
      pt1=SelectedMuo[i]->pt();
    }
  }
  if(hasAtLeastOneHighPtMuo){
    for(unsigned int i=0; i<SelectedMuo.size(); i++){
      if(SelectedMuo[i]==Muon1) continue;
      foundMuonMM=true;
      if(SelectedMuo[i]->pt()>pt2){
	Muon2=SelectedMuo[i];
	pt2=SelectedMuo[i]->pt();
      }
    }
  }
  if(foundMuonMM){
    if(Muon1->pt()>Muon2->pt()){
      SelectedMuon1=Muon1;
      SelectedMuon2=Muon2;
    }else{
      SelectedMuon1=Muon2;
      SelectedMuon2=Muon1;
    }
  }
}


void Analyzer::SelectTrackerMuon(edm::Handle<pat::MuonCollection> muoH, vector<pat::MuonCollection::const_iterator> & SelectedMuo, reco::Vertex primaryVertex){
  for(pat::MuonCollection::const_iterator muon = muoH->begin(); muon != muoH->end(); ++muon) {
    if(muon->pt()<10) continue;
    if(abs(muon->eta())>2.4) continue;
    if(abs(muon->phi())>3.2) continue;
    if(!(muon->isTrackerMuon())) continue;
    if(muon->numberOfMatches()<=1) continue;
    if(fabs(muon->muonBestTrack()->dz(primaryVertex.position()))>=0.5) continue;
    if(fabs(muon->dB())>=0.2 ) continue;
    if(muon->innerTrack()->hitPattern().numberOfValidPixelHits()<=0) continue;
    if(muon->innerTrack()->hitPattern().trackerLayersWithMeasurement()<=5) continue;
    if((muon->muonBestTrack()->ptError()/muon->muonBestTrack()->pt())>0.3) continue;
    SelectedMuo.push_back(muon);
  }
}

float Analyzer::MuonPFIso(pat::MuonCollection::const_iterator muon, bool highpt){
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


float Analyzer::ElectronPFIso(pat::ElectronCollection::const_iterator electron, float rho){
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

float Analyzer::MuonCorrPFIso(pat::MuonCollection::const_iterator muon, bool highpt){
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


float Analyzer::ElectronCorrPFIso(pat::ElectronCollection::const_iterator electron, float rho){
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

float Analyzer::MuonDETIso(pat::MuonCollection::const_iterator SelectedMuo, vector<pat::MuonCollection::const_iterator> SelectedHighptMuo){
  float isovar = SelectedMuo->trackIso()/SelectedMuo->pt();
  for(unsigned int j = 0; j< SelectedHighptMuo.size();++j){
    if(SelectedMuo==SelectedHighptMuo[j]) continue;
    double dR = ROOT::Math::VectorUtil::DeltaR(SelectedMuo->p4(),SelectedHighptMuo[j]->p4());
    if(dR < 0.3) isovar = isovar - ((SelectedHighptMuo[j]->track()->pt())/SelectedMuo->pt());
  }
  return isovar;
}

void Analyzer::BtagVeto(edm::Handle<pat::JetCollection> ak5jetCands, int & nbtagsL, int & nbtagsM, int & nbtagsT,
			pat::JetCollection::const_iterator SelectedJet){
  for(pat::JetCollection::const_iterator ak5 = ak5jetCands->begin(); ak5 != ak5jetCands->end(); ++ak5) {
    if(ROOT::Math::VectorUtil::DeltaR(ak5->p4(),SelectedJet->p4())<0.8) continue;
    double discCSV = ak5->bDiscriminator("combinedSecondaryVertexBJetTags");
    if(discCSV>0.244) nbtagsL++; //loose working point
    if(discCSV>0.679) nbtagsM++; //medium working point
    if(discCSV>0.898) nbtagsT++; //tight working point
  }
}

void Analyzer::Efficiency(float & genEvent, bool isData, Handle<vector<reco::GenParticle> > genParts){
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

void Analyzer::svfit(edm::Handle<pat::METCollection> metRaw, edm::Handle<pat::METCollection> met, LorentzVector SelectedLep1, 
		     LorentzVector SelectedLep2, TLorentzVector PrunedJet, float & MassSVFit, float & XMassSVFit, float & dRJetZSVFit, float & ptSVFit, 
		     pat::JetCollection::const_iterator SelectedJet, int category){
  TMatrixD covMET(2, 2); // PFMET significance matrix
  covMET[0][0] = (metRaw->front() ).getSignificanceMatrix()(0,0);
  covMET[1][0] = (metRaw->front() ).getSignificanceMatrix()(1,0);
  covMET[0][1] = (metRaw->front() ).getSignificanceMatrix()(0,1);
  covMET[1][1] = (metRaw->front() ).getSignificanceMatrix()(1,1);
  std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
  if(category==0 || category==1 || category==2){
    measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, SelectedLep1));
    measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, SelectedLep2));
  }
  if(category==3 || category==4){
    measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay, SelectedLep1));
    measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, SelectedLep2));
  }
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

void Analyzer::FillTree(int category, TTree *Tree, pat::JetCollection::const_iterator SelectedJet, pat::TauCollection::const_iterator SelectedTau,
			pat::MuonCollection::const_iterator SelectedMuo, pat::MuonCollection::const_iterator SelectedMuo1, pat::MuonCollection::const_iterator SelectedMuo2,
			pat::ElectronCollection::const_iterator SelectedEle, pat::ElectronCollection::const_iterator SelectedEle1, 
			pat::ElectronCollection::const_iterator SelectedEle2, edm::Handle<reco::VertexCollection> vertices,
			edm::Handle<pat::METCollection> metRaw, edm::Handle<pat::METCollection> met, edm::Handle<pat::METCollection> uncorrmet,
			float prunedMass, int nbtagsL, int nbtagsM, int nbtagsT, bool isFired_HLT, bool isFired_HLT_PFJet320, bool isFired_HLT_HT650,
			double MyWeight, float genEvent, float rho, int EleMuo, int MuoMuo, int EleEle, int MuoTau, int EleTau){
  
  math::PtEtaPhiELorentzVector lep1;
  math::PtEtaPhiELorentzVector lep2;
  LorentzVector lep1SVFit;
  LorentzVector lep2SVFit;
  float lep1Charge=0;  float lep2Charge=0;
  float lep1PFIso=100; float lep1CorrPFIso=100;
  float lep2PFIso=100; float lep2CorrPFIso=100;
  if(category==0){
    lep1 = SelectedEle->p4();
    lep2 = SelectedMuo->p4();
    lep1SVFit = SelectedEle->p4();
    lep2SVFit = SelectedMuo->p4();
    lep1Charge=SelectedEle->charge();
    lep2Charge=SelectedMuo->charge();
    lep1PFIso=ElectronPFIso(SelectedEle,rho);
    lep1CorrPFIso=ElectronCorrPFIso(SelectedEle,rho);
    lep2PFIso=MuonPFIso(SelectedMuo,true);
    lep2CorrPFIso=MuonCorrPFIso(SelectedMuo,true);
  }else if(category==1){
    lep1 = SelectedMuo1->p4();
    lep2 = SelectedMuo2->p4();
    lep1SVFit = SelectedMuo1->p4();
    lep2SVFit = SelectedMuo2->p4();
    lep1Charge=SelectedMuo1->charge();
    lep2Charge=SelectedMuo2->charge();
    lep1PFIso=MuonPFIso(SelectedMuo1,true);
    lep1CorrPFIso=MuonCorrPFIso(SelectedMuo1,true);
    lep2PFIso=MuonPFIso(SelectedMuo2,true);
    lep2CorrPFIso=MuonCorrPFIso(SelectedMuo2,true);
  }else if(category==2){
    lep1 = SelectedEle1->p4();
    lep2 = SelectedEle2->p4();
    lep1SVFit = SelectedEle1->p4();
    lep2SVFit = SelectedEle2->p4();
    lep1Charge=SelectedEle1->charge();
    lep2Charge=SelectedEle2->charge();
    lep1PFIso=ElectronPFIso(SelectedEle1,rho);
    lep1CorrPFIso=ElectronCorrPFIso(SelectedEle1,rho);
    lep2PFIso=ElectronPFIso(SelectedEle2,rho);
    lep2CorrPFIso=ElectronCorrPFIso(SelectedEle2,rho);
  }else if(category==3){
    lep1 = SelectedTau->p4();
    lep2 = SelectedMuo->p4();
    lep1SVFit = SelectedTau->p4();
    lep2SVFit = SelectedMuo->p4();
    lep1Charge=SelectedTau->charge();
    lep2Charge=SelectedMuo->charge();
    lep1PFIso=-10;
    lep1CorrPFIso=-10;
    lep2PFIso=MuonPFIso(SelectedMuo,true);
    lep2CorrPFIso=MuonCorrPFIso(SelectedMuo,true);
  }else if(category==4){
    lep1 = SelectedTau->p4();
    lep2 = SelectedEle->p4();
    lep1SVFit = SelectedTau->p4();
    lep2SVFit = SelectedEle->p4();
    lep1Charge=SelectedTau->charge();
    lep2Charge=SelectedEle->charge();
    lep1PFIso=-10;
    lep1CorrPFIso=-10;
    lep2PFIso=ElectronPFIso(SelectedEle,rho);
    lep2CorrPFIso=ElectronCorrPFIso(SelectedEle,rho);
  }

  //SVFIT
  math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedJet->pt(),SelectedJet->eta(),SelectedJet->phi(),prunedMass);
  TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());
  float MassSVFit = -1;
  float XMassSVFit = -1;
  float dRJetZSVFit = -1;
  float ptSVFit = -1;
  svfit(metRaw, met, lep1SVFit, lep2SVFit, PrunedJet, MassSVFit, XMassSVFit, dRJetZSVFit, ptSVFit, SelectedJet, category);
  //COLLINEAR APPROXIMATION
  TLorentzVector CATauTau; bool CA = false;
  float a = (lep1.py()*met->begin()->px()-lep1.px()*met->begin()->py())/
    (lep2.px()*lep1.py()-lep2.py()*lep1.px());
  float b = (lep2.py()*met->begin()->px()-lep2.px()*met->begin()->py())/
    (lep1.px()*lep2.py()-lep1.py()*lep2.px());
  if(((1+a)*(1+b))>0) {
    CATauTau.SetPxPyPzE((1+a)*lep2.px()+(1+b)*lep1.px(),
			(1+a)*lep2.py()+(1+b)*lep1.py(),
			(1+a)*lep2.pz()+(1+b)*lep1.pz(),
			(1+a)*lep2.energy()+(1+b)*lep1.energy());
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
  math::PtEtaPhiELorentzVector dilep; dilep = lep1+lep2;
  math::PtEtaPhiELorentzVector dilepmet; dilepmet = lep1+lep2+met->begin()->p4();
  math::PtEtaPhiELorentzVector dilepjet; dilepjet = lep1+lep2+PrunedJet_prov;
  math::PtEtaPhiELorentzVector dilepmetjet; dilepmetjet = lep1+lep2+met->begin()->p4()+PrunedJet_prov;
  m_jetPt=SelectedJet->pt();
  m_jetEta=SelectedJet->eta();
  m_jetMass=prunedMass;
  m_jetSubjettiness=SelectedJet->userFloat("tau2")/SelectedJet->userFloat("tau1");
  m_lep1Pt=lep1.pt();
  m_lep1Eta=lep1.eta();
  m_lep1Charge=lep1Charge;
  m_lep1PFIso=lep1PFIso;
  m_lep1CorrPFIso=lep1CorrPFIso;
  m_lep2Pt=lep2.pt();
  m_lep2Eta=lep2.eta();
  m_lep2Charge=lep2Charge;
  m_lep2PFIso=lep2PFIso;
  m_lep2CorrPFIso=lep2CorrPFIso;
  m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet->p4());
  m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedJet->p4());
  m_dRJetLep2=ROOT::Math::VectorUtil::DeltaR(lep2,SelectedJet->p4());
  m_dRJetLep1=ROOT::Math::VectorUtil::DeltaR(lep1,SelectedJet->p4());
  m_dPhiLep1Met=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),lep1);
  m_dRLep1Met=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),lep1);
  m_dRLep1Lep2=ROOT::Math::VectorUtil::DeltaR(lep2,lep1);
  m_dPhiLep2Met=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),lep2);
  m_dRLep2Met=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),lep2);
  m_dRZZVis=ROOT::Math::VectorUtil::DeltaR(dilep,PrunedJet);
  m_dRZZEff=ROOT::Math::VectorUtil::DeltaR(dilepmet,PrunedJet);
  m_dRZZSvFit=dRJetZSVFit;
  m_dRZZCA=dRJetZCA;
  m_charge=lep1Charge*lep2Charge;
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
  m_sideband=(int)(prunedMass<70);
  m_NVertices=vertices->size();
  m_PUWeight=MyWeight;
  m_metPx=met->begin()->px();
  m_metPy=met->begin()->py();
  m_NeventsTOT=NeventsTOT_;
  m_xsec=xsec_;
  m_lumi=lumi_;
  m_weight=xsec_*lumi_/NeventsTOT_;
  m_genEvent = genEvent;
  m_EleMuo=EleMuo;
  m_MuoMuo=MuoMuo;
  m_EleEle=EleEle;
  m_MuoTau=MuoTau;
  m_EleTau=EleTau;
  Tree->Fill();
}

//define this as a plug-in
DEFINE_FWK_MODULE(Analyzer);
