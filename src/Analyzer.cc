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
#include "CMGTools/External/interface/PileupJetIdentifier.h"

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
                 pat::JetCollection::const_iterator & SelectedJet, float & prunedMass, float massMin, float massMax, bool tau21DOWN);
  void SelectTau(edm::Handle<pat::TauCollection> tauHandle, pat::JetCollection::const_iterator SelectedJet, 
		 std::vector<pat::TauCollection::const_iterator> & SelectedTau, bool foundJet);
  void SelectMuon(edm::Handle<pat::MuonCollection> muoH, pat::JetCollection::const_iterator SelectedJet,
		  std::vector<pat::MuonCollection::const_iterator> & SelectedMuon, bool foundJet, reco::Vertex primaryVertex, bool fully);
  void SelectElectron(edm::Handle<pat::ElectronCollection> eleH, pat::JetCollection::const_iterator SelectedJet,
		      std::vector<pat::ElectronCollection::const_iterator> & SelectedElectron, bool foundJet, reco::Vertex primaryVertex, float rho, bool fully);
  void SelectTrackerMuon(edm::Handle<pat::MuonCollection> muoH, std::vector<pat::MuonCollection::const_iterator> & SelectedMuo, reco::Vertex primaryVertex);



  void SelectEM(pat::ElectronCollection::const_iterator & SelectedElectron, pat::MuonCollection::const_iterator & SelectedMuon, bool & foundEM,
		std::vector<pat::ElectronCollection::const_iterator> SelectedElectrons, std::vector<pat::MuonCollection::const_iterator> SelectedMuons);
  void SelectMM(pat::MuonCollection::const_iterator & SelectedMuon1, pat::MuonCollection::const_iterator & SelectedMuon2, bool & foundMM, 
		std::vector<pat::MuonCollection::const_iterator> SelectedTrackerMuo, pat::JetCollection::const_iterator SelectedJet, bool foundJet, 
		reco::Vertex primaryVertex);
  void SelectEE(pat::ElectronCollection::const_iterator & SelectedElectron1, pat::ElectronCollection::const_iterator & SelectedElectron2, bool & foundEE,
		std::vector<pat::ElectronCollection::const_iterator> SelectedElectrons);
  void SelectMT(pat::MuonCollection::const_iterator & SelectedMuon, pat::TauCollection::const_iterator & SelectedTau, bool & foundMT,
		std::vector<pat::MuonCollection::const_iterator> SelectedMuons, std::vector<pat::TauCollection::const_iterator> SelectedTaus);
  void SelectET(pat::ElectronCollection::const_iterator & SelectedElectron, pat::TauCollection::const_iterator & SelectedTau, bool & foundET,
		std::vector<pat::ElectronCollection::const_iterator> SelectedElectrons, std::vector<pat::TauCollection::const_iterator> SelectedTaus);



  void Efficiency(float & genEvent, bool isData, edm::Handle<std::vector<reco::GenParticle> > genParts);
  void svfit(edm::Handle<pat::METCollection> metRaw, edm::Handle<pat::METCollection> met, LorentzVector SelectedTau, LorentzVector SelectedMuon, TLorentzVector PrunedJet,
	     float & MassSVFit, float & XMassSVFit, float & dRJetZSVFit, float & ptSVFit, pat::JetCollection::const_iterator SelectedJet, int category);
  
  float MuonPFIso(pat::MuonCollection::const_iterator muon, bool highpt);
  float ElectronPFIso(pat::ElectronCollection::const_iterator electron, float rho);
  float MuonCorrPFIso(pat::MuonCollection::const_iterator muon, bool highpt);
  float ElectronCorrPFIso(pat::ElectronCollection::const_iterator electron, float rho);
  float MuonDETIso(pat::MuonCollection::const_iterator SelectedMuo, std::vector<pat::MuonCollection::const_iterator> SelectedHighptMuo);
  void BtagVeto(int & njet1, int & nbtagsL1, int & nbtagsM1, int & nbtagsT1, int & njet2, int & nbtagsL2, int & nbtagsM2, int & nbtagsT2,
		int & njet3, int & nbtagsL3, int & nbtagsM3, int & nbtagsT3, int & njet4, int & nbtagsL4, int & nbtagsM4, int & nbtagsT4,
		pat::JetCollection::const_iterator SelectedJet, math::PtEtaPhiELorentzVector lep1, math::PtEtaPhiELorentzVector lep2, const edm::Event& iEvent);
  void FillTree(int category, TTree *Tree, pat::JetCollection::const_iterator SelectedJet, pat::TauCollection::const_iterator SelectedTau,
		pat::MuonCollection::const_iterator SelectedMuo, pat::MuonCollection::const_iterator SelectedMuo1, pat::MuonCollection::const_iterator SelectedMuo2,
		std::vector<pat::MuonCollection::const_iterator> SelectedTrackerMuo,
		pat::ElectronCollection::const_iterator SelectedEle, pat::ElectronCollection::const_iterator SelectedEle1, 
		pat::ElectronCollection::const_iterator SelectedEle2, edm::Handle<reco::VertexCollection> vertices,
		edm::Handle<pat::METCollection> metRaw, edm::Handle<pat::METCollection> met, edm::Handle<pat::METCollection> uncorrmet,
		float prunedMass, bool isFired_HLT, bool isFired_HLT_PFJet320, bool isFired_HLT_HT650,
		double MyWeight, float genEvent, float rho, int EleMuo, int MuoMuo, int EleEle, int MuoTau, int EleTau, const edm::Event& iEvent);
  
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
  TTree *TreeSB3EleMuo;
  TTree *TreeSB3MuoMuo;
  TTree *TreeSB3EleEle;
  TTree *TreeSB3MuoTau;
  TTree *TreeSB3EleTau;
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
  float m_lep1DETIso;
  float m_lep2Pt;
  float m_lep2Eta;
  float m_lep2Charge;
  float m_lep2PFIso;
  float m_lep2CorrPFIso;
  float m_lep2DETIso;
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
  int m_njet1;
  int m_nbtagsL1;
  int m_nbtagsM1;
  int m_nbtagsT1;
  int m_njet2;
  int m_nbtagsL2;
  int m_nbtagsM2;
  int m_nbtagsT2;
  int m_njet3;
  int m_nbtagsL3;
  int m_nbtagsM3;
  int m_nbtagsT3;
  int m_njet4;
  int m_nbtagsL4;
  int m_nbtagsM4;
  int m_nbtagsT4;
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
  TFile *da_=new TFile ("/data06/users/spiezia/EXO/CMSSW_5_3_13/src/Analyzer/EDBRTauAnalyzer/data/MyDataPileupHistogram_True.root");
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

    
  //--------------------------------------------- SELECT OBJECTS ---------------------------------------------//

  //JET SELECTION
  pat::JetCollection::const_iterator SelectedJet;
  pat::JetCollection::const_iterator SelectedSB1Jet;
  pat::JetCollection::const_iterator SelectedSB2Jet;
  pat::JetCollection::const_iterator SelectedSB3Jet;
  float prunedMass   =-9999; bool foundJet   =false;   
  float prunedMassSB1=-9999; bool foundSB1Jet=false;
  float prunedMassSB2=-9999; bool foundSB2Jet=false;
  float prunedMassSB3=-9999; bool foundSB3Jet=false;
  SelectJet(CA8JetswithQjets, CA8JetsPruned, foundJet,    SelectedJet,    prunedMass,    70,  110,   true);
  SelectJet(CA8JetswithQjets, CA8JetsPruned, foundSB1Jet, SelectedSB1Jet, prunedMassSB1, 20,  70,    true);
  SelectJet(CA8JetswithQjets, CA8JetsPruned, foundSB2Jet, SelectedSB2Jet, prunedMassSB2, 110, 99999, true);
  SelectJet(CA8JetswithQjets, CA8JetsPruned, foundSB3Jet, SelectedSB3Jet, prunedMassSB3, 20,  99999, false);

  //TAU SELECTION
  std::vector<pat::TauCollection::const_iterator> SelectedTausMT;
  std::vector<pat::TauCollection::const_iterator> SelectedSB1TausMT;
  std::vector<pat::TauCollection::const_iterator> SelectedSB2TausMT;
  std::vector<pat::TauCollection::const_iterator> SelectedSB3TausMT;
  std::vector<pat::TauCollection::const_iterator> SelectedTausET;
  std::vector<pat::TauCollection::const_iterator> SelectedSB1TausET;
  std::vector<pat::TauCollection::const_iterator> SelectedSB2TausET;
  std::vector<pat::TauCollection::const_iterator> SelectedSB3TausET;
  SelectTau(tauMuTauHandle, SelectedJet,    SelectedTausMT,    foundJet);
  SelectTau(tauMuTauHandle, SelectedSB1Jet, SelectedSB1TausMT, foundSB1Jet);
  SelectTau(tauMuTauHandle, SelectedSB2Jet, SelectedSB2TausMT, foundSB2Jet); 
  SelectTau(tauMuTauHandle, SelectedSB3Jet, SelectedSB3TausMT, foundSB3Jet); 
  SelectTau(tauElTauHandle, SelectedJet,    SelectedTausET,    foundJet);
  SelectTau(tauElTauHandle, SelectedSB1Jet, SelectedSB1TausET, foundSB1Jet);
  SelectTau(tauElTauHandle, SelectedSB2Jet, SelectedSB2TausET, foundSB2Jet); 
  SelectTau(tauElTauHandle, SelectedSB3Jet, SelectedSB3TausET, foundSB3Jet); 

  //MUON SELECTION
  std::vector<pat::MuonCollection::const_iterator> SelectedMuonsEM;
  std::vector<pat::MuonCollection::const_iterator> SelectedSB1MuonsEM;
  std::vector<pat::MuonCollection::const_iterator> SelectedSB2MuonsEM;
  std::vector<pat::MuonCollection::const_iterator> SelectedSB3MuonsEM;
  std::vector<pat::MuonCollection::const_iterator> SelectedMuonsMT;
  std::vector<pat::MuonCollection::const_iterator> SelectedSB1MuonsMT;
  std::vector<pat::MuonCollection::const_iterator> SelectedSB2MuonsMT;
  std::vector<pat::MuonCollection::const_iterator> SelectedSB3MuonsMT;
  SelectMuon(muoH, SelectedJet,    SelectedMuonsEM,    foundJet,    primaryVertex, true);
  SelectMuon(muoH, SelectedSB1Jet, SelectedSB1MuonsEM, foundSB1Jet, primaryVertex, true);
  SelectMuon(muoH, SelectedSB2Jet, SelectedSB2MuonsEM, foundSB2Jet, primaryVertex, true);
  SelectMuon(muoH, SelectedSB3Jet, SelectedSB3MuonsEM, foundSB3Jet, primaryVertex, true);
  SelectMuon(muoH, SelectedJet,    SelectedMuonsMT,    foundJet,    primaryVertex, false);
  SelectMuon(muoH, SelectedSB1Jet, SelectedSB1MuonsMT, foundSB1Jet, primaryVertex, false);
  SelectMuon(muoH, SelectedSB2Jet, SelectedSB2MuonsMT, foundSB2Jet, primaryVertex, false);
  SelectMuon(muoH, SelectedSB3Jet, SelectedSB3MuonsMT, foundSB3Jet, primaryVertex, false);

  //ELECTRON SELECTION
  std::vector<pat::ElectronCollection::const_iterator> SelectedElectronsEEandEM;
  std::vector<pat::ElectronCollection::const_iterator> SelectedSB1ElectronsEEandEM;
  std::vector<pat::ElectronCollection::const_iterator> SelectedSB2ElectronsEEandEM;
  std::vector<pat::ElectronCollection::const_iterator> SelectedSB3ElectronsEEandEM;
  std::vector<pat::ElectronCollection::const_iterator> SelectedElectronsET;
  std::vector<pat::ElectronCollection::const_iterator> SelectedSB1ElectronsET;
  std::vector<pat::ElectronCollection::const_iterator> SelectedSB2ElectronsET;
  std::vector<pat::ElectronCollection::const_iterator> SelectedSB3ElectronsET;
  SelectElectron(eleH,   SelectedJet,    SelectedElectronsEEandEM,    foundJet,    primaryVertex, rho, true);
  SelectElectron(eleH,   SelectedSB1Jet, SelectedSB1ElectronsEEandEM, foundSB1Jet, primaryVertex, rho, true);
  SelectElectron(eleH,   SelectedSB2Jet, SelectedSB2ElectronsEEandEM, foundSB2Jet, primaryVertex, rho, true);
  SelectElectron(eleH,   SelectedSB3Jet, SelectedSB3ElectronsEEandEM, foundSB3Jet, primaryVertex, rho, true);
  SelectElectron(eleETH, SelectedJet,    SelectedElectronsET,         foundJet,    primaryVertex, rho, false);
  SelectElectron(eleETH, SelectedSB1Jet, SelectedSB1ElectronsET,      foundSB1Jet, primaryVertex, rho, false);
  SelectElectron(eleETH, SelectedSB2Jet, SelectedSB2ElectronsET,      foundSB2Jet, primaryVertex, rho, false);
  SelectElectron(eleETH, SelectedSB3Jet, SelectedSB3ElectronsET,      foundSB3Jet, primaryVertex, rho, false);


  //--------------------------------------------- SELECT EVENTS ---------------------------------------------//
  //SELECT ELEMUO EVENTS - SR
  pat::ElectronCollection::const_iterator SelectedElectronEM;
  pat::MuonCollection::const_iterator SelectedMuonEM;
  bool foundEM=false;
  SelectEM(SelectedElectronEM, SelectedMuonEM, foundEM, SelectedElectronsEEandEM, SelectedMuonsEM);
  //SELECT ELEMUO EVENTS - SB1
  pat::ElectronCollection::const_iterator SelectedSB1ElectronEM;
  pat::MuonCollection::const_iterator SelectedSB1MuonEM;
  bool foundSB1EM=false;
  SelectEM(SelectedSB1ElectronEM, SelectedSB1MuonEM, foundSB1EM, SelectedSB1ElectronsEEandEM, SelectedSB1MuonsEM);
  //SELECT ELEMUO EVENTS - SB2
  pat::ElectronCollection::const_iterator SelectedSB2ElectronEM;
  pat::MuonCollection::const_iterator SelectedSB2MuonEM;
  bool foundSB2EM=false;
  SelectEM(SelectedSB2ElectronEM, SelectedSB2MuonEM, foundSB2EM, SelectedSB2ElectronsEEandEM, SelectedSB2MuonsEM);
  //SELECT ELEMUO EVENTS - SB3
  pat::ElectronCollection::const_iterator SelectedSB3ElectronEM;
  pat::MuonCollection::const_iterator SelectedSB3MuonEM;
  bool foundSB3EM=false;
  SelectEM(SelectedSB3ElectronEM, SelectedSB3MuonEM, foundSB3EM, SelectedSB3ElectronsEEandEM, SelectedSB3MuonsEM);

  //SELECT MUOMUO EVENTS - SR
  vector<pat::MuonCollection::const_iterator> SelectedTrackerMuo;
  SelectTrackerMuon(muoH, SelectedTrackerMuo, primaryVertex);
  bool foundMM=false; 
  pat::MuonCollection::const_iterator SelectedMuon1MM;
  pat::MuonCollection::const_iterator SelectedMuon2MM;
  SelectMM(SelectedMuon1MM, SelectedMuon2MM, foundMM, SelectedTrackerMuo, SelectedJet, foundJet, primaryVertex);
  //SELECT MUOMUO EVENTS - SB1
  bool foundSB1MM=false;
  pat::MuonCollection::const_iterator SelectedSB1Muon1MM;
  pat::MuonCollection::const_iterator SelectedSB1Muon2MM;
  SelectMM(SelectedSB1Muon1MM, SelectedSB1Muon2MM, foundSB1MM, SelectedTrackerMuo, SelectedSB1Jet, foundSB1Jet, primaryVertex);
  //SELECT MUOMUO EVENTS - SB2
  bool foundSB2MM=false;
  pat::MuonCollection::const_iterator SelectedSB2Muon1MM;
  pat::MuonCollection::const_iterator SelectedSB2Muon2MM;
  SelectMM(SelectedSB2Muon1MM, SelectedSB2Muon2MM, foundSB2MM, SelectedTrackerMuo, SelectedSB2Jet, foundSB2Jet, primaryVertex);
  //SELECT MUOMUO EVENTS - SB3
  bool foundSB3MM=false;
  pat::MuonCollection::const_iterator SelectedSB3Muon1MM;
  pat::MuonCollection::const_iterator SelectedSB3Muon2MM;
  SelectMM(SelectedSB3Muon1MM, SelectedSB3Muon2MM, foundSB3MM, SelectedTrackerMuo, SelectedSB3Jet, foundSB3Jet, primaryVertex);

  //SELECT ELEELE EVENTS - SR
  pat::ElectronCollection::const_iterator SelectedElectron1EE;
  pat::ElectronCollection::const_iterator SelectedElectron2EE;
  bool foundEE=false;
  SelectEE(SelectedElectron1EE, SelectedElectron2EE, foundEE, SelectedElectronsEEandEM);
  //SELECT ELEELE EVENTS - SB1
  pat::ElectronCollection::const_iterator SelectedSB1Electron1EE;
  pat::ElectronCollection::const_iterator SelectedSB1Electron2EE;
  bool foundSB1EE=false;
  SelectEE(SelectedSB1Electron1EE, SelectedSB1Electron2EE, foundSB1EE, SelectedSB1ElectronsEEandEM);
  //SELECT ELEELE EVENTS - SB2
  pat::ElectronCollection::const_iterator SelectedSB2Electron1EE;
  pat::ElectronCollection::const_iterator SelectedSB2Electron2EE;
  bool foundSB2EE=false;
  SelectEE(SelectedSB2Electron1EE, SelectedSB2Electron2EE, foundSB2EE, SelectedSB2ElectronsEEandEM);
  //SELECT ELEELE EVENTS - SB3
  pat::ElectronCollection::const_iterator SelectedSB3Electron1EE;
  pat::ElectronCollection::const_iterator SelectedSB3Electron2EE;
  bool foundSB3EE=false;
  SelectEE(SelectedSB3Electron1EE, SelectedSB3Electron2EE, foundSB3EE, SelectedSB3ElectronsEEandEM);

  //SELECT MUOTAU EVENTS - SR
  pat::MuonCollection::const_iterator SelectedMuonMT;
  pat::TauCollection::const_iterator SelectedTauMT;
  bool foundMT=false;
  SelectMT(SelectedMuonMT, SelectedTauMT, foundMT, SelectedMuonsMT, SelectedTausMT);
  //SELECT MUOTAU EVENTS - SB1
  pat::MuonCollection::const_iterator SelectedSB1MuonMT;
  pat::TauCollection::const_iterator SelectedSB1TauMT;
  bool foundSB1MT=false;
  SelectMT(SelectedSB1MuonMT, SelectedSB1TauMT, foundSB1MT, SelectedSB1MuonsMT, SelectedSB1TausMT);
  //SELECT MUOTAU EVENTS - SB2
  pat::MuonCollection::const_iterator SelectedSB2MuonMT;
  pat::TauCollection::const_iterator SelectedSB2TauMT;
  bool foundSB2MT=false;
  SelectMT(SelectedSB2MuonMT, SelectedSB2TauMT, foundSB2MT, SelectedSB2MuonsMT, SelectedSB2TausMT);
  //SELECT MUOTAU EVENTS - SB3
  pat::MuonCollection::const_iterator SelectedSB3MuonMT;
  pat::TauCollection::const_iterator SelectedSB3TauMT;
  bool foundSB3MT=false;
  SelectMT(SelectedSB3MuonMT, SelectedSB3TauMT, foundSB3MT, SelectedSB3MuonsMT, SelectedSB3TausMT);

  //SELECT ELETAU EVENTS - SR
  pat::ElectronCollection::const_iterator SelectedElectronET;
  pat::TauCollection::const_iterator SelectedTauET;
  bool foundET=false;
  SelectET(SelectedElectronET, SelectedTauET, foundET, SelectedElectronsET, SelectedTausET);
  //SELECT ELETAU EVENTS - SB1
  pat::ElectronCollection::const_iterator SelectedSB1ElectronET;
  pat::TauCollection::const_iterator SelectedSB1TauET;
  bool foundSB1ET=false;
  SelectET(SelectedSB1ElectronET, SelectedSB1TauET, foundSB1ET, SelectedSB1ElectronsET, SelectedSB1TausET);
  //SELECT ELETAU EVENTS - SB2
  pat::ElectronCollection::const_iterator SelectedSB2ElectronET;
  pat::TauCollection::const_iterator SelectedSB2TauET;
  bool foundSB2ET=false;
  SelectET(SelectedSB2ElectronET, SelectedSB2TauET, foundSB2ET, SelectedSB2ElectronsET, SelectedSB2TausET);
  //SELECT ELETAU EVENTS - SB3
  pat::ElectronCollection::const_iterator SelectedSB3ElectronET;
  pat::TauCollection::const_iterator SelectedSB3TauET;
  bool foundSB3ET=false;
  SelectET(SelectedSB3ElectronET, SelectedSB3TauET, foundSB3ET, SelectedSB3ElectronsET, SelectedSB3TausET);


  //--------------------------------------------- FILL TREES ---------------------------------------------//
  pat::TauCollection::const_iterator      SelectedTauFake;
  pat::MuonCollection::const_iterator     SelectedMuonFake;
  pat::MuonCollection::const_iterator     SelectedMuon1Fake;
  pat::MuonCollection::const_iterator     SelectedMuon2Fake;
  pat::ElectronCollection::const_iterator SelectedElectronFake;
  pat::ElectronCollection::const_iterator SelectedElectron1Fake;
  pat::ElectronCollection::const_iterator SelectedElectron2Fake;

  int EleMuo = 0; int MuoMuo = 0; int EleEle = 0; int MuoTau = 0; int EleTau = 0;
  if(foundJet && foundEM) EleMuo=1;
  if(foundJet && foundMM) MuoMuo=1;
  if(foundJet && foundEE) EleEle=1;
  if(foundJet && foundMT) MuoTau=1;
  if(foundJet && foundET) EleTau=1;
  //ELE-MUO ANALYSIS - SR
  if(foundJet && foundEM){
    cout<<"EleMuo "<<iEvent.id().event()<<"; muon pt "<<SelectedMuonEM->pt()<<"; ele pt "<<SelectedElectronEM->pt()<<
      "; jet pt "<<SelectedJet->pt()<<"; jet mass "<<prunedMass<<"; met "<<met->begin()->pt()<<endl;
    FillTree(0,TreeEleMuo,SelectedJet,SelectedTauFake,SelectedMuonEM,SelectedMuon1Fake,SelectedMuon2Fake,SelectedTrackerMuo,SelectedElectronEM,SelectedElectron1Fake,
	     SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMass, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuo, MuoMuo, EleEle, MuoTau, EleTau, iEvent);
  }
  //MUO-MUO ANALYSIS - SR
  if(foundJet && foundMM){
    cout<<"MuoMuo "<<iEvent.id().event()<<"; muon1 pt "<<SelectedMuon1MM->pt()<<"; muon2 pt "<<SelectedMuon2MM->pt()<<
      "; jet pt "<<SelectedJet->pt()<<"; jet mass "<<prunedMass<<"; met "<<met->begin()->pt()<<endl;
    FillTree(1,TreeMuoMuo,SelectedJet,SelectedTauFake,SelectedMuonFake,SelectedMuon1MM,SelectedMuon2MM,SelectedTrackerMuo,SelectedElectronFake,SelectedElectron1Fake,
	     SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMass, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuo, MuoMuo, EleEle, MuoTau, EleTau, iEvent);
  }
  //ELE-ELE ANALYSIS - SR
  if(foundJet && foundEE){
    cout<<"EleEle "<<iEvent.id().event()<<"; ele1 pt "<<SelectedElectron1EE->pt()<<"; ele2 pt "<<SelectedElectron2EE->pt()<<
      "; jet pt "<<SelectedJet->pt()<<"; jet mass "<<prunedMass<<"; met "<<met->begin()->pt()<<endl;
    FillTree(2,TreeEleEle,SelectedJet,SelectedTauFake,SelectedMuonFake,SelectedMuon1Fake,SelectedMuon2Fake,SelectedTrackerMuo,SelectedElectronFake,SelectedElectron1EE,
	     SelectedElectron2EE, vertices, metRaw, met, uncorrmet, prunedMass, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuo, MuoMuo, EleEle, MuoTau, EleTau, iEvent);
  }
  //MUO-TAU ANALYSIS - SR
  if(foundJet && foundMT){
    cout<<"MuoTau "<<iEvent.id().event()<<"; tau pt "<<SelectedTauMT->pt()<<"; muon pt "<<SelectedMuonMT->pt()<<
      "; jet pt "<<SelectedJet->pt()<<"; jet mass "<<prunedMass<<"; met "<<met->begin()->pt()<<endl;
    FillTree(3,TreeMuoTau,SelectedJet,SelectedTauMT,SelectedMuonMT,SelectedMuon1Fake,SelectedMuon2Fake,SelectedTrackerMuo,SelectedElectronFake,SelectedElectron1Fake,
	     SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMass, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuo, MuoMuo, EleEle, MuoTau, EleTau, iEvent);
  }
  //ELE-TAU ANALYSIS - SR
  if(foundJet && foundET){
    cout<<"EleTau "<<iEvent.id().event()<<"; tau pt "<<SelectedTauET->pt()<<"; electron pt "<<SelectedElectronET->pt()<<
      "; jet pt "<<SelectedJet->pt()<<"; jet mass "<<prunedMass<<"; met "<<met->begin()->pt()<<endl;
    FillTree(4,TreeEleTau,SelectedJet,SelectedTauET,SelectedMuonFake,SelectedMuon1Fake,SelectedMuon2Fake,SelectedTrackerMuo,SelectedElectronET,SelectedElectron1Fake,
	     SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMass, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuo, MuoMuo, EleEle, MuoTau, EleTau, iEvent);
  }

  int EleMuoSB1 = 0; int MuoMuoSB1 = 0; int EleEleSB1 = 0; int MuoTauSB1 = 0; int EleTauSB1 = 0;
  if(foundSB1Jet && foundSB1EM) EleMuoSB1=1;
  if(foundSB1Jet && foundSB1MM) MuoMuoSB1=1;
  if(foundSB1Jet && foundSB1EE) EleEleSB1=1;
  if(foundSB1Jet && foundSB1MT) MuoTauSB1=1;
  if(foundSB1Jet && foundSB1ET) EleTauSB1=1;
  //ELE-MUO ANALYSIS - SB1
  if(foundSB1Jet && foundSB1EM){
    FillTree(0,TreeSB1EleMuo,SelectedSB1Jet,SelectedTauFake,SelectedSB1MuonEM,SelectedMuon1Fake,SelectedMuon2Fake,SelectedTrackerMuo,SelectedSB1ElectronEM,
	     SelectedElectron1Fake,SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMassSB1, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB1, MuoMuoSB1, EleEleSB1, MuoTauSB1, EleTauSB1, iEvent);
  }
  //MUO-MUO ANALYSIS - SB1
  if(foundSB1Jet && foundSB1MM){
    FillTree(1,TreeSB1MuoMuo,SelectedSB1Jet,SelectedTauFake,SelectedMuonFake,SelectedSB1Muon1MM,SelectedSB1Muon2MM,SelectedTrackerMuo,SelectedElectronFake,
	     SelectedElectron1Fake,SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMassSB1, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB1, MuoMuoSB1, EleEleSB1, MuoTauSB1, EleTauSB1, iEvent);
  }
  //ELE-ELE ANALYSIS - SB1
  if(foundSB1Jet && foundSB1EE){
    FillTree(2,TreeSB1EleEle,SelectedSB1Jet,SelectedTauFake,SelectedMuonFake,SelectedMuon1Fake,SelectedMuon2Fake,SelectedTrackerMuo,SelectedElectronFake,
	     SelectedSB1Electron1EE,SelectedSB1Electron2EE, vertices, metRaw, met, uncorrmet, prunedMassSB1, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB1, MuoMuoSB1, EleEleSB1, MuoTauSB1, EleTauSB1, iEvent);
  }
  //MUO-TAU ANALYSIS - SB1
  if(foundSB1Jet && foundSB1MT){
    FillTree(3,TreeSB1MuoTau,SelectedSB1Jet,SelectedSB1TauMT,SelectedSB1MuonMT,SelectedMuon1Fake,SelectedMuon2Fake,SelectedTrackerMuo,SelectedElectronFake,
	     SelectedElectron1Fake,SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMassSB1, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB1, MuoMuoSB1, EleEleSB1, MuoTauSB1, EleTauSB1, iEvent);
  }
  //ELE-TAU ANALYSIS - SB1
  if(foundSB1Jet && foundSB1ET){
    FillTree(4,TreeSB1EleTau,SelectedSB1Jet,SelectedSB1TauET,SelectedMuonFake,SelectedMuon1Fake,SelectedMuon2Fake,SelectedTrackerMuo,SelectedSB1ElectronET,
	     SelectedElectron1Fake,SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMassSB1, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB1, MuoMuoSB1, EleEleSB1, MuoTauSB1, EleTauSB1, iEvent);
  }

  int EleMuoSB2 = 0; int MuoMuoSB2 = 0; int EleEleSB2 = 0; int MuoTauSB2 = 0; int EleTauSB2 = 0;
  if(foundSB2Jet && foundSB2EM) EleMuoSB2=1;
  if(foundSB2Jet && foundSB2MM) MuoMuoSB2=1;
  if(foundSB2Jet && foundSB2EE) EleEleSB2=1;
  if(foundSB2Jet && foundSB2MT) MuoTauSB2=1;
  if(foundSB2Jet && foundSB2ET) EleTauSB2=1;
  //ELE-MUO ANALYSIS - SB2
  if(foundSB2Jet && foundSB2EM){
    FillTree(0,TreeSB2EleMuo,SelectedSB2Jet,SelectedTauFake,SelectedSB2MuonEM,SelectedMuon1Fake,SelectedMuon2Fake,SelectedTrackerMuo,SelectedSB2ElectronEM,
	     SelectedElectron1Fake,SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMassSB2, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB2, MuoMuoSB2, EleEleSB2, MuoTauSB2, EleTauSB2, iEvent);
  }
  //MUO-MUO ANALYSIS - SB2
  if(foundSB2Jet && foundSB2MM){
    FillTree(1,TreeSB2MuoMuo,SelectedSB2Jet,SelectedTauFake,SelectedMuonFake,SelectedSB2Muon1MM,SelectedSB2Muon2MM,SelectedTrackerMuo,SelectedElectronFake,
	     SelectedElectron1Fake,SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMassSB2, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB2, MuoMuoSB2, EleEleSB2, MuoTauSB2, EleTauSB2, iEvent);
  }
  //ELE-ELE ANALYSIS - SB2
  if(foundSB2Jet && foundSB2EE){
    FillTree(2,TreeSB2EleEle,SelectedSB2Jet,SelectedTauFake,SelectedMuonFake,SelectedMuon1Fake,SelectedMuon2Fake,SelectedTrackerMuo,SelectedElectronFake,
	     SelectedSB2Electron1EE,SelectedSB2Electron2EE, vertices, metRaw, met, uncorrmet, prunedMassSB2, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB2, MuoMuoSB2, EleEleSB2, MuoTauSB2, EleTauSB2, iEvent);
  }
  //MUO-TAU ANALYSIS - SB2
  if(foundSB2Jet && foundSB2MT){
    FillTree(3,TreeSB2MuoTau,SelectedSB2Jet,SelectedSB2TauMT,SelectedSB2MuonMT,SelectedMuon1Fake,SelectedMuon2Fake,SelectedTrackerMuo,SelectedElectronFake,
	     SelectedElectron1Fake,SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMassSB2, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB2, MuoMuoSB2, EleEleSB2, MuoTauSB2, EleTauSB2, iEvent);
  }
  //ELE-TAU ANALYSIS - SB2
  if(foundSB2Jet && foundSB2ET){
    FillTree(4,TreeSB2EleTau,SelectedSB2Jet,SelectedSB2TauET,SelectedMuonFake,SelectedMuon1Fake,SelectedMuon2Fake,SelectedTrackerMuo,SelectedSB2ElectronET,
	     SelectedElectron1Fake,SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMassSB2, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB2, MuoMuoSB2, EleEleSB2, MuoTauSB2, EleTauSB2, iEvent);
  }
  
  int EleMuoSB3 = 0; int MuoMuoSB3 = 0; int EleEleSB3 = 0; int MuoTauSB3 = 0; int EleTauSB3 = 0;
  if(foundSB3Jet && foundSB3EM) EleMuoSB3=1;
  if(foundSB3Jet && foundSB3MM) MuoMuoSB3=1;
  if(foundSB3Jet && foundSB3EE) EleEleSB3=1;
  if(foundSB3Jet && foundSB3MT) MuoTauSB3=1;
  if(foundSB3Jet && foundSB3ET) EleTauSB3=1;
  //ELE-MUO ANALYSIS - SB3
  if(foundSB3Jet && foundSB3EM){
    FillTree(0,TreeSB3EleMuo,SelectedSB3Jet,SelectedTauFake,SelectedSB3MuonEM,SelectedMuon1Fake,SelectedMuon2Fake,SelectedTrackerMuo,SelectedSB3ElectronEM,
	     SelectedElectron1Fake,SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMassSB3, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB3, MuoMuoSB3, EleEleSB3, MuoTauSB3, EleTauSB3, iEvent);
  }
  //MUO-MUO ANALYSIS - SB3
  if(foundSB3Jet && foundSB3MM){
    FillTree(1,TreeSB3MuoMuo,SelectedSB3Jet,SelectedTauFake,SelectedMuonFake,SelectedSB3Muon1MM,SelectedSB3Muon2MM,SelectedTrackerMuo,SelectedElectronFake,
	     SelectedElectron1Fake,SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMassSB3, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB3, MuoMuoSB3, EleEleSB3, MuoTauSB3, EleTauSB3, iEvent);
  }
  //ELE-ELE ANALYSIS - SB3
  if(foundSB3Jet && foundSB3EE){
    FillTree(2,TreeSB3EleEle,SelectedSB3Jet,SelectedTauFake,SelectedMuonFake,SelectedMuon1Fake,SelectedMuon2Fake,SelectedTrackerMuo,SelectedElectronFake,
	     SelectedSB3Electron1EE,SelectedSB3Electron2EE, vertices, metRaw, met, uncorrmet, prunedMassSB3, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB3, MuoMuoSB3, EleEleSB3, MuoTauSB3, EleTauSB3, iEvent);
  }
  //MUO-TAU ANALYSIS - SB3
  if(foundSB3Jet && foundSB3MT){
    FillTree(3,TreeSB3MuoTau,SelectedSB3Jet,SelectedSB3TauMT,SelectedSB3MuonMT,SelectedMuon1Fake,SelectedMuon2Fake,SelectedTrackerMuo,SelectedElectronFake,
	     SelectedElectron1Fake,SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMassSB3, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB3, MuoMuoSB3, EleEleSB3, MuoTauSB3, EleTauSB3, iEvent);
  }
  //ELE-TAU ANALYSIS - SB3
  if(foundSB3Jet && foundSB3ET){
    FillTree(4,TreeSB3EleTau,SelectedSB3Jet,SelectedSB3TauET,SelectedMuonFake,SelectedMuon1Fake,SelectedMuon2Fake,SelectedTrackerMuo,SelectedSB3ElectronET,
	     SelectedElectron1Fake,SelectedElectron2Fake, vertices, metRaw, met, uncorrmet, prunedMassSB3, isFired_HLT, isFired_HLT_PFJet320, 
	     isFired_HLT_HT650, MyWeight, genEvent, rho, EleMuoSB3, MuoMuoSB3, EleEleSB3, MuoTauSB3, EleTauSB3, iEvent);
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
  TreeEleMuo->Branch("lep1DETIso", &m_lep1DETIso, "lep1DETIso/f");
  TreeEleMuo->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeEleMuo->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeEleMuo->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeEleMuo->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeEleMuo->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeEleMuo->Branch("lep2DETIso", &m_lep2DETIso, "lep2DETIso/f");
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
  TreeEleMuo->Branch("njet1", &m_njet1, "njet1/i");
  TreeEleMuo->Branch("nbtagsL1", &m_nbtagsL1, "nbtagsL1/i");
  TreeEleMuo->Branch("nbtagsM1", &m_nbtagsM1, "nbtagsM1/i");
  TreeEleMuo->Branch("nbtagsT1", &m_nbtagsT1, "nbtagsT1/i");
  TreeEleMuo->Branch("njet2", &m_njet2, "njet2/i");
  TreeEleMuo->Branch("nbtagsL2", &m_nbtagsL2, "nbtagsL2/i");
  TreeEleMuo->Branch("nbtagsM2", &m_nbtagsM2, "nbtagsM2/i");
  TreeEleMuo->Branch("nbtagsT2", &m_nbtagsT2, "nbtagsT2/i");
  TreeEleMuo->Branch("njet3", &m_njet3, "njet3/i");
  TreeEleMuo->Branch("nbtagsL3", &m_nbtagsL3, "nbtagsL3/i");
  TreeEleMuo->Branch("nbtagsM3", &m_nbtagsM3, "nbtagsM3/i");
  TreeEleMuo->Branch("nbtagsT3", &m_nbtagsT3, "nbtagsT3/i");
  TreeEleMuo->Branch("njet4", &m_njet4, "njet4/i");
  TreeEleMuo->Branch("nbtagsL4", &m_nbtagsL4, "nbtagsL4/i");
  TreeEleMuo->Branch("nbtagsM4", &m_nbtagsM4, "nbtagsM4/i");
  TreeEleMuo->Branch("nbtagsT4", &m_nbtagsT4, "nbtagsT4/i");
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
  TreeMuoMuo->Branch("lep1DETIso", &m_lep1DETIso, "lep1DETIso/f");
  TreeMuoMuo->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeMuoMuo->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeMuoMuo->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeMuoMuo->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeMuoMuo->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeMuoMuo->Branch("lep2DETIso", &m_lep2DETIso, "lep2DETIso/f");
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
  TreeMuoMuo->Branch("njet1", &m_njet1, "njet1/i");
  TreeMuoMuo->Branch("nbtagsL1", &m_nbtagsL1, "nbtagsL1/i");
  TreeMuoMuo->Branch("nbtagsM1", &m_nbtagsM1, "nbtagsM1/i");
  TreeMuoMuo->Branch("nbtagsT1", &m_nbtagsT1, "nbtagsT1/i");
  TreeMuoMuo->Branch("njet2", &m_njet2, "njet2/i");
  TreeMuoMuo->Branch("nbtagsL2", &m_nbtagsL2, "nbtagsL2/i");
  TreeMuoMuo->Branch("nbtagsM2", &m_nbtagsM2, "nbtagsM2/i");
  TreeMuoMuo->Branch("nbtagsT2", &m_nbtagsT2, "nbtagsT2/i");
  TreeMuoMuo->Branch("njet3", &m_njet3, "njet3/i");
  TreeMuoMuo->Branch("nbtagsL3", &m_nbtagsL3, "nbtagsL3/i");
  TreeMuoMuo->Branch("nbtagsM3", &m_nbtagsM3, "nbtagsM3/i");
  TreeMuoMuo->Branch("nbtagsT3", &m_nbtagsT3, "nbtagsT3/i");
  TreeMuoMuo->Branch("njet4", &m_njet4, "njet4/i");
  TreeMuoMuo->Branch("nbtagsL4", &m_nbtagsL4, "nbtagsL4/i");
  TreeMuoMuo->Branch("nbtagsM4", &m_nbtagsM4, "nbtagsM4/i");
  TreeMuoMuo->Branch("nbtagsT4", &m_nbtagsT4, "nbtagsT4/i");
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
  TreeEleEle->Branch("lep1DETIso", &m_lep1DETIso, "lep1DETIso/f");
  TreeEleEle->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeEleEle->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeEleEle->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeEleEle->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeEleEle->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeEleEle->Branch("lep2DETIso", &m_lep2DETIso, "lep2DETIso/f");
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
  TreeEleEle->Branch("njet1", &m_njet1, "njet1/i");
  TreeEleEle->Branch("nbtagsL1", &m_nbtagsL1, "nbtagsL1/i");
  TreeEleEle->Branch("nbtagsM1", &m_nbtagsM1, "nbtagsM1/i");
  TreeEleEle->Branch("nbtagsT1", &m_nbtagsT1, "nbtagsT1/i");
  TreeEleEle->Branch("njet2", &m_njet2, "njet2/i");
  TreeEleEle->Branch("nbtagsL2", &m_nbtagsL2, "nbtagsL2/i");
  TreeEleEle->Branch("nbtagsM2", &m_nbtagsM2, "nbtagsM2/i");
  TreeEleEle->Branch("nbtagsT2", &m_nbtagsT2, "nbtagsT2/i");
  TreeEleEle->Branch("njet3", &m_njet3, "njet3/i");
  TreeEleEle->Branch("nbtagsL3", &m_nbtagsL3, "nbtagsL3/i");
  TreeEleEle->Branch("nbtagsM3", &m_nbtagsM3, "nbtagsM3/i");
  TreeEleEle->Branch("nbtagsT3", &m_nbtagsT3, "nbtagsT3/i");
  TreeEleEle->Branch("njet4", &m_njet4, "njet4/i");
  TreeEleEle->Branch("nbtagsL4", &m_nbtagsL4, "nbtagsL4/i");
  TreeEleEle->Branch("nbtagsM4", &m_nbtagsM4, "nbtagsM4/i");
  TreeEleEle->Branch("nbtagsT4", &m_nbtagsT4, "nbtagsT4/i");
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
  TreeMuoTau->Branch("lep1DETIso", &m_lep1DETIso, "lep1DETIso/f");
  TreeMuoTau->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeMuoTau->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeMuoTau->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeMuoTau->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeMuoTau->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeMuoTau->Branch("lep2DETIso", &m_lep2DETIso, "lep2DETIso/f");
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
  TreeMuoTau->Branch("njet1", &m_njet1, "njet1/i");
  TreeMuoTau->Branch("nbtagsL1", &m_nbtagsL1, "nbtagsL1/i");
  TreeMuoTau->Branch("nbtagsM1", &m_nbtagsM1, "nbtagsM1/i");
  TreeMuoTau->Branch("nbtagsT1", &m_nbtagsT1, "nbtagsT1/i");
  TreeMuoTau->Branch("njet2", &m_njet2, "njet2/i");
  TreeMuoTau->Branch("nbtagsL2", &m_nbtagsL2, "nbtagsL2/i");
  TreeMuoTau->Branch("nbtagsM2", &m_nbtagsM2, "nbtagsM2/i");
  TreeMuoTau->Branch("nbtagsT2", &m_nbtagsT2, "nbtagsT2/i");
  TreeMuoTau->Branch("njet3", &m_njet3, "njet3/i");
  TreeMuoTau->Branch("nbtagsL3", &m_nbtagsL3, "nbtagsL3/i");
  TreeMuoTau->Branch("nbtagsM3", &m_nbtagsM3, "nbtagsM3/i");
  TreeMuoTau->Branch("nbtagsT3", &m_nbtagsT3, "nbtagsT3/i");
  TreeMuoTau->Branch("njet4", &m_njet4, "njet4/i");
  TreeMuoTau->Branch("nbtagsL4", &m_nbtagsL4, "nbtagsL4/i");
  TreeMuoTau->Branch("nbtagsM4", &m_nbtagsM4, "nbtagsM4/i");
  TreeMuoTau->Branch("nbtagsT4", &m_nbtagsT4, "nbtagsT4/i");
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
  TreeEleTau->Branch("lep1DETIso", &m_lep1DETIso, "lep1DETIso/f");
  TreeEleTau->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeEleTau->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeEleTau->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeEleTau->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeEleTau->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeEleTau->Branch("lep2DETIso", &m_lep2DETIso, "lep2DETIso/f");
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
  TreeEleTau->Branch("njet1", &m_njet1, "njet1/i");
  TreeEleTau->Branch("nbtagsL1", &m_nbtagsL1, "nbtagsL1/i");
  TreeEleTau->Branch("nbtagsM1", &m_nbtagsM1, "nbtagsM1/i");
  TreeEleTau->Branch("nbtagsT1", &m_nbtagsT1, "nbtagsT1/i");
  TreeEleTau->Branch("njet2", &m_njet2, "njet2/i");
  TreeEleTau->Branch("nbtagsL2", &m_nbtagsL2, "nbtagsL2/i");
  TreeEleTau->Branch("nbtagsM2", &m_nbtagsM2, "nbtagsM2/i");
  TreeEleTau->Branch("nbtagsT2", &m_nbtagsT2, "nbtagsT2/i");
  TreeEleTau->Branch("njet3", &m_njet3, "njet3/i");
  TreeEleTau->Branch("nbtagsL3", &m_nbtagsL3, "nbtagsL3/i");
  TreeEleTau->Branch("nbtagsM3", &m_nbtagsM3, "nbtagsM3/i");
  TreeEleTau->Branch("nbtagsT3", &m_nbtagsT3, "nbtagsT3/i");
  TreeEleTau->Branch("njet4", &m_njet4, "njet4/i");
  TreeEleTau->Branch("nbtagsL4", &m_nbtagsL4, "nbtagsL4/i");
  TreeEleTau->Branch("nbtagsM4", &m_nbtagsM4, "nbtagsM4/i");
  TreeEleTau->Branch("nbtagsT4", &m_nbtagsT4, "nbtagsT4/i");
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
  TreeSB1EleMuo->Branch("lep1DETIso", &m_lep1DETIso, "lep1DETIso/f");
  TreeSB1EleMuo->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB1EleMuo->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB1EleMuo->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB1EleMuo->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB1EleMuo->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeSB1EleMuo->Branch("lep2DETIso", &m_lep2DETIso, "lep2DETIso/f");
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
  TreeSB1EleMuo->Branch("njet1", &m_njet1, "njet1/i");
  TreeSB1EleMuo->Branch("nbtagsL1", &m_nbtagsL1, "nbtagsL1/i");
  TreeSB1EleMuo->Branch("nbtagsM1", &m_nbtagsM1, "nbtagsM1/i");
  TreeSB1EleMuo->Branch("nbtagsT1", &m_nbtagsT1, "nbtagsT1/i");
  TreeSB1EleMuo->Branch("njet2", &m_njet2, "njet2/i");
  TreeSB1EleMuo->Branch("nbtagsL2", &m_nbtagsL2, "nbtagsL2/i");
  TreeSB1EleMuo->Branch("nbtagsM2", &m_nbtagsM2, "nbtagsM2/i");
  TreeSB1EleMuo->Branch("nbtagsT2", &m_nbtagsT2, "nbtagsT2/i");
  TreeSB1EleMuo->Branch("njet3", &m_njet3, "njet3/i");
  TreeSB1EleMuo->Branch("nbtagsL3", &m_nbtagsL3, "nbtagsL3/i");
  TreeSB1EleMuo->Branch("nbtagsM3", &m_nbtagsM3, "nbtagsM3/i");
  TreeSB1EleMuo->Branch("nbtagsT3", &m_nbtagsT3, "nbtagsT3/i");
  TreeSB1EleMuo->Branch("njet4", &m_njet4, "njet4/i");
  TreeSB1EleMuo->Branch("nbtagsL4", &m_nbtagsL4, "nbtagsL4/i");
  TreeSB1EleMuo->Branch("nbtagsM4", &m_nbtagsM4, "nbtagsM4/i");
  TreeSB1EleMuo->Branch("nbtagsT4", &m_nbtagsT4, "nbtagsT4/i");
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
  TreeSB1MuoMuo->Branch("lep1DETIso", &m_lep1DETIso, "lep1DETIso/f");
  TreeSB1MuoMuo->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB1MuoMuo->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB1MuoMuo->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB1MuoMuo->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB1MuoMuo->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeSB1MuoMuo->Branch("lep2DETIso", &m_lep2DETIso, "lep2DETIso/f");
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
  TreeSB1MuoMuo->Branch("njet1", &m_njet1, "njet1/i");
  TreeSB1MuoMuo->Branch("nbtagsL1", &m_nbtagsL1, "nbtagsL1/i");
  TreeSB1MuoMuo->Branch("nbtagsM1", &m_nbtagsM1, "nbtagsM1/i");
  TreeSB1MuoMuo->Branch("nbtagsT1", &m_nbtagsT1, "nbtagsT1/i");
  TreeSB1MuoMuo->Branch("njet2", &m_njet2, "njet2/i");
  TreeSB1MuoMuo->Branch("nbtagsL2", &m_nbtagsL2, "nbtagsL2/i");
  TreeSB1MuoMuo->Branch("nbtagsM2", &m_nbtagsM2, "nbtagsM2/i");
  TreeSB1MuoMuo->Branch("nbtagsT2", &m_nbtagsT2, "nbtagsT2/i");
  TreeSB1MuoMuo->Branch("njet3", &m_njet3, "njet3/i");
  TreeSB1MuoMuo->Branch("nbtagsL3", &m_nbtagsL3, "nbtagsL3/i");
  TreeSB1MuoMuo->Branch("nbtagsM3", &m_nbtagsM3, "nbtagsM3/i");
  TreeSB1MuoMuo->Branch("nbtagsT3", &m_nbtagsT3, "nbtagsT3/i");
  TreeSB1MuoMuo->Branch("njet4", &m_njet4, "njet4/i");
  TreeSB1MuoMuo->Branch("nbtagsL4", &m_nbtagsL4, "nbtagsL4/i");
  TreeSB1MuoMuo->Branch("nbtagsM4", &m_nbtagsM4, "nbtagsM4/i");
  TreeSB1MuoMuo->Branch("nbtagsT4", &m_nbtagsT4, "nbtagsT4/i");
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
  TreeSB1EleEle->Branch("lep1DETIso", &m_lep1DETIso, "lep1DETIso/f");
  TreeSB1EleEle->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB1EleEle->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB1EleEle->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB1EleEle->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB1EleEle->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeSB1EleEle->Branch("lep2DETIso", &m_lep2DETIso, "lep2DETIso/f");
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
  TreeSB1EleEle->Branch("njet1", &m_njet1, "njet1/i");
  TreeSB1EleEle->Branch("nbtagsL1", &m_nbtagsL1, "nbtagsL1/i");
  TreeSB1EleEle->Branch("nbtagsM1", &m_nbtagsM1, "nbtagsM1/i");
  TreeSB1EleEle->Branch("nbtagsT1", &m_nbtagsT1, "nbtagsT1/i");
  TreeSB1EleEle->Branch("njet2", &m_njet2, "njet2/i");
  TreeSB1EleEle->Branch("nbtagsL2", &m_nbtagsL2, "nbtagsL2/i");
  TreeSB1EleEle->Branch("nbtagsM2", &m_nbtagsM2, "nbtagsM2/i");
  TreeSB1EleEle->Branch("nbtagsT2", &m_nbtagsT2, "nbtagsT2/i");
  TreeSB1EleEle->Branch("njet3", &m_njet3, "njet3/i");
  TreeSB1EleEle->Branch("nbtagsL3", &m_nbtagsL3, "nbtagsL3/i");
  TreeSB1EleEle->Branch("nbtagsM3", &m_nbtagsM3, "nbtagsM3/i");
  TreeSB1EleEle->Branch("nbtagsT3", &m_nbtagsT3, "nbtagsT3/i");
  TreeSB1EleEle->Branch("njet4", &m_njet4, "njet4/i");
  TreeSB1EleEle->Branch("nbtagsL4", &m_nbtagsL4, "nbtagsL4/i");
  TreeSB1EleEle->Branch("nbtagsM4", &m_nbtagsM4, "nbtagsM4/i");
  TreeSB1EleEle->Branch("nbtagsT4", &m_nbtagsT4, "nbtagsT4/i");
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
  TreeSB1MuoTau->Branch("lep1DETIso", &m_lep1DETIso, "lep1DETIso/f");
  TreeSB1MuoTau->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB1MuoTau->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB1MuoTau->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB1MuoTau->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB1MuoTau->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeSB1MuoTau->Branch("lep2DETIso", &m_lep2DETIso, "lep2DETIso/f");
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
  TreeSB1MuoTau->Branch("njet1", &m_njet1, "njet1/i");
  TreeSB1MuoTau->Branch("nbtagsL1", &m_nbtagsL1, "nbtagsL1/i");
  TreeSB1MuoTau->Branch("nbtagsM1", &m_nbtagsM1, "nbtagsM1/i");
  TreeSB1MuoTau->Branch("nbtagsT1", &m_nbtagsT1, "nbtagsT1/i");
  TreeSB1MuoTau->Branch("njet2", &m_njet2, "njet2/i");
  TreeSB1MuoTau->Branch("nbtagsL2", &m_nbtagsL2, "nbtagsL2/i");
  TreeSB1MuoTau->Branch("nbtagsM2", &m_nbtagsM2, "nbtagsM2/i");
  TreeSB1MuoTau->Branch("nbtagsT2", &m_nbtagsT2, "nbtagsT2/i");
  TreeSB1MuoTau->Branch("njet3", &m_njet3, "njet3/i");
  TreeSB1MuoTau->Branch("nbtagsL3", &m_nbtagsL3, "nbtagsL3/i");
  TreeSB1MuoTau->Branch("nbtagsM3", &m_nbtagsM3, "nbtagsM3/i");
  TreeSB1MuoTau->Branch("nbtagsT3", &m_nbtagsT3, "nbtagsT3/i");
  TreeSB1MuoTau->Branch("njet4", &m_njet4, "njet4/i");
  TreeSB1MuoTau->Branch("nbtagsL4", &m_nbtagsL4, "nbtagsL4/i");
  TreeSB1MuoTau->Branch("nbtagsM4", &m_nbtagsM4, "nbtagsM4/i");
  TreeSB1MuoTau->Branch("nbtagsT4", &m_nbtagsT4, "nbtagsT4/i");
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
  TreeSB1EleTau->Branch("lep1DETIso", &m_lep1DETIso, "lep1DETIso/f");
  TreeSB1EleTau->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB1EleTau->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB1EleTau->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB1EleTau->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB1EleTau->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeSB1EleTau->Branch("lep2DETIso", &m_lep2DETIso, "lep2DETIso/f");
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
  TreeSB1EleTau->Branch("njet1", &m_njet1, "njet1/i");
  TreeSB1EleTau->Branch("nbtagsL1", &m_nbtagsL1, "nbtagsL1/i");
  TreeSB1EleTau->Branch("nbtagsM1", &m_nbtagsM1, "nbtagsM1/i");
  TreeSB1EleTau->Branch("nbtagsT1", &m_nbtagsT1, "nbtagsT1/i");
  TreeSB1EleTau->Branch("njet2", &m_njet2, "njet2/i");
  TreeSB1EleTau->Branch("nbtagsL2", &m_nbtagsL2, "nbtagsL2/i");
  TreeSB1EleTau->Branch("nbtagsM2", &m_nbtagsM2, "nbtagsM2/i");
  TreeSB1EleTau->Branch("nbtagsT2", &m_nbtagsT2, "nbtagsT2/i");
  TreeSB1EleTau->Branch("njet3", &m_njet3, "njet3/i");
  TreeSB1EleTau->Branch("nbtagsL3", &m_nbtagsL3, "nbtagsL3/i");
  TreeSB1EleTau->Branch("nbtagsM3", &m_nbtagsM3, "nbtagsM3/i");
  TreeSB1EleTau->Branch("nbtagsT3", &m_nbtagsT3, "nbtagsT3/i");
  TreeSB1EleTau->Branch("njet4", &m_njet4, "njet4/i");
  TreeSB1EleTau->Branch("nbtagsL4", &m_nbtagsL4, "nbtagsL4/i");
  TreeSB1EleTau->Branch("nbtagsM4", &m_nbtagsM4, "nbtagsM4/i");
  TreeSB1EleTau->Branch("nbtagsT4", &m_nbtagsT4, "nbtagsT4/i");
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
  TreeSB2EleMuo->Branch("lep1DETIso", &m_lep1DETIso, "lep1DETIso/f");
  TreeSB2EleMuo->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB2EleMuo->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB2EleMuo->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB2EleMuo->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB2EleMuo->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeSB2EleMuo->Branch("lep2DETIso", &m_lep2DETIso, "lep2DETIso/f");
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
  TreeSB2EleMuo->Branch("njet1", &m_njet1, "njet1/i");
  TreeSB2EleMuo->Branch("nbtagsL1", &m_nbtagsL1, "nbtagsL1/i");
  TreeSB2EleMuo->Branch("nbtagsM1", &m_nbtagsM1, "nbtagsM1/i");
  TreeSB2EleMuo->Branch("nbtagsT1", &m_nbtagsT1, "nbtagsT1/i");
  TreeSB2EleMuo->Branch("njet2", &m_njet2, "njet2/i");
  TreeSB2EleMuo->Branch("nbtagsL2", &m_nbtagsL2, "nbtagsL2/i");
  TreeSB2EleMuo->Branch("nbtagsM2", &m_nbtagsM2, "nbtagsM2/i");
  TreeSB2EleMuo->Branch("nbtagsT2", &m_nbtagsT2, "nbtagsT2/i");
  TreeSB2EleMuo->Branch("njet3", &m_njet3, "njet3/i");
  TreeSB2EleMuo->Branch("nbtagsL3", &m_nbtagsL3, "nbtagsL3/i");
  TreeSB2EleMuo->Branch("nbtagsM3", &m_nbtagsM3, "nbtagsM3/i");
  TreeSB2EleMuo->Branch("nbtagsT3", &m_nbtagsT3, "nbtagsT3/i");
  TreeSB2EleMuo->Branch("njet4", &m_njet4, "njet4/i");
  TreeSB2EleMuo->Branch("nbtagsL4", &m_nbtagsL4, "nbtagsL4/i");
  TreeSB2EleMuo->Branch("nbtagsM4", &m_nbtagsM4, "nbtagsM4/i");
  TreeSB2EleMuo->Branch("nbtagsT4", &m_nbtagsT4, "nbtagsT4/i");
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
  TreeSB2EleMuo->Branch("lep1DETIso", &m_lep1DETIso, "lep1DETIso/f");
  TreeSB2MuoMuo->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB2MuoMuo->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB2MuoMuo->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB2MuoMuo->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB2MuoMuo->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeSB2MuoMuo->Branch("lep2DETIso", &m_lep2DETIso, "lep2DETIso/f");
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
  TreeSB2MuoMuo->Branch("njet1", &m_njet1, "njet1/i");
  TreeSB2MuoMuo->Branch("nbtagsL1", &m_nbtagsL1, "nbtagsL1/i");
  TreeSB2MuoMuo->Branch("nbtagsM1", &m_nbtagsM1, "nbtagsM1/i");
  TreeSB2MuoMuo->Branch("nbtagsT1", &m_nbtagsT1, "nbtagsT1/i");
  TreeSB2MuoMuo->Branch("njet2", &m_njet2, "njet2/i");
  TreeSB2MuoMuo->Branch("nbtagsL2", &m_nbtagsL2, "nbtagsL2/i");
  TreeSB2MuoMuo->Branch("nbtagsM2", &m_nbtagsM2, "nbtagsM2/i");
  TreeSB2MuoMuo->Branch("nbtagsT2", &m_nbtagsT2, "nbtagsT2/i");
  TreeSB2MuoMuo->Branch("njet3", &m_njet3, "njet3/i");
  TreeSB2MuoMuo->Branch("nbtagsL3", &m_nbtagsL3, "nbtagsL3/i");
  TreeSB2MuoMuo->Branch("nbtagsM3", &m_nbtagsM3, "nbtagsM3/i");
  TreeSB2MuoMuo->Branch("nbtagsT3", &m_nbtagsT3, "nbtagsT3/i");
  TreeSB2MuoMuo->Branch("njet4", &m_njet4, "njet4/i");
  TreeSB2MuoMuo->Branch("nbtagsL4", &m_nbtagsL4, "nbtagsL4/i");
  TreeSB2MuoMuo->Branch("nbtagsM4", &m_nbtagsM4, "nbtagsM4/i");
  TreeSB2MuoMuo->Branch("nbtagsT4", &m_nbtagsT4, "nbtagsT4/i");
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
  TreeSB2EleEle->Branch("lep1DETIso", &m_lep1DETIso, "lep1DETIso/f");
  TreeSB2EleEle->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB2EleEle->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB2EleEle->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB2EleEle->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB2EleEle->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeSB2EleEle->Branch("lep2DETIso", &m_lep2DETIso, "lep2DETIso/f");
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
  TreeSB2EleEle->Branch("njet1", &m_njet1, "njet1/i");
  TreeSB2EleEle->Branch("nbtagsL1", &m_nbtagsL1, "nbtagsL1/i");
  TreeSB2EleEle->Branch("nbtagsM1", &m_nbtagsM1, "nbtagsM1/i");
  TreeSB2EleEle->Branch("nbtagsT1", &m_nbtagsT1, "nbtagsT1/i");
  TreeSB2EleEle->Branch("njet2", &m_njet2, "njet2/i");
  TreeSB2EleEle->Branch("nbtagsL2", &m_nbtagsL2, "nbtagsL2/i");
  TreeSB2EleEle->Branch("nbtagsM2", &m_nbtagsM2, "nbtagsM2/i");
  TreeSB2EleEle->Branch("nbtagsT2", &m_nbtagsT2, "nbtagsT2/i");
  TreeSB2EleEle->Branch("njet3", &m_njet3, "njet3/i");
  TreeSB2EleEle->Branch("nbtagsL3", &m_nbtagsL3, "nbtagsL3/i");
  TreeSB2EleEle->Branch("nbtagsM3", &m_nbtagsM3, "nbtagsM3/i");
  TreeSB2EleEle->Branch("nbtagsT3", &m_nbtagsT3, "nbtagsT3/i");
  TreeSB2EleEle->Branch("njet4", &m_njet4, "njet4/i");
  TreeSB2EleEle->Branch("nbtagsL4", &m_nbtagsL4, "nbtagsL4/i");
  TreeSB2EleEle->Branch("nbtagsM4", &m_nbtagsM4, "nbtagsM4/i");
  TreeSB2EleEle->Branch("nbtagsT4", &m_nbtagsT4, "nbtagsT4/i");
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
  TreeSB2MuoTau->Branch("lep1DETIso", &m_lep1DETIso, "lep1DETIso/f");
  TreeSB2MuoTau->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB2MuoTau->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB2MuoTau->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB2MuoTau->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB2MuoTau->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeSB2MuoTau->Branch("lep2DETIso", &m_lep2DETIso, "lep2DETIso/f");
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
  TreeSB2MuoTau->Branch("njet1", &m_njet1, "njet1/i");
  TreeSB2MuoTau->Branch("nbtagsL1", &m_nbtagsL1, "nbtagsL1/i");
  TreeSB2MuoTau->Branch("nbtagsM1", &m_nbtagsM1, "nbtagsM1/i");
  TreeSB2MuoTau->Branch("nbtagsT1", &m_nbtagsT1, "nbtagsT1/i");
  TreeSB2MuoTau->Branch("njet2", &m_njet2, "njet2/i");
  TreeSB2MuoTau->Branch("nbtagsL2", &m_nbtagsL2, "nbtagsL2/i");
  TreeSB2MuoTau->Branch("nbtagsM2", &m_nbtagsM2, "nbtagsM2/i");
  TreeSB2MuoTau->Branch("nbtagsT2", &m_nbtagsT2, "nbtagsT2/i");
  TreeSB2MuoTau->Branch("njet3", &m_njet3, "njet3/i");
  TreeSB2MuoTau->Branch("nbtagsL3", &m_nbtagsL3, "nbtagsL3/i");
  TreeSB2MuoTau->Branch("nbtagsM3", &m_nbtagsM3, "nbtagsM3/i");
  TreeSB2MuoTau->Branch("nbtagsT3", &m_nbtagsT3, "nbtagsT3/i");
  TreeSB2MuoTau->Branch("njet4", &m_njet4, "njet4/i");
  TreeSB2MuoTau->Branch("nbtagsL4", &m_nbtagsL4, "nbtagsL4/i");
  TreeSB2MuoTau->Branch("nbtagsM4", &m_nbtagsM4, "nbtagsM4/i");
  TreeSB2MuoTau->Branch("nbtagsT4", &m_nbtagsT4, "nbtagsT4/i");
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
  TreeSB2EleTau->Branch("lep1DETIso", &m_lep1DETIso, "lep1DETIso/f");
  TreeSB2EleTau->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB2EleTau->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB2EleTau->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB2EleTau->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB2EleTau->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeSB2EleTau->Branch("lep2DETIso", &m_lep2DETIso, "lep2DETIso/f");
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
  TreeSB2EleTau->Branch("njet1", &m_njet1, "njet1/i");
  TreeSB2EleTau->Branch("nbtagsL1", &m_nbtagsL1, "nbtagsL1/i");
  TreeSB2EleTau->Branch("nbtagsM1", &m_nbtagsM1, "nbtagsM1/i");
  TreeSB2EleTau->Branch("nbtagsT1", &m_nbtagsT1, "nbtagsT1/i");
  TreeSB2EleTau->Branch("njet2", &m_njet2, "njet2/i");
  TreeSB2EleTau->Branch("nbtagsL2", &m_nbtagsL2, "nbtagsL2/i");
  TreeSB2EleTau->Branch("nbtagsM2", &m_nbtagsM2, "nbtagsM2/i");
  TreeSB2EleTau->Branch("nbtagsT2", &m_nbtagsT2, "nbtagsT2/i");
  TreeSB2EleTau->Branch("njet3", &m_njet3, "njet3/i");
  TreeSB2EleTau->Branch("nbtagsL3", &m_nbtagsL3, "nbtagsL3/i");
  TreeSB2EleTau->Branch("nbtagsM3", &m_nbtagsM3, "nbtagsM3/i");
  TreeSB2EleTau->Branch("nbtagsT3", &m_nbtagsT3, "nbtagsT3/i");
  TreeSB2EleTau->Branch("njet4", &m_njet4, "njet4/i");
  TreeSB2EleTau->Branch("nbtagsL4", &m_nbtagsL4, "nbtagsL4/i");
  TreeSB2EleTau->Branch("nbtagsM4", &m_nbtagsM4, "nbtagsM4/i");
  TreeSB2EleTau->Branch("nbtagsT4", &m_nbtagsT4, "nbtagsT4/i");
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

  TreeSB3EleMuo = fs->make<TTree>("TreeSB3EleMuo", "TreeSB3EleMuo");
  TreeSB3EleMuo->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB3EleMuo->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB3EleMuo->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB3EleMuo->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB3EleMuo->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB3EleMuo->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB3EleMuo->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeSB3EleMuo->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeSB3EleMuo->Branch("dPhiLep1Met", &m_dPhiLep1Met, "dPhiLep1Met/f");
  TreeSB3EleMuo->Branch("dRLep1Met", &m_dRLep1Met, "dRLep1Met/f");
  TreeSB3EleMuo->Branch("dRLep1Lep2", &m_dRLep1Lep2, "dRLep1Lep2/f");
  TreeSB3EleMuo->Branch("dPhiLep2Met", &m_dPhiLep2Met, "dPhiLep2Met/f");
  TreeSB3EleMuo->Branch("dRLep2Met", &m_dRLep2Met, "dRLep2Met/f");
  TreeSB3EleMuo->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB3EleMuo->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB3EleMuo->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB3EleMuo->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB3EleMuo->Branch("lep1Pt", &m_lep1Pt, "lep1Pt/f");
  TreeSB3EleMuo->Branch("lep1Eta", &m_lep1Eta, "lep1Eta/f");
  TreeSB3EleMuo->Branch("lep1Charge", &m_lep1Charge, "lep1Charge/f");
  TreeSB3EleMuo->Branch("lep1PFIso", &m_lep1PFIso, "lep1PFIso/f");
  TreeSB3EleMuo->Branch("lep1CorrPFIso", &m_lep1CorrPFIso, "lep1CorrPFIso/f");
  TreeSB3EleMuo->Branch("lep1DETIso", &m_lep1DETIso, "lep1DETIso/f");
  TreeSB3EleMuo->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB3EleMuo->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB3EleMuo->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB3EleMuo->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB3EleMuo->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeSB3EleMuo->Branch("lep2DETIso", &m_lep2DETIso, "lep2DETIso/f");
  TreeSB3EleMuo->Branch("charge", &m_charge, "charge/f");
  TreeSB3EleMuo->Branch("met", &m_met, "met/f");
  TreeSB3EleMuo->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeSB3EleMuo->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeSB3EleMuo->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  TreeSB3EleMuo->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeSB3EleMuo->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeSB3EleMuo->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeSB3EleMuo->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeSB3EleMuo->Branch("PtSvfit", &m_PtSvfit, "PtSvfit/f");
  TreeSB3EleMuo->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeSB3EleMuo->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeSB3EleMuo->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeSB3EleMuo->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeSB3EleMuo->Branch("njet1", &m_njet1, "njet1/i");
  TreeSB3EleMuo->Branch("nbtagsL1", &m_nbtagsL1, "nbtagsL1/i");
  TreeSB3EleMuo->Branch("nbtagsM1", &m_nbtagsM1, "nbtagsM1/i");
  TreeSB3EleMuo->Branch("nbtagsT1", &m_nbtagsT1, "nbtagsT1/i");
  TreeSB3EleMuo->Branch("njet2", &m_njet2, "njet2/i");
  TreeSB3EleMuo->Branch("nbtagsL2", &m_nbtagsL2, "nbtagsL2/i");
  TreeSB3EleMuo->Branch("nbtagsM2", &m_nbtagsM2, "nbtagsM2/i");
  TreeSB3EleMuo->Branch("nbtagsT2", &m_nbtagsT2, "nbtagsT2/i");
  TreeSB3EleMuo->Branch("njet3", &m_njet3, "njet3/i");
  TreeSB3EleMuo->Branch("nbtagsL3", &m_nbtagsL3, "nbtagsL3/i");
  TreeSB3EleMuo->Branch("nbtagsM3", &m_nbtagsM3, "nbtagsM3/i");
  TreeSB3EleMuo->Branch("nbtagsT3", &m_nbtagsT3, "nbtagsT3/i");
  TreeSB3EleMuo->Branch("njet4", &m_njet4, "njet4/i");
  TreeSB3EleMuo->Branch("nbtagsL4", &m_nbtagsL4, "nbtagsL4/i");
  TreeSB3EleMuo->Branch("nbtagsM4", &m_nbtagsM4, "nbtagsM4/i");
  TreeSB3EleMuo->Branch("nbtagsT4", &m_nbtagsT4, "nbtagsT4/i");
  TreeSB3EleMuo->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeSB3EleMuo->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeSB3EleMuo->Branch("trigger", &m_trigger, "trigger/i");
  TreeSB3EleMuo->Branch("sideband", &m_sideband, "sideband/i");
  TreeSB3EleMuo->Branch("NVertices", &m_NVertices, "NVertices/i");
  TreeSB3EleMuo->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeSB3EleMuo->Branch("metPx", &m_metPx, "metPx/f");
  TreeSB3EleMuo->Branch("metPy", &m_metPy, "metPy/f");
  TreeSB3EleMuo->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeSB3EleMuo->Branch("xsec", &m_xsec, "xsec/d");
  TreeSB3EleMuo->Branch("lumi", &m_lumi, "lumi/d");
  TreeSB3EleMuo->Branch("weight", &m_weight, "weight/f");
  TreeSB3EleMuo->Branch("genEvent", &m_genEvent, "genEvent/f");
  TreeSB3EleMuo->Branch("EleMuo", &m_EleMuo, "EleMuo/i");
  TreeSB3EleMuo->Branch("MuoMuo", &m_MuoMuo, "MuoMuo/i");
  TreeSB3EleMuo->Branch("EleEle", &m_EleEle, "EleEle/i");
  TreeSB3EleMuo->Branch("MuoTau", &m_MuoTau, "MuoTau/i");
  TreeSB3EleMuo->Branch("EleTau", &m_EleTau, "EleTau/i");

  TreeSB3MuoMuo = fs->make<TTree>("TreeSB3MuoMuo", "TreeSB3MuoMuo");
  TreeSB3MuoMuo->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB3MuoMuo->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB3MuoMuo->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB3MuoMuo->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB3MuoMuo->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB3MuoMuo->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB3MuoMuo->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeSB3MuoMuo->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeSB3MuoMuo->Branch("dPhiLep1Met", &m_dPhiLep1Met, "dPhiLep1Met/f");
  TreeSB3MuoMuo->Branch("dRLep1Met", &m_dRLep1Met, "dRLep1Met/f");
  TreeSB3MuoMuo->Branch("dRLep1Lep2", &m_dRLep1Lep2, "dRLep1Lep2/f");
  TreeSB3MuoMuo->Branch("dPhiLep2Met", &m_dPhiLep2Met, "dPhiLep2Met/f");
  TreeSB3MuoMuo->Branch("dRLep2Met", &m_dRLep2Met, "dRLep2Met/f");
  TreeSB3MuoMuo->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB3MuoMuo->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB3MuoMuo->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB3MuoMuo->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB3MuoMuo->Branch("lep1Pt", &m_lep1Pt, "lep1Pt/f");
  TreeSB3MuoMuo->Branch("lep1Eta", &m_lep1Eta, "lep1Eta/f");
  TreeSB3MuoMuo->Branch("lep1Charge", &m_lep1Charge, "lep1Charge/f");
  TreeSB3MuoMuo->Branch("lep1PFIso", &m_lep1PFIso, "lep1PFIso/f");
  TreeSB3MuoMuo->Branch("lep1CorrPFIso", &m_lep1CorrPFIso, "lep1CorrPFIso/f");
  TreeSB3EleMuo->Branch("lep1DETIso", &m_lep1DETIso, "lep1DETIso/f");
  TreeSB3MuoMuo->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB3MuoMuo->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB3MuoMuo->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB3MuoMuo->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB3MuoMuo->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeSB3MuoMuo->Branch("lep2DETIso", &m_lep2DETIso, "lep2DETIso/f");
  TreeSB3MuoMuo->Branch("charge", &m_charge, "charge/f");
  TreeSB3MuoMuo->Branch("met", &m_met, "met/f");
  TreeSB3MuoMuo->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeSB3MuoMuo->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeSB3MuoMuo->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  TreeSB3MuoMuo->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeSB3MuoMuo->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeSB3MuoMuo->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeSB3MuoMuo->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeSB3MuoMuo->Branch("PtSvfit", &m_PtSvfit, "PtSvfit/f");
  TreeSB3MuoMuo->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeSB3MuoMuo->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeSB3MuoMuo->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeSB3MuoMuo->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeSB3MuoMuo->Branch("njet1", &m_njet1, "njet1/i");
  TreeSB3MuoMuo->Branch("nbtagsL1", &m_nbtagsL1, "nbtagsL1/i");
  TreeSB3MuoMuo->Branch("nbtagsM1", &m_nbtagsM1, "nbtagsM1/i");
  TreeSB3MuoMuo->Branch("nbtagsT1", &m_nbtagsT1, "nbtagsT1/i");
  TreeSB3MuoMuo->Branch("njet2", &m_njet2, "njet2/i");
  TreeSB3MuoMuo->Branch("nbtagsL2", &m_nbtagsL2, "nbtagsL2/i");
  TreeSB3MuoMuo->Branch("nbtagsM2", &m_nbtagsM2, "nbtagsM2/i");
  TreeSB3MuoMuo->Branch("nbtagsT2", &m_nbtagsT2, "nbtagsT2/i");
  TreeSB3MuoMuo->Branch("njet3", &m_njet3, "njet3/i");
  TreeSB3MuoMuo->Branch("nbtagsL3", &m_nbtagsL3, "nbtagsL3/i");
  TreeSB3MuoMuo->Branch("nbtagsM3", &m_nbtagsM3, "nbtagsM3/i");
  TreeSB3MuoMuo->Branch("nbtagsT3", &m_nbtagsT3, "nbtagsT3/i");
  TreeSB3MuoMuo->Branch("njet4", &m_njet4, "njet4/i");
  TreeSB3MuoMuo->Branch("nbtagsL4", &m_nbtagsL4, "nbtagsL4/i");
  TreeSB3MuoMuo->Branch("nbtagsM4", &m_nbtagsM4, "nbtagsM4/i");
  TreeSB3MuoMuo->Branch("nbtagsT4", &m_nbtagsT4, "nbtagsT4/i");
  TreeSB3MuoMuo->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeSB3MuoMuo->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeSB3MuoMuo->Branch("trigger", &m_trigger, "trigger/i");
  TreeSB3MuoMuo->Branch("sideband", &m_sideband, "sideband/i");
  TreeSB3MuoMuo->Branch("NVertices", &m_NVertices, "NVertices/i");
  TreeSB3MuoMuo->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeSB3MuoMuo->Branch("metPx", &m_metPx, "metPx/f");
  TreeSB3MuoMuo->Branch("metPy", &m_metPy, "metPy/f");
  TreeSB3MuoMuo->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeSB3MuoMuo->Branch("xsec", &m_xsec, "xsec/d");
  TreeSB3MuoMuo->Branch("lumi", &m_lumi, "lumi/d");
  TreeSB3MuoMuo->Branch("weight", &m_weight, "weight/f");
  TreeSB3MuoMuo->Branch("genEvent", &m_genEvent, "genEvent/f");
  TreeSB3MuoMuo->Branch("EleMuo", &m_EleMuo, "EleMuo/i");
  TreeSB3MuoMuo->Branch("MuoMuo", &m_MuoMuo, "MuoMuo/i");
  TreeSB3MuoMuo->Branch("EleEle", &m_EleEle, "EleEle/i");
  TreeSB3MuoMuo->Branch("MuoTau", &m_MuoTau, "MuoTau/i");
  TreeSB3MuoMuo->Branch("EleTau", &m_EleTau, "EleTau/i");

  TreeSB3EleEle = fs->make<TTree>("TreeSB3EleEle", "TreeSB3EleEle");
  TreeSB3EleEle->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB3EleEle->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB3EleEle->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB3EleEle->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB3EleEle->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB3EleEle->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB3EleEle->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeSB3EleEle->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeSB3EleEle->Branch("dPhiLep1Met", &m_dPhiLep1Met, "dPhiLep1Met/f");
  TreeSB3EleEle->Branch("dRLep1Met", &m_dRLep1Met, "dRLep1Met/f");
  TreeSB3EleEle->Branch("dRLep1Lep2", &m_dRLep1Lep2, "dRLep1Lep2/f");
  TreeSB3EleEle->Branch("dPhiLep2Met", &m_dPhiLep2Met, "dPhiLep2Met/f");
  TreeSB3EleEle->Branch("dRLep2Met", &m_dRLep2Met, "dRLep2Met/f");
  TreeSB3EleEle->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB3EleEle->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB3EleEle->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB3EleEle->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB3EleEle->Branch("lep1Pt", &m_lep1Pt, "lep1Pt/f");
  TreeSB3EleEle->Branch("lep1Eta", &m_lep1Eta, "lep1Eta/f");
  TreeSB3EleEle->Branch("lep1Charge", &m_lep1Charge, "lep1Charge/f");
  TreeSB3EleEle->Branch("lep1PFIso", &m_lep1PFIso, "lep1PFIso/f");
  TreeSB3EleEle->Branch("lep1CorrPFIso", &m_lep1CorrPFIso, "lep1CorrPFIso/f");
  TreeSB3EleEle->Branch("lep1DETIso", &m_lep1DETIso, "lep1DETIso/f");
  TreeSB3EleEle->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB3EleEle->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB3EleEle->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB3EleEle->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB3EleEle->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeSB3EleEle->Branch("lep2DETIso", &m_lep2DETIso, "lep2DETIso/f");
  TreeSB3EleEle->Branch("charge", &m_charge, "charge/f");
  TreeSB3EleEle->Branch("met", &m_met, "met/f");
  TreeSB3EleEle->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeSB3EleEle->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeSB3EleEle->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  TreeSB3EleEle->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeSB3EleEle->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeSB3EleEle->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeSB3EleEle->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeSB3EleEle->Branch("PtSvfit", &m_PtSvfit, "PtSvfit/f");
  TreeSB3EleEle->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeSB3EleEle->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeSB3EleEle->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeSB3EleEle->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeSB3EleEle->Branch("njet1", &m_njet1, "njet1/i");
  TreeSB3EleEle->Branch("nbtagsL1", &m_nbtagsL1, "nbtagsL1/i");
  TreeSB3EleEle->Branch("nbtagsM1", &m_nbtagsM1, "nbtagsM1/i");
  TreeSB3EleEle->Branch("nbtagsT1", &m_nbtagsT1, "nbtagsT1/i");
  TreeSB3EleEle->Branch("njet2", &m_njet2, "njet2/i");
  TreeSB3EleEle->Branch("nbtagsL2", &m_nbtagsL2, "nbtagsL2/i");
  TreeSB3EleEle->Branch("nbtagsM2", &m_nbtagsM2, "nbtagsM2/i");
  TreeSB3EleEle->Branch("nbtagsT2", &m_nbtagsT2, "nbtagsT2/i");
  TreeSB3EleEle->Branch("njet3", &m_njet3, "njet3/i");
  TreeSB3EleEle->Branch("nbtagsL3", &m_nbtagsL3, "nbtagsL3/i");
  TreeSB3EleEle->Branch("nbtagsM3", &m_nbtagsM3, "nbtagsM3/i");
  TreeSB3EleEle->Branch("nbtagsT3", &m_nbtagsT3, "nbtagsT3/i");
  TreeSB3EleEle->Branch("njet4", &m_njet4, "njet4/i");
  TreeSB3EleEle->Branch("nbtagsL4", &m_nbtagsL4, "nbtagsL4/i");
  TreeSB3EleEle->Branch("nbtagsM4", &m_nbtagsM4, "nbtagsM4/i");
  TreeSB3EleEle->Branch("nbtagsT4", &m_nbtagsT4, "nbtagsT4/i");
  TreeSB3EleEle->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeSB3EleEle->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeSB3EleEle->Branch("trigger", &m_trigger, "trigger/i");
  TreeSB3EleEle->Branch("sideband", &m_sideband, "sideband/i");
  TreeSB3EleEle->Branch("NVertices", &m_NVertices, "NVertices/i");
  TreeSB3EleEle->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeSB3EleEle->Branch("metPx", &m_metPx, "metPx/f");
  TreeSB3EleEle->Branch("metPy", &m_metPy, "metPy/f");
  TreeSB3EleEle->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeSB3EleEle->Branch("xsec", &m_xsec, "xsec/d");
  TreeSB3EleEle->Branch("lumi", &m_lumi, "lumi/d");
  TreeSB3EleEle->Branch("weight", &m_weight, "weight/f");
  TreeSB3EleEle->Branch("genEvent", &m_genEvent, "genEvent/f");
  TreeSB3EleEle->Branch("EleMuo", &m_EleMuo, "EleMuo/i");
  TreeSB3EleEle->Branch("MuoMuo", &m_MuoMuo, "MuoMuo/i");
  TreeSB3EleEle->Branch("EleEle", &m_EleEle, "EleEle/i");
  TreeSB3EleEle->Branch("MuoTau", &m_MuoTau, "MuoTau/i");
  TreeSB3EleEle->Branch("EleTau", &m_EleTau, "EleTau/i");

  TreeSB3MuoTau = fs->make<TTree>("TreeSB3MuoTau", "TreeSB3MuoTau");
  TreeSB3MuoTau->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB3MuoTau->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB3MuoTau->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB3MuoTau->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB3MuoTau->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB3MuoTau->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB3MuoTau->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeSB3MuoTau->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeSB3MuoTau->Branch("dPhiLep1Met", &m_dPhiLep1Met, "dPhiLep1Met/f");
  TreeSB3MuoTau->Branch("dRLep1Met", &m_dRLep1Met, "dRLep1Met/f");
  TreeSB3MuoTau->Branch("dRLep1Lep2", &m_dRLep1Lep2, "dRLep1Lep2/f");
  TreeSB3MuoTau->Branch("dPhiLep2Met", &m_dPhiLep2Met, "dPhiLep2Met/f");
  TreeSB3MuoTau->Branch("dRLep2Met", &m_dRLep2Met, "dRLep2Met/f");
  TreeSB3MuoTau->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB3MuoTau->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB3MuoTau->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB3MuoTau->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB3MuoTau->Branch("lep1Pt", &m_lep1Pt, "lep1Pt/f");
  TreeSB3MuoTau->Branch("lep1Eta", &m_lep1Eta, "lep1Eta/f");
  TreeSB3MuoTau->Branch("lep1Charge", &m_lep1Charge, "lep1Charge/f");
  TreeSB3MuoTau->Branch("lep1PFIso", &m_lep1PFIso, "lep1PFIso/f");
  TreeSB3MuoTau->Branch("lep1CorrPFIso", &m_lep1CorrPFIso, "lep1CorrPFIso/f");
  TreeSB3MuoTau->Branch("lep1DETIso", &m_lep1DETIso, "lep1DETIso/f");
  TreeSB3MuoTau->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB3MuoTau->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB3MuoTau->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB3MuoTau->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB3MuoTau->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeSB3MuoTau->Branch("lep2DETIso", &m_lep2DETIso, "lep2DETIso/f");
  TreeSB3MuoTau->Branch("charge", &m_charge, "charge/f");
  TreeSB3MuoTau->Branch("met", &m_met, "met/f");
  TreeSB3MuoTau->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeSB3MuoTau->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeSB3MuoTau->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  TreeSB3MuoTau->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeSB3MuoTau->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeSB3MuoTau->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeSB3MuoTau->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeSB3MuoTau->Branch("PtSvfit", &m_PtSvfit, "PtSvfit/f");
  TreeSB3MuoTau->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeSB3MuoTau->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeSB3MuoTau->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeSB3MuoTau->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeSB3MuoTau->Branch("njet1", &m_njet1, "njet1/i");
  TreeSB3MuoTau->Branch("nbtagsL1", &m_nbtagsL1, "nbtagsL1/i");
  TreeSB3MuoTau->Branch("nbtagsM1", &m_nbtagsM1, "nbtagsM1/i");
  TreeSB3MuoTau->Branch("nbtagsT1", &m_nbtagsT1, "nbtagsT1/i");
  TreeSB3MuoTau->Branch("njet2", &m_njet2, "njet2/i");
  TreeSB3MuoTau->Branch("nbtagsL2", &m_nbtagsL2, "nbtagsL2/i");
  TreeSB3MuoTau->Branch("nbtagsM2", &m_nbtagsM2, "nbtagsM2/i");
  TreeSB3MuoTau->Branch("nbtagsT2", &m_nbtagsT2, "nbtagsT2/i");
  TreeSB3MuoTau->Branch("njet3", &m_njet3, "njet3/i");
  TreeSB3MuoTau->Branch("nbtagsL3", &m_nbtagsL3, "nbtagsL3/i");
  TreeSB3MuoTau->Branch("nbtagsM3", &m_nbtagsM3, "nbtagsM3/i");
  TreeSB3MuoTau->Branch("nbtagsT3", &m_nbtagsT3, "nbtagsT3/i");
  TreeSB3MuoTau->Branch("njet4", &m_njet4, "njet4/i");
  TreeSB3MuoTau->Branch("nbtagsL4", &m_nbtagsL4, "nbtagsL4/i");
  TreeSB3MuoTau->Branch("nbtagsM4", &m_nbtagsM4, "nbtagsM4/i");
  TreeSB3MuoTau->Branch("nbtagsT4", &m_nbtagsT4, "nbtagsT4/i");
  TreeSB3MuoTau->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeSB3MuoTau->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeSB3MuoTau->Branch("trigger", &m_trigger, "trigger/i");
  TreeSB3MuoTau->Branch("sideband", &m_sideband, "sideband/i");
  TreeSB3MuoTau->Branch("NVertices", &m_NVertices, "NVertices/i");
  TreeSB3MuoTau->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeSB3MuoTau->Branch("metPx", &m_metPx, "metPx/f");
  TreeSB3MuoTau->Branch("metPy", &m_metPy, "metPy/f");
  TreeSB3MuoTau->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeSB3MuoTau->Branch("xsec", &m_xsec, "xsec/d");
  TreeSB3MuoTau->Branch("lumi", &m_lumi, "lumi/d");
  TreeSB3MuoTau->Branch("weight", &m_weight, "weight/f");
  TreeSB3MuoTau->Branch("genEvent", &m_genEvent, "genEvent/f");
  TreeSB3MuoTau->Branch("EleMuo", &m_EleMuo, "EleMuo/i");
  TreeSB3MuoTau->Branch("MuoMuo", &m_MuoMuo, "MuoMuo/i");
  TreeSB3MuoTau->Branch("EleEle", &m_EleEle, "EleEle/i");
  TreeSB3MuoTau->Branch("MuoTau", &m_MuoTau, "MuoTau/i");
  TreeSB3MuoTau->Branch("EleTau", &m_EleTau, "EleTau/i");

  TreeSB3EleTau = fs->make<TTree>("TreeSB3EleTau", "TreeSB3EleTau");
  TreeSB3EleTau->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB3EleTau->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB3EleTau->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB3EleTau->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB3EleTau->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB3EleTau->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB3EleTau->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeSB3EleTau->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeSB3EleTau->Branch("dPhiLep1Met", &m_dPhiLep1Met, "dPhiLep1Met/f");
  TreeSB3EleTau->Branch("dRLep1Met", &m_dRLep1Met, "dRLep1Met/f");
  TreeSB3EleTau->Branch("dRLep1Lep2", &m_dRLep1Lep2, "dRLep1Lep2/f");
  TreeSB3EleTau->Branch("dPhiLep2Met", &m_dPhiLep2Met, "dPhiLep2Met/f");
  TreeSB3EleTau->Branch("dRLep2Met", &m_dRLep2Met, "dRLep2Met/f");
  TreeSB3EleTau->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB3EleTau->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB3EleTau->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB3EleTau->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB3EleTau->Branch("lep1Pt", &m_lep1Pt, "lep1Pt/f");
  TreeSB3EleTau->Branch("lep1Eta", &m_lep1Eta, "lep1Eta/f");
  TreeSB3EleTau->Branch("lep1Charge", &m_lep1Charge, "lep1Charge/f");
  TreeSB3EleTau->Branch("lep1PFIso", &m_lep1PFIso, "lep1PFIso/f");
  TreeSB3EleTau->Branch("lep1CorrPFIso", &m_lep1CorrPFIso, "lep1CorrPFIso/f");
  TreeSB3EleTau->Branch("lep1DETIso", &m_lep1DETIso, "lep1DETIso/f");
  TreeSB3EleTau->Branch("lep2Pt", &m_lep2Pt, "lep2Pt/f");
  TreeSB3EleTau->Branch("lep2Eta", &m_lep2Eta, "lep2Eta/f");
  TreeSB3EleTau->Branch("lep2Charge", &m_lep2Charge, "lep2Charge/f");
  TreeSB3EleTau->Branch("lep2PFIso", &m_lep2PFIso, "lep2PFIso/f");
  TreeSB3EleTau->Branch("lep2CorrPFIso", &m_lep2CorrPFIso, "lep2CorrPFIso/f");
  TreeSB3EleTau->Branch("lep2DETIso", &m_lep2DETIso, "lep2DETIso/f");
  TreeSB3EleTau->Branch("charge", &m_charge, "charge/f");
  TreeSB3EleTau->Branch("met", &m_met, "met/f");
  TreeSB3EleTau->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeSB3EleTau->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeSB3EleTau->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  TreeSB3EleTau->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeSB3EleTau->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeSB3EleTau->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeSB3EleTau->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeSB3EleTau->Branch("PtSvfit", &m_PtSvfit, "PtSvfit/f");
  TreeSB3EleTau->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeSB3EleTau->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeSB3EleTau->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeSB3EleTau->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeSB3EleTau->Branch("njet1", &m_njet1, "njet1/i");
  TreeSB3EleTau->Branch("nbtagsL1", &m_nbtagsL1, "nbtagsL1/i");
  TreeSB3EleTau->Branch("nbtagsM1", &m_nbtagsM1, "nbtagsM1/i");
  TreeSB3EleTau->Branch("nbtagsT1", &m_nbtagsT1, "nbtagsT1/i");
  TreeSB3EleTau->Branch("njet2", &m_njet2, "njet2/i");
  TreeSB3EleTau->Branch("nbtagsL2", &m_nbtagsL2, "nbtagsL2/i");
  TreeSB3EleTau->Branch("nbtagsM2", &m_nbtagsM2, "nbtagsM2/i");
  TreeSB3EleTau->Branch("nbtagsT2", &m_nbtagsT2, "nbtagsT2/i");
  TreeSB3EleTau->Branch("njet3", &m_njet3, "njet3/i");
  TreeSB3EleTau->Branch("nbtagsL3", &m_nbtagsL3, "nbtagsL3/i");
  TreeSB3EleTau->Branch("nbtagsM3", &m_nbtagsM3, "nbtagsM3/i");
  TreeSB3EleTau->Branch("nbtagsT3", &m_nbtagsT3, "nbtagsT3/i");
  TreeSB3EleTau->Branch("njet4", &m_njet4, "njet4/i");
  TreeSB3EleTau->Branch("nbtagsL4", &m_nbtagsL4, "nbtagsL4/i");
  TreeSB3EleTau->Branch("nbtagsM4", &m_nbtagsM4, "nbtagsM4/i");
  TreeSB3EleTau->Branch("nbtagsT4", &m_nbtagsT4, "nbtagsT4/i");
  TreeSB3EleTau->Branch("trigger320", &m_trigger320, "trigger320/i");
  TreeSB3EleTau->Branch("trigger650", &m_trigger650, "trigger650/i");
  TreeSB3EleTau->Branch("trigger", &m_trigger, "trigger/i");
  TreeSB3EleTau->Branch("sideband", &m_sideband, "sideband/i");
  TreeSB3EleTau->Branch("NVertices", &m_NVertices, "NVertices/i");
  TreeSB3EleTau->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeSB3EleTau->Branch("metPx", &m_metPx, "metPx/f");
  TreeSB3EleTau->Branch("metPy", &m_metPy, "metPy/f");
  TreeSB3EleTau->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  TreeSB3EleTau->Branch("xsec", &m_xsec, "xsec/d");
  TreeSB3EleTau->Branch("lumi", &m_lumi, "lumi/d");
  TreeSB3EleTau->Branch("weight", &m_weight, "weight/f");
  TreeSB3EleTau->Branch("genEvent", &m_genEvent, "genEvent/f");
  TreeSB3EleTau->Branch("EleMuo", &m_EleMuo, "EleMuo/i");
  TreeSB3EleTau->Branch("MuoMuo", &m_MuoMuo, "MuoMuo/i");
  TreeSB3EleTau->Branch("EleEle", &m_EleEle, "EleEle/i");
  TreeSB3EleTau->Branch("MuoTau", &m_MuoTau, "MuoTau/i");
  TreeSB3EleTau->Branch("EleTau", &m_EleTau, "EleTau/i");
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
			 bool & foundJet, pat::JetCollection::const_iterator & SelectedJet,
			 float & prunedMass, float massMin, float massMax, bool tau21DOWN){

  float ptZ=-99;
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
    if(fabs(jet->eta())>1.0 && fabs(jet->eta())<1.5 && jet->neutralMultiplicity()!=0){ if(jet->chargedMultiplicity()/jet->neutralMultiplicity()>2) continue;}
    if(jet->nConstituents()<=1) continue;
    if(jet->pt()<400) continue;
    if(fabs(jet->eta())>2.4) continue;
    if(!(mass>massMin && mass<massMax))  continue;
    if(tau21DOWN){  if(jet->userFloat("tau2")/jet->userFloat("tau1")>0.75) continue;}
    if(!tau21DOWN){ if(jet->userFloat("tau2")/jet->userFloat("tau1")<0.75) continue;}
    foundJet=true;
    if(jet->pt()>ptZ){
      prunedMass=mass;
      ptZ=jet->pt();
      SelectedJet=jet;
    }
  }
}

void Analyzer::SelectTau(edm::Handle<pat::TauCollection> tauHandle,
			 pat::JetCollection::const_iterator SelectedJet,
			 std::vector<pat::TauCollection::const_iterator> & SelectedTau,
			 bool foundJet){
  
  for (pat::TauCollection::const_iterator patTau = tauHandle->begin(); patTau != tauHandle->end(); ++patTau ) {
    if(patTau->pt()<20) continue;
    if(fabs(patTau->eta())>2.3) continue;
    if(patTau->tauID("decayModeFindingNewDMs")<0.5) continue;
    if(patTau->tauID("againstMuonLoose")<0.5) continue;
    if(patTau->tauID("againstElectronLoose")<0.5) continue;
    if(patTau->tauID("byVLooseIsolationMVA3newDMwoLT")<0.5) continue;
    if(foundJet){
      if(ROOT::Math::VectorUtil::DeltaR(patTau->p4(),SelectedJet->p4())>0.8) SelectedTau.push_back(patTau);
    }
  }
}

void Analyzer::SelectMuon(edm::Handle<pat::MuonCollection> muoH,
			  pat::JetCollection::const_iterator SelectedJet,
			  std::vector<pat::MuonCollection::const_iterator> & SelectedMuon,
			  bool foundJet, reco::Vertex primaryVertex, bool fully){
  for(pat::MuonCollection::const_iterator muon = muoH->begin(); muon != muoH->end(); ++muon) {
    if(!(muon->isGlobalMuon())) continue;
    reco::TrackRef cktTrack = (muon::tevOptimized(*muon, 200, 17., 40., 0.25)).first;
    if(cktTrack->pt()<10) continue;
    if(fabs(cktTrack->eta())>2.4) continue;
    if(fabs(cktTrack->phi())>3.2) continue;
    if((cktTrack->ptError()/cktTrack->pt())>0.3) continue;
    if(muon->globalTrack()->hitPattern().numberOfValidMuonHits()<=0) continue;
    if(muon->numberOfMatches()<=1) continue;
    if(fabs(cktTrack->dxy(primaryVertex.position()))>=0.2) continue;
    if(fabs(cktTrack->dz( primaryVertex.position()))>=0.5) continue;
    if(muon->innerTrack()->hitPattern().numberOfValidPixelHits()<=0) continue;
    if(muon->innerTrack()->hitPattern().trackerLayersWithMeasurement()<=5) continue;
    if(fully){if(MuonPFIso(muon, true)>0.2) continue;}
    else     {if(MuonCorrPFIso(muon, true)>0.2) continue;}
    if(foundJet) {
      if(ROOT::Math::VectorUtil::DeltaR(muon->p4(),SelectedJet->p4())>0.8) SelectedMuon.push_back(muon);
    }
  }
}

void Analyzer::SelectElectron(edm::Handle<pat::ElectronCollection> eleH,
			      pat::JetCollection::const_iterator SelectedJet,
			      std::vector<pat::ElectronCollection::const_iterator> & SelectedElectron,
			      bool foundJet, reco::Vertex primaryVertex, float rho, bool fully){
  bool passEle=false;
  for(pat::ElectronCollection::const_iterator electron = eleH->begin(); electron != eleH->end(); ++electron) {
    if(electron->pt()<10) continue;
    if(fully){
      if(ElectronPFIso(electron,rho)>0.1) continue;
      if(electron->pt()<20){if(ElectronPFIso(electron,rho)>0.07) continue;}
    } else {
      if(ElectronCorrPFIso(electron,rho)>0.1) continue;
      if(electron->pt()<20){if(ElectronCorrPFIso(electron,rho)>0.07) continue;}
    }
    if(fabs(electron->superCluster()->eta())<=1.479){
      if(fabs(electron->deltaEtaSuperClusterTrackAtVtx())>=0.004) continue;
      if(fabs(electron->deltaPhiSuperClusterTrackAtVtx())>=0.030) continue;
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
      if(fabs(electron->deltaEtaSuperClusterTrackAtVtx())>=0.005) continue;
      if(fabs(electron->deltaPhiSuperClusterTrackAtVtx())>=0.020) continue;
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
      if(ROOT::Math::VectorUtil::DeltaR(electron->p4(),SelectedJet->p4())>0.8) SelectedElectron.push_back(electron);
    }
  }
}

void Analyzer::SelectEM(pat::ElectronCollection::const_iterator & SelectedElectron, pat::MuonCollection::const_iterator & SelectedMuon, bool & foundEM,
			std::vector<pat::ElectronCollection::const_iterator> SelectedElectrons, std::vector<pat::MuonCollection::const_iterator> SelectedMuons){
  float pt1=-99; float pt2=-99;
  for(unsigned int i=0; i<SelectedElectrons.size(); i++){
    for(unsigned int j=0; j<SelectedMuons.size(); j++){
      math::PtEtaPhiELorentzVector dilep; dilep = SelectedElectrons[i]->p4() + SelectedMuons[j]->p4();
      if(ROOT::Math::VectorUtil::DeltaR(SelectedElectrons[i]->p4(),SelectedMuons[j]->p4())>0.1 && dilep.mass()>10){
	foundEM=true;
	if(SelectedElectrons[i]->pt()>pt1){
	  pt1=SelectedElectrons[i]->pt();
	  SelectedElectron=SelectedElectrons[i];
	}
	if(SelectedMuons[j]->pt()>pt2){
	  pt2=SelectedMuons[j]->pt();
	  SelectedMuon=SelectedMuons[j];
	}
      }
    }
  }
}

void Analyzer::SelectMM(pat::MuonCollection::const_iterator & SelectedMuon1, 
			    pat::MuonCollection::const_iterator & SelectedMuon2, 
			    bool & foundMM, std::vector<pat::MuonCollection::const_iterator> SelectedTrackerMuo, 
			    pat::JetCollection::const_iterator SelectedJet, bool foundJet, reco::Vertex primaryVertex){
  vector<pat::MuonCollection::const_iterator> SelectedMuo;
  for(unsigned int i=0; i<SelectedTrackerMuo.size(); i++){
    if(!foundJet) continue;
    if((MuonDETIso(SelectedTrackerMuo[i],SelectedTrackerMuo))>0.1) continue;
    if(ROOT::Math::VectorUtil::DeltaR(SelectedTrackerMuo[i]->p4(),SelectedJet->p4())<0.8) continue;
    SelectedMuo.push_back(SelectedTrackerMuo[i]);
  }
  vector<pat::MuonCollection::const_iterator> SelectedGlobalMuo;
  for(unsigned int i=0; i<SelectedMuo.size(); i++){
    if(!(SelectedMuo[i]->isGlobalMuon())) continue;
    reco::TrackRef cktTrack = (muon::tevOptimized(*(SelectedMuo[i]), 200, 17., 40., 0.25)).first;
    if(cktTrack->pt()<10) continue;
    if(fabs(cktTrack->eta())>2.4) continue;
    if(fabs(cktTrack->phi())>3.2) continue;
    if((cktTrack->ptError()/cktTrack->pt())>0.3) continue;
    if(SelectedMuo[i]->globalTrack()->hitPattern().numberOfValidMuonHits()<=0) continue;
    if(SelectedMuo[i]->numberOfMatches()<=1) continue;
    if(fabs(cktTrack->dxy(primaryVertex.position()))>=0.2) continue;
    if(fabs(cktTrack->dz( primaryVertex.position()))>=0.5) continue;
    if(SelectedMuo[i]->innerTrack()->hitPattern().numberOfValidPixelHits()<=0) continue;
    if(SelectedMuo[i]->innerTrack()->hitPattern().trackerLayersWithMeasurement()<=5) continue;
    SelectedGlobalMuo.push_back(SelectedMuo[i]);
  }
  float pt=-99;
  pat::MuonCollection::const_iterator Muon1;
  pat::MuonCollection::const_iterator Muon2;
  for(unsigned int i=0; i<SelectedGlobalMuo.size(); i++){
    for(unsigned int j=0; j<SelectedMuo.size(); j++){
      math::PtEtaPhiELorentzVector dilep; dilep = SelectedGlobalMuo[i]->p4() + SelectedMuo[j]->p4();
      if(ROOT::Math::VectorUtil::DeltaR(SelectedGlobalMuo[i]->p4(),SelectedMuo[j]->p4())>0.1 && dilep.mass()>10){
	foundMM=true;
	if(pt<SelectedGlobalMuo[i]->pt()+SelectedMuo[j]->pt()){
	  pt=SelectedGlobalMuo[i]->pt()+SelectedMuo[j]->pt();
	  Muon1=SelectedGlobalMuo[i];
	  Muon2=SelectedMuo[j];
	}
      }
    }
  }
  if(foundMM){
    if(Muon1->pt()>Muon2->pt()){
      SelectedMuon1=Muon1;
      SelectedMuon2=Muon2;
    }else{
      SelectedMuon1=Muon2;
      SelectedMuon2=Muon1;
    }
  }
}

void Analyzer::SelectEE(pat::ElectronCollection::const_iterator & SelectedElectron1, pat::ElectronCollection::const_iterator & SelectedElectron2, bool & foundEE,
			std::vector<pat::ElectronCollection::const_iterator> SelectedElectrons){
  float pt=-99;
  pat::ElectronCollection::const_iterator Electron1;
  pat::ElectronCollection::const_iterator Electron2;
  for(unsigned int i=0; i<SelectedElectrons.size(); i++){
    for(unsigned int j=0; j<SelectedElectrons.size(); j++){
      math::PtEtaPhiELorentzVector dilep; dilep = SelectedElectrons[i]->p4() + SelectedElectrons[j]->p4();
      if(ROOT::Math::VectorUtil::DeltaR(SelectedElectrons[i]->p4(),SelectedElectrons[j]->p4())>0.1 && dilep.mass()>10){
	foundEE=true;
	if(pt<SelectedElectrons[i]->pt()+SelectedElectrons[j]->pt()){
	  pt=SelectedElectrons[i]->pt()+SelectedElectrons[j]->pt();
	  Electron1=SelectedElectrons[i];
	  Electron2=SelectedElectrons[j];
	}
      }
    }
  }
  if(foundEE){
    if(Electron1->pt()>Electron2->pt()){
      SelectedElectron1=Electron1;
      SelectedElectron2=Electron2;
    }else{
      SelectedElectron1=Electron2;
      SelectedElectron2=Electron1;
    }
  }
}

void Analyzer::SelectMT(pat::MuonCollection::const_iterator & SelectedMuon, pat::TauCollection::const_iterator & SelectedTau, bool & foundMT,
			std::vector<pat::MuonCollection::const_iterator> SelectedMuons, std::vector<pat::TauCollection::const_iterator> SelectedTaus){
  float pt1=-99; float pt2=-99;
  for(unsigned int i=0; i<SelectedMuons.size(); i++){
    for(unsigned int j=0; j<SelectedTaus.size(); j++){
      math::PtEtaPhiELorentzVector dilep; dilep = SelectedMuons[i]->p4() + SelectedTaus[j]->p4();
      if(ROOT::Math::VectorUtil::DeltaR(SelectedMuons[i]->p4(),SelectedTaus[j]->p4())>0.1 && dilep.mass()>10){
	foundMT=true;
	if(SelectedMuons[i]->pt()>pt1){
	  pt1=SelectedMuons[i]->pt();
	  SelectedMuon=SelectedMuons[i];
	}
	if(SelectedTaus[j]->pt()>pt2){
	  pt2=SelectedTaus[j]->pt();
	  SelectedTau=SelectedTaus[j];
	}
      }
    }
  }
}

void Analyzer::SelectET(pat::ElectronCollection::const_iterator & SelectedElectron, pat::TauCollection::const_iterator & SelectedTau, bool & foundET,
			std::vector<pat::ElectronCollection::const_iterator> SelectedElectrons, std::vector<pat::TauCollection::const_iterator> SelectedTaus){
  float pt1=-99; float pt2=-99;
  for(unsigned int i=0; i<SelectedElectrons.size(); i++){
    for(unsigned int j=0; j<SelectedTaus.size(); j++){
      math::PtEtaPhiELorentzVector dilep; dilep = SelectedElectrons[i]->p4() + SelectedTaus[j]->p4();
      if(ROOT::Math::VectorUtil::DeltaR(SelectedElectrons[i]->p4(),SelectedTaus[j]->p4())>0.1 && dilep.mass()>10){
	foundET=true;
	if(SelectedElectrons[i]->pt()>pt1){
	  pt1=SelectedElectrons[i]->pt();
	  SelectedElectron=SelectedElectrons[i];
	}
	if(SelectedTaus[j]->pt()>pt2){
	  pt2=SelectedTaus[j]->pt();
	  SelectedTau=SelectedTaus[j];
	}
      }
    }
  }
}


void Analyzer::SelectTrackerMuon(edm::Handle<pat::MuonCollection> muoH, vector<pat::MuonCollection::const_iterator> & SelectedMuo, reco::Vertex primaryVertex){
  for(pat::MuonCollection::const_iterator muon = muoH->begin(); muon != muoH->end(); ++muon) {
    if(muon->pt()<10) continue;
    if(fabs(muon->eta())>2.4) continue;
    if(fabs(muon->phi())>3.2) continue;
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

void Analyzer::BtagVeto(int & njet1, int & nbtagsL1, int & nbtagsM1, int & nbtagsT1, int & njet2, int & nbtagsL2, int & nbtagsM2, int & nbtagsT2,
			int & njet3, int & nbtagsL3, int & nbtagsM3, int & nbtagsT3, int & njet4, int & nbtagsL4, int & nbtagsM4, int & nbtagsT4,
			pat::JetCollection::const_iterator SelectedJet, math::PtEtaPhiELorentzVector lep1, math::PtEtaPhiELorentzVector lep2, const edm::Event& iEvent){

  //edm::Handle<pat::JetCollection> ak5jetCands;
  edm::Handle<edm::View<pat::Jet> > ak5jetCands;
  iEvent.getByLabel(ak5JetColl_,ak5jetCands);
  edm::Handle<ValueMap<float> > puJetIdMVA;
  iEvent.getByLabel("puJetMvaAK5CHS","fullDiscriminant", puJetIdMVA);
  edm::Handle<ValueMap<int> > puJetIdFlag;
  iEvent.getByLabel("puJetMvaAK5CHS","fullId", puJetIdFlag);

  for ( unsigned int i=0; i<ak5jetCands->size(); ++i ) {
    const pat::Jet & ak5 = ak5jetCands->at(i);
    //for(pat::JetCollection::const_iterator ak5 = ak5jetCands->begin(); ak5 != ak5jetCands->end(); ++ak5) {
    if(ak5.pt()<20) continue;
    if(fabs(ak5.eta())>2.4) continue;
    if(ROOT::Math::VectorUtil::DeltaR(ak5.p4(),SelectedJet->p4())<0.8) continue;
    if(ROOT::Math::VectorUtil::DeltaR(ak5.p4(),lep1)<0.5) continue;
    if(ROOT::Math::VectorUtil::DeltaR(ak5.p4(),lep2)<0.5) continue;
    //float mva   = (*puJetIdMVA)[ak5jetCands->refAt(i)];
    int idflag = (*puJetIdFlag)[ak5jetCands->refAt(i)];
    double discCSV = ak5.bDiscriminator("combinedSecondaryVertexBJetTags");
    njet1++;
    if(discCSV>0.244) nbtagsL1++; //loose working point
    if(discCSV>0.679) nbtagsM1++; //medium working point
    if(discCSV>0.898) nbtagsT1++; //tight working point
    if(!(PileupJetIdentifier::passJetId(idflag, PileupJetIdentifier::kLoose))){
      njet2++;
      if(discCSV>0.244) nbtagsL2++;
      if(discCSV>0.679) nbtagsM2++;
      if(discCSV>0.898) nbtagsT2++;
    }
    if(!(PileupJetIdentifier::passJetId(idflag, PileupJetIdentifier::kMedium))){
      njet3++;
      if(discCSV>0.244) nbtagsL3++;
      if(discCSV>0.679) nbtagsM3++;
      if(discCSV>0.898) nbtagsT3++;
    }
    if(!(PileupJetIdentifier::passJetId(idflag, PileupJetIdentifier::kTight))){
      njet4++;
      if(discCSV>0.244) nbtagsL4++;
      if(discCSV>0.679) nbtagsM4++;
      if(discCSV>0.898) nbtagsT4++;
    }
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
			std::vector<pat::MuonCollection::const_iterator> SelectedTrackerMuo,
			pat::ElectronCollection::const_iterator SelectedEle, pat::ElectronCollection::const_iterator SelectedEle1, 
			pat::ElectronCollection::const_iterator SelectedEle2, edm::Handle<reco::VertexCollection> vertices,
			edm::Handle<pat::METCollection> metRaw, edm::Handle<pat::METCollection> met, edm::Handle<pat::METCollection> uncorrmet,
			float prunedMass, bool isFired_HLT, bool isFired_HLT_PFJet320, bool isFired_HLT_HT650,
			double MyWeight, float genEvent, float rho, int EleMuo, int MuoMuo, int EleEle, int MuoTau, int EleTau, const edm::Event& iEvent){
  
  math::PtEtaPhiELorentzVector lep1;
  math::PtEtaPhiELorentzVector lep2;
  LorentzVector lep1SVFit;
  LorentzVector lep2SVFit;
  float lep1Charge=0;  float lep2Charge=0;
  float lep1PFIso=100; float lep1CorrPFIso=100; float lep1DETIso=100;
  float lep2PFIso=100; float lep2CorrPFIso=100;	float lep2DETIso=100;
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
    lep2DETIso=MuonDETIso(SelectedMuo,SelectedTrackerMuo);
  }else if(category==1){
    lep1 = SelectedMuo1->p4();
    lep2 = SelectedMuo2->p4();
    lep1SVFit = SelectedMuo1->p4();
    lep2SVFit = SelectedMuo2->p4();
    lep1Charge=SelectedMuo1->charge();
    lep2Charge=SelectedMuo2->charge();
    lep1PFIso=MuonPFIso(SelectedMuo1,true);
    lep1CorrPFIso=MuonCorrPFIso(SelectedMuo1,true);
    lep1DETIso=MuonDETIso(SelectedMuo1,SelectedTrackerMuo);
    lep2PFIso=MuonPFIso(SelectedMuo2,true);
    lep2CorrPFIso=MuonCorrPFIso(SelectedMuo2,true);
    lep2DETIso=MuonDETIso(SelectedMuo2,SelectedTrackerMuo);
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
    lep2DETIso=MuonDETIso(SelectedMuo,SelectedTrackerMuo);
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

  //BTAG VETO
  int njet1=0; int nbtagsL1=0; int nbtagsM1=0; int nbtagsT1=0; int njet2=0; int nbtagsL2=0; int nbtagsM2=0; int nbtagsT2=0;
  int njet3=0; int nbtagsL3=0; int nbtagsM3=0; int nbtagsT3=0; int njet4=0; int nbtagsL4=0; int nbtagsM4=0; int nbtagsT4=0;
  BtagVeto(njet1, nbtagsL1, nbtagsM1, nbtagsT1, njet2, nbtagsL2, nbtagsM2, nbtagsT2, 
	   njet3, nbtagsL3, nbtagsM3, nbtagsT3, njet4, nbtagsL4, nbtagsM4, nbtagsT4, SelectedJet, lep1, lep2, iEvent);

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
  if(((1+a)*(1+b))>0 && ((1+a)*(1+b))<9999999999) {
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

  if(prunedMass>70 && prunedMass<110){
    //if(category==0) cout<<"EleMuo "<<iEvent.id().event()<<"; MassSVFit "<<MassSVFit<<"; ptSVFit "<<ptSVFit<<"; XMassSVFit "<<XMassSVFit<<endl;
    //if(category==1) cout<<"MuoMuo "<<iEvent.id().event()<<"; MassSVFit "<<MassSVFit<<"; ptSVFit "<<ptSVFit<<"; XMassSVFit "<<XMassSVFit<<endl;
    //if(category==2) cout<<"EleEle "<<iEvent.id().event()<<"; MassSVFit "<<MassSVFit<<"; ptSVFit "<<ptSVFit<<"; XMassSVFit "<<XMassSVFit<<endl;
    //if(category==3) cout<<"MuoTau "<<iEvent.id().event()<<"; MassSVFit "<<MassSVFit<<"; ptSVFit "<<ptSVFit<<"; XMassSVFit "<<XMassSVFit<<endl;
    //if(category==4) cout<<"EleTau "<<iEvent.id().event()<<"; MassSVFit "<<MassSVFit<<"; ptSVFit "<<ptSVFit<<"; XMassSVFit "<<XMassSVFit<<endl;
  }

  m_jetPt=SelectedJet->pt();
  m_jetEta=SelectedJet->eta();
  m_jetMass=prunedMass;
  m_jetSubjettiness=SelectedJet->userFloat("tau2")/SelectedJet->userFloat("tau1");
  m_lep1Pt=lep1.pt();
  m_lep1Eta=lep1.eta();
  m_lep1Charge=lep1Charge;
  m_lep1PFIso=lep1PFIso;
  m_lep1CorrPFIso=lep1CorrPFIso;
  m_lep1DETIso=lep1DETIso;
  m_lep2Pt=lep2.pt();
  m_lep2Eta=lep2.eta();
  m_lep2Charge=lep2Charge;
  m_lep2PFIso=lep2PFIso;
  m_lep2CorrPFIso=lep2CorrPFIso;
  m_lep2DETIso=lep2DETIso;
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
  m_njet1=njet1;
  m_nbtagsL1=nbtagsL1;
  m_nbtagsM1=nbtagsM1;
  m_nbtagsT1=nbtagsT1;
  m_njet2=njet2;
  m_nbtagsL2=nbtagsL2;
  m_nbtagsM2=nbtagsM2;
  m_nbtagsT2=nbtagsT2;
  m_njet3=njet3;
  m_nbtagsL3=nbtagsL3;
  m_nbtagsM3=nbtagsM3;
  m_nbtagsT3=nbtagsT3;
  m_njet4=njet4;
  m_nbtagsL4=nbtagsL4;
  m_nbtagsM4=nbtagsM4;
  m_nbtagsT4=nbtagsT4;
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
