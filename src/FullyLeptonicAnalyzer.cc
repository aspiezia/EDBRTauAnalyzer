// -*- C++ -*-
//
// Package: FullyLeptonicAnalyzer
// Class: FullyLeptonicAnalyzer
//
/**\class FullyLeptonicAnalyzer FullyLeptonicAnalyzer.cc ExoDiBosonResonances/EDBRTauAnalyzer/src/FullyLeptonicAnalyzer.cc

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

class FullyLeptonicAnalyzer : public edm::EDAnalyzer {
public:
  explicit FullyLeptonicAnalyzer(const edm::ParameterSet&);
  ~FullyLeptonicAnalyzer();
  
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
  void SelectJet(edm::Handle<pat::JetCollection> CA8JetswithQjets, edm::Handle<pat::JetCollection> CA8JetsPruned,
                 std::vector<pat::JetCollection::const_iterator> & SelectedJet, std::vector<float> & SelectedPrunedMass, int & Njet);
  void SelectElectron(edm::Handle<pat::ElectronCollection> eleH, std::vector<pat::ElectronCollection::const_iterator> & SelectedEle, int & Nele, float rho);
  void SelectElectronCutBased(edm::Handle<pat::ElectronCollection> eleH, std::vector<pat::ElectronCollection::const_iterator> & SelectedEle, int & Nele, float rho, reco::Vertex primaryVertex);
  void SelectTightMuon( edm::Handle<pat::MuonCollection> muoH, std::vector<pat::MuonCollection::const_iterator> & SelectedMuo, int & Nmuo, reco::Vertex primaryVertex);
  void SelectHighptMuon( edm::Handle<pat::MuonCollection> muoH, std::vector<pat::MuonCollection::const_iterator> & SelectedMuo, int & Nmuo, reco::Vertex primaryVertex);
  void SelectTrackerMuon(edm::Handle<pat::MuonCollection> muoH, std::vector<pat::MuonCollection::const_iterator> & SelectedMuo, int & Nmuo, reco::Vertex primaryVertex);
  void SelectTrackerGlobalID(std::vector<pat::MuonCollection::const_iterator> SelectedMuo, reco::Vertex primaryVertex, bool & hasAtLeastOneHighPtMuo);
  float ElectronPFIso(pat::ElectronCollection::const_iterator electron, float rho);
  bool ElectronDETIso(pat::ElectronCollection::const_iterator electron, float rho);
  float MuonPFIso(pat::MuonCollection::const_iterator muon, bool highpt);
  float MuonDETIso(pat::MuonCollection::const_iterator SelectedMuo, std::vector<pat::MuonCollection::const_iterator> SelectedHighptMuo);
  void BtagVeto(edm::Handle<pat::JetCollection> ak5jetCands, int & nbtagsL, int & nbtagsM, int & nbtagsT, pat::JetCollection::const_iterator SelectedJet);
  float SVFitMass(edm::Handle<pat::METCollection> met, edm::Handle<pat::METCollection> metRaw, TLorentzVector & SVFitTauTau, LorentzVector lep1, LorentzVector lep2);
  void CollinearApproximation(edm::Handle<pat::METCollection> met, TLorentzVector & CATauTau, TLorentzVector PrunedJet, LorentzVector lep1, LorentzVector lep2, bool & CA);
 
  TH1D* Nevents;

  //TREE
  TTree *TreeEleEle;                     TTree *TreeEleMuo;		      TTree *TreeMuoMuo;
  float m_jetPtEleEle;			 float m_jetPtEleMuo;		      float m_jetPtMuoMuo;
  float m_jetEtaEleEle;			 float m_jetEtaEleMuo;		      float m_jetEtaMuoMuo;
  float m_jetMassEleEle;		 float m_jetMassEleMuo;	      	      float m_jetMassMuoMuo;
  float m_jetSubjettinessEleEle;	 float m_jetSubjettinessEleMuo;       float m_jetSubjettinessMuoMuo;
  float m_dRJetLep1EleEle;		 float m_dRJetLep1EleMuo;	      float m_dRJetLep1MuoMuo;
  float m_dRJetLep2EleEle;		 float m_dRJetLep2EleMuo;	      float m_dRJetLep2MuoMuo;
  float m_dRJetMetEleEle;		 float m_dRJetMetEleMuo;	      float m_dRJetMetMuoMuo;
  float m_dPhiJetMetEleEle;		 float m_dPhiJetMetEleMuo;	      float m_dPhiJetMetMuoMuo;
  float m_dRLepMet1EleEle;		 float m_dRLepMet1EleMuo;	      float m_dRLepMet1MuoMuo;
  float m_dPhiLepMet1EleEle;		 float m_dPhiLepMet1EleMuo;	      float m_dPhiLepMet1MuoMuo;
  float m_dRLepMet2EleEle;               float m_dRLepMet2EleMuo;             float m_dRLepMet2MuoMuo;
  float m_dPhiLepMet2EleEle;     	 float m_dPhiLepMet2EleMuo;           float m_dPhiLepMet2MuoMuo;
  float m_dRLepLepEleEle;                float m_dRLepLepEleMuo;              float m_dRLepLepMuoMuo;
  float m_dRZZVisEleEle;                 float m_dRZZVisEleMuo;		      float m_dRZZVisMuoMuo;
  float m_dRZZEffEleEle;		 float m_dRZZEffEleMuo;		      float m_dRZZEffMuoMuo;
  float m_dRZZSvFitEleEle;		 float m_dRZZSvFitEleMuo;	      float m_dRZZSvFitMuoMuo;
  float m_dRZZCAEleEle;                  float m_dRZZCAEleMuo;                float m_dRZZCAMuoMuo;
  float m_LepPt1EleEle;	        	 float m_LepPt1EleMuo;	              float m_LepPt1MuoMuo;
  float m_LepEta1EleEle;		 float m_LepEta1EleMuo;	              float m_LepEta1MuoMuo;
  float m_LepPFIso1EleEle;		 float m_LepPFIso1EleMuo;	      float m_LepPFIso1MuoMuo;
  float m_LepDetIso1EleEle;		 float m_LepDetIso1EleMuo;	      float m_LepDetIso1MuoMuo;
  float m_LepPt2EleEle;	          	 float m_LepPt2EleMuo;	              float m_LepPt2MuoMuo;
  float m_LepEta2EleEle;		 float m_LepEta2EleMuo;	              float m_LepEta2MuoMuo;
  float m_LepPFIso2EleEle;      	 float m_LepPFIso2EleMuo;             float m_LepPFIso2MuoMuo;
  float m_LepDetIso2EleEle;     	 float m_LepDetIso2EleMuo;            float m_LepDetIso2MuoMuo;
  float m_MassVisEleEle;		 float m_MassVisEleMuo;	      	      float m_MassVisMuoMuo;
  float m_MassEffEleEle;		 float m_MassEffEleMuo;	      	      float m_MassEffMuoMuo;
  float m_MassSvfitEleEle;		 float m_MassSvfitEleMuo;	      float m_MassSvfitMuoMuo;
  float m_MassCAEleEle;			 float m_MassCAEleMuo;		      float m_MassCAMuoMuo;
  float m_XMassVisEleEle;		 float m_XMassVisEleMuo;	      float m_XMassVisMuoMuo;
  float m_XMassEffEleEle;		 float m_XMassEffEleMuo;	      float m_XMassEffMuoMuo;
  float m_XMassSVFitEleEle;		 float m_XMassSVFitEleMuo;	      float m_XMassSVFitMuoMuo;
  float m_XMassCAEleEle;              	 float m_XMassCAEleMuo;               float m_XMassCAMuoMuo;
  float m_metEleEle;                     float m_metEleMuo;		      float m_metMuoMuo;
  int m_nbtagsLEleEle;			 int m_nbtagsLEleMuo;		      int m_nbtagsLMuoMuo;
  int m_nbtagsMEleEle;			 int m_nbtagsMEleMuo;		      int m_nbtagsMMuoMuo;
  int m_nbtagsTEleEle;                   int m_nbtagsTEleMuo;                 int m_nbtagsTMuoMuo;
  int m_trigger650EleEle;                int m_trigger650EleMuo;              int m_trigger650MuoMuo;
  int m_trigger320EleEle;                int m_trigger320EleMuo;              int m_trigger320MuoMuo;
  int m_NVerticesEleEle;                 int m_NVerticesEleMuo;		      int m_NVerticesMuoMuo;
  float m_PUWeightEleEle;		 float m_PUWeightEleMuo;	      float m_PUWeightMuoMuo;
  float m_metPxEleEle;			 float m_metPxEleMuo;		      float m_metPxMuoMuo;
  float m_metPyEleEle;			 float m_metPyEleMuo;		      float m_metPyMuoMuo;
  float m_metEtaEleEle;			 float m_metEtaEleMuo;		      float m_metEtaMuoMuo;
  float m_metPhiEleEle;                  float m_metPhiEleMuo;                float m_metPhiMuoMuo;
  
  edm::LumiReWeighting LumiWeights_;
  bool isData; 
  bool sideband; 
  edm::InputTag vtxColl_;
  edm::InputTag jetColl_;
  edm::InputTag jetPrunedColl_;
  edm::InputTag metColl_;
  edm::InputTag electronColl_;
  edm::InputTag muonColl_;
  edm::InputTag metRawColl_;
  edm::InputTag ak5JetColl_;

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
FullyLeptonicAnalyzer::FullyLeptonicAnalyzer(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
  isData = iConfig.getUntrackedParameter<bool>("isData_");
  sideband = iConfig.getUntrackedParameter<bool>("sideband_");
  vtxColl_ = iConfig.getParameter<edm::InputTag>("vtxColl"); 
  jetColl_ = iConfig.getParameter<edm::InputTag>("jetColl"); 
  jetPrunedColl_ = iConfig.getParameter<edm::InputTag>("jetPrunedColl"); 
  metColl_ = iConfig.getParameter<edm::InputTag>("metColl"); 
  electronColl_ = iConfig.getParameter<edm::InputTag>("electronColl"); 
  muonColl_ = iConfig.getParameter<edm::InputTag>("muonColl"); 
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


FullyLeptonicAnalyzer::~FullyLeptonicAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event ------------
void
FullyLeptonicAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  edm::Handle<pat::ElectronCollection> eleH;
  iEvent.getByLabel(electronColl_, eleH);

  edm::Handle<pat::MuonCollection> muoH;
  iEvent.getByLabel(muonColl_, muoH);

  edm::Handle<pat::METCollection> met;
  iEvent.getByLabel(metColl_, met);

  edm::Handle<pat::METCollection> metRaw;
  iEvent.getByLabel(metRawColl_, metRaw);

  edm::Handle<pat::JetCollection> ak5jetCands;
  iEvent.getByLabel(ak5JetColl_,ak5jetCands);

  edm::Handle<double> rhoHandle;
  iEvent.getByLabel("kt6PFJets", "rho", rhoHandle);
  float rho = *(rhoHandle.product());



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

  //SELECT THE FAT JET WITH THE HIGHTEST pt 
  int Njet = 0;
  vector<pat::JetCollection::const_iterator> SelectedJet;
  vector<float> SelectedPrunedMass;
  SelectJet(CA8JetswithQjets, CA8JetsPruned, SelectedJet, SelectedPrunedMass, Njet);
  float jetpt=-99; int jetInd = -1;
  for(unsigned int j=0; j<SelectedJet.size(); j++){
    if(SelectedJet[j]->pt()>jetpt) {
      jetInd = j;
      jetpt = SelectedJet[j]->pt();
    }
  }

  //SELECT THE ELCTRONS AND APPLY THE CUTS BETWEEN THE JET AND THE ELECTRONS
  int Nele = 0;
  vector<pat::ElectronCollection::const_iterator> SelectedInitialEle;
  vector<pat::ElectronCollection::const_iterator> SelectedEle;
  SelectElectronCutBased(eleH, SelectedInitialEle, Nele, rho, primaryVertex);
  for(unsigned int i=0; i<SelectedInitialEle.size(); i++){
    if(jetInd==-1) continue;
    if(ElectronPFIso(SelectedInitialEle[i],rho)>0.1) continue;
    if(ROOT::Math::VectorUtil::DeltaR(SelectedInitialEle[i]->p4(),SelectedJet[jetInd]->p4())<0.8) continue;
    if(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedInitialEle[i]->p4())>1) continue;
    SelectedEle.push_back(SelectedInitialEle[i]);
  }
   
  //SELECT THE MUONS AND APPLY THE CUTS BETWEEN THE JET AND THE MUONS
  int Nmuo = 0;
  vector<pat::MuonCollection::const_iterator> SelectedTrackerMuo;
  vector<pat::MuonCollection::const_iterator> SelectedHighptMuo;
  vector<pat::MuonCollection::const_iterator> SelectedMuo;
  SelectHighptMuon( muoH, SelectedHighptMuo, Nmuo, primaryVertex);
  SelectTrackerMuon(muoH, SelectedTrackerMuo, Nmuo, primaryVertex);
  for(unsigned int i=0; i<SelectedTrackerMuo.size(); i++){
    if(jetInd==-1) continue;
    if(ROOT::Math::VectorUtil::DeltaR(SelectedTrackerMuo[i]->p4(),SelectedJet[jetInd]->p4())<0.8) continue;
    if(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedTrackerMuo[i]->p4())>1) continue;
    if((MuonDETIso(SelectedTrackerMuo[i],SelectedTrackerMuo))>0.2) continue;
    SelectedMuo.push_back(SelectedTrackerMuo[i]);
  }
  bool hasAtLeastOneHighPtMuo = false;
  SelectTrackerGlobalID(SelectedMuo, primaryVertex, hasAtLeastOneHighPtMuo);
   
  //BTAG VETO
  int nbtagsL=0; int nbtagsM=0; int nbtagsT=0;
  if(jetInd!=-1) BtagVeto(ak5jetCands, nbtagsL, nbtagsM, nbtagsT, SelectedJet[jetInd]);

  //DI-ELECTRON SELECTION
  if(SelectedEle.size()==2 && SelectedJet.size()>0){
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedJet[jetInd]->pt(),SelectedJet[jetInd]->eta(),SelectedJet[jetInd]->phi(),SelectedPrunedMass[jetInd]);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());

    //SVFIT
    TLorentzVector SVFitTauTau;
    float MassSVFit = SVFitMass(met, metRaw, SVFitTauTau, SelectedEle[0]->p4(), SelectedEle[1]->p4());
    float XMassSVFit = (SVFitTauTau+PrunedJet).M();
    float dRJetZSVFit = ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,PrunedJet);
 
    //COLLINEAR APPROXIMATION
    TLorentzVector CATauTau; bool CA = false;
    CollinearApproximation(met, CATauTau, PrunedJet, SelectedEle[0]->p4(), SelectedEle[1]->p4(), CA);
    float MassCA = 0;
    float XmassCA = 0.;
    float dRJetZCA = 0.;
    if(CA){
      MassCA = CATauTau.M();
      XmassCA = (CATauTau+PrunedJet).M();
      dRJetZCA = ROOT::Math::VectorUtil::DeltaR(CATauTau,PrunedJet);
    }

    pat::ElectronCollection::const_iterator SelectedLep1;
    pat::ElectronCollection::const_iterator SelectedLep2;
    if(SelectedEle[0]->pt()>SelectedEle[1]->pt()){
      SelectedLep1=SelectedEle[0];
      SelectedLep2=SelectedEle[1];
    } else {
      SelectedLep1=SelectedEle[1];
      SelectedLep2=SelectedEle[0];
    }
    math::PtEtaPhiELorentzVector dilep; dilep = SelectedLep1->p4()+SelectedLep2->p4();
    math::PtEtaPhiELorentzVector dilepmet; dilepmet = SelectedLep1->p4()+SelectedLep2->p4()+met->begin()->p4();
    math::PtEtaPhiELorentzVector dilepjet; dilepjet = SelectedLep1->p4()+SelectedLep2->p4()+PrunedJet_prov;
    math::PtEtaPhiELorentzVector dilepmetjet; dilepmetjet = SelectedLep1->p4()+SelectedLep2->p4()+met->begin()->p4()+PrunedJet_prov;

    //PLOT - TREE
    m_jetPtEleEle=SelectedJet[jetInd]->pt();
    m_jetEtaEleEle=SelectedJet[jetInd]->eta();
    m_jetMassEleEle=SelectedPrunedMass[jetInd];
    m_jetSubjettinessEleEle=SelectedJet[jetInd]->userFloat("tau2")/SelectedJet[jetInd]->userFloat("tau1");
    m_dRJetLep1EleEle=ROOT::Math::VectorUtil::DeltaR(SelectedLep1->p4(),SelectedJet[jetInd]->p4());
    m_dRJetLep2EleEle=ROOT::Math::VectorUtil::DeltaR(SelectedLep2->p4(),SelectedJet[jetInd]->p4());
    m_dRJetMetEleEle=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedJet[jetInd]->p4());
    m_dPhiJetMetEleEle=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet[jetInd]->p4());
    m_dRLepMet1EleEle=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedLep1->p4());
    m_dPhiLepMet1EleEle=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedLep1->p4());
    m_dRLepMet2EleEle=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedLep2->p4());
    m_dPhiLepMet2EleEle=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedLep2->p4());
    m_dRLepLepEleEle=ROOT::Math::VectorUtil::DeltaR(SelectedLep1->p4(),SelectedLep2->p4());
    m_dRZZVisEleEle=ROOT::Math::VectorUtil::DeltaR(dilep,PrunedJet);
    m_dRZZEffEleEle=ROOT::Math::VectorUtil::DeltaR(dilepmet,PrunedJet);
    m_dRZZSvFitEleEle=dRJetZSVFit;
    m_dRZZCAEleEle=dRJetZCA;
    m_LepPt1EleEle=SelectedLep1->pt();
    m_LepEta1EleEle=SelectedLep1->eta();
    m_LepPFIso1EleEle=ElectronPFIso(SelectedLep1,rho);
    m_LepDetIso1EleEle=(SelectedLep1->userIso(0)+SelectedLep1->userIso(1)+SelectedLep1->userIso(2))/SelectedLep1->pt();
    m_LepPt2EleEle=SelectedLep2->pt();
    m_LepEta2EleEle=SelectedLep2->eta();
    m_LepPFIso2EleEle=ElectronPFIso(SelectedLep2,rho);
    m_LepDetIso2EleEle=(SelectedLep2->userIso(0)+SelectedLep2->userIso(1)+SelectedLep2->userIso(2))/SelectedLep2->pt();
    m_MassVisEleEle=dilep.mass();
    m_MassEffEleEle=dilepmet.mass();
    m_MassSvfitEleEle=MassSVFit;
    m_MassCAEleEle=MassCA;
    m_XMassVisEleEle=dilepjet.mass();
    m_XMassEffEleEle=dilepmetjet.mass();
    m_XMassSVFitEleEle=XMassSVFit;
    m_XMassCAEleEle=XmassCA;
    m_metEleEle=met->begin()->pt();
    m_nbtagsLEleEle=nbtagsL;
    m_nbtagsMEleEle=nbtagsM;
    m_nbtagsTEleEle=nbtagsT;
    m_trigger320EleEle=(int)isFired_HLT_PFJet320;
    m_trigger650EleEle=(int)isFired_HLT_HT650; 
    m_NVerticesEleEle=vertices->size();
    m_PUWeightEleEle=MyWeight;
    m_metPxEleEle=met->begin()->px();
    m_metPyEleEle=met->begin()->py();
    m_metEtaEleEle=met->begin()->eta();
    m_metPhiEleEle=met->begin()->phi();
    TreeEleEle->Fill();
  }



  //DI-MUON SELECTION
  if(SelectedMuo.size()==2 && hasAtLeastOneHighPtMuo == true && SelectedJet.size()>0){
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedJet[jetInd]->pt(),SelectedJet[jetInd]->eta(),SelectedJet[jetInd]->phi(),SelectedPrunedMass[jetInd]);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());

    //SVFIT
    TLorentzVector SVFitTauTau;
    float MassSVFit = SVFitMass(met, metRaw, SVFitTauTau, SelectedMuo[0]->p4(), SelectedMuo[1]->p4());
    float XMassSVFit = (SVFitTauTau+PrunedJet).M();
    float dRJetZSVFit = ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,PrunedJet);

    //COLLINEAR APPROXIMATION
    TLorentzVector CATauTau; bool CA = false;
    CollinearApproximation(met, CATauTau, PrunedJet, SelectedMuo[0]->p4(), SelectedMuo[1]->p4(), CA);
    float MassCA = 0;
    float XmassCA = 0.;
    float dRJetZCA = 0.;
    if(CA){
      MassCA = CATauTau.M();
      XmassCA = (CATauTau+PrunedJet).M();
      dRJetZCA = ROOT::Math::VectorUtil::DeltaR(CATauTau,PrunedJet);
    }

    pat::MuonCollection::const_iterator SelectedLep1;
    pat::MuonCollection::const_iterator SelectedLep2;
    if(SelectedMuo[0]->pt()>SelectedMuo[1]->pt()){
      SelectedLep1=SelectedMuo[0];
      SelectedLep2=SelectedMuo[1];
    } else {
      SelectedLep1=SelectedMuo[1];
      SelectedLep2=SelectedMuo[0];
    }
    math::PtEtaPhiELorentzVector dilep; dilep = SelectedLep1->p4()+SelectedLep2->p4();
    math::PtEtaPhiELorentzVector dilepmet; dilepmet = SelectedLep1->p4()+SelectedLep2->p4()+met->begin()->p4();
    math::PtEtaPhiELorentzVector dilepjet; dilepjet = SelectedLep1->p4()+SelectedLep2->p4()+PrunedJet_prov;
    math::PtEtaPhiELorentzVector dilepmetjet; dilepmetjet = SelectedLep1->p4()+SelectedLep2->p4()+met->begin()->p4()+PrunedJet_prov;

    //PLOT - TREE
    m_jetPtMuoMuo=SelectedJet[jetInd]->pt();
    m_jetEtaMuoMuo=SelectedJet[jetInd]->eta();
    m_jetMassMuoMuo=SelectedPrunedMass[jetInd];
    m_jetSubjettinessMuoMuo=SelectedJet[jetInd]->userFloat("tau2")/SelectedJet[jetInd]->userFloat("tau1");
    m_dRJetLep1MuoMuo=ROOT::Math::VectorUtil::DeltaR(SelectedLep1->p4(),SelectedJet[jetInd]->p4());
    m_dRJetLep2MuoMuo=ROOT::Math::VectorUtil::DeltaR(SelectedLep2->p4(),SelectedJet[jetInd]->p4());
    m_dRJetMetMuoMuo=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedJet[jetInd]->p4());
    m_dPhiJetMetMuoMuo=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet[jetInd]->p4());
    m_dRLepMet1MuoMuo=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedLep1->p4());
    m_dPhiLepMet1MuoMuo=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedLep1->p4());
    m_dRLepMet2MuoMuo=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedLep2->p4());
    m_dPhiLepMet2MuoMuo=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedLep2->p4());
    m_dRLepLepMuoMuo=ROOT::Math::VectorUtil::DeltaR(SelectedLep1->p4(),SelectedLep2->p4());
    m_dRZZVisMuoMuo=ROOT::Math::VectorUtil::DeltaR(dilep,PrunedJet);
    m_dRZZEffMuoMuo=ROOT::Math::VectorUtil::DeltaR(dilepmet,PrunedJet);
    m_dRZZSvFitMuoMuo=dRJetZSVFit;
    m_dRZZCAMuoMuo=dRJetZCA;
    m_LepPt1MuoMuo=SelectedLep1->pt();
    m_LepEta1MuoMuo=SelectedLep1->eta();
    m_LepPFIso1MuoMuo=MuonPFIso(SelectedLep1,false);
    m_LepDetIso1MuoMuo=(MuonDETIso(SelectedLep1, SelectedMuo));
    m_LepPt2MuoMuo=SelectedLep2->pt();
    m_LepEta2MuoMuo=SelectedLep2->eta();
    m_LepPFIso2MuoMuo=MuonPFIso(SelectedLep2,false);
    m_LepDetIso2MuoMuo=(MuonDETIso(SelectedLep2, SelectedMuo));
    m_MassVisMuoMuo=dilep.mass();
    m_MassEffMuoMuo=dilepmet.mass();
    m_MassSvfitMuoMuo=MassSVFit;
    m_MassCAMuoMuo=MassCA;
    m_XMassVisMuoMuo=dilepjet.mass();
    m_XMassEffMuoMuo=dilepmetjet.mass();
    m_XMassSVFitMuoMuo=XMassSVFit;
    m_XMassCAMuoMuo=XmassCA;
    m_metMuoMuo=met->begin()->pt();
    m_nbtagsLMuoMuo=nbtagsL;
    m_nbtagsMMuoMuo=nbtagsM;
    m_nbtagsTMuoMuo=nbtagsT;
    m_trigger320MuoMuo=(int)isFired_HLT_PFJet320;
    m_trigger650MuoMuo=(int)isFired_HLT_HT650;
    m_NVerticesMuoMuo=vertices->size();
    m_PUWeightMuoMuo=MyWeight;
    m_metPxMuoMuo=met->begin()->px();
    m_metPyMuoMuo=met->begin()->py();
    m_metEtaMuoMuo=met->begin()->eta();
    m_metPhiMuoMuo=met->begin()->phi();
    TreeMuoMuo->Fill();
  }



  //ELECTRON-MUON SELECTION
  SelectedEle.clear(); SelectedMuo.clear();
  for(unsigned int i=0; i<SelectedInitialEle.size(); i++){
    if(jetInd==-1) continue;
    if(ROOT::Math::VectorUtil::DeltaR(SelectedInitialEle[i]->p4(),SelectedJet[jetInd]->p4())<0.8) continue;
    if(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedInitialEle[i]->p4())>1) continue;
    if(ElectronPFIso(SelectedInitialEle[i],rho)>0.1) continue;
    SelectedEle.push_back(SelectedInitialEle[i]);
  }
  for(unsigned int i=0; i<SelectedHighptMuo.size(); i++){
    if(jetInd==-1) continue;
    if(ROOT::Math::VectorUtil::DeltaR(SelectedHighptMuo[i]->p4(),SelectedJet[jetInd]->p4())<0.8) continue;
    if(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedHighptMuo[i]->p4())>1) continue;
    if(MuonPFIso(SelectedHighptMuo[i], true)>0.2) continue;
    SelectedMuo.push_back(SelectedHighptMuo[i]);
  }

  if(SelectedEle.size()==1 && SelectedMuo.size()==1 && SelectedJet.size()>0){
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedJet[jetInd]->pt(),SelectedJet[jetInd]->eta(),SelectedJet[jetInd]->phi(),SelectedPrunedMass[jetInd]);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());

    //SVFIT
    TLorentzVector SVFitTauTau;
    float MassSVFit = SVFitMass(met, metRaw, SVFitTauTau, SelectedEle[0]->p4(), SelectedMuo[0]->p4());
    float XMassSVFit = (SVFitTauTau+PrunedJet).M();
    float dRJetZSVFit = ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,PrunedJet);

    //COLLINEAR APPROXIMATION
    TLorentzVector CATauTau; bool CA = false;
    CollinearApproximation(met, CATauTau, PrunedJet, SelectedEle[0]->p4(), SelectedMuo[0]->p4(), CA);
    float MassCA = 0;
    float XmassCA = 0.;
    float dRJetZCA = 0.;
    if(CA){
      MassCA = CATauTau.M();
      XmassCA = (CATauTau+PrunedJet).M();
      dRJetZCA = ROOT::Math::VectorUtil::DeltaR(CATauTau,PrunedJet);
    }

    pat::ElectronCollection::const_iterator SelectedLep1;
    pat::MuonCollection::const_iterator SelectedLep2;
    SelectedLep1=SelectedEle[0];
    SelectedLep2=SelectedMuo[0];
    math::PtEtaPhiELorentzVector dilep; dilep = SelectedLep1->p4()+SelectedLep2->p4();
    math::PtEtaPhiELorentzVector dilepmet; dilepmet = SelectedLep1->p4()+SelectedLep2->p4()+met->begin()->p4();
    math::PtEtaPhiELorentzVector dilepjet; dilepjet = SelectedLep1->p4()+SelectedLep2->p4()+PrunedJet_prov;
    math::PtEtaPhiELorentzVector dilepmetjet; dilepmetjet = SelectedLep1->p4()+SelectedLep2->p4()+met->begin()->p4()+PrunedJet_prov;

    //PLOT - TREE
    m_jetPtEleMuo=SelectedJet[jetInd]->pt();
    m_jetEtaEleMuo=SelectedJet[jetInd]->eta();
    m_jetMassEleMuo=SelectedPrunedMass[jetInd];
    m_jetSubjettinessEleMuo=SelectedJet[jetInd]->userFloat("tau2")/SelectedJet[jetInd]->userFloat("tau1");
    m_dRJetLep1EleMuo=ROOT::Math::VectorUtil::DeltaR(SelectedLep1->p4(),SelectedJet[jetInd]->p4());
    m_dRJetLep2EleMuo=ROOT::Math::VectorUtil::DeltaR(SelectedLep2->p4(),SelectedJet[jetInd]->p4());
    m_dRJetMetEleMuo=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedJet[jetInd]->p4());
    m_dPhiJetMetEleMuo=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet[jetInd]->p4());
    m_dRLepMet1EleMuo=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedLep1->p4());
    m_dPhiLepMet1EleMuo=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedLep1->p4());
    m_dRLepMet2EleMuo=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedLep2->p4());
    m_dPhiLepMet2EleMuo=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedLep2->p4());
    m_dRLepLepEleMuo=ROOT::Math::VectorUtil::DeltaR(SelectedLep1->p4(),SelectedLep2->p4());
    m_dRZZVisEleMuo=ROOT::Math::VectorUtil::DeltaR(dilep,PrunedJet);
    m_dRZZEffEleMuo=ROOT::Math::VectorUtil::DeltaR(dilepmet,PrunedJet);
    m_dRZZSvFitEleMuo=dRJetZSVFit;
    m_dRZZCAEleMuo=dRJetZCA;
    m_LepPt1EleMuo=SelectedLep1->pt();
    m_LepEta1EleMuo=SelectedLep1->eta();
    m_LepPFIso1EleMuo=ElectronPFIso(SelectedLep1,rho);
    m_LepDetIso1EleMuo=(SelectedLep1->userIso(0)+SelectedLep1->userIso(1)+SelectedLep1->userIso(2))/SelectedLep1->pt();
    m_LepPt2EleMuo=SelectedLep2->pt();
    m_LepEta2EleMuo=SelectedLep2->eta();
    m_LepPFIso2EleMuo=MuonPFIso(SelectedLep2,false);
    m_LepDetIso2EleMuo=(MuonDETIso(SelectedLep2, SelectedMuo));
    m_MassVisEleMuo=dilep.mass();
    m_MassEffEleMuo=dilepmet.mass();
    m_MassSvfitEleMuo=MassSVFit;
    m_MassCAEleMuo=MassCA;
    m_XMassVisEleMuo=dilepjet.mass();
    m_XMassEffEleMuo=dilepmetjet.mass();
    m_XMassSVFitEleMuo=XMassSVFit;
    m_XMassCAEleMuo=XmassCA;
    m_metEleMuo=met->begin()->pt();
    m_nbtagsLEleMuo=nbtagsL;
    m_nbtagsMEleMuo=nbtagsM;
    m_nbtagsTEleMuo=nbtagsT;
    m_trigger320EleMuo=(int)isFired_HLT_PFJet320;
    m_trigger650EleMuo=(int)isFired_HLT_HT650;
    m_NVerticesEleMuo=vertices->size();
    m_PUWeightEleMuo=MyWeight;
    m_metPxEleMuo=met->begin()->px();
    m_metPyEleMuo=met->begin()->py();
    m_metEtaEleMuo=met->begin()->eta();
    m_metPhiEleMuo=met->begin()->phi();
    TreeEleMuo->Fill();
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
FullyLeptonicAnalyzer::beginJob()
{
  Service<TFileService> fs;
  Nevents = fs->make<TH1D>("Nevents", "Nevents", 3, -0.5, 2.5);

  TreeEleEle = fs->make<TTree>("TreeEleEle", "TreeEleEle");
  TreeEleEle->Branch("jetPtEleEle", &m_jetPtEleEle, "jetPtEleEle/f");
  TreeEleEle->Branch("jetEtaEleEle", &m_jetEtaEleEle, "jetEtaEleEle/f");
  TreeEleEle->Branch("jetMassEleEle", &m_jetMassEleEle, "jetMassEleEle/f");
  TreeEleEle->Branch("jetSubjettinessEleEle", &m_jetSubjettinessEleEle, "jetSubjettinessEleEle/f");
  TreeEleEle->Branch("dRJetLep1EleEle", &m_dRJetLep1EleEle, "dRJetLep1EleEle/f");
  TreeEleEle->Branch("dRJetLep2EleEle", &m_dRJetLep2EleEle, "dRJetLep2EleEle/f");
  TreeEleEle->Branch("dRJetMetEleEle", &m_dRJetMetEleEle, "dRJetMetEleEle/f");
  TreeEleEle->Branch("dPhiJetMetEleEle", &m_dPhiJetMetEleEle, "dPhiJetMetEleEle/f");
  TreeEleEle->Branch("dRLepMet1EleEle", &m_dRLepMet1EleEle, "dRLepMet1EleEle/f");
  TreeEleEle->Branch("dPhiLepMet1EleEle", &m_dPhiLepMet1EleEle, "dPhiLepMet1EleEle/f");
  TreeEleEle->Branch("dRLepMet2EleEle", &m_dRLepMet2EleEle, "dRLepMet2EleEle/f");
  TreeEleEle->Branch("dPhiLepMet2EleEle", &m_dPhiLepMet2EleEle, "dPhiLepMet2EleEle/f");
  TreeEleEle->Branch("dRLepLepEleEle", &m_dRLepLepEleEle, "dRLepLepEleEle/f");
  TreeEleEle->Branch("LepPt1EleEle", &m_LepPt1EleEle, "LepPt1EleEle/f");
  TreeEleEle->Branch("LepEta1EleEle", &m_LepEta1EleEle, "LepEta1EleEle/f");
  TreeEleEle->Branch("LepPFIso1EleEle", &m_LepPFIso1EleEle, "LepPFIso1EleEle/f");
  TreeEleEle->Branch("LepDetIso1EleEle", &m_LepDetIso1EleEle, "LepDetIso1EleEle/f");
  TreeEleEle->Branch("LepPt2EleEle", &m_LepPt2EleEle, "LepPt2EleEle/f");
  TreeEleEle->Branch("LepEta2EleEle", &m_LepEta2EleEle, "LepEta2EleEle/f");
  TreeEleEle->Branch("LepPFIso2EleEle", &m_LepPFIso2EleEle, "LepPFIso2EleEle/f");
  TreeEleEle->Branch("LepDetIso2EleEle", &m_LepDetIso2EleEle, "LepDetIso2EleEle/f");
  TreeEleEle->Branch("MassVisEleEle", &m_MassVisEleEle, "MassVisEleEle/f");
  TreeEleEle->Branch("MassEffEleEle", &m_MassEffEleEle, "MassEffEleEle/f");
  TreeEleEle->Branch("MassSvfitEleEle", &m_MassSvfitEleEle, "MassSvfitEleEle/f");
  TreeEleEle->Branch("MassCAEleEle", &m_MassCAEleEle, "MassCAEleEle/f");
  TreeEleEle->Branch("XMassVisEleEle", &m_XMassVisEleEle, "XMassVisEleEle/f");
  TreeEleEle->Branch("XMassEffEleEle", &m_XMassEffEleEle, "XMassEffEleEle/f");
  TreeEleEle->Branch("XMassSVFitEleEle", &m_XMassSVFitEleEle, "XMassSVFitEleEle/f");
  TreeEleEle->Branch("XMassCAEleEle", &m_XMassCAEleEle, "XMassCAEleEle/f");
  TreeEleEle->Branch("dRZZVisEleEle", &m_dRZZVisEleEle, "dRZZVisEleEle/f");
  TreeEleEle->Branch("dRZZEffEleEle", &m_dRZZEffEleEle, "dRZZEffEleEle/f");
  TreeEleEle->Branch("dRZZSvFitEleEle", &m_dRZZSvFitEleEle, "dRZZSvFitEleEle/f");
  TreeEleEle->Branch("dRZZCAEleEle", &m_dRZZCAEleEle, "dRZZCAEleEle/f");
  TreeEleEle->Branch("metEleEle", &m_metEleEle, "metEleEle/f");
  TreeEleEle->Branch("nbtagsLEleEle", &m_nbtagsLEleEle, "nbtagsLEleEle/i");
  TreeEleEle->Branch("nbtagsMEleEle", &m_nbtagsMEleEle, "nbtagsMEleEle/i");
  TreeEleEle->Branch("nbtagsTEleEle", &m_nbtagsTEleEle, "nbtagsTEleEle/i");
  TreeEleEle->Branch("trigger320EleEle", &m_trigger320EleEle, "trigger320EleEle/i");
  TreeEleEle->Branch("trigger650EleEle", &m_trigger650EleEle, "trigger650EleEle/i");
  TreeEleEle->Branch("NVerticesEleEle", &m_NVerticesEleEle, "NVerticesEleEle/i");
  TreeEleEle->Branch("PUWeightEleEle", &m_PUWeightEleEle, "PUWeightEleEle/f");
  TreeEleEle->Branch("metPxEleEle", &m_metPxEleEle, "metPxEleEle/f");
  TreeEleEle->Branch("metPyEleEle", &m_metPyEleEle, "metPyEleEle/f");
  TreeEleEle->Branch("metEtaEleEle", &m_metEtaEleEle, "metEtaEleEle/f");
  TreeEleEle->Branch("metPhiEleEle", &m_metPhiEleEle, "metPhiEleEle/f");

  TreeMuoMuo = fs->make<TTree>("TreeMuoMuo", "TreeMuoMuo");
  TreeMuoMuo->Branch("jetPtMuoMuo", &m_jetPtMuoMuo, "jetPtMuoMuo/f");
  TreeMuoMuo->Branch("jetEtaMuoMuo", &m_jetEtaMuoMuo, "jetEtaMuoMuo/f");
  TreeMuoMuo->Branch("jetMassMuoMuo", &m_jetMassMuoMuo, "jetMassMuoMuo/f");
  TreeMuoMuo->Branch("jetSubjettinessMuoMuo", &m_jetSubjettinessMuoMuo, "jetSubjettinessMuoMuo/f");
  TreeMuoMuo->Branch("dRJetLep1MuoMuo", &m_dRJetLep1MuoMuo, "dRJetLep1MuoMuo/f");
  TreeMuoMuo->Branch("dRJetLep2MuoMuo", &m_dRJetLep2MuoMuo, "dRJetLep2MuoMuo/f");
  TreeMuoMuo->Branch("dRJetMetMuoMuo", &m_dRJetMetMuoMuo, "dRJetMetMuoMuo/f");
  TreeMuoMuo->Branch("dPhiJetMetMuoMuo", &m_dPhiJetMetMuoMuo, "dPhiJetMetMuoMuo/f");
  TreeMuoMuo->Branch("dRLepMet1MuoMuo", &m_dRLepMet1MuoMuo, "dRLepMet1MuoMuo/f");
  TreeMuoMuo->Branch("dPhiLepMet1MuoMuo", &m_dPhiLepMet1MuoMuo, "dPhiLepMet1MuoMuo/f");
  TreeMuoMuo->Branch("dRLepMet2MuoMuo", &m_dRLepMet2MuoMuo, "dRLepMet2MuoMuo/f");
  TreeMuoMuo->Branch("dPhiLepMet2MuoMuo", &m_dPhiLepMet2MuoMuo, "dPhiLepMet2MuoMuo/f");
  TreeMuoMuo->Branch("dRLepLepMuoMuo", &m_dRLepLepMuoMuo, "dRLepLepMuoMuo/f");
  TreeMuoMuo->Branch("LepPt1MuoMuo", &m_LepPt1MuoMuo, "LepPt1MuoMuo/f");
  TreeMuoMuo->Branch("LepEta1MuoMuo", &m_LepEta1MuoMuo, "LepEta1MuoMuo/f");
  TreeMuoMuo->Branch("LepPFIso1MuoMuo", &m_LepPFIso1MuoMuo, "LepPFIso1MuoMuo/f");
  TreeMuoMuo->Branch("LepDetIso1MuoMuo", &m_LepDetIso1MuoMuo, "LepDetIso1MuoMuo/f");
  TreeMuoMuo->Branch("LepPt2MuoMuo", &m_LepPt2MuoMuo, "LepPt2MuoMuo/f");
  TreeMuoMuo->Branch("LepEta2MuoMuo", &m_LepEta2MuoMuo, "LepEta2MuoMuo/f");
  TreeMuoMuo->Branch("LepPFIso2MuoMuo", &m_LepPFIso2MuoMuo, "LepPFIso2MuoMuo/f");
  TreeMuoMuo->Branch("LepDetIso2MuoMuo", &m_LepDetIso2MuoMuo, "LepDetIso2MuoMuo/f");
  TreeMuoMuo->Branch("MassVisMuoMuo", &m_MassVisMuoMuo, "MassVisMuoMuo/f");
  TreeMuoMuo->Branch("MassEffMuoMuo", &m_MassEffMuoMuo, "MassEffMuoMuo/f");
  TreeMuoMuo->Branch("MassSvfitMuoMuo", &m_MassSvfitMuoMuo, "MassSvfitMuoMuo/f");
  TreeMuoMuo->Branch("MassCAMuoMuo", &m_MassCAMuoMuo, "MassCAMuoMuo/f");
  TreeMuoMuo->Branch("XMassVisMuoMuo", &m_XMassVisMuoMuo, "XMassVisMuoMuo/f");
  TreeMuoMuo->Branch("XMassEffMuoMuo", &m_XMassEffMuoMuo, "XMassEffMuoMuo/f");
  TreeMuoMuo->Branch("XMassSVFitMuoMuo", &m_XMassSVFitMuoMuo, "XMassSVFitMuoMuo/f");
  TreeMuoMuo->Branch("XMassCAMuoMuo", &m_XMassCAMuoMuo, "XMassCAMuoMuo/f");
  TreeMuoMuo->Branch("dRZZVisMuoMuo", &m_dRZZVisMuoMuo, "dRZZVisMuoMuo/f");
  TreeMuoMuo->Branch("dRZZEffMuoMuo", &m_dRZZEffMuoMuo, "dRZZEffMuoMuo/f");
  TreeMuoMuo->Branch("dRZZSvFitMuoMuo", &m_dRZZSvFitMuoMuo, "dRZZSvFitMuoMuo/f");
  TreeMuoMuo->Branch("dRZZCAMuoMuo", &m_dRZZCAMuoMuo, "dRZZCAMuoMuo/f");
  TreeMuoMuo->Branch("metMuoMuo", &m_metMuoMuo, "metMuoMuo/f");
  TreeMuoMuo->Branch("nbtagsLMuoMuo", &m_nbtagsLMuoMuo, "nbtagsLMuoMuo/i");
  TreeMuoMuo->Branch("nbtagsMMuoMuo", &m_nbtagsMMuoMuo, "nbtagsMMuoMuo/i");
  TreeMuoMuo->Branch("nbtagsTMuoMuo", &m_nbtagsTMuoMuo, "nbtagsTMuoMuo/i");
  TreeMuoMuo->Branch("trigger320MuoMuo", &m_trigger320MuoMuo, "trigger320MuoMuo/i");
  TreeMuoMuo->Branch("trigger650MuoMuo", &m_trigger650MuoMuo, "trigger650MuoMuo/i");
  TreeMuoMuo->Branch("NVerticesMuoMuo", &m_NVerticesMuoMuo, "NVerticesMuoMuo/i");
  TreeMuoMuo->Branch("PUWeightMuoMuo", &m_PUWeightMuoMuo, "PUWeightMuoMuo/f");
  TreeMuoMuo->Branch("metPxMuoMuo", &m_metPxMuoMuo, "metPxMuoMuo/f");
  TreeMuoMuo->Branch("metPyMuoMuo", &m_metPyMuoMuo, "metPyMuoMuo/f");
  TreeMuoMuo->Branch("metEtaMuoMuo", &m_metEtaMuoMuo, "metEtaMuoMuo/f");
  TreeMuoMuo->Branch("metPhiMuoMuo", &m_metPhiMuoMuo, "metPhiMuoMuo/f");

  TreeEleMuo = fs->make<TTree>("TreeEleMuo", "TreeEleMuo");
  TreeEleMuo->Branch("jetPtEleMuo", &m_jetPtEleMuo, "jetPtEleMuo/f");
  TreeEleMuo->Branch("jetEtaEleMuo", &m_jetEtaEleMuo, "jetEtaEleMuo/f");
  TreeEleMuo->Branch("jetMassEleMuo", &m_jetMassEleMuo, "jetMassEleMuo/f");
  TreeEleMuo->Branch("jetSubjettinessEleMuo", &m_jetSubjettinessEleMuo, "jetSubjettinessEleMuo/f");
  TreeEleMuo->Branch("dRJetLep1EleMuo", &m_dRJetLep1EleMuo, "dRJetLep1EleMuo/f");
  TreeEleMuo->Branch("dRJetLep2EleMuo", &m_dRJetLep2EleMuo, "dRJetLep2EleMuo/f");
  TreeEleMuo->Branch("dRJetMetEleMuo", &m_dRJetMetEleMuo, "dRJetMetEleMuo/f");
  TreeEleMuo->Branch("dPhiJetMetEleMuo", &m_dPhiJetMetEleMuo, "dPhiJetMetEleMuo/f");
  TreeEleMuo->Branch("dRLepMet1EleMuo", &m_dRLepMet1EleMuo, "dRLepMet1EleMuo/f");
  TreeEleMuo->Branch("dPhiLepMet1EleMuo", &m_dPhiLepMet1EleMuo, "dPhiLepMet1EleMuo/f");
  TreeEleMuo->Branch("dRLepMet2EleMuo", &m_dRLepMet2EleMuo, "dRLepMet2EleMuo/f");
  TreeEleMuo->Branch("dPhiLepMet2EleMuo", &m_dPhiLepMet2EleMuo, "dPhiLepMet2EleMuo/f");
  TreeEleMuo->Branch("dRLepLepEleMuo", &m_dRLepLepEleMuo, "dRLepLepEleMuo/f");
  TreeEleMuo->Branch("LepPt1EleMuo", &m_LepPt1EleMuo, "LepPt1EleMuo/f");
  TreeEleMuo->Branch("LepEta1EleMuo", &m_LepEta1EleMuo, "LepEta1EleMuo/f");
  TreeEleMuo->Branch("LepPFIso1EleMuo", &m_LepPFIso1EleMuo, "LepPFIso1EleMuo/f");
  TreeEleMuo->Branch("LepDetIso1EleMuo", &m_LepDetIso1EleMuo, "LepDetIso1EleMuo/f");
  TreeEleMuo->Branch("LepPt2EleMuo", &m_LepPt2EleMuo, "LepPt2EleMuo/f");
  TreeEleMuo->Branch("LepEta2EleMuo", &m_LepEta2EleMuo, "LepEta2EleMuo/f");
  TreeEleMuo->Branch("LepPFIso2EleMuo", &m_LepPFIso2EleMuo, "LepPFIso2EleMuo/f");
  TreeEleMuo->Branch("LepDetIso2EleMuo", &m_LepDetIso2EleMuo, "LepDetIso2EleMuo/f");
  TreeEleMuo->Branch("MassVisEleMuo", &m_MassVisEleMuo, "MassVisEleMuo/f");
  TreeEleMuo->Branch("MassEffEleMuo", &m_MassEffEleMuo, "MassEffEleMuo/f");
  TreeEleMuo->Branch("MassSvfitEleMuo", &m_MassSvfitEleMuo, "MassSvfitEleMuo/f");
  TreeEleMuo->Branch("MassCAEleMuo", &m_MassCAEleMuo, "MassCAEleMuo/f");
  TreeEleMuo->Branch("XMassVisEleMuo", &m_XMassVisEleMuo, "XMassVisEleMuo/f");
  TreeEleMuo->Branch("XMassEffEleMuo", &m_XMassEffEleMuo, "XMassEffEleMuo/f");
  TreeEleMuo->Branch("XMassSVFitEleMuo", &m_XMassSVFitEleMuo, "XMassSVFitEleMuo/f");
  TreeEleMuo->Branch("XMassCAEleMuo", &m_XMassCAEleMuo, "XMassCAEleMuo/f");
  TreeEleMuo->Branch("dRZZVisEleMuo", &m_dRZZVisEleMuo, "dRZZVisEleMuo/f");
  TreeEleMuo->Branch("dRZZEffEleMuo", &m_dRZZEffEleMuo, "dRZZEffEleMuo/f");
  TreeEleMuo->Branch("dRZZSvFitEleMuo", &m_dRZZSvFitEleMuo, "dRZZSvFitEleMuo/f");
  TreeEleMuo->Branch("dRZZCAEleMuo", &m_dRZZCAEleMuo, "dRZZCAEleMuo/f");
  TreeEleMuo->Branch("metEleMuo", &m_metEleMuo, "metEleMuo/f");
  TreeEleMuo->Branch("nbtagsLEleMuo", &m_nbtagsLEleMuo, "nbtagsLEleMuo/i");
  TreeEleMuo->Branch("nbtagsMEleMuo", &m_nbtagsMEleMuo, "nbtagsMEleMuo/i");
  TreeEleMuo->Branch("nbtagsTEleMuo", &m_nbtagsTEleMuo, "nbtagsTEleMuo/i");
  TreeEleMuo->Branch("trigger320EleMuo", &m_trigger320EleMuo, "trigger320EleMuo/i");
  TreeEleMuo->Branch("trigger650EleMuo", &m_trigger650EleMuo, "trigger650EleMuo/i");
  TreeEleMuo->Branch("NVerticesEleMuo", &m_NVerticesEleMuo, "NVerticesEleMuo/i");
  TreeEleMuo->Branch("PUWeightEleMuo", &m_PUWeightEleMuo, "PUWeightEleMuo/f");
  TreeEleMuo->Branch("metPxEleMuo", &m_metPxEleMuo, "metPxEleMuo/f");
  TreeEleMuo->Branch("metPyEleMuo", &m_metPyEleMuo, "metPyEleMuo/f");
  TreeEleMuo->Branch("metEtaEleMuo", &m_metEtaEleMuo, "metEtaEleMuo/f");
  TreeEleMuo->Branch("metPhiEleMuo", &m_metPhiEleMuo, "metPhiEleMuo/f");
}

// ------------ method called once each job just after ending the event loop ------------
void
FullyLeptonicAnalyzer::endJob()
{
}

// ------------ method called when starting to processes a run ------------
void
FullyLeptonicAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run ------------
void
FullyLeptonicAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block ------------
void
FullyLeptonicAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block ------------
void
FullyLeptonicAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module ------------
void
FullyLeptonicAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void FullyLeptonicAnalyzer::SelectJet(edm::Handle<pat::JetCollection> CA8JetswithQjets,
				      edm::Handle<pat::JetCollection> CA8JetsPruned,
				      vector<pat::JetCollection::const_iterator> & SelectedJet,
				      vector<float> & SelectedPrunedMass, int & Njet){

  SelectedJet.clear(); SelectedPrunedMass.clear();
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
    if(sideband) {if(!(mass>20 && mass<70))  continue;}
    else         {if(!(mass>70 && mass<110)) continue;}
    if(jet->userFloat("tau2")/jet->userFloat("tau1")>0.75) continue;
    Njet = Njet + 1;
    SelectedJet.push_back(jet);
    SelectedPrunedMass.push_back(mass);
  }
}


void FullyLeptonicAnalyzer::SelectElectron(edm::Handle<pat::ElectronCollection> eleH, vector<pat::ElectronCollection::const_iterator> & SelectedEle, int & Nele, float rho){
  for(pat::ElectronCollection::const_iterator electron = eleH->begin(); electron != eleH->end(); ++electron) {
    if(electron->pt()<40) continue;
    if(!(abs(electron->eta())<1.4442 || (abs(electron->eta())>1.5666 && abs(electron->eta())<2.5))) continue;
    if(abs(electron->phi())>3.2) continue;
    if(electron->userInt("HEEPId")!=0) continue;
    Nele = Nele + 1;
    SelectedEle.push_back(electron);
  }
}


void FullyLeptonicAnalyzer::SelectElectronCutBased(edm::Handle<pat::ElectronCollection> eleH, vector<pat::ElectronCollection::const_iterator> & SelectedEle, int & Nele, float rho, reco::Vertex primaryVertex){
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
    Nele = Nele + 1;
    SelectedEle.push_back(electron);
  }
}


void FullyLeptonicAnalyzer::SelectTightMuon(edm::Handle<pat::MuonCollection> muoH, vector<pat::MuonCollection::const_iterator> & SelectedMuo, int & Nmuo, reco::Vertex primaryVertex){
  for(pat::MuonCollection::const_iterator muon = muoH->begin(); muon != muoH->end(); ++muon) {
    if(muon->pt()<10) continue;
    if(abs(muon->eta())>2.4) continue;
    if(abs(muon->phi())>3.2) continue;
    if(!(muon->isGlobalMuon())) continue;
    if(!(muon->isPFMuon())) continue;
    if(muon->globalTrack()->normalizedChi2()>=10) continue;
    if(muon->globalTrack()->hitPattern().numberOfValidMuonHits()<=0) continue;
    if(muon->numberOfMatches()<=1) continue;
    if(fabs(muon->dB())>=0.2 ) continue;
    if(fabs(muon->muonBestTrack()->dz(primaryVertex.position()))>=0.5) continue;
    if(muon->innerTrack()->hitPattern().numberOfValidPixelHits()<=0) continue;
    if(muon->innerTrack()->hitPattern().trackerLayersWithMeasurement()<=5) continue;
    Nmuo = Nmuo + 1;
    SelectedMuo.push_back(muon);
  }
}



void FullyLeptonicAnalyzer::SelectHighptMuon(edm::Handle<pat::MuonCollection> muoH, vector<pat::MuonCollection::const_iterator> & SelectedMuo, int & Nmuo, reco::Vertex primaryVertex){
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
    SelectedMuo.push_back(muon);
    Nmuo = Nmuo + 1;
  }
}


void FullyLeptonicAnalyzer::SelectTrackerMuon(edm::Handle<pat::MuonCollection> muoH, vector<pat::MuonCollection::const_iterator> & SelectedMuo,int & Nmuo, reco::Vertex primaryVertex){
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
    Nmuo = Nmuo + 1;
  }
}

void FullyLeptonicAnalyzer::SelectTrackerGlobalID(vector<pat::MuonCollection::const_iterator> SelectedMuo, reco::Vertex primaryVertex, bool & hasAtLeastOneHighPtMuo){
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


float FullyLeptonicAnalyzer::ElectronPFIso(pat::ElectronCollection::const_iterator electron, float rho){
  float chargedHadronIso = electron->chargedHadronIso();
  float neutralHadronIso = electron->neutralHadronIso();
  float photonIso = electron->photonIso();
  float thiseta = electron->superCluster()->eta();
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

float FullyLeptonicAnalyzer::MuonPFIso(pat::MuonCollection::const_iterator muon, bool highpt){
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

float FullyLeptonicAnalyzer::MuonDETIso(pat::MuonCollection::const_iterator SelectedMuo, vector<pat::MuonCollection::const_iterator> SelectedHighptMuo){
  float isovar = SelectedMuo->trackIso()/SelectedMuo->pt();
  for(unsigned int j = 0; j< SelectedHighptMuo.size();++j){
    if(SelectedMuo==SelectedHighptMuo[j]) continue;
    double dR = ROOT::Math::VectorUtil::DeltaR(SelectedMuo->p4(),SelectedHighptMuo[j]->p4());
    if(dR < 0.3) isovar = isovar - ((SelectedHighptMuo[j]->track()->pt())/SelectedMuo->pt());
  }
  return isovar;
}

bool FullyLeptonicAnalyzer::ElectronDETIso(pat::ElectronCollection::const_iterator electron, float rho){
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

void FullyLeptonicAnalyzer::BtagVeto(edm::Handle<pat::JetCollection> ak5jetCands, int & nbtagsL, int & nbtagsM, int & nbtagsT,
				     pat::JetCollection::const_iterator SelectedJet){
  for(pat::JetCollection::const_iterator ak5 = ak5jetCands->begin(); ak5 != ak5jetCands->end(); ++ak5) {
    if(ROOT::Math::VectorUtil::DeltaR(ak5->p4(),SelectedJet->p4())<0.8) continue;
    double discCSV = ak5->bDiscriminator("combinedSecondaryVertexBJetTags");
    if(discCSV>0.244) nbtagsL++; //loose working point
    if(discCSV>0.679) nbtagsM++; //medium working point
    if(discCSV>0.898) nbtagsT++; //tight working point
  }
}

float FullyLeptonicAnalyzer::SVFitMass(edm::Handle<pat::METCollection> met, edm::Handle<pat::METCollection> metRaw, TLorentzVector & SVFitTauTau,
				       LorentzVector lep1, LorentzVector lep2){
  TMatrixD covMET(2, 2); // PFMET significance matrix
  covMET[0][0] = (metRaw->front() ).getSignificanceMatrix()(0,0);
  covMET[1][0] = (metRaw->front() ).getSignificanceMatrix()(1,0);
  covMET[0][1] = (metRaw->front() ).getSignificanceMatrix()(0,1);
  covMET[1][1] = (metRaw->front() ).getSignificanceMatrix()(1,1);
  std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, lep1));
  measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, lep2));
  NSVfitStandaloneAlgorithm algo(measuredTauLeptons, met->begin()->momentum(), covMET, 0);
  algo.addLogM(false);
  algo.integrateMarkovChain();
  SVFitTauTau.SetPtEtaPhiM(algo.pt(), algo.eta(), algo.phi(), algo.getMass());
  double MassSVFit = algo.getMass();
  return MassSVFit;
}


void FullyLeptonicAnalyzer::CollinearApproximation(edm::Handle<pat::METCollection> met, TLorentzVector & CATauTau, TLorentzVector PrunedJet, LorentzVector lep1, LorentzVector lep2, bool & CA){
  float a = (lep2.py()*met->begin()->px()-lep2.px()*met->begin()->py())/
    (lep1.px()*lep2.py()-lep1.py()*lep2.px());
  float b = (lep1.py()*met->begin()->px()-lep1.px()*met->begin()->py())/
    (lep2.px()*lep1.py()-lep2.py()*lep1.px());
  if(((1+a)*(1+b))>0) {
    CATauTau.SetPxPyPzE((1+a)*lep1.px()+(1+b)*lep2.px(),
			(1+a)*lep1.py()+(1+b)*lep2.py(),
			(1+a)*lep1.pz()+(1+b)*lep2.pz(),
			(1+a)*lep1.energy()+(1+b)*lep2.energy());
    CA=true;
  }
  else CA = false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(FullyLeptonicAnalyzer);
