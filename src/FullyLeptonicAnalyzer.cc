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
                 std::vector<pat::JetCollection::const_iterator> & SelectedJet, std::vector<float> & SelectedPrunedMass, int & Njet, float massMin, float massMax);
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

  TTree *TreeSignalEff;
  float m_genEvent;

  //TREE
  TTree *TreeEleEle;
  TTree *TreeEleMuo;
  TTree *TreeMuoMuo;
  TTree *TreeSB1EleEle;
  TTree *TreeSB1EleMuo;
  TTree *TreeSB1MuoMuo;
  TTree *TreeSB2EleEle;
  TTree *TreeSB2EleMuo;
  TTree *TreeSB2MuoMuo;
  float m_jetPt;
  float m_jetEta;
  float m_jetMass;
  float m_jetSubjettiness;
  float m_dRJetLep1;
  float m_dRJetLep2;
  float m_dRJetMet;
  float m_dPhiJetMet;
  float m_dRLepMet1;
  float m_dPhiLepMet1;
  float m_dRLepMet2;
  float m_dPhiLepMet2;
  float m_dRLepLep;
  float m_dRZZVis;
  float m_dRZZEff;
  float m_dRZZSvFit;
  float m_dRZZCA;
  float m_lepPt1;
  float m_lepEta1;
  float m_lepPFIso1;
  float m_lepDetIso1;
  float m_lepCharge1;
  float m_lepPt2;
  float m_lepEta2;
  float m_lepPFIso2;
  float m_lepDetIso2;
  float m_lepCharge2;
  int m_lepNumberEle;
  int m_lepNumberMuo;
  float m_charge;
  float m_MassVis;
  float m_MassEff;
  float m_MassSvfit;
  float m_MassCA;
  float m_XMassVis;
  float m_XMassEff;
  float m_XMassSVFit;
  float m_XMassCA;
  float m_met;
  float m_metPhi;
  float m_uncorrmet;
  float m_uncorrmetPhi;
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
  edm::InputTag metColl_;
  edm::InputTag electronColl_;
  edm::InputTag muonColl_;
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
FullyLeptonicAnalyzer::FullyLeptonicAnalyzer(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
  isData = iConfig.getUntrackedParameter<bool>("isData_");
  vtxColl_ = iConfig.getParameter<edm::InputTag>("vtxColl"); 
  jetColl_ = iConfig.getParameter<edm::InputTag>("jetColl"); 
  jetPrunedColl_ = iConfig.getParameter<edm::InputTag>("jetPrunedColl"); 
  metColl_ = iConfig.getParameter<edm::InputTag>("metColl"); 
  electronColl_ = iConfig.getParameter<edm::InputTag>("electronColl"); 
  muonColl_ = iConfig.getParameter<edm::InputTag>("muonColl"); 
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

  //SELECT THE FAT JET WITH THE HIGHTEST pt in SR
  int Njet = 0;
  vector<pat::JetCollection::const_iterator> SelectedJet;
  vector<float> SelectedPrunedMass;
  SelectJet(CA8JetswithQjets, CA8JetsPruned, SelectedJet, SelectedPrunedMass, Njet, 70., 110.);
  float jetpt=-99; int jetInd = -1;
  for(unsigned int j=0; j<SelectedJet.size(); j++){
    if(SelectedJet[j]->pt()>jetpt) {
      jetInd = j;
      jetpt = SelectedJet[j]->pt();
    }
  }

  //SELECT THE FAT JET WITH THE HIGHTEST pt in SB1
  int NSB1jet = 0;
  vector<pat::JetCollection::const_iterator> SelectedSB1Jet;
  vector<float> SelectedSB1PrunedMass;
  SelectJet(CA8JetswithQjets, CA8JetsPruned, SelectedSB1Jet, SelectedSB1PrunedMass, NSB1jet, 20., 70.);
  float jetSB1pt=-99; int jetSB1Ind = -1;
  for(unsigned int j=0; j<SelectedSB1Jet.size(); j++){
    if(SelectedSB1Jet[j]->pt()>jetSB1pt) {
      jetSB1Ind = j;
      jetSB1pt = SelectedSB1Jet[j]->pt();
    }
  }

  //SELECT THE FAT JET WITH THE HIGHTEST pt in SB2
  int NSB2jet = 0;
  vector<pat::JetCollection::const_iterator> SelectedSB2Jet;
  vector<float> SelectedSB2PrunedMass;
  SelectJet(CA8JetswithQjets, CA8JetsPruned, SelectedSB2Jet, SelectedSB2PrunedMass, NSB2jet, 110., 999999.);
  float jetSB2pt=-99; int jetSB2Ind = -1;
  for(unsigned int j=0; j<SelectedSB2Jet.size(); j++){
    if(SelectedSB2Jet[j]->pt()>jetSB2pt) {
      jetSB2Ind = j;
      jetSB2pt = SelectedSB2Jet[j]->pt();
    }
  }

  //SELECT THE ELECTRONS AND APPLY THE CUTS BETWEEN THE JET AND THE ELECTRONS in SR
  int Nele = 0;
  vector<pat::ElectronCollection::const_iterator> SelectedInitialEle;
  vector<pat::ElectronCollection::const_iterator> SelectedEle;
  SelectElectronCutBased(eleH, SelectedInitialEle, Nele, rho, primaryVertex);
  for(unsigned int i=0; i<SelectedInitialEle.size(); i++){
    if(jetInd==-1) continue;
    if(ElectronPFIso(SelectedInitialEle[i],rho)>0.1) continue;
    if(ROOT::Math::VectorUtil::DeltaR(SelectedInitialEle[i]->p4(),SelectedJet[jetInd]->p4())<0.8) continue;
    SelectedEle.push_back(SelectedInitialEle[i]);
  }

  //SELECT THE ELECTRONS AND APPLY THE CUTS BETWEEN THE JET AND THE ELECTRONS in SB1
  int NSB1ele = 0;
  vector<pat::ElectronCollection::const_iterator> SelectedSB1Ele;
  for(unsigned int i=0; i<SelectedInitialEle.size(); i++){
    if(jetSB1Ind==-1) continue;
    if(ElectronPFIso(SelectedInitialEle[i],rho)>0.1) continue;
    if(ROOT::Math::VectorUtil::DeltaR(SelectedInitialEle[i]->p4(),SelectedSB1Jet[jetSB1Ind]->p4())<0.8) continue;
    SelectedSB1Ele.push_back(SelectedInitialEle[i]);
  }

  //SELECT THE ELECTRONS AND APPLY THE CUTS BETWEEN THE JET AND THE ELECTRONS in SB2
  int NSB2ele = 0;
  vector<pat::ElectronCollection::const_iterator> SelectedSB2Ele;
  for(unsigned int i=0; i<SelectedInitialEle.size(); i++){
    if(jetSB2Ind==-1) continue;
    if(ElectronPFIso(SelectedInitialEle[i],rho)>0.1) continue;
    if(ROOT::Math::VectorUtil::DeltaR(SelectedInitialEle[i]->p4(),SelectedSB2Jet[jetSB2Ind]->p4())<0.8) continue;
    SelectedSB2Ele.push_back(SelectedInitialEle[i]);
  }
   
  //SELECT THE MUONS AND APPLY THE CUTS BETWEEN THE JET AND THE MUONS in SR
  int Nmuo = 0;
  vector<pat::MuonCollection::const_iterator> SelectedTrackerMuo;
  vector<pat::MuonCollection::const_iterator> SelectedHighptMuo;
  vector<pat::MuonCollection::const_iterator> SelectedMuo;
  SelectHighptMuon( muoH, SelectedHighptMuo, Nmuo, primaryVertex);
  SelectTrackerMuon(muoH, SelectedTrackerMuo, Nmuo, primaryVertex);
  for(unsigned int i=0; i<SelectedTrackerMuo.size(); i++){
    if(jetInd==-1) continue;
    if((MuonDETIso(SelectedTrackerMuo[i],SelectedTrackerMuo))>0.2) continue;
    if(ROOT::Math::VectorUtil::DeltaR(SelectedTrackerMuo[i]->p4(),SelectedJet[jetInd]->p4())<0.8) continue;
    SelectedMuo.push_back(SelectedTrackerMuo[i]);
  }
  bool hasAtLeastOneHighPtMuo = false;
  SelectTrackerGlobalID(SelectedMuo, primaryVertex, hasAtLeastOneHighPtMuo);
   
  //SELECT THE MUONS AND APPLY THE CUTS BETWEEN THE JET AND THE MUONS in SB1
  int NSB1muo = 0;
  vector<pat::MuonCollection::const_iterator> SelectedSB1Muo;
  for(unsigned int i=0; i<SelectedTrackerMuo.size(); i++){
    if(jetSB1Ind==-1) continue;
    if((MuonDETIso(SelectedTrackerMuo[i],SelectedTrackerMuo))>0.2) continue;
    if(ROOT::Math::VectorUtil::DeltaR(SelectedTrackerMuo[i]->p4(),SelectedSB1Jet[jetSB1Ind]->p4())<0.8) continue;
    SelectedSB1Muo.push_back(SelectedTrackerMuo[i]);
  }
  bool hasAtLeastOneHighPtMuoSB1 = false;
  SelectTrackerGlobalID(SelectedSB1Muo, primaryVertex, hasAtLeastOneHighPtMuoSB1);
   
  //SELECT THE MUONS AND APPLY THE CUTS BETWEEN THE JET AND THE MUONS in SB2
  int NSB2muo = 0;
  vector<pat::MuonCollection::const_iterator> SelectedSB2Muo;
  for(unsigned int i=0; i<SelectedTrackerMuo.size(); i++){
    if(jetSB2Ind==-1) continue;
    if((MuonDETIso(SelectedTrackerMuo[i],SelectedTrackerMuo))>0.2) continue;
    if(ROOT::Math::VectorUtil::DeltaR(SelectedTrackerMuo[i]->p4(),SelectedSB2Jet[jetSB2Ind]->p4())<0.8) continue;
    SelectedSB2Muo.push_back(SelectedTrackerMuo[i]);
  }
  bool hasAtLeastOneHighPtMuoSB2 = false;
  SelectTrackerGlobalID(SelectedSB2Muo, primaryVertex, hasAtLeastOneHighPtMuoSB2);
   
  //BTAG VETO
  int nbtagsL=0; int nbtagsM=0; int nbtagsT=0;
  if(jetInd!=-1) BtagVeto(ak5jetCands, nbtagsL, nbtagsM, nbtagsT, SelectedJet[jetInd]);
  int nSB1btagsL=0; int nSB1btagsM=0; int nSB1btagsT=0;
  if(jetSB1Ind!=-1) BtagVeto(ak5jetCands, nbtagsL, nbtagsM, nbtagsT, SelectedSB1Jet[jetSB1Ind]);
  int nSB2btagsL=0; int nSB2btagsM=0; int nSB2btagsT=0;
  if(jetSB2Ind!=-1) BtagVeto(ak5jetCands, nbtagsL, nbtagsM, nbtagsT, SelectedSB2Jet[jetSB2Ind]);

  //DI-ELECTRON SELECTION - SR
  if(SelectedEle.size()>1 && SelectedJet.size()>0){
    cout<<iEvent.id().event()<<"; electron1 pt "<<SelectedEle[0]->pt()<<"; electron2 pt "<<SelectedEle[1]->pt()<<"; jet pt "<<SelectedJet[jetInd]->pt()<<"; met "<<met->begin()->pt()<<endl;
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedJet[jetInd]->pt(),SelectedJet[jetInd]->eta(),SelectedJet[jetInd]->phi(),SelectedPrunedMass[jetInd]);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());
    //SVFIT
    TLorentzVector SVFitTauTau;
    float MassSVFit = SVFitMass(met, metRaw, SVFitTauTau, SelectedEle[0]->p4(), SelectedEle[1]->p4());
    float XMassSVFit = -1;
    float dRJetZSVFit = -1;
    if(SVFitTauTau.Pt()>0){
      XMassSVFit = (SVFitTauTau+PrunedJet).M();
      dRJetZSVFit = ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,PrunedJet);
    } 
    //COLLINEAR APPROXIMATION
    TLorentzVector CATauTau; bool CA = false;
    CollinearApproximation(met, CATauTau, PrunedJet, SelectedEle[0]->p4(), SelectedEle[1]->p4(), CA);
    float MassCA = -1;
    float XmassCA = -1.;
    float dRJetZCA = -1.;
    if(CA){
      MassCA = CATauTau.M();
      XmassCA = (CATauTau+PrunedJet).M();
      dRJetZCA = ROOT::Math::VectorUtil::DeltaR(CATauTau,PrunedJet);
    }
    pat::ElectronCollection::const_iterator SelectedLep1;
    pat::ElectronCollection::const_iterator SelectedLep2;
    float leppt=-99; unsigned int lepInd1 = -1; unsigned int lepInd2 = -1;
    for(unsigned int j=0; j<SelectedEle.size(); j++){
      if(SelectedEle[j]->pt()>leppt) {
	lepInd1 = j;
	leppt = SelectedEle[j]->pt();
      }
    }
    leppt=-99;
    for(unsigned int j=0; j<SelectedEle.size(); j++){
      if(j==lepInd1) continue;
      if(SelectedEle[j]->pt()>leppt) {
	lepInd2 = j;
	leppt = SelectedEle[j]->pt();
      }
    }
    SelectedLep1=SelectedEle[lepInd1];
    SelectedLep2=SelectedEle[lepInd2];
    math::PtEtaPhiELorentzVector dilep; dilep = SelectedLep1->p4()+SelectedLep2->p4();
    math::PtEtaPhiELorentzVector dilepmet; dilepmet = SelectedLep1->p4()+SelectedLep2->p4()+met->begin()->p4();
    math::PtEtaPhiELorentzVector dilepjet; dilepjet = SelectedLep1->p4()+SelectedLep2->p4()+PrunedJet_prov;
    math::PtEtaPhiELorentzVector dilepmetjet; dilepmetjet = SelectedLep1->p4()+SelectedLep2->p4()+met->begin()->p4()+PrunedJet_prov;
    //PLOT - TREE
    m_jetPt=SelectedJet[jetInd]->pt();
    m_jetEta=SelectedJet[jetInd]->eta();
    m_jetMass=SelectedPrunedMass[jetInd];
    m_jetSubjettiness=SelectedJet[jetInd]->userFloat("tau2")/SelectedJet[jetInd]->userFloat("tau1");
    m_dRJetLep1=ROOT::Math::VectorUtil::DeltaR(SelectedLep1->p4(),SelectedJet[jetInd]->p4());
    m_dRJetLep2=ROOT::Math::VectorUtil::DeltaR(SelectedLep2->p4(),SelectedJet[jetInd]->p4());
    m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedJet[jetInd]->p4());
    m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet[jetInd]->p4());
    m_dRLepMet1=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedLep1->p4());
    m_dPhiLepMet1=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedLep1->p4());
    m_dRLepMet2=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedLep2->p4());
    m_dPhiLepMet2=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedLep2->p4());
    m_dRLepLep=ROOT::Math::VectorUtil::DeltaR(SelectedLep1->p4(),SelectedLep2->p4());
    m_dRZZVis=ROOT::Math::VectorUtil::DeltaR(dilep,PrunedJet);
    m_dRZZEff=ROOT::Math::VectorUtil::DeltaR(dilepmet,PrunedJet);
    m_dRZZSvFit=dRJetZSVFit;
    m_dRZZCA=dRJetZCA;
    m_lepPt1=SelectedLep1->pt();
    m_lepEta1=SelectedLep1->eta();
    m_lepPFIso1=ElectronPFIso(SelectedLep1,rho);
    m_lepDetIso1=(SelectedLep1->userIso(0)+SelectedLep1->userIso(1)+SelectedLep1->userIso(2))/SelectedLep1->pt();
    m_lepCharge1=SelectedLep1->charge();
    m_lepPt2=SelectedLep2->pt();
    m_lepEta2=SelectedLep2->eta();
    m_lepPFIso2=ElectronPFIso(SelectedLep2,rho);
    m_lepDetIso2=(SelectedLep2->userIso(0)+SelectedLep2->userIso(1)+SelectedLep2->userIso(2))/SelectedLep2->pt();
    m_lepCharge2=SelectedLep2->charge();
    m_lepNumberEle=SelectedEle.size();
    m_lepNumberMuo=SelectedMuo.size();
    m_charge=SelectedLep1->charge()*SelectedLep2->charge();
    m_MassVis=dilep.mass();
    m_MassEff=dilepmet.mass();
    m_MassSvfit=MassSVFit;
    m_MassCA=MassCA;
    m_XMassVis=dilepjet.mass();
    m_XMassEff=dilepmetjet.mass();
    m_XMassSVFit=XMassSVFit;
    m_XMassCA=XmassCA;
    m_met=met->begin()->pt();
    m_metPhi=met->begin()->phi();
    m_uncorrmet=uncorrmet->begin()->pt();
    m_uncorrmetPhi=uncorrmet->begin()->phi();
    m_nbtagsL=nbtagsL;
    m_nbtagsM=nbtagsM;
    m_nbtagsT=nbtagsT;
    m_trigger320=(int)isFired_HLT_PFJet320;
    m_trigger650=(int)isFired_HLT_HT650; 
    m_trigger=(int)isFired_HLT; 
    m_sideband=(int)(SelectedPrunedMass[jetInd]<70);
    m_NVertices=vertices->size();
    m_PUWeight=MyWeight;
    m_metPx=met->begin()->px();
    m_metPy=met->begin()->py();
    m_NeventsTOT=NeventsTOT_;
    m_xsec=xsec_;
    m_lumi=lumi_;
    m_weight=xsec_*lumi_/NeventsTOT_;
    m_genEvent = genEvent;
    TreeEleEle->Fill();
  }

  //DI-MUON SELECTION - SR
  if(SelectedMuo.size()>1 && hasAtLeastOneHighPtMuo == true && SelectedJet.size()>0){
    cout<<iEvent.id().event()<<"; muon1 pt "<<SelectedMuo[0]->pt()<<"; muon2 pt "<<SelectedMuo[1]->pt()<<"; jet pt "<<SelectedJet[jetInd]->pt()<<"; met "<<met->begin()->pt()<<endl;
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedJet[jetInd]->pt(),SelectedJet[jetInd]->eta(),SelectedJet[jetInd]->phi(),SelectedPrunedMass[jetInd]);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());
    //SVFIT
    TLorentzVector SVFitTauTau;
    float MassSVFit = SVFitMass(met, metRaw, SVFitTauTau, SelectedMuo[0]->p4(), SelectedMuo[1]->p4());
    float XMassSVFit = -1;
    float dRJetZSVFit = -1;
    if(SVFitTauTau.Pt()>0){
      XMassSVFit = (SVFitTauTau+PrunedJet).M();
      dRJetZSVFit = ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,PrunedJet);
    }
    //COLLINEAR APPROXIMATION
    TLorentzVector CATauTau; bool CA = false;
    CollinearApproximation(met, CATauTau, PrunedJet, SelectedMuo[0]->p4(), SelectedMuo[1]->p4(), CA);
    float MassCA = -1;
    float XmassCA = -1.;
    float dRJetZCA = -1.;
    if(CA){
      MassCA = CATauTau.M();
      XmassCA = (CATauTau+PrunedJet).M();
      dRJetZCA = ROOT::Math::VectorUtil::DeltaR(CATauTau,PrunedJet);
    }
    pat::MuonCollection::const_iterator SelectedLep1;
    pat::MuonCollection::const_iterator SelectedLep2;
    float leppt=-99; unsigned int lepInd1 = -1; unsigned int lepInd2 = -1;
    for(unsigned int j=0; j<SelectedMuo.size(); j++){
      if(SelectedMuo[j]->pt()>leppt) {
	lepInd1 = j;
	leppt = SelectedMuo[j]->pt();
      }
    }
    leppt=-99;
    for(unsigned int j=0; j<SelectedMuo.size(); j++){
      if(j==lepInd1) continue;
      if(SelectedMuo[j]->pt()>leppt) {
	lepInd2 = j;
	leppt = SelectedMuo[j]->pt();
      }
    }
    SelectedLep1=SelectedMuo[lepInd1];
    SelectedLep2=SelectedMuo[lepInd2];
    math::PtEtaPhiELorentzVector dilep; dilep = SelectedLep1->p4()+SelectedLep2->p4();
    math::PtEtaPhiELorentzVector dilepmet; dilepmet = SelectedLep1->p4()+SelectedLep2->p4()+met->begin()->p4();
    math::PtEtaPhiELorentzVector dilepjet; dilepjet = SelectedLep1->p4()+SelectedLep2->p4()+PrunedJet_prov;
    math::PtEtaPhiELorentzVector dilepmetjet; dilepmetjet = SelectedLep1->p4()+SelectedLep2->p4()+met->begin()->p4()+PrunedJet_prov;
    //PLOT - TREE
    m_jetPt=SelectedJet[jetInd]->pt();
    m_jetEta=SelectedJet[jetInd]->eta();
    m_jetMass=SelectedPrunedMass[jetInd];
    m_jetSubjettiness=SelectedJet[jetInd]->userFloat("tau2")/SelectedJet[jetInd]->userFloat("tau1");
    m_dRJetLep1=ROOT::Math::VectorUtil::DeltaR(SelectedLep1->p4(),SelectedJet[jetInd]->p4());
    m_dRJetLep2=ROOT::Math::VectorUtil::DeltaR(SelectedLep2->p4(),SelectedJet[jetInd]->p4());
    m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedJet[jetInd]->p4());
    m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet[jetInd]->p4());
    m_dRLepMet1=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedLep1->p4());
    m_dPhiLepMet1=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedLep1->p4());
    m_dRLepMet2=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedLep2->p4());
    m_dPhiLepMet2=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedLep2->p4());
    m_dRLepLep=ROOT::Math::VectorUtil::DeltaR(SelectedLep1->p4(),SelectedLep2->p4());
    m_dRZZVis=ROOT::Math::VectorUtil::DeltaR(dilep,PrunedJet);
    m_dRZZEff=ROOT::Math::VectorUtil::DeltaR(dilepmet,PrunedJet);
    m_dRZZSvFit=dRJetZSVFit;
    m_dRZZCA=dRJetZCA;
    m_lepPt1=SelectedLep1->pt();
    m_lepEta1=SelectedLep1->eta();
    m_lepPFIso1=MuonPFIso(SelectedLep1,false);
    m_lepDetIso1=(MuonDETIso(SelectedLep1, SelectedMuo));
    m_lepCharge1=SelectedLep1->charge();
    m_lepPt2=SelectedLep2->pt();
    m_lepEta2=SelectedLep2->eta();
    m_lepPFIso2=MuonPFIso(SelectedLep2,false);
    m_lepDetIso2=(MuonDETIso(SelectedLep2, SelectedMuo));
    m_lepCharge2=SelectedLep2->charge();
    m_lepNumberEle=SelectedEle.size();
    m_lepNumberMuo=SelectedMuo.size();
    m_charge=SelectedLep1->charge()*SelectedLep2->charge();
    m_MassVis=dilep.mass();
    m_MassEff=dilepmet.mass();
    m_MassSvfit=MassSVFit;
    m_MassCA=MassCA;
    m_XMassVis=dilepjet.mass();
    m_XMassEff=dilepmetjet.mass();
    m_XMassSVFit=XMassSVFit;
    m_XMassCA=XmassCA;
    m_met=met->begin()->pt();
    m_metPhi=met->begin()->phi();
    m_uncorrmet=uncorrmet->begin()->pt();
    m_uncorrmetPhi=uncorrmet->begin()->phi();
    m_nbtagsL=nbtagsL;
    m_nbtagsM=nbtagsM;
    m_nbtagsT=nbtagsT;
    m_trigger320=(int)isFired_HLT_PFJet320;
    m_trigger650=(int)isFired_HLT_HT650;
    m_trigger=(int)isFired_HLT; 
    m_sideband=(int)(SelectedPrunedMass[jetInd]<70);
    m_NVertices=vertices->size();
    m_PUWeight=MyWeight;
    m_metPx=met->begin()->px();
    m_metPy=met->begin()->py();
    m_NeventsTOT=NeventsTOT_;
    m_xsec=xsec_;
    m_lumi=lumi_;
    m_weight=xsec_*lumi_/NeventsTOT_;
    m_genEvent = genEvent;
    TreeMuoMuo->Fill();
  }

  //ELECTRON-MUON SELECTION - SR
  SelectedEle.clear(); SelectedMuo.clear();
  for(unsigned int i=0; i<SelectedInitialEle.size(); i++){
    if(jetInd==-1) continue;
    if(ROOT::Math::VectorUtil::DeltaR(SelectedInitialEle[i]->p4(),SelectedJet[jetInd]->p4())<0.8) continue;
    //if(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedInitialEle[i]->p4())>1) continue;
    if(ElectronPFIso(SelectedInitialEle[i],rho)>0.1) continue;
    SelectedEle.push_back(SelectedInitialEle[i]);
  }
  for(unsigned int i=0; i<SelectedHighptMuo.size(); i++){
    if(jetInd==-1) continue;
    if(ROOT::Math::VectorUtil::DeltaR(SelectedHighptMuo[i]->p4(),SelectedJet[jetInd]->p4())<0.8) continue;
    //if(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedHighptMuo[i]->p4())>1) continue;
    if(MuonPFIso(SelectedHighptMuo[i], true)>0.2) continue;
    SelectedMuo.push_back(SelectedHighptMuo[i]);
  }
  if(SelectedEle.size()>0 && SelectedMuo.size()>0 && SelectedJet.size()>0){
    cout<<iEvent.id().event()<<"; muon pt "<<SelectedMuo[0]->pt()<<"; electron pt "<<SelectedEle[0]->pt()<<"; jet pt "<<SelectedJet[jetInd]->pt()<<"; met "<<met->begin()->pt()<<endl;
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedJet[jetInd]->pt(),SelectedJet[jetInd]->eta(),SelectedJet[jetInd]->phi(),SelectedPrunedMass[jetInd]);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());
    //SVFIT
    TLorentzVector SVFitTauTau;
    float MassSVFit = SVFitMass(met, metRaw, SVFitTauTau, SelectedEle[0]->p4(), SelectedMuo[0]->p4());
    float XMassSVFit = -1;
    float dRJetZSVFit = -1;
    if(SVFitTauTau.Pt()>0){
      XMassSVFit = (SVFitTauTau+PrunedJet).M();
      dRJetZSVFit = ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,PrunedJet);
    }
    //COLLINEAR APPROXIMATION
    TLorentzVector CATauTau; bool CA = false;
    CollinearApproximation(met, CATauTau, PrunedJet, SelectedEle[0]->p4(), SelectedMuo[0]->p4(), CA);
    float MassCA = -1;
    float XmassCA = -1.;
    float dRJetZCA = -1.;
    if(CA){
      MassCA = CATauTau.M();
      XmassCA = (CATauTau+PrunedJet).M();
      dRJetZCA = ROOT::Math::VectorUtil::DeltaR(CATauTau,PrunedJet);
    }
    pat::ElectronCollection::const_iterator SelectedLep1;
    pat::MuonCollection::const_iterator SelectedLep2;
    float leppt=-99; int lepInd1 = -1; int lepInd2 = -1;
    for(unsigned int j=0; j<SelectedEle.size(); j++){
      if(SelectedEle[j]->pt()>leppt) {
	lepInd1 = j;
	leppt = SelectedEle[j]->pt();
      }
    }
    leppt=-99;
    for(unsigned int j=0; j<SelectedMuo.size(); j++){
      if(SelectedMuo[j]->pt()>leppt) {
	lepInd2 = j;
	leppt = SelectedMuo[j]->pt();
      }
    }
    SelectedLep1=SelectedEle[lepInd1];
    SelectedLep2=SelectedMuo[lepInd2];
    math::PtEtaPhiELorentzVector dilep; dilep = SelectedLep1->p4()+SelectedLep2->p4();
    math::PtEtaPhiELorentzVector dilepmet; dilepmet = SelectedLep1->p4()+SelectedLep2->p4()+met->begin()->p4();
    math::PtEtaPhiELorentzVector dilepjet; dilepjet = SelectedLep1->p4()+SelectedLep2->p4()+PrunedJet_prov;
    math::PtEtaPhiELorentzVector dilepmetjet; dilepmetjet = SelectedLep1->p4()+SelectedLep2->p4()+met->begin()->p4()+PrunedJet_prov;
    //PLOT - TREE
    m_jetPt=SelectedJet[jetInd]->pt();
    m_jetEta=SelectedJet[jetInd]->eta();
    m_jetMass=SelectedPrunedMass[jetInd];
    m_jetSubjettiness=SelectedJet[jetInd]->userFloat("tau2")/SelectedJet[jetInd]->userFloat("tau1");
    m_dRJetLep1=ROOT::Math::VectorUtil::DeltaR(SelectedLep1->p4(),SelectedJet[jetInd]->p4());
    m_dRJetLep2=ROOT::Math::VectorUtil::DeltaR(SelectedLep2->p4(),SelectedJet[jetInd]->p4());
    m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedJet[jetInd]->p4());
    m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet[jetInd]->p4());
    m_dRLepMet1=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedLep1->p4());
    m_dPhiLepMet1=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedLep1->p4());
    m_dRLepMet2=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedLep2->p4());
    m_dPhiLepMet2=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedLep2->p4());
    m_dRLepLep=ROOT::Math::VectorUtil::DeltaR(SelectedLep1->p4(),SelectedLep2->p4());
    m_dRZZVis=ROOT::Math::VectorUtil::DeltaR(dilep,PrunedJet);
    m_dRZZEff=ROOT::Math::VectorUtil::DeltaR(dilepmet,PrunedJet);
    m_dRZZSvFit=dRJetZSVFit;
    m_dRZZCA=dRJetZCA;
    m_lepPt1=SelectedLep1->pt();
    m_lepEta1=SelectedLep1->eta();
    m_lepPFIso1=ElectronPFIso(SelectedLep1,rho);
    m_lepDetIso1=(SelectedLep1->userIso(0)+SelectedLep1->userIso(1)+SelectedLep1->userIso(2))/SelectedLep1->pt();
    m_lepCharge1=SelectedLep1->charge();
    m_lepPt2=SelectedLep2->pt();
    m_lepEta2=SelectedLep2->eta();
    m_lepPFIso2=MuonPFIso(SelectedLep2,false);
    m_lepDetIso2=(MuonDETIso(SelectedLep2, SelectedMuo));
    m_lepCharge2=SelectedLep2->charge();
    m_lepNumberEle=SelectedEle.size();
    m_lepNumberMuo=SelectedMuo.size();
    m_charge=SelectedLep1->charge()*SelectedLep2->charge();
    m_MassVis=dilep.mass();
    m_MassEff=dilepmet.mass();
    m_MassSvfit=MassSVFit;
    m_MassCA=MassCA;
    m_XMassVis=dilepjet.mass();
    m_XMassEff=dilepmetjet.mass();
    m_XMassSVFit=XMassSVFit;
    m_XMassCA=XmassCA;
    m_met=met->begin()->pt();
    m_metPhi=met->begin()->phi();
    m_uncorrmet=uncorrmet->begin()->pt();
    m_uncorrmetPhi=uncorrmet->begin()->phi();
    m_nbtagsL=nbtagsL;
    m_nbtagsM=nbtagsM;
    m_nbtagsT=nbtagsT;
    m_trigger320=(int)isFired_HLT_PFJet320;
    m_trigger650=(int)isFired_HLT_HT650;
    m_trigger=(int)isFired_HLT; 
    m_sideband=(int)(SelectedPrunedMass[jetInd]<70);
    m_NVertices=vertices->size();
    m_PUWeight=MyWeight;
    m_metPx=met->begin()->px();
    m_metPy=met->begin()->py();
    m_NeventsTOT=NeventsTOT_;
    m_xsec=xsec_;
    m_lumi=lumi_;
    m_weight=xsec_*lumi_/NeventsTOT_;
    m_genEvent = genEvent;
    TreeEleMuo->Fill();
  }

  //DI-ELECTRON SELECTION - SB1
  if(SelectedSB1Ele.size()>1 && SelectedSB1Jet.size()>0){
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedSB1Jet[jetSB1Ind]->pt(),SelectedSB1Jet[jetSB1Ind]->eta(),SelectedSB1Jet[jetSB1Ind]->phi(),
						SelectedSB1PrunedMass[jetSB1Ind]);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());
    //SVFIT
    TLorentzVector SVFitTauTau;
    float MassSVFit = SVFitMass(met, metRaw, SVFitTauTau, SelectedSB1Ele[0]->p4(), SelectedSB1Ele[1]->p4());
    float XMassSVFit = -1;
    float dRJetZSVFit = -1;
    if(SVFitTauTau.Pt()>0){
      XMassSVFit = (SVFitTauTau+PrunedJet).M();
      dRJetZSVFit = ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,PrunedJet);
    } 
    //COLLINEAR APPROXIMATION
    TLorentzVector CATauTau; bool CA = false;
    CollinearApproximation(met, CATauTau, PrunedJet, SelectedSB1Ele[0]->p4(), SelectedSB1Ele[1]->p4(), CA);
    float MassCA = -1;
    float XmassCA = -1.;
    float dRJetZCA = -1.;
    if(CA){
      MassCA = CATauTau.M();
      XmassCA = (CATauTau+PrunedJet).M();
      dRJetZCA = ROOT::Math::VectorUtil::DeltaR(CATauTau,PrunedJet);
    }
    pat::ElectronCollection::const_iterator SelectedSB1Lep1;
    pat::ElectronCollection::const_iterator SelectedSB1Lep2;
    float leppt=-99; unsigned int lepInd1 = -1; unsigned int lepInd2 = -1;
    for(unsigned int j=0; j<SelectedSB1Ele.size(); j++){
      if(SelectedSB1Ele[j]->pt()>leppt) {
	lepInd1 = j;
	leppt = SelectedSB1Ele[j]->pt();
      }
    }
    leppt=-99;
    for(unsigned int j=0; j<SelectedSB1Ele.size(); j++){
      if(j==lepInd1) continue;
      if(SelectedSB1Ele[j]->pt()>leppt) {
	lepInd2 = j;
	leppt = SelectedSB1Ele[j]->pt();
      }
    }
    SelectedSB1Lep1=SelectedSB1Ele[lepInd1];
    SelectedSB1Lep2=SelectedSB1Ele[lepInd2];
    math::PtEtaPhiELorentzVector dilep; dilep = SelectedSB1Lep1->p4()+SelectedSB1Lep2->p4();
    math::PtEtaPhiELorentzVector dilepmet; dilepmet = SelectedSB1Lep1->p4()+SelectedSB1Lep2->p4()+met->begin()->p4();
    math::PtEtaPhiELorentzVector dilepjet; dilepjet = SelectedSB1Lep1->p4()+SelectedSB1Lep2->p4()+PrunedJet_prov;
    math::PtEtaPhiELorentzVector dilepmetjet; dilepmetjet = SelectedSB1Lep1->p4()+SelectedSB1Lep2->p4()+met->begin()->p4()+PrunedJet_prov;
    //PLOT - TREE
    m_jetPt=SelectedSB1Jet[jetSB1Ind]->pt();
    m_jetEta=SelectedSB1Jet[jetSB1Ind]->eta();
    m_jetMass=SelectedSB1PrunedMass[jetSB1Ind];
    m_jetSubjettiness=SelectedSB1Jet[jetSB1Ind]->userFloat("tau2")/SelectedSB1Jet[jetSB1Ind]->userFloat("tau1");
    m_dRJetLep1=ROOT::Math::VectorUtil::DeltaR(SelectedSB1Lep1->p4(),SelectedSB1Jet[jetSB1Ind]->p4());
    m_dRJetLep2=ROOT::Math::VectorUtil::DeltaR(SelectedSB1Lep2->p4(),SelectedSB1Jet[jetSB1Ind]->p4());
    m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB1Jet[jetSB1Ind]->p4());
    m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB1Jet[jetSB1Ind]->p4());
    m_dRLepMet1=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB1Lep1->p4());
    m_dPhiLepMet1=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB1Lep1->p4());
    m_dRLepMet2=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB1Lep2->p4());
    m_dPhiLepMet2=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB1Lep2->p4());
    m_dRLepLep=ROOT::Math::VectorUtil::DeltaR(SelectedSB1Lep1->p4(),SelectedSB1Lep2->p4());
    m_dRZZVis=ROOT::Math::VectorUtil::DeltaR(dilep,PrunedJet);
    m_dRZZEff=ROOT::Math::VectorUtil::DeltaR(dilepmet,PrunedJet);
    m_dRZZSvFit=dRJetZSVFit;
    m_dRZZCA=dRJetZCA;
    m_lepPt1=SelectedSB1Lep1->pt();
    m_lepEta1=SelectedSB1Lep1->eta();
    m_lepPFIso1=ElectronPFIso(SelectedSB1Lep1,rho);
    m_lepDetIso1=(SelectedSB1Lep1->userIso(0)+SelectedSB1Lep1->userIso(1)+SelectedSB1Lep1->userIso(2))/SelectedSB1Lep1->pt();
    m_lepCharge1=SelectedSB1Lep1->charge();
    m_lepPt2=SelectedSB1Lep2->pt();
    m_lepEta2=SelectedSB1Lep2->eta();
    m_lepPFIso2=ElectronPFIso(SelectedSB1Lep2,rho);
    m_lepDetIso2=(SelectedSB1Lep2->userIso(0)+SelectedSB1Lep2->userIso(1)+SelectedSB1Lep2->userIso(2))/SelectedSB1Lep2->pt();
    m_lepCharge2=SelectedSB1Lep2->charge();
    m_lepNumberEle=SelectedSB1Ele.size();
    m_lepNumberMuo=SelectedSB1Muo.size();
    m_charge=SelectedSB1Lep1->charge()*SelectedSB1Lep2->charge();
    m_MassVis=dilep.mass();
    m_MassEff=dilepmet.mass();
    m_MassSvfit=MassSVFit;
    m_MassCA=MassCA;
    m_XMassVis=dilepjet.mass();
    m_XMassEff=dilepmetjet.mass();
    m_XMassSVFit=XMassSVFit;
    m_XMassCA=XmassCA;
    m_met=met->begin()->pt();
    m_metPhi=met->begin()->phi();
    m_uncorrmet=uncorrmet->begin()->pt();
    m_uncorrmetPhi=uncorrmet->begin()->phi();
    m_nbtagsL=nSB1btagsL;
    m_nbtagsM=nSB1btagsM;
    m_nbtagsT=nSB1btagsT;
    m_trigger320=(int)isFired_HLT_PFJet320;
    m_trigger650=(int)isFired_HLT_HT650; 
    m_trigger=(int)isFired_HLT; 
    m_sideband=(int)(SelectedSB1PrunedMass[jetSB1Ind]<70);
    m_NVertices=vertices->size();
    m_PUWeight=MyWeight;
    m_metPx=met->begin()->px();
    m_metPy=met->begin()->py();
    m_NeventsTOT=NeventsTOT_;
    m_xsec=xsec_;
    m_lumi=lumi_;
    m_weight=xsec_*lumi_/NeventsTOT_;
    m_genEvent = genEvent;
    TreeSB1EleEle->Fill();
  }

  //DI-MUON SELECTION - SB1
  if(SelectedSB1Muo.size()>1 && hasAtLeastOneHighPtMuoSB1 == true && SelectedSB1Jet.size()>0){
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedSB1Jet[jetSB1Ind]->pt(),SelectedSB1Jet[jetSB1Ind]->eta(),SelectedSB1Jet[jetSB1Ind]->phi(),
						SelectedSB1PrunedMass[jetSB1Ind]);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());
    //SVFIT
    TLorentzVector SVFitTauTau;
    float MassSVFit = SVFitMass(met, metRaw, SVFitTauTau, SelectedSB1Muo[0]->p4(), SelectedSB1Muo[1]->p4());
    float XMassSVFit = -1;
    float dRJetZSVFit = -1;
    if(SVFitTauTau.Pt()>0){
      XMassSVFit = (SVFitTauTau+PrunedJet).M();
      dRJetZSVFit = ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,PrunedJet);
    }
    //COLLINEAR APPROXIMATION
    TLorentzVector CATauTau; bool CA = false;
    CollinearApproximation(met, CATauTau, PrunedJet, SelectedSB1Muo[0]->p4(), SelectedSB1Muo[1]->p4(), CA);
    float MassCA = -1;
    float XmassCA = -1.;
    float dRJetZCA = -1.;
    if(CA){
      MassCA = CATauTau.M();
      XmassCA = (CATauTau+PrunedJet).M();
      dRJetZCA = ROOT::Math::VectorUtil::DeltaR(CATauTau,PrunedJet);
    }
    pat::MuonCollection::const_iterator SelectedSB1Lep1;
    pat::MuonCollection::const_iterator SelectedSB1Lep2;
    float leppt=-99; unsigned int lepInd1 = -1; unsigned int lepInd2 = -1;
    for(unsigned int j=0; j<SelectedSB1Muo.size(); j++){
      if(SelectedSB1Muo[j]->pt()>leppt) {
	lepInd1 = j;
	leppt = SelectedSB1Muo[j]->pt();
      }
    }
    leppt=-99;
    for(unsigned int j=0; j<SelectedSB1Muo.size(); j++){
      if(j==lepInd1) continue;
      if(SelectedSB1Muo[j]->pt()>leppt) {
	lepInd2 = j;
	leppt = SelectedSB1Muo[j]->pt();
      }
    }
    SelectedSB1Lep1=SelectedSB1Muo[lepInd1];
    SelectedSB1Lep2=SelectedSB1Muo[lepInd2];
    math::PtEtaPhiELorentzVector dilep; dilep = SelectedSB1Lep1->p4()+SelectedSB1Lep2->p4();
    math::PtEtaPhiELorentzVector dilepmet; dilepmet = SelectedSB1Lep1->p4()+SelectedSB1Lep2->p4()+met->begin()->p4();
    math::PtEtaPhiELorentzVector dilepjet; dilepjet = SelectedSB1Lep1->p4()+SelectedSB1Lep2->p4()+PrunedJet_prov;
    math::PtEtaPhiELorentzVector dilepmetjet; dilepmetjet = SelectedSB1Lep1->p4()+SelectedSB1Lep2->p4()+met->begin()->p4()+PrunedJet_prov;
    //PLOT - TREE
    m_jetPt=SelectedSB1Jet[jetSB1Ind]->pt();
    m_jetEta=SelectedSB1Jet[jetSB1Ind]->eta();
    m_jetMass=SelectedSB1PrunedMass[jetSB1Ind];
    m_jetSubjettiness=SelectedSB1Jet[jetSB1Ind]->userFloat("tau2")/SelectedSB1Jet[jetSB1Ind]->userFloat("tau1");
    m_dRJetLep1=ROOT::Math::VectorUtil::DeltaR(SelectedSB1Lep1->p4(),SelectedSB1Jet[jetSB1Ind]->p4());
    m_dRJetLep2=ROOT::Math::VectorUtil::DeltaR(SelectedSB1Lep2->p4(),SelectedSB1Jet[jetSB1Ind]->p4());
    m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB1Jet[jetSB1Ind]->p4());
    m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB1Jet[jetSB1Ind]->p4());
    m_dRLepMet1=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB1Lep1->p4());
    m_dPhiLepMet1=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB1Lep1->p4());
    m_dRLepMet2=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB1Lep2->p4());
    m_dPhiLepMet2=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB1Lep2->p4());
    m_dRLepLep=ROOT::Math::VectorUtil::DeltaR(SelectedSB1Lep1->p4(),SelectedSB1Lep2->p4());
    m_dRZZVis=ROOT::Math::VectorUtil::DeltaR(dilep,PrunedJet);
    m_dRZZEff=ROOT::Math::VectorUtil::DeltaR(dilepmet,PrunedJet);
    m_dRZZSvFit=dRJetZSVFit;
    m_dRZZCA=dRJetZCA;
    m_lepPt1=SelectedSB1Lep1->pt();
    m_lepEta1=SelectedSB1Lep1->eta();
    m_lepPFIso1=MuonPFIso(SelectedSB1Lep1,false);
    m_lepDetIso1=(MuonDETIso(SelectedSB1Lep1, SelectedSB1Muo));
    m_lepCharge1=SelectedSB1Lep1->charge();
    m_lepPt2=SelectedSB1Lep2->pt();
    m_lepEta2=SelectedSB1Lep2->eta();
    m_lepPFIso2=MuonPFIso(SelectedSB1Lep2,false);
    m_lepDetIso2=(MuonDETIso(SelectedSB1Lep2, SelectedSB1Muo));
    m_lepCharge2=SelectedSB1Lep2->charge();
    m_lepNumberEle=SelectedSB1Ele.size();
    m_lepNumberMuo=SelectedSB1Muo.size();
    m_charge=SelectedSB1Lep1->charge()*SelectedSB1Lep2->charge();
    m_MassVis=dilep.mass();
    m_MassEff=dilepmet.mass();
    m_MassSvfit=MassSVFit;
    m_MassCA=MassCA;
    m_XMassVis=dilepjet.mass();
    m_XMassEff=dilepmetjet.mass();
    m_XMassSVFit=XMassSVFit;
    m_XMassCA=XmassCA;
    m_met=met->begin()->pt();
    m_metPhi=met->begin()->phi();
    m_uncorrmet=uncorrmet->begin()->pt();
    m_uncorrmetPhi=uncorrmet->begin()->phi();
    m_nbtagsL=nSB1btagsL;
    m_nbtagsM=nSB1btagsM;
    m_nbtagsT=nSB1btagsT;
    m_trigger320=(int)isFired_HLT_PFJet320;
    m_trigger650=(int)isFired_HLT_HT650;
    m_trigger=(int)isFired_HLT; 
    m_sideband=(int)(SelectedSB1PrunedMass[jetSB1Ind]<70);
    m_NVertices=vertices->size();
    m_PUWeight=MyWeight;
    m_metPx=met->begin()->px();
    m_metPy=met->begin()->py();
    m_NeventsTOT=NeventsTOT_;
    m_xsec=xsec_;
    m_lumi=lumi_;
    m_weight=xsec_*lumi_/NeventsTOT_;
    m_genEvent = genEvent;
    TreeSB1MuoMuo->Fill();
  }

  //ELECTRON-MUON SELECTION - SB1
  SelectedSB1Ele.clear(); SelectedSB1Muo.clear();
  for(unsigned int i=0; i<SelectedInitialEle.size(); i++){
    if(jetSB1Ind==-1) continue;
    if(ROOT::Math::VectorUtil::DeltaR(SelectedInitialEle[i]->p4(),SelectedSB1Jet[jetSB1Ind]->p4())<0.8) continue;
    //if(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedInitialEle[i]->p4())>1) continue;
    if(ElectronPFIso(SelectedInitialEle[i],rho)>0.1) continue;
    SelectedSB1Ele.push_back(SelectedInitialEle[i]);
  }
  for(unsigned int i=0; i<SelectedHighptMuo.size(); i++){
    if(jetSB1Ind==-1) continue;
    if(ROOT::Math::VectorUtil::DeltaR(SelectedHighptMuo[i]->p4(),SelectedSB1Jet[jetSB1Ind]->p4())<0.8) continue;
    //if(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedHighptMuo[i]->p4())>1) continue;
    if(MuonPFIso(SelectedHighptMuo[i], true)>0.2) continue;
    SelectedSB1Muo.push_back(SelectedHighptMuo[i]);
  }
  if(SelectedSB1Ele.size()>0 && SelectedSB1Muo.size()>0 && SelectedSB1Jet.size()>0){
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedSB1Jet[jetSB1Ind]->pt(),SelectedSB1Jet[jetSB1Ind]->eta(),SelectedSB1Jet[jetSB1Ind]->phi(),
						SelectedSB1PrunedMass[jetSB1Ind]);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());
    //SVFIT
    TLorentzVector SVFitTauTau;
    float MassSVFit = SVFitMass(met, metRaw, SVFitTauTau, SelectedSB1Ele[0]->p4(), SelectedSB1Muo[0]->p4());
    float XMassSVFit = -1;
    float dRJetZSVFit = -1;
    if(SVFitTauTau.Pt()>0){
      XMassSVFit = (SVFitTauTau+PrunedJet).M();
      dRJetZSVFit = ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,PrunedJet);
    }
    //COLLINEAR APPROXIMATION
    TLorentzVector CATauTau; bool CA = false;
    CollinearApproximation(met, CATauTau, PrunedJet, SelectedSB1Ele[0]->p4(), SelectedSB1Muo[0]->p4(), CA);
    float MassCA = -1;
    float XmassCA = -1.;
    float dRJetZCA = -1.;
    if(CA){
      MassCA = CATauTau.M();
      XmassCA = (CATauTau+PrunedJet).M();
      dRJetZCA = ROOT::Math::VectorUtil::DeltaR(CATauTau,PrunedJet);
    }
    pat::ElectronCollection::const_iterator SelectedSB1Lep1;
    pat::MuonCollection::const_iterator SelectedSB1Lep2;
    float leppt=-99; int lepInd1 = -1; int lepInd2 = -1;
    for(unsigned int j=0; j<SelectedSB1Ele.size(); j++){
      if(SelectedSB1Ele[j]->pt()>leppt) {
	lepInd1 = j;
	leppt = SelectedSB1Ele[j]->pt();
      }
    }
    leppt=-99;
    for(unsigned int j=0; j<SelectedSB1Muo.size(); j++){
      if(SelectedSB1Muo[j]->pt()>leppt) {
	lepInd2 = j;
	leppt = SelectedSB1Muo[j]->pt();
      }
    }
    SelectedSB1Lep1=SelectedSB1Ele[lepInd1];
    SelectedSB1Lep2=SelectedSB1Muo[lepInd2];
    math::PtEtaPhiELorentzVector dilep; dilep = SelectedSB1Lep1->p4()+SelectedSB1Lep2->p4();
    math::PtEtaPhiELorentzVector dilepmet; dilepmet = SelectedSB1Lep1->p4()+SelectedSB1Lep2->p4()+met->begin()->p4();
    math::PtEtaPhiELorentzVector dilepjet; dilepjet = SelectedSB1Lep1->p4()+SelectedSB1Lep2->p4()+PrunedJet_prov;
    math::PtEtaPhiELorentzVector dilepmetjet; dilepmetjet = SelectedSB1Lep1->p4()+SelectedSB1Lep2->p4()+met->begin()->p4()+PrunedJet_prov;
    //PLOT - TREE
    m_jetPt=SelectedSB1Jet[jetSB1Ind]->pt();
    m_jetEta=SelectedSB1Jet[jetSB1Ind]->eta();
    m_jetMass=SelectedSB1PrunedMass[jetSB1Ind];
    m_jetSubjettiness=SelectedSB1Jet[jetSB1Ind]->userFloat("tau2")/SelectedSB1Jet[jetSB1Ind]->userFloat("tau1");
    m_dRJetLep1=ROOT::Math::VectorUtil::DeltaR(SelectedSB1Lep1->p4(),SelectedSB1Jet[jetSB1Ind]->p4());
    m_dRJetLep2=ROOT::Math::VectorUtil::DeltaR(SelectedSB1Lep2->p4(),SelectedSB1Jet[jetSB1Ind]->p4());
    m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB1Jet[jetSB1Ind]->p4());
    m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB1Jet[jetSB1Ind]->p4());
    m_dRLepMet1=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB1Lep1->p4());
    m_dPhiLepMet1=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB1Lep1->p4());
    m_dRLepMet2=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB1Lep2->p4());
    m_dPhiLepMet2=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB1Lep2->p4());
    m_dRLepLep=ROOT::Math::VectorUtil::DeltaR(SelectedSB1Lep1->p4(),SelectedSB1Lep2->p4());
    m_dRZZVis=ROOT::Math::VectorUtil::DeltaR(dilep,PrunedJet);
    m_dRZZEff=ROOT::Math::VectorUtil::DeltaR(dilepmet,PrunedJet);
    m_dRZZSvFit=dRJetZSVFit;
    m_dRZZCA=dRJetZCA;
    m_lepPt1=SelectedSB1Lep1->pt();
    m_lepEta1=SelectedSB1Lep1->eta();
    m_lepPFIso1=ElectronPFIso(SelectedSB1Lep1,rho);
    m_lepDetIso1=(SelectedSB1Lep1->userIso(0)+SelectedSB1Lep1->userIso(1)+SelectedSB1Lep1->userIso(2))/SelectedSB1Lep1->pt();
    m_lepCharge1=SelectedSB1Lep1->charge();
    m_lepPt2=SelectedSB1Lep2->pt();
    m_lepEta2=SelectedSB1Lep2->eta();
    m_lepPFIso2=MuonPFIso(SelectedSB1Lep2,false);
    m_lepDetIso2=(MuonDETIso(SelectedSB1Lep2, SelectedSB1Muo));
    m_lepCharge2=SelectedSB1Lep2->charge();
    m_lepNumberEle=SelectedSB1Ele.size();
    m_lepNumberMuo=SelectedSB1Muo.size();
    m_charge=SelectedSB1Lep1->charge()*SelectedSB1Lep2->charge();
    m_MassVis=dilep.mass();
    m_MassEff=dilepmet.mass();
    m_MassSvfit=MassSVFit;
    m_MassCA=MassCA;
    m_XMassVis=dilepjet.mass();
    m_XMassEff=dilepmetjet.mass();
    m_XMassSVFit=XMassSVFit;
    m_XMassCA=XmassCA;
    m_met=met->begin()->pt();
    m_metPhi=met->begin()->phi();
    m_uncorrmet=uncorrmet->begin()->pt();
    m_uncorrmetPhi=uncorrmet->begin()->phi();
    m_nbtagsL=nSB1btagsL;
    m_nbtagsM=nSB1btagsM;
    m_nbtagsT=nSB1btagsT;
    m_trigger320=(int)isFired_HLT_PFJet320;
    m_trigger650=(int)isFired_HLT_HT650;
    m_trigger=(int)isFired_HLT; 
    m_sideband=(int)(SelectedSB1PrunedMass[jetSB1Ind]<70);
    m_NVertices=vertices->size();
    m_PUWeight=MyWeight;
    m_metPx=met->begin()->px();
    m_metPy=met->begin()->py();
    m_NeventsTOT=NeventsTOT_;
    m_xsec=xsec_;
    m_lumi=lumi_;
    m_weight=xsec_*lumi_/NeventsTOT_;
    m_genEvent = genEvent;
    TreeSB1EleMuo->Fill();
  }

  //DI-ELECTRON SELECTION - SB2
  if(SelectedSB2Ele.size()>1 && SelectedSB2Jet.size()>0){
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedSB2Jet[jetSB2Ind]->pt(),SelectedSB2Jet[jetSB2Ind]->eta(),SelectedSB2Jet[jetSB2Ind]->phi(),
						SelectedSB2PrunedMass[jetSB2Ind]);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());
    //SVFIT
    TLorentzVector SVFitTauTau;
    float MassSVFit = SVFitMass(met, metRaw, SVFitTauTau, SelectedSB2Ele[0]->p4(), SelectedSB2Ele[1]->p4());
    float XMassSVFit = -1;
    float dRJetZSVFit = -1;
    if(SVFitTauTau.Pt()>0){
      XMassSVFit = (SVFitTauTau+PrunedJet).M();
      dRJetZSVFit = ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,PrunedJet);
    } 
    //COLLINEAR APPROXIMATION
    TLorentzVector CATauTau; bool CA = false;
    CollinearApproximation(met, CATauTau, PrunedJet, SelectedSB2Ele[0]->p4(), SelectedSB2Ele[1]->p4(), CA);
    float MassCA = -1;
    float XmassCA = -1.;
    float dRJetZCA = -1.;
    if(CA){
      MassCA = CATauTau.M();
      XmassCA = (CATauTau+PrunedJet).M();
      dRJetZCA = ROOT::Math::VectorUtil::DeltaR(CATauTau,PrunedJet);
    }
    pat::ElectronCollection::const_iterator SelectedSB2Lep1;
    pat::ElectronCollection::const_iterator SelectedSB2Lep2;
    float leppt=-99; unsigned int lepInd1 = -1; unsigned int lepInd2 = -1;
    for(unsigned int j=0; j<SelectedSB2Ele.size(); j++){
      if(SelectedSB2Ele[j]->pt()>leppt) {
	lepInd1 = j;
	leppt = SelectedSB2Ele[j]->pt();
      }
    }
    leppt=-99;
    for(unsigned int j=0; j<SelectedSB2Ele.size(); j++){
      if(j==lepInd1) continue;
      if(SelectedSB2Ele[j]->pt()>leppt) {
	lepInd2 = j;
	leppt = SelectedSB2Ele[j]->pt();
      }
    }
    SelectedSB2Lep1=SelectedSB2Ele[lepInd1];
    SelectedSB2Lep2=SelectedSB2Ele[lepInd2];
    math::PtEtaPhiELorentzVector dilep; dilep = SelectedSB2Lep1->p4()+SelectedSB2Lep2->p4();
    math::PtEtaPhiELorentzVector dilepmet; dilepmet = SelectedSB2Lep1->p4()+SelectedSB2Lep2->p4()+met->begin()->p4();
    math::PtEtaPhiELorentzVector dilepjet; dilepjet = SelectedSB2Lep1->p4()+SelectedSB2Lep2->p4()+PrunedJet_prov;
    math::PtEtaPhiELorentzVector dilepmetjet; dilepmetjet = SelectedSB2Lep1->p4()+SelectedSB2Lep2->p4()+met->begin()->p4()+PrunedJet_prov;
    //PLOT - TREE
    m_jetPt=SelectedSB2Jet[jetSB2Ind]->pt();
    m_jetEta=SelectedSB2Jet[jetSB2Ind]->eta();
    m_jetMass=SelectedSB2PrunedMass[jetSB2Ind];
    m_jetSubjettiness=SelectedSB2Jet[jetSB2Ind]->userFloat("tau2")/SelectedSB2Jet[jetSB2Ind]->userFloat("tau1");
    m_dRJetLep1=ROOT::Math::VectorUtil::DeltaR(SelectedSB2Lep1->p4(),SelectedSB2Jet[jetSB2Ind]->p4());
    m_dRJetLep2=ROOT::Math::VectorUtil::DeltaR(SelectedSB2Lep2->p4(),SelectedSB2Jet[jetSB2Ind]->p4());
    m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB2Jet[jetSB2Ind]->p4());
    m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB2Jet[jetSB2Ind]->p4());
    m_dRLepMet1=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB2Lep1->p4());
    m_dPhiLepMet1=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB2Lep1->p4());
    m_dRLepMet2=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB2Lep2->p4());
    m_dPhiLepMet2=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB2Lep2->p4());
    m_dRLepLep=ROOT::Math::VectorUtil::DeltaR(SelectedSB2Lep1->p4(),SelectedSB2Lep2->p4());
    m_dRZZVis=ROOT::Math::VectorUtil::DeltaR(dilep,PrunedJet);
    m_dRZZEff=ROOT::Math::VectorUtil::DeltaR(dilepmet,PrunedJet);
    m_dRZZSvFit=dRJetZSVFit;
    m_dRZZCA=dRJetZCA;
    m_lepPt1=SelectedSB2Lep1->pt();
    m_lepEta1=SelectedSB2Lep1->eta();
    m_lepPFIso1=ElectronPFIso(SelectedSB2Lep1,rho);
    m_lepDetIso1=(SelectedSB2Lep1->userIso(0)+SelectedSB2Lep1->userIso(1)+SelectedSB2Lep1->userIso(2))/SelectedSB2Lep1->pt();
    m_lepCharge1=SelectedSB2Lep1->charge();
    m_lepPt2=SelectedSB2Lep2->pt();
    m_lepEta2=SelectedSB2Lep2->eta();
    m_lepPFIso2=ElectronPFIso(SelectedSB2Lep2,rho);
    m_lepDetIso2=(SelectedSB2Lep2->userIso(0)+SelectedSB2Lep2->userIso(1)+SelectedSB2Lep2->userIso(2))/SelectedSB2Lep2->pt();
    m_lepCharge2=SelectedSB2Lep2->charge();
    m_lepNumberEle=SelectedSB2Ele.size();
    m_lepNumberMuo=SelectedSB2Muo.size();
    m_charge=SelectedSB2Lep1->charge()*SelectedSB2Lep2->charge();
    m_MassVis=dilep.mass();
    m_MassEff=dilepmet.mass();
    m_MassSvfit=MassSVFit;
    m_MassCA=MassCA;
    m_XMassVis=dilepjet.mass();
    m_XMassEff=dilepmetjet.mass();
    m_XMassSVFit=XMassSVFit;
    m_XMassCA=XmassCA;
    m_met=met->begin()->pt();
    m_metPhi=met->begin()->phi();
    m_uncorrmet=uncorrmet->begin()->pt();
    m_uncorrmetPhi=uncorrmet->begin()->phi();
    m_nbtagsL=nSB2btagsL;
    m_nbtagsM=nSB2btagsM;
    m_nbtagsT=nSB2btagsT;
    m_trigger320=(int)isFired_HLT_PFJet320;
    m_trigger650=(int)isFired_HLT_HT650; 
    m_trigger=(int)isFired_HLT; 
    m_sideband=(int)(SelectedSB2PrunedMass[jetSB2Ind]<70);
    m_NVertices=vertices->size();
    m_PUWeight=MyWeight;
    m_metPx=met->begin()->px();
    m_metPy=met->begin()->py();
    m_NeventsTOT=NeventsTOT_;
    m_xsec=xsec_;
    m_lumi=lumi_;
    m_weight=xsec_*lumi_/NeventsTOT_;
    m_genEvent = genEvent;
    TreeSB2EleEle->Fill();
  }

  //DI-MUON SELECTION - SB2
  if(SelectedSB2Muo.size()>1 && hasAtLeastOneHighPtMuoSB2 == true && SelectedSB2Jet.size()>0){
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedSB2Jet[jetSB2Ind]->pt(),SelectedSB2Jet[jetSB2Ind]->eta(),SelectedSB2Jet[jetSB2Ind]->phi(),
						SelectedSB2PrunedMass[jetSB2Ind]);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());
    //SVFIT
    TLorentzVector SVFitTauTau;
    float MassSVFit = SVFitMass(met, metRaw, SVFitTauTau, SelectedSB2Muo[0]->p4(), SelectedSB2Muo[1]->p4());
    float XMassSVFit = -1;
    float dRJetZSVFit = -1;
    if(SVFitTauTau.Pt()>0){
      XMassSVFit = (SVFitTauTau+PrunedJet).M();
      dRJetZSVFit = ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,PrunedJet);
    }
    //COLLINEAR APPROXIMATION
    TLorentzVector CATauTau; bool CA = false;
    CollinearApproximation(met, CATauTau, PrunedJet, SelectedSB2Muo[0]->p4(), SelectedSB2Muo[1]->p4(), CA);
    float MassCA = -1;
    float XmassCA = -1.;
    float dRJetZCA = -1.;
    if(CA){
      MassCA = CATauTau.M();
      XmassCA = (CATauTau+PrunedJet).M();
      dRJetZCA = ROOT::Math::VectorUtil::DeltaR(CATauTau,PrunedJet);
    }
    pat::MuonCollection::const_iterator SelectedSB2Lep1;
    pat::MuonCollection::const_iterator SelectedSB2Lep2;
    float leppt=-99; unsigned int lepInd1 = -1; unsigned int lepInd2 = -1;
    for(unsigned int j=0; j<SelectedSB2Muo.size(); j++){
      if(SelectedSB2Muo[j]->pt()>leppt) {
	lepInd1 = j;
	leppt = SelectedSB2Muo[j]->pt();
      }
    }
    leppt=-99;
    for(unsigned int j=0; j<SelectedSB2Muo.size(); j++){
      if(j==lepInd1) continue;
      if(SelectedSB2Muo[j]->pt()>leppt) {
	lepInd2 = j;
	leppt = SelectedSB2Muo[j]->pt();
      }
    }
    SelectedSB2Lep1=SelectedSB2Muo[lepInd1];
    SelectedSB2Lep2=SelectedSB2Muo[lepInd2];
    math::PtEtaPhiELorentzVector dilep; dilep = SelectedSB2Lep1->p4()+SelectedSB2Lep2->p4();
    math::PtEtaPhiELorentzVector dilepmet; dilepmet = SelectedSB2Lep1->p4()+SelectedSB2Lep2->p4()+met->begin()->p4();
    math::PtEtaPhiELorentzVector dilepjet; dilepjet = SelectedSB2Lep1->p4()+SelectedSB2Lep2->p4()+PrunedJet_prov;
    math::PtEtaPhiELorentzVector dilepmetjet; dilepmetjet = SelectedSB2Lep1->p4()+SelectedSB2Lep2->p4()+met->begin()->p4()+PrunedJet_prov;
    //PLOT - TREE
    m_jetPt=SelectedSB2Jet[jetSB2Ind]->pt();
    m_jetEta=SelectedSB2Jet[jetSB2Ind]->eta();
    m_jetMass=SelectedSB2PrunedMass[jetSB2Ind];
    m_jetSubjettiness=SelectedSB2Jet[jetSB2Ind]->userFloat("tau2")/SelectedSB2Jet[jetSB2Ind]->userFloat("tau1");
    m_dRJetLep1=ROOT::Math::VectorUtil::DeltaR(SelectedSB2Lep1->p4(),SelectedSB2Jet[jetSB2Ind]->p4());
    m_dRJetLep2=ROOT::Math::VectorUtil::DeltaR(SelectedSB2Lep2->p4(),SelectedSB2Jet[jetSB2Ind]->p4());
    m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB2Jet[jetSB2Ind]->p4());
    m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB2Jet[jetSB2Ind]->p4());
    m_dRLepMet1=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB2Lep1->p4());
    m_dPhiLepMet1=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB2Lep1->p4());
    m_dRLepMet2=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB2Lep2->p4());
    m_dPhiLepMet2=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB2Lep2->p4());
    m_dRLepLep=ROOT::Math::VectorUtil::DeltaR(SelectedSB2Lep1->p4(),SelectedSB2Lep2->p4());
    m_dRZZVis=ROOT::Math::VectorUtil::DeltaR(dilep,PrunedJet);
    m_dRZZEff=ROOT::Math::VectorUtil::DeltaR(dilepmet,PrunedJet);
    m_dRZZSvFit=dRJetZSVFit;
    m_dRZZCA=dRJetZCA;
    m_lepPt1=SelectedSB2Lep1->pt();
    m_lepEta1=SelectedSB2Lep1->eta();
    m_lepPFIso1=MuonPFIso(SelectedSB2Lep1,false);
    m_lepDetIso1=(MuonDETIso(SelectedSB2Lep1, SelectedSB2Muo));
    m_lepCharge1=SelectedSB2Lep1->charge();
    m_lepPt2=SelectedSB2Lep2->pt();
    m_lepEta2=SelectedSB2Lep2->eta();
    m_lepPFIso2=MuonPFIso(SelectedSB2Lep2,false);
    m_lepDetIso2=(MuonDETIso(SelectedSB2Lep2, SelectedSB2Muo));
    m_lepCharge2=SelectedSB2Lep2->charge();
    m_lepNumberEle=SelectedSB2Ele.size();
    m_lepNumberMuo=SelectedSB2Muo.size();
    m_charge=SelectedSB2Lep1->charge()*SelectedSB2Lep2->charge();
    m_MassVis=dilep.mass();
    m_MassEff=dilepmet.mass();
    m_MassSvfit=MassSVFit;
    m_MassCA=MassCA;
    m_XMassVis=dilepjet.mass();
    m_XMassEff=dilepmetjet.mass();
    m_XMassSVFit=XMassSVFit;
    m_XMassCA=XmassCA;
    m_met=met->begin()->pt();
    m_metPhi=met->begin()->phi();
    m_uncorrmet=uncorrmet->begin()->pt();
    m_uncorrmetPhi=uncorrmet->begin()->phi();
    m_nbtagsL=nSB2btagsL;
    m_nbtagsM=nSB2btagsM;
    m_nbtagsT=nSB2btagsT;
    m_trigger320=(int)isFired_HLT_PFJet320;
    m_trigger650=(int)isFired_HLT_HT650;
    m_trigger=(int)isFired_HLT; 
    m_sideband=(int)(SelectedSB2PrunedMass[jetSB2Ind]<70);
    m_NVertices=vertices->size();
    m_PUWeight=MyWeight;
    m_metPx=met->begin()->px();
    m_metPy=met->begin()->py();
    m_NeventsTOT=NeventsTOT_;
    m_xsec=xsec_;
    m_lumi=lumi_;
    m_weight=xsec_*lumi_/NeventsTOT_;
    m_genEvent = genEvent;
    TreeSB2MuoMuo->Fill();
  }

  //ELECTRON-MUON SELECTION - SB2
  SelectedSB2Ele.clear(); SelectedSB2Muo.clear();
  for(unsigned int i=0; i<SelectedInitialEle.size(); i++){
    if(jetSB2Ind==-1) continue;
    if(ROOT::Math::VectorUtil::DeltaR(SelectedInitialEle[i]->p4(),SelectedSB2Jet[jetSB2Ind]->p4())<0.8) continue;
    //if(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedInitialEle[i]->p4())>1) continue;
    if(ElectronPFIso(SelectedInitialEle[i],rho)>0.1) continue;
    SelectedSB2Ele.push_back(SelectedInitialEle[i]);
  }
  for(unsigned int i=0; i<SelectedHighptMuo.size(); i++){
    if(jetSB2Ind==-1) continue;
    if(ROOT::Math::VectorUtil::DeltaR(SelectedHighptMuo[i]->p4(),SelectedSB2Jet[jetSB2Ind]->p4())<0.8) continue;
    //if(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedHighptMuo[i]->p4())>1) continue;
    if(MuonPFIso(SelectedHighptMuo[i], true)>0.2) continue;
    SelectedSB2Muo.push_back(SelectedHighptMuo[i]);
  }
  if(SelectedSB2Ele.size()>0 && SelectedSB2Muo.size()>0 && SelectedSB2Jet.size()>0){
    math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedSB2Jet[jetSB2Ind]->pt(),SelectedSB2Jet[jetSB2Ind]->eta(),SelectedSB2Jet[jetSB2Ind]->phi(),
						SelectedSB2PrunedMass[jetSB2Ind]);
    TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());
    //SVFIT
    TLorentzVector SVFitTauTau;
    float MassSVFit = SVFitMass(met, metRaw, SVFitTauTau, SelectedSB2Ele[0]->p4(), SelectedSB2Muo[0]->p4());
    float XMassSVFit = -1;
    float dRJetZSVFit = -1;
    if(SVFitTauTau.Pt()>0){
      XMassSVFit = (SVFitTauTau+PrunedJet).M();
      dRJetZSVFit = ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,PrunedJet);
    }
    //COLLINEAR APPROXIMATION
    TLorentzVector CATauTau; bool CA = false;
    CollinearApproximation(met, CATauTau, PrunedJet, SelectedSB2Ele[0]->p4(), SelectedSB2Muo[0]->p4(), CA);
    float MassCA = -1;
    float XmassCA = -1.;
    float dRJetZCA = -1.;
    if(CA){
      MassCA = CATauTau.M();
      XmassCA = (CATauTau+PrunedJet).M();
      dRJetZCA = ROOT::Math::VectorUtil::DeltaR(CATauTau,PrunedJet);
    }
    pat::ElectronCollection::const_iterator SelectedSB2Lep1;
    pat::MuonCollection::const_iterator SelectedSB2Lep2;
    float leppt=-99; int lepInd1 = -1; int lepInd2 = -1;
    for(unsigned int j=0; j<SelectedSB2Ele.size(); j++){
      if(SelectedSB2Ele[j]->pt()>leppt) {
	lepInd1 = j;
	leppt = SelectedSB2Ele[j]->pt();
      }
    }
    leppt=-99;
    for(unsigned int j=0; j<SelectedSB2Muo.size(); j++){
      if(SelectedSB2Muo[j]->pt()>leppt) {
	lepInd2 = j;
	leppt = SelectedSB2Muo[j]->pt();
      }
    }
    SelectedSB2Lep1=SelectedSB2Ele[lepInd1];
    SelectedSB2Lep2=SelectedSB2Muo[lepInd2];
    math::PtEtaPhiELorentzVector dilep; dilep = SelectedSB2Lep1->p4()+SelectedSB2Lep2->p4();
    math::PtEtaPhiELorentzVector dilepmet; dilepmet = SelectedSB2Lep1->p4()+SelectedSB2Lep2->p4()+met->begin()->p4();
    math::PtEtaPhiELorentzVector dilepjet; dilepjet = SelectedSB2Lep1->p4()+SelectedSB2Lep2->p4()+PrunedJet_prov;
    math::PtEtaPhiELorentzVector dilepmetjet; dilepmetjet = SelectedSB2Lep1->p4()+SelectedSB2Lep2->p4()+met->begin()->p4()+PrunedJet_prov;
    //PLOT - TREE
    m_jetPt=SelectedSB2Jet[jetSB2Ind]->pt();
    m_jetEta=SelectedSB2Jet[jetSB2Ind]->eta();
    m_jetMass=SelectedSB2PrunedMass[jetSB2Ind];
    m_jetSubjettiness=SelectedSB2Jet[jetSB2Ind]->userFloat("tau2")/SelectedSB2Jet[jetSB2Ind]->userFloat("tau1");
    m_dRJetLep1=ROOT::Math::VectorUtil::DeltaR(SelectedSB2Lep1->p4(),SelectedSB2Jet[jetSB2Ind]->p4());
    m_dRJetLep2=ROOT::Math::VectorUtil::DeltaR(SelectedSB2Lep2->p4(),SelectedSB2Jet[jetSB2Ind]->p4());
    m_dRJetMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB2Jet[jetSB2Ind]->p4());
    m_dPhiJetMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB2Jet[jetSB2Ind]->p4());
    m_dRLepMet1=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB2Lep1->p4());
    m_dPhiLepMet1=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB2Lep1->p4());
    m_dRLepMet2=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedSB2Lep2->p4());
    m_dPhiLepMet2=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedSB2Lep2->p4());
    m_dRLepLep=ROOT::Math::VectorUtil::DeltaR(SelectedSB2Lep1->p4(),SelectedSB2Lep2->p4());
    m_dRZZVis=ROOT::Math::VectorUtil::DeltaR(dilep,PrunedJet);
    m_dRZZEff=ROOT::Math::VectorUtil::DeltaR(dilepmet,PrunedJet);
    m_dRZZSvFit=dRJetZSVFit;
    m_dRZZCA=dRJetZCA;
    m_lepPt1=SelectedSB2Lep1->pt();
    m_lepEta1=SelectedSB2Lep1->eta();
    m_lepPFIso1=ElectronPFIso(SelectedSB2Lep1,rho);
    m_lepDetIso1=(SelectedSB2Lep1->userIso(0)+SelectedSB2Lep1->userIso(1)+SelectedSB2Lep1->userIso(2))/SelectedSB2Lep1->pt();
    m_lepCharge1=SelectedSB2Lep1->charge();
    m_lepPt2=SelectedSB2Lep2->pt();
    m_lepEta2=SelectedSB2Lep2->eta();
    m_lepPFIso2=MuonPFIso(SelectedSB2Lep2,false);
    m_lepDetIso2=(MuonDETIso(SelectedSB2Lep2, SelectedSB2Muo));
    m_lepCharge2=SelectedSB2Lep2->charge();
    m_lepNumberEle=SelectedSB2Ele.size();
    m_lepNumberMuo=SelectedSB2Muo.size();
    m_charge=SelectedSB2Lep1->charge()*SelectedSB2Lep2->charge();
    m_MassVis=dilep.mass();
    m_MassEff=dilepmet.mass();
    m_MassSvfit=MassSVFit;
    m_MassCA=MassCA;
    m_XMassVis=dilepjet.mass();
    m_XMassEff=dilepmetjet.mass();
    m_XMassSVFit=XMassSVFit;
    m_XMassCA=XmassCA;
    m_met=met->begin()->pt();
    m_metPhi=met->begin()->phi();
    m_uncorrmet=uncorrmet->begin()->pt();
    m_uncorrmetPhi=uncorrmet->begin()->phi();
    m_nbtagsL=nSB2btagsL;
    m_nbtagsM=nSB2btagsM;
    m_nbtagsT=nSB2btagsT;
    m_trigger320=(int)isFired_HLT_PFJet320;
    m_trigger650=(int)isFired_HLT_HT650;
    m_trigger=(int)isFired_HLT; 
    m_sideband=(int)(SelectedSB2PrunedMass[jetSB2Ind]<70);
    m_NVertices=vertices->size();
    m_PUWeight=MyWeight;
    m_metPx=met->begin()->px();
    m_metPy=met->begin()->py();
    m_NeventsTOT=NeventsTOT_;
    m_xsec=xsec_;
    m_lumi=lumi_;
    m_weight=xsec_*lumi_/NeventsTOT_;
    m_genEvent = genEvent;
    TreeSB2EleMuo->Fill();
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

  TreeSignalEff = fs->make<TTree>("TreeSignalEff", "TreeSignalEff");
  TreeSignalEff->Branch("genEvent", &m_genEvent, "genEvent/f");

  TreeEleEle = fs->make<TTree>("TreeEleEle", "TreeEleEle");
  TreeEleEle->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeEleEle->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeEleEle->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeEleEle->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeEleEle->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeEleEle->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeEleEle->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeEleEle->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeEleEle->Branch("dRLepMet1", &m_dRLepMet1, "dRLepMet1/f");
  TreeEleEle->Branch("dPhiLepMet1", &m_dPhiLepMet1, "dPhiLepMet1/f");
  TreeEleEle->Branch("dRLepMet2", &m_dRLepMet2, "dRLepMet2/f");
  TreeEleEle->Branch("dPhiLepMet2", &m_dPhiLepMet2, "dPhiLepMet2/f");
  TreeEleEle->Branch("dRLepLep", &m_dRLepLep, "dRLepLep/f");
  TreeEleEle->Branch("lepPt1", &m_lepPt1, "lepPt1/f");
  TreeEleEle->Branch("lepEta1", &m_lepEta1, "lepEta1/f");
  TreeEleEle->Branch("lepPFIso1", &m_lepPFIso1, "lepPFIso1/f");
  TreeEleEle->Branch("lepDetIso1", &m_lepDetIso1, "lepDetIso1/f");
  TreeEleEle->Branch("lepCharge1", &m_lepCharge1, "lepCharge1/f");
  TreeEleEle->Branch("lepPt2", &m_lepPt2, "lepPt2/f");
  TreeEleEle->Branch("lepEta2", &m_lepEta2, "lepEta2/f");
  TreeEleEle->Branch("lepPFIso2", &m_lepPFIso2, "lepPFIso2/f");
  TreeEleEle->Branch("lepDetIso2", &m_lepDetIso2, "lepDetIso2/f");
  TreeEleEle->Branch("lepCharge2", &m_lepCharge2, "lepCharge2/f");
  TreeEleEle->Branch("lepNumberEle", &m_lepNumberEle, "lepNumberEle/i");
  TreeEleEle->Branch("lepNumberMuo", &m_lepNumberMuo, "lepNumberMuo/i");
  TreeEleEle->Branch("charge", &m_charge, "charge/f");
  TreeEleEle->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeEleEle->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeEleEle->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeEleEle->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeEleEle->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeEleEle->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeEleEle->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeEleEle->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeEleEle->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeEleEle->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeEleEle->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeEleEle->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeEleEle->Branch("met", &m_met, "met/f");
  TreeEleEle->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeEleEle->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeEleEle->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
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

  TreeMuoMuo = fs->make<TTree>("TreeMuoMuo", "TreeMuoMuo");
  TreeMuoMuo->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeMuoMuo->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeMuoMuo->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeMuoMuo->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeMuoMuo->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeMuoMuo->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeMuoMuo->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeMuoMuo->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeMuoMuo->Branch("dRLepMet1", &m_dRLepMet1, "dRLepMet1/f");
  TreeMuoMuo->Branch("dPhiLepMet1", &m_dPhiLepMet1, "dPhiLepMet1/f");
  TreeMuoMuo->Branch("dRLepMet2", &m_dRLepMet2, "dRLepMet2/f");
  TreeMuoMuo->Branch("dPhiLepMet2", &m_dPhiLepMet2, "dPhiLepMet2/f");
  TreeMuoMuo->Branch("dRLepLep", &m_dRLepLep, "dRLepLep/f");
  TreeMuoMuo->Branch("lepPt1", &m_lepPt1, "lepPt1/f");
  TreeMuoMuo->Branch("lepEta1", &m_lepEta1, "lepEta1/f");
  TreeMuoMuo->Branch("lepPFIso1", &m_lepPFIso1, "lepPFIso1/f");
  TreeMuoMuo->Branch("lepDetIso1", &m_lepDetIso1, "lepDetIso1/f");
  TreeMuoMuo->Branch("lepCharge1", &m_lepCharge1, "lepCharge1/f");
  TreeMuoMuo->Branch("lepPt2", &m_lepPt2, "lepPt2/f");
  TreeMuoMuo->Branch("lepEta2", &m_lepEta2, "lepEta2/f");
  TreeMuoMuo->Branch("lepPFIso2", &m_lepPFIso2, "lepPFIso2/f");
  TreeMuoMuo->Branch("lepDetIso2", &m_lepDetIso2, "lepDetIso2/f");
  TreeMuoMuo->Branch("lepCharge2", &m_lepCharge2, "lepCharge2/f");
  TreeMuoMuo->Branch("lepNumberEle", &m_lepNumberEle, "lepNumberEle/i");
  TreeMuoMuo->Branch("lepNumberMuo", &m_lepNumberMuo, "lepNumberMuo/i");
  TreeMuoMuo->Branch("charge", &m_charge, "charge/f");
  TreeMuoMuo->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeMuoMuo->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeMuoMuo->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeMuoMuo->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeMuoMuo->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeMuoMuo->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeMuoMuo->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeMuoMuo->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeMuoMuo->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeMuoMuo->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeMuoMuo->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeMuoMuo->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeMuoMuo->Branch("met", &m_met, "met/f");
  TreeMuoMuo->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeMuoMuo->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeMuoMuo->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
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

  TreeEleMuo = fs->make<TTree>("TreeEleMuo", "TreeEleMuo");
  TreeEleMuo->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeEleMuo->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeEleMuo->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeEleMuo->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeEleMuo->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeEleMuo->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeEleMuo->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeEleMuo->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeEleMuo->Branch("dRLepMet1", &m_dRLepMet1, "dRLepMet1/f");
  TreeEleMuo->Branch("dPhiLepMet1", &m_dPhiLepMet1, "dPhiLepMet1/f");
  TreeEleMuo->Branch("dRLepMet2", &m_dRLepMet2, "dRLepMet2/f");
  TreeEleMuo->Branch("dPhiLepMet2", &m_dPhiLepMet2, "dPhiLepMet2/f");
  TreeEleMuo->Branch("dRLepLep", &m_dRLepLep, "dRLepLep/f");
  TreeEleMuo->Branch("lepPt1", &m_lepPt1, "lepPt1/f");
  TreeEleMuo->Branch("lepEta1", &m_lepEta1, "lepEta1/f");
  TreeEleMuo->Branch("lepPFIso1", &m_lepPFIso1, "lepPFIso1/f");
  TreeEleMuo->Branch("lepDetIso1", &m_lepDetIso1, "lepDetIso1/f");
  TreeEleMuo->Branch("lepCharge1", &m_lepCharge1, "lepCharge1/f");
  TreeEleMuo->Branch("lepPt2", &m_lepPt2, "lepPt2/f");
  TreeEleMuo->Branch("lepEta2", &m_lepEta2, "lepEta2/f");
  TreeEleMuo->Branch("lepPFIso2", &m_lepPFIso2, "lepPFIso2/f");
  TreeEleMuo->Branch("lepDetIso2", &m_lepDetIso2, "lepDetIso2/f");
  TreeEleMuo->Branch("lepCharge2", &m_lepCharge2, "lepCharge2/f");
  TreeEleMuo->Branch("lepNumberEle", &m_lepNumberEle, "lepNumberEle/i");
  TreeEleMuo->Branch("lepNumberMuo", &m_lepNumberMuo, "lepNumberMuo/i");
  TreeEleMuo->Branch("charge", &m_charge, "charge/f");
  TreeEleMuo->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeEleMuo->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeEleMuo->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeEleMuo->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeEleMuo->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeEleMuo->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeEleMuo->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeEleMuo->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeEleMuo->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeEleMuo->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeEleMuo->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeEleMuo->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeEleMuo->Branch("met", &m_met, "met/f");
  TreeEleMuo->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeEleMuo->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeEleMuo->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
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




  TreeSB1EleEle = fs->make<TTree>("TreeSB1EleEle", "TreeSB1EleEle");
  TreeSB1EleEle->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB1EleEle->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB1EleEle->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB1EleEle->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB1EleEle->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeSB1EleEle->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeSB1EleEle->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB1EleEle->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB1EleEle->Branch("dRLepMet1", &m_dRLepMet1, "dRLepMet1/f");
  TreeSB1EleEle->Branch("dPhiLepMet1", &m_dPhiLepMet1, "dPhiLepMet1/f");
  TreeSB1EleEle->Branch("dRLepMet2", &m_dRLepMet2, "dRLepMet2/f");
  TreeSB1EleEle->Branch("dPhiLepMet2", &m_dPhiLepMet2, "dPhiLepMet2/f");
  TreeSB1EleEle->Branch("dRLepLep", &m_dRLepLep, "dRLepLep/f");
  TreeSB1EleEle->Branch("lepPt1", &m_lepPt1, "lepPt1/f");
  TreeSB1EleEle->Branch("lepEta1", &m_lepEta1, "lepEta1/f");
  TreeSB1EleEle->Branch("lepPFIso1", &m_lepPFIso1, "lepPFIso1/f");
  TreeSB1EleEle->Branch("lepDetIso1", &m_lepDetIso1, "lepDetIso1/f");
  TreeSB1EleEle->Branch("lepCharge1", &m_lepCharge1, "lepCharge1/f");
  TreeSB1EleEle->Branch("lepPt2", &m_lepPt2, "lepPt2/f");
  TreeSB1EleEle->Branch("lepEta2", &m_lepEta2, "lepEta2/f");
  TreeSB1EleEle->Branch("lepPFIso2", &m_lepPFIso2, "lepPFIso2/f");
  TreeSB1EleEle->Branch("lepDetIso2", &m_lepDetIso2, "lepDetIso2/f");
  TreeSB1EleEle->Branch("lepCharge2", &m_lepCharge2, "lepCharge2/f");
  TreeSB1EleEle->Branch("lepNumberEle", &m_lepNumberEle, "lepNumberEle/i");
  TreeSB1EleEle->Branch("lepNumberMuo", &m_lepNumberMuo, "lepNumberMuo/i");
  TreeSB1EleEle->Branch("charge", &m_charge, "charge/f");
  TreeSB1EleEle->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeSB1EleEle->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeSB1EleEle->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeSB1EleEle->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeSB1EleEle->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeSB1EleEle->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeSB1EleEle->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeSB1EleEle->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeSB1EleEle->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB1EleEle->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB1EleEle->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB1EleEle->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB1EleEle->Branch("met", &m_met, "met/f");
  TreeSB1EleEle->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeSB1EleEle->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeSB1EleEle->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
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

  TreeSB1MuoMuo = fs->make<TTree>("TreeSB1MuoMuo", "TreeSB1MuoMuo");
  TreeSB1MuoMuo->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB1MuoMuo->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB1MuoMuo->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB1MuoMuo->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB1MuoMuo->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeSB1MuoMuo->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeSB1MuoMuo->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB1MuoMuo->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB1MuoMuo->Branch("dRLepMet1", &m_dRLepMet1, "dRLepMet1/f");
  TreeSB1MuoMuo->Branch("dPhiLepMet1", &m_dPhiLepMet1, "dPhiLepMet1/f");
  TreeSB1MuoMuo->Branch("dRLepMet2", &m_dRLepMet2, "dRLepMet2/f");
  TreeSB1MuoMuo->Branch("dPhiLepMet2", &m_dPhiLepMet2, "dPhiLepMet2/f");
  TreeSB1MuoMuo->Branch("dRLepLep", &m_dRLepLep, "dRLepLep/f");
  TreeSB1MuoMuo->Branch("lepPt1", &m_lepPt1, "lepPt1/f");
  TreeSB1MuoMuo->Branch("lepEta1", &m_lepEta1, "lepEta1/f");
  TreeSB1MuoMuo->Branch("lepPFIso1", &m_lepPFIso1, "lepPFIso1/f");
  TreeSB1MuoMuo->Branch("lepDetIso1", &m_lepDetIso1, "lepDetIso1/f");
  TreeSB1MuoMuo->Branch("lepCharge1", &m_lepCharge1, "lepCharge1/f");
  TreeSB1MuoMuo->Branch("lepPt2", &m_lepPt2, "lepPt2/f");
  TreeSB1MuoMuo->Branch("lepEta2", &m_lepEta2, "lepEta2/f");
  TreeSB1MuoMuo->Branch("lepPFIso2", &m_lepPFIso2, "lepPFIso2/f");
  TreeSB1MuoMuo->Branch("lepDetIso2", &m_lepDetIso2, "lepDetIso2/f");
  TreeSB1MuoMuo->Branch("lepCharge2", &m_lepCharge2, "lepCharge2/f");
  TreeSB1MuoMuo->Branch("lepNumberEle", &m_lepNumberEle, "lepNumberEle/i");
  TreeSB1MuoMuo->Branch("lepNumberMuo", &m_lepNumberMuo, "lepNumberMuo/i");
  TreeSB1MuoMuo->Branch("charge", &m_charge, "charge/f");
  TreeSB1MuoMuo->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeSB1MuoMuo->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeSB1MuoMuo->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeSB1MuoMuo->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeSB1MuoMuo->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeSB1MuoMuo->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeSB1MuoMuo->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeSB1MuoMuo->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeSB1MuoMuo->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB1MuoMuo->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB1MuoMuo->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB1MuoMuo->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB1MuoMuo->Branch("met", &m_met, "met/f");
  TreeSB1MuoMuo->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeSB1MuoMuo->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeSB1MuoMuo->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
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

  TreeSB1EleMuo = fs->make<TTree>("TreeSB1EleMuo", "TreeSB1EleMuo");
  TreeSB1EleMuo->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB1EleMuo->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB1EleMuo->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB1EleMuo->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB1EleMuo->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeSB1EleMuo->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeSB1EleMuo->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB1EleMuo->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB1EleMuo->Branch("dRLepMet1", &m_dRLepMet1, "dRLepMet1/f");
  TreeSB1EleMuo->Branch("dPhiLepMet1", &m_dPhiLepMet1, "dPhiLepMet1/f");
  TreeSB1EleMuo->Branch("dRLepMet2", &m_dRLepMet2, "dRLepMet2/f");
  TreeSB1EleMuo->Branch("dPhiLepMet2", &m_dPhiLepMet2, "dPhiLepMet2/f");
  TreeSB1EleMuo->Branch("dRLepLep", &m_dRLepLep, "dRLepLep/f");
  TreeSB1EleMuo->Branch("lepPt1", &m_lepPt1, "lepPt1/f");
  TreeSB1EleMuo->Branch("lepEta1", &m_lepEta1, "lepEta1/f");
  TreeSB1EleMuo->Branch("lepPFIso1", &m_lepPFIso1, "lepPFIso1/f");
  TreeSB1EleMuo->Branch("lepDetIso1", &m_lepDetIso1, "lepDetIso1/f");
  TreeSB1EleMuo->Branch("lepCharge1", &m_lepCharge1, "lepCharge1/f");
  TreeSB1EleMuo->Branch("lepPt2", &m_lepPt2, "lepPt2/f");
  TreeSB1EleMuo->Branch("lepEta2", &m_lepEta2, "lepEta2/f");
  TreeSB1EleMuo->Branch("lepPFIso2", &m_lepPFIso2, "lepPFIso2/f");
  TreeSB1EleMuo->Branch("lepDetIso2", &m_lepDetIso2, "lepDetIso2/f");
  TreeSB1EleMuo->Branch("lepCharge2", &m_lepCharge2, "lepCharge2/f");
  TreeSB1EleMuo->Branch("lepNumberEle", &m_lepNumberEle, "lepNumberEle/i");
  TreeSB1EleMuo->Branch("lepNumberMuo", &m_lepNumberMuo, "lepNumberMuo/i");
  TreeSB1EleMuo->Branch("charge", &m_charge, "charge/f");
  TreeSB1EleMuo->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeSB1EleMuo->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeSB1EleMuo->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeSB1EleMuo->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeSB1EleMuo->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeSB1EleMuo->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeSB1EleMuo->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeSB1EleMuo->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeSB1EleMuo->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB1EleMuo->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB1EleMuo->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB1EleMuo->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB1EleMuo->Branch("met", &m_met, "met/f");
  TreeSB1EleMuo->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeSB1EleMuo->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeSB1EleMuo->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
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




  TreeSB2EleEle = fs->make<TTree>("TreeSB2EleEle", "TreeSB2EleEle");
  TreeSB2EleEle->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB2EleEle->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB2EleEle->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB2EleEle->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB2EleEle->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeSB2EleEle->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeSB2EleEle->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB2EleEle->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB2EleEle->Branch("dRLepMet1", &m_dRLepMet1, "dRLepMet1/f");
  TreeSB2EleEle->Branch("dPhiLepMet1", &m_dPhiLepMet1, "dPhiLepMet1/f");
  TreeSB2EleEle->Branch("dRLepMet2", &m_dRLepMet2, "dRLepMet2/f");
  TreeSB2EleEle->Branch("dPhiLepMet2", &m_dPhiLepMet2, "dPhiLepMet2/f");
  TreeSB2EleEle->Branch("dRLepLep", &m_dRLepLep, "dRLepLep/f");
  TreeSB2EleEle->Branch("lepPt1", &m_lepPt1, "lepPt1/f");
  TreeSB2EleEle->Branch("lepEta1", &m_lepEta1, "lepEta1/f");
  TreeSB2EleEle->Branch("lepPFIso1", &m_lepPFIso1, "lepPFIso1/f");
  TreeSB2EleEle->Branch("lepDetIso1", &m_lepDetIso1, "lepDetIso1/f");
  TreeSB2EleEle->Branch("lepCharge1", &m_lepCharge1, "lepCharge1/f");
  TreeSB2EleEle->Branch("lepPt2", &m_lepPt2, "lepPt2/f");
  TreeSB2EleEle->Branch("lepEta2", &m_lepEta2, "lepEta2/f");
  TreeSB2EleEle->Branch("lepPFIso2", &m_lepPFIso2, "lepPFIso2/f");
  TreeSB2EleEle->Branch("lepDetIso2", &m_lepDetIso2, "lepDetIso2/f");
  TreeSB2EleEle->Branch("lepCharge2", &m_lepCharge2, "lepCharge2/f");
  TreeSB2EleEle->Branch("lepNumberEle", &m_lepNumberEle, "lepNumberEle/i");
  TreeSB2EleEle->Branch("lepNumberMuo", &m_lepNumberMuo, "lepNumberMuo/i");
  TreeSB2EleEle->Branch("charge", &m_charge, "charge/f");
  TreeSB2EleEle->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeSB2EleEle->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeSB2EleEle->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeSB2EleEle->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeSB2EleEle->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeSB2EleEle->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeSB2EleEle->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeSB2EleEle->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeSB2EleEle->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB2EleEle->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB2EleEle->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB2EleEle->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB2EleEle->Branch("met", &m_met, "met/f");
  TreeSB2EleEle->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeSB2EleEle->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeSB2EleEle->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
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

  TreeSB2MuoMuo = fs->make<TTree>("TreeSB2MuoMuo", "TreeSB2MuoMuo");
  TreeSB2MuoMuo->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB2MuoMuo->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB2MuoMuo->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB2MuoMuo->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB2MuoMuo->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeSB2MuoMuo->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeSB2MuoMuo->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB2MuoMuo->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB2MuoMuo->Branch("dRLepMet1", &m_dRLepMet1, "dRLepMet1/f");
  TreeSB2MuoMuo->Branch("dPhiLepMet1", &m_dPhiLepMet1, "dPhiLepMet1/f");
  TreeSB2MuoMuo->Branch("dRLepMet2", &m_dRLepMet2, "dRLepMet2/f");
  TreeSB2MuoMuo->Branch("dPhiLepMet2", &m_dPhiLepMet2, "dPhiLepMet2/f");
  TreeSB2MuoMuo->Branch("dRLepLep", &m_dRLepLep, "dRLepLep/f");
  TreeSB2MuoMuo->Branch("lepPt1", &m_lepPt1, "lepPt1/f");
  TreeSB2MuoMuo->Branch("lepEta1", &m_lepEta1, "lepEta1/f");
  TreeSB2MuoMuo->Branch("lepPFIso1", &m_lepPFIso1, "lepPFIso1/f");
  TreeSB2MuoMuo->Branch("lepDetIso1", &m_lepDetIso1, "lepDetIso1/f");
  TreeSB2MuoMuo->Branch("lepCharge1", &m_lepCharge1, "lepCharge1/f");
  TreeSB2MuoMuo->Branch("lepPt2", &m_lepPt2, "lepPt2/f");
  TreeSB2MuoMuo->Branch("lepEta2", &m_lepEta2, "lepEta2/f");
  TreeSB2MuoMuo->Branch("lepPFIso2", &m_lepPFIso2, "lepPFIso2/f");
  TreeSB2MuoMuo->Branch("lepDetIso2", &m_lepDetIso2, "lepDetIso2/f");
  TreeSB2MuoMuo->Branch("lepCharge2", &m_lepCharge2, "lepCharge2/f");
  TreeSB2MuoMuo->Branch("lepNumberEle", &m_lepNumberEle, "lepNumberEle/i");
  TreeSB2MuoMuo->Branch("lepNumberMuo", &m_lepNumberMuo, "lepNumberMuo/i");
  TreeSB2MuoMuo->Branch("charge", &m_charge, "charge/f");
  TreeSB2MuoMuo->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeSB2MuoMuo->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeSB2MuoMuo->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeSB2MuoMuo->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeSB2MuoMuo->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeSB2MuoMuo->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeSB2MuoMuo->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeSB2MuoMuo->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeSB2MuoMuo->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB2MuoMuo->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB2MuoMuo->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB2MuoMuo->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB2MuoMuo->Branch("met", &m_met, "met/f");
  TreeSB2MuoMuo->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeSB2MuoMuo->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeSB2MuoMuo->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
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

  TreeSB2EleMuo = fs->make<TTree>("TreeSB2EleMuo", "TreeSB2EleMuo");
  TreeSB2EleMuo->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeSB2EleMuo->Branch("jetEta", &m_jetEta, "jetEta/f");
  TreeSB2EleMuo->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeSB2EleMuo->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  TreeSB2EleMuo->Branch("dRJetLep1", &m_dRJetLep1, "dRJetLep1/f");
  TreeSB2EleMuo->Branch("dRJetLep2", &m_dRJetLep2, "dRJetLep2/f");
  TreeSB2EleMuo->Branch("dRJetMet", &m_dRJetMet, "dRJetMet/f");
  TreeSB2EleMuo->Branch("dPhiJetMet", &m_dPhiJetMet, "dPhiJetMet/f");
  TreeSB2EleMuo->Branch("dRLepMet1", &m_dRLepMet1, "dRLepMet1/f");
  TreeSB2EleMuo->Branch("dPhiLepMet1", &m_dPhiLepMet1, "dPhiLepMet1/f");
  TreeSB2EleMuo->Branch("dRLepMet2", &m_dRLepMet2, "dRLepMet2/f");
  TreeSB2EleMuo->Branch("dPhiLepMet2", &m_dPhiLepMet2, "dPhiLepMet2/f");
  TreeSB2EleMuo->Branch("dRLepLep", &m_dRLepLep, "dRLepLep/f");
  TreeSB2EleMuo->Branch("lepPt1", &m_lepPt1, "lepPt1/f");
  TreeSB2EleMuo->Branch("lepEta1", &m_lepEta1, "lepEta1/f");
  TreeSB2EleMuo->Branch("lepPFIso1", &m_lepPFIso1, "lepPFIso1/f");
  TreeSB2EleMuo->Branch("lepDetIso1", &m_lepDetIso1, "lepDetIso1/f");
  TreeSB2EleMuo->Branch("lepCharge1", &m_lepCharge1, "lepCharge1/f");
  TreeSB2EleMuo->Branch("lepPt2", &m_lepPt2, "lepPt2/f");
  TreeSB2EleMuo->Branch("lepEta2", &m_lepEta2, "lepEta2/f");
  TreeSB2EleMuo->Branch("lepPFIso2", &m_lepPFIso2, "lepPFIso2/f");
  TreeSB2EleMuo->Branch("lepDetIso2", &m_lepDetIso2, "lepDetIso2/f");
  TreeSB2EleMuo->Branch("lepCharge2", &m_lepCharge2, "lepCharge2/f");
  TreeSB2EleMuo->Branch("lepNumberEle", &m_lepNumberEle, "lepNumberEle/i");
  TreeSB2EleMuo->Branch("lepNumberMuo", &m_lepNumberMuo, "lepNumberMuo/i");
  TreeSB2EleMuo->Branch("charge", &m_charge, "charge/f");
  TreeSB2EleMuo->Branch("MassVis", &m_MassVis, "MassVis/f");
  TreeSB2EleMuo->Branch("MassEff", &m_MassEff, "MassEff/f");
  TreeSB2EleMuo->Branch("MassSvfit", &m_MassSvfit, "MassSvfit/f");
  TreeSB2EleMuo->Branch("MassCA", &m_MassCA, "MassCA/f");
  TreeSB2EleMuo->Branch("XMassVis", &m_XMassVis, "XMassVis/f");
  TreeSB2EleMuo->Branch("XMassEff", &m_XMassEff, "XMassEff/f");
  TreeSB2EleMuo->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");
  TreeSB2EleMuo->Branch("XMassCA", &m_XMassCA, "XMassCA/f");
  TreeSB2EleMuo->Branch("dRZZVis", &m_dRZZVis, "dRZZVis/f");
  TreeSB2EleMuo->Branch("dRZZEff", &m_dRZZEff, "dRZZEff/f");
  TreeSB2EleMuo->Branch("dRZZSvFit", &m_dRZZSvFit, "dRZZSvFit/f");
  TreeSB2EleMuo->Branch("dRZZCA", &m_dRZZCA, "dRZZCA/f");
  TreeSB2EleMuo->Branch("met", &m_met, "met/f");
  TreeSB2EleMuo->Branch("metPhi", &m_metPhi, "metPhi/f");
  TreeSB2EleMuo->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  TreeSB2EleMuo->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
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
				      vector<float> & SelectedPrunedMass, int & Njet,
				      float massMin, float massMax){

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
    if(!(mass>massMin && mass<massMax))  continue;
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
