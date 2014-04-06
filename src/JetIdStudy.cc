// -*- C++ -*-
//
// Package:    JetIdStudy
// Class:      JetIdStudy
// 
/**\class JetIdStudy JetIdStudy.cc ExoDiBosonResonances/EDBRTauAnalyzer/src/JetIdStudy.cc

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

//new inclusion
#include "Math/VectorUtil.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TauAnalysis/CandidateTools/interface/NSVfitStandaloneAlgorithm.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//
// class declaration
//

class JetIdStudy : public edm::EDAnalyzer {
public:
  explicit JetIdStudy(const edm::ParameterSet&);
  ~JetIdStudy();
  
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
  float m_jetPt;
  float m_jetMass;
  float m_jetEta;
  float m_jetSubjettiness;
  float m_jetMuonEF;
  float m_jetPhotonEF;
  float m_jetChargedEmEF;
  float m_jetNeutralHEF;
  float m_jetChargedHEF;
  int m_cutPt;
  int m_cutMass;
  int m_cutEta;
  int m_cutSubjettiness;
  int m_cutMuonEF;
  int m_cutPhotonEF;
  int m_cutChargedEmEF;
  int m_cutNeutralHEF;
  int m_cutChargedHEF;

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
JetIdStudy::JetIdStudy(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
}


JetIdStudy::~JetIdStudy()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
JetIdStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel("primaryVertexFilter", vertices);
  reco::Vertex primaryVertex;
  primaryVertex = vertices->at(0);

  edm::Handle<pat::JetCollection> CA8JetswithQjets;
  iEvent.getByLabel("selectedPatJetsCA8CHSwithQJetsForBoostedTaus", CA8JetswithQjets);
  edm::Handle<pat::JetCollection> CA8JetsPruned;
  iEvent.getByLabel("selectedPatJetsCA8CHSprunedForBoostedTaus", CA8JetsPruned);

  Handle<vector<reco::GenParticle> > genParts;
  iEvent.getByLabel("genParticles", genParts);

  vector<math::PtEtaPhiELorentzVector> genquark;
  for(size_t ngenPart=0; ngenPart<genParts->size(); ngenPart++){
    const reco::GenParticle & genPart = (*genParts)[ngenPart];
    if(abs(genPart.pdgId())>0 && abs(genPart.pdgId())<7){
      const reco::Candidate * mom = genPart.mother();
      if(mom->pdgId()==23){
	math::PtEtaPhiELorentzVector gen_prov; gen_prov=genPart.p4();
	genquark.push_back(gen_prov);
      }
    }
  }

  bool matching=false;
  float dRmatch = 9999.; pat::JetCollection::const_iterator SelectedJet;
  for(pat::JetCollection::const_iterator jet = CA8JetswithQjets->begin(); jet != CA8JetswithQjets->end(); ++jet) {
    for(unsigned int i=0; i<genquark.size(); i++){
      if(ROOT::Math::VectorUtil::DeltaR(jet->p4(),genquark[i])<dRmatch && ROOT::Math::VectorUtil::DeltaR(jet->p4(),genquark[i])<0.8){
	dRmatch = ROOT::Math::VectorUtil::DeltaR(jet->p4(),genquark[i]);
	SelectedJet = jet;
	matching=true;
      }
    }
  }

  if(matching){
    float dRmin = 9999.; float mass = 0.;
    for(pat::JetCollection::const_iterator jetPruned = CA8JetsPruned->begin(); jetPruned != CA8JetsPruned->end(); ++jetPruned) {
      float dRtmp = ROOT::Math::VectorUtil::DeltaR(SelectedJet->p4(),jetPruned->p4());
      if(dRtmp<dRmin && dRtmp<0.7 ){//matching failed if greater than jet radius
	dRmin=dRtmp;
	mass=jetPruned->mass();
      }
    }
    
    m_jetPt = SelectedJet->pt();
    m_jetMass = mass;
    m_jetEta = SelectedJet->eta();
    m_jetSubjettiness = SelectedJet->userFloat("tau2")/SelectedJet->userFloat("tau1");
    m_jetMuonEF = SelectedJet->muonEnergyFraction();
    m_jetPhotonEF = SelectedJet->photonEnergyFraction();
    m_jetChargedEmEF = SelectedJet->chargedEmEnergyFraction();
    m_jetNeutralHEF = SelectedJet->neutralHadronEnergyFraction();
    m_jetChargedHEF = SelectedJet->chargedHadronEnergyFraction();
    m_cutPt = (int)(SelectedJet->pt()>350);
    m_cutMass = (int)(mass>70 && mass<110);
    m_cutEta = (int)(abs(SelectedJet->eta())<2.4);
    m_cutSubjettiness = (int)((SelectedJet->userFloat("tau2")/SelectedJet->userFloat("tau1"))<0.75);
    m_cutMuonEF = (int)(SelectedJet->muonEnergyFraction()<0.99);
    m_cutPhotonEF = (int)(SelectedJet->photonEnergyFraction()<0.99);
    m_cutChargedEmEF = (int)(SelectedJet->chargedEmEnergyFraction()<0.99);
    m_cutNeutralHEF = (int)(SelectedJet->neutralHadronEnergyFraction()<0.99);
    m_cutChargedHEF = (int)(SelectedJet->chargedHadronEnergyFraction()>0.0);
    Tree->Fill();
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
JetIdStudy::beginJob()
{
  Service<TFileService> fs;

  Tree = fs->make<TTree>("Tree", "Tree");
  Tree->Branch("jetPt", &m_jetPt, "jetPt/f");
  Tree->Branch("jetMass", &m_jetMass, "jetMass/f");
  Tree->Branch("jetEta", &m_jetEta, "jetEta/f");
  Tree->Branch("jetSubjettiness", &m_jetSubjettiness, "jetSubjettiness/f");
  Tree->Branch("jetMuonEF", &m_jetMuonEF, "jetMuonEF/f");
  Tree->Branch("jetPhotonEF", &m_jetPhotonEF, "jetPhotonEF/f");
  Tree->Branch("jetPhotonEF", &m_jetPhotonEF, "jetPhotonEF/f");
  Tree->Branch("jetNeutralHEF", &m_jetNeutralHEF, "jetNeutralHEF/f");
  Tree->Branch("jetChargedHEF", &m_jetChargedHEF, "jetChargedHEF/f");
  Tree->Branch("cutPt", &m_cutPt, "cutPt/i");
  Tree->Branch("cutMass", &m_cutMass, "cutMass/i");
  Tree->Branch("cutEta", &m_cutEta, "cutEta/i");
  Tree->Branch("cutSubjettiness", &m_cutSubjettiness, "cutSubjettiness/i");
  Tree->Branch("cutMuonEF", &m_cutMuonEF, "cutMuonEF/i");
  Tree->Branch("cutPhotonEF", &m_cutPhotonEF, "cutPhotonEF/i");
  Tree->Branch("cutPhotonEF", &m_cutPhotonEF, "cutPhotonEF/i");
  Tree->Branch("cutNeutralHEF", &m_cutNeutralHEF, "cutNeutralHEF/i");
  Tree->Branch("cutChargedHEF", &m_cutChargedHEF, "cutChargedHEF/i");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetIdStudy::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
JetIdStudy::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
JetIdStudy::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
JetIdStudy::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
JetIdStudy::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetIdStudy::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(JetIdStudy);
