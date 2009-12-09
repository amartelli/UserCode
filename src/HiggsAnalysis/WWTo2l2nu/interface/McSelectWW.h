#ifndef McSelectWW_h
#define McSelectWW_h

/*!
  class for selecting  W bosons from the generated event.
  Electrons from the W decay are also accesible.
*/

#include <CLHEP/Vector/LorentzVector.h>
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "HiggsAnalysis/WWTo2l2nu/interface/TreeContent.h"


class GenParticleCandidate;



class McSelectWW {

public:
  typedef math::XYZTLorentzVector LorentzVector;
  McSelectWW (const edm::Event& e, edm::InputTag mcTruthCollection_);


  //! "Real W" is the closer to the true W mass  
  const reco::Candidate * getRealW() { return realW_; }
  const reco::Candidate * getRealZ() { return realZ_; }   

  const reco::Candidate * getLeptonFromW() { return leptonFromW_;}
  const reco::Candidate * getNuFromW() { return nuFromW_;}
  const reco::Candidate * getLeptonFromZ() { return leptonFromZ_;}
  const reco::Candidate * getCLeptonFromZ() { return cleptonFromZ_;}

  const reco::Candidate * getLeptonFromWAftPh () { return leptonFromWAftPh_;}
  const reco::Candidate * getNuFromWAftPh () { return nuFromWAftPh_;}
  const reco::Candidate * getLeptonFromZAftPh () { return leptonFromZAftPh_;}
  const reco::Candidate * getCLeptonFromZAftPh () { return cleptonFromZAftPh_;}
  
  //! For internal consistency! Verify if there are two W bosons in the event
  bool isValid() {return validity_;}

  bool isElectrons3() { return Electrons3_; }
  bool isElectrons2Muon1() { return Electrons2Muon1_; }
  bool isMuons3() { return Muons3_; }
  bool isMuons2Electron1() { return Muons2Electron1_; }

 // return WW combined system
  LorentzVector getWZSystem() {return WZSystem_;}

  void SelectFillTree(TreeContent& myTreeVariables_x);

  void printWWEvent();
  
  const static float WMASSPDG;
  const static float ZMASSPDG;
  

private:
  
  void Init();
  void InitAfterPhotos();


  //const HepMC::GenEvent *generatorEvent_;

  bool validity_;
 

  const reco::Candidate* realW_;
  const reco::Candidate* realZ_;
  
  // electrons before PHOTOS
  const reco::Candidate* leptonFromW_;
  const reco::Candidate* nuFromW_;
  const reco::Candidate* leptonFromZ_;
  const reco::Candidate* cleptonFromZ_;
  
  // electrons after PHOTOS
  const reco::Candidate* leptonFromWAftPh_;
  const reco::Candidate* nuFromWAftPh_;
  const reco::Candidate* leptonFromZAftPh_;
  const reco::Candidate* cleptonFromZAftPh_;
  int da1_rW_pdgID_; 
  int da2_rW_pdgID_;
  int da1_rZ_pdgID_;
  int da2_rZ_pdgID_;

  bool Electrons3_;
  bool Electrons2Muon1_;
  bool Muons3_;
  bool Muons2Electron1_;


  edm::Handle<reco::GenParticleCollection> genCandidates_;
  LorentzVector WZSystem_;
  
};

#endif

