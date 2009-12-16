#ifndef WWMUONSELECTOR
#define WWMUONSELECTOR

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

class WWMuonSelector {
 public:
  explicit WWMuonSelector(const edm::ParameterSet&);
  ~WWMuonSelector();

  typedef reco::MuonCollection collection;
  typedef std::vector<reco::MuonRef>::const_iterator const_iterator;

  const_iterator begin() const {return selected_.begin();}
  const_iterator end() const {return selected_.end();}
  
  void select(edm::Handle<reco::MuonCollection>,const edm::Event&, 
	      const edm::EventSetup&);

 private:

  // ----------member data ---------------------------
 
  std::vector<reco::MuonRef> selected_;
  edm::InputTag muonLabel_;
  float muonPtMinBarrel_;
  float muonPtMinEndcap_;
  float muonPMinEndcap_;
  float muonEtaMax_;
};

#endif
   

