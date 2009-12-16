#ifndef WWMUONISOLATOR
#define WWMUONISOLATOR

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"


class WWMuonIsolator{
 public:
  explicit WWMuonIsolator(const edm::ParameterSet&);
  ~WWMuonIsolator();

  // Collections to be selected
  typedef reco::MuonCollection collection;
  typedef std::vector<reco::MuonRef>::const_iterator const_iterator;

  //define iterators with above typedef
  const_iterator begin() const{return selected_.begin();}
  const_iterator end() const{return selected_.end();}


  void select(edm::Handle<reco::MuonCollection>,const edm::Event&, 
	      const edm::EventSetup&);
 private:
  
  // ----------member data ---------------------------
  
  edm::InputTag theMuonLabel;
  edm::InputTag theTrackerIsoDepositLabel;
  edm::InputTag theEcalIsoDepositLabel;
  edm::InputTag theHcalIsoDepositLabel;
  bool doRefCheck_;
  edm::InputTag selectedMuonsRefLabel_;
  std::vector<reco::MuonRef> selected_;
  double theTrackIsolCut;
  double theCaloIsolCut;

  TH1F* m_SumPt_over_Pt_MuonTk;
  TH1F* m_SumPt_MuonTk;
  TH1F* m_SumPt_MuonCalo;

};


#endif
