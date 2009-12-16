#ifndef WWPRESELECTIONMARKER
#define WWPRESELECTIONMARKER

#include <memory>
#include <math.h>
#include <vector>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"





#include "TH1F.h"

class WWPreselectionMarker : public edm::EDFilter {
 public:
  explicit WWPreselectionMarker(const edm::ParameterSet&);
  ~WWPreselectionMarker();
  
 private:
  void beginJob(const edm::EventSetup&) ;
  bool filter(edm::Event&, const edm::EventSetup&);
  void endJob() ;
  
  // ----------member data ---------------------------
  edm::InputTag muonslabel_;
  edm::InputTag electronslabel_;
  edm::InputTag calometlabel_;
  //  edm::InputTag jetslabel_;

  double leptonPtMinMin_;
  double leptonPtMaxMin_;
  double leptonEtaMax_;
  //  double leptonChargeCombination_;
  double metMin_;
  double invMassMin_;  
  int selectedEvents[7];

  TH1F* m_totalEvents;
  TH1F* m_passedEvents;
  TH1F* m_filterEfficiency;
};

#endif

