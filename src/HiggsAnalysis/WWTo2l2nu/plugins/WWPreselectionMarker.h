#ifndef WWPRESELECTIONMARKER
#define WWPRESELECTIONMARKER

#include <memory>
#include <math.h>
#include <vector>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"



class WWPreselectionMarker : public edm::EDProducer {
public:
  explicit WWPreselectionMarker(const edm::ParameterSet&);
  ~WWPreselectionMarker();

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  edm::InputTag muonslabel_;
  edm::InputTag electronslabel_;
  edm::InputTag calometlabel_;
  edm::InputTag jetslabel_;

  double leptonPtMinMin_;
  double leptonPtMaxMin_;
  double leptonEtaMax_;
  double leptonChargeCombination_;
  double metMin_;
  double invMassMin_;  
  int selectedEvents[7];
};

#endif

