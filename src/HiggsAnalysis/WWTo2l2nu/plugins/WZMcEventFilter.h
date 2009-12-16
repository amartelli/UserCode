#ifndef WZMcEventFilter_h
#define WZMcEventFilter_h


#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"


class WZMcEventFilter : public edm::EDFilter {

public:
  explicit WZMcEventFilter (const edm::ParameterSet&);
  ~WZMcEventFilter();

private:
  
  bool filter (edm::Event& , const edm::EventSetup&);
  void beginJob(const edm::EventSetup&) ;
  void endJob() ;

  std::string mctype_;
  edm::InputTag mcTruthCollection_;
  //  edm::Handle<reco::GenParticleCollection> genCandidates_;

  const reco::Candidate* realW_;
  const reco::Candidate* realZ_;
  
  int da1_rW_pdgID_; 
  int da2_rW_pdgID_;
  int da1_rZ_pdgID_;
  int da2_rZ_pdgID_;

  TH1F* m_MCtotalEvents;
  TH1F* m_MCpassedEvents;
  TH1F* m_MCfilterEfficiency;
  
};

#endif

