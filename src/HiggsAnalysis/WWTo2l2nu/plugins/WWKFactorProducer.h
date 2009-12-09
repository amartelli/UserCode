#ifndef WWKFACTORPRODUCER
#define WWKFACTORPRODUCER

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HepMC/WeightContainer.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "HiggsAnalysis/WWTo2l2nu/plugins/WWKFactorList.h"
#include "TH1D.h"
#include "TFile.h"

#include <vector>

//
// class decleration
//

class WWKFactorProducer : public edm::EDProducer {
   public:
      explicit WWKFactorProducer(const edm::ParameterSet&);
      ~WWKFactorProducer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
     
  std::string inputFilename_;
  int  processID_;
  std::vector<int> altProcessID_;
  WWKfactorList* pt_histo_;
  bool debug_;
  // use NNLO for alternative Kfactors?
  bool useNNLO_;
};

#endif

