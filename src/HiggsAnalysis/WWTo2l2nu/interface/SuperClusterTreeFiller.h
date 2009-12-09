#ifndef SuperClusterTreeFiller_h
#define SuperClusterTreeFiller_h

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "HiggsAnalysis/WWTo2l2nu/interface/TreeContent.h"
#include <TTree.h>

using namespace cms;
using namespace edm;
using namespace reco;



class SuperClusterTreeFiller {

public:

  // Constructors

  // Dump everything
  SuperClusterTreeFiller(edm::InputTag collectionTag,
			 const edm::Event&, const edm::EventSetup&,
			 int maxSC = 500);


  // Destructor
  virtual ~SuperClusterTreeFiller();

  // Operators

  // run number and all of that --- to implement

  virtual void SClusterFillTree(TreeContent& myTreeVars_);


protected:
  

  virtual void writeSCInfo(const SuperCluster *cand, TreeContent& myTreeVars_);

  //  virtual void treeSCInfo(const std::string colPrefix, const std::string colSuffix);
  
  // Friends
  edm::InputTag collectionTag_;
  const edm::Event& iEvent_;
  const edm::EventSetup& iSetup_;

  int maxSC_;
  std::string *trkIndexName_;
};

#endif // SuperClusterTreeFiller_h
