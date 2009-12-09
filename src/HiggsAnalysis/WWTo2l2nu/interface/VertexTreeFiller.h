#ifndef VertexTreeFiller_h
#define VertexTreeFiller_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "HiggsAnalysis/WWTo2l2nu/interface/TreeContent.h"
#include <TTree.h>

using namespace cms;
using namespace edm;
using namespace reco;


class VertexTreeFiller{

public:

  //! Constructors
  VertexTreeFiller(edm::InputTag vtxcollectionTag,
		   const edm::Event& iEvent, const edm::EventSetup& iSetup);

  //! Destructor
  virtual ~VertexTreeFiller();

  void VertexFillTree (TreeContent& myTreeVariables);

private:

  edm::InputTag vtxcollectionTag_;
  const edm::Event& iEvent_;
  const edm::EventSetup& iSetup_;



};

#endif // VertexTreeFiller_h
