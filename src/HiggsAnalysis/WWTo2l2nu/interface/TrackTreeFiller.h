#ifndef TrackTreeFiller_h
#define TrackTreeFiller_h

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "HiggsAnalysis/WWTo2l2nu/interface/TreeContent.h"

//#include "HiggsAnalysis/WWTo2l2nu/interface/CandidateTreeFiller.h"
#include "DataFormats/Math/interface/Point3D.h"
#include <TTree.h>

using namespace cms;
using namespace edm;
using namespace reco;




class TrackTreeFiller {

 public:

  // Constructors

  // Dump everything
  TrackTreeFiller(edm::InputTag vertexCollection,
		  edm::InputTag collectionTag,
		  const edm::Event&, const edm::EventSetup& );


  // Destructor
  virtual ~TrackTreeFiller();



  /// Find Primary Vertex
  void findPrimaryVertex();

  // Operators


  void TrackFillTree(TreeContent& myTreeVariables);

  

 private:

  void writeCandInfo(const Candidate *cand, TreeContent& myTreeVariables_);  
  void writeTrkInfo(const Candidate *cand, TrackRef trkRef, TreeContent& myTreeVariables_);


  edm::InputTag vertexCollection_;
  edm::InputTag collectionTag_;

  const edm::Event& iEvent_;
  const edm::EventSetup& iSetup_;

  // Primary Vertex in point format
  float x0, y0, z0;
};

#endif // TrackTreeFiller_h
