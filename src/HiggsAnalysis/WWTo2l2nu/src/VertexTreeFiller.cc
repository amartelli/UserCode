// system include files
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
//#include "DataFormats/JetReco/interface/CaloJetCollection.h"
//#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"


#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "HiggsAnalysis/WWTo2l2nu/interface/TreeContent.h"
//#include "HiggsAnalysis/WWTo2l2nu/interface/EleIDTreeFiller.h"
#include "HiggsAnalysis/WWTo2l2nu/interface/VertexTreeFiller.h"

#include "FWCore/Framework/interface/TriggerNames.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include <TTree.h>
#include <string>

using namespace edm;
using namespace reco;





VertexTreeFiller::VertexTreeFiller(edm::InputTag vtxcollectionTag,
				   const edm::Event& iEvent, const edm::EventSetup& iSetup):
  vtxcollectionTag_(vtxcollectionTag), iEvent_(iEvent), iSetup_(iSetup)
{}


// ---------------------------------------------------------------


VertexTreeFiller::~VertexTreeFiller() {}


// ---------------------------------------------------------------
void VertexTreeFiller::VertexFillTree(TreeContent& myTreeVariables)
{
  int nVtx = 0;

  Handle<reco::VertexCollection> primaryVertex;
  try { iEvent_.getByLabel(vtxcollectionTag_, primaryVertex); }
  catch(cms::Exception& ex ) {edm::LogWarning("CmsMuonFiller") << "Can't get candidate collection: " << vtxcollectionTag_; }

  if(primaryVertex->size() > 0) {
    for(VertexCollection::const_iterator v = primaryVertex->begin();
	v != primaryVertex->end(); ++v){
      float SumPt = 0.0;
      if((*v).tracksSize() > 0){
	std::vector<TrackBaseRef >::const_iterator t;
	for( t = (*v).tracks_begin(); t != (*v).tracks_end(); t++){
	  if((**t).charge() < -1 || (**t).charge() > 1){
	    //illegal charge
	  } else {
	    SumPt += (**t).pt();
	  }
	}
	myTreeVariables.PVx->push_back((*v).x());
	myTreeVariables.PVy->push_back((*v).y());
	myTreeVariables.PVz->push_back((*v).z());
	myTreeVariables.PVErrx->push_back((*v).xError());
	myTreeVariables.PVErry->push_back((*v).yError());
	myTreeVariables.PVErrz->push_back((*v).zError());
	myTreeVariables.SumPt->push_back(SumPt);
	myTreeVariables.ndof->push_back((*v).ndof());
	myTreeVariables.chi2->push_back((*v).chi2());
	nVtx++;
      }
    }
  } else {
    myTreeVariables.PVx->push_back(-1.);
    myTreeVariables.PVy->push_back(-1.);
    myTreeVariables.PVz->push_back(-1.);
    myTreeVariables.PVErrx->push_back(-1.);
    myTreeVariables.PVErry->push_back(-1.);
    myTreeVariables.PVErrz->push_back(-1.);
    myTreeVariables.SumPt->push_back(-1.);
    myTreeVariables.ndof->push_back(-1.);
    myTreeVariables.chi2->push_back(-1.);
  }

  
}

