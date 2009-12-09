// system include files
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"

//#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "HiggsAnalysis/WWTo2l2nu/interface/TreeContent.h"
//#include "HiggsAnalysis/WWTo2l2nu/interface/CandidateTreeFiller.h"
#include "HiggsAnalysis/WWTo2l2nu/interface/TrackTreeFiller.h"
#include "DataFormats/Math/interface/Point3D.h"

#include <TTree.h>

#include <string>

using namespace edm;
using namespace reco;




TrackTreeFiller::TrackTreeFiller(edm::InputTag vertexCollection,
				 edm::InputTag collectionTag,
				 const edm::Event& iEvent, const edm::EventSetup& iSetup):
  vertexCollection_(vertexCollection), collectionTag_(collectionTag), iEvent_(iEvent), iSetup_(iSetup)
{
  x0 = 0.;
  y0 = 0.;
  z0 = 0.;
}


TrackTreeFiller::~TrackTreeFiller() {}

//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out


void TrackTreeFiller::findPrimaryVertex(){

  edm::Handle< reco::VertexCollection>  primaryVertex  ;
  try { iEvent_.getByLabel(vertexCollection_, primaryVertex); }
  catch ( cms::Exception& ex ) { edm::LogWarning("TrackTreeFiller") << "Can't get candidate collection: " 
								    << vertexCollection_; }
  
  if(primaryVertex->size() < 1) { // there is no vertex in the event
    x0 = 0.;
    y0 = 0.;
    z0 = 0.;
  } 
  else {
    float MaxSumPt = -1.;
    VertexCollection::const_iterator vMax = primaryVertex->begin();

    // calculate the vertex pT 
    for(VertexCollection::const_iterator v = primaryVertex->begin();
	v != primaryVertex->end(); ++v){
      float SumPt = 0.0;
      if((*v).tracksSize() > 0){
	std::vector<TrackBaseRef >::const_iterator t;
	for( t = (*v).tracks_begin(); t != (*v).tracks_end(); t++){
	  if((**t).charge() < -1 || (**t).charge() > 1){ /*illegal charge */ } 
	  else { SumPt += (**t).pt(); }
	}
      }
      
      if(SumPt > MaxSumPt) {
	MaxSumPt = SumPt;
	vMax  = v;
      } 
    }

    x0 = vMax->x();
    y0 = vMax->y();
    z0 = vMax->z();
  }
}


void TrackTreeFiller::TrackFillTree(TreeContent& myTreeVariables_){

  findPrimaryVertex();

  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent_.getByLabel(collectionTag_, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("TrackTreeFiller") << "Can't get candidate collection: " 
								    << collectionTag_; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();

  if(collection) {

    myTreeVariables_.track_ncand = collection->size();

    edm::View<reco::Candidate>::const_iterator cand;
    for(cand=collection->begin(); cand!=collection->end(); ++cand) {
      // fill basic kinematics

      writeCandInfo(&(*cand), myTreeVariables_);
      
      // fill tracks extra informations
      TrackRef trkRef = cand->get<TrackRef>();
      writeTrkInfo(&(*cand), trkRef, myTreeVariables_);
    }
  }
  else {
    myTreeVariables_.track_ncand = collection->size();
  }
  
}


void TrackTreeFiller::writeCandInfo(const Candidate *cand, TreeContent& myTreeVariables_)
{

  myTreeVariables_.track_charge -> push_back((int)cand->charge());
  myTreeVariables_.track_energy -> push_back(cand->energy());
  myTreeVariables_.track_pt -> push_back(cand->pt());
  myTreeVariables_.track_et -> push_back(cand->et());
  myTreeVariables_.track_momentum -> push_back(cand->p());
  myTreeVariables_.track_momentumX -> push_back(cand->px());
  myTreeVariables_.track_momentumY -> push_back(cand->py());
  myTreeVariables_.track_momentumZ -> push_back(cand->pz());
  myTreeVariables_.track_vertexX -> push_back(cand->vx());
  myTreeVariables_.track_vertexY -> push_back(cand->vy());
  myTreeVariables_.track_vertexZ -> push_back(cand->vz());
  myTreeVariables_.track_theta -> push_back(cand->theta());
  myTreeVariables_.track_eta -> push_back(cand->eta());
  myTreeVariables_.track_phi -> push_back(cand->phi());
  myTreeVariables_.track_x ->push_back(cand->momentum().x());
  myTreeVariables_.track_y ->push_back(cand->momentum().y());
  myTreeVariables_.track_z ->push_back(cand->momentum().z());
  myTreeVariables_.track_mass ->push_back(cand-> mass());
  myTreeVariables_.track_mt ->push_back(cand-> mt());
  myTreeVariables_.track_pdgId ->push_back(cand->pdgId());
  myTreeVariables_.track_nDau ->push_back(cand->numberOfDaughters());
}



void TrackTreeFiller::writeTrkInfo(const Candidate *cand, TrackRef trkRef,
				   TreeContent& myTreeVariables_) {
  
  if(&trkRef) {
  
    
    // Find the vertex the track belongs to
    Handle<reco::VertexCollection> primaryVertex;
    try { iEvent_.getByLabel(vertexCollection_, primaryVertex); }
    catch ( cms::Exception& ex ) { edm::LogWarning("TrackTreeFiller") << "Can't get candidate collection: " 
								      << vertexCollection_; }
    
    int iVtx = -1;
    int counter = 0;
    double weight = 0.;
    // if(saveVtxTrk_) {
      if(primaryVertex->size() >0 ) { // there is at least one vertex in the event
	for(VertexCollection::const_iterator v = primaryVertex->begin();
	    v != primaryVertex->end(); ++v){
	  double tmpw = v->trackWeight(trkRef);
	  if(tmpw > weight) {
	    if(weight > 0) edm::LogWarning("TrackTreeFiller") << "I found this track in two vertices!!!!!!" ;
	    weight = tmpw;
	    iVtx = counter;
	  }
	  counter++;
	}
      }
      //  }
    
    // vertex information
    myTreeVariables_.track_vtxIndex->push_back(iVtx);
    myTreeVariables_.track_vtxWeight->push_back(weight);
  

    // Inner Tracker information
    myTreeVariables_.track_pxAtInner->push_back(trkRef->innerMomentum().x());
    myTreeVariables_.track_pyAtInner->push_back(trkRef->innerMomentum().y());
    myTreeVariables_.track_pzAtInner->push_back(trkRef->innerMomentum().z());
      
    myTreeVariables_.track_xAtInner->push_back(trkRef->innerPosition().x());
    myTreeVariables_.track_yAtInner->push_back(trkRef->innerPosition().y());
    myTreeVariables_.track_zAtInner->push_back(trkRef->innerPosition().z());
      
    // Outer Tracker information
    myTreeVariables_.track_pxAtOuter->push_back(trkRef->outerMomentum().x());
    myTreeVariables_.track_pyAtOuter->push_back(trkRef->outerMomentum().y());
    myTreeVariables_.track_pzAtOuter->push_back(trkRef->outerMomentum().z());
      
    myTreeVariables_.track_xAtOuter->push_back(trkRef->outerPosition().x());
    myTreeVariables_.track_yAtOuter->push_back(trkRef->outerPosition().y());
    myTreeVariables_.track_zAtOuter->push_back(trkRef->outerPosition().z());

   
    // track quality
    myTreeVariables_.track_ValidHits->push_back(trkRef->numberOfValidHits());
    myTreeVariables_.track_LostHits ->push_back(trkRef->numberOfLostHits());
    myTreeVariables_.track_NormalizedChi2->push_back(trkRef->normalizedChi2());
    myTreeVariables_.track_recHitsSize->push_back(trkRef->recHitsSize());

    /// vtx position
    myTreeVariables_.track_Vx ->push_back(trkRef->vx());
    myTreeVariables_.track_Vy ->push_back(trkRef->vy());
    myTreeVariables_.track_Vz ->push_back(trkRef->vz());

    // distance w.r.t. (0,0,0)
    myTreeVariables_.track_Dxy->push_back(trkRef->dxy());
    myTreeVariables_.track_D0 ->push_back(trkRef->d0());
    myTreeVariables_.track_Dsz->push_back(trkRef->dsz());
    myTreeVariables_.track_Dz->push_back(trkRef->dz());

    myTreeVariables_.track_DxyError->push_back(trkRef->dxyError());
    myTreeVariables_.track_D0Error ->push_back(trkRef->d0Error());
    myTreeVariables_.track_DszError->push_back(trkRef->dszError());
    myTreeVariables_.track_DzError ->push_back(trkRef->dzError());

    // distance w.r.t. primary vertex
    myTreeVariables_.track_DxyPV->push_back(trkRef->dxy(math::XYZPoint(x0,y0,z0)));
    myTreeVariables_.track_DszPV->push_back(trkRef->dsz(math::XYZPoint(x0,y0,z0)));
    myTreeVariables_.track_DzPV->push_back(trkRef->dz(math::XYZPoint(x0,y0,z0)));
     
  }
}

