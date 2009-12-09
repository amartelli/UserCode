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

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "HiggsAnalysis/WWTo2l2nu/interface/TreeContent.h"
#include "HiggsAnalysis/WWTo2l2nu/interface/SuperClusterTreeFiller.h"

#include <TTree.h>

#include <string>

using namespace edm;
using namespace reco;



SuperClusterTreeFiller::SuperClusterTreeFiller(edm::InputTag collectionTag,
					       const edm::Event& iEvent, const edm::EventSetup& iSetup,
					       int maxSC):
  collectionTag_(collectionTag), iEvent_(iEvent), iSetup_(iSetup)
{
  trkIndexName_ = new std::string("n");
  maxSC_ = maxSC;
}


SuperClusterTreeFiller::~SuperClusterTreeFiller() 
{}


// --------------------------------------------------------------------- //


void SuperClusterTreeFiller::SClusterFillTree(TreeContent& myTreeVars_)
{
  Handle<SuperClusterCollection> collectionHandle;
  try { iEvent_.getByLabel(collectionTag_, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("SuperClusterTreeFiller") << "Can't get SC Collection: " << collectionTag_; }
  const SuperClusterCollection *collection = collectionHandle.product();

  if(collection) 
    {
      if((int)collection->size() > maxSC_)
	{
	  edm::LogError("SuperClusterTreeFiller") << "Track length " << collection->size() 
						  << " is too long for declared max length for tree "
						  << maxSC_ 
						  << ". Collection will be truncated ";
	}
      
      myTreeVars_.nSC = collection->size();
      
      SuperClusterCollection::const_iterator cand;
  
      for(cand=collection->begin(); cand!=collection->end(); cand++) 
	{
	  // fill basic kinematics
	  writeSCInfo(&(*cand), myTreeVars_);
	}
    }
  else 
    {
      myTreeVars_.nSC = collection->size();
    }
  
  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the 
  // tree 
  
}






void SuperClusterTreeFiller::writeSCInfo(const SuperCluster *cand, 
					 TreeContent& myTreeVars_) 
{
  myTreeVars_.nBC->push_back((int)cand->clustersSize());
  
  int ncry = 0;
 
  for(reco::CaloCluster_iterator bcItr = cand->clustersBegin(); bcItr!=cand->clustersEnd(); ++bcItr){         
    ncry += (*bcItr)->size();                                                                       
  }  
  myTreeVars_.nCrystals->push_back(ncry);
  
  myTreeVars_.rawEnergy->push_back((float)cand->rawEnergy());
  myTreeVars_.energy->push_back((float)cand->energy());
  myTreeVars_.eta->push_back((float)cand->position().eta());
  myTreeVars_.phi->push_back((float)cand->position().phi());
}


