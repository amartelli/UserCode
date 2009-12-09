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
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "HiggsAnalysis/WWTo2l2nu/interface/TreeContent.h"
//#include "HiggsAnalysis/WWTo2l2nu/interface/CandidateTreeFiller.h"
#include "HiggsAnalysis/WWTo2l2nu/interface/MuonTreeFiller.h"

#include <TTree.h>

#include <string>

using namespace edm;
using namespace reco;



MuonTreeFiller::MuonTreeFiller(edm::InputTag collectionTag,
			       const edm::Event& iEvent, const edm::EventSetup& iSetup):
  collectionTag_(collectionTag), iEvent_(iEvent), iSetup_(iSetup)
{}

//--------------
// Destructor --
//--------------

MuonTreeFiller::~MuonTreeFiller() {}


//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out



void MuonTreeFiller::writeCollectionToTree(TreeContent& myTreeVariables_){

  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent_.getByLabel(collectionTag_, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("MuonTreeFiller") << "Can't get candidate collection: " 
								   << collectionTag_; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();
  
  if(collection) {
    
    myTreeVariables_.muon_ncand = collection->size();

    edm::View<reco::Candidate>::const_iterator cand;
    for(cand = collection->begin(); cand != collection->end(); ++cand) {
      
      //Fill candidate info
      writeCandInfo(&(*cand), myTreeVariables_);
      
      // fill muon extra information
      const reco::Muon *muon = dynamic_cast< const reco::Muon *> ( &(*cand));
      writeMuonInfo(&(*muon), myTreeVariables_);

      // fill tracks extra informations (only if the muon has a tracker track)
      writeTrkInfo(&(*cand),&(*muon), myTreeVariables_);

    }
  }
  else {
    myTreeVariables_.muon_ncand = 0;
  }
 
}





//void MuonTreeFiller::writeCandInfo(const Candidate *cand, const Muon *muon, TreeContent& myTreeVariables_) {

void MuonTreeFiller::writeCandInfo(const Candidate *cand, TreeContent& myTreeVariables_) 
{

  myTreeVariables_.muon_charge -> push_back((int)cand->charge());
  myTreeVariables_.muon_energy -> push_back(cand->energy());
  myTreeVariables_.muon_pt -> push_back(cand->pt());
  myTreeVariables_.muon_et -> push_back(cand->et());
  myTreeVariables_.muon_momentum -> push_back(cand->p());
  myTreeVariables_.muon_momentumX -> push_back(cand->px());
  myTreeVariables_.muon_momentumY -> push_back(cand->py());
  myTreeVariables_.muon_momentumZ -> push_back(cand->pz());
  myTreeVariables_.muon_vertexX -> push_back(cand->vx());
  myTreeVariables_.muon_vertexY -> push_back(cand->vy());
  myTreeVariables_.muon_vertexZ -> push_back(cand->vz());
  myTreeVariables_.muon_theta -> push_back(cand->theta());
  myTreeVariables_.muon_eta -> push_back(cand->eta());
  myTreeVariables_.muon_phi -> push_back(cand->phi());
  myTreeVariables_.muon_x ->push_back(cand->momentum().x());
  myTreeVariables_.muon_y ->push_back(cand->momentum().y());
  myTreeVariables_.muon_z ->push_back(cand->momentum().z());
  myTreeVariables_.muon_mass ->push_back(cand-> mass()); 
  myTreeVariables_.muon_mt ->push_back(cand-> mt()); 
  myTreeVariables_.muon_pdgId ->push_back(cand->pdgId()); 
  myTreeVariables_.muon_nDau ->push_back(cand->numberOfDaughters());
}





void MuonTreeFiller::writeMuonInfo(const Muon *muon, TreeContent& myTreeVariables_) {
  if(&muon) {

    myTreeVariables_.muon_isGlobal->push_back(muon->isGlobalMuon());
    myTreeVariables_.muon_isTracker->push_back(muon->isTrackerMuon());
    myTreeVariables_.muon_isStandAlone->push_back(muon->isStandAloneMuon());
    myTreeVariables_.muon_isCalo->push_back(muon->isCaloMuon());

    // default isolation variables 0.3
    MuonIsolation Iso03  = muon->isolationR03();
    myTreeVariables_.muon_sumPt03->push_back(Iso03.sumPt);
    myTreeVariables_.muon_emEt03->push_back(Iso03.emEt);
    myTreeVariables_.muon_hadEt03->push_back(Iso03.hadEt);
    myTreeVariables_.muon_hoEt03->push_back(Iso03.hoEt);
    myTreeVariables_.muon_nTrk03->push_back(Iso03.nTracks);
    myTreeVariables_.muon_nJets03->push_back(Iso03.nJets);

    // default isolation variables 0.5
    MuonIsolation Iso05  = muon->isolationR05();
    myTreeVariables_.muon_sumPt05->push_back(Iso05.sumPt);
    myTreeVariables_.muon_emEt05->push_back(Iso05.emEt);
    myTreeVariables_.muon_hadEt05->push_back(Iso05.hadEt);
    myTreeVariables_.muon_hoEt05->push_back(Iso05.hoEt);
    myTreeVariables_.muon_nTrk05->push_back(Iso05.nTracks);
    myTreeVariables_.muon_nJets05->push_back(Iso05.nJets);

    // Expected deposits in CALO
    myTreeVariables_.muon_EcalExpDepo->push_back(muon->calEnergy().em);
    myTreeVariables_.muon_HcalExpDepo->push_back(muon->calEnergy().had);
    myTreeVariables_.muon_HoExpDepo->push_back(muon->calEnergy().ho);
    myTreeVariables_.muon_emS9->push_back(muon->calEnergy().emS9);
    myTreeVariables_.muon_hadS9->push_back(muon->calEnergy().hadS9);
    myTreeVariables_.muon_hoS9->push_back(muon->calEnergy().hoS9);
    myTreeVariables_.muon_CaloComp->push_back(muon->caloCompatibility());

  } else {
    
    // default isolation variables 0.3
    myTreeVariables_.muon_sumPt03->push_back(-1.);
    myTreeVariables_.muon_emEt03->push_back(-1.);
    myTreeVariables_.muon_hadEt03->push_back(-1.);
    myTreeVariables_.muon_hoEt03->push_back(-1.);
    myTreeVariables_.muon_nTrk03->push_back(-1.);
    myTreeVariables_.muon_nJets03->push_back(-1.);
    // default isolation variables 0.5
     myTreeVariables_.muon_sumPt05->push_back(-1.);
    myTreeVariables_.muon_emEt05->push_back(-1.);
    myTreeVariables_.muon_hadEt05->push_back(-1.);
    myTreeVariables_.muon_hoEt05->push_back(-1.);
    myTreeVariables_.muon_nTrk05->push_back(-1.);
    myTreeVariables_.muon_nJets05->push_back(-1.);

    // Expected deposits in CALO
    myTreeVariables_.muon_EcalExpDepo->push_back(-1.);
    myTreeVariables_.muon_HcalExpDepo->push_back(-1.);
    myTreeVariables_.muon_HoExpDepo->push_back(-1.);
    myTreeVariables_.muon_emS9->push_back(-1.);
    myTreeVariables_.muon_hadS9->push_back(-1.);
    myTreeVariables_.muon_hoS9->push_back(-1.);
    myTreeVariables_.muon_CaloComp->push_back(-1.);
  }
}


 
void MuonTreeFiller::writeTrkInfo(const Candidate *cand, const Muon *muon, TreeContent& myTreeVariables_) {

  
  TrackRef trkRef;
  bool hasTrackerTrack = false;

  if( & muon ) {  
    if ( muon->track().isNonnull() ) {
      hasTrackerTrack = true;
      trkRef = cand->get<TrackRef>();
    }
  }

  if( hasTrackerTrack && &trkRef!=0 ) {
    
    myTreeVariables_.muon_pxAtInner->push_back(trkRef->innerMomentum().x());
    myTreeVariables_.muon_pyAtInner->push_back(trkRef->innerMomentum().y());
    myTreeVariables_.muon_pzAtInner->push_back(trkRef->innerMomentum().z());
    
    myTreeVariables_.muon_xAtInner->push_back(trkRef->innerPosition().x());
    myTreeVariables_.muon_yAtInner->push_back(trkRef->innerPosition().y());
    myTreeVariables_.muon_zAtInner->push_back(trkRef->innerPosition().z());
    
    myTreeVariables_.muon_pxAtOuter->push_back(trkRef->outerMomentum().x());
    myTreeVariables_.muon_pyAtOuter->push_back(trkRef->outerMomentum().y());
    myTreeVariables_.muon_pzAtOuter->push_back(trkRef->outerMomentum().z());
      
    myTreeVariables_.muon_xAtOuter->push_back(trkRef->outerPosition().x());
    myTreeVariables_.muon_yAtOuter->push_back(trkRef->outerPosition().y());
    myTreeVariables_.muon_zAtOuter->push_back(trkRef->outerPosition().z());

    myTreeVariables_.muon_TrackVx ->push_back(trkRef->vx());
    myTreeVariables_.muon_TrackVy ->push_back(trkRef->vy());
    myTreeVariables_.muon_TrackVz ->push_back(trkRef->vz());
  }

     /*
    else {
      myTreeVariables_.muon_pxAtInner->push_back( -1.0 );
      myTreeVariables_.muon_pyAtInner->push_back( -1.0 );
      myTreeVariables_.muon_pzAtInner->push_back( -1.0 );
      
      myTreeVariables_.muon_xAtInner->push_back( -1.0 );
      myTreeVariables_.muon_yAtInner->push_back( -1.0 );
      myTreeVariables_.muon_zAtInner->push_back( -1.0 );
      
      myTreeVariables_.muon_pxAtOuter->push_back( -1.0 );
      myTreeVariables_.muon_pyAtOuter->push_back( -1.0 );
      myTreeVariables_.muon_pzAtOuter->push_back( -1.0 );
      
      myTreeVariables_.muon_xAtOuter->push_back( -1.0 );
      myTreeVariables_.muon_yAtOuter->push_back( -1.0 );
      myTreeVariables_.muon_zAtOuter->push_back( -1.0 );
  
    myTreeVariables_.muon_pxAtInner->push_back(-1.);
    myTreeVariables_.muon_pyAtInner->push_back(-1.);
    myTreeVariables_.muon_pzAtInner->push_back(-1.);

    myTreeVariables_.muon_xAtInner->push_back(-1.);
    myTreeVariables_.muon_yAtInner->push_back(-1.);
    myTreeVariables_.muon_zAtInner->push_back(-1.);

    myTreeVariables_.muon_pxAtOuter->push_back(-1.);
    myTreeVariables_.muon_pyAtOuter->push_back(-1.);
    myTreeVariables_.muon_pzAtOuter->push_back(-1.);

    myTreeVariables_.muon_xAtOuter->push_back(-1.);
    myTreeVariables_.muon_yAtOuter->push_back(-1.);
    myTreeVariables_.muon_zAtOuter->push_back(-1.);

    myTreeVariables_.muon_TrackNormalizedChi2->push_back(-1.);

    myTreeVariables_.muon_TrackDxy->push_back(-1.);
    myTreeVariables_.muon_TrackD0 ->push_back(-1.);
    myTreeVariables_.muon_TrackDsz->push_back(-1.);
    myTreeVariables_.muon_TrackDz ->push_back(-1.);

    myTreeVariables_.muon_TrackDxyError->push_back(-1.);
    myTreeVariables_.muon_TrackD0Error ->push_back(-1.);
    myTreeVariables_.muon_TrackDszError->push_back(-1.);
    myTreeVariables_.muon_TrackDzError ->push_back(-1.);

    myTreeVariables_.muon_TrackValidHits->push_back(-1.);					       
    myTreeVariables_.muon_TrackLostHits ->push_back(-1.);

    myTreeVariables_.muon_TrackVx ->push_back(-1.);
    myTreeVariables_.muon_TrackVy ->push_back(-1.);
    myTreeVariables_.muon_TrackVz ->push_back(-1.);
  }

     */

}

