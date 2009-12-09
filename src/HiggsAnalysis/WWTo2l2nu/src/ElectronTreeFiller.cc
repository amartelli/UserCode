// system include files
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//#include "DataFormats/METReco/interface/CaloMET.h"
//#include "DataFormats/METReco/interface/CaloMETCollection.h"
//#include "DataFormats/METReco/interface/GenMET.h"
//#include "DataFormats/METReco/interface/GenMETCollection.h"

//#include "DataFormats/JetReco/interface/CaloJet.h"
//#include "DataFormats/JetReco/interface/CaloJetCollection.h"
//#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloID.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/DetId/interface/DetId.h"
//#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "HiggsAnalysis/WWTo2l2nu/interface/TreeContent.h"
//#include "HiggsAnalysis/WWTo2l2nu/interface/EleIDTreeFiller.h"
//#include "HiggsAnalysis/WWTo2l2nu/interface/CandidateTreeFiller.h"
#include "HiggsAnalysis/WWTo2l2nu/interface/ElectronTreeFiller.h"

#include <TTree.h>

#include <string>

using namespace edm;
using namespace reco;



ElectronTreeFiller::ElectronTreeFiller(edm::InputTag recoEleCollectionTag,
                                       const edm::Event& iEvent, const edm::EventSetup& iSetup):
  recoEleCollectionTag_(recoEleCollectionTag), iEvent_(iEvent), iSetup_(iSetup)
{}


ElectronTreeFiller::~ElectronTreeFiller() {}


// -------------------------------------------------------- //

void ElectronTreeFiller::writeCollectionToTree(TreeContent& myTreeVariables_){
  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent_.getByLabel(recoEleCollectionTag_, collectionHandle); }
  catch ( cms::Exception& ex )
    { edm::LogWarning("ElectronTreeFiller") << "Can't get electron candidate collection: " << recoEleCollectionTag_; }
  const edm::View<reco::Candidate>* collection = collectionHandle.product();
  
  if(collection) {
    myTreeVariables_.ele_ncand = collection->size();
    
    edm::View<reco::Candidate>::const_iterator cand;
    for(cand=collection->begin(); cand!=collection->end(); ++cand) {

      writeCandInfo(&(*cand), myTreeVariables_);

      // fill Cluster + EleID variables + Isolation
      const SuperClusterRef sclusRef = cand->get<SuperClusterRef>();
      // const GsfElectronRef electronRef = cand->get<ElectronRef>();
      const reco::GsfElectron* electronRef = dynamic_cast< const reco::GsfElectron *> ( &(*cand));
      
      writeEcalInfo(&(*cand), sclusRef, (electronRef), myTreeVariables_ );

      // fill (GSF) Track Adapter
      const GsfTrackRef trkRef = cand->get<GsfTrackRef>();
      writeTrkInfo(&(*cand), electronRef, trkRef, myTreeVariables_);
      

      //      eleIDFiller.writeCollectionToTree(myTreeVariables_);

    }
    
  } // end if(collection)
  else {
    myTreeVariables_.ele_ncand = 0;
  }
  
}

void ElectronTreeFiller::writeCandInfo(const Candidate *cand, TreeContent& myTreeVariables_)
{
  
  myTreeVariables_.ele_charge -> push_back((int)cand->charge());
  myTreeVariables_.ele_energy -> push_back(cand->energy());
  myTreeVariables_.ele_pt -> push_back(cand->pt());
  myTreeVariables_.ele_et -> push_back(cand->et());
  myTreeVariables_.ele_momentum -> push_back(cand->p());
  myTreeVariables_.ele_momentumX -> push_back(cand->px());
  myTreeVariables_.ele_momentumY -> push_back(cand->py());
  myTreeVariables_.ele_momentumZ -> push_back(cand->pz());
  myTreeVariables_.ele_vertexX -> push_back(cand->vx());
  myTreeVariables_.ele_vertexY -> push_back(cand->vy());
  myTreeVariables_.ele_vertexZ -> push_back(cand->vz());
  myTreeVariables_.ele_theta -> push_back(cand->theta());
  myTreeVariables_.ele_eta -> push_back(cand->eta());
  myTreeVariables_.ele_phi -> push_back(cand->phi());
  myTreeVariables_.ele_x ->push_back(cand->momentum().x());
  myTreeVariables_.ele_y ->push_back(cand->momentum().y());
  myTreeVariables_.ele_z ->push_back(cand->momentum().z());
  myTreeVariables_.ele_mass ->push_back(cand-> mass());
  myTreeVariables_.ele_mt ->push_back(cand-> mt());
  myTreeVariables_.ele_pdgId ->push_back(cand->pdgId());
  myTreeVariables_.ele_nDau ->push_back(cand->numberOfDaughters());
}



void ElectronTreeFiller::writeEcalInfo(const Candidate *cand, SuperClusterRef sclusRef,
				       const GsfElectron* electronRef, TreeContent& myTreeVariables_) 
{

  if(&sclusRef) {
    // Cluster related variables
    myTreeVariables_.ele_nClu->push_back(sclusRef->clustersSize());

    myTreeVariables_.ele_ScEraw->push_back(sclusRef->rawEnergy());
    myTreeVariables_.ele_ScEcal->push_back(sclusRef->energy());
    
    myTreeVariables_.ele_ScEta->push_back(sclusRef->eta());
    myTreeVariables_.ele_ScPhi->push_back(sclusRef->phi());
  }

  //EleID
  //&&&&&   // accessors

  myTreeVariables_.ele_covEtaEta -> push_back(electronRef->sigmaEtaEta());
  myTreeVariables_.ele_coviEtaiEta -> push_back(electronRef->sigmaIetaIeta());

  myTreeVariables_.ele_e1x5->push_back(electronRef->e1x5());
  myTreeVariables_.ele_e2x5Max->push_back(electronRef->e2x5Max());
  myTreeVariables_.ele_e5x5->push_back(electronRef->e5x5());
  myTreeVariables_.ele_HoE->push_back(electronRef->hcalOverEcal());

  myTreeVariables_.ele_ScEoP->push_back(electronRef->eSuperClusterOverP());
  myTreeVariables_.ele_SeedEoP->push_back(electronRef->eSeedClusterOverP());
  myTreeVariables_.ele_SeedEoPout->push_back(electronRef->eSeedClusterOverPout());
  myTreeVariables_.ele_eEleClusteroPout->push_back(electronRef->eEleClusterOverPout());

  myTreeVariables_.ele_deltaEtaSuperClusterTrackAtVtx -> push_back(electronRef->deltaEtaSuperClusterTrackAtVtx());
  myTreeVariables_.ele_deltaEtaSeedClusterTrackAtCalo  -> push_back(electronRef->deltaEtaSeedClusterTrackAtCalo() );
  myTreeVariables_.ele_deltaEtaEleClusterTrackAtCalo  -> push_back(electronRef->deltaEtaEleClusterTrackAtCalo() );
  myTreeVariables_.ele_deltaPhiSuperClusterTrackAtVtx  -> push_back(electronRef->deltaPhiSuperClusterTrackAtVtx() );
  myTreeVariables_.ele_deltaPhiSeedClusterTrackAtCalo  -> push_back(electronRef->deltaPhiSeedClusterTrackAtCalo() );
  myTreeVariables_.ele_deltaPhiEleClusterTrackAtCalo  -> push_back(electronRef->deltaPhiEleClusterTrackAtCalo() );



  //Ele Isolation
  //03
  myTreeVariables_.ele_dr03TkSumPt -> push_back(electronRef->dr03TkSumPt());
  myTreeVariables_.ele_dr03EcalRecHitSumEt -> push_back(electronRef->dr03EcalRecHitSumEt());
  myTreeVariables_.ele_dr03HcalDepth1TowerSumEt -> push_back(electronRef->dr03HcalDepth1TowerSumEt());
  myTreeVariables_.ele_dr03HcalDepth2TowerSumEt -> push_back(electronRef->dr03HcalDepth2TowerSumEt());
  myTreeVariables_.ele_dr03HcalTowerSumEt -> push_back(electronRef->dr03HcalTowerSumEt());

  //04
  myTreeVariables_.ele_dr04TkSumPt -> push_back(electronRef->dr04TkSumPt());
  myTreeVariables_.ele_dr04EcalRecHitSumEt -> push_back(electronRef->dr04EcalRecHitSumEt());
  myTreeVariables_.ele_dr04HcalDepth1TowerSumEt -> push_back(electronRef->dr04HcalDepth1TowerSumEt());
  myTreeVariables_.ele_dr04HcalDepth2TowerSumEt -> push_back(electronRef->dr04HcalDepth2TowerSumEt());
  myTreeVariables_.ele_dr04HcalTowerSumEt -> push_back(electronRef->dr04HcalTowerSumEt());

}


void ElectronTreeFiller::writeTrkInfo(const Candidate* cand, const GsfElectron* electronRef, 
				      const GsfTrackRef trkRef, TreeContent& myTreeVariables_) 
{

  if(&trkRef) {

    myTreeVariables_.ele_trackPositionAtVtxX -> push_back(electronRef->trackPositionAtVtx().x());
    myTreeVariables_.ele_trackPositionAtVtxY -> push_back(electronRef->trackPositionAtVtx().y());
    myTreeVariables_.ele_trackPositionAtVtxZ -> push_back(electronRef->trackPositionAtVtx().z());

    myTreeVariables_.ele_trackPositionAtCaloX -> push_back(electronRef->trackPositionAtCalo().x());
    myTreeVariables_.ele_trackPositionAtCaloY -> push_back(electronRef->trackPositionAtCalo().y());
    myTreeVariables_.ele_trackPositionAtCaloZ -> push_back(electronRef->trackPositionAtCalo().z());


    myTreeVariables_.ele_trackMomentumAtVtxX -> push_back(electronRef->trackMomentumAtVtx().x());
    myTreeVariables_.ele_trackMomentumAtVtxY -> push_back(electronRef->trackMomentumAtVtx().y());
    myTreeVariables_.ele_trackMomentumAtVtxZ -> push_back(electronRef->trackMomentumAtVtx().z());

    myTreeVariables_.ele_trackMomentumAtCaloX -> push_back(electronRef->trackMomentumAtCalo().x());
    myTreeVariables_.ele_trackMomentumAtCaloY -> push_back(electronRef->trackMomentumAtCalo().y());
    myTreeVariables_.ele_trackMomentumAtCaloZ -> push_back(electronRef->trackMomentumAtCalo().z());

    myTreeVariables_.ele_VtxX -> push_back(trkRef->vertex().x());
    myTreeVariables_.ele_VtxY -> push_back(trkRef->vertex().y());
    myTreeVariables_.ele_VtxZ -> push_back(trkRef->vertex().z());

  }

  /*
  else {

    myTreeVariables_.ele_trackPositionAtVtxX -> push_back();
    myTreeVariables_.ele_trackPositionAtVtxY -> push_back(electronRef.ele_trackPositionAtVtx().y());
    myTreeVariables_.ele_trackPositionAtVtxZ -> push_back(electronRef.ele_trackPositionAtVtx().z());

    myTreeVariables_.ele_trackPositionAtCaloX -> push_back(electronRef.ele_trackPositionAtCalo().x());
    myTreeVariables_.ele_trackPositionAtCaloY -> push_back(electronRef.ele_trackPositionAtCalo().y());
    myTreeVariables_.ele_trackPositionAtCaloZ -> push_back(electronRef.ele_trackPositionAtCalo().z());


    myTreeVariables_.ele_trackMomentumAtVtxX -> push_back(electronRef.ele_trackMomentumAtVtx().x());
    myTreeVariables_.ele_trackMomentumAtVtxY -> push_back(electronRef.ele_trackMomentumAtVtx().y());
    myTreeVariables_.ele_trackMomentumAtVtxZ -> push_back(electronRef.ele_trackMomentumAtVtx().z());

    myTreeVariables_.ele_trackMomentumAtCaloX -> push_back(electronRef.ele_trackMomentumAtCalo().x());
    myTreeVariables_.ele_trackMomentumAtCaloY -> push_back(electronRef.ele_trackMomentumAtCalo().y());
    myTreeVariables_.ele_trackMomentumAtCaloZ -> push_back(electronRef.ele_trackMomentumAtCalo().z());


  }
  */
}
