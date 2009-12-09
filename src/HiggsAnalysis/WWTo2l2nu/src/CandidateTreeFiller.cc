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

//#include "DataFormats/JetReco/interface/CaloJet.h"
//#include "DataFormats/JetReco/interface/CaloJetCollection.h"
//#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/OverlapChecker.h"

#include "HiggsAnalysis/WWTo2l2nu/interface/TreeContent.h"
//#include "HiggsAnalysis/WWTo2l2nu/interface/EleIDTreeFiller.h"
#include "HiggsAnalysis/WWTo2l2nu/interface/CandidateTreeFiller.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


#include <TTree.h>

#include <string>

using namespace edm;
using namespace reco;



CandidateTreeFiller::CandidateTreeFiller(edm::InputTag recoCollectionTag,
                                         const edm::Event& iEvent, const edm::EventSetup& iSetup):
  recoCollectionTag_(recoCollectionTag), iEvent_(iEvent), iSetup_(iSetup)
{
  Initialise();
}


CandidateTreeFiller::~CandidateTreeFiller() {}


// -------------------------------------------------- //

void CandidateTreeFiller::writeCollectionToTree(TreeContent& myTreeVariables_, std::string name){

  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent_.getByLabel(recoCollectionTag_, collectionHandle); }
  catch ( cms::Exception& ex ) {edm::LogWarning("CandidateTreeFiller") << "Can't get candidate collection: " 
								       << recoCollectionTag_; }
  
  const edm::View<reco::Candidate> *collection = collectionHandle.product();
  
  if(collection) {
    ncand_ = collection->size();

    edm::View<reco::Candidate>::const_iterator cand;
    for(cand=collection->begin(); cand!=collection->end(); ++cand) {
      // fill basic kinematics
      writeCandInfo(&(*cand), ncand_);

      if(name == "met"){
	setTreeBranches(myTreeVariables_.met_charge, myTreeVariables_.met_energy, myTreeVariables_.met_et,
 			myTreeVariables_.met_momentum, myTreeVariables_.met_vertexX,
 			myTreeVariables_.met_vertexY, myTreeVariables_.met_vertexZ,
 			myTreeVariables_.met_theta, myTreeVariables_.met_eta, myTreeVariables_.met_phi,
 			myTreeVariables_.met_x, myTreeVariables_.met_y, myTreeVariables_.met_z,
 			myTreeVariables_.met_mass, myTreeVariables_.met_mt,
 			myTreeVariables_.met_pdgId, myTreeVariables_.met_nDau);
      }

      else if(name == "genmet"){
	setTreeBranches(myTreeVariables_.genmet_charge, myTreeVariables_.genmet_energy, myTreeVariables_.genmet_et,
			myTreeVariables_.genmet_momentum, myTreeVariables_.genmet_vertexX,
			myTreeVariables_.genmet_vertexY, myTreeVariables_.genmet_vertexZ,
			myTreeVariables_.genmet_theta, myTreeVariables_.genmet_eta, myTreeVariables_.genmet_phi,
			myTreeVariables_.genmet_x, myTreeVariables_.genmet_y, myTreeVariables_.genmet_z,
			myTreeVariables_.genmet_mass, myTreeVariables_.genmet_mt,
			myTreeVariables_.genmet_pdgId, myTreeVariables_.genmet_nDau);
      }

    }
  }
  else {
    ncand_ = 0;
  }
  
  if(name == "met") myTreeVariables_.met_ncand = ncand_;
  else if(name == "genmet") myTreeVariables_.genmet_ncand = ncand_;
}


void CandidateTreeFiller::writeCandInfo(const Candidate* cand, int ncand){ 

  cand_ = cand;		
  charge_ = ( (int)cand_->charge() );
  energy_ = (cand_->energy());
  pt_ = (cand_->pt());
  et_ = (cand_->et());
  momentum_ = (cand_->p());
  theta_ = (cand_->theta());
  eta_ = (cand_->eta());
  phi_ = (cand_->phi());
  x_ = (cand_->momentum().x());
  y_ = (cand_->momentum().y());
  z_ = (cand_->momentum().z());
  vertexX_ = (cand_->vx());
  vertexY_ = (cand_->vy());
  vertexZ_ = (cand_->vz());
  mass_ = (cand_->mass());
  mt_ = (cand_->mt());
  pdgId_ = (cand_->pdgId());
  nDau_ = (cand_->numberOfDaughters());

}



void CandidateTreeFiller::setTreeBranches(std::vector<int>* charge, std::vector<float>* energy, std::vector<float>* et, 
					  std::vector<float>* momentum, std::vector<float>* vertexX,  
					  std::vector<float>* vertexY, std::vector<float>* vertexZ, 
					  std::vector<float>* theta, std::vector<float>* eta, std::vector<float>* phi, 
					  std::vector<float>* x, std::vector<float>* y, std::vector<float>* z, 
					  std::vector<float>* mass, std::vector<float>* mt, 
					  std::vector<int>* pdgId, std::vector<int>* nDau)
{

  charge->push_back(charge_);
  energy->push_back(energy_);
  et->push_back(et_);
  momentum->push_back(momentum_);
  vertexX->push_back(vertexX_);
  vertexY->push_back(vertexY_);
  vertexZ->push_back(vertexZ_);
  theta->push_back(theta_);
  eta->push_back(eta_);
  phi->push_back(phi_);
  x->push_back(x_);
  y->push_back(y_);
  z->push_back(z_);
  mass->push_back(mass_);
  mt->push_back(mt_);
  pdgId->push_back(pdgId_);
  nDau->push_back(nDau_);
  
}



void CandidateTreeFiller::Initialise() {
  charge_ = 0;
  energy_ = 0;
  et_ = 0;
  momentum_ = 0;
  theta_ = 0;
  eta_ = 0;
  phi_ = 0;
  x_ = 0;
  y_ = 0;
  z_ = 0;
  vertexX_ = 0;
  vertexY_ = 0;
  vertexZ_ = 0;
  mass_ = 0;
  mt_ = 0;
  pdgId_ = 0;
  ncand_ = 0;
  nDau_ = 0;
}

