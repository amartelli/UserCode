#ifndef ElectronTreeFiller_h
#define ElectronTreeFiller_h

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

/*
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
*/

#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
//#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
//#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
//#include "DataFormats/EgammaCandidates/interface/Electron.h"


#include "HiggsAnalysis/WWTo2l2nu/interface/TreeContent.h"
//#include "HiggsAnalysis/WWTo2l2nu/interface/CandidateTreeFiller.h"
#include <TTree.h>

using namespace cms;
using namespace edm;
using namespace reco;



class ElectronTreeFiller {

 public:

  //! Dump everything                                                                                                                                
  ElectronTreeFiller(edm::InputTag recoCollection,
                     const edm::Event&, const edm::EventSetup&);


  //! Dump  everything if fatTree is true and less informations otherwise                                                                            
    
  //! Destructor
  virtual ~ElectronTreeFiller();

  // ----------------------------- //
  
  void writeCollectionToTree(TreeContent& myTreeVariables_);


 private:

  void writeCandInfo(const Candidate *cand, TreeContent& myTreeVariables_);
  void writeEcalInfo(const Candidate *cand, const SuperClusterRef sclusRef,
		     const GsfElectron* electronRef, TreeContent& myTreeVariables_);

  void writeTrkInfo(const Candidate* cand, const GsfElectron* electronRef, 
		    const GsfTrackRef trkRef, TreeContent& myTreeVariables_);

  edm::InputTag recoEleCollectionTag_;
  const edm::Event& iEvent_;
  const edm::EventSetup& iSetup_;

  
};

#endif // ElectronTreeFiller_h
