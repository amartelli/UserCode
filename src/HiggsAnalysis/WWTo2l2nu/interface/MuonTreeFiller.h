#ifndef MuonTreeFiller_h
#define MuonTreeFiller_h

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "TrackingTools/TrackAssociator/interface/CachedTrajectory.h"
#include "TrackingTools/TrackAssociator/interface/CaloDetIdAssociator.h"
#include "TrackingTools/TrackAssociator/interface/EcalDetIdAssociator.h"
#include "TrackingTools/TrackAssociator/interface/MuonDetIdAssociator.h"
#include "TrackingTools/TrackAssociator/interface/HcalDetIdAssociator.h"
#include "TrackingTools/TrackAssociator/interface/HODetIdAssociator.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "HiggsAnalysis/WWTo2l2nu/interface/TreeContent.h"
//#include "HiggsAnalysis/WWTo2l2nu/interface/CandidateTreeFiller.h"
#include <TTree.h>

using namespace cms;
using namespace edm;
using namespace reco;


class MuonTreeFiller {

 public:

  // Dump everything

  MuonTreeFiller(edm::InputTag collectionTag,
		 const edm::Event&, const edm::EventSetup&);


  
  // Destructor
  virtual ~MuonTreeFiller();


  //! write the muon related informations for the given collection

  void writeCollectionToTree(TreeContent& myTreeVariables_);

  
 private:

  void writeCandInfo(const Candidate *cand, TreeContent& myTreeVariables_);

  void writeTrkInfo(const Candidate *cand, const Muon *muon, TreeContent& myTreeVariables_);
  void writeMuonInfo(const Muon *muon, TreeContent& myTreeVariables_);

  edm::InputTag collectionTag_;
  const edm::Event& iEvent_;
  const edm::EventSetup& iSetup_;

  // Geometry
  EcalDetIdAssociator ecalDetIdAssociator_;
  HcalDetIdAssociator hcalDetIdAssociator_;
  HODetIdAssociator   hoDetIdAssociator_;
  CaloDetIdAssociator caloDetIdAssociator_;
  MuonDetIdAssociator muonDetIdAssociator_;


  CachedTrajectory cachedTrajectory_;
  edm::ESHandle<MagneticField> bField;
  edm::ESHandle<CaloGeometry> theCaloGeometry_;
  edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry_;

};

#endif // MuonTreeFiller_h
