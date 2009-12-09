#include <FWCore/PluginManager/interface/ModuleDef.h>
#include <FWCore/Framework/interface/MakerMacros.h>



#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "PhysicsTools/UtilAlgos/interface/ObjectSelector.h"
#include "PhysicsTools/UtilAlgos/interface/SortCollectionSelector.h"
#include "PhysicsTools/UtilAlgos/interface/PtMinSelector.h"
#include "PhysicsTools/UtilAlgos/interface/SingleElementCollectionSelector.h"
#include "PhysicsTools/UtilAlgos/interface/ObjectCountFilter.h"
#include "PhysicsTools/UtilAlgos/interface/SingleObjectSelector.h"

#include "PhysicsTools/UtilAlgos/interface/ObjectSelector.h"
#include "HiggsAnalysis/WWTo2l2nu/plugins/WWEleAmbiguityResolve.h"
#include "HiggsAnalysis/WWTo2l2nu/plugins/WWTreeProducer.h"


//#include "PhysicsTools/UtilAlgos/interface/SingleObjectRefVectorSelector.h"
#include "DataFormats/Common/interface/RefVector.h" 

#include "HiggsAnalysis/WWTo2l2nu/plugins/WWElectronSelector.h"
#include "HiggsAnalysis/WWTo2l2nu/plugins/WWMuonSelector.h"
#include "HiggsAnalysis/WWTo2l2nu/plugins/WWElectronIsolator.h"
#include "HiggsAnalysis/WWTo2l2nu/plugins/WWMuonIsolator.h"
//#include "HiggsAnalysis/WWTo2l2nu/plugins/WWJetCleaner.h"

typedef ObjectSelector<WWElectronSelector> WWElectronSelection;
typedef ObjectSelector<WWElectronSelector, edm::RefVector<reco::GsfElectronCollection> > WWElectronSelectionRef;
typedef ObjectSelector<WWMuonSelector> WWMuonSelection;
typedef ObjectSelector<WWMuonSelector, edm::RefVector<reco::MuonCollection> > WWMuonSelectionRef;

typedef ObjectSelector<WWElectronIsolator> WWElectronIsolation;
typedef ObjectSelector<WWElectronIsolator, edm::RefVector<reco::GsfElectronCollection> > WWElectronIsolationRef;
typedef ObjectSelector<WWMuonIsolator> WWMuonIsolation;
typedef ObjectSelector<WWMuonIsolator, edm::RefVector<reco::MuonCollection> > WWMuonIsolationRef;
//Jet cleaning
//typedef ObjectSelector<WWJetCleaner> WWJetCleaning;
//typedef ObjectSelector<WWJetCleaner, edm::RefVector<reco::CaloJetCollection> > WWJetCleaningRef;

typedef ObjectSelector< WWEleAmbiguityResolve, reco::GsfElectronRefVector > AmbResolver ;




DEFINE_SEAL_MODULE () ;
DEFINE_ANOTHER_FWK_MODULE(WWElectronSelection);
DEFINE_ANOTHER_FWK_MODULE(WWElectronSelectionRef);
DEFINE_ANOTHER_FWK_MODULE(WWMuonSelection);
DEFINE_ANOTHER_FWK_MODULE(WWMuonSelectionRef);

DEFINE_ANOTHER_FWK_MODULE(WWElectronIsolation);
DEFINE_ANOTHER_FWK_MODULE(WWElectronIsolationRef);
DEFINE_ANOTHER_FWK_MODULE(WWMuonIsolation);
DEFINE_ANOTHER_FWK_MODULE(WWMuonIsolationRef);
//DEFINE_ANOTHER_FWK_MODULE(WWJetCleaning);
//DEFINE_ANOTHER_FWK_MODULE(WWJetCleaningRef) ;


DEFINE_ANOTHER_FWK_MODULE (AmbResolver) ;
DEFINE_ANOTHER_FWK_MODULE (WWTreeProducer);
