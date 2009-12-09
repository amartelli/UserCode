// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "AnalysisDataFormats/Egamma/interface/ElectronID.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "HiggsAnalysis/WWTo2l2nu/plugins/WWElectronSelector.h"

#include <string.h>

WWElectronSelector::WWElectronSelector(const edm::ParameterSet& iConfig)
{
  electronIdCutsLabel_ = iConfig.getParameter<edm::InputTag>("electronIdCutsLabel");
  //  electronIdLikelihoodLabel_ = iConfig.getParameter<edm::InputTag>("electronIdLikelihoodLabel"); 
  useCuts_ = iConfig.getParameter<bool>("useCuts");
  //  likelihoodThreshold_ = iConfig.getUntrackedParameter<double>("likelihoodThreshold",-1.); //0.5 
}


WWElectronSelector::~WWElectronSelector()
{}


void
WWElectronSelector::select (edm::Handle<reco::GsfElectronCollection> electrons,
                             const edm::Event& iEvent,
		             const edm::EventSetup& iEventSetup)
{

  using namespace edm;
  using namespace reco;
  selected_.clear();

  std::vector<edm::Handle<edm::ValueMap<float> > > eIDValueMap(2);
  
  if( iEvent.getByLabel( electronIdCutsLabel_ , eIDValueMap[0] )){


//   if( iEvent.getByLabel( electronIdCutsLabel_ , eIDValueMap[0] )  &&
//       iEvent.getByLabel( electronIdLikelihoodLabel_ , eIDValueMap[1] ) ) {

    const edm::ValueMap<float> & eIdmapCuts = * eIDValueMap[0] ;
    //    const edm::ValueMap<float> & eIdmapLikelihood = * eIDValueMap[1] ;

    // Loop over electrons
    for (unsigned int i = 0; i < electrons->size(); i++) {	  
      Ref<reco::GsfElectronCollection> electronRef(electrons,i);
      bool eleid = false;

      if( useCuts_ ) {
	if( eIdmapCuts[electronRef] > -1. )   // 0.
	  eleid = true;
      }
      //else {
      //	float likelihood = eIdmapLikelihood[electronRef];
      //	if( likelihood > likelihoodThreshold_ ) eleid = true;
      //      }
      
      if (eleid==true)
	selected_.push_back (electronRef);
      
    }
  } else {
    LogWarning("WWElectronSelector") << electronIdCutsLabel_ << " or " << electronIdLikelihoodLabel_ << " not available";
  }

}

