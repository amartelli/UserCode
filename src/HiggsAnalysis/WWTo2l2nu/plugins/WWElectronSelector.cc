// system include files
#include <memory>

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
  elecPtMin_   = iConfig.getParameter<double>("electronPtMin");
  elecEtaMax_  = iConfig.getParameter<double>("electronEtaMax");
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

  edm::Handle< edm::ValueMap<float> >  eIDValueMap;  
  if( iEvent.getByLabel( electronIdCutsLabel_ , eIDValueMap )){

  const edm::ValueMap<float> & eIdmapCuts = * eIDValueMap ;
  
  // Loop over electrons
  for (unsigned int i = 0; i < electrons->size(); i++) {	  
    Ref<reco::GsfElectronCollection> electronRef(electrons,i);
    if (electronRef->pt() >= elecPtMin_ && fabs(electronRef->eta()) < elecEtaMax_){

      bool eleid = false;
      if( eIdmapCuts[electronRef] > -1 ) eleid = true; //0
      
      if (eleid == true) selected_.push_back (electronRef);
    }
  }
  } else {
    LogWarning("WWElectronSelector") << electronIdCutsLabel_ << " not available";
  }
  
}

