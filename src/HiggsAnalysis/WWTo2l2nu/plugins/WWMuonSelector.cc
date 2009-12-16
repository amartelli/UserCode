#include <memory>
#include <iostream>
#include "HiggsAnalysis/WWTo2l2nu/plugins/WWMuonSelector.h"
#include "DataFormats/MuonReco/interface/Muon.h"


WWMuonSelector::WWMuonSelector(const edm::ParameterSet& pset)
{
  muonPtMinBarrel_  = pset.getParameter<double>("muonPtMinBarrel");
  muonPtMinEndcap_  = pset.getParameter<double>("muonPtMinEndcap");
  muonPMinEndcap_   = pset.getParameter<double>("muonPMinEndcap");
  muonEtaMax_       = pset.getParameter<double>("muonEtaMax");
}


WWMuonSelector::~WWMuonSelector()
{}


void WWMuonSelector::select(edm::Handle<reco::MuonCollection> muons,
                        const edm::Event& iEvent,
                        const edm::EventSetup& iEventSetup)
{
  using namespace edm;
  using namespace reco;
  selected_.clear();

  //  std::cout << " WWMuonSelector::select " << muons->size() << std::endl;


  for(unsigned int i=0; i<muons->size(); ++i)
    {
      bool isEndcap = false;
      bool isBarrel = false;
      
      if ( fabs( (*muons)[i].eta() ) > muonEtaMax_ ) continue;
      if ( fabs( (*muons)[i].eta() ) < 1.1 ) {
	isBarrel = true;
      } else {
	isEndcap = true;
      }

      if ( isEndcap && ( (*muons)[i].pt() < muonPtMinEndcap_ ||
			 (*muons)[i].p() < muonPMinEndcap_ )) continue;
      if ( isBarrel && (*muons)[i].pt() < muonPtMinBarrel_ ) continue;
      
      Ref<MuonCollection> muonRef(muons,i);
      selected_.push_back(muonRef);
    }
 
}

 
//define this as a plug-in
//DEFINE_FWK_MODULE(WWMuonSelector);
