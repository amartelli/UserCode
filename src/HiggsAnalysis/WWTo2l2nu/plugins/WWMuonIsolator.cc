// system include files
#include <memory>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"

#include "HiggsAnalysis/WWTo2l2nu/plugins/WWMuonIsolator.h"


WWMuonIsolator::WWMuonIsolator(const edm::ParameterSet& iConfig)
{
  using namespace edm;  

  selectedMuonsRefLabel_=iConfig.getParameter<InputTag>("SelectedMuonRefCollectionLabel");

  theTrackIsolCut   = iConfig.getParameter<double>("trackIsolCut"); 
  theCaloIsolCut    = iConfig.getParameter<double>("caloIsolCut"); 
  doRefCheck_=iConfig.getParameter<bool>("doRefCheck");


  // *****************  // nota: per aggiungere le altre var di isolamento
  // *****************  //  ex in ECAL in HCAl in TK o combined
  // *****************  //


  //produces<reco::MuonCollection>();
}

WWMuonIsolator::~WWMuonIsolator()
{}

void WWMuonIsolator::select(edm::Handle<reco::MuonCollection> muons,
			    const edm::Event& iEvent,
			    const edm::EventSetup& iEventSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace std;

  selected_.clear();
  double muonTrackerDeposit30 = -9999;  
  double muonEcalDeposit30    = -9999; 
  double muonHcalDeposit30    = -9999;

  Handle<RefVector<MuonCollection> >muonsRef;
  if(doRefCheck_==true)
    iEvent.getByLabel(selectedMuonsRefLabel_,muonsRef);


  bool muTrackIsol=false;
  bool muCaloIsol=false;

  for (unsigned i = 0;i<muons->size(); i++){ 
    muTrackIsol = false;
    muCaloIsol = false;
      
    muonTrackerDeposit30= (*muons)[i].isolationR03().sumPt;   
    muonEcalDeposit30   = (*muons)[i].isolationR03().emEt;   
    muonHcalDeposit30   = (*muons)[i].isolationR03().hadEt; 
      
    // *** tracker isolation cut
    if(muonTrackerDeposit30 < theTrackIsolCut){ muTrackIsol = true;}
      
    // *** calo isolation cut 
    if( muonEcalDeposit30+muonHcalDeposit30 < theCaloIsolCut){ muCaloIsol = true;}

    Ref<MuonCollection>muonRAWRef(muons,i);
    bool selected=true;
    if(doRefCheck_==true)
      if (find(muonsRef->begin(), muonsRef->end(),muonRAWRef)==muonsRef->end())
	{
	  selected=false;
	  //cout<<"Isolated muon without ID"<<endl;
	}
      
    if(muTrackIsol==true && muCaloIsol == true && selected==true)
      {
	selected_.push_back(muonRAWRef);
      }
  } // loop over muons

}
//define this as a plug-in
//DEFINE_FWK_MODULE(WWMuonIsolator);
