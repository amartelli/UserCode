#include "DataFormats/EgammaCandidates/interface/GsfElectronIsoCollection.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"
#include "HiggsAnalysis/WWTo2l2nu/plugins/WWElectronIsolator.h"


WWElectronIsolator::WWElectronIsolator(const edm::ParameterSet& iConfig)
{
  selectedElectronsRefLabel_ = iConfig.getParameter<edm::InputTag>("SelectedElectronRefCollectionLabel");
  trackIsolationProducer_    = iConfig.getParameter<edm::InputTag>("TrackIsolationProducerLabel");
  doRefCheck_		     = iConfig.getParameter<bool>("doRefCheck");
  theTrackIsolCut_	     = iConfig.getParameter<double>("TrackIsolCut");
}


WWElectronIsolator::~WWElectronIsolator()
{
}

void WWElectronIsolator::select(edm::Handle<reco::GsfElectronCollection> electrons,
                                 const edm::Event& iEvent,
		                 const edm::EventSetup& iEventSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace std;

  selected_.clear();
  Handle<RefVector<GsfElectronCollection> > electronsRef;
  //  std::cout << "BBBBBBBBBBBBBBBBBBBBBBBBBBBBB  sono in  WWElectronIsolator::select " << std::endl;

  edm::Handle< edm::ValueMap<double> > tkIsolationHandle;
  try { iEvent.getByLabel(trackIsolationProducer_, tkIsolationHandle); }
  catch ( cms::Exception& ex ) { printf("Can't get tracker isolation product\n"); }

  const edm::ValueMap<double>& tkIsolationVal = *tkIsolationHandle;

  if(doRefCheck_ == true)
    iEvent.getByLabel(selectedElectronsRefLabel_,electronsRef);

  for(unsigned i =0; i<electrons->size(); i++) {
    
    Ref<reco::GsfElectronCollection> electronRAWRef(electrons,i);
    double sumPtOverEt = tkIsolationVal[electronRAWRef];
     
    bool selected = true;
    if(doRefCheck_ == true) {
      if (find(electronsRef->begin(), electronsRef->end(),electronRAWRef)==electronsRef->end()) {
	selected = false;
      }
    }

    //    std::cout << "sumPtOverEt =  " << sumPtOverEt << "  mentre theTrackIsolCut_ =  " << theTrackIsolCut_ << std::endl; 

    if (sumPtOverEt < theTrackIsolCut_ && selected == true)
      selected_.push_back(electronRAWRef);
  }
}
