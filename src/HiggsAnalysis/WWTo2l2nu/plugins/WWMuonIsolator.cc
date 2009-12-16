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

  selectedMuonsRefLabel_ = iConfig.getParameter<InputTag>("SelectedMuonRefCollectionLabel");

  theTrackIsolCut   = iConfig.getParameter<double>("trackIsolCut"); 
  theCaloIsolCut    = iConfig.getParameter<double>("caloIsolCut"); 
  doRefCheck_ = iConfig.getParameter<bool>("doRefCheck");

  edm::Service<TFileService> fs;

  m_SumPt_over_Pt_MuonTk = fs -> make<TH1F>("SumPt_over_Pt_MuonTk", "SumPt_over_Pt_MuonTk", 1000,  0., 100.);
  m_SumPt_MuonTk = fs -> make<TH1F>("SumPt_MuonTk", "SumPt_MuonTk", 1000,  0., 100.);
  m_SumPt_MuonCalo = fs -> make<TH1F>("SumPt_MuonCalo", "SumPt_MuonCalo", 1000,  0., 100.);
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
  double muonTrackerDeposit30_overPt = -9999;

  Handle<RefVector<MuonCollection> >muonsRef;
  if(doRefCheck_ == true)
    iEvent.getByLabel(selectedMuonsRefLabel_,muonsRef);


  bool muTrackIsol = false;
  bool muCaloIsol = false;

  //  std::cout << " WWMuonIsolator::select " << muons->size() << std::endl;


  for (unsigned i = 0; i<muons->size(); ++i){ 
    muTrackIsol = false;
    muCaloIsol = false;
      
    muonTrackerDeposit30 = (*muons)[i].isolationR03().sumPt;   
    muonEcalDeposit30   = (*muons)[i].isolationR03().emEt;   
    muonHcalDeposit30   = (*muons)[i].isolationR03().hadEt; 
    muonTrackerDeposit30_overPt = muonTrackerDeposit30 / (*muons)[i].pt() ;

    m_SumPt_over_Pt_MuonTk->Fill(muonTrackerDeposit30_overPt);
    m_SumPt_MuonTk->Fill(muonTrackerDeposit30);
    m_SumPt_MuonTk->Fill(muonEcalDeposit30+muonHcalDeposit30);

    // *** tracker isolation cut
    if(muonTrackerDeposit30_overPt < theTrackIsolCut){ muTrackIsol = true;}
    //    else std::cout << " }}}}}}}}}}}}}muonTrackerDeposit30_overPt " << muonTrackerDeposit30_overPt << std::endl;
    // *** calo isolation cut 
    if( muonEcalDeposit30 + muonHcalDeposit30 < theCaloIsolCut){ muCaloIsol = true;}
    //    else std::cout << " }}}}}}}}}}}}}muonEcalDeposit30 + muonHcalDeposit30 " << muonEcalDeposit30 + muonHcalDeposit30 << std::endl;

    Ref<MuonCollection> muonRAWRef(muons,i);
    bool selected = true;
    if(doRefCheck_ == true)
      if (find(muonsRef->begin(), muonsRef->end(),muonRAWRef) == muonsRef->end())
	{
	  selected = false;
	  //cout<<"Isolated muon without ID"<<endl;
	}
      
    //    if(muTrackIsol == true && muCaloIsol == true && selected == true)
    if(muTrackIsol == true && muCaloIsol == true )
      {
	selected_.push_back(muonRAWRef);
      }
  } // loop over muons

}
//define this as a plug-in
//DEFINE_FWK_MODULE(WWMuonIsolator);
