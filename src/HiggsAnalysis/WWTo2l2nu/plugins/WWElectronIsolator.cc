//#include "HiggsAnalysis/WWTo2l2nu/plugins/EleTrackerIsolationAlgo.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"
#include "HiggsAnalysis/WWTo2l2nu/plugins/WWElectronIsolator.h"


WWElectronIsolator::WWElectronIsolator(const edm::ParameterSet& iConfig)
{
  theTrackIsolCut_     = iConfig.getParameter<double>("TrackIsolCut");
  theAbsTrackIsolCut_     = iConfig.getParameter<double>("AbsTrackIsolCut");
  theCaloIsolCut_     = iConfig.getParameter<double>("CaloIsolCut");
  theECALIsolCut_     = iConfig.getParameter<double>("ECALIsolCut");
  //  absolute_                  = iConfig.getParameter<bool>("absolute");

  edm::Service<TFileService> fs;

  m_SumPt_over_Pt_EleTk = fs -> make<TH1F>("SumPt_over_Pt_EleTk", "SumPt_over_Pt_EleTk", 1000,  0., 100.);
  m_SumPt_EleTk = fs -> make<TH1F>("SumPt_EleTk", "SumPt_EleTk", 1000,  0., 100.);
  m_SumPt_EleCalo = fs -> make<TH1F>("SumPt_EleCalo", "SumPt_EleCalo", 1000,  0., 100.);
  m_SumPt_EleEcal = fs -> make<TH1F>("SumPt_EleEcal", "SumPt_EleEcal", 1000,  0., 100.);
}


WWElectronIsolator::~WWElectronIsolator()
{}

void WWElectronIsolator::select(edm::Handle<reco::GsfElectronCollection> electrons,
                                 const edm::Event& iEvent,
		                 const edm::EventSetup& iEventSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace std;

//   if ( !absolute_ && theTrackIsolCut_ > 1.0 ) {
//     LogWarning("HWWElectronIsolator") << "Requested relative electron tracker isolation with cut: " 
//                                       << theTrackIsolCut_ << " > 1.0. Possible misconfiguration...";
//   }

  selected_.clear();
  
  for(unsigned i =0; i<electrons->size(); ++i) {
    
    Ref<reco::GsfElectronCollection> electronRef(electrons,i);
    double pt = electronRef->pt();
    double sumPt = electronRef->dr03TkSumPt();
    double sumPt_o_Pt = electronRef->dr03TkSumPt() / pt;
    double sumEcal = electronRef->dr03EcalRecHitSumEt();
    double sumHcal = electronRef->dr03HcalTowerSumEt();

    m_SumPt_EleTk ->Fill(sumPt);
    m_SumPt_over_Pt_EleTk->Fill(sumPt_o_Pt);
    m_SumPt_EleCalo ->Fill(sumEcal+sumHcal);
    m_SumPt_EleEcal->Fill(sumEcal);

    if ( sumPt > theAbsTrackIsolCut_ ) continue;
    if ( sumPt_o_Pt > theTrackIsolCut_ ) continue;
    if ( sumEcal > theECALIsolCut_ ) continue;
    if ( sumHcal > theCaloIsolCut_ ) continue;


    selected_.push_back(electronRef);

   


  }
}
