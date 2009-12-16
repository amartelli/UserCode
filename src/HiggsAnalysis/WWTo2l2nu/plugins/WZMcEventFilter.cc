#include "HiggsAnalysis/WWTo2l2nu/plugins/WZMcEventFilter.h"
#include <iostream>

//#include <map>
#include <string>
#include <vector>

using namespace std;

WZMcEventFilter::WZMcEventFilter(const edm::ParameterSet& iConfig):
  mctype_(iConfig.getUntrackedParameter<std::string>( "mctype")),
  mcTruthCollection_(iConfig.getParameter<edm::InputTag>("mcTruthCollection"))
{
  edm::Service<TFileService> fs;
  
  m_MCtotalEvents = fs -> make<TH1F>("MCtotalEvents", "MCtotalEvents", 1,  0., 1.);
  m_MCpassedEvents = fs -> make<TH1F>("MCpassedEvents", "MCpassedEvents", 1,  0., 1.);
  m_MCfilterEfficiency = fs -> make<TH1F>("MCfilterEfficiency", "MCfilterEfficiency", 1,  0., 1.);
}


WZMcEventFilter::~WZMcEventFilter()
{}

bool WZMcEventFilter::filter(edm::Event& e, const edm::EventSetup& iSetup) 
{
  
  if(mctype_ != "signal") return true;
  
  m_MCtotalEvents -> Fill(0.5);

  // 1. make a loop untill both W and Z are found
  // 2. then follow the decay tree
  
  vector<const reco::Candidate*> tmpW;
  vector<const reco::Candidate*> tmpZ;
  vector<const reco::Candidate*> Wdecays;
  vector<const reco::Candidate*> Zdecays;
  

  edm::Handle<reco::GenParticleCollection> genCandidates_;
  e.getByLabel(mcTruthCollection_, genCandidates_);
  if(! (genCandidates_.isValid ()) )
    {
      std::cerr << "WZMcEventFilter::Warning: " << genCandidates_ << " not available" << std::endl;
      return false;
    }
  

  reco::GenParticleCollection::const_iterator p = genCandidates_->begin();
  
  while (p != genCandidates_->end()){
    if (abs(p->pdgId()) == 24 && p->status() == 3 ) tmpW.push_back(&(*p));
    if (abs(p->pdgId()) == 23 && p->status() == 3 ) tmpZ.push_back(&(*p));
    p++;
  }
  
  if (tmpW.size() != 1 || tmpZ.size() != 1 ) {
    std::cout << "WZMcEventFilter::Filter(): no WZ bosons in the event!" << endl;
    return false;
  } 
  
  realW_ = tmpW.at(0);
  realZ_ = tmpZ.at(0);
  
  vector<const reco::Candidate*> da_realW, da_realZ;  

  for (unsigned int bb = 0; bb < realW_->numberOfDaughters(); ++bb){
    if (realW_->daughter(bb)->status() != 2 ) da_realW.push_back(realW_->daughter(bb));
  }
  for (unsigned int bb = 0; bb < realZ_->numberOfDaughters(); ++bb){
    if (realZ_->daughter(bb)->status() != 2 )  da_realZ.push_back(realZ_->daughter(bb));
  }

  if ( (da_realW.size() != 2) || (da_realZ.size() != 2) ){
    std::cout <<"WZMcEventFilter::Filter(): not a WZ -> 4 body decay!"<<endl;
    return false;
  }

  int da1_rW_pdgID = da_realW.at(0) ->pdgId();
  int da2_rW_pdgID = da_realW.at(1) ->pdgId();
  int da1_rZ_pdgID = da_realZ.at(0) ->pdgId();
  int da2_rZ_pdgID = da_realZ.at(1) ->pdgId();
  
  if(abs(da1_rW_pdgID) ==  15 || abs(da1_rW_pdgID) == 16 || abs(da2_rW_pdgID) == 15 || abs(da2_rW_pdgID) == 16 || 
     abs(da1_rZ_pdgID) ==  15 || abs(da1_rZ_pdgID) == 16 || abs(da2_rZ_pdgID) == 15 || abs(da2_rZ_pdgID) == 16) return false;

  m_MCpassedEvents -> Fill(0.5);

  int nTotalEvents = static_cast<int>(m_MCtotalEvents -> GetBinContent(1));
  int nPassedEvents = static_cast<int>(m_MCpassedEvents -> GetBinContent(1));

  m_MCfilterEfficiency -> SetBinContent(1, 1.* nPassedEvents/nTotalEvents);
  

  return true;
}



void WZMcEventFilter::beginJob(const edm::EventSetup&)
{}


void WZMcEventFilter::endJob() 
{}

//define this as a plug-in                                                                                                                           
DEFINE_FWK_MODULE(WZMcEventFilter);
