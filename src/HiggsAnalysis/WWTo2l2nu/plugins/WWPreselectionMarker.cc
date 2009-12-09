#include "HiggsAnalysis/WWTo2l2nu/plugins/WWPreselectionMarker.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include <iostream>



WWPreselectionMarker::WWPreselectionMarker(const edm::ParameterSet& iConfig)
{
	using namespace edm;
	using namespace reco;

	
	muonslabel_=iConfig.getParameter<InputTag>( "MuonLabel");
	electronslabel_=iConfig.getParameter<InputTag>( "ElectronLabel");
	calometlabel_=iConfig.getParameter<InputTag>( "CaloMetLabel");
	jetslabel_=iConfig.getParameter<InputTag>( "JetLabel");
	
	leptonPtMinMin_=iConfig.getParameter<double>("LeptonPtMinMin");        
	leptonPtMaxMin_=iConfig.getParameter<double>("LeptonPtMaxMin");        
	leptonEtaMax_=iConfig.getParameter<double>("LeptonEtaMax");      
	leptonChargeCombination_=iConfig.getParameter<double>("LeptonChargeCombination");  
	metMin_=iConfig.getParameter<double>("MetMin");  
 	invMassMin_=iConfig.getParameter<double>("InvMassMin");      
	produces<bool>();
}


WWPreselectionMarker::~WWPreselectionMarker()
{}


void WWPreselectionMarker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{       
        using namespace edm;
        using namespace reco;
	using namespace std;

        Handle<MuonCollection> muons;
        Handle<GsfElectronCollection> electrons;
        Handle<CaloMETCollection> met;
        Handle<CaloJetCollection> jets;		
        std::auto_ptr<bool> pOut(new bool);
        *pOut=true;
        
        iEvent.getByLabel(muonslabel_, muons);
        iEvent.getByLabel(electronslabel_, electrons);
        iEvent.getByLabel(calometlabel_, met);
        iEvent.getByLabel(jetslabel_, jets);
        std::vector<const RecoCandidate*> leptons;
        std::vector<int> pdgID_reco;


	//  ----- Lepton selection 
        GsfElectronCollection::const_iterator electron;
        for(electron=electrons->begin(); electron != electrons->end(); ++electron)
        {
	  if(electron->pt() >= leptonPtMinMin_ && fabs(electron->eta()) <= leptonEtaMax_)
	     {
                        const RecoCandidate *lepton=&(*electron);
                        leptons.push_back(lepton);
			int charge = 11;
			if(lepton->charge() == 1 ) charge = -11;
			pdgID_reco.push_back(charge);
	     }
        }
	
        MuonCollection::const_iterator muon;
        for(muon=muons->begin(); muon != muons->end(); ++muon)
	  {
	    if(muon->isGlobalMuon() != 1) continue;
	    if(muon->pt() >= leptonPtMinMin_ && fabs(muon->eta()) <= leptonEtaMax_)
	      {
                        const RecoCandidate *lepton=&(*muon);
                        leptons.push_back(lepton);
			int charge = 13;
			if(lepton->charge() == 1 ) charge = -13;
			pdgID_reco.push_back(charge);
	      }
	  }
	
	// event counter
	selectedEvents[0]++;
  	// find hardest lepton
        if(leptons.size() < 3)
                *pOut=false;
        else
	// 3 good leptons found
        { 

	  selectedEvents[1]++;

	  const RecoCandidate *lepton1 = leptons[0];
	  const RecoCandidate *lepton2 = leptons[1];
	  const RecoCandidate *lepton3 = leptons[2];

	  int pdgID_reco1 = pdgID_reco.at(0);
	  int pdgID_reco2 = pdgID_reco.at(1);
	  int pdgID_reco3 = pdgID_reco.at(2);

	  if(leptons.size() > 3) {

	    std::vector<const RecoCandidate*>::const_iterator lepton;
	    int pdgID_reco_it = 0;
	    for(lepton=leptons.begin(); lepton != leptons.end(); ++lepton){
	      
	      //choose hardest leptons as lepton1 and lepton2 and lepton3 (change the softest one)

	      //change the last
	      if((*lepton)->pt() > lepton1->pt() && (*lepton)->pt() > lepton2->pt() && (*lepton)->pt() > lepton3->pt()) {
		// 1 2 3
		if(lepton1->pt() > lepton2->pt() && lepton1->pt() > lepton3->pt() && lepton2->pt() > lepton3->pt() ) 
		  { lepton3 = (*lepton); pdgID_reco3 = pdgID_reco.at(pdgID_reco_it); } 
		//1 3 2 
		else if(lepton1->pt() > lepton2->pt() && lepton1->pt() > lepton3->pt() && lepton3->pt() > lepton2->pt() ) 
		  { lepton2 = (*lepton); pdgID_reco2 = pdgID_reco.at(pdgID_reco_it); } 
		//2 1 3 
		else if(lepton2->pt() > lepton1->pt() && lepton2->pt() > lepton3->pt() && lepton1->pt() > lepton3->pt() ) 
		  { lepton3 = (*lepton); pdgID_reco.at(2) = pdgID_reco.at(pdgID_reco_it); } 
		//2 3 1
		else if(lepton2->pt() > lepton1->pt() && lepton2->pt() > lepton3->pt() && lepton3->pt() > lepton1->pt() ) 
		  { lepton1 = (*lepton); pdgID_reco1 = pdgID_reco.at(pdgID_reco_it); } 
	      }

	      // change 1 or 2
	      else if( (*lepton)->pt() > lepton1->pt() && (*lepton)->pt() > lepton2->pt() && lepton3->pt() > (*lepton)->pt()){
		// 1 2 change 2
		if(lepton1->pt() > lepton2->pt() ) { lepton2 = (*lepton); pdgID_reco2 = pdgID_reco.at(pdgID_reco_it); } 
		else { lepton1 = (*lepton); pdgID_reco1 = pdgID_reco.at(pdgID_reco_it); }   // change 1
	      }

	      // change 1 or 3
	      else if( (*lepton)->pt() > lepton1->pt() && (*lepton)->pt() > lepton3->pt() && lepton2->pt() > (*lepton)->pt()){
		// 1 3 change 3
		if(lepton1->pt() > lepton3->pt() ) { lepton3 = (*lepton); pdgID_reco3 = pdgID_reco.at(pdgID_reco_it); } 
		else { lepton1 = (*lepton); pdgID_reco1 = pdgID_reco.at(pdgID_reco_it); }   // change 1
	      }

	      // change 2 or 3
	      else if( (*lepton)->pt() > lepton2->pt() && (*lepton)->pt() > lepton3->pt() && lepton1->pt() > (*lepton)->pt()){
		// 2 3 change 3
		if(lepton2->pt() > lepton3->pt() ) { lepton3 = (*lepton); pdgID_reco3 = pdgID_reco.at(pdgID_reco_it); } 
		else { lepton2 = (*lepton); pdgID_reco2 = pdgID_reco.at(pdgID_reco_it); }   // change 2
	      }

	      // change 1
	      else if( (*lepton)->pt() > lepton1->pt() && lepton2->pt() > (*lepton)->pt() && lepton3->pt() > (*lepton)->pt())
		{ lepton1 = (*lepton); pdgID_reco1 = pdgID_reco.at(pdgID_reco_it); } 

	      // change 2
	      else if( (*lepton)->pt() > lepton2->pt() && lepton1->pt() > (*lepton)->pt() && lepton3->pt() > (*lepton)->pt())
		{ lepton2 = (*lepton); pdgID_reco2 = pdgID_reco.at(pdgID_reco_it); } 

	      // change 3
	      else if( (*lepton)->pt() > lepton3->pt() && lepton1->pt() > (*lepton)->pt() && lepton2->pt() > (*lepton)->pt())
		{ lepton3 = (*lepton); pdgID_reco3 = pdgID_reco.at(pdgID_reco_it); } 

	      ++pdgID_reco_it;
	    }
	  }


	  // one lepton with pt>15 GeV
	  if(lepton1->pt() < leptonPtMaxMin_  && lepton2->pt() < leptonPtMaxMin_ && lepton3->pt() < leptonPtMaxMin_) 
	    *pOut=false;
	  else selectedEvents[2]++;

	  // unlike signed leptons
          if(fabs(pdgID_reco1+pdgID_reco2+pdgID_reco3) == 13 || fabs(pdgID_reco1+pdgID_reco2+pdgID_reco3) == 11)
	    selectedEvents[3]++;
	  else  *pOut=false;                
	  
	  //met cut  
	  reco::CaloMETCollection::const_iterator met_iter;
	  for(met_iter=met->begin(); met_iter != met->end(); ++met_iter)
	    {
	      if(met_iter->pt() < metMin_)
		*pOut=false;
	      else selectedEvents[4]++;
	    }

	  // invariant mass cut                
	  double entot = 0.;
	  double pxtot = 0.;
	  double pytot = 0.;
	  double pztot = 0.;
	  double inmass = 0.;
	  double inmass_sqrt = 0.;
	  if(pdgID_reco1+pdgID_reco2 == 0 ){
	    entot = lepton1->energy()+lepton2->energy();
	    pxtot = lepton1->px()+lepton2->px();
	    pytot = lepton1->py()+lepton2->py();
	    pztot = lepton1->pz()+lepton2->pz();
	    inmass = entot*entot - pxtot*pxtot - pytot*pytot - pztot*pztot;
	    inmass_sqrt=sqrt(inmass);                
	  }
	  else if(pdgID_reco1+pdgID_reco3 == 0 ){
	    entot = lepton1->energy()+lepton3->energy();
	    pxtot = lepton1->px()+lepton3->px();
	    pytot = lepton1->py()+lepton3->py();
	    pztot = lepton1->pz()+lepton3->pz();
	    inmass = entot*entot - pxtot*pxtot - pytot*pytot - pztot*pztot;
	    inmass_sqrt=sqrt(inmass);                
	  }
	  else if(pdgID_reco2+pdgID_reco3 == 0 ){
	    entot = lepton2->energy()+lepton3->energy();
	    pxtot = lepton2->px()+lepton3->px();
	    pytot = lepton2->py()+lepton3->py();
	    pztot = lepton2->pz()+lepton3->pz();
	    inmass = entot*entot - pxtot*pxtot - pytot*pytot - pztot*pztot;
	    inmass_sqrt=sqrt(inmass);                
	  }

	  if(inmass_sqrt<invMassMin_)*pOut=false;
	  else selectedEvents[5]++;    
	}
        if(*pOut==true) selectedEvents[6]++;
        iEvent.put(pOut);       
}


void WWPreselectionMarker::beginJob(const edm::EventSetup&)
{
	for(int i=0;i<7;i++)
	  selectedEvents[i]=0;	
}


void WWPreselectionMarker::endJob() {
  std::string cutnames[7];
  cutnames[0]="# # tot. events          : ";
  cutnames[1]="# of Leptons > 2         : ";
  cutnames[2]="one lepton with pt>15GeV : ";
  cutnames[3]="ZLepton opposite charge  : ";
  cutnames[4]="MET     > 10             : ";
  cutnames[5]="Inv Mass  > 10           : ";
  cutnames[6]="All cuts                 : ";
  std::cerr<<" =================================== "<<std::endl;
  for(int i=0; i<7; ++i) 
	std::cerr<<" #"<<i <<" "<< cutnames[i] << selectedEvents[i]<<std::endl;
  std::cerr<<" =================================== "<<std::endl;}

//define this as a plug-in
DEFINE_FWK_MODULE(WWPreselectionMarker);

