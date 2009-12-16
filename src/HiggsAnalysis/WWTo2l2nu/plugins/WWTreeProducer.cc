// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Run.h"

#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "HiggsAnalysis/WWTo2l2nu/interface/TreeContent.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "HiggsAnalysis/WWTo2l2nu/interface/RunInfoTreeFiller.h"
#include "HiggsAnalysis/WWTo2l2nu/interface/ElectronTreeFiller.h"
#include "HiggsAnalysis/WWTo2l2nu/interface/SuperClusterTreeFiller.h"
#include "HiggsAnalysis/WWTo2l2nu/interface/MuonTreeFiller.h"
#include "HiggsAnalysis/WWTo2l2nu/interface/TrackTreeFiller.h"
#include "HiggsAnalysis/WWTo2l2nu/interface/VertexTreeFiller.h"
#include "HiggsAnalysis/WWTo2l2nu/interface/CandidateTreeFiller.h"

#include "HiggsAnalysis/WWTo2l2nu/interface/McSelectWW.h"

#include "HiggsAnalysis/WWTo2l2nu/plugins/WWTreeProducer.h"



WWTreeProducer::WWTreeProducer(const edm::ParameterSet& iConfig)
{
  nameFile_      = iConfig.getUntrackedParameter<std::string>("nameFile", "RootOutput.root");
  nameTree_      = iConfig.getUntrackedParameter<std::string>("nameTree", "BaseTree");
  
  // Candidate Collections
  dumpElectrons_      = iConfig.getUntrackedParameter<bool>("dumpElectrons", false);
  dumpSCs_            = iConfig.getUntrackedParameter<bool>("dumpSCs", false);
  dumpMuons_          = iConfig.getUntrackedParameter<bool>("dumpMuons", false);
  dumpTracks_         = iConfig.getUntrackedParameter<bool>("dumpTracks", false);
  dumpVertices_       = iConfig.getUntrackedParameter<bool>("dumpVertices", false);
  dumpMet_            = iConfig.getUntrackedParameter<bool>("dumpMet", false);
  dumpGenMet_         = iConfig.getUntrackedParameter<bool>("dumpGenMet", false);

  // data run informations
  dumpRunInfo_ = iConfig.getUntrackedParameter<bool>("dumpRunInfo",false);
  dumpPreselInfo_     = iConfig.getUntrackedParameter<bool>("dumpPreselInfo", false);

  mctype_                  = iConfig.getUntrackedParameter<std::string>("mctype", "bau");
  electronCollection_      = iConfig.getParameter<edm::InputTag>("electronCollection");
  ecalBarrelSCCollection_  = iConfig.getParameter<edm::InputTag>("ecalBarrelSCCollection");
  ecalEndcapSCCollection_  = iConfig.getParameter<edm::InputTag>("ecalEndcapSCCollection");
  muonCollection_          = iConfig.getParameter<edm::InputTag>("muonCollection");
  trackCollection_         = iConfig.getParameter<edm::InputTag>("trackCollection");
  vertexCollection_        = iConfig.getParameter<edm::InputTag>("vertexCollection");
  metCollection_           = iConfig.getParameter<edm::InputTag>("metCollection");
  genMetCollection_        = iConfig.getParameter<edm::InputTag>("genMetCollection");
  mcTruthCollection_       = iConfig.getParameter<edm::InputTag>("mcTruthCollection");
}


WWTreeProducer::~WWTreeProducer() { }



// ------------ method called to for each event  ------------
void WWTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 if(nEvents_ % 10000 == 0) std::cerr << ">>>>>> WWTreeProducer::Analyze >>>>>>  " << nEvents_ << std::endl;
  ++nEvents_;

  using namespace edm;



  // Fill Tree
  initializeBranches(tree_, myTreeVariables_);

  bool eventOk = false;

  MCInfo(iEvent, iSetup, eventOk);
  if(eventOk == false) return; 

  if(dumpRunInfo_){
    /// fill the run info (run number, event, ...)
    RunInfoTreeFiller runTree (iEvent);
    runTree.RunFillTree(myTreeVariables_);
  }


//   //  std::cout << "Fill the Preselected "<<std::endl;
//   // fill preselection output
//   if (dumpPreselInfo_) {
//     Handle<bool> selected;
//     iEvent.getByLabel("preselectionMarker", selected );
//     bool isSelected = *selected;
//     myTreeVariables_.evtPresel = isSelected;
//   }
  

  //  std::cout << "Fill the Electrons "<<std::endl;
  // fill Electrons block
  if(dumpElectrons_) {
    ElectronTreeFiller electronTree(electronCollection_, iEvent, iSetup);     
    electronTree.writeCollectionToTree(myTreeVariables_);
  }


  //  std::cout << "Fill the SC "<<std::endl;
  // fill SC block
  if (dumpSCs_) {
    SuperClusterTreeFiller treeFillBarrel(ecalBarrelSCCollection_, iEvent, iSetup, 100);
    treeFillBarrel.SClusterFillTree(myTreeVariables_);

    SuperClusterTreeFiller treeFillEndcap(ecalEndcapSCCollection_, iEvent, iSetup, 100);
    treeFillEndcap.SClusterFillTree(myTreeVariables_);
  }

  //  fill muons block
  if(dumpMuons_) {
    MuonTreeFiller muonTree(muonCollection_, iEvent, iSetup);
    muonTree.writeCollectionToTree(myTreeVariables_);

  }



   //  std::cout << "Fill the Tracks "<<std::endl;
   // fill track block
   if(dumpTracks_) {
     TrackTreeFiller treeFiller(vertexCollection_, trackCollection_, iEvent, iSetup);
     treeFiller.TrackFillTree(myTreeVariables_);
   }


  //  std::cout << "Fill the PVertex "<<std::endl;
  //fill Primary Vertex and associated tracks
  if(dumpVertices_){
    VertexTreeFiller treeFill(vertexCollection_, iEvent, iSetup);
    treeFill.VertexFillTree(myTreeVariables_);
  }
  
  
  //  std::cout << "Fill the MET "<<std::endl;
  // fill MET block
  if(dumpMet_) {
    CandidateTreeFiller treeRecoFill(metCollection_, iEvent, iSetup);
    
    std::string localname = "met";
    treeRecoFill.writeCollectionToTree(myTreeVariables_, localname);

  }
  
  //  std::cout << "Fill the MCMET "<<std::endl;
  // dump generated MET
  if(dumpGenMet_) {
    CandidateTreeFiller treeGenFill(genMetCollection_, iEvent, iSetup);
    
    std::string localname = "genmet";
    treeGenFill.writeCollectionToTree(myTreeVariables_, localname);
    }
    
  
   //  std::cout << " Filling the Tree " << std::endl;
  tree_->Fill();
  //  std::cout << " Tree fillato " << std::endl;

    
}



// ------------ method called once each job just before starting event loop  ------------
void WWTreeProducer::beginJob(const edm::EventSetup&) {

  nEvents_ = 0;
  is3e_ = 0;
  is3mu_ = 0;
  is2e1mu_ = 0;
  is2mu1e_ = 0;

  // Create File
  fileOut_ = TFile::Open(nameFile_.c_str(), "RECREATE");

  // Initialize Tree
  //  tree_ = new TTree ( "WZAnalysis","WZAnalysis" );
  tree_ = new TTree ( nameTree_.c_str(), nameTree_.c_str() );
  setBranches (tree_, myTreeVariables_);

  jevt_ = 1;
}



// ------------ method called once each job just after ending the event loop  ------------
void  WWTreeProducer::endJob() {
  fileOut_->cd();
  tree_->Write () ;
  fileOut_->Close();
  std::cout << " Total of 3e events = " << is3e_ << std::endl; 
  std::cout << " Total of 2e1mu events = " << is2e1mu_ << std::endl; 
  std::cout << " Total of 3mu events = " << is3mu_ << std::endl; 
  std::cout << " Total of 2mu1e events = " << is2mu1e_ << std::endl; 

}


// ------------- MCInfo

void WWTreeProducer::MCInfo(const edm::Event& e, const edm::EventSetup& iSetup, bool& eventOk)
{
  // *** monte carlo info ***
  bool is3e = false;
  char mc[200];
  strcpy(mc, mctype_.c_str());
  //compare two strings up to 6 charact
  if (strncmp(mc, "signal", 6) == 0) {
    McSelectWW selectorMcWW(e, mcTruthCollection_);
    std::string eventId_local = "signal";
    myTreeVariables_.mcEventId->push_back(eventId_local);

    if(selectorMcWW.isElectrons3()) {
      std::string WWevent_local = "3Ele";
      myTreeVariables_.mcWWevent->push_back(WWevent_local);
      selectorMcWW.SelectFillTree(myTreeVariables_);
      is3e = true;
      ++is3e_;
      eventOk = true;
    }

    else if(selectorMcWW.isElectrons2Muon1()) {
      std::string WWevent_local = "2Ele1Mu";
      myTreeVariables_.mcWWevent->push_back(WWevent_local);    
      selectorMcWW.SelectFillTree(myTreeVariables_);
      ++is2e1mu_;
      eventOk = true;
    }

    else if(selectorMcWW.isMuons2Electron1()) {
      std::string WWevent_local = "2Mu1Ele";
      myTreeVariables_.mcWWevent->push_back(WWevent_local);     
      selectorMcWW.SelectFillTree(myTreeVariables_);
      ++is2mu1e_;
      eventOk = true;
    }

    else if(selectorMcWW.isMuons3()) {
      std::string WWevent_local = "3Mu";
      myTreeVariables_.mcWWevent->push_back(WWevent_local);     
      selectorMcWW.SelectFillTree(myTreeVariables_);
      ++is3mu_;
      eventOk = true;
    }
    else {
      std::string WWevent_local = "Other";
      myTreeVariables_.mcWWevent->push_back(WWevent_local);
      //      selectorMcWW.SelectFillTree(myTreeVariables_);
      eventOk = true;
    }
  }
  else {
    eventOk = true;
    std::string WWevent_local = "BG";
    myTreeVariables_.mcWWevent->push_back(WWevent_local);
  }
}

