// -*- C++ -*-
//
// Package:    WWKFactorProducer
// Class:      WWKFactorProducer
//
/**\class WWKFactorProducer WWKFactorProducer.cc HiggsAnalysis/WWKFactorProducer/src/WWKFactorProducer.cc
 
Description: <one line class summary>
 
Implementation:
<Notes on implementation>
*/
//
// Original Author:  Joanna Weng
//         Created:  Fri Feb  1 15:30:42 CET 2008
// $Id: WWKFactorProducer.cc,v 1.1 2009/04/08 12:46:34 amartell Exp $
//
//

// system include files
#include <memory>
#include <iostream>

// user include files
#include "HiggsAnalysis/WWTo2l2nu/plugins/WWKFactorProducer.h"


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "HepMC/WeightContainer.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "TH1D.h"
#include "TFile.h"



WWKFactorProducer::WWKFactorProducer(const edm::ParameterSet& iConfig): pt_histo_(0),debug_(0)
{
  produces<double>();
  inputFilename_=iConfig.getUntrackedParameter<std::string>("inputFilename","dummy.root");  
  processID_ = iConfig.getUntrackedParameter<int>("ProcessID",0);  
  debug_ = iConfig.getUntrackedParameter<bool>("Debug",0);
  useNNLO_ = iConfig.getUntrackedParameter<bool>("UseNNLO",false);
  edm::FileInPath path_inputFilename(inputFilename_.c_str());
  pt_histo_ = new  WWKfactorList("KFactorList", path_inputFilename.fullPath().c_str() );

  // add possibility of alternative processes (VBF)
  std::vector<int> defpid1;
  defpid1.push_back(-1);
  altProcessID_ = iConfig.getUntrackedParameter<std::vector<int> >("AltProcessID",defpid1);
}


WWKFactorProducer::~WWKFactorProducer()
{
}

// ------------ method called to for each event  ------------
void WWKFactorProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)

{
  using namespace edm;
  // get HepMC::GenEvent ...
  Handle<HepMCProduct> evt_h;
  iEvent.getByType(evt_h);
  
  // get process id;  Does not work -> Ask Guillelmo?
  // Handle< int > genProcessID;
  //  iEvent.getByLabel( "genEventProcID", genProcessID );
  // int ThisprocessID = *genProcessID;
   // if ( ThisprocessID ==  processID){

  // Event weight to be stored in event
  std::auto_ptr<double> pweight(new double);
  *pweight=1.;

  HepMC::GenEvent* evt = new  HepMC::GenEvent(*(evt_h->GetEvent()));   
  if (debug_)std::cout << " Process Id " << evt->signal_process_id()   << std::endl; 
  // Gluon Fusion found
  if ( processID_ == evt->signal_process_id() ){
    // look for Higgs Boson and determine its pt
    std::vector<HepMC::GenParticle*> higgs;
    for(HepMC::GenEvent::particle_iterator it = evt->particles_begin(); it != evt->particles_end(); ++it){
      if(abs((*it)->pdg_id())==25)
	higgs.push_back(*it);
    }
    //Print the complete table
    //if (debug_) std::cout << (*pt_histo_);
    math::XYZTLorentzVector tot_momentum(higgs[0]->momentum());
    //calculate bin size
    double binsize = pt_histo_->GetXaxis()->GetXmax() /pt_histo_->GetNbinsX();
    double higgspt = tot_momentum.pt();
    // which bin ?
    int bin=int ((higgspt /binsize)) + 1 ;
    // overflow protection: use last entry
    if(bin>pt_histo_->GetNbinsX()) bin=pt_histo_->GetNbinsX();
    if (debug_){
    std::cout <<" Bin Size "<< binsize <<std::endl;
    std::cout <<" Higgs Pt "<< higgspt <<std::endl;
    std::cout <<" Bin  "<< bin <<std::endl;
    std::cout <<" KFactor "<<   pt_histo_->GetBinContent(bin) <<std::endl;
    }
    // get KFactor
    *pweight=  pt_histo_->GetBinContent(bin);
  } else {
    bool isAlter=false;
    for(unsigned int j=0; j<altProcessID_.size(); ++j)
      if(evt->signal_process_id()==altProcessID_[j]) isAlter=true;
    if(isAlter)  {
      if(!useNNLO_)
	*pweight=  pt_histo_->GetAlterKfactor();
      else
	*pweight=  pt_histo_->GetAlterNNLOKfactor();
      if(debug_) std::cout <<" KFactor "<< *pweight <<std::endl;
    }
  }
  iEvent.put(pweight);
  delete evt;  
}

// ------------ method called once each job just before starting event loop  ------------
void
WWKFactorProducer::beginJob(const edm::EventSetup&)
{}

// ------------ method called once each job just after ending the event loop  ------------
void
WWKFactorProducer::endJob()
{
  if (pt_histo_) delete(pt_histo_);
}


//define this as a plug-in
DEFINE_FWK_MODULE(WWKFactorProducer);

