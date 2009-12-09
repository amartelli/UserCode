#ifndef WWELECTRONISOLATOR
#define WWELECTRONISOLATOR

#include <memory>
#include <vector>
#include <math.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"


class WWElectronIsolator{
 public:
   explicit WWElectronIsolator(const edm::ParameterSet&);
   ~WWElectronIsolator();

   // Collections to be selected
   //   typedef reco::PixelMatchGsfElectronCollection collection;
   typedef reco::GsfElectronCollection collection;
   //   typedef std::vector<reco::PixelMatchGsfElectronRef> ::const_iterator const_iterator;
   typedef std::vector<reco::GsfElectronRef> ::const_iterator const_iterator;

   //define iterators with above typedef
   const_iterator begin () const { return selected_.begin () ; }
   const_iterator end () const { return  selected_.end () ; }

   //   void select (edm::Handle<reco::PixelMatchGsfElectronCollection>,
   void select (edm::Handle<reco::GsfElectronCollection>,
                const edm::Event&, 
                const edm::EventSetup&) ;
	
 private:	
   //   std::vector<reco::PixelMatchGsfElectronRef> selected_;
   std::vector<reco::GsfElectronRef> selected_;
   edm::InputTag tracksLabel_;
   edm::InputTag electronsLabel_;
   edm::InputTag trackIsolationProducer_;
   bool doRefCheck_;
   edm::InputTag selectedElectronsRefLabel_;
   double theTrackIsolCut_; 
};


#endif
