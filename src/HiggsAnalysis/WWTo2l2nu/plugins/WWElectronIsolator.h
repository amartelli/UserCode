#ifndef WWELECTRONISOLATOR
#define WWELECTRONISOLATOR

#include <memory>
#include <vector>
#include <math.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"

class WWElectronIsolator{
 public:
   explicit WWElectronIsolator(const edm::ParameterSet&);
   ~WWElectronIsolator();

   typedef reco::GsfElectronCollection collection;
   typedef std::vector<reco::GsfElectronRef> ::const_iterator const_iterator;

   const_iterator begin () const { return selected_.begin () ; }
   const_iterator end () const { return  selected_.end () ; }

   void select (edm::Handle<reco::GsfElectronCollection>,
                const edm::Event&, 
                const edm::EventSetup&) ;
	
 private:	
   std::vector<reco::GsfElectronRef> selected_;
   double theTrackIsolCut_; 
   double theAbsTrackIsolCut_; 
   double theCaloIsolCut_; 
   double theECALIsolCut_; 
   //   bool absolute_;

   TH1F* m_SumPt_over_Pt_EleTk;
   TH1F* m_SumPt_EleTk;
   TH1F* m_SumPt_EleCalo;
   TH1F* m_SumPt_EleEcal;

};


#endif
