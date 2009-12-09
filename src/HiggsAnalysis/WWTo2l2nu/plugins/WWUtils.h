#ifndef WWUtils_h
#define WWUtils_h
// -*- C++ -*-
//
// Package:    VBFProcessFilter
// Class:      VBFProcessFilter
// 
/* 
 
 Description: filter events based on the Pythia ProcessID and the Pt_hat
 Implementation: inherits from generic EDFilter
 
 */
//
// $Id: WWUtils.h,v 1.1 2009/03/21 11:03:24 amartell Exp $
//
//
// system include files
#include <memory>

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "TLorentzVector.h"

typedef reco::CaloJetCollection::const_iterator WWjetIt ;

void 
setMomentum (TLorentzVector & myvector, 
             const reco::Candidate & gen) ;

std::pair<WWjetIt,WWjetIt>	
findTagJets (WWjetIt begin, WWjetIt end,
             double jetPtMin, double jetEtaMax) ;

double deltaPhi (double phi1, double phi2) ;

#endif

