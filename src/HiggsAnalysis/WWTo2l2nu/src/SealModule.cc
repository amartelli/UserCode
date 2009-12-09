#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/eventsetupdata_registration_macro.h"
#include "HiggsAnalysis/WWTo2l2nu/interface/HtoWWEleCloner.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE (HtoWWEleCloner) ;
