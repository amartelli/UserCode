#ifndef RunInfoTreeFiller_h
#define RunInfoTreeFiller_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "HiggsAnalysis/WWTo2l2nu/interface/TreeContent.h"


class RunInfoTreeFiller {

public:
  
  RunInfoTreeFiller(const edm::Event&);
  virtual ~RunInfoTreeFiller();

  void RunFillTree(TreeContent& myTreeVariables_);

protected:

  const edm::Event& iEvent_;

};

#endif // RunInfoTreeFiller_h
