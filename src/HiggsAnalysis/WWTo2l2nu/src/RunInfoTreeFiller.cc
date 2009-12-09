#include "HiggsAnalysis/WWTo2l2nu/interface/RunInfoTreeFiller.h"

using namespace edm;
using namespace std;

RunInfoTreeFiller::RunInfoTreeFiller(const edm::Event& iEvent):
  iEvent_(iEvent) {}

RunInfoTreeFiller::~RunInfoTreeFiller() {}

void RunInfoTreeFiller::RunFillTree(TreeContent& myTreeVariables_)
{ 
  myTreeVariables_.runN = iEvent_.id().run();
  myTreeVariables_.evtN = iEvent_.id().event();
}

