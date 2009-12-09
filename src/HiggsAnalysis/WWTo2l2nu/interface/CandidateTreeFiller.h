
#ifndef CandidateTreeFiller_h
#define CandidateTreeFiller_h

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"

#include "HiggsAnalysis/WWTo2l2nu/interface/TreeContent.h"
#include <TTree.h>

using namespace cms;
using namespace edm;
using namespace reco;



class CandidateTreeFiller {

 public:

  //! Dump everything                                                                                                                                
  CandidateTreeFiller(edm::InputTag recoCollectionTag,
                      const edm::Event& iEvent, const edm::EventSetup& iSetup);

  //! Destructor                                                                                                                                     
  virtual ~CandidateTreeFiller();

  // --------------------------------------------- //

  virtual void writeCollectionToTree(TreeContent& myTreeVariables, std::string name);

  virtual void setTreeBranches(std::vector<int>* charge, std::vector<float>* energy, std::vector<float>* et, 
			       std::vector<float>* momentum, std::vector<float>* vertexX,  
			       std::vector<float>* vertexY, std::vector<float>* vertexZ, 
			       std::vector<float>* theta, std::vector<float>* eta, std::vector<float>* phi, 
			       std::vector<float>* x, std::vector<float>* y, std::vector<float>* z, 
			       std::vector<float>* mass, std::vector<float>* mt, 
			       std::vector<int>* pdgId, std::vector<int>* nDau);
  

   

 protected:

  virtual void Initialise();
  
  virtual void writeCandInfo(const Candidate *cand, int ncand);
  
  //  virtual void writeMcMatchInfo(const edm::View<reco::Candidate> *recoCollection,
  //				const edm::View<reco::Candidate> *genCollection,
  //			        std::vector<int>* mcIndex );


  // Friends
  // std::vector< const edm::View<reco::Candidate>* > daugCollectionList_;

  edm::InputTag recoCollectionTag_;
  const edm::Event& iEvent_;
  const edm::EventSetup& iSetup_;
  edm::InputTag genCollectionTag_;
  
  //  edm::InputTag matchMap_;
  
  //  bool saveCand_;
  //  bool doMcMatch_;
  
  const Candidate* cand_;


  // Tree Variables

  int charge_;
  float energy_, pt_, et_, momentum_, momentumX_, momentumY_, momentumZ_;
  float vertexX_, vertexY_, vertexZ_;
  float theta_, eta_, phi_;
  float x_, y_, z_;
  float mass_, mt_;
  int pdgId_;
  int nDau_;

  int ncand_;
  
};

#endif // CandidateTreeFiller_h
