#ifndef WWTreeProducer_h
#define WWTreeProducer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "HiggsAnalysis/WWTo2l2nu/interface/TreeContent.h"
#include <TFile.h>

class WWTreeProducer : public edm::EDAnalyzer {
 public:
  explicit WWTreeProducer(const edm::ParameterSet&);
  ~WWTreeProducer();


 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void MCInfo(const edm::Event&, const edm::EventSetup&, bool& eventOk);

private:
  
  int is3e_;
  int is3mu_;
  int is2e1mu_;
  int is2mu1e_;
  int nEvents_;

  //! name 
  std::string nameFile_;
  std::string nameTree_;

  bool dumpRunInfo_;
  bool dumpPreselInfo_;

  bool dumpElectrons_;
  bool dumpSCs_;
  bool dumpMuons_;
  bool dumpTracks_;
  bool dumpVertices_;
  bool dumpMet_, dumpGenMet_;

  std::string mctype_;
  edm::InputTag electronCollection_, muonCollection_;
  edm::InputTag metCollection_, genMetCollection_;
  edm::InputTag vertexCollection_;
  edm::InputTag trackCollection_;
  edm::InputTag ecalBarrelSCCollection_, ecalEndcapSCCollection_;
  edm::InputTag mcTruthCollection_;

  int jevt_;


  //! ROOT file with the plain ROOT tree inside
  TFile* fileOut_;


  TreeContent myTreeVariables_ ;
  TTree* tree_ ;


};
#endif // WWTreeProducer_h

