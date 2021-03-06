#ifndef WWEleAmbiguityResolve_h
#define WWEleAmbiguityResolve_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <functional>
#include <vector>
#include <map>


class WWEleAmbiguityResolve{

 public:

  typedef const reco::GsfElectron * electron ;
  typedef reco::GsfElectronCollection collection ;
  typedef reco::GsfElectronRefVector container;
  typedef container::const_iterator const_iterator ;
  
  //! ctor
  WWEleAmbiguityResolve (const edm::ParameterSet& conf) ;
  //!dtor
  ~WWEleAmbiguityResolve () ;

  //! iterator to the begin of the selected collection
  const_iterator begin () const { return m_selected.begin () ; }
  
  //! iterator to the end of the selected collection
  const_iterator end () const { return m_selected.end () ; }

  //! the actual selector method 
  void select (edm::Handle<collection>, const edm::Event&, const edm::EventSetup& ) ;
     

 private:

  //! find the non ambiguous electrons, initialise the ambiguity map
  void Init(const edm::Handle<collection> & input);
  
  //! ambiguity resolution
  void ResolveByEoverP(const edm::Handle<collection> & input);

  //! the selected collection
  container m_selected ;
 
  //! map between the ambiguous electrons
  std::vector<std::pair<unsigned int, unsigned int> > ambEle;

  //! if doRefCheck, select only among the elements present in the reducedElectron collection
  bool m_doRefCheck;
  edm::InputTag m_reducedElectronsRefCollectionLabel;
  edm::Handle< edm::RefVector<collection> > m_reduced;

};  

#endif

