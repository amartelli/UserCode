// my includes
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackExtraBase.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "HiggsAnalysis/WWTo2l2nu/plugins/WWEleAmbiguityResolve.h"

//CLHEP
#include <CLHEP/Vector/LorentzVector.h>

#include <iostream>

WWEleAmbiguityResolve::WWEleAmbiguityResolve (const edm::ParameterSet& conf)
{
  m_reducedElectronsRefCollectionLabel = conf.getParameter<edm::InputTag>("reducedElectronsRefCollectionLabel");
  m_doRefCheck                         = conf.getParameter<bool>("doRefCheck");

  //  std::cout << " WWEleAmbiguityResolve::WWEleAmbiguityResolve " << std::endl;
}


// ------------------------------------------------------------


void 
WWEleAmbiguityResolve::select (edm::Handle<collection> input, 
				const edm::Event& evt,  const edm::EventSetup& evtSetup ) 
{

  //  std::cout << " WWEleAmbiguityResolve::select "<< std::endl;
  m_selected.clear() ;
  ambEle.clear();

  if ( m_doRefCheck ) 
    evt.getByLabel(m_reducedElectronsRefCollectionLabel, m_reduced);

  edm::LogInfo("WWEleAmbiguityResolve") << "WWEleAmbiguityResolve starting, original collection size = " 
					 << input->size();
  Init(input);
  
  ResolveByEoverP(input);

  edm::LogInfo("WWEleAmbiguityResolve") << "WWEleAmbiguityResolve ending, final collection size = " 
					 << m_selected.size();

  return ;
}       


// ------------------------------------------------------------


WWEleAmbiguityResolve::~WWEleAmbiguityResolve()
{}


void WWEleAmbiguityResolve::Init(const edm::Handle<collection> & input) {

  //  std::cout << " WWEleAmbiguityResolve::Init "<< std::endl;

  // This method is doing the check for ambiguity, where ambiguity is
  // 2 electron candidates sharing the same cluster or the same track.
  // This method fills : 
  // --- resolvedEle collections with the no-ambigus electrons
  
  int no_ambiguity;
  std::vector<unsigned int> amb_index;

  for (unsigned int ii=0; ii < input->size(); ii++){  

    no_ambiguity=0;
 
    // first electron
    const reco::GsfElectron & thisEle1 = (*input)[ii];
    reco::GsfTrackRef thisTrack1  = (thisEle1).gsfTrack();
    reco::SuperClusterRef thisSc1 = (thisEle1).superCluster();

    // 2nd loop over all the electron candidates
    for (unsigned int jj=ii+1; jj<input->size(); jj++){

      // second electron
      const reco::GsfElectron & thisEle2 = (*input)[jj];
      reco::GsfTrackRef thisTrack2  = (thisEle2).gsfTrack();
      reco::SuperClusterRef thisSc2 = (thisEle2).superCluster();

      // sharing sc or track
      if ( (thisSc1 == thisSc2) || (thisTrack1 == thisTrack2) ){      
        amb_index.push_back(jj);
        ambEle.push_back(std::pair<unsigned int, unsigned int>(ii, jj));
      }
      else { no_ambiguity++; }
    }
    
    bool test=true;
    for (unsigned int iii=0; iii<amb_index.size(); iii++){
      if (amb_index[iii] == ii ) test=false;
    }
    if ((unsigned int)no_ambiguity == (unsigned int)input->size()-1-ii && test==true) {
      reco::GsfElectronRef electronRef(input, ii);
      edm::LogInfo("WWEleAmbiguityResolve") << "easy, not ambiguous at all. SuperClusterRef->energy() = " 
					     << electronRef->superCluster()->energy();
      bool presentInReducedCollection=true;
      if ( m_doRefCheck ) {
	if (find(m_reduced->begin(), m_reduced->end(), electronRef) == m_reduced->end())
	  presentInReducedCollection=false;
      }
      
      if ( presentInReducedCollection )
	m_selected.push_back(electronRef);
    }
  }
}


void WWEleAmbiguityResolve::ResolveByEoverP(const edm::Handle<collection> & input)
{
  //  std::cout << " WWEleAmbiguityResolve::ResolveByEoverP "<< std::endl;

  // This method resolves the ambiguities using E(supercluster) over P , 
  // the closest to one is declared as the good candidate.
  // It uses the map associating the ambigus electron candidates
  // It continues filling:
  // --- resolvedEle collections with the chosen electrons candidates 

  if ( input->size() !=0 ) {

    std::vector<std::pair<unsigned int, unsigned int> >::iterator it;

    for ( it=ambEle.begin(); it<ambEle.end(); it++){
      
      // loop over electrons which have one or more than one ambiguity
      int bestEleId=it->first;
      int multiAmbEleId=(int)it->first;
      while((int)it->first==multiAmbEleId && it<ambEle.end()) {
	//	std::cout << "multiAmbEleId = " << multiAmbEleId << "\tbestEleId = " << bestEleId << std::endl;
	const reco::GsfElectron & bestEle = (*input)[bestEleId];
	const reco::GsfElectron & compEle = (*input)[it->second];
	edm::LogInfo("WWEleAmbiguityResolve") << "compEle.superCluster()->energy() = " 
					       << compEle.superCluster()->energy() 
					       << "\tbestEle.superCluster()->energy() = " 
					       << bestEle.superCluster()->energy();
	if(fabs(compEle.eSuperClusterOverP()-1) <= fabs(bestEle.eSuperClusterOverP()-1)) bestEleId=it->second;
	it++;
      }
      edm::LogInfo("WWEleAmbiguityResolve") << "looping on the ambiguous electrons: best ele index = " 
					     <<  bestEleId << std::endl;
      reco::GsfElectronRef electronRef(input, bestEleId);
      
      bool presentInReducedCollection=true;
      if ( m_doRefCheck ) {
	if (find(m_reduced->begin(), m_reduced->end(), electronRef) == m_reduced->end())
	  presentInReducedCollection=false;
      }
      
      if ( presentInReducedCollection )
	m_selected.push_back(electronRef);
    }
    
  }

}


