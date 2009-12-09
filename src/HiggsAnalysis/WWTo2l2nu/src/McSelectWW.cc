#include "HiggsAnalysis/WWTo2l2nu/interface/McSelectWW.h"

#include <map>
#include <set>
#include <vector>

using namespace std;

McSelectWW::McSelectWW(const edm::Event& e, edm::InputTag mcTruthCollection_)
{
  validity_ = false;
  Electrons3_ = false;
  Electrons2Muon1_ = false;
  Muons3_ = false;
  Muons2Electron1_ = false;
  
  e.getByLabel(mcTruthCollection_, genCandidates_);
  if(! (genCandidates_.isValid ()) )
    {
      std::cerr << "McSelectWW::Warning: " << genCandidates_ << " not available" << std::endl;
      return;
    }
  
  Init();
  //  InitAfterPhotos();
}

void McSelectWW::Init() {
  
  // 1. make a loop untill both W and Z are found
  // 2. then follow the decay tree
  
  vector<const reco::Candidate*> tmpW;
  vector<const reco::Candidate*> tmpZ;
  vector<const reco::Candidate*> Wdecays;
  vector<const reco::Candidate*> Zdecays;


  reco::GenParticleCollection::const_iterator p = genCandidates_->begin();

  while (p != genCandidates_->end()){
    //    if (p->status() == 3 ) Wdecays.push_back(&(*p));
    if (abs(p->pdgId()) == 24 && p->status() == 3 ) tmpW.push_back(&(*p));
    if (abs(p->pdgId()) == 23 && p->status() == 3 ) tmpZ.push_back(&(*p));
    //    if (tmpW.size() == 2 ) break;
    p++;
  }
  
  if (tmpW.size() != 1 || tmpZ.size() != 1 ) {
    std::cout << "McSelectWW::Init(): no WZ bosons in the event!" << endl;
    return;
  } 
  else 
    {
      validity_ = true;        
      
      realW_ = tmpW.at(0);
      realZ_ = tmpZ.at(0);
    
      vector<const reco::Candidate*> da_realW, da_realZ;  

      for (unsigned int bb = 0; bb < realW_->numberOfDaughters(); ++bb){
	if (realW_->daughter(bb)->status() !=2 )  da_realW.push_back(realW_->daughter(bb));
	// != 2 avoids gamma brehm
	//	if (realW_->daughter(bb)->status() == 3)  da_realW.push_back(realW_->daughter(bb));
      }
      for (unsigned int bb = 0; bb < realZ_->numberOfDaughters(); ++bb){
	if (realZ_->daughter(bb)->status() !=2 )  da_realZ.push_back(realZ_->daughter(bb));
	// ???
	//	if (virtualW_->daughter(bb)->status() == 3)  da_virtualW.push_back(virtualW_->daughter(bb));
      }
      if ( (da_realW.size() != 2) || (da_realZ.size() != 2) ){
	std::cout <<"McSelectWW::Init(): not a WZ -> 4 body decay!"<<endl;
	return;
      }


      WZSystem_ = tmpW[0]->p4() + tmpZ[0]->p4();
      Wdecays.push_back(da_realW[0]);
      Wdecays.push_back(da_realW[1]);
      Zdecays.push_back(da_realZ[0]);
      Zdecays.push_back(da_realZ[1]);
      
   

      int da1_rW_pdgID = da_realW.at(0) ->pdgId();
      int da2_rW_pdgID = da_realW.at(1) ->pdgId();
      int da1_rZ_pdgID = da_realZ.at(0) ->pdgId();
      int da2_rZ_pdgID = da_realZ.at(1) ->pdgId();
    
      // e nu e+ e-    
      if ( (abs(da1_rW_pdgID) == 11) && (abs(da2_rW_pdgID) == 12) && (da1_rZ_pdgID == -11) && (da2_rZ_pdgID == 11) ) {
	Electrons3_ = true;
	leptonFromW_ = da_realW.at(0); 
	nuFromW_ = da_realW.at(1); 
	leptonFromZ_ = da_realZ.at(1); 
	cleptonFromZ_ = da_realZ.at(0);
 
	da1_rW_pdgID_ = da_realW.at(0) ->pdgId();
	da2_rW_pdgID_ = da_realW.at(1) ->pdgId();
	da1_rZ_pdgID_ = da_realZ.at(1) ->pdgId();
	da2_rZ_pdgID_ = da_realZ.at(0) ->pdgId();
      }

      // nu e e+ e-
      else if ( (abs(da1_rW_pdgID) == 12) && (abs(da2_rW_pdgID) == 11) && (da1_rZ_pdgID == -11) && (da2_rZ_pdgID == 11) ) {
	Electrons3_ = true;
	leptonFromW_ = da_realW.at(1); 
	nuFromW_ = da_realW.at(0); 
	leptonFromZ_ = da_realZ.at(1); 
	cleptonFromZ_ = da_realZ.at(0); 

	da1_rW_pdgID_ = da_realW.at(1) ->pdgId();
	da2_rW_pdgID_ = da_realW.at(0) ->pdgId();
	da1_rZ_pdgID_ = da_realZ.at(1) ->pdgId();
	da2_rZ_pdgID_ = da_realZ.at(0) ->pdgId();
      }


      // e nu e- e+    
      else if ( (abs(da1_rW_pdgID) == 11) && (abs(da2_rW_pdgID) == 12) && (da1_rZ_pdgID == 11) && (da2_rZ_pdgID == -11) ) {
	Electrons3_ = true;
	leptonFromW_ = da_realW.at(0); 
	nuFromW_ = da_realW.at(1); 
	leptonFromZ_ = da_realZ.at(0); 
	cleptonFromZ_ = da_realZ.at(1); 

	da1_rW_pdgID_ = da_realW.at(0) ->pdgId();
	da2_rW_pdgID_ = da_realW.at(1) ->pdgId();
	da1_rZ_pdgID_ = da_realZ.at(0) ->pdgId();
	da2_rZ_pdgID_ = da_realZ.at(1) ->pdgId();
      }

      // nu e e- e+
      else if ( (abs(da1_rW_pdgID) == 12) && (abs(da2_rW_pdgID) == 11) && (da1_rZ_pdgID == 11) && (da2_rZ_pdgID == -11) ) {
	Electrons3_ = true;
	leptonFromW_ = da_realW.at(1); 
	nuFromW_ = da_realW.at(0); 
	leptonFromZ_ = da_realZ.at(0); 
	cleptonFromZ_ = da_realZ.at(1); 

	da1_rW_pdgID_ = da_realW.at(1) ->pdgId();
	da2_rW_pdgID_ = da_realW.at(0) ->pdgId();
	da1_rZ_pdgID_ = da_realZ.at(0) ->pdgId();
	da2_rZ_pdgID_ = da_realZ.at(1) ->pdgId();
      }


      //mu nu e- e+
      else if ( (abs(da1_rW_pdgID) == 13) && (abs(da2_rW_pdgID) == 14) && (da1_rZ_pdgID == 11) && (da2_rZ_pdgID == -11) ) {
	Electrons2Muon1_ = true;
	leptonFromW_ = da_realW.at(0); 
	nuFromW_ = da_realW.at(1); 
	leptonFromZ_ = da_realZ.at(0); 
	cleptonFromZ_ = da_realZ.at(1); 

	da1_rW_pdgID_ = da_realW.at(0) ->pdgId();
	da2_rW_pdgID_ = da_realW.at(1) ->pdgId();
	da1_rZ_pdgID_ = da_realZ.at(0) ->pdgId();
	da2_rZ_pdgID_ = da_realZ.at(1) ->pdgId();
      }


      //nu mu e- e+
      else if ( (abs(da1_rW_pdgID) == 14) && (abs(da2_rW_pdgID) == 13) && (da1_rZ_pdgID == 11) && (da2_rZ_pdgID == -11) ) {
	Electrons2Muon1_ = true;
	leptonFromW_ = da_realW.at(1); 
	nuFromW_ = da_realW.at(0); 
	leptonFromZ_ = da_realZ.at(0); 
	cleptonFromZ_ = da_realZ.at(1);
 
	da1_rW_pdgID_ = da_realW.at(1) ->pdgId();
	da2_rW_pdgID_ = da_realW.at(0) ->pdgId();
	da1_rZ_pdgID_ = da_realZ.at(0) ->pdgId();
	da2_rZ_pdgID_ = da_realZ.at(1) ->pdgId();
      }


      //mu nu e+ e-
      else if ( (abs(da1_rW_pdgID) == 13) && (abs(da2_rW_pdgID) == 14) && (da1_rZ_pdgID == -11) && (da2_rZ_pdgID == 11) ) {
	Electrons2Muon1_ = true;
	leptonFromW_ = da_realW.at(0); 
	nuFromW_ = da_realW.at(1); 
	leptonFromZ_ = da_realZ.at(1); 
	cleptonFromZ_ = da_realZ.at(0); 
     
	da1_rW_pdgID_ = da_realW.at(0) ->pdgId();
	da2_rW_pdgID_ = da_realW.at(1) ->pdgId();
	da1_rZ_pdgID_ = da_realZ.at(1) ->pdgId();
	da2_rZ_pdgID_ = da_realZ.at(0) ->pdgId();
      }

      //nu mu e+ e-
      else if ( (abs(da1_rW_pdgID) == 14) && (abs(da2_rW_pdgID) == 13) && (da1_rZ_pdgID == -11) && (da2_rZ_pdgID == 11) ) {
	Electrons2Muon1_ = true;
	leptonFromW_ = da_realW.at(1); 
	nuFromW_ = da_realW.at(0); 
	leptonFromZ_ = da_realZ.at(1); 
	cleptonFromZ_ = da_realZ.at(0); 

	da1_rW_pdgID_ = da_realW.at(1) ->pdgId();
	da2_rW_pdgID_ = da_realW.at(0) ->pdgId();
	da1_rZ_pdgID_ = da_realZ.at(1) ->pdgId();
	da2_rZ_pdgID_ = da_realZ.at(0) ->pdgId();
      }

      /////////

      //mu nu mu+ mu-
      else if ( (abs(da1_rW_pdgID) == 13) && (abs(da2_rW_pdgID) == 14) && (da1_rZ_pdgID == -11) && (da2_rZ_pdgID == 11) ) {
	Muons3_ = true;
	leptonFromW_ = da_realW.at(0); 
	nuFromW_ = da_realW.at(1); 
	leptonFromZ_ = da_realZ.at(1); 
	cleptonFromZ_ = da_realZ.at(0); 

	da1_rW_pdgID_ = da_realW.at(0) ->pdgId();
	da2_rW_pdgID_ = da_realW.at(1) ->pdgId();
	da1_rZ_pdgID_ = da_realZ.at(1) ->pdgId();
	da2_rZ_pdgID_ = da_realZ.at(0) ->pdgId();
      }


      //nu mu mu+ mu-
      else if ( (abs(da1_rW_pdgID) == 14) && (abs(da2_rW_pdgID) == 13) && (da1_rZ_pdgID == -13) && (da2_rZ_pdgID == 13) ) {
	Muons3_ = true;
	leptonFromW_ = da_realW.at(1); 
	nuFromW_ = da_realW.at(0); 
	leptonFromZ_ = da_realZ.at(1); 
	cleptonFromZ_ = da_realZ.at(0); 

	da1_rW_pdgID_ = da_realW.at(1) ->pdgId();
	da2_rW_pdgID_ = da_realW.at(0) ->pdgId();
	da1_rZ_pdgID_ = da_realZ.at(1) ->pdgId();
	da2_rZ_pdgID_ = da_realZ.at(0) ->pdgId();
      }


      //mu nu mu- mu+
      else if ( (abs(da1_rW_pdgID) == 13) && (abs(da2_rW_pdgID) == 14) && (da1_rZ_pdgID == 13) && (da2_rZ_pdgID == -13) ) {
	Muons3_ = true;
	leptonFromW_ = da_realW.at(0); 
	nuFromW_ = da_realW.at(1); 
	leptonFromZ_ = da_realZ.at(0); 
	cleptonFromZ_ = da_realZ.at(1); 

	da1_rW_pdgID_ = da_realW.at(0) ->pdgId();
	da2_rW_pdgID_ = da_realW.at(1) ->pdgId();
	da1_rZ_pdgID_ = da_realZ.at(0) ->pdgId();
	da2_rZ_pdgID_ = da_realZ.at(1) ->pdgId();
      }


      //nu mu mu+ mu-
      else if ( (abs(da1_rW_pdgID) == 14) && (abs(da2_rW_pdgID) == 13) && (da1_rZ_pdgID == -13) && (da2_rZ_pdgID == 13) ) {
	Muons3_ = true;
	leptonFromW_ = da_realW.at(1); 
	nuFromW_ = da_realW.at(0); 
	leptonFromZ_ = da_realZ.at(1); 
	cleptonFromZ_ = da_realZ.at(0); 

	da1_rW_pdgID_ = da_realW.at(1) ->pdgId();
	da2_rW_pdgID_ = da_realW.at(0) ->pdgId();
	da1_rZ_pdgID_ = da_realZ.at(1) ->pdgId();
	da2_rZ_pdgID_ = da_realZ.at(0) ->pdgId();
      }


      //nu e mu+ mu-
      else if ( (abs(da1_rW_pdgID) == 12) && (abs(da2_rW_pdgID) == 11) && (da1_rZ_pdgID == -13) && (da2_rZ_pdgID == 13) ) {
	Muons2Electron1_ = true;
	leptonFromW_ = da_realW.at(1); 
	nuFromW_ = da_realW.at(0); 
	leptonFromZ_ = da_realZ.at(1); 
	cleptonFromZ_ = da_realZ.at(0); 

	da1_rW_pdgID_ = da_realW.at(1) ->pdgId();
	da2_rW_pdgID_ = da_realW.at(0) ->pdgId();
	da1_rZ_pdgID_ = da_realZ.at(1) ->pdgId();
	da2_rZ_pdgID_ = da_realZ.at(0) ->pdgId();
      }

      //nu e mu- mu+
      else if ( (abs(da1_rW_pdgID) == 12) && (abs(da2_rW_pdgID) == 11) && (da1_rZ_pdgID == 13) && (da2_rZ_pdgID == -13) ) {
	Muons2Electron1_ = true;
	leptonFromW_ = da_realW.at(1); 
	nuFromW_ = da_realW.at(0); 
	leptonFromZ_ = da_realZ.at(0); 
	cleptonFromZ_ = da_realZ.at(1); 

	da1_rW_pdgID_ = da_realW.at(1) ->pdgId();
	da2_rW_pdgID_ = da_realW.at(0) ->pdgId();
	da1_rZ_pdgID_ = da_realZ.at(0) ->pdgId();
	da2_rZ_pdgID_ = da_realZ.at(1) ->pdgId();
      }


      //e nu mu+ mu-
      else if ( (abs(da1_rW_pdgID) == 11) && (abs(da2_rW_pdgID) == 12) && (da1_rZ_pdgID == -13) && (da2_rZ_pdgID == 13) ) {
	Muons2Electron1_ = true;
	leptonFromW_ = da_realW.at(0); 
	nuFromW_ = da_realW.at(1); 
	leptonFromZ_ = da_realZ.at(1); 
	cleptonFromZ_ = da_realZ.at(0); 

	da1_rW_pdgID_ = da_realW.at(0) ->pdgId();
	da2_rW_pdgID_ = da_realW.at(1) ->pdgId();
	da1_rZ_pdgID_ = da_realZ.at(1) ->pdgId();
	da2_rZ_pdgID_ = da_realZ.at(0) ->pdgId();
      }

      //e nu mu- mu+
      else if ( (abs(da1_rW_pdgID) == 11) && (abs(da2_rW_pdgID) == 12) && (da1_rZ_pdgID == 13) && (da2_rZ_pdgID == -13) ) {
	Muons2Electron1_ = true;
	leptonFromW_ = da_realW.at(0); 
	nuFromW_ = da_realW.at(1); 
	leptonFromZ_ = da_realZ.at(0); 
	cleptonFromZ_ = da_realZ.at(1); 

	da1_rW_pdgID_ = da_realW.at(0) ->pdgId();
	da2_rW_pdgID_ = da_realW.at(1) ->pdgId();
	da1_rZ_pdgID_ = da_realZ.at(0) ->pdgId();
	da2_rZ_pdgID_ = da_realZ.at(1) ->pdgId();
      }


      else {
	cout << "SelectMcWW::Init(): strange WZ decay  event!" << endl;
	cout << "ID of particles from WZ decay: "<< Wdecays[0]->pdgId()  << " | " << Wdecays[1]->pdgId() 
	     << " | " << Zdecays[0]->pdgId() << " | " <<  Zdecays[1]->pdgId() << endl; 
      }	
      
    } 
}



/*
void McSelectWW::InitAfterPhotos() {

  // start from W and W, make a loop through their dauthers and decide which is electron and which is positron
  if( realW_ && realZ_){
    const reco::Candidate *tmpGenParticle;
    for (unsigned int i = 0 ; i <realW_->numberOfDaughters() ; i++) {
      tmpGenParticle = realW_->daughter(i);
      if (abs(tmpGenParticle->pdgId()) == 11) {
	bool loopend = false;
	unsigned j = 0;
	while ( j < tmpGenParticle->numberOfDaughters() && !loopend){
	  if (tmpGenParticle->daughter(j)->pdgId()==tmpGenParticle->pdgId()) {
	    if (tmpGenParticle->pdgId() ==  11) electronFromRealWAftPh_ = tmpGenParticle->daughter(j);
	    if (tmpGenParticle->pdgId() == -11) nuFromRealWAftPh_ = tmpGenParticle->daughter(j);
	    loopend=true;
	  }
	  ++j;
	}
      }
    } 
    

    for (unsigned int i = 0 ; i <virtualW_->numberOfDaughters() ; i++) { 
      tmpGenParticle = virtualW_->daughter(i);
      if (abs(tmpGenParticle->pdgId()) == 11) {
	bool loopend=false;
	unsigned j=0;
	while ( j < tmpGenParticle->numberOfDaughters() && !loopend){
	  if (tmpGenParticle->daughter(j)->pdgId()==tmpGenParticle->pdgId()) {
	    if (tmpGenParticle->pdgId() ==  11) electronFromVirtualWAftPh_ = tmpGenParticle->daughter(j);
	    if (tmpGenParticle->pdgId() == -11) nuFromVirtualWAftPh_ = tmpGenParticle->daughter(j);
	    loopend=true;
	  }
	  ++j;
	}
      }
    }
  }
}

*/


// -------- FillTree

void McSelectWW::SelectFillTree(TreeContent& myTreeVariables_){

  LorentzVector leptonsWWOrig[4];
  // electronsWWOrigPh[4];

  leptonsWWOrig[0]= leptonFromW_->p4();
  leptonsWWOrig[1]= nuFromW_->p4();
  leptonsWWOrig[2]= leptonFromZ_->p4();
  leptonsWWOrig[3]= cleptonFromZ_->p4() ;

//   electronsWWOrigPh[0]= electronFromRealWAftPh_->p4();
//   electronsWWOrigPh[1]= nuFromRealWAftPh_->p4();
//   electronsWWOrigPh[2]= electronFromVirtualWAftPh_->p4();
//   electronsWWOrigPh[3]= nuFromVirtualWAftPh_->p4() ;

  for(int imc=0; imc<4; ++imc){

  //   std::cout << " ciao 1 " << electronsWWOrig[imc].px() << std::endl;
//     std::cout << " ciao 2 " << electronsWWOrig[imc].py() << std::endl;
//     std::cout << " ciao 2 " << electronsWWOrig[imc].pz() << std::endl;
//     std::cout << " ciao 4 " << electronsWWOrig[imc].e() << std::endl;


    myTreeVariables_.mcEl_px->push_back(leptonsWWOrig[imc].px());
    myTreeVariables_.mcEl_py->push_back(leptonsWWOrig[imc].py());
    myTreeVariables_.mcEl_pz->push_back(leptonsWWOrig[imc].pz());
    myTreeVariables_.mcEl_e->push_back(leptonsWWOrig[imc].e());



//     std::cout << " ciao 1 " << electronsWWOrigPh[imc].px() << std::endl;
//     std::cout << " ciao 2 " << electronsWWOrigPh[imc].py() << std::endl;
//     std::cout << " ciao 2 " << electronsWWOrigPh[imc].pz() << std::endl;
//     std::cout << " ciao 4 " << electronsWWOrigPh[imc].e() << std::endl;



//     myTreeVariables_.mcElph_px->push_back(electronsWWOrigPh[imc].px());
//     myTreeVariables_.mcElph_py->push_back(electronsWWOrigPh[imc].py());
//     myTreeVariables_.mcElph_pz->push_back(electronsWWOrigPh[imc].pz());
//     myTreeVariables_.mcElph_e->push_back(electronsWWOrigPh[imc].e());

  }

  myTreeVariables_.mcEl_pdgID->push_back(da1_rW_pdgID_);
  myTreeVariables_.mcEl_pdgID->push_back(da2_rW_pdgID_);
  myTreeVariables_.mcEl_pdgID->push_back(da1_rZ_pdgID_);
  myTreeVariables_.mcEl_pdgID->push_back(da2_rZ_pdgID_);
}


// --------------------- 

void McSelectWW::printWWEvent(){
  cout << "---------------------------------------------------------------------------------------------------" << endl;
  cout << "----------------: (px,py,pz,E) | eta | phi | pT | m -----------------------------------------------" << endl;
  cout << " Lepton from W: " <<  leptonFromW_->p4()  << " | " 
       << leptonFromW_->eta()  << " | "  
       << leptonFromW_->phi() << " | " 
       << leptonFromW_->pt() << " | " 
       << leptonFromW_->mass() << endl;
    
  cout << " Nu from W: " <<  nuFromW_->p4() << " | " 
       << nuFromW_->eta() << " | "   
       << nuFromW_->phi() << " | " 
       << nuFromW_->pt()  << " | " 
       << nuFromW_->mass() << endl;
       
  cout << "Lepton from Z: " <<  leptonFromZ_->p4() << " | " 
       << leptonFromZ_->eta()  << " | " 
       << leptonFromZ_->phi() << " | "  
       << leptonFromZ_->pt() << " | " 
       << leptonFromZ_->mass() << endl;
 
  cout << "CLepton from Z: " <<  cleptonFromZ_->p4() << " | " 
       << cleptonFromZ_->eta() << " | " 
       << cleptonFromZ_->phi() << " | "  
       << cleptonFromZ_->pt() << " | " 
       << cleptonFromZ_->mass() << endl;
     
    
  cout << "    real W boson: " << realW_->p4() << " | " << realW_->eta() << " | " <<  realW_->phi() 
       << " | " <<  realW_->pt() << " | " <<  realW_->mass() << endl;
  cout << " real Z boson: " << realZ_->p4() << " | " <<  realZ_->eta() << " | " <<  realZ_->phi() 
       << " | " <<  realZ_->pt() << " | " <<  realZ_->mass() << endl;
  cout << "---------------------------------------------------------------------------------------------------" << endl;
}


