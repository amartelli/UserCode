#ifndef TreeContent_h
#define TreeContent_h

#include "TChain.h" 
#include <vector>
#include <string>



struct TreeContent
{
  unsigned int runN;
  unsigned int evtN;

  //MC INFO
  std::vector<float>* mcEl_px;
  std::vector<float>* mcEl_py; 
  std::vector<float>* mcEl_pz; 
  std::vector<float>* mcEl_e;
  std::vector<int>* mcEl_pdgID;

  std::vector<std::string>* mcEventId;
  std::vector<std::string>* mcWWevent;


  //Selected
  bool evtPresel;


  //Electrons ----- start
  //Ele tracks
  std::vector<float> *ele_trackPositionAtVtxX, *ele_trackPositionAtVtxY, *ele_trackPositionAtVtxZ;
  std::vector<float> *ele_trackPositionAtCaloX, *ele_trackPositionAtCaloY, *ele_trackPositionAtCaloZ;

  std::vector<float> *ele_trackMomentumAtVtxX, *ele_trackMomentumAtVtxY, *ele_trackMomentumAtVtxZ;
  std::vector<float> *ele_trackMomentumAtCaloX, *ele_trackMomentumAtCaloY, *ele_trackMomentumAtCaloZ;

  std::vector<float> *ele_VtxX, *ele_VtxY, *ele_VtxZ;

  //SCref
  std::vector<float> *ele_ScEcal, *ele_ScEraw, *ele_ScEta, *ele_ScPhi;
  std::vector<int> *ele_nClu;

  //Ele ID
  std::vector<float> *ele_covEtaEta, *ele_coviEtaiEta;
  std::vector<float> *ele_e1x5, *ele_e2x5Max, *ele_e5x5;
  std::vector<float> *ele_HoE;

  std::vector<float> *ele_ScEoP, *ele_SeedEoP, *ele_SeedEoPout, *ele_eEleClusteroPout;
  std::vector<float> *ele_deltaEtaSuperClusterTrackAtVtx, *ele_deltaEtaSeedClusterTrackAtCalo;
  std::vector<float> *ele_deltaEtaEleClusterTrackAtCalo, *ele_deltaPhiSuperClusterTrackAtVtx;
  std::vector<float> *ele_deltaPhiSeedClusterTrackAtCalo, *ele_deltaPhiEleClusterTrackAtCalo;

  //Ele Iso
  //03
  std::vector<float> *ele_dr03TkSumPt, *ele_dr03EcalRecHitSumEt;
  std::vector<float> *ele_dr03HcalDepth1TowerSumEt, *ele_dr03HcalDepth2TowerSumEt, *ele_dr03HcalTowerSumEt;
  //04
  std::vector<float> *ele_dr04TkSumPt, *ele_dr04EcalRecHitSumEt;
  std::vector<float> *ele_dr04HcalDepth1TowerSumEt, *ele_dr04HcalDepth2TowerSumEt, *ele_dr04HcalTowerSumEt;
  //  std::vector<float> *ele_lat, *ele_phiLat, *ele_etaLat, *ele_a20, *ele_a42;

  // Ele Cand
  std::vector<int> *ele_charge;
  std::vector<float> *ele_energy, *ele_pt, *ele_et, *ele_momentum, *ele_momentumX, *ele_momentumY, *ele_momentumZ;
  std::vector<float> *ele_vertexX, *ele_vertexY, *ele_vertexZ;
  std::vector<float> *ele_theta, *ele_eta, *ele_phi;
  std::vector<float> *ele_x, *ele_y, *ele_z;
  std::vector<float> *ele_mass, *ele_mt;
  std::vector<int> *ele_pdgId;
  std::vector<int> *ele_nDau;
  int ele_ncand;
  //Electrons ----- end

  // SuperCluster ---- start
  int nSC;
  std::vector<int>  *nBC, *nCrystals;
  std::vector<float> *rawEnergy, *energy, *eta, *phi;
  // SuperCluster ---- end


  //Muons  ----- start

  //Muon track
  std::vector<float> *muon_pxAtInner, *muon_pyAtInner, *muon_pzAtInner;
  std::vector<float> *muon_xAtInner, *muon_yAtInner, *muon_zAtInner;
  std::vector<float> *muon_pxAtOuter, *muon_pyAtOuter, *muon_pzAtOuter;
  std::vector<float> *muon_xAtOuter, *muon_yAtOuter, *muon_zAtOuter;
  std::vector<float> *muon_TrackVx, *muon_TrackVy, *muon_TrackVz;

  //Muon Iso & altro
  std::vector<int> *muon_isGlobal, *muon_isTracker, *muon_isStandAlone, *muon_isCalo;

  std::vector<float> *muon_sumPt03, *muon_emEt03, *muon_hadEt03, *muon_hoEt03, *muon_nTrk03, *muon_nJets03;
  std::vector<float> *muon_sumPt05, *muon_emEt05, *muon_hadEt05, *muon_hoEt05, *muon_nTrk05, *muon_nJets05;
  
  std::vector<float> *muon_EcalExpDepo, *muon_HcalExpDepo, *muon_HoExpDepo;
  std::vector<float> *muon_emS9, *muon_hadS9, *muon_hoS9, *muon_CaloComp;

  // Muon Cand
  std::vector<int> *muon_charge;
  std::vector<float> *muon_energy, *muon_pt, *muon_et, *muon_momentum, *muon_momentumX, *muon_momentumY, *muon_momentumZ;
  std::vector<float> *muon_vertexX, *muon_vertexY, *muon_vertexZ;
  std::vector<float> *muon_theta, *muon_eta, *muon_phi;
  std::vector<float> *muon_x, *muon_y, *muon_z;
  std::vector<float> *muon_mass, *muon_mt;
  std::vector<int> *muon_pdgId;
  std::vector<int> *muon_nDau;
  int muon_ncand;
  //Muons  ----- end



  //Tracks ---- start
  std::vector<int> *track_vtxIndex;
  std::vector<float> *track_vtxWeight;

  std::vector<float> *track_pxAtInner, *track_pyAtInner, *track_pzAtInner;
  std::vector<float> *track_xAtInner, *track_yAtInner, *track_zAtInner;
  std::vector<float> *track_pxAtOuter, *track_pyAtOuter, *track_pzAtOuter;
  std::vector<float> *track_xAtOuter, *track_yAtOuter, *track_zAtOuter;


  //trk quality
  std::vector<float> *track_ValidHits, *track_LostHits;
  std::vector<float> *track_NormalizedChi2, *track_recHitsSize;

  std::vector<float> *track_Vx, *track_Vy, *track_Vz;

  std::vector<float> *track_Dxy, *track_D0, *track_Dsz, *track_Dz;
  std::vector<float> *track_DxyError, *track_D0Error, *track_DszError, *track_DzError;

  std::vector<float> *track_DxyPV, /**track_D0PV,*/ *track_DszPV, *track_DzPV;
  //  std::vector<float> *trackDxyErrorPV, *trackD0ErrorPV, *trackDszErrorPV, *trackDzErrorPV;

  // Track Cand
  std::vector<int> *track_charge;
  std::vector<float> *track_energy, *track_pt, *track_et, *track_momentum;
  std::vector<float> *track_momentumX, *track_momentumY, *track_momentumZ;
  std::vector<float> *track_vertexX, *track_vertexY, *track_vertexZ;
  std::vector<float> *track_theta, *track_eta, *track_phi;
  std::vector<float> *track_x, *track_y, *track_z;
  std::vector<float> *track_mass, *track_mt;
  std::vector<int> *track_pdgId;
  std::vector<int> *track_nDau;
  int track_ncand;
  // Track ---- end

  //Vertex
  std::vector<float> *PVx, *PVy, *PVz;
  std::vector<float> *PVErrx, *PVErry, *PVErrz;
  std::vector<float> *SumPt, *ndof, *chi2;  

  //MET
  std::vector<int> *met_charge;
  std::vector<float> *met_energy, *met_et, *met_momentum;
  std::vector<float> *met_vertexX, *met_vertexY, *met_vertexZ;
  std::vector<float> *met_theta, *met_eta, *met_phi;
  std::vector<float> *met_x, *met_y, *met_z;
  std::vector<float> *met_mass, *met_mt;
  std::vector<int> *met_pdgId;
  std::vector<int> *met_nDau;
  int met_ncand;

  //GENMET 
  std::vector<int> *genmet_charge;
  std::vector<float> *genmet_energy, *genmet_et, *genmet_momentum;
  std::vector<float> *genmet_vertexX, *genmet_vertexY, *genmet_vertexZ;
  std::vector<float> *genmet_theta, *genmet_eta, *genmet_phi;
  std::vector<float> *genmet_x, *genmet_y, *genmet_z;
  std::vector<float> *genmet_mass, *genmet_mt;
  std::vector<int> *genmet_pdgId;
  std::vector<int> *genmet_nDau;
  int genmet_ncand;



  //PFMET
  /*  std::vector<int> *pfmet_charge;
  std::vector<float> *pfmet_energy, *pfmet_et, *pfmet_momentum;
  std::vector<float> *pfmet_vertexX, *pfmet_vertexY, *pfmet_vertexZ;
  std::vector<float> *pfmet_theta, *pfmet_eta, *pfmet_phi;
  std::vector<float> *pfmet_x, *pfmet_y, *pfmet_z;
  std::vector<float> *pfmet_mass, *pfmet_mt;
  std::vector<int> *pfmet_pdgId;
  std::vector<int> *pfmet_nDau;
  int pfmet_ncand;
  */
  
  //CORRJET
  /* std::vector<float>* corrjet_alpha;
  std::vector<float>* corrjet_emFrac;
  std::vector<float>* corrjet_hadFrac;
  std::vector<float>* corrjet_flavourId;
  */

  //CORR JETCANDIDATE         
  /*  std::vector<int> *corrjet_charge;
  std::vector<float> *corrjet_energy, *corrjet_et, *corrjet_momentum;
  std::vector<float> *corrjet_vertexX, *corrjet_vertexY, *corrjet_vertexZ;
  std::vector<float> *corrjet_theta, *corrjet_eta, *corrjet_phi;
  std::vector<float> *corrjet_x, *corrjet_y, *corrjet_z;
  std::vector<float> *corrjet_mass, *corrjet_mt;
  std::vector<int> *corrjet_pdgId;
  std::vector<int> *corrjet_nDau;
  int corrjet_ncand;
  */

  //JET
  /*  std::vector<float>* jet_alpha;
  std::vector<float>* jet_emFrac;
  std::vector<float>* jet_hadFrac;
  std::vector<float>* jet_flavourId;
  */

  //JETCANDIDATE         

  /*  std::vector<int> *jet_charge;
  std::vector<float> *jet_energy, *jet_et, *jet_momentum;
  std::vector<float> *jet_vertexX, *jet_vertexY, *jet_vertexZ;
  std::vector<float> *jet_theta, *jet_eta, *jet_phi;
  std::vector<float> *jet_x, *jet_y, *jet_z;
  std::vector<float> *jet_mass, *jet_mt;
  std::vector<int> *jet_pdgId;
  std::vector<int> *jet_nDau;
  int jet_ncand;
  */

  //PFCORRJET

  /*  std::vector<float>* pfcorrjet_alpha;
  std::vector<float>* pfcorrjet_emFrac;
  std::vector<float>* pfcorrjet_hadFrac;
  std::vector<float>* pfcorrjet_flavourId;
  */

  //PFCORRJETCANDIDATE         
  /*  std::vector<int> *pfcorrjet_charge;
  std::vector<float> *pfcorrjet_energy, *pfcorrjet_et, *pfcorrjet_momentum;
  std::vector<float> *pfcorrjet_vertexX, *pfcorrjet_vertexY, *pfcorrjet_vertexZ;
  std::vector<float> *pfcorrjet_theta, *pfcorrjet_eta, *pfcorrjet_phi;
  std::vector<float> *pfcorrjet_x, *pfcorrjet_y, *pfcorrjet_z;
  std::vector<float> *pfcorrjet_mass, *pfcorrjet_mt;
  std::vector<int> *pfcorrjet_pdgId;
  std::vector<int> *pfcorrjet_nDau;
  int pfcorrjet_ncand;
  */


  //GENJET
  /*  std::vector<float>* genjet_alpha;
  std::vector<float>* genjet_emFrac;
  std::vector<float>* genjet_hadFrac;
  std::vector<float>* genjet_flavourId;
  */

  //GENJETCANDIDATE         
  /*  std::vector<int> *genjet_charge;
  std::vector<float> *genjet_energy, *genjet_et, *genjet_momentum;
  std::vector<float> *genjet_vertexX, *genjet_vertexY, *genjet_vertexZ;
  std::vector<float> *genjet_theta, *genjet_eta, *genjet_phi;
  std::vector<float> *genjet_x, *genjet_y, *genjet_z;
  std::vector<float> *genjet_mass, *genjet_mt;
  std::vector<int> *genjet_pdgId;
  std::vector<int> *genjet_nDau;
  int genjet_ncand;
  */

};




// ------------------------------------------------------------------------
//! branch addresses settings

void setBranchAddresses(TTree* chain, TreeContent& treeVars);






// ------------------------------------------------------------------------
//! create branches for a tree

void setBranches(TTree* chain, TreeContent& treeVars);






// ------------------------------------------------------------------------
//! initialize branches

void initializeBranches(TTree* chain, TreeContent& treeVars);



#endif

