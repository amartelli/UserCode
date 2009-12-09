#include "HiggsAnalysis/WWTo2l2nu/interface/TreeContent.h"

#include <iostream>



void setBranchAddresses(TTree* chain, TreeContent& treeVars)
{
  // RunINFO
  chain->SetBranchAddress("runN", &treeVars.runN);
  chain->SetBranchAddress("evtN", &treeVars.evtN);

  //MC INFO
  treeVars.mcEl_px = new std::vector<float>;
  treeVars.mcEl_py = new std::vector<float>;
  treeVars.mcEl_pz = new std::vector<float>;
  treeVars.mcEl_e = new std::vector<float>;
  treeVars.mcEl_pdgID = new std::vector<int>;

  chain->SetBranchAddress("mcEl_px", &treeVars.mcEl_px);
  chain->SetBranchAddress("mcEl_py", &treeVars.mcEl_py);
  chain->SetBranchAddress("mcEl_pz", &treeVars.mcEl_pz);
  chain->SetBranchAddress("mcEl_e", &treeVars.mcEl_e);
  chain->SetBranchAddress("mcEl_pdgID", &treeVars.mcEl_pdgID);

  treeVars.mcEventId = new std::vector<std::string>;
  treeVars.mcWWevent = new std::vector<std::string>;

  chain->SetBranchAddress("mcEventId", &treeVars.mcEventId);
  chain->SetBranchAddress("mcWWevent", &treeVars.mcWWevent);


  //PRESELECTED events
  chain->SetBranchAddress("evtPresel",  &treeVars.evtPresel);

  //Electrons ---- start
  // ---------- Ele track
  treeVars.ele_trackPositionAtVtxX = new std::vector<float>;
  treeVars.ele_trackPositionAtVtxY = new std::vector<float>;
  treeVars.ele_trackPositionAtVtxZ = new std::vector<float>;

  treeVars.ele_trackPositionAtCaloX = new std::vector<float>;
  treeVars.ele_trackPositionAtCaloY = new std::vector<float>;
  treeVars.ele_trackPositionAtCaloZ = new std::vector<float>;

  treeVars.ele_trackMomentumAtVtxX = new std::vector<float>;
  treeVars.ele_trackMomentumAtVtxY = new std::vector<float>;
  treeVars.ele_trackMomentumAtVtxZ = new std::vector<float>;

  treeVars.ele_trackMomentumAtCaloX = new std::vector<float>;
  treeVars.ele_trackMomentumAtCaloY = new std::vector<float>;
  treeVars.ele_trackMomentumAtCaloZ = new std::vector<float>;

  treeVars.ele_VtxX = new std::vector<float>;
  treeVars.ele_VtxY = new std::vector<float>;
  treeVars.ele_VtxZ = new std::vector<float>;

  // --------  Ele Sc
  treeVars.ele_ScEcal = new std::vector<float>;
  treeVars.ele_ScEraw = new std::vector<float>;
  treeVars.ele_ScEta = new std::vector<float>;
  treeVars.ele_ScPhi = new std::vector<float>;
  treeVars.ele_nClu = new std::vector<int>;

  // ---------- Ele ID
                                                                                                                   
  treeVars.ele_covEtaEta = new std::vector<float>;
  treeVars.ele_coviEtaiEta = new std::vector<float>;
  treeVars.ele_e1x5 = new std::vector<float>;
  treeVars.ele_e2x5Max = new std::vector<float>;
  treeVars.ele_e5x5 = new std::vector<float>;
  treeVars.ele_HoE = new std::vector<float>;

  treeVars.ele_ScEoP = new std::vector<float>;
  treeVars.ele_SeedEoP = new std::vector<float>;
  treeVars.ele_SeedEoPout = new std::vector<float>;
  treeVars.ele_eEleClusteroPout = new std::vector<float>;
  treeVars.ele_deltaEtaSuperClusterTrackAtVtx = new std::vector<float>;
  treeVars.ele_deltaEtaSeedClusterTrackAtCalo = new std::vector<float>;
  treeVars.ele_deltaEtaEleClusterTrackAtCalo = new std::vector<float>;
  treeVars.ele_deltaPhiSuperClusterTrackAtVtx = new std::vector<float>;
  treeVars.ele_deltaPhiSeedClusterTrackAtCalo = new std::vector<float>;
  treeVars.ele_deltaPhiEleClusterTrackAtCalo = new std::vector<float>;

  // ------------ Ele Iso
  //--03
  treeVars.ele_dr03TkSumPt = new std::vector<float>;
  treeVars.ele_dr03EcalRecHitSumEt = new std::vector<float>;
  treeVars.ele_dr03HcalDepth1TowerSumEt = new std::vector<float>;
  treeVars.ele_dr03HcalDepth2TowerSumEt = new std::vector<float>;
  treeVars.ele_dr03HcalTowerSumEt = new std::vector<float>;
  //--04  
  treeVars.ele_dr04TkSumPt = new std::vector<float>;
  treeVars.ele_dr04EcalRecHitSumEt = new std::vector<float>;
  treeVars.ele_dr04HcalDepth1TowerSumEt = new std::vector<float>;
  treeVars.ele_dr04HcalDepth2TowerSumEt = new std::vector<float>;
  treeVars.ele_dr04HcalTowerSumEt = new std::vector<float>;

  // ----- Ele Cand
  treeVars.ele_charge = new std::vector<int>;
  treeVars.ele_energy = new std::vector<float>;
  treeVars.ele_pt = new std::vector<float>;
  treeVars.ele_et = new std::vector<float>;
  treeVars.ele_momentum = new std::vector<float>;
  treeVars.ele_momentumX = new std::vector<float>;
  treeVars.ele_momentumY = new std::vector<float>;
  treeVars.ele_momentumZ = new std::vector<float>;
  treeVars.ele_theta = new std::vector<float>;
  treeVars.ele_eta = new std::vector<float>;
  treeVars.ele_phi = new std::vector<float>;
  treeVars.ele_x = new std::vector<float>;
  treeVars.ele_y = new std::vector<float>;
  treeVars.ele_z = new std::vector<float>;
  treeVars.ele_vertexX = new std::vector<float>;
  treeVars.ele_vertexY = new std::vector<float>;
  treeVars.ele_vertexZ = new std::vector<float>;
  treeVars.ele_mass = new std::vector<float>;
  treeVars.ele_mt = new std::vector<float>;
  treeVars.ele_pdgId = new std::vector<int>;
  treeVars.ele_nDau = new std::vector<int>;

    
  chain->SetBranchAddress("ele_trackPositionAtVtxX", &treeVars.ele_trackPositionAtVtxX);
  chain->SetBranchAddress("ele_trackPositionAtVtxY", &treeVars.ele_trackPositionAtVtxY);
  chain->SetBranchAddress("ele_trackPositionAtVtxZ", &treeVars.ele_trackPositionAtVtxZ);
  chain->SetBranchAddress("ele_trackPositionAtCaloX", &treeVars.ele_trackPositionAtCaloX);
  chain->SetBranchAddress("ele_trackPositionAtCaloY", &treeVars.ele_trackPositionAtCaloY);
  chain->SetBranchAddress("ele_trackPositionAtCaloZ", &treeVars.ele_trackPositionAtCaloZ);
  chain->SetBranchAddress("ele_trackMomentumAtVtxX", &treeVars.ele_trackMomentumAtVtxX);
  chain->SetBranchAddress("ele_trackMomentumAtVtxY", &treeVars.ele_trackMomentumAtVtxY);
  chain->SetBranchAddress("ele_trackMomentumAtVtxZ", &treeVars.ele_trackMomentumAtVtxZ);
  chain->SetBranchAddress("ele_trackMomentumAtCaloX", &treeVars.ele_trackMomentumAtCaloX);
  chain->SetBranchAddress("ele_trackMomentumAtCaloY", &treeVars.ele_trackMomentumAtCaloY);
  chain->SetBranchAddress("ele_trackMomentumAtCaloZ", &treeVars.ele_trackMomentumAtCaloZ);
  chain->SetBranchAddress("ele_VtxX", &treeVars.ele_VtxX);
  chain->SetBranchAddress("ele_VtxY", &treeVars.ele_VtxY);
  chain->SetBranchAddress("ele_VtxZ", &treeVars.ele_VtxZ);

  chain->SetBranchAddress("ele_ScEcal", &treeVars.ele_ScEcal);
  chain->SetBranchAddress("ele_ScEraw", &treeVars.ele_ScEraw);
  chain->SetBranchAddress("ele_ScEta", &treeVars.ele_ScEta);
  chain->SetBranchAddress("ele_ScPhi", &treeVars.ele_ScPhi);
  chain->SetBranchAddress("ele_nClu", &treeVars.ele_nClu);

  chain->SetBranchAddress("ele_covEtaEta", &treeVars.ele_covEtaEta);
  chain->SetBranchAddress("ele_coviEtaiEta", &treeVars.ele_coviEtaiEta);
  chain->SetBranchAddress("ele_e1x5", &treeVars.ele_e1x5);
  chain->SetBranchAddress("ele_e2x5Max", &treeVars.ele_e2x5Max);
  chain->SetBranchAddress("ele_e5x5", &treeVars.ele_e5x5);
  chain->SetBranchAddress("ele_HoE", &treeVars.ele_HoE);

  chain->SetBranchAddress("ele_ScEoP", &treeVars.ele_ScEoP);
  chain->SetBranchAddress("ele_SeedEoP", &treeVars.ele_SeedEoP);
  chain->SetBranchAddress("ele_eEleClusteroPout", &treeVars.ele_eEleClusteroPout);
  chain->SetBranchAddress("ele_deltaEtaSuperClusterTrackAtVtx", &treeVars.ele_deltaEtaSuperClusterTrackAtVtx);
  chain->SetBranchAddress("ele_deltaEtaSeedClusterTrackAtCalo", &treeVars.ele_deltaEtaSeedClusterTrackAtCalo);
  chain->SetBranchAddress("ele_deltaEtaEleClusterTrackAtCalo", &treeVars.ele_deltaEtaEleClusterTrackAtCalo);
  chain->SetBranchAddress("ele_deltaPhiSuperClusterTrackAtVtx", &treeVars.ele_deltaPhiSuperClusterTrackAtVtx);
  chain->SetBranchAddress("ele_deltaPhiSeedClusterTrackAtCalo", &treeVars.ele_deltaPhiSeedClusterTrackAtCalo);
  chain->SetBranchAddress("ele_deltaPhiEleClusterTrackAtCalo", &treeVars.ele_deltaPhiEleClusterTrackAtCalo);


  chain->SetBranchAddress("ele_dr03TkSumPt", &treeVars.ele_dr03TkSumPt);
  chain->SetBranchAddress("ele_dr03EcalRecHitSumEt", &treeVars.ele_dr03EcalRecHitSumEt);
  chain->SetBranchAddress("ele_dr03HcalDepth1TowerSumEt", &treeVars.ele_dr03HcalDepth1TowerSumEt);
  chain->SetBranchAddress("ele_dr03HcalDepth2TowerSumEt", &treeVars.ele_dr03HcalDepth2TowerSumEt);
  chain->SetBranchAddress("ele_dr03HcalTowerSumEt", &treeVars.ele_dr03HcalTowerSumEt);

  chain->SetBranchAddress("ele_dr04TkSumPt", &treeVars.ele_dr04TkSumPt);
  chain->SetBranchAddress("ele_dr04EcalRecHitSumEt", &treeVars.ele_dr04EcalRecHitSumEt);
  chain->SetBranchAddress("ele_dr04HcalDepth1TowerSumEt", &treeVars.ele_dr04HcalDepth1TowerSumEt);
  chain->SetBranchAddress("ele_dr04HcalDepth2TowerSumEt", &treeVars.ele_dr04HcalDepth2TowerSumEt);
  chain->SetBranchAddress("ele_dr04HcalTowerSumEt", &treeVars.ele_dr04HcalTowerSumEt);

  chain->SetBranchAddress("ele_charge", &treeVars.ele_charge);
  chain->SetBranchAddress("ele_energy", &treeVars.ele_energy);
  chain->SetBranchAddress("ele_pt", &treeVars.ele_pt);
  chain->SetBranchAddress("ele_et", &treeVars.ele_et);
  chain->SetBranchAddress("ele_momentum", &treeVars.ele_momentum);
  chain->SetBranchAddress("ele_momentumX", &treeVars.ele_momentumX);
  chain->SetBranchAddress("ele_momentumY", &treeVars.ele_momentumY);
  chain->SetBranchAddress("ele_momentumZ", &treeVars.ele_momentumZ);
  chain->SetBranchAddress("ele_theta", &treeVars.ele_theta);
  chain->SetBranchAddress("ele_eta", &treeVars.ele_eta);
  chain->SetBranchAddress("ele_phi", &treeVars.ele_phi);  
  chain->SetBranchAddress("ele_x", &treeVars.ele_x);
  chain->SetBranchAddress("ele_y", &treeVars.ele_y);
  chain->SetBranchAddress("ele_z", &treeVars.ele_z);
  chain->SetBranchAddress("ele_vertexX", &treeVars.ele_vertexX);
  chain->SetBranchAddress("ele_vertexY", &treeVars.ele_vertexY);
  chain->SetBranchAddress("ele_vertexZ", &treeVars.ele_vertexZ);
  chain->SetBranchAddress("ele_mass", &treeVars.ele_mass);
  chain->SetBranchAddress("ele_mt", &treeVars.ele_mt);
  chain->SetBranchAddress("ele_pdgId", &treeVars.ele_pdgId);
  chain->SetBranchAddress("ele_ncand", &treeVars.ele_ncand);
  chain->SetBranchAddress("ele_nDau", &treeVars.ele_nDau);
  // Electron ---------- end

  //SClusters ---- start
  treeVars.nBC = new std::vector<int>;
  treeVars.nCrystals = new std::vector<int>;
  treeVars.rawEnergy = new std::vector<float>;
  treeVars.energy = new std::vector<float>;
  treeVars.eta = new std::vector<float>;
  treeVars.phi = new std::vector<float>;

  chain->SetBranchAddress("nBC", &treeVars.nBC);
  chain->SetBranchAddress("nCrystals", &treeVars.nCrystals);
  chain->SetBranchAddress("rawEnergy", &treeVars.rawEnergy);
  chain->SetBranchAddress("energy", &treeVars.energy);
  chain->SetBranchAddress("eta", &treeVars.eta);
  chain->SetBranchAddress("phi", &treeVars.phi);
  //SClusters ---- end

  //Muons --------- start

  //Muon track
  treeVars.muon_pxAtInner = new std::vector<float>;
  treeVars.muon_pyAtInner = new std::vector<float>;
  treeVars.muon_pzAtInner = new std::vector<float>;
  treeVars.muon_xAtInner = new std::vector<float>;
  treeVars.muon_yAtInner = new std::vector<float>;
  treeVars.muon_zAtInner = new std::vector<float>;
  treeVars.muon_pxAtOuter = new std::vector<float>;
  treeVars.muon_pyAtOuter = new std::vector<float>;
  treeVars.muon_pzAtOuter = new std::vector<float>;
  treeVars.muon_xAtOuter = new std::vector<float>;
  treeVars.muon_yAtOuter = new std::vector<float>;
  treeVars.muon_zAtOuter = new std::vector<float>;
  treeVars.muon_TrackVx = new std::vector<float>;
  treeVars.muon_TrackVy = new std::vector<float>;
  treeVars.muon_TrackVz = new std::vector<float>;

  treeVars.muon_isGlobal = new std::vector<int>;
  treeVars.muon_isTracker = new std::vector<int>;
  treeVars.muon_isStandAlone = new std::vector<int>;
  treeVars.muon_isCalo = new std::vector<int>;
  
  treeVars.muon_sumPt03 = new std::vector<float>;
  treeVars.muon_emEt03 = new std::vector<float>;
  treeVars.muon_hadEt03 = new std::vector<float>;
  treeVars.muon_hoEt03 = new std::vector<float>;
  treeVars.muon_nTrk03 = new std::vector<float>;
  treeVars.muon_nJets03 = new std::vector<float>;
  treeVars.muon_sumPt05 = new std::vector<float>;
  treeVars.muon_emEt05 = new std::vector<float>;
  treeVars.muon_hadEt05 = new std::vector<float>;
  treeVars.muon_hoEt05 = new std::vector<float>;
  treeVars.muon_nTrk05 = new std::vector<float>;
  treeVars.muon_nJets05 = new std::vector<float>;

  treeVars.muon_EcalExpDepo = new std::vector<float>;
  treeVars.muon_HcalExpDepo = new std::vector<float>;
  treeVars.muon_HoExpDepo = new std::vector<float>;
  treeVars.muon_emS9 = new std::vector<float>;
  treeVars.muon_hadS9 = new std::vector<float>;
  treeVars.muon_hoS9 = new std::vector<float>;
  treeVars.muon_CaloComp = new std::vector<float>;

  //Muon cand
  treeVars.muon_charge = new std::vector<int>;
  treeVars.muon_energy = new std::vector<float>;
  treeVars.muon_pt = new std::vector<float>;
  treeVars.muon_et = new std::vector<float>;
  treeVars.muon_momentum = new std::vector<float>;
  treeVars.muon_momentumX = new std::vector<float>;
  treeVars.muon_momentumY = new std::vector<float>;
  treeVars.muon_momentumZ = new std::vector<float>;
  treeVars.muon_theta = new std::vector<float>;
  treeVars.muon_eta = new std::vector<float>;
  treeVars.muon_phi = new std::vector<float>;
  treeVars.muon_x = new std::vector<float>;
  treeVars.muon_y = new std::vector<float>;
  treeVars.muon_z = new std::vector<float>;
  treeVars.muon_vertexX = new std::vector<float>;
  treeVars.muon_vertexY = new std::vector<float>;
  treeVars.muon_vertexZ = new std::vector<float>;
  treeVars.muon_mass = new std::vector<float>;
  treeVars.muon_mt = new std::vector<float>;
  treeVars.muon_pdgId = new std::vector<int>;
  treeVars.muon_nDau = new std::vector<int>;

  chain->SetBranchAddress("muon_pxAtInner", &treeVars.muon_pxAtInner);
  chain->SetBranchAddress("muon_pyAtInner", &treeVars.muon_pyAtInner);
  chain->SetBranchAddress("muon_pzAtInner", &treeVars.muon_pzAtInner);
  chain->SetBranchAddress("muon_xAtInner", &treeVars.muon_xAtInner);
  chain->SetBranchAddress("muon_yAtInner", &treeVars.muon_yAtInner);
  chain->SetBranchAddress("muon_zAtInner", &treeVars.muon_zAtInner);
  chain->SetBranchAddress("muon_pxAtOuter", &treeVars.muon_pxAtOuter);
  chain->SetBranchAddress("muon_pyAtOuter", &treeVars.muon_pyAtOuter);
  chain->SetBranchAddress("muon_pzAtOuter", &treeVars.muon_pzAtOuter);
  chain->SetBranchAddress("muon_xAtOuter", &treeVars.muon_xAtOuter);
  chain->SetBranchAddress("muon_yAtOuter", &treeVars.muon_yAtOuter);
  chain->SetBranchAddress("muon_zAtOuter", &treeVars.muon_zAtOuter);
  chain->SetBranchAddress("muon_TrackVx", &treeVars.muon_TrackVx);
  chain->SetBranchAddress("muon_TrackVy", &treeVars.muon_TrackVy);
  chain->SetBranchAddress("muon_TrackVz", &treeVars.muon_TrackVz);


  chain->SetBranchAddress("muon_isGlobal", &treeVars.muon_isGlobal);
  chain->SetBranchAddress("muon_isTracker", &treeVars.muon_isTracker);
  chain->SetBranchAddress("muon_isStandAlone", &treeVars.muon_isStandAlone);
  chain->SetBranchAddress("muon_isCalo", &treeVars.muon_isCalo);

  chain->SetBranchAddress("muon_sumPt03", &treeVars.muon_sumPt03);
  chain->SetBranchAddress("muon_emEt03", &treeVars.muon_sumPt03);
  chain->SetBranchAddress("muon_hadEt03", &treeVars.muon_hadEt03);
  chain->SetBranchAddress("muon_hoEt03", &treeVars.muon_hoEt03);
  chain->SetBranchAddress("muon_nTrk03", &treeVars.muon_nTrk03);
  chain->SetBranchAddress("muon_nJets03", &treeVars.muon_nJets03);

  chain->SetBranchAddress("muon_sumPt05", &treeVars.muon_sumPt05);
  chain->SetBranchAddress("muon_emEt05", &treeVars.muon_emEt05);
  chain->SetBranchAddress("muon_hadEt05", &treeVars.muon_hadEt05);
  chain->SetBranchAddress("muon_hoEt05",  &treeVars.muon_hoEt05);
  chain->SetBranchAddress("muon_nTrk05", &treeVars.muon_nTrk05);
  chain->SetBranchAddress("muon_nJets05", &treeVars.muon_nJets05);

  chain->SetBranchAddress("muon_EcalExpDepo", &treeVars.muon_EcalExpDepo);
  chain->SetBranchAddress("muon_HcalExpDepo", &treeVars.muon_HcalExpDepo);
  chain->SetBranchAddress("muon_HoExpDepo", &treeVars.muon_HoExpDepo);
  chain->SetBranchAddress("muon_emS9", &treeVars.muon_emS9);
  chain->SetBranchAddress("muon_hadS9", &treeVars.muon_hadS9);
  chain->SetBranchAddress("muon_hoS9", &treeVars.muon_hoS9);
  chain->SetBranchAddress("muon_CaloComp", &treeVars.muon_CaloComp);


  chain->SetBranchAddress("muon_charge", &treeVars.muon_charge);
  chain->SetBranchAddress("muon_energy", &treeVars.muon_energy);
  chain->SetBranchAddress("muon_pt", &treeVars.muon_pt);
  chain->SetBranchAddress("muon_et", &treeVars.muon_et);
  chain->SetBranchAddress("muon_momentum", &treeVars.muon_momentum);
  chain->SetBranchAddress("muon_momentumX", &treeVars.muon_momentumX);
  chain->SetBranchAddress("muon_momentumY", &treeVars.muon_momentumY);
  chain->SetBranchAddress("muon_momentumZ", &treeVars.muon_momentumZ);
  chain->SetBranchAddress("muon_theta", &treeVars.muon_theta);
  chain->SetBranchAddress("muon_eta", &treeVars.muon_eta);
  chain->SetBranchAddress("muon_phi", &treeVars.muon_phi);  
  chain->SetBranchAddress("muon_x", &treeVars.muon_x);
  chain->SetBranchAddress("muon_y", &treeVars.muon_y);
  chain->SetBranchAddress("muon_z", &treeVars.muon_z);
  chain->SetBranchAddress("muon_vertexX", &treeVars.muon_vertexX);
  chain->SetBranchAddress("muon_vertexY", &treeVars.muon_vertexY);
  chain->SetBranchAddress("muon_vertexZ", &treeVars.muon_vertexZ);
  chain->SetBranchAddress("muon_mass", &treeVars.muon_mass);
  chain->SetBranchAddress("muon_mt", &treeVars.muon_mt);
  chain->SetBranchAddress("muon_pdgId", &treeVars.muon_pdgId);
  chain->SetBranchAddress("muon_ncand", &treeVars.muon_ncand);
  chain->SetBranchAddress("muon_nDau", &treeVars.muon_nDau);
  //Muons ------ end

  //Tracks --- start
  treeVars.track_vtxIndex  = new std::vector<int>;
  treeVars.track_vtxWeight = new std::vector<float>;

  treeVars.track_pxAtInner = new std::vector<float>;
  treeVars.track_pyAtInner = new std::vector<float>;
  treeVars.track_pzAtInner = new std::vector<float>;
  treeVars.track_xAtInner = new std::vector<float>;
  treeVars.track_yAtInner = new std::vector<float>;
  treeVars.track_zAtInner = new std::vector<float>;
  treeVars.track_pxAtOuter = new std::vector<float>;
  treeVars.track_pyAtOuter = new std::vector<float>;
  treeVars.track_pzAtOuter = new std::vector<float>;
  treeVars.track_xAtOuter = new std::vector<float>;
  treeVars.track_yAtOuter = new std::vector<float>;
  treeVars.track_zAtOuter = new std::vector<float>;

  treeVars.track_ValidHits  = new std::vector<float>;
  treeVars.track_LostHits  = new std::vector<float>;
  treeVars.track_NormalizedChi2 = new std::vector<float>;
  treeVars.track_recHitsSize  = new std::vector<float>;

  treeVars.track_Vx = new std::vector<float>;
  treeVars.track_Vy = new std::vector<float>;
  treeVars.track_Vz = new std::vector<float>;

  treeVars.track_Dxy = new std::vector<float>;
  treeVars.track_D0 = new std::vector<float>;
  treeVars.track_Dsz = new std::vector<float>;
  treeVars.track_Dz = new std::vector<float>;

  treeVars.track_DxyError = new std::vector<float>;
  treeVars.track_D0Error = new std::vector<float>;
  treeVars.track_DszError = new std::vector<float>;
  treeVars.track_DzError = new std::vector<float>;

  treeVars.track_DxyPV = new std::vector<float>;
  treeVars.track_DszPV = new std::vector<float>;
  treeVars.track_DzPV = new std::vector<float>;


  //Track Cand
  treeVars.track_charge = new std::vector<int>;
  treeVars.track_energy = new std::vector<float>;
  treeVars.track_pt = new std::vector<float>;
  treeVars.track_et = new std::vector<float>;
  treeVars.track_momentum = new std::vector<float>;
  treeVars.track_momentumX = new std::vector<float>;
  treeVars.track_momentumY = new std::vector<float>;
  treeVars.track_momentumZ = new std::vector<float>;
  treeVars.track_theta = new std::vector<float>;
  treeVars.track_eta = new std::vector<float>;
  treeVars.track_phi = new std::vector<float>;
  treeVars.track_x = new std::vector<float>;
  treeVars.track_y = new std::vector<float>;
  treeVars.track_z = new std::vector<float>;
  treeVars.track_vertexX = new std::vector<float>;
  treeVars.track_vertexY = new std::vector<float>;
  treeVars.track_vertexZ = new std::vector<float>;
  treeVars.track_mass = new std::vector<float>;
  treeVars.track_mt = new std::vector<float>;
  treeVars.track_pdgId = new std::vector<int>;
  treeVars.track_nDau = new std::vector<int>;

  chain->SetBranchAddress("track_vtxIndex", &treeVars.track_vtxIndex);
  chain->SetBranchAddress("track_vtxWeight", &treeVars.track_vtxWeight);

  chain->SetBranchAddress("track_pxAtInner", &treeVars.track_pxAtInner);
  chain->SetBranchAddress("track_pyAtInner", &treeVars.track_pyAtInner);
  chain->SetBranchAddress("track_pzAtInner", &treeVars.track_pzAtInner);
  chain->SetBranchAddress("track_xAtInner", &treeVars.track_xAtInner);
  chain->SetBranchAddress("track_yAtInner", &treeVars.track_yAtInner);
  chain->SetBranchAddress("track_zAtInner", &treeVars.track_zAtInner);
  chain->SetBranchAddress("track_pxAtOuter", &treeVars.track_pxAtOuter);
  chain->SetBranchAddress("track_pyAtOuter", &treeVars.track_pyAtOuter);
  chain->SetBranchAddress("track_pzAtOuter", &treeVars.track_pzAtOuter);
  chain->SetBranchAddress("track_xAtOuter", &treeVars.track_xAtOuter);
  chain->SetBranchAddress("track_yAtOuter", &treeVars.track_yAtOuter);
  chain->SetBranchAddress("track_zAtOuter", &treeVars.track_zAtOuter);


  chain->SetBranchAddress("track_ValidHits", &treeVars.track_ValidHits);
  chain->SetBranchAddress("track_LostHits", &treeVars.track_LostHits);
  chain->SetBranchAddress("track_NormalizedChi2", &treeVars.track_NormalizedChi2);
  chain->SetBranchAddress("track_recHitsSize", &treeVars.track_recHitsSize);

  chain->SetBranchAddress("track_Vx", &treeVars.track_Vx);
  chain->SetBranchAddress("track_Vy", &treeVars.track_Vy);
  chain->SetBranchAddress("track_Vz", &treeVars.track_Vz);


  chain->SetBranchAddress("track_Dxy", &treeVars.track_Dxy);
  chain->SetBranchAddress("track_D0", &treeVars.track_D0);
  chain->SetBranchAddress("track_Dsz", &treeVars.track_Dsz);
  chain->SetBranchAddress("track_Dz", &treeVars.track_Dz);
  chain->SetBranchAddress("track_DxyError", &treeVars.track_DxyError);
  chain->SetBranchAddress("track_D0Error", &treeVars.track_D0Error);
  chain->SetBranchAddress("track_DszError", &treeVars.track_DszError);
  chain->SetBranchAddress("track_DzError", &treeVars.track_DzError);
  chain->SetBranchAddress("track_DxyPV", &treeVars.track_DxyPV);
  chain->SetBranchAddress("track_DszPV", &treeVars.track_DszPV);
  chain->SetBranchAddress("track_DzPV", &treeVars.track_DzPV);


  chain->SetBranchAddress("track_charge", &treeVars.track_charge);
  chain->SetBranchAddress("track_energy", &treeVars.track_energy);
  chain->SetBranchAddress("track_pt", &treeVars.track_pt);
  chain->SetBranchAddress("track_et", &treeVars.track_et);
  chain->SetBranchAddress("track_momentum", &treeVars.track_momentum);
  chain->SetBranchAddress("track_momentumX", &treeVars.track_momentumX);
  chain->SetBranchAddress("track_momentumY", &treeVars.track_momentumY);
  chain->SetBranchAddress("track_momentumZ", &treeVars.track_momentumZ);
  chain->SetBranchAddress("track_theta", &treeVars.track_theta);
  chain->SetBranchAddress("track_eta", &treeVars.track_eta);
  chain->SetBranchAddress("track_phi", &treeVars.track_phi);
  chain->SetBranchAddress("track_x", &treeVars.track_x);
  chain->SetBranchAddress("track_y", &treeVars.track_y);
  chain->SetBranchAddress("track_z", &treeVars.track_z);
  chain->SetBranchAddress("track_vertexX", &treeVars.track_vertexX);
  chain->SetBranchAddress("track_vertexY", &treeVars.track_vertexY);
  chain->SetBranchAddress("track_vertexZ", &treeVars.track_vertexZ);
  chain->SetBranchAddress("track_mass", &treeVars.track_mass);
  chain->SetBranchAddress("track_mt", &treeVars.track_mt);
  chain->SetBranchAddress("track_pdgId", &treeVars.track_pdgId);
  chain->SetBranchAddress("track_ncand", &treeVars.track_ncand);
  chain->SetBranchAddress("track_nDau", &treeVars.track_nDau);
  //Tracks ---- end

  //Vertex --- start
  treeVars.PVx = new std::vector<float>;
  treeVars.PVy = new std::vector<float>;
  treeVars.PVz = new std::vector<float>;

  treeVars.PVErrx = new std::vector<float>;
  treeVars.PVErry = new std::vector<float>;
  treeVars.PVErrz = new std::vector<float>;

  treeVars.SumPt = new std::vector<float>;
  treeVars.ndof = new std::vector<float>;
  treeVars.chi2 = new std::vector<float>;

  chain->SetBranchAddress("PVx", &treeVars.PVx);
  chain->SetBranchAddress("PVy", &treeVars.PVy);
  chain->SetBranchAddress("PVz", &treeVars.PVz);
  chain->SetBranchAddress("PVErrx", &treeVars.PVErrx);
  chain->SetBranchAddress("PVErry", &treeVars.PVErry);
  chain->SetBranchAddress("PVErrz", &treeVars.PVErrz);

  chain->SetBranchAddress("SumPt", &treeVars.SumPt);
  chain->SetBranchAddress("ndof", &treeVars.ndof);
  chain->SetBranchAddress("chi2", &treeVars.chi2);
  //Vertex --- end

  // missing transvers energy ---- start

  //MET
  treeVars.met_charge = new std::vector<int>;
  treeVars.met_energy = new std::vector<float>;
  treeVars.met_et = new std::vector<float>;
  treeVars.met_momentum = new std::vector<float>;
  treeVars.met_theta = new std::vector<float>;
  treeVars.met_eta = new std::vector<float>;
  treeVars.met_phi = new std::vector<float>;
  treeVars.met_x = new std::vector<float>;
  treeVars.met_y = new std::vector<float>;
  treeVars.met_z = new std::vector<float>;
  treeVars.met_vertexX = new std::vector<float>;
  treeVars.met_vertexY = new std::vector<float>;
  treeVars.met_vertexZ = new std::vector<float>;
  treeVars.met_mass = new std::vector<float>;
  treeVars.met_mt = new std::vector<float>;
  treeVars.met_pdgId = new std::vector<int>;
  treeVars.met_nDau = new std::vector<int>;

  chain->SetBranchAddress("met_charge", &treeVars.met_charge);
  chain->SetBranchAddress("met_energy", &treeVars.met_energy);
  chain->SetBranchAddress("met_et", &treeVars.met_et);
  chain->SetBranchAddress("met_momentum", &treeVars.met_momentum);
  chain->SetBranchAddress("met_theta", &treeVars.met_theta);
  chain->SetBranchAddress("met_eta", &treeVars.met_eta);
  chain->SetBranchAddress("met_phi", &treeVars.met_phi);
  chain->SetBranchAddress("met_x", &treeVars.met_x);
  chain->SetBranchAddress("met_y", &treeVars.met_y);
  chain->SetBranchAddress("met_z", &treeVars.met_z);
  chain->SetBranchAddress("met_vertexX", &treeVars.met_vertexX);
  chain->SetBranchAddress("met_vertexY", &treeVars.met_vertexY);
  chain->SetBranchAddress("met_vertexZ", &treeVars.met_vertexZ);
  chain->SetBranchAddress("met_mass", &treeVars.met_mass);
  chain->SetBranchAddress("met_mt", &treeVars.met_mt);
  chain->SetBranchAddress("met_pdgId", &treeVars.met_pdgId);
  chain->SetBranchAddress("met_ncand", &treeVars.met_ncand);
  chain->SetBranchAddress("met_nDau", &treeVars.met_nDau);

  //GENMET
  treeVars.genmet_charge = new std::vector<int>;
  treeVars.genmet_energy = new std::vector<float>;
  treeVars.genmet_et = new std::vector<float>;
  treeVars.genmet_momentum = new std::vector<float>;
  treeVars.genmet_theta = new std::vector<float>;
  treeVars.genmet_eta = new std::vector<float>;
  treeVars.genmet_phi = new std::vector<float>;
  treeVars.genmet_x = new std::vector<float>;
  treeVars.genmet_y = new std::vector<float>;
  treeVars.genmet_z = new std::vector<float>;
  treeVars.genmet_vertexX = new std::vector<float>;
  treeVars.genmet_vertexY = new std::vector<float>;
  treeVars.genmet_vertexZ = new std::vector<float>;
  treeVars.genmet_mass = new std::vector<float>;
  treeVars.genmet_mt = new std::vector<float>;
  treeVars.genmet_pdgId = new std::vector<int>;
  treeVars.genmet_nDau = new std::vector<int>;


  chain->SetBranchAddress("genmet_charge", &treeVars.genmet_charge);
  chain->SetBranchAddress("genmet_energy", &treeVars.genmet_energy);
  chain->SetBranchAddress("genmet_et", &treeVars.genmet_et);
  chain->SetBranchAddress("genmet_momentum", &treeVars.genmet_momentum);
  chain->SetBranchAddress("genmet_theta", &treeVars.genmet_theta);
  chain->SetBranchAddress("genmet_eta", &treeVars.genmet_eta);
  chain->SetBranchAddress("genmet_phi", &treeVars.genmet_phi);
  chain->SetBranchAddress("genmet_x", &treeVars.genmet_x);
  chain->SetBranchAddress("genmet_y", &treeVars.genmet_y);
  chain->SetBranchAddress("genmet_z", &treeVars.genmet_z);
  chain->SetBranchAddress("genmet_vertexX", &treeVars.genmet_vertexX);
  chain->SetBranchAddress("genmet_vertexY", &treeVars.genmet_vertexY);
  chain->SetBranchAddress("genmet_vertexZ", &treeVars.genmet_vertexZ);
  chain->SetBranchAddress("genmet_mass", &treeVars.genmet_mass);
  chain->SetBranchAddress("genmet_mt", &treeVars.genmet_mt);
  chain->SetBranchAddress("genmet_pdgId", &treeVars.genmet_pdgId);
  chain->SetBranchAddress("genmet_ncand", &treeVars.genmet_ncand);
  chain->SetBranchAddress("genmet_nDau", &treeVars.genmet_nDau);

  // missing transvers energy ---- end

}



void setBranches(TTree* chain, TreeContent& treeVars)
{
  // RunINFO 
  chain->Branch("runN", &treeVars.runN, "runN/i");
  chain->Branch("evtN", &treeVars.evtN, "evtN/i");

  //MC INFO
  treeVars.mcEl_px = new std::vector<float> ;
  treeVars.mcEl_py = new std::vector<float> ;
  treeVars.mcEl_pz = new std::vector<float> ;
  treeVars.mcEl_e = new std::vector<float> ;
  treeVars.mcEl_pdgID = new std::vector<int> ;

  chain->Branch("mcEl_px", "std::vector<float>", &treeVars.mcEl_px);
  chain->Branch("mcEl_py", "std::vector<float>", &treeVars.mcEl_py);
  chain->Branch("mcEl_pz", "std::vector<float>", &treeVars.mcEl_pz);
  chain->Branch("mcEl_e", "std::vector<float>", &treeVars.mcEl_e);
  chain->Branch("mcEl_pdgID", "std::vector<int>", &treeVars.mcEl_pdgID);

  treeVars.mcEventId = new std::vector<std::string>;
  treeVars.mcWWevent = new std::vector<std::string>;

  chain->Branch("mcEventId", "std::vector<std::string>", &treeVars.mcEventId);
  chain->Branch("mcWWevent", "std::vector<std::string>", &treeVars.mcWWevent);


  //PRESELECTED events  
  chain->Branch("evtPresel",  &treeVars.evtPresel, "evtPresel/b");


  //Electrons ---- start
  // ---------- Ele track
  treeVars.ele_trackPositionAtVtxX = new std::vector<float>;
  treeVars.ele_trackPositionAtVtxY = new std::vector<float>;
  treeVars.ele_trackPositionAtVtxZ = new std::vector<float>;

  treeVars.ele_trackPositionAtCaloX = new std::vector<float>;
  treeVars.ele_trackPositionAtCaloY = new std::vector<float>;
  treeVars.ele_trackPositionAtCaloZ = new std::vector<float>;

  treeVars.ele_trackMomentumAtVtxX = new std::vector<float>;
  treeVars.ele_trackMomentumAtVtxY = new std::vector<float>;
  treeVars.ele_trackMomentumAtVtxZ = new std::vector<float>;

  treeVars.ele_trackMomentumAtCaloX = new std::vector<float>;
  treeVars.ele_trackMomentumAtCaloY = new std::vector<float>;
  treeVars.ele_trackMomentumAtCaloZ = new std::vector<float>;

  treeVars.ele_VtxX = new std::vector<float>;
  treeVars.ele_VtxY = new std::vector<float>;
  treeVars.ele_VtxZ = new std::vector<float>;

  // --------  Ele Sc
  treeVars.ele_ScEcal = new std::vector<float>;
  treeVars.ele_ScEraw = new std::vector<float>;
  treeVars.ele_ScEta = new std::vector<float>;
  treeVars.ele_ScPhi = new std::vector<float>;
  treeVars.ele_nClu = new std::vector<int>;

  // ---------- Ele ID
      
  treeVars.ele_covEtaEta = new std::vector<float>;
  treeVars.ele_coviEtaiEta = new std::vector<float>;
  treeVars.ele_e1x5 = new std::vector<float>;
  treeVars.ele_e2x5Max = new std::vector<float>;
  treeVars.ele_e5x5 = new std::vector<float>;
  treeVars.ele_HoE = new std::vector<float>;

  treeVars.ele_ScEoP = new std::vector<float>;
  treeVars.ele_SeedEoP = new std::vector<float>;
  treeVars.ele_SeedEoPout = new std::vector<float>;
  treeVars.ele_eEleClusteroPout = new std::vector<float>;
  treeVars.ele_deltaEtaSuperClusterTrackAtVtx = new std::vector<float>;
  treeVars.ele_deltaEtaSeedClusterTrackAtCalo = new std::vector<float>;
  treeVars.ele_deltaEtaEleClusterTrackAtCalo = new std::vector<float>;
  treeVars.ele_deltaPhiSuperClusterTrackAtVtx = new std::vector<float>;
  treeVars.ele_deltaPhiSeedClusterTrackAtCalo = new std::vector<float>;
  treeVars.ele_deltaPhiEleClusterTrackAtCalo = new std::vector<float>;

  // ------------ Ele Iso
  //--03
  treeVars.ele_dr03TkSumPt = new std::vector<float>;
  treeVars.ele_dr03EcalRecHitSumEt = new std::vector<float>;
  treeVars.ele_dr03HcalDepth1TowerSumEt = new std::vector<float>;
  treeVars.ele_dr03HcalDepth2TowerSumEt = new std::vector<float>;
  treeVars.ele_dr03HcalTowerSumEt = new std::vector<float>;
  //--04  
  treeVars.ele_dr04TkSumPt = new std::vector<float>;
  treeVars.ele_dr04EcalRecHitSumEt = new std::vector<float>;
  treeVars.ele_dr04HcalDepth1TowerSumEt = new std::vector<float>;
  treeVars.ele_dr04HcalDepth2TowerSumEt = new std::vector<float>;
  treeVars.ele_dr04HcalTowerSumEt = new std::vector<float>;

  // ----- Ele Cand
  treeVars.ele_charge = new std::vector<int>;
  treeVars.ele_energy = new std::vector<float>;
  treeVars.ele_pt = new std::vector<float>;
  treeVars.ele_et = new std::vector<float>;
  treeVars.ele_momentum = new std::vector<float>;
  treeVars.ele_momentumX = new std::vector<float>;
  treeVars.ele_momentumY = new std::vector<float>;
  treeVars.ele_momentumZ = new std::vector<float>;
  treeVars.ele_theta = new std::vector<float>;
  treeVars.ele_eta = new std::vector<float>;
  treeVars.ele_phi = new std::vector<float>;
  treeVars.ele_x = new std::vector<float>;
  treeVars.ele_y = new std::vector<float>;
  treeVars.ele_z = new std::vector<float>;
  treeVars.ele_vertexX = new std::vector<float>;
  treeVars.ele_vertexY = new std::vector<float>;
  treeVars.ele_vertexZ = new std::vector<float>;
  treeVars.ele_mass = new std::vector<float>;
  treeVars.ele_mt = new std::vector<float>;
  treeVars.ele_pdgId = new std::vector<int>;
  treeVars.ele_nDau = new std::vector<int>;


  chain->Branch("ele_trackPositionAtVtxX", "std::vector<float>", &treeVars.ele_trackPositionAtVtxX);
  chain->Branch("ele_trackPositionAtVtxY", "std::vector<float>", &treeVars.ele_trackPositionAtVtxY);
  chain->Branch("ele_trackPositionAtVtxZ", "std::vector<float>", &treeVars.ele_trackPositionAtVtxZ);
  chain->Branch("ele_trackPositionAtCaloX", "std::vector<float>", &treeVars.ele_trackPositionAtCaloX);
  chain->Branch("ele_trackPositionAtCaloY", "std::vector<float>", &treeVars.ele_trackPositionAtCaloY);
  chain->Branch("ele_trackPositionAtCaloZ", "std::vector<float>", &treeVars.ele_trackPositionAtCaloZ);
  chain->Branch("ele_trackMomentumAtVtxX", "std::vector<float>", &treeVars.ele_trackMomentumAtVtxX);
  chain->Branch("ele_trackMomentumAtVtxY", "std::vector<float>", &treeVars.ele_trackMomentumAtVtxY);
  chain->Branch("ele_trackMomentumAtVtxZ", "std::vector<float>", &treeVars.ele_trackMomentumAtVtxZ);
  chain->Branch("ele_trackMomentumAtCaloX", "std::vector<float>", &treeVars.ele_trackMomentumAtCaloX);
  chain->Branch("ele_trackMomentumAtCaloY", "std::vector<float>", &treeVars.ele_trackMomentumAtCaloY);
  chain->Branch("ele_trackMomentumAtCaloZ", "std::vector<float>", &treeVars.ele_trackMomentumAtCaloZ);
  chain->Branch("ele_VtxX", "std::vector<float>", &treeVars.ele_VtxX);
  chain->Branch("ele_VtxY", "std::vector<float>", &treeVars.ele_VtxY);
  chain->Branch("ele_VtxZ", "std::vector<float>", &treeVars.ele_VtxZ);

  chain->Branch("ele_ScEcal", "std::vector<float>", &treeVars.ele_ScEcal);
  chain->Branch("ele_ScEraw", "std::vector<float>", &treeVars.ele_ScEraw);
  chain->Branch("ele_ScEta", "std::vector<float>", &treeVars.ele_ScEta);
  chain->Branch("ele_ScPhi", "std::vector<float>", &treeVars.ele_ScPhi);
  chain->Branch("ele_nClu", "std::vector<int>", &treeVars.ele_nClu);

  chain->Branch("ele_covEtaEta", "std::vector<float>", &treeVars.ele_covEtaEta);
  chain->Branch("ele_coviEtaiEta", "std::vector<float>", &treeVars.ele_coviEtaiEta);
  chain->Branch("ele_e1x5", "std::vector<float>", &treeVars.ele_e1x5);
  chain->Branch("ele_e2x5Max", "std::vector<float>", &treeVars.ele_e2x5Max);
  chain->Branch("ele_e5x5", "std::vector<float>", &treeVars.ele_e5x5);
  chain->Branch("ele_HoE", "std::vector<float>", &treeVars.ele_HoE);

  chain->Branch("ele_ScEoP", "std::vector<float>", &treeVars.ele_ScEoP);
  chain->Branch("ele_SeedEoP", "std::vector<float>", &treeVars.ele_SeedEoP);
  chain->Branch("ele_eEleClusteroPout", "std::vector<float>", &treeVars.ele_eEleClusteroPout);
  chain->Branch("ele_deltaEtaSuperClusterTrackAtVtx", "std::vector<float>", &treeVars.ele_deltaEtaSuperClusterTrackAtVtx);
  chain->Branch("ele_deltaEtaSeedClusterTrackAtCalo", "std::vector<float>", &treeVars.ele_deltaEtaSeedClusterTrackAtCalo);
  chain->Branch("ele_deltaEtaEleClusterTrackAtCalo", "std::vector<float>", &treeVars.ele_deltaEtaEleClusterTrackAtCalo);
  chain->Branch("ele_deltaPhiSuperClusterTrackAtVtx", "std::vector<float>", &treeVars.ele_deltaPhiSuperClusterTrackAtVtx);
  chain->Branch("ele_deltaPhiSeedClusterTrackAtCalo", "std::vector<float>", &treeVars.ele_deltaPhiSeedClusterTrackAtCalo);
  chain->Branch("ele_deltaPhiEleClusterTrackAtCalo", "std::vector<float>", &treeVars.ele_deltaPhiEleClusterTrackAtCalo);


  chain->Branch("ele_dr03TkSumPt", "std::vector<float>", &treeVars.ele_dr03TkSumPt);
  chain->Branch("ele_dr03EcalRecHitSumEt", "std::vector<float>", &treeVars.ele_dr03EcalRecHitSumEt);
  chain->Branch("ele_dr03HcalDepth1TowerSumEt", "std::vector<float>", &treeVars.ele_dr03HcalDepth1TowerSumEt);
  chain->Branch("ele_dr03HcalDepth2TowerSumEt", "std::vector<float>", &treeVars.ele_dr03HcalDepth2TowerSumEt);
  chain->Branch("ele_dr03HcalTowerSumEt", "std::vector<float>", &treeVars.ele_dr03HcalTowerSumEt);

  chain->Branch("ele_dr04TkSumPt", "std::vector<float>", &treeVars.ele_dr04TkSumPt);
  chain->Branch("ele_dr04EcalRecHitSumEt", "std::vector<float>", &treeVars.ele_dr04EcalRecHitSumEt);
  chain->Branch("ele_dr04HcalDepth1TowerSumEt", "std::vector<float>", &treeVars.ele_dr04HcalDepth1TowerSumEt);
  chain->Branch("ele_dr04HcalDepth2TowerSumEt", "std::vector<float>", &treeVars.ele_dr04HcalDepth2TowerSumEt);
  chain->Branch("ele_dr04HcalTowerSumEt", "std::vector<float>", &treeVars.ele_dr04HcalTowerSumEt);

  chain->Branch("ele_charge", "std::vector<int>", &treeVars.ele_charge);
  chain->Branch("ele_energy", "std::vector<float>", &treeVars.ele_energy);
  chain->Branch("ele_pt", "std::vector<float>", &treeVars.ele_pt);
  chain->Branch("ele_et", "std::vector<float>", &treeVars.ele_et);
  chain->Branch("ele_momentum", "std::vector<float>", &treeVars.ele_momentum);
  chain->Branch("ele_momentumX", "std::vector<float>", &treeVars.ele_momentumX);
  chain->Branch("ele_momentumY", "std::vector<float>", &treeVars.ele_momentumY);
  chain->Branch("ele_momentumZ", "std::vector<float>", &treeVars.ele_momentumZ);
  chain->Branch("ele_theta", "std::vector<float>", &treeVars.ele_theta);
  chain->Branch("ele_eta", "std::vector<float>", &treeVars.ele_eta);
  chain->Branch("ele_phi", "std::vector<float>", &treeVars.ele_phi);  
  chain->Branch("ele_x", "std::vector<float>", &treeVars.ele_x);
  chain->Branch("ele_y", "std::vector<float>", &treeVars.ele_y);
  chain->Branch("ele_z", "std::vector<float>", &treeVars.ele_z);
  chain->Branch("ele_vertexX", "std::vector<float>", &treeVars.ele_vertexX);
  chain->Branch("ele_vertexY", "std::vector<float>", &treeVars.ele_vertexY);
  chain->Branch("ele_vertexZ", "std::vector<float>", &treeVars.ele_vertexZ);
  chain->Branch("ele_mass", "std::vector<float>", &treeVars.ele_mass);
  chain->Branch("ele_mt", "std::vector<float>", &treeVars.ele_mt);
  chain->Branch("ele_pdgId", "std::vector<int>", &treeVars.ele_pdgId);
  chain->Branch("ele_ncand", &treeVars.ele_ncand, "ele_ncand/I");
  chain->Branch("ele_nDau", "std::vector<int>", &treeVars.ele_nDau);

  // Electron ---------- end

  //SClusters ---- start
  treeVars.nBC = new std::vector<int>;
  treeVars.nCrystals = new std::vector<int>;
  treeVars.rawEnergy = new std::vector<float>;
  treeVars.energy = new std::vector<float>;
  treeVars.eta = new std::vector<float>;
  treeVars.phi = new std::vector<float>;

  chain->Branch("nBC", "std::vector<int>", &treeVars.nBC);
  chain->Branch("nCrystals", "std::vector<int>", &treeVars.nCrystals);
  chain->Branch("rawEnergy", "std::vector<float>", &treeVars.rawEnergy);
  chain->Branch("energy", "std::vector<float>", &treeVars.energy);
  chain->Branch("eta", "std::vector<float>", &treeVars.eta);
  chain->Branch("phi", "std::vector<float>", &treeVars.phi);
  //SClusters ---- end


  //Muons --------- start

  //Muon track
  treeVars.muon_pxAtInner = new std::vector<float>;
  treeVars.muon_pyAtInner = new std::vector<float>;
  treeVars.muon_pzAtInner = new std::vector<float>;
  treeVars.muon_xAtInner = new std::vector<float>;
  treeVars.muon_yAtInner = new std::vector<float>;
  treeVars.muon_zAtInner = new std::vector<float>;
  treeVars.muon_pxAtOuter = new std::vector<float>;
  treeVars.muon_pyAtOuter = new std::vector<float>;
  treeVars.muon_pzAtOuter = new std::vector<float>;
  treeVars.muon_xAtOuter = new std::vector<float>;
  treeVars.muon_yAtOuter = new std::vector<float>;
  treeVars.muon_zAtOuter = new std::vector<float>;
  treeVars.muon_TrackVx = new std::vector<float>;
  treeVars.muon_TrackVy = new std::vector<float>;
  treeVars.muon_TrackVz = new std::vector<float>;

  treeVars.muon_isGlobal = new std::vector<int>;
  treeVars.muon_isTracker = new std::vector<int>;
  treeVars.muon_isStandAlone = new std::vector<int>;
  treeVars.muon_isCalo = new std::vector<int>;
  
  treeVars.muon_sumPt03 = new std::vector<float>;
  treeVars.muon_emEt03 = new std::vector<float>;
  treeVars.muon_hadEt03 = new std::vector<float>;
  treeVars.muon_hoEt03 = new std::vector<float>;
  treeVars.muon_nTrk03 = new std::vector<float>;
  treeVars.muon_nJets03 = new std::vector<float>;
  treeVars.muon_sumPt05 = new std::vector<float>;
  treeVars.muon_emEt05 = new std::vector<float>;
  treeVars.muon_hadEt05 = new std::vector<float>;
  treeVars.muon_hoEt05 = new std::vector<float>;
  treeVars.muon_nTrk05 = new std::vector<float>;
  treeVars.muon_nJets05 = new std::vector<float>;

  treeVars.muon_EcalExpDepo = new std::vector<float>;
  treeVars.muon_HcalExpDepo = new std::vector<float>;
  treeVars.muon_HoExpDepo = new std::vector<float>;
  treeVars.muon_emS9 = new std::vector<float>;
  treeVars.muon_hadS9 = new std::vector<float>;
  treeVars.muon_hoS9 = new std::vector<float>;
  treeVars.muon_CaloComp = new std::vector<float>;

  //Muon cand
  treeVars.muon_charge = new std::vector<int>;
  treeVars.muon_energy = new std::vector<float>;
  treeVars.muon_pt = new std::vector<float>;
  treeVars.muon_et = new std::vector<float>;
  treeVars.muon_momentum = new std::vector<float>;
  treeVars.muon_momentumX = new std::vector<float>;
  treeVars.muon_momentumY = new std::vector<float>;
  treeVars.muon_momentumZ = new std::vector<float>;
  treeVars.muon_theta = new std::vector<float>;
  treeVars.muon_eta = new std::vector<float>;
  treeVars.muon_phi = new std::vector<float>;
  treeVars.muon_x = new std::vector<float>;
  treeVars.muon_y = new std::vector<float>;
  treeVars.muon_z = new std::vector<float>;
  treeVars.muon_vertexX = new std::vector<float>;
  treeVars.muon_vertexY = new std::vector<float>;
  treeVars.muon_vertexZ = new std::vector<float>;
  treeVars.muon_mass = new std::vector<float>;
  treeVars.muon_mt = new std::vector<float>;
  treeVars.muon_pdgId = new std::vector<int>;
  treeVars.muon_nDau = new std::vector<int>;

  chain->Branch("muon_pxAtInner", "std::vector<float>", &treeVars.muon_pxAtInner);
  chain->Branch("muon_pyAtInner", "std::vector<float>", &treeVars.muon_pyAtInner);
  chain->Branch("muon_pzAtInner", "std::vector<float>", &treeVars.muon_pzAtInner);
  chain->Branch("muon_xAtInner", "std::vector<float>", &treeVars.muon_xAtInner);
  chain->Branch("muon_yAtInner", "std::vector<float>", &treeVars.muon_yAtInner);
  chain->Branch("muon_zAtInner", "std::vector<float>", &treeVars.muon_zAtInner);
  chain->Branch("muon_pxAtOuter", "std::vector<float>", &treeVars.muon_pxAtOuter);
  chain->Branch("muon_pyAtOuter", "std::vector<float>", &treeVars.muon_pyAtOuter);
  chain->Branch("muon_pzAtOuter", "std::vector<float>", &treeVars.muon_pzAtOuter);
  chain->Branch("muon_xAtOuter", "std::vector<float>", &treeVars.muon_xAtOuter);
  chain->Branch("muon_yAtOuter", "std::vector<float>", &treeVars.muon_yAtOuter);
  chain->Branch("muon_zAtOuter", "std::vector<float>", &treeVars.muon_zAtOuter);
  chain->Branch("muon_TrackVx", "std::vector<float>", &treeVars.muon_TrackVx);
  chain->Branch("muon_TrackVy", "std::vector<float>", &treeVars.muon_TrackVy);
  chain->Branch("muon_TrackVz", "std::vector<float>", &treeVars.muon_TrackVz);

  chain->Branch("muon_isGlobal", "std::vector<int>", &treeVars.muon_isGlobal);
  chain->Branch("muon_isTracker", "std::vector<int>",  &treeVars.muon_isTracker);
  chain->Branch("muon_isStandAlone", "std::vector<int>", &treeVars.muon_isStandAlone);
  chain->Branch("muon_isCalo", "std::vector<int>", &treeVars.muon_isCalo);

  chain->Branch("muon_sumPt03", "std::vector<float>", &treeVars.muon_sumPt03);
  chain->Branch("muon_emEt03", "std::vector<float>", &treeVars.muon_sumPt03);
  chain->Branch("muon_hadEt03", "std::vector<float>", &treeVars.muon_hadEt03);
  chain->Branch("muon_hoEt03", "std::vector<float>", &treeVars.muon_hoEt03);
  chain->Branch("muon_nTrk03", "std::vector<float>", &treeVars.muon_nTrk03);
  chain->Branch("muon_nJets03", "std::vector<float>", &treeVars.muon_nJets03);

  chain->Branch("muon_sumPt05", "std::vector<float>", &treeVars.muon_sumPt05);
  chain->Branch("muon_emEt05", "std::vector<float>", &treeVars.muon_emEt05);
  chain->Branch("muon_hadEt05", "std::vector<float>", &treeVars.muon_hadEt05);
  chain->Branch("muon_hoEt05", "std::vector<float>", &treeVars.muon_hoEt05);
  chain->Branch("muon_nTrk05", "std::vector<float>", &treeVars.muon_nTrk05);
  chain->Branch("muon_nJets05", "std::vector<float>", &treeVars.muon_nJets05);

  chain->Branch("muon_EcalExpDepo", "std::vector<float>", &treeVars.muon_EcalExpDepo);
  chain->Branch("muon_HcalExpDepo", "std::vector<float>", &treeVars.muon_HcalExpDepo);
  chain->Branch("muon_HoExpDepo", "std::vector<float>", &treeVars.muon_HoExpDepo);
  chain->Branch("muon_emS9", "std::vector<float>", &treeVars.muon_emS9);
  chain->Branch("muon_hadS9", "std::vector<float>", &treeVars.muon_hadS9);
  chain->Branch("muon_hoS9", "std::vector<float>", &treeVars.muon_hoS9);
  chain->Branch("muon_CaloComp", "std::vector<float>", &treeVars.muon_CaloComp);


  chain->Branch("muon_charge", "std::vector<int>", &treeVars.muon_charge);
  chain->Branch("muon_energy", "std::vector<float>",  &treeVars.muon_energy);
  chain->Branch("muon_pt",  "std::vector<float>", &treeVars.muon_pt);
  chain->Branch("muon_et",  "std::vector<float>", &treeVars.muon_et);
  chain->Branch("muon_momentum",  "std::vector<float>", &treeVars.muon_momentum);
  chain->Branch("muon_momentumX", "std::vector<float>",  &treeVars.muon_momentumX);
  chain->Branch("muon_momentumY", "std::vector<float>",  &treeVars.muon_momentumY);
  chain->Branch("muon_momentumZ", "std::vector<float>",  &treeVars.muon_momentumZ);
  chain->Branch("muon_theta", "std::vector<float>",  &treeVars.muon_theta);
  chain->Branch("muon_eta", "std::vector<float>",  &treeVars.muon_eta);
  chain->Branch("muon_phi", "std::vector<float>",  &treeVars.muon_phi);  
  chain->Branch("muon_x", "std::vector<float>",  &treeVars.muon_x);
  chain->Branch("muon_y", "std::vector<float>",  &treeVars.muon_y);
  chain->Branch("muon_z", "std::vector<float>",  &treeVars.muon_z);
  chain->Branch("muon_vertexX", "std::vector<float>",  &treeVars.muon_vertexX);
  chain->Branch("muon_vertexY", "std::vector<float>",  &treeVars.muon_vertexY);
  chain->Branch("muon_vertexZ", "std::vector<float>",  &treeVars.muon_vertexZ);
  chain->Branch("muon_mass", "std::vector<float>",  &treeVars.muon_mass);
  chain->Branch("muon_mt", "std::vector<float>",  &treeVars.muon_mt);
  chain->Branch("muon_pdgId", "std::vector<int>",  &treeVars.muon_pdgId);
  chain->Branch("muon_ncand", &treeVars.muon_ncand, "muon_ncand/I");
  chain->Branch("muon_nDau", "std::vector<int>", &treeVars.muon_nDau);

  //Muons ------ end

  //Tracks --- start
  treeVars.track_vtxIndex  = new std::vector<int>;
  treeVars.track_vtxWeight = new std::vector<float>;

  treeVars.track_pxAtInner = new std::vector<float>;
  treeVars.track_pyAtInner = new std::vector<float>;
  treeVars.track_pzAtInner = new std::vector<float>;
  treeVars.track_xAtInner = new std::vector<float>;
  treeVars.track_yAtInner = new std::vector<float>;
  treeVars.track_zAtInner = new std::vector<float>;
  treeVars.track_pxAtOuter = new std::vector<float>;
  treeVars.track_pyAtOuter = new std::vector<float>;
  treeVars.track_pzAtOuter = new std::vector<float>;
  treeVars.track_xAtOuter = new std::vector<float>;
  treeVars.track_yAtOuter = new std::vector<float>;
  treeVars.track_zAtOuter = new std::vector<float>;

  treeVars.track_ValidHits  = new std::vector<float>;
  treeVars.track_LostHits  = new std::vector<float>;
  treeVars.track_NormalizedChi2 = new std::vector<float>;
  treeVars.track_recHitsSize  = new std::vector<float>;

  treeVars.track_Vx = new std::vector<float>;
  treeVars.track_Vy = new std::vector<float>;
  treeVars.track_Vz = new std::vector<float>;

  treeVars.track_Dxy = new std::vector<float>;
  treeVars.track_D0 = new std::vector<float>;
  treeVars.track_Dsz = new std::vector<float>;
  treeVars.track_Dz = new std::vector<float>;

  treeVars.track_DxyError = new std::vector<float>;
  treeVars.track_D0Error = new std::vector<float>;
  treeVars.track_DszError = new std::vector<float>;
  treeVars.track_DzError = new std::vector<float>;

  treeVars.track_DxyPV = new std::vector<float>;
  treeVars.track_DszPV = new std::vector<float>;
  treeVars.track_DzPV = new std::vector<float>;


  //Track Cand                                                                                                               
  treeVars.track_charge = new std::vector<int>;
  treeVars.track_energy = new std::vector<float>;
  treeVars.track_pt = new std::vector<float>;
  treeVars.track_et = new std::vector<float>;
  treeVars.track_momentum = new std::vector<float>;
  treeVars.track_momentumX = new std::vector<float>;
  treeVars.track_momentumY = new std::vector<float>;
  treeVars.track_momentumZ = new std::vector<float>;
  treeVars.track_theta = new std::vector<float>;
  treeVars.track_eta = new std::vector<float>;
  treeVars.track_phi = new std::vector<float>;
  treeVars.track_x = new std::vector<float>;
  treeVars.track_y = new std::vector<float>;
  treeVars.track_z = new std::vector<float>;
  treeVars.track_vertexX = new std::vector<float>;
  treeVars.track_vertexY = new std::vector<float>;
  treeVars.track_vertexZ = new std::vector<float>;
  treeVars.track_mass = new std::vector<float>;
  treeVars.track_mt = new std::vector<float>;
  treeVars.track_pdgId = new std::vector<int>;
  treeVars.track_nDau = new std::vector<int>;

  chain->Branch("track_vtxIndex", "std::vector<int>", &treeVars.track_vtxIndex);
  chain->Branch("track_vtxWeight", "std::vector<float>", &treeVars.track_vtxWeight);

  chain->Branch("track_pxAtInner", "std::vector<float>", &treeVars.track_pxAtInner);
  chain->Branch("track_pyAtInner", "std::vector<float>", &treeVars.track_pyAtInner);
  chain->Branch("track_pzAtInner", "std::vector<float>", &treeVars.track_pzAtInner);
  chain->Branch("track_xAtInner", "std::vector<float>", &treeVars.track_xAtInner);
  chain->Branch("track_yAtInner", "std::vector<float>", &treeVars.track_yAtInner);
  chain->Branch("track_zAtInner", "std::vector<float>", &treeVars.track_zAtInner);
  chain->Branch("track_pxAtOuter", "std::vector<float>", &treeVars.track_pxAtOuter);
  chain->Branch("track_pyAtOuter", "std::vector<float>", &treeVars.track_pyAtOuter);
  chain->Branch("track_pzAtOuter", "std::vector<float>", &treeVars.track_pzAtOuter);
  chain->Branch("track_xAtOuter", "std::vector<float>", &treeVars.track_xAtOuter);
  chain->Branch("track_yAtOuter", "std::vector<float>", &treeVars.track_yAtOuter);
  chain->Branch("track_zAtOuter", "std::vector<float>", &treeVars.track_zAtOuter);

  chain->Branch("track_ValidHits", "std::vector<float>", &treeVars.track_ValidHits);
  chain->Branch("track_LostHits", "std::vector<float>", &treeVars.track_LostHits);
  chain->Branch("track_NormalizedChi2", "std::vector<float>", &treeVars.track_NormalizedChi2);
  chain->Branch("track_recHitsSize", "std::vector<float>", &treeVars.track_recHitsSize);

  chain->Branch("track_Vx", "std::vector<float>", &treeVars.track_Vx);
  chain->Branch("track_Vy", "std::vector<float>", &treeVars.track_Vy);
  chain->Branch("track_Vz", "std::vector<float>", &treeVars.track_Vz);

  chain->Branch("track_Dxy", "std::vector<float>", &treeVars.track_Dxy);
  chain->Branch("track_D0", "std::vector<float>", &treeVars.track_D0);
  chain->Branch("track_Dsz", "std::vector<float>", &treeVars.track_Dsz);
  chain->Branch("track_Dz", "std::vector<float>", &treeVars.track_Dz);
  chain->Branch("track_DxyError", "std::vector<float>", &treeVars.track_DxyError);
  chain->Branch("track_D0Error", "std::vector<float>", &treeVars.track_D0Error);
  chain->Branch("track_DszError", "std::vector<float>", &treeVars.track_DszError);
  chain->Branch("track_DzPV", "std::vector<float>", &treeVars.track_DzPV);

  chain->Branch("track_charge", "std::vector<int>", &treeVars.track_charge);
  chain->Branch("track_energy", "std::vector<float>", &treeVars.track_energy);
  chain->Branch("track_pt", "std::vector<float>", &treeVars.track_pt);
  chain->Branch("track_et", "std::vector<float>", &treeVars.track_et);
  chain->Branch("track_momentum", "std::vector<float>", &treeVars.track_momentum);
  chain->Branch("track_momentumX", "std::vector<float>", &treeVars.track_momentumX);
  chain->Branch("track_momentumY", "std::vector<float>", &treeVars.track_momentumY);
  chain->Branch("track_momentumZ", "std::vector<float>", &treeVars.track_momentumZ);
  chain->Branch("track_theta", "std::vector<float>", &treeVars.track_theta);
  chain->Branch("track_eta", "std::vector<float>", &treeVars.track_eta);
  chain->Branch("track_phi", "std::vector<float>", &treeVars.track_phi);
  chain->Branch("track_x", "std::vector<float>", &treeVars.track_x);
  chain->Branch("track_y", "std::vector<float>", &treeVars.track_y);
  chain->Branch("track_z", "std::vector<float>", &treeVars.track_z);

  chain->Branch("track_vertexX", "std::vector<float>", &treeVars.track_vertexX);
  chain->Branch("track_vertexY", "std::vector<float>", &treeVars.track_vertexY);
  chain->Branch("track_vertexZ", "std::vector<float>", &treeVars.track_vertexZ);
  chain->Branch("track_mass", "std::vector<float>", &treeVars.track_mass);
  chain->Branch("track_mt", "std::vector<float>", &treeVars.track_mt);
  chain->Branch("track_pdgId", "std::vector<int>", &treeVars.track_pdgId);
  chain->Branch("track_ncand", &treeVars.track_ncand, "track_ncand/I");
  chain->Branch("track_nDau", "std::vector<int>", &treeVars.track_nDau);
  //Tracks ---- end                                                                                                          

  //Vertex --- start


  //Vertex
  treeVars.PVx = new std::vector<float>;
  treeVars.PVy = new std::vector<float>;
  treeVars.PVz = new std::vector<float>;

  treeVars.PVErrx = new std::vector<float>;
  treeVars.PVErry = new std::vector<float>;
  treeVars.PVErrz = new std::vector<float>;

  treeVars.SumPt = new std::vector<float>;
  treeVars.ndof = new std::vector<float>;
  treeVars.chi2 = new std::vector<float>;

  chain->Branch("PVx", "std::vector<float>", &treeVars.PVx);
  chain->Branch("PVy", "std::vector<float>", &treeVars.PVy);
  chain->Branch("PVz", "std::vector<float>", &treeVars.PVz);
  chain->Branch("PVErrx", "std::vector<float>", &treeVars.PVErrx);
  chain->Branch("PVErry", "std::vector<float>", &treeVars.PVErry);
  chain->Branch("PVErrz", "std::vector<float>", &treeVars.PVErrz);

  chain->Branch("SumPt", "std::vector<float>", &treeVars.SumPt);
  chain->Branch("ndof", "std::vector<float>", &treeVars.ndof);
  chain->Branch("chi2", "std::vector<float>", &treeVars.chi2);

  //MET
  treeVars.met_charge = new std::vector<int>;
  treeVars.met_energy = new std::vector<float>;
  treeVars.met_et = new std::vector<float>;
  treeVars.met_momentum = new std::vector<float>;
  treeVars.met_theta = new std::vector<float>;
  treeVars.met_eta = new std::vector<float>;
  treeVars.met_phi = new std::vector<float>;
  treeVars.met_x = new std::vector<float>;
  treeVars.met_y = new std::vector<float>;
  treeVars.met_z = new std::vector<float>;
  treeVars.met_vertexX = new std::vector<float>;
  treeVars.met_vertexY = new std::vector<float>;
  treeVars.met_vertexZ = new std::vector<float>;
  treeVars.met_mass = new std::vector<float>;
  treeVars.met_mt = new std::vector<float>;
  treeVars.met_pdgId = new std::vector<int>;
  treeVars.met_nDau = new std::vector<int>;

  chain->Branch("met_charge", "std::vector<int>", &treeVars.met_charge);
  chain->Branch("met_energy", "std::vector<float>", &treeVars.met_energy);
  chain->Branch("met_et", "std::vector<float>", &treeVars.met_et);
  chain->Branch("met_momentum", "std::vector<float>", &treeVars.met_momentum);
  chain->Branch("met_theta", "std::vector<float>", &treeVars.met_theta);
  chain->Branch("met_eta", "std::vector<float>", &treeVars.met_eta);
  chain->Branch("met_phi", "std::vector<float>", &treeVars.met_phi);
  chain->Branch("met_x", "std::vector<float>", &treeVars.met_x);
  chain->Branch("met_y", "std::vector<float>", &treeVars.met_y);
  chain->Branch("met_z", "std::vector<float>", &treeVars.met_z);
  chain->Branch("met_vertexX", "std::vector<float>", &treeVars.met_vertexX);
  chain->Branch("met_vertexY", "std::vector<float>", &treeVars.met_vertexY);
  chain->Branch("met_vertexZ", "std::vector<float>", &treeVars.met_vertexZ);
  chain->Branch("met_mass", "std::vector<float>", &treeVars.met_mass);
  chain->Branch("met_mt", "std::vector<float>", &treeVars.met_mt);
  chain->Branch("met_pdgId", "std::vector<int>", &treeVars.met_pdgId);
  chain->Branch("met_ncand", &treeVars.met_ncand, "met_ncand/I");
  chain->Branch("met_nDau", "std::vector<int>", &treeVars.met_nDau);

  //GENMET
  treeVars.genmet_charge = new std::vector<int>;
  treeVars.genmet_energy = new std::vector<float>;
  treeVars.genmet_et = new std::vector<float>;
  treeVars.genmet_momentum = new std::vector<float>;
  treeVars.genmet_theta = new std::vector<float>;
  treeVars.genmet_eta = new std::vector<float>;
  treeVars.genmet_phi = new std::vector<float>;
  treeVars.genmet_x = new std::vector<float>;
  treeVars.genmet_y = new std::vector<float>;
  treeVars.genmet_z = new std::vector<float>;
  treeVars.genmet_vertexX = new std::vector<float>;
  treeVars.genmet_vertexY = new std::vector<float>;
  treeVars.genmet_vertexZ = new std::vector<float>;
  treeVars.genmet_mass = new std::vector<float>;
  treeVars.genmet_mt = new std::vector<float>;
  treeVars.genmet_pdgId = new std::vector<int>;
  treeVars.genmet_nDau = new std::vector<int>;

  chain->Branch("genmet_charge", "std::vector<int>", &treeVars.genmet_charge);
  chain->Branch("genmet_energy", "std::vector<float>", &treeVars.genmet_energy);
  chain->Branch("genmet_et", "std::vector<float>", &treeVars.genmet_et);
  chain->Branch("genmet_momentum", "std::vector<float>", &treeVars.genmet_momentum);
  chain->Branch("genmet_theta", "std::vector<float>", &treeVars.genmet_theta);
  chain->Branch("genmet_eta", "std::vector<float>", &treeVars.genmet_eta);
  chain->Branch("genmet_phi", "std::vector<float>", &treeVars.genmet_phi);
  chain->Branch("genmet_x", "std::vector<float>", &treeVars.genmet_x);
  chain->Branch("genmet_y", "std::vector<float>", &treeVars.genmet_y);
  chain->Branch("genmet_z", "std::vector<float>", &treeVars.genmet_z);
  chain->Branch("genmet_vertexX", "std::vector<float>", &treeVars.genmet_vertexX);
  chain->Branch("genmet_vertexY", "std::vector<float>", &treeVars.genmet_vertexY);
  chain->Branch("genmet_vertexZ", "std::vector<float>", &treeVars.genmet_vertexZ);
  chain->Branch("genmet_mass", "std::vector<float>", &treeVars.genmet_mass);
  chain->Branch("genmet_mt", "std::vector<float>", &treeVars.genmet_mt);
  chain->Branch("genmet_pdgId", "std::vector<int>", &treeVars.genmet_pdgId);
  chain->Branch("genmet_ncand", &treeVars.genmet_ncand, "genmet_ncand/I");
  chain->Branch("genmet_nDau", "std::vector<int>", &treeVars.genmet_nDau);
}




void initializeBranches(TTree* chain, TreeContent& treeVars)
{
  treeVars.runN = 0;
  treeVars.evtN = 0; 

  treeVars.mcEl_px->clear();
  treeVars.mcEl_py->clear();
  treeVars.mcEl_pz->clear();
  treeVars.mcEl_e->clear();
  treeVars.mcEl_pdgID->clear();
  
  treeVars.mcEventId->clear();
  treeVars.mcWWevent->clear();

  //Preselected
  treeVars.evtPresel = false;

  //Electrons ---- start                                                                                                     
  // ---------- Ele track                                                                                                    
  treeVars.ele_trackPositionAtVtxX->clear();     
  treeVars.ele_trackPositionAtVtxY->clear();     
  treeVars.ele_trackPositionAtVtxZ->clear();     

  treeVars.ele_trackPositionAtCaloX->clear();    
  treeVars.ele_trackPositionAtCaloY->clear();    
  treeVars.ele_trackPositionAtCaloZ->clear();    

  treeVars.ele_trackMomentumAtVtxX->clear();     
  treeVars.ele_trackMomentumAtVtxY->clear();     
  treeVars.ele_trackMomentumAtVtxZ->clear();     

  treeVars.ele_trackMomentumAtCaloX->clear();    
  treeVars.ele_trackMomentumAtCaloY->clear();
  treeVars.ele_trackMomentumAtCaloZ->clear();     

  treeVars.ele_VtxX->clear();  
  treeVars.ele_VtxY->clear();  
  treeVars.ele_VtxZ->clear();  

  // --------  Ele Sc                                                                                                        
  treeVars.ele_ScEcal->clear();   
  treeVars.ele_ScEraw->clear();   
  treeVars.ele_ScEta->clear();    
  treeVars.ele_ScPhi->clear();    
  treeVars.ele_nClu->clear();     

  // ---------- Ele ID                                                                                                       

  treeVars.ele_covEtaEta->clear();   
  treeVars.ele_coviEtaiEta->clear(); 
  treeVars.ele_e1x5->clear();    
  treeVars.ele_e2x5Max->clear(); 
  treeVars.ele_e5x5->clear();    
  treeVars.ele_HoE->clear();     

  treeVars.ele_ScEoP->clear();      
  treeVars.ele_SeedEoP->clear();    
  treeVars.ele_SeedEoPout->clear(); 
  treeVars.ele_eEleClusteroPout->clear();    
  treeVars.ele_deltaEtaSuperClusterTrackAtVtx->clear();      
  treeVars.ele_deltaEtaSeedClusterTrackAtCalo->clear();      
  treeVars.ele_deltaEtaEleClusterTrackAtCalo->clear();       
  treeVars.ele_deltaPhiSuperClusterTrackAtVtx->clear();
  treeVars.ele_deltaPhiSeedClusterTrackAtCalo->clear();     
  treeVars.ele_deltaPhiEleClusterTrackAtCalo->clear();      

  // ------------ Ele Iso                                                                                                    
  //--03                                                                                                                     
  treeVars.ele_dr03TkSumPt->clear();      
  treeVars.ele_dr03EcalRecHitSumEt->clear();      
  treeVars.ele_dr03HcalDepth1TowerSumEt->clear();     
  treeVars.ele_dr03HcalDepth2TowerSumEt->clear();  
  treeVars.ele_dr03HcalTowerSumEt->clear();    
  //--04                                                                                                                     
  treeVars.ele_dr04TkSumPt->clear();          
  treeVars.ele_dr04EcalRecHitSumEt->clear();     
  treeVars.ele_dr04HcalDepth1TowerSumEt->clear();    
  treeVars.ele_dr04HcalDepth2TowerSumEt->clear();    
  treeVars.ele_dr04HcalTowerSumEt->clear();          

  // ----- Ele Cand                                                                                                          
  treeVars.ele_charge->clear();     
  treeVars.ele_energy->clear();     
  treeVars.ele_pt->clear();         
  treeVars.ele_et->clear();         
  treeVars.ele_momentum->clear();     
  treeVars.ele_momentumX->clear();    
  treeVars.ele_momentumY->clear();    
  treeVars.ele_momentumZ->clear();    
  treeVars.ele_theta->clear();        
  treeVars.ele_eta->clear();          
  treeVars.ele_phi->clear();          
  treeVars.ele_x->clear();        
  treeVars.ele_y->clear();       
  treeVars.ele_z->clear();       
  treeVars.ele_vertexX->clear();     
  treeVars.ele_vertexY->clear();     
  treeVars.ele_vertexZ->clear();     
  treeVars.ele_mass->clear();        
  treeVars.ele_mt->clear();          
  treeVars.ele_pdgId->clear();       
  treeVars.ele_nDau->clear();        
  // Electron ---------- end                                                                                                 

  //SClusters ---- start                                                                                                     
  treeVars.nBC->clear();         
  treeVars.nCrystals->clear();   
  treeVars.rawEnergy->clear();    
  treeVars.energy->clear();       
  treeVars.eta->clear();    
  treeVars.phi->clear();    
  //SClusters ---- end                                                                                                       

  //Muons --------- start                                                                                                    

  //Muon track   
  treeVars.muon_pxAtInner->clear();   
  treeVars.muon_pyAtInner->clear();   
  treeVars.muon_pzAtInner->clear();  
  treeVars.muon_xAtInner->clear();   
  treeVars.muon_yAtInner->clear();   
  treeVars.muon_zAtInner->clear();   
  treeVars.muon_pxAtOuter->clear();  
  treeVars.muon_pyAtOuter->clear();  
  treeVars.muon_pzAtOuter->clear();  
  treeVars.muon_xAtOuter->clear();   
  treeVars.muon_yAtOuter->clear();   
  treeVars.muon_zAtOuter->clear();   
  treeVars.muon_TrackVx->clear();    
  treeVars.muon_TrackVy->clear();    
  treeVars.muon_TrackVz->clear();    

  treeVars.muon_isGlobal->clear();   
  treeVars.muon_isTracker->clear();  
  treeVars.muon_isStandAlone->clear();    
  treeVars.muon_isCalo->clear();     

  treeVars.muon_sumPt03->clear();    
  treeVars.muon_emEt03->clear();     
  treeVars.muon_hadEt03->clear();    
  treeVars.muon_hoEt03->clear();     
  treeVars.muon_nTrk03->clear();     
  treeVars.muon_nJets03->clear();    
  treeVars.muon_sumPt05->clear();    
  treeVars.muon_emEt05->clear();     
  treeVars.muon_hadEt05->clear();    
  treeVars.muon_hoEt05->clear();     
  treeVars.muon_nTrk05->clear();     
  treeVars.muon_nJets05->clear();    

  treeVars.muon_EcalExpDepo->clear();    
  treeVars.muon_HcalExpDepo->clear();    
  treeVars.muon_HoExpDepo->clear();      
  treeVars.muon_emS9->clear();           
  treeVars.muon_hadS9->clear();          
  treeVars.muon_hoS9->clear();      
  treeVars.muon_CaloComp->clear();  

  //Muon cand                                                                                                                
  treeVars.muon_charge->clear();      
  treeVars.muon_energy->clear();      
  treeVars.muon_pt->clear();          
  treeVars.muon_et->clear();      
  treeVars.muon_momentum->clear();    
  treeVars.muon_momentumX->clear();   
  treeVars.muon_momentumY->clear();   
  treeVars.muon_momentumZ->clear();   
  treeVars.muon_theta->clear();       
  treeVars.muon_eta->clear();         
  treeVars.muon_phi->clear();    
  treeVars.muon_x->clear();      
  treeVars.muon_y->clear();      
  treeVars.muon_z->clear();      
  treeVars.muon_vertexX->clear();
  treeVars.muon_vertexY->clear();    
  treeVars.muon_vertexZ->clear();    
  treeVars.muon_mass->clear();       
  treeVars.muon_mt->clear();         
  treeVars.muon_pdgId->clear();      
  treeVars.muon_nDau->clear();       
  //Muons ------ end                        

  //Tracks --- start                                                                                                         
  treeVars.track_vtxIndex ->clear();    
  treeVars.track_vtxWeight->clear();    

  treeVars.track_pxAtInner->clear();    
  treeVars.track_pyAtInner->clear();    
  treeVars.track_pzAtInner->clear();    
  treeVars.track_xAtInner->clear();     
  treeVars.track_yAtInner->clear();     
  treeVars.track_zAtInner->clear();     
  treeVars.track_pxAtOuter->clear();    
  treeVars.track_pyAtOuter->clear();    
  treeVars.track_pzAtOuter->clear();    
  treeVars.track_xAtOuter->clear();     
  treeVars.track_yAtOuter->clear();     
  treeVars.track_zAtOuter->clear();     

  treeVars.track_ValidHits ->clear();   
  treeVars.track_LostHits ->clear();    
  treeVars.track_NormalizedChi2->clear();      
  treeVars.track_recHitsSize ->clear();        

  treeVars.track_Vx->clear();   
  treeVars.track_Vy->clear();   
  treeVars.track_Vz->clear();   

  treeVars.track_Dxy->clear();  
  treeVars.track_D0->clear();   
  treeVars.track_Dsz->clear();  
  treeVars.track_Dz->clear();   

  treeVars.track_DxyError->clear();        
  treeVars.track_D0Error->clear();     
  treeVars.track_DszError->clear();    
  treeVars.track_DzError->clear();     

  treeVars.track_DxyPV->clear();       
  treeVars.track_DszPV->clear();    
  treeVars.track_DzPV->clear();     
  //Track Cand                                                                                                               
  treeVars.track_charge->clear();   
  treeVars.track_energy->clear();   
  treeVars.track_pt->clear();       
  treeVars.track_et->clear();       
  treeVars.track_momentum->clear();     
  treeVars.track_momentumX->clear();    
  treeVars.track_momentumY->clear();    
  treeVars.track_momentumZ->clear();    
  treeVars.track_theta->clear();   
  treeVars.track_eta->clear();     
  treeVars.track_phi->clear();     
  treeVars.track_x->clear();       
  treeVars.track_y->clear();       
  treeVars.track_z->clear();       
  treeVars.track_vertexX->clear();    
  treeVars.track_vertexY->clear();    
  treeVars.track_vertexZ->clear();    

  treeVars.track_mass->clear();       
  treeVars.track_mt->clear();         
  treeVars.track_pdgId->clear();      
  treeVars.track_nDau->clear();       
  //Tracks ---- end                                                                                                          

  //Vertex --- start                                                                                                         
  treeVars.PVx->clear();      
  treeVars.PVy->clear();      
  treeVars.PVz->clear();      

  treeVars.PVErrx->clear();     
  treeVars.PVErry->clear();     
  treeVars.PVErrz->clear();     

  treeVars.SumPt->clear();      
  treeVars.ndof->clear();       
  treeVars.chi2->clear();       
  //Vertex --- end                                                                                                           

  // missing transvers energy ---- start                                                                                     

  //MET                                                                                                                      
  treeVars.met_charge->clear();       
  treeVars.met_energy->clear();       
  treeVars.met_et->clear();           
  treeVars.met_momentum->clear();     
  treeVars.met_theta->clear();        
  treeVars.met_eta->clear();          
  treeVars.met_phi->clear();      
  treeVars.met_x->clear();        
  treeVars.met_y->clear();        
  treeVars.met_z->clear();        
  treeVars.met_vertexX->clear();     
  treeVars.met_vertexY->clear();     
  treeVars.met_vertexZ->clear();     
  treeVars.met_mass->clear();      
  treeVars.met_mt->clear();        
  treeVars.met_pdgId->clear();     
  treeVars.met_nDau->clear();      


  //GENMET                                                                                                                   
  treeVars.genmet_charge->clear();      
  treeVars.genmet_energy->clear();      
  treeVars.genmet_et->clear();        
  treeVars.genmet_momentum->clear();     
  treeVars.genmet_theta->clear();      
  treeVars.genmet_eta->clear();    
  treeVars.genmet_phi->clear();     
  treeVars.genmet_x->clear();      
  treeVars.genmet_y->clear();     
  treeVars.genmet_z->clear();    
  treeVars.genmet_vertexX->clear();       
  treeVars.genmet_vertexY->clear();       
  treeVars.genmet_vertexZ->clear();    
  treeVars.genmet_mass->clear();       
  treeVars.genmet_mt->clear();      
  treeVars.genmet_pdgId->clear();  
  treeVars.genmet_nDau->clear();  

  // missing transvers energy ---- end
}

//  treeVars.muon_ncand = 0;
