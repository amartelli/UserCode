import FWCore.ParameterSet.Config as cms

preselectionMarker = cms.EDFilter("WWPreselectionMarker",
                                  ElectronLabel = cms.InputTag("isolatedElectrons"),
                                  MuonLabel = cms.InputTag("isolatedMuons"),
                                  JetLabel = cms.InputTag("iterativeCone5CaloJets"),
                                  CaloMetLabel = cms.InputTag("muonCorrectedMET"),
                                  
                                  LeptonPtMinMin = cms.double(10.),   
                                  LeptonPtMaxMin = cms.double(15.0),
                                  LeptonEtaMax = cms.double(2.5),
                                  LeptonChargeCombination = cms.double(-1),
                                  MetMin = cms.double(0.),  #30
                                  InvMassMin = cms.double(12.0)   
                                  )


