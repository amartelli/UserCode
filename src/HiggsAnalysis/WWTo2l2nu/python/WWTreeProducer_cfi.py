import FWCore.ParameterSet.Config as cms

treeProducer = cms.EDFilter("WWTreeProducer",
                          nameFile = cms.untracked.string('WZ_3l-ReducedTree_prova.root'),
                          nameTree = cms.untracked.string('WZAnalysis'),

                          dumpElectrons = cms.untracked.bool(True),
                          dumpSCs = cms.untracked.bool(True),
                          dumpMuons = cms.untracked.bool(True),
                          dumpTracks = cms.untracked.bool(True),
                          dumpVertices = cms.untracked.bool(True),
                          dumpMet = cms.untracked.bool(True),
                          dumpGenMet = cms.untracked.bool(True),

                          dumpRunInfo = cms.untracked.bool(True),
                          dumpPreselInfo = cms.untracked.bool(True),

                            
                          mctype = cms.untracked.string('signal'),
                          electronCollection = cms.InputTag("ambiguityResolvedElectrons"),
                          ecalBarrelSCCollection = cms.InputTag("correctedHybridSuperClusters"),
                          ecalEndcapSCCollection = cms.InputTag("multi5x5SuperClusters","multi5x5EndcapSuperClusters"),
                          muonCollection = cms.InputTag("isolatedMuons"), # preselection
                          trackCollection = cms.InputTag("trackCandidates"),
                          vertexCollection = cms.InputTag("offlinePrimaryVertices"),
                          metCollection = cms.InputTag("muonCorrectedMET"), # preselection
                          genMetCollection = cms.InputTag("genMetTrue"),    #("genMet"),
                          mcTruthCollection = cms.InputTag("genParticles"),
                            )

