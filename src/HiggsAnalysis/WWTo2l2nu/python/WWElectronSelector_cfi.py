import FWCore.ParameterSet.Config as cms

selectedElectrons = cms.EDFilter("WWElectronSelection",
                                 src = cms.InputTag("gsfElectrons"),
                                 electronIdCutsLabel = cms.InputTag("egammaIDCutsLoose"),
                                 electronIdLikelihoodLabel = cms.InputTag("egammaIDLikelihood"),
                                 useCuts = cms.bool(True)
#                                 likelihoodThreshold_ = cms.double(-1.)   #0.5
                                 )

