import FWCore.ParameterSet.Config as cms

selectedElectronsRef = cms.EDFilter("WWElectronSelectionRef",
                                    src = cms.InputTag("gsfElectrons"),
                                    electronIdCutsLabel = cms.InputTag("egammaIDCutsLoose"),
                                    electronIdLikelihoodLabel = cms.InputTag("egammaIDLikelihood"),
                                    useCuts = cms.bool(True)
                                    # likelihoodThreshold_ = cms.double(0.5)
                                    )

