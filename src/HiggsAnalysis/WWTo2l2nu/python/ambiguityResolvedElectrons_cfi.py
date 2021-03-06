import FWCore.ParameterSet.Config as cms

ambiguityResolvedElectrons = cms.EDFilter("AmbResolver",
                                          src = cms.InputTag("gsfElectrons"),
                                          filter = cms.bool(False),
                                          reducedElectronsRefCollectionLabel = cms.InputTag("isolatedElectronsRef"),
                                          doRefCheck = cms.bool(True)
                                          )

