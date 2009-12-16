import FWCore.ParameterSet.Config as cms

selectedElectrons = cms.EDFilter("WWElectronSelection",
#                                 src = cms.InputTag("ambiguityResolvedElectrons"),
                                 src = cms.InputTag("gsfElectrons"),
                                 electronIdCutsLabel = cms.InputTag("eidLoose"),
                                 electronEtaMax     = cms.double(2.5),
                                 electronPtMin      = cms.double(5.0),
                                 )

