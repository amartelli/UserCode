import FWCore.ParameterSet.Config as cms

isolatedElectronsRef = cms.EDFilter("WWElectronIsolationRef",
                                    src = cms.InputTag("gsfElectrons"),
                                    #src = cms.InputTag("selectedElectrons"),
                                    SelectedElectronRefCollectionLabel = cms.InputTag("selectedElectronsRef"),
                                    TrackIsolationProducerLabel = cms.InputTag("egammaTrackerIsolationLoose"),
                                    TrackIsolCut = cms.double(0.1),
                                    doRefCheck = cms.bool(True)
                                    )

