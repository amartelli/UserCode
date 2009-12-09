import FWCore.ParameterSet.Config as cms

isolatedElectrons = cms.EDFilter("WWElectronIsolation",
                                 src = cms.InputTag("gsfElectrons"),
                                 SelectedElectronRefCollectionLabel = cms.InputTag("selectedElectronsRef"),
                                 TrackIsolationProducerLabel = cms.InputTag("egammaTrackerIsolationLoose"),
                                 TrackIsolCut = cms.double(20.),  #0.1
                                 doRefCheck = cms.bool(True)
                                 )

