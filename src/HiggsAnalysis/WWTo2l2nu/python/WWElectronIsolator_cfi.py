import FWCore.ParameterSet.Config as cms

isolatedElectrons = cms.EDFilter("WWElectronIsolation",
                                 src = cms.InputTag("selectedElectrons"),
                                 TrackIsolCut = cms.double(100.0), # GeV: absolute isolation     ? 20
                                 AbsTrackIsolCut = cms.double(100.0), # GeV: absolute isolation     ? 6
                                 CaloIsolCut = cms.double(100.0), # GeV: absolute isolation     ? 10
                                 ECALIsolCut = cms.double(100.0), # GeV: absolute isolation     ? 10
#                                 absolute = cms.bool(True)
                                 )

