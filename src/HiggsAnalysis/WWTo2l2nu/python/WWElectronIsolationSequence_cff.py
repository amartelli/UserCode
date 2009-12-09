import FWCore.ParameterSet.Config as cms

import RecoEgamma.EgammaIsolationAlgos.electronTrackIsolationScone_cfi
egammaTrackerIsolationLoose = RecoEgamma.EgammaIsolationAlgos.electronTrackIsolationScone_cfi.electronTrackIsolationScone.clone()
egammaTrackerIsolationLoose.intRadius = 0.015
egammaTrackerIsolationLoose.extRadius = 0.4
egammaTrackerIsolationLoose.maxVtxDist = 0.2

WWElectronIsolationSequence = cms.Sequence( egammaTrackerIsolationLoose )


