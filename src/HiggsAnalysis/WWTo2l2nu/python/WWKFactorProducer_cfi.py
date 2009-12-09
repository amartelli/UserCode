import FWCore.ParameterSet.Config as cms

KFactorProducer = cms.EDProducer("WWKFactorProducer",
                                 inputFilename = cms.untracked.string('HiggsAnalysis/WWTo2l2nu/data/160_10TeV.dat'),
                                 # PYTHIA Process to apply KFactor
                                 ProcessID = cms.untracked.int32(102),
                                 AltProcessID = cms.untracked.vint32(123, 124),
                                 Debug =cms.untracked.bool(False),
                                 UseNNLO = cms.untracked.bool(False)
                                 )

