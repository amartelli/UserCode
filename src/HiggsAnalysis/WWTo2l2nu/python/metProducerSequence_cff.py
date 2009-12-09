import FWCore.ParameterSet.Config as cms

#from RecoMET.Configuration.GenMETParticles_cff import *
#from RecoMET.Configuration.RecoGenMET_cff import *
from JetMETCorrections.Type1MET.MetType1Corrections_cff import *

from RecoMET.Configuration.RecoPFMET_cff import *
#from RecoMET.METProducers.PFMET_cfi import *

pfMET = cms.EDProducer("PFMET",
                       PFCandidates = cms.InputTag("particleFlow"),
                       verbose = cms.untracked.bool(True)
                       )

metSequence = cms.Sequence (pfMET)
