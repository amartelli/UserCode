import FWCore.ParameterSet.Config as cms

selectedMuonsRef = cms.EDFilter("WWMuonSelectionRef",
                                src = cms.InputTag("muons"),
                                muonPtMinEndcap = cms.double(3.0),
                                muonPMinEndcap = cms.double(9.0),
                                muonPtMinBarrel = cms.double(5.0),
                                muonEtaMax = cms.double(2.4),
                                )

