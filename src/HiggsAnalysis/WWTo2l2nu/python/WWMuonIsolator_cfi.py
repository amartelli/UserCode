import FWCore.ParameterSet.Config as cms

isolatedMuons = cms.EDFilter("WWMuonIsolation",
                             src = cms.InputTag("selectedMuons"),
                             SelectedMuonRefCollectionLabel = cms.InputTag("selectedMuonsRef"),

                             # isolation parameters
                             hcalIsoDepositLabel = cms.InputTag("muGlobalIsoDepositCalByAssociatorTowers","hcal"),
                             ecalIsoDepositLabel = cms.InputTag("muGlobalIsoDepositCalByAssociatorTowers","ecal"),
                             trackerIsoDepositLabel = cms.InputTag("muGlobalIsoDepositCtfTk"),

                             trackIsolCut = cms.double(100.), #10
                             caloIsolCut = cms.double(100.), #10
                             doRefCheck = cms.bool(True),
                             )


