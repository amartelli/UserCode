import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.WWTo2l2nu.WWMetCorrector_cfi import *
from HiggsAnalysis.WWTo2l2nu.WWMuonIsolator_cfi import *
from HiggsAnalysis.WWTo2l2nu.WWMuonIsolatorRef_cfi import *
from HiggsAnalysis.WWTo2l2nu.WWMuonSelector_cfi import *
from HiggsAnalysis.WWTo2l2nu.WWMuonSelectorRef_cfi import *
from HiggsAnalysis.WWTo2l2nu.WWElectronIdSequence_cff import *
from HiggsAnalysis.WWTo2l2nu.WWElectronIsolationSequence_cff import *
from HiggsAnalysis.WWTo2l2nu.WWElectronIsolator_cfi import *
from HiggsAnalysis.WWTo2l2nu.WWElectronIsolatorRef_cfi import *
from HiggsAnalysis.WWTo2l2nu.WWElectronSelector_cfi import *
from HiggsAnalysis.WWTo2l2nu.WWElectronSelectorRef_cfi import *
from HiggsAnalysis.WWTo2l2nu.WWPreselectionMarker_cfi import *

WWTo2l2nuPreselectionSequence = cms.Sequence(
    muonCorrectedMET *
    selectedMuons *
    selectedMuonsRef *
    isolatedMuons *
    isolatedMuonsRef *
    WWElectronIdSequence *
    WWElectronIsolationSequence *
    selectedElectrons *
    selectedElectronsRef *
    isolatedElectrons *
    isolatedElectronsRef * 
    preselectionMarker
    )
    

