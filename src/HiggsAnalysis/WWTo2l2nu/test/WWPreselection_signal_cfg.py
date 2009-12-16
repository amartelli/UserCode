import FWCore.ParameterSet.Config as cms

process = cms.Process('WWPreselection')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryPilot2_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')

# --- common preselection code ---
#process.load("HiggsAnalysis.WWTo2l2nu.WZMcEventFilter_cfi")
process.load("HiggsAnalysis.WWTo2l2nu.WWPreselectionSequence_cff")
process.McEventFilter.mctype = 'signal'

# --- tree producer ---
process.load("HiggsAnalysis.WWTo2l2nu.WWTreeProducer_cfi")
process.treeProducer.nameFile = 'WZ_3l-ReducedTree_04_12_09.root'
process.treeProducer.nameTree = 'WZAnalysis'
process.treeProducer.mctype = 'signal'
process.treeProducer.dumpPreselInfo = True
process.treeProducer.dumpTracks = False
process.treeProducer.dumpVertices = True



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            debugFlag = cms.untracked.bool(True),
                            debugVebosity = cms.untracked.uint32(10),
                            fileNames = cms.untracked.vstring(

    'file:/tmp/amartell/signal/EC84ED85-E8B1-DE11-A097-00E08134051C.root'
    )
                            )

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("efficiencies_signal.root"),
    closeFileFast = cms.untracked.bool(True)
    )


#process.out = cms.OutputModule("PoolOutputModule",
#                               verbose = cms.untracked.bool(False),
#                               fileName = cms.untracked.string('WZ_3l-ReducedTree_02_12_09.root'),
#                               outputCommands = cms.untracked.vstring(
#    'keep *'
#    )
#                               )



#process.p = cms.Path ( #process.KFactorProducer *
#                       process.WWTo2l2nuPreselectionSequence 
#                      #process.eleIsolationSequence *
#                      # process.ambiguityResolvedElectrons 
#                      #process.trackCandidates 
#			)



process.o = cms.EndPath ( process.WWTo2l2nuPreselectionSequence* process.treeProducer)

