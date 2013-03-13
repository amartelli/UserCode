# Auto generated configuration file
# using: 
# Revision: 1.381.2.11 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/Generator/python/DYToLL_M_50_TuneZ2star_8TeV_pythia6_tauola_cff.py -s GEN:ProductionFilterSequence,SIM,DIGI,L1,DIGI2RAW,HLT,RAW2DIGI,RECO --conditions auto:mc --datatier RECOSIM --eventcontent RAWSIM -n 10 --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('reSIMtoAOD')

completeApproach = True
effectiveApproach = True

from SimG4Core.Application.g4SimHits_cfi import g4SimHits
if completeApproach:
    g4SimHits.ECalSD.AgeingWithSlopeLY = cms.untracked.bool(True)
    g4SimHits.ECalSD.InstLuminosity = cms.double(10e34)
    g4SimHits.ECalSD.DelivLuminosity = cms.double(500.)
    effectiveApproach = False

from SimCalorimetry.EcalSimProducers.ecalDigiParameters_cff import *
if not effectiveApproach:
    ecal_digi_parameters.UseLCcorrection = cms.untracked.bool(False)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
#process.load('SimGeneral.MixingModule.mix_2012_Summer_50ns_PoissonOOTPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

if completeApproach:
    process.GlobalTag.toGet = cms.VPSet(
        ## noise D
        cms.PSet(record = cms.string("EcalPedestalsRcd"),
                 tag = cms.string("EcalPedestals_TL500_IL1E34_mc"),
                 connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL")
                 ),
        ## laser D
        cms.PSet(record = cms.string("EcalLaserAPDPNRatiosRcd"),
                 tag = cms.string("EcalLaserAPDPNRatios_TL500_IL1E34_mc"),
                 connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL")
                 ),
        #    cms.PSet(record = cms.string("EcalLaserAPDPNRatiosRefRcd"),
        #             tag = cms.string("EcalLaserAPDPNRatiosRef_mc"),
        #             connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_FROM21X")
        #             ),
        cms.PSet(record = cms.string('EcalLaserAlphasRcd'),
                 tag = cms.string('EcalLaserAlphas_EB_sic1_btcp1_EE_sic1_btcp1'),
                 connect = cms.untracked.string('frontier://FrontierPrep/CMS_COND_ECAL')
                 ),
        #VPT aging
        cms.PSet(record = cms.string('EcalLinearCorrectionsRcd'),
                 tag = cms.string('EcalLinearCorrections_mc'),
                 connect = cms.untracked.string('frontier://FrontierPrep/CMS_COND_ECAL')
                 )
        )
    
if effectiveApproach:
    process.GlobalTag.toGet = cms.VPSet(
        ## noise D
        cms.PSet(record = cms.string("EcalPedestalsRcd"),
                 tag = cms.string("EcalPedestals_TL500_IL1E34_mc"),
                 connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL")
                 ),
        ## laser D
        cms.PSet(record = cms.string("EcalLaserAPDPNRatiosRcd"),
                 tag = cms.string("EcalLaserAPDPNRatios_TL500_IL1E34_mc"),
                 connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL")
                 ),
        #    cms.PSet(record = cms.string("EcalLaserAPDPNRatiosRefRcd"),
        #             tag = cms.string("EcalLaserAPDPNRatiosRef_mc"),
        #             connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_FROM21X")
        #             ),
        cms.PSet(record = cms.string('EcalLaserAlphasRcd'),
                 tag = cms.string('EcalLaserAlphas_EB_sic1_btcp1_EE_sic1_btcp1'),
                 connect = cms.untracked.string('frontier://FrontierPrep/CMS_COND_ECAL')
                 ),
        #VPT aging
        cms.PSet(record = cms.string('EcalLinearCorrectionsRcd'),
                 tag = cms.string('EcalLinearCorrections_mc'),
                 connect = cms.untracked.string('frontier://FrontierPrep/CMS_COND_ECAL')
                 ),
        # IC + constant term smearing
        cms.PSet(record = cms.string('EcalIntercalibConstantsMCRcd'),
                 tag = cms.string('EcalIntercalibConstants_TL500_IL1E34_mc'),
                 connect = cms.untracked.string('frontier://FrontierPrep/CMS_COND_ECAL')
                 )
        )

if not effectiveApproach and not completeApproach:
    process.GlobalTag.toGet = cms.VPSet(
        #VPT aging
        cms.PSet(record = cms.string('EcalLinearCorrectionsRcd'),
                 tag = cms.string('EcalLinearCorrections_mc'),
                 connect = cms.untracked.string('frontier://FrontierPrep/CMS_COND_ECAL')
                 )
        )
    
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5)
)

process.options = cms.untracked.PSet(
)

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'root://eoscms.cern.ch//eos/cms/store/caf/user/amartell/Simulation/GEN_SIM/DYToEE_M20/4E50EFAA-2662-E211-839E-001EC9AA9F13.root'
    )
                            )

# Output definition

if completeApproach:
    process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
                                            outputCommands = process.AODSIMEventContent.outputCommands,
                                            fileName = cms.untracked.string('DYToEE_reSIM_noPU_TL500_CompleteApproach.root'),
                                            )
if effectiveApproach:
    process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
                                            outputCommands = process.AODSIMEventContent.outputCommands,
                                            fileName = cms.untracked.string('DYToEE_reSIM_noPU_TL500_EffectiveApproach.root'),
                                            )

if not effectiveApproach and not completeApproach:
    process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
                                            outputCommands = process.AODSIMEventContent.outputCommands,
                                            fileName = cms.untracked.string('DYToEE_reSIM_noPU_TL500_StandardApproach.root'),
                                            )

# Additional output definition

# Path and EndPath definitions
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.AODSIMoutput_step = cms.EndPath(process.AODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
#process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend([process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.AODSIMoutput_step])

