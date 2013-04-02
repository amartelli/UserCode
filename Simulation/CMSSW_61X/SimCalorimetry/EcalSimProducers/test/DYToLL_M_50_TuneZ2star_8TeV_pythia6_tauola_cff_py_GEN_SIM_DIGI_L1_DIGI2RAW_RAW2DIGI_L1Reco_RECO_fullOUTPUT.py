# Auto generated configuration file
# using: 
# Revision: 1.11 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Configuration/Generator/python/DYToLL_M_50_TuneZ2star_8TeV_pythia6_tauola_cff.py -s GEN:ProductionFilterSequence,SIM,DIGI,L1,DIGI2RAW,RAW2DIGI,L1Reco,RECO --conditions START61_V11::All --datatier AODSIM --eventcontent AODSIM -n 10 --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

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
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)


## load the DB tags consistent with the scenarios
if completeApproach:
    process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(record = cms.string("EcalPedestalsRcd"),
                 tag = cms.string("EcalPedestals_TL500_IL1E34_mc"),
                 connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL")
                 ),
        cms.PSet(record = cms.string("EcalLaserAPDPNRatiosRcd"),
                 tag = cms.string("EcalLaserAPDPNRatios_TL500_IL1E34_mc"),
                 connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL")
                 ),
        cms.PSet(record = cms.string('EcalLaserAlphasRcd'),
                 tag = cms.string('EcalLaserAlphas_EB_sic1_btcp1_EE_sic1_btcp1'),
                 connect = cms.untracked.string('frontier://FrontierPrep/CMS_COND_ECAL')
                 )
        )
    
if effectiveApproach:
    process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(record = cms.string("EcalPedestalsRcd"),
                 tag = cms.string("EcalPedestals_TL500_IL1E34_mc"),
                 connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL")
                 ),
        cms.PSet(record = cms.string("EcalLaserAPDPNRatiosRcd"),
                 tag = cms.string("EcalLaserAPDPNRatios_TL500_IL1E34_mc"),
                 connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL")
                 ),
        cms.PSet(record = cms.string('EcalLaserAlphasRcd'),
                 tag = cms.string('EcalLaserAlphas_EB_sic1_btcp1_EE_sic1_btcp1'),
                 connect = cms.untracked.string('frontier://FrontierPrep/CMS_COND_ECAL')
                 ),
        cms.PSet(record = cms.string('EcalIntercalibConstantsMCRcd'),
                 tag = cms.string('EcalIntercalibConstants_TL500_IL1E34_mc'),
                 connect = cms.untracked.string('frontier://FrontierPrep/CMS_COND_ECAL')
                 )
        )


# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('Configuration/Generator/python/DYToLL_M_50_TuneZ2star_8TeV_pythia6_tauola_cff.py nevts:10'),
    name = cms.untracked.string('Applications')
)

# Output definition
process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
                                        compressionLevel = cms.untracked.int32(4),
                                        compressionAlgorithm = cms.untracked.string('LZMA'),
                                        eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
                                        outputCommands = process.AODSIMEventContent.outputCommands,
                                        fileName = cms.untracked.string('DYToLL_M_50_TuneZ2star_8TeV_pythia6_tauola_cff_py_standard.root'),
                                        dataset = cms.untracked.PSet(
                                            filterName = cms.untracked.string(''),
                                            dataTier = cms.untracked.string('AODSIM')
                                            ),
                                        SelectEvents = cms.untracked.PSet(
                                            SelectEvents = cms.vstring('generation_step')
                                            )
                                        )


# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'START61_V11::All', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    ExternalDecays = cms.PSet(
        Tauola = cms.untracked.PSet(
            UseTauolaPolarization = cms.bool(True),
            InputCards = cms.PSet(
                mdtau = cms.int32(0),
                pjak2 = cms.int32(0),
                pjak1 = cms.int32(0)
            )
        ),
        parameterSets = cms.vstring('Tauola')
    ),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(8000.0),
    crossSection = cms.untracked.double(762.0),
    UseExternalGenerators = cms.untracked.bool(True),
    PythiaParameters = cms.PSet(
        pythiaUESettings = cms.vstring('MSTU(21)=1     ! Check on possible errors during program execution', 
            'MSTJ(22)=2     ! Decay those unstable particles', 
            'PARJ(71)=10 .  ! for which ctau  10 mm', 
            'MSTP(33)=0     ! no K factors in hard cross sections', 
            'MSTP(2)=1      ! which order running alphaS', 
            'MSTP(51)=10042 ! structure function chosen (external PDF CTEQ6L1)', 
            'MSTP(52)=2     ! work with LHAPDF', 
            'PARP(82)=1.921 ! pt cutoff for multiparton interactions', 
            'PARP(89)=1800. ! sqrts for which PARP82 is set', 
            'PARP(90)=0.227 ! Multiple interactions: rescaling power', 
            'MSTP(95)=6     ! CR (color reconnection parameters)', 
            'PARP(77)=1.016 ! CR', 
            'PARP(78)=0.538 ! CR', 
            'PARP(80)=0.1   ! Prob. colored parton from BBR', 
            'PARP(83)=0.356 ! Multiple interactions: matter distribution parameter', 
            'PARP(84)=0.651 ! Multiple interactions: matter distribution parameter', 
            'PARP(62)=1.025 ! ISR cutoff', 
            'MSTP(91)=1     ! Gaussian primordial kT', 
            'PARP(93)=10.0  ! primordial kT-max', 
            'MSTP(81)=21    ! multiple parton interactions 1 is Pythia default', 
            'MSTP(82)=4     ! Defines the multi-parton model'),
        processParameters = cms.vstring('MSEL=0            !User defined processes', 
            'MSUB(1)=1         !Incl Z0/gamma* production', 
            'MSTP(43)=3        !Both Z0 and gamma*', 
            'MDME(174,1)=0     !Z decay into d dbar', 
            'MDME(175,1)=0     !Z decay into u ubar', 
            'MDME(176,1)=0     !Z decay into s sbar', 
            'MDME(177,1)=0     !Z decay into c cbar', 
            'MDME(178,1)=0     !Z decay into b bbar', 
            'MDME(179,1)=0     !Z decay into t tbar', 
            'MDME(182,1)=1     !Z decay into e- e+', 
            'MDME(183,1)=0     !Z decay into nu_e nu_ebar', 
            'MDME(184,1)=1     !Z decay into mu- mu+', 
            'MDME(185,1)=0     !Z decay into nu_mu nu_mubar', 
            'MDME(186,1)=1     !Z decay into tau- tau+', 
            'MDME(187,1)=0     !Z decay into nu_tau nu_taubar', 
            'CKIN(1)=50.       !Minimum sqrt(s_hat) value (=Z mass)'),
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters')
    )
)


process.ProductionFilterSequence = cms.Sequence(process.generator)


# Path and EndPath definitions
process.generation_step = cms.Path(process.ProductionFilterSequence)
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
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.AODSIMoutput_step)

# Additional output definition
process.AODSIMoutput.outputCommands.append('keep *_rawDataCollector_*_*')
process.AODSIMoutput.outputCommands.append('keep *_simEcalUnsuppressedDigis_*_*')
process.AODSIMoutput.outputCommands.append('keep *_*_ebDigis_*')
process.AODSIMoutput.outputCommands.append('keep *_*_eeDigis_*')
process.AODSIMoutput.outputCommands.append('keep *_*_EcalRecHitsEB_*')
process.AODSIMoutput.outputCommands.append('keep *_*_EcalRecHitsEE_*')
process.AODSIMoutput.outputCommands.append('keep *_*_EcalRecHitsES_*')
process.AODSIMoutput.outputCommands.append('keep *_simEcalPreshowerDigis_*_*')
process.AODSIMoutput.outputCommands.append('keep *_ecalPreshowerDigis_*_*')

if completeApproach:
    process.AODSIMoutput.outputCommands.fileName = cms.untracked.string('DYToLL_M_50_TuneZ2star_8TeV_pythia6_tauola_cff_py_complete.root')
if effectiveApproach:
    process.AODSIMoutput.outputCommands.fileName = cms.untracked.string('DYToLL_M_50_TuneZ2star_8TeV_pythia6_tauola_cff_py_effective.root')



