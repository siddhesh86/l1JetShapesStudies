# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: l1Ntuple -s RAW2DIGI --python_filename=l1ntuple_maker_2024_mc_13_3_0.py -n 100 --no_output --era=Run3 --mc --conditions=133X_mcRun3_2024_realistic_v8 --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulMCFromRAWSimHcalTP --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleRAWEMUGEN_MC --customise=L1Trigger/Configuration/customiseSettings.L1TSettingsToCaloParams_2023_v0_4 --filein=/store/mc/Run3Winter24Digi/SingleNeutrino_Pt-2To20-gun/GEN-SIM-RAW/133X_mcRun3_2024_realistic_v8-v2/2540000/038bda40-23b3-4038-a546-6397626ae3e2.root --no_exec

# Config file from Elena on 11/03/2024
# Layer1 SFs from Elena on 11/03/2024:
# https://indico.cern.ch/event/1392887/contributions/5854592/attachments/2817136/4918395/2024_3_11_L1TDPGNews.pdf#page=7
# HCAL response corrections are applied in deriving layer1 SFs
# https://raw.githubusercontent.com/elenavernazza/CaloL1CalibrationProducer/master/caloParams/caloParams_2023_JAX_ECAL_19_HCAL_51_Phys_newCalib_cfi.py


import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('RAW2DIGI',Run3)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('/store/mc/Run3Winter24Digi/SingleNeutrino_Pt-2To20-gun/GEN-SIM-RAW/133X_mcRun3_2024_realistic_v8-v2/2540000/038bda40-23b3-4038-a546-6397626ae3e2.root'),
    fileNames = cms.untracked.vstring('/store/mc/Run3Winter24Digi/QCD_PT-15to7000_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/FlatPU0to80_133X_mcRun3_2024_realistic_v8-v3/80000/0007e1fa-9d03-4902-a1e7-ececc95e54cd.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    TryToContinue = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToCallForTryToContinue = cms.untracked.vstring(),
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('l1Ntuple nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '133X_mcRun3_2024_realistic_v8', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.endjob_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

## Implement HCAL response correction manually following https://docs.google.com/document/d/1T0ileOTXzM7kgJx0V_dbJzcWr_JNbyi1V5mvssnigAI/edit#heading=h.txprel27aytz
CONDDIR="Debug/HcalDebug/test/conditions"
process.es_prefer = cms.ESPrefer('HcalTextCalibrations','es_ascii')
process.es_ascii = cms.ESSource('HcalTextCalibrations',
   input = cms.VPSet(
      cms.PSet(
          object = cms.string('RespCorrs'),
          #file   = cms.FileInPath('HCALConditions/RespCorrs/HcalRespCorrs_2023_v3.0_data.txt') # cms.FileInPath(CONDDIR+'/RespCorrs/HcalRespCorrs_2023_v3.0_data.txt')
          #file   = cms.FileInPath(CONDDIR+'/RespCorrs/HcalRespCorrs_2023_v3.0_data.txt')
          file   = cms.FileInPath('L1Trigger/L1TCalorimeter/data/HcalRespCorrs_2023_v3.0_data.txt')
      ),
      cms.PSet(
          object = cms.string('Gains'),
          #file   = cms.FileInPath('HCALConditions/Gains/HcalGains_2023_v2.0_data.txt') # cms.FileInPath(CONDDIR+'/Gains/HcalGains_2023_v2.0_data.txt')
          #file   = cms.FileInPath(CONDDIR+'/Gains/HcalGains_2023_v2.0_data.txt')
          file   = cms.FileInPath('L1Trigger/L1TCalorimeter/data/HcalGains_2023_v2.0_data.txt')
      ),
   )
)
## ---------


# Automatic addition of the customisation function from L1Trigger.Configuration.customiseReEmul
from L1Trigger.Configuration.customiseReEmul import L1TReEmulMCFromRAWSimHcalTP 

#call to customisation function L1TReEmulMCFromRAWSimHcalTP imported from L1Trigger.Configuration.customiseReEmul
process = L1TReEmulMCFromRAWSimHcalTP(process)

# Automatic addition of the customisation function from L1Trigger.L1TNtuples.customiseL1Ntuple
from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleRAWEMUGEN_MC 

#call to customisation function L1NtupleRAWEMUGEN_MC imported from L1Trigger.L1TNtuples.customiseL1Ntuple
process = L1NtupleRAWEMUGEN_MC(process)

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseSettings
#from L1Trigger.Configuration.customiseSettings import L1TSettingsToCaloParams_2023_v0_4
#from L1Trigger.Configuration.customiseSettings import L1TSettingsToCaloParams_2023_v0_4_L1SFOlivier20240226
#from L1Trigger.Configuration.customiseSettings import L1TSettingsToCaloParams_2023_v0_4_L1SFLLR20240304
from L1Trigger.Configuration.customiseSettings import L1TSettingsToCaloParams_2023_v0_4_L1SFLLR20240311woZSHF

#call to customisation function L1TSettingsToCaloParams_2023_v0_4 imported from L1Trigger.Configuration.customiseSettings
#process = L1TSettingsToCaloParams_2023_v0_4(process)
#process = L1TSettingsToCaloParams_2023_v0_4_L1SFOlivier20240226(process)
#process = L1TSettingsToCaloParams_2023_v0_4_L1SFLLR20240304(process)
process = L1TSettingsToCaloParams_2023_v0_4_L1SFLLR20240311woZSHF(process)

# End of customisation functions


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
