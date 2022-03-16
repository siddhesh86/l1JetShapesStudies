# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: l1Ntuple -s RAW2DIGI --python_filename=l1ntuple_maker_run2_2018_data.py -n 100 --no_output --era=Run2_2018 --data --conditions=112X_dataRun2_v7 --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulFromRAWsimHcalTP --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleAODRAWEMUCalo --customise=L1Trigger/Configuration/customiseSettings.L1TSettingsToCaloParams_2018_v1_3 --customise_commands=process.HcalTPGCoderULUT.LUTGenerationMode=cms.bool(False) --filein=root://cms-xrd-global.cern.ch//store/data/Run2018D/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/321/457/00000/FEB94F25-1AA6-E811-BC22-FA163E177DDE.root --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

process = cms.Process('RAW2DIGI',Run2_2018)


# use command-line options
import FWCore.ParameterSet.VarParsing as VP
args = VP.VarParsing('analysis')

args.register('maxEvt',  -1,  VP.VarParsing.multiplicity.singleton, VP.VarParsing.varType.int,   'Number of events to process (-1 for all)')
args.register('prtEvt', 10000,  VP.VarParsing.multiplicity.singleton, VP.VarParsing.varType.int,   'Print out every Nth event')
args.register('nVtxMin',   0,  VP.VarParsing.multiplicity.singleton, VP.VarParsing.varType.int,   'Minimum # of reconstructed vertices (pileup)')
args.register('nVtxMax', 999,  VP.VarParsing.multiplicity.singleton, VP.VarParsing.varType.int,   'Maximum # of reconstructed vertices (pileup)')
args.register('HCALPFA', 'PFA2',  VP.VarParsing.multiplicity.singleton, VP.VarParsing.varType.string,   'HCAL Pulse Filtering Algorithm name. PFA2, PFA1p, PFA1 ')

args.parseArguments()
print '\n\nProcessing test_HCAL_OOT_PUS_data.py: maxEvt: %d, prtEvt: %d, nVtsMin: %d, nVtsMax: %d, HCALPFA: %s \n' % (args.maxEvt, args.prtEvt,args.nVtxMin,args.nVtxMax,args.HCALPFA)


# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(args.maxEvt),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)
process.MessageLogger.cerr.FwkReport.reportEvery = args.prtEvt

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/data/Run2018D/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/321/457/00000/FEB94F25-1AA6-E811-BC22-FA163E177DDE.root'),
    secondaryFileNames = cms.untracked.vstring()
)


## HCAL OOT PU subtraction
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/HcalPileupMitigation r_19
if args.HCALPFA == "PFA2":
    print "Setting PFA2 settings"
    process.load("SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff")
    
    process.simHcalTriggerPrimitiveDigis.numberOfFilterPresamplesHBQIE11 = 0
    process.simHcalTriggerPrimitiveDigis.numberOfFilterPresamplesHEQIE11 = 0
    process.simHcalTriggerPrimitiveDigis.weightsQIE11 = {
        "ieta1" :  [1.0, 1.0],
        "ieta2" :  [1.0, 1.0],
        "ieta3" :  [1.0, 1.0],
        "ieta4" :  [1.0, 1.0],
        "ieta5" :  [1.0, 1.0],
        "ieta6" :  [1.0, 1.0],
        "ieta7" :  [1.0, 1.0],
        "ieta8" :  [1.0, 1.0],
        "ieta9" :  [1.0, 1.0],
        "ieta10" : [1.0, 1.0],
        "ieta11" : [1.0, 1.0],
        "ieta12" : [1.0, 1.0],
        "ieta13" : [1.0, 1.0],
        "ieta14" : [1.0, 1.0],
        "ieta15" : [1.0, 1.0],
        "ieta16" : [1.0, 1.0],
        "ieta17" : [1.0, 1.0],
        "ieta18" : [1.0, 1.0],
        "ieta19" : [1.0, 1.0],
        "ieta20" : [1.0, 1.0],
        "ieta21" : [1.0, 1.0],
        "ieta22" : [1.0, 1.0],
        "ieta23" : [1.0, 1.0],
        "ieta24" : [1.0, 1.0],
        "ieta25" : [1.0, 1.0],
        "ieta26" : [1.0, 1.0],
        "ieta27" : [1.0, 1.0],
        "ieta28" : [1.0, 1.0]
    }
elif args.HCALPFA == "PFA1p":
    print "Setting PFA1p settings"    
    process.load("SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff")
    
    process.simHcalTriggerPrimitiveDigis.numberOfFilterPresamplesHBQIE11 = 1
    process.simHcalTriggerPrimitiveDigis.numberOfFilterPresamplesHEQIE11 = 1
    process.simHcalTriggerPrimitiveDigis.weightsQIE11 = {
        "ieta1" :  [-0.47, 1.0],
        "ieta2" :  [-0.47, 1.0],
        "ieta3" :  [-0.47, 1.0],
        "ieta4" :  [-0.47, 1.0],
        "ieta5" :  [-0.47, 1.0],
        "ieta6" :  [-0.47, 1.0],
        "ieta7" :  [-0.47, 1.0],
        "ieta8" :  [-0.47, 1.0],
        "ieta9" :  [-0.47, 1.0],
        "ieta10" : [-0.47, 1.0],
        "ieta11" : [-0.47, 1.0],
        "ieta12" : [-0.47, 1.0],
        "ieta13" : [-0.47, 1.0],
        "ieta14" : [-0.47, 1.0],
        "ieta15" : [-0.47, 1.0],
        "ieta16" : [-0.47, 1.0],
        "ieta17" : [-0.47, 1.0],
        "ieta18" : [-0.47, 1.0],
        "ieta19" : [-0.47, 1.0],
        "ieta20" : [-0.47, 1.0],
        "ieta21" : [-0.43, 1.0],
        "ieta22" : [-0.43, 1.0],
        "ieta23" : [-0.43, 1.0],
        "ieta24" : [-0.43, 1.0],
        "ieta25" : [-0.43, 1.0],
        "ieta26" : [-0.43, 1.0],
        "ieta27" : [-0.43, 1.0],
        "ieta28" : [-0.43, 1.0]
    }
    
    process.HcalTPGCoderULUT.contain1TSHB = True
    process.HcalTPGCoderULUT.contain1TSHE = True
    
    # Pick one of the pairs of lines below based on the intended scenario for running
    process.HcalTPGCoderULUT.containPhaseNSHB = 0.0 # For Run2 2018 Data
    process.HcalTPGCoderULUT.containPhaseNSHE = 0.0 # For Run2 2018 Data
    
elif args.HCALPFA == "PFA1":
    print "Setting PFA1 settings"
    process.load("SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff")
    
    process.simHcalTriggerPrimitiveDigis.numberOfFilterPresamplesHBQIE11 = 1
    process.simHcalTriggerPrimitiveDigis.numberOfFilterPresamplesHEQIE11 = 1
    process.simHcalTriggerPrimitiveDigis.weightsQIE11 = {
        "ieta1" :  [0.0, 1.0],
        "ieta2" :  [0.0, 1.0],
        "ieta3" :  [0.0, 1.0],
        "ieta4" :  [0.0, 1.0],
        "ieta5" :  [0.0, 1.0],
        "ieta6" :  [0.0, 1.0],
        "ieta7" :  [0.0, 1.0],
        "ieta8" :  [0.0, 1.0],
        "ieta9" :  [0.0, 1.0],
        "ieta10" : [0.0, 1.0],
        "ieta11" : [0.0, 1.0],
        "ieta12" : [0.0, 1.0],
        "ieta13" : [0.0, 1.0],
        "ieta14" : [0.0, 1.0],
        "ieta15" : [0.0, 1.0],
        "ieta16" : [0.0, 1.0],
        "ieta17" : [0.0, 1.0],
        "ieta18" : [0.0, 1.0],
        "ieta19" : [0.0, 1.0],
        "ieta20" : [0.0, 1.0],
        "ieta21" : [0.0, 1.0],
        "ieta22" : [0.0, 1.0],
        "ieta23" : [0.0, 1.0],
        "ieta24" : [0.0, 1.0],
        "ieta25" : [0.0, 1.0],
        "ieta26" : [0.0, 1.0],
        "ieta27" : [0.0, 1.0],
        "ieta28" : [0.0, 1.0]
    }
    
    process.HcalTPGCoderULUT.contain1TSHB = True
    process.HcalTPGCoderULUT.contain1TSHE = True
    
    # Pick one of the pairs of lines below based on the intended scenario for running
    process.HcalTPGCoderULUT.containPhaseNSHB = 0.0 # For Run2 2018 Data
    process.HcalTPGCoderULUT.containPhaseNSHE = 0.0 # For Run2 2018 Data


process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(1)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
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
process.GlobalTag = GlobalTag(process.GlobalTag, '112X_dataRun2_v7', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.endjob_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseReEmul
from L1Trigger.Configuration.customiseReEmul import L1TReEmulFromRAWsimHcalTP 

#call to customisation function L1TReEmulFromRAWsimHcalTP imported from L1Trigger.Configuration.customiseReEmul
process = L1TReEmulFromRAWsimHcalTP(process)

# Automatic addition of the customisation function from L1Trigger.L1TNtuples.customiseL1Ntuple
from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleAODRAWEMUCalo 

#call to customisation function L1NtupleAODRAWEMUCalo imported from L1Trigger.L1TNtuples.customiseL1Ntuple
process = L1NtupleAODRAWEMUCalo(process)

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseSettings
from L1Trigger.Configuration.customiseSettings import L1TSettingsToCaloParams_2018_v1_3 

#call to customisation function L1TSettingsToCaloParams_2018_v1_3 imported from L1Trigger.Configuration.customiseSettings
process = L1TSettingsToCaloParams_2018_v1_3(process)


# Add filter on number of reconstructed vertices (pileup)
process.goodVertex = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("isValid & !isFake"),
    filter = cms.bool(True)
)
process.countVertices = cms.EDFilter("VertexCountFilter",
    src = cms.InputTag("goodVertex"),
    minNumber = cms.uint32(args.nVtxMin),
    maxNumber = cms.uint32(args.nVtxMax)
)
process.nGoodVerticesFilterSequence = cms.Sequence(process.goodVertex*process.countVertices)

print("\n# Pre-final L1TReEmul sequence:  ")
print("# {0}".format(process.L1TReEmul))
print("# {0}".format(process.schedule))


## Add primary vertex filter to *EVERY* path in schedule
## Imitating PrefireVetoFilter in L1TNtuples/python/customiseL1Ntuple.py
for path in process.schedule:
    if str(path) == str(process.endjob_step): continue  ## Don't add filter to endjob_step
    path.insert(0,process.nGoodVerticesFilterSequence)

# End of customisation functions


print("\n# Final L1TReEmul sequence:  ")
print("# {0}".format(process.L1TReEmul))
print("# {0}".format(process.schedule))
# for path in process.schedule:
#     print ''
#     print("# {0}".format(path))

# Customisation from command line

process.HcalTPGCoderULUT.LUTGenerationMode=cms.bool(False)
#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion



## Customize output ROOT file name
#subprocess.call(['mkdir', 'output'])
out_name = 'L1Ntuple_HCAL_TP_OOT_PUS'
if args.nVtxMin > 0:   out_name += '_nVtxMin_%d' % args.nVtxMin
if args.nVtxMax < 999: out_name += '_nVtxMax_%d' % args.nVtxMax
#out_name += '_nSamp_%d' % args.samples
#out_name += '_nPresamp_%d' % args.presamps
#out_name += '_HB_' +(str(WGT_HB )[1:-1]).replace(', ', '_').replace('-','n').replace('.','p')
#out_name += '_HE1_'+(str(WGT_HE1)[1:-1]).replace(', ', '_').replace('-','n').replace('.','p')
#out_name += '_HE2_'+(str(WGT_HE2)[1:-1]).replace(', ','_').replace('-','n').replace('.','p')
out_name += '_%s' % args.HCALPFA 
if args.maxEvt == -1:
    out_name += '_all'
else:
    out_name += '_%dk' % (args.maxEvt/1000)
out_name += '.root'
print '\nWill output root file %s' % out_name
 
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(out_name)
)


