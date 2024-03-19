
'''
crab submit -c cragConfig.py

crab status -d <dir name>
'''




import CRABClient
from CRABClient.UserUtilities import config

config = config()

#config.General.requestName = 
#config.General.workArea = 'L1TNtuple_JEC2024v0_13_3_0_L1SFLLR20240304'
config.General.workArea = 'L1TNtuple_JEC2024v0_13_3_0_L1SFLLR20240311wZSHF4p5GeV' 
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'l1ntuple_maker_2023_data_13_3_0_L1SFLLR20240304.py'
config.JobType.psetName = 'l1ntuple_maker_2023_data_13_3_0_L1SFLLR20240311wZSHF4p5GeV.py'
#config.JobType.pyCfgParams = ['maxEvt=-1', 'prtEvt=10000', 'nVtxMin=50', 'HCALPFA=%s' % (scheme)] 
#config.JobType.outputFiles = ['L1Ntuple_HCAL.root']

#config.Data.inputDataset = '/Muon0/Run2023B-ZMu-PromptReco-v1/RAW-RECO' # '/Muon/Run2022G-ZMu-PromptReco-v1/RAW-RECO' #'/Muon/Run2022G-PromptReco-v1/AOD'
config.Data.inputDataset = '/Muon0/Run2023D-ZMu-PromptReco-v1/RAW-RECO'
#config.Data.inputDataset = '/Muon0/Run2023D-ZMu-PromptReco-v2/RAW-RECO'
#config.Data.secondaryInputDataset = '/Muon/Run2022G-v1/RAW' #/Muon/Run2022G-ZMu-PromptReco-v1/RAW-RECO'

config.General.requestName = config.Data.inputDataset.split('/')[2]
config.Data.outputDatasetTag = '%s_%s' % (config.General.workArea, config.General.requestName)

config.Data.ignoreLocality = False
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased' # 'Automatic' #'LumiBased' 'FileBased'
config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json' # 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json' #'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions22/Cert_Collisions2022_eraG_362433_362760_Golden.json' 
#config.Data.runRange = '362695,362696'
#config.Data.runRange = '366727-367079' Run2023B
config.Data.runRange = '369844-370790'
config.Data.unitsPerJob = 20
#config.Data.totalUnits = 30

#        config.Data.outLFNDirBase = '/store/user/ssawant/'
#config.Data.useParent = True 

#        config.Site.storageSite = 'T2_IN_TIFR' # Choose your site.  

config.Site.storageSite = 'T2_CH_CERN' # Choose your site
config.Data.outLFNDirBase = '/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/'
