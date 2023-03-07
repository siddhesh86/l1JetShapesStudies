
'''
crab submit -c cragConfig.py

crab status -d <dir name>
'''




import CRABClient
from CRABClient.UserUtilities import config

config = config()

#config.General.requestName = 
config.General.workArea = 'L1TNtuple_13_0_0_pre4' 
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'l1ntuple_maker_2022_data.py'
#config.JobType.pyCfgParams = ['maxEvt=-1', 'prtEvt=10000', 'nVtxMin=50', 'HCALPFA=%s' % (scheme)] 
#config.JobType.outputFiles = ['L1Ntuple_HCAL.root']

config.Data.inputDataset = '/Muon/Run2022G-ZMu-PromptReco-v1/RAW-RECO' #'/Muon/Run2022G-PromptReco-v1/AOD'
#config.Data.secondaryInputDataset = '/Muon/Run2022G-v1/RAW' #/Muon/Run2022G-ZMu-PromptReco-v1/RAW-RECO'

config.General.requestName = config.Data.inputDataset.split('/')[2]
config.Data.outputDatasetTag = '%s_%s' % (config.General.workArea, config.General.requestName)

config.Data.ignoreLocality = False
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic' # 'Automatic' #'LumiBased' 'FileBased'
config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions22/Cert_Collisions2022_eraG_362433_362760_Golden.json'
config.Data.runRange = '362695,362696'
config.Data.unitsPerJob = 500
#config.Data.totalUnits = 30

#        config.Data.outLFNDirBase = '/store/user/ssawant/'
#config.Data.useParent = True 

#        config.Site.storageSite = 'T2_IN_TIFR' # Choose your site.  

config.Site.storageSite = 'T2_CH_CERN' # Choose your site
config.Data.outLFNDirBase = '/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/'
