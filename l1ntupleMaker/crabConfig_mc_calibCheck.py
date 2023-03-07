
'''
crab submit -c cragConfig.py

crab status 
'''




import CRABClient
from CRABClient.UserUtilities import config

config = config()

#config.General.requestName = 
config.General.workArea = 'L1TNtuple_forL1JetL2Calib_12_6_0_pre1_JECLUT2022v5' 
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'l1ntuple_maker_run3_mc_check_calib.py' #'l1ntuple_maker_run3_mc_calib.py'
#config.JobType.pyCfgParams = ['maxEvt=-1', 'prtEvt=10000', 'nVtxMin=50', 'HCALPFA=%s' % (scheme)] 
#config.JobType.outputFiles = ['L1Ntuple_HCAL.root']

config.Data.inputDataset = '/QCD_Pt15to7000_TuneCP5_13p6TeV-pythia8/Run3Winter22DR-L1TPU0to99FEVT_castor_122X_mcRun3_2021_realistic_v9-v2/AODSIM'
config.Data.secondaryInputDataset = '/QCD_Pt15to7000_TuneCP5_13p6TeV-pythia8/Run3Winter22DR-L1TPU0to99FEVT_castor_122X_mcRun3_2021_realistic_v9-v2/GEN-SIM-DIGI-RAW'

config.General.requestName = config.Data.inputDataset.split('/')[2]
#config.Data.outputDatasetTag = 'Run3Summer21DR-FlatPU0to80FEVT_castor_120X_mcRun3_2021_realistic_v6-v1_'
config.Data.outputDatasetTag = '%s_%s' % (config.General.workArea, config.General.requestName)

config.Data.ignoreLocality = False
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic' # 'FileBased' # 'Automatic' #'LumiBased' 
#config.Data.unitsPerJob = 1
#config.Data.totalUnits = 30

#        config.Data.outLFNDirBase = '/store/user/ssawant/'
#config.Data.useParent = True 

#        config.Site.storageSite = 'T2_IN_TIFR' # Choose your site.  

config.Site.storageSite = 'T2_CH_CERN' # Choose your site
config.Data.outLFNDirBase = '/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/'
