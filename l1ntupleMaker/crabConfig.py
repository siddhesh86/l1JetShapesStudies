
'''
crab submit -c cragConfig.py

crab status 
'''




import CRABClient
from CRABClient.UserUtilities import config

config = config()

#config.General.requestName = 
config.General.workArea = 'L1TNtuple_forHCalL2Calib_PFA1p' 
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'l1ntuple_maker_run3_mc_PFA1p.py'
#config.JobType.pyCfgParams = ['maxEvt=-1', 'prtEvt=10000', 'nVtxMin=50', 'HCALPFA=%s' % (scheme)] 
#config.JobType.outputFiles = ['L1Ntuple_HCAL.root']

config.Data.inputDataset = '/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/Run3Summer21DR-FlatPU0to80FEVT_castor_120X_mcRun3_2021_realistic_v6-v1/AODSIM'
config.Data.secondaryInputDataset = '/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/Run3Summer21DR-FlatPU0to80FEVT_castor_120X_mcRun3_2021_realistic_v6-v1/GEN-SIM-DIGI-RAW'

config.General.requestName = config.Data.inputDataset.split('/')[2]
#config.Data.outputDatasetTag = 'Run3Summer21DR-FlatPU0to80FEVT_castor_120X_mcRun3_2021_realistic_v6-v1_'
config.Data.outputDatasetTag = '%s_%s' % (config.General.workArea, config.General.requestName)

config.Data.ignoreLocality = False
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased' # 'Automatic' #'LumiBased' 
config.Data.unitsPerJob = 1
#config.Data.totalUnits = 30

#        config.Data.outLFNDirBase = '/store/user/ssawant/'
#config.Data.useParent = True 

#        config.Site.storageSite = 'T2_IN_TIFR' # Choose your site.  

config.Site.storageSite = 'T2_CH_CERN' # Choose your site
config.Data.outLFNDirBase = '/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/'
