
'''
crab submit -c cragConfig.py

crab status 
'''




import CRABClient
from CRABClient.UserUtilities import config

config = config()

#config.General.requestName = 
config.General.workArea = 'L1TNtuple_forL1JetL2Calib_13_1_0_pre2_HBZS0p5' 
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'l1ntuple_maker_run3_mc_13_1_0_pre2_HBZS0p5.py' 
#config.JobType.pyCfgParams = ['maxEvt=-1', 'prtEvt=10000', 'nVtxMin=50', 'HCALPFA=%s' % (scheme)] 
#config.JobType.outputFiles = ['L1Ntuple_HCAL.root']

#config.Data.inputDataset = '/QCD_PT-15to7000_TuneCP5_13p6TeV_pythia8/Run3Winter23Digi-FlatPU0to80_126X_mcRun3_2023_forPU65_v1-v1/GEN-SIM-RAW' # '/QCD_Pt15to7000_TuneCP5_13p6TeV-pythia8/Run3Winter22DR-L1TPU0to99FEVT_castor_122X_mcRun3_2021_realistic_v9-v2/GEN-SIM-DIGI-RAW' #'/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/Run3Summer21DR-FlatPU0to80FEVT_castor_120X_mcRun3_2021_realistic_v6-v1/AODSIM'
#config.Data.secondaryInputDataset = '/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/Run3Summer21DR-FlatPU0to80FEVT_castor_120X_mcRun3_2021_realistic_v6-v1/GEN-SIM-DIGI-RAW'

#config.Data.inputDataset = '/SinglePionGun_E0p2to200/Run3Winter23Digi-NoPU_126X_mcRun3_2023_forPU65_v1-v2/GEN-SIM-RAW'
#config.Data.inputDataset = '/SinglePionGun_E200to500/Run3Winter23Digi-NoPU_126X_mcRun3_2023_forPU65_v1-v2/GEN-SIM-RAW'
#config.Data.inputDataset = '/SinglePhoton_Pt-0To200_gun/Run3Winter23Digi-EpsilonPU_126X_mcRun3_2023_forPU65_v1-v1/GEN-SIM-RAW'

config.Data.inputDataset = '/QCD_Pt15to7000_TuneCP5_13p6TeV-pythia8/Run3Winter22DR-L1TPU0to99FEVT_castor_122X_mcRun3_2021_realistic_v9-v2/GEN-SIM-DIGI-RAW'


config.General.requestName = config.Data.inputDataset.split('/')[2]
#config.Data.outputDatasetTag = 'Run3Summer21DR-FlatPU0to80FEVT_castor_120X_mcRun3_2021_realistic_v6-v1_'
config.Data.outputDatasetTag = '%s_%s' % (config.General.workArea, config.General.requestName)

config.Data.ignoreLocality = False
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased' # 'FileBased' # 'Automatic' #'LumiBased' 
config.Data.unitsPerJob = 10
#config.Data.totalUnits = 30

#        config.Data.outLFNDirBase = '/store/user/ssawant/'
#config.Data.useParent = True 

#        config.Site.storageSite = 'T2_IN_TIFR' # Choose your site.  

config.Site.storageSite = 'T2_CH_CERN' # Choose your site
config.Data.outLFNDirBase = '/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/'
