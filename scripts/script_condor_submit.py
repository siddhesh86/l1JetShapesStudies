#!/bin/env python
####!/usr/bin/env python

'''

''' 



from collections import OrderedDict as OD
import json
import os
import subprocess
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import time
import datetime
import sys
import glob


runLocally = False

commands = {
    #"PFA2":  "time python L1T_JetMET_Res.py def   /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/TTToSemiLeptonic_TuneCP5_14TeV-powheg-pythia8/L1TNtuple_HCal_OOT_PUS_mc_effi_PFA2_wCaloOption_Run3Winter20DRPremixMiniAOD-110X_mcRun3_2021_realistic_v6-v2/210510_221640/0000/ ",
    #"PFA1p": "time python L1T_JetMET_Res.py PFA1p /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/TTToSemiLeptonic_TuneCP5_14TeV-powheg-pythia8/L1TNtuple_HCal_OOT_PUS_mc_effi_PFA1p_wCaloOption_Run3Winter20DRPremixMiniAOD-110X_mcRun3_2021_realistic_v6-v2/210510_222636/0000/ "

	# Run3 ttbar After wLUTGenFalse_PFA1pRun3ContainPhaseNSm2
    #"PFA2":  "time python L1T_JetMET_Res.py def   /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/TTToSemiLeptonic_TuneCP5_14TeV-powheg-pythia8/L1TNtuple_HCal_OOT_PUS_mc_effi_PFA2_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_Run3Winter20DRPremixMiniAOD-110X_mcRun3_2021_realistic_v6-v2/210729_133650/0000/ ",
	#"PFA1p": "time python L1T_JetMET_Res.py PFA1p /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/TTToSemiLeptonic_TuneCP5_14TeV-powheg-pythia8/L1TNtuple_HCal_OOT_PUS_mc_effi_PFA1p_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_Run3Winter20DRPremixMiniAOD-110X_mcRun3_2021_realistic_v6-v2/210729_133844/0000/ ",


    # 2018 data with larger statistics
    #"PFA2":  "time python L1T_JetMET_Res.py def   /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/ZeroBias/L1TNtuple_HCal_OOT_PUS_data_PFA2_wCaloOption_Run2018D-PromptReco-v2/210707_122906/0000/  ",
    #"PFA1p": "time python L1T_JetMET_Res.py PFA1p /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/ZeroBias/L1TNtuple_HCal_OOT_PUS_data_PFA1p_wCaloOption_Run2018D-PromptReco-v2/210707_122636/0000/  "

    # 2018 SingleMu data
    #"PFA2":  "time python L1T_JetMET_Res.py def   /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA2_wCaloOption_Run2018D-ZMu-12Nov2019_UL2018-v4/210714_114235/0000/  ",
    #"PFA1p": "time python L1T_JetMET_Res.py PFA1p /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA1p_wCaloOption_Run2018D-ZMu-12Nov2019_UL2018-v4/210714_114456/0000/  ",

    # 2018 SingleMu data: 2018 D era
    #"PFA2":  "time python L1T_JetMET_Res.py def   /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA2_wCaloOption_Run2018D-ZMu-12Nov2019_UL2018-v4/210725_101056/*/   ",
    #"PFA1p": "time python L1T_JetMET_Res.py PFA1p /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA1p_wCaloOption_Run2018D-ZMu-12Nov2019_UL2018-v4/210725_100240/*/  ",
    # Reading '*' in directory name has some issue in the python code, which is yet to resolve. So run one directory for a time being.
    #"PFA2":  "time python L1T_JetMET_Res.py def   /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA2_wCaloOption_Run2018D-ZMu-12Nov2019_UL2018-v4/210725_101056/0000/   ",
    #"PFA1p": "time python L1T_JetMET_Res.py PFA1p /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA1p_wCaloOption_Run2018D-ZMu-12Nov2019_UL2018-v4/210725_100240/0000/  ",
	

    # 2018 SingleMu data: 2018 D era. After wLUTGenFalse_PFA1pRun3ContainPhaseNSm2
    #"PFA2":  "time python L1T_JetMET_Res.py def   '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA2_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_Run2018D-ZMu-12Nov2019_UL2018-v4/210729_132736/00*/' ",
    #"PFA1p": "time python L1T_JetMET_Res.py PFA1p '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA1p_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_Run2018D-ZMu-12Nov2019_UL2018-v4/210729_133100/00*/' ",
    #"PFA2":  "time python L1T_JetMET_Res.py def   /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA2_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_Run2018D-ZMu-12Nov2019_UL2018-v4/210729_132736/0000/ ",
    #"PFA1p": "time python L1T_JetMET_Res.py PFA1p /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA1p_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_Run2018D-ZMu-12Nov2019_UL2018-v4/210729_133100/0000/ ",

    # 2018 SingleMu data: 2018 D era. After wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxlt25  nVts < 25
    #"PFA2":  "time python L1T_JetMET_Res.py def   '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA2_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxlt25_Run2018D-ZMu-12Nov2019_UL2018-v4/210809_215415/00*/' ",


    ### 2021/09/14
    # 2018 SingleMu data:  After wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxgt50 nVtx > 50
    "PFA2_2018_SingleMu_nVtxgt50":  "time python L1T_JetMET_Res.py def   '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA2_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxgt50_Run2018*/210914*/00*/'    2018_SingleMu_nVtxgt50 ",
    #"PFA1p": "time python L1T_JetMET_Res.py PFA1p '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA1p_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxgt50_Run2018*/210914*/00*/'   nVtxgt50 ",

    # 2018 SingleMu data:  After wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxlt25 nVtx < 25
    "PFA2_2018_SingleMu_nVtxlt25":  "time python L1T_JetMET_Res.py def   '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA2_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxlt25_Run2018*/210914*/00*/'    2018_SingleMu_nVtxlt25 ",
    #"PFA1p": "time python L1T_JetMET_Res.py PFA1p '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA1p_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxlt25_Run2018*/210914*/00*/'   nVtxlt25 "


    # 2018 ZeroBias data:  After wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxgt50 nVtx > 50
    #"PFA2_2018_ZeroBias_nVtxgt50":  "time python L1T_JetMET_Res.py def   '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/ZeroBias/L1TNtuple_HCal_OOT_PUS_data_PFA2_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxgt50_Run2018*/210914*/00*/'    2018_ZeroBias_nVtxgt50 ",
    #"PFA1p": "time python L1T_JetMET_Res.py PFA1p '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA1p_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxgt50_Run2018*/210914*/00*/'   nVtxgt50 ",

    # 2018 ZeroBias data:  After wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxlt25 nVtx < 25
    #"PFA2_2018_ZeroBias_nVtxlt25":  "time python L1T_JetMET_Res.py def   '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/ZeroBias/L1TNtuple_HCal_OOT_PUS_data_PFA2_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxlt25_Run2018*/210914*/00*/'    2018_ZeroBias_nVtxlt25 "
    #"PFA1p": "time python L1T_JetMET_Res.py PFA1p '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA1p_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxlt25_Run2018*/210914*/00*/'   nVtxlt25 "


    
}

nSplitsInput = 60



pwd = os.getcwd()
print "pwd: ",pwd

for PUS, command0 in commands.items():
    for kSplit in range(nSplitsInput):
        M_quantulesEvts = kSplit
        N_parts         = nSplitsInput
        runName         = "%d_%d" % (N_parts, M_quantulesEvts)

        condor_exec_file = 'condor_exec_HHInference_%s_%s.sh' % (PUS, runName) 

        with open(condor_exec_file, 'w') as f:
            f.write("#!/bin/bash  \n\n")

            if "t3storage3" in pwd:
                f.write("export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch \n")
                f.write("export SCRAM_ARCH=slc6_amd64_gcc700  \n")
                f.write("source /cvmfs/cms.cern.ch/cmsset_default.sh \n\n")
                #f.write("cd ")
                f.write("export X509_USER_PROXY=/home/ssawant/x509up_u56558 \n")
            elif "cern" in pwd:
                f.write("source /afs/cern.ch/user/s/ssawant/.bashrc  \n")

                #f.write("unset PYTHONPATH \n")
                #f.write("unset PYTHONHOME \n")

                #f.write("export PYTHONHOME=/usr/local  \n")

                '''
                # https://twiki.cern.ch/twiki/bin/view/Main/HomerWolfeCMSSWAndGDB : Do you see "ImportError: No module named site"
                Run scram tool info python, you'll see
                ...
                PYTHON_BASE=/afs/cern.ch/cms/slc5_amd64_gcc472/external/python/2.7.3-cms4
                ...
                or whatever. set PYTHONHOME equal to this value. 

                '''
                
                
                '''
                f.write("export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch \n")
                #f.write("export SCRAM_ARCH=slc7_amd64_gcc700   \n")
                f.write("export SCRAM_ARCH=slc7_amd64_gcc900    \n")
                f.write("source /cvmfs/cms.cern.ch/cmsset_default.sh \n\n")
                f.write("export X509_USER_PROXY=/afs/cern.ch/user/s/ssawant/x509up_u108989  \n")
                '''

                f.write("export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch \n")
                f.write("source /cvmfs/cms.cern.ch/cmsset_default.sh \n\n")                
                f.write("export SCRAM_ARCH=slc7_amd64_gcc900    \n")
                f.write("CMS_VERSION=\"CMSSW_11_2_0\" \n")
                f.write("scramv1 project CMSSW ${CMS_VERSION} \n")
                f.write("cd ${CMS_VERSION} \n")
                f.write("eval `scram runtime -sh` \n")
                f.write("export X509_USER_PROXY=/afs/cern.ch/user/s/ssawant/x509up_u108989  \n")
                
                #f.write("export PYTHONHOME=/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/python/2.7.14-omkpbe4  \n")
                f.write("export PYTHONHOME=/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/python/2.7.15-bcolbf2   \n")

                
                f.write("echo PYTHONHOME: $PYTHONHOME  \n")
                

            #f.write("cd /afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/CMSSW_11_2_0/src/   \n")
            #f.write("cmsenv   \n")
            #f.write("cd /afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/myStudies/run_1/mimicMLJetRec/Run3_MC/tmp/  \n")
            f.write("cd %s  \n" % (pwd))
            f.write("eval \n")
            f.write("%s %d %d  \n" % (command0, N_parts, M_quantulesEvts))


        condor_submit_file = 'condor_submit_HHInference_%s_%s.sh' % (PUS, runName)
        with open(condor_submit_file, 'w') as f:
            f.write("universe = vanilla \n")
            f.write("executable = %s \n" % condor_exec_file)
            f.write("getenv = TRUE \n")
            f.write("log = HHInference_%s_%s.log \n" % (PUS, runName))
            f.write("output = HHInference_%s_%s.out \n" % (PUS, runName))
            f.write("error = HHInference_%s_%s.error \n" % (PUS, runName))
            f.write("notification = never \n")
            f.write("should_transfer_files = YES \n")
            f.write("when_to_transfer_output = ON_EXIT \n")
            #f.write("+JobFlavour = \"workday\" \n")
            #f.write("+JobFlavour = \"espresso\" \n") # 20 mins
            #f.write("+JobFlavour = \"microcentury\" \n") # 1 hours
            #f.write("+JobFlavour = \"longlunch\" \n") # 2 hours
            #f.write("+JobFlavour = \"workday\" \n") # 8 hours
            f.write("+JobFlavour = \"tomorrow\" \n") # 1 day
            f.write("queue \n")


        os.system("chmod a+x %s" % condor_exec_file)
        os.system("chmod a+x %s" % condor_submit_file)


        cmd1 = "condor_submit %s" % condor_submit_file
        print "Now:  %s " % cmd1
        if not runLocally:
            os.system(cmd1)
        else:
            os.system("time source ./%s &" % (condor_exec_file))

os.system("wait ")
print "\n\n\nJobs done\n\n"
