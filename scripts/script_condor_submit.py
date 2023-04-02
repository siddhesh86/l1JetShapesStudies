#!/bin/env python
####!/usr/bin/env python

'''

''' 



from collections import OrderedDict as OD
import json
import os
import subprocess
#import ROOT
#ROOT.PyConfig.IgnoreCommandLineOptions = True
import time
import datetime
import sys
import glob


commands = {
    ### l1ntuples w/ jetL2CalibSyed
    #"PFA1p_Run3_QCD_Pt15to7000_l1NtupleChunkyDonut":  "time python L1T_HCALL2Calib_stage1.py  --HcalPUS PFA1p   --l1ntuple '/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v5_20220311/run_1/l1ntuples_wJetL2Calib_test/L1Ntuple_HCAL_TP_OOT_PUS_PFA1p_wJetL2CalibSyed_nEvts100.root'   --PUrangeTag nVtxAll  --l1MatchGen --l1NtupleChunkyDonut ",
    #"PFA1p_Run3_QCD_Pt15to7000_l1NtupleChunkyDonut":  "time python L1T_HCALL2Calib_stage1.py  --HcalPUS PFA1p   --l1ntuple '/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v5_20220311/run_1/l1ntuples_wJetL2Calib_test/L1Ntuple_HCAL_TP_OOT_PUS_PFA1p_wJetL2CalibSyed_nEvts20k.root'   --PUrangeTag nVtxAll  --l1MatchGen --l1NtupleChunkyDonut ",
    #"PFA1p_Run3_QCD_Pt15to7000_l1NtupleChunkyDonut":  "time python L1T_HCALL2Calib_stage1.py  --HcalPUS PFA1p   --l1ntuple '/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v5_20220311/run_1/l1ntuples_wJetL2Calib_test/L1Ntuple_HCAL_TP_OOT_PUS_PFA1p_wJetL2CalibSyed_nEvts10k_v6p5.root'   --PUrangeTag nVtxAll  --l1MatchGen --l1NtupleChunkyDonut ",
    # 12_3_0_pre6
    #"PFA1p_Run3_QCD_Pt15to7000_l1NtupleChunkyDonut":  "time python L1T_HCALL2Calib_stage1.py  --HcalPUS PFA1p   --l1ntuple '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/QCD_Pt15to7000_TuneCP5_13p6TeV-pythia8/L1TNtuple_forHCalL2Calib_PFA1p_wHCALL1Run2Scheme_Run3Winter22DR-L1TPU0to99FEVT_castor_122X_mcRun3_2021_realistic_v9-v2/220625_085154/00*/L1Ntuple_*.root'   --PUrangeTag nVtxAll  --l1MatchGen --l1NtupleChunkyDonut ",
    #"PFA1p_Run3_QCD_Pt15to7000_l1NtupleChunkyDonut":  "time python L1T_HCALL2Calib_stage1.py  --HcalPUS PFA1p   --l1ntuple '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/QCD_Pt15to7000_TuneCP5_13p6TeV-pythia8/L1TNtuple_forHCalL2Calib_PFA1p_wHCALL1IterativeScheme_Run3Winter22DR-L1TPU0to99FEVT_castor_122X_mcRun3_2021_realistic_v9-v2/220625_091308/00*/L1Ntuple_*.root'   --PUrangeTag nVtxAll  --l1MatchGen --l1NtupleChunkyDonut ",

    # 12_6_0_pre1 
    #"PFA1p_Run3_QCD_Pt15to7000_l1NtupleChunkyDonut":  "time python L1T_HCALL2Calib_stage1.py  --HcalPUS PFA1p   --l1ntuple '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/QCD_Pt15to7000_TuneCP5_13p6TeV-pythia8/L1TNtuple_forL1JetL2Calib_12_6_0_pre1_Run3Winter22DR-L1TPU0to99FEVT_castor_122X_mcRun3_2021_realistic_v9-v2/220922_062812/0000/L1Ntuple_*.root'   --PUrangeTag nVtxAll  --l1MatchGen --l1NtupleChunkyDonut ", # for MLInputs

    #"PFA1p_Run3_QCD_Pt15to7000_l1NtupleChunkyDonut_JEC2022v4_test":  "time python L1T_HCALL2Calib_stage1.py  --HcalPUS PFA1p   --l1ntuple '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/QCD_Pt15to7000_TuneCP5_13p6TeV-pythia8/L1TNtuple_forL1JetL2Calib_12_6_0_pre1_JECLUT2022v4_Run3Winter22DR-L1TPU0to99FEVT_castor_122X_mcRun3_2021_realistic_v9-v2/220928_161945/0000/L1Ntuple_101.root'   --PUrangeTag nVtxAll  --l1MatchGen --l1NtupleChunkyDonut ",
    #"PFA1p_Run3_QCD_Pt15to7000_l1NtupleChunkyDonut_JEC2022v4_test":  "time python L1T_HCALL2Calib_stage1.py  --HcalPUS PFA1p   --l1ntuple '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/QCD_Pt15to7000_TuneCP5_13p6TeV-pythia8/L1TNtuple_forL1JetL2Calib_12_6_0_pre1_JECLUT2022v4_Run3Winter22DR-L1TPU0to99FEVT_castor_122X_mcRun3_2021_realistic_v9-v2/220928_161945/0000/L1Ntuple_101.root'   --PUrangeTag nVtxAll  --l1MatchGen --l1NtupleChunkyDonut ", 
    #"PFA1p_Run3_QCD_Pt15to7000_l1NtupleChunkyDonut_JEC2022v4_test":  "time python L1T_HCALL2Calib_stage1.py  --HcalPUS PFA1p   --l1ntuple '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/QCD_Pt15to7000_TuneCP5_13p6TeV-pythia8/L1TNtuple_forL1JetL2Calib_12_6_0_pre1_JECLUT2022v4_Run3Winter22DR-L1TPU0to99FEVT_castor_122X_mcRun3_2021_realistic_v9-v2/220928_161945/0000/L1Ntuple_*.root'   --PUrangeTag nVtxAll  --l1MatchGen --l1NtupleChunkyDonut ",
    #"PFA1p_Run3_QCD_Pt15to7000_l1NtupleChunkyDonut_JEC2022v5_test":  "time python L1T_HCALL2Calib_stage1.py  --HcalPUS PFA1p   --l1ntuple '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/QCD_Pt15to7000_TuneCP5_13p6TeV-pythia8/L1TNtuple_forL1JetL2Calib_12_6_0_pre1_JECLUT2022v5_Run3Winter22DR-L1TPU0to99FEVT_castor_122X_mcRun3_2021_realistic_v9-v2/221007_181344/0000/L1Ntuple_*.root'   --PUrangeTag nVtxAll  --l1MatchGen --l1NtupleChunkyDonut ",

    #"2022_Default_l1NtupleChunkyDonut_JEC2022v7":  "time python3 L1T_HCALL2Calib_stage1.py  --HcalPUS PFA1p   --l1ntuple '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/Muon/L1TNtuple_13_0_0_pre4_Run2022G-ZMu-PromptReco-v1/230221_161135/0000/L1Ntuple_*.root'   --PUrangeTag nVtxAll  --l1MatchOffline --l1NtupleChunkyDonut ",
    #"2022_Default_l1NtupleChunkyDonut_JEC2022v7":  "time python L1T_HCALL2Calib_stage1.py  --HcalPUS PFA1p   --l1ntuple '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/Muon/L1TNtuple_13_0_0_pre4_Run2022G-ZMu-PromptReco-v1/230221_161135/0000/L1Ntuple_1-1.root'   --PUrangeTag nVtxAll  --l1MatchOffline --l1NtupleChunkyDonut ",
    #"2022_Default_l1NtupleChunkyDonut_JEC2022v7":  "time python3 L1T_HCALL2Calib_stage1.py  --HcalPUS PFA1p   --l1ntuple '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/Muon/L1TNtuple_13_0_0_pre4_Run2022G-ZMu-PromptReco-v1/230221_161135/0000/L1Ntuple_1-1.root'   --PUrangeTag nVtxAll  --l1MatchOffline --l1NtupleChunkyDonut ",

    # 13_1_0_pre2_HBZS0p5
    #"2023_MC_l1NtupleChunkyDonut_JEC2023v0":  "time python3 L1T_HCALL2Calib_stage1.py  --HcalPUS PFA1p   --l1ntuple '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/QCD_PT-15to7000_TuneCP5_13p6TeV_pythia8/L1TNtuple_forL1JetL2Calib_13_1_0_pre2_HBZS0p5_Run3Winter23Digi-FlatPU0to80_126X_mcRun3_2023_forPU65_v1-v1/230328_213104/0000/L1Ntuple_*'   --PUrangeTag nVtxAll  --l1MatchGen  --l1NtupleChunkyDonut    ",
    #"2023_MC_l1NtupleChunkyDonut_JEC2023v0":  "time python3 L1T_HCALL2Calib_stage1.py  --HcalPUS PFA1p   --l1ntuple '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/QCD_PT-15to7000_TuneCP5_13p6TeV_pythia8/L1TNtuple_forL1JetL2Calib_13_1_0_pre2_HBZS0p5_Run3Winter23Digi-FlatPU0to80_126X_mcRun3_2023_forPU65_v1-v1/230328_213104/0000/L1Ntuple_982.root'   --PUrangeTag nVtxAll  --l1MatchGen  --l1NtupleChunkyDonut  ",
    "2022_Data_l1NtupleChunkyDonut_JEC2023v0":  "time python3 L1T_HCALL2Calib_stage1.py  --HcalPUS PFA1p   --l1ntuple '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/Muon/L1TNtuple_forL1JetL2Calib_13_1_0_pre2_HBZS0p5_Run2022G-ZMu-PromptReco-v1/230329_100000/0000/L1Ntuple_*.root'   --PUrangeTag nVtxAll  --l1MatchOffline  --l1NtupleChunkyDonut    ",
    #"2022_Data_l1NtupleChunkyDonut_JEC2023v0":  "time python3 L1T_HCALL2Calib_stage1.py  --HcalPUS PFA1p   --l1ntuple '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/Muon/L1TNtuple_forL1JetL2Calib_13_1_0_pre2_HBZS0p5_Run2022G-ZMu-PromptReco-v1/230329_100000/0000/L1Ntuple_9.root'   --PUrangeTag nVtxAll  --l1MatchOffline  --l1NtupleChunkyDonut    ",
} 

runLocally = False 
nSplitsInput = 1 if runLocally else 50



pwd = os.getcwd()
print("pwd: ",pwd)

for PUS, command0 in commands.items():
    for kSplit in range(nSplitsInput):
        M_quantulesEvts = kSplit
        N_parts         = nSplitsInput
        runName         = "%d_%d" % (N_parts, M_quantulesEvts)

        condor_exec_file = 'condor_exec_L1TStudies_%s_%s.sh' % (PUS, runName) 

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
                '''
                f.write("export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch \n")
                f.write("source /cvmfs/cms.cern.ch/cmsset_default.sh \n\n")                
                f.write("export SCRAM_ARCH=slc7_amd64_gcc900    \n")
                f.write("CMS_VERSION=\"CMSSW_12_3_0_pre1\" \n")
                f.write("scramv1 project CMSSW ${CMS_VERSION} \n")
                f.write("cd ${CMS_VERSION} \n")
                f.write("eval `scram runtime -sh` \n")
                f.write("export X509_USER_PROXY=/afs/cern.ch/user/s/ssawant/x509up_u108989  \n")
                
                #f.write("export PYTHONHOME=/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/python/2.7.14-omkpbe4  \n")
                f.write("export PYTHONHOME=/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/python/2.7.15-bcolbf2   \n")
                
                
                f.write("export PYTHONHOME=/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/python/2.7.15-bcolbf2   \n")

                f.write("echo PYTHONHOME: $PYTHONHOME  \n")
                '''

            f.write("source /cvmfs/cms.cern.ch/cmsset_default.sh  \n")

            '''
            f.write("echo PATH: $PATH  \n")
           
            f.write("export PYTHONHOME=/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/python/2.7.15-bcolbf2   \n")
            f.write("echo PYTHONHOME: $PYTHONHOME  \n")
            f.write("export PYTHON_BASE=/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/python/2.7.15-bcolbf2   \n")
            f.write("echo PYTHON_BASE: $PYTHON_BASE  \n")
            '''
            #f.write("cd /afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/CMSSW_11_2_0/src/   \n")
            #f.write("cd /afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v5_20220311/CMSSW_12_6_0_pre1/CMSSW_12_6_0_pre1/src/    \n")
            '''
            f.write("echo CMSSW_BASE: $CMSSW_BASE  \n")
            f.write("echo PATH: $PATH  \n")
            f.write("echo PYTHONHOME: $PYTHONHOME  \n")
            f.write('export PYTHONPATH="$PYTHONPATH:/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/CMSSW_11_2_0/src:/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/CMSSW_11_2_0/lib/slc7_amd64_gcc900"   \n')
            f.write("echo PYTHONPATH: $PYTHONPATH  \n")
            '''
            
            f.write("cd /afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v6_20230213/CMSSW_13_0_0_pre4/src   \n")
            f.write("pwd   \n")
            f.write("cmsenv   \n")
            f.write('eval `scramv1 runtime -sh` \n') # which is alias for cmsenv
            f.write("echo ROOTSYSY $ROOTSYS   \n")
            f.write("export PATH=$PATH:$ROOTSYS/bin   \n")
            f.write("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib   \n")
            f.write("export PYTHONPATH=$PYTHONPATH:$ROOTSYS/lib   \n")
                                   
            
            #f.write("cd /afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/myStudies/run_1/mimicMLJetRec/Run3_MC/tmp/  \n")
            f.write("cd %s  \n" % (pwd))
            f.write("pwd   \n")
            f.write("eval \n")
            #f.write("%s %d %d  \n" % (command0, N_parts, M_quantulesEvts))
            f.write("%s  --N_parts %d  --M_quantilesIpFilesSet %d  \n" % (command0, N_parts, M_quantulesEvts))


        condor_submit_file = 'condor_submit_L1TStudies_%s_%s.sh' % (PUS, runName)
        with open(condor_submit_file, 'w') as f:
            f.write("universe = vanilla \n")
            f.write("executable = %s \n" % condor_exec_file)
            f.write("getenv = TRUE \n")
            f.write("log = L1TStudies_%s_%s.log \n" % (PUS, runName))
            f.write("output = L1TStudies_%s_%s.out \n" % (PUS, runName))
            f.write("error = L1TStudies_%s_%s.error \n" % (PUS, runName))
            f.write("notification = never \n")
            f.write("should_transfer_files = YES \n")
            f.write("when_to_transfer_output = ON_EXIT \n")
            #f.write("+JobFlavour = \"workday\" \n")
            #f.write("+JobFlavour = \"espresso\" \n") # 20 mins
            #f.write("+JobFlavour = \"microcentury\" \n") # 1 hours
            #f.write("+JobFlavour = \"longlunch\" \n") # 2 hours
            f.write("+JobFlavour = \"workday\" \n") # 8 hours
            #f.write("+JobFlavour = \"tomorrow\" \n") # 1 day
            f.write("queue \n")


        os.system("chmod a+x %s" % condor_exec_file)
        os.system("chmod a+x %s" % condor_submit_file)


        cmd1 = "condor_submit %s" % condor_submit_file
        print("Now:  %s " % cmd1)
        if not runLocally:
            os.system(cmd1)
        else:
            os.system("time source ./%s &" % (condor_exec_file))

os.system("wait ")
print("\n\n\nJobs done\n\n")
