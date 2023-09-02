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
    #"PFA1p_Run3_QCD_Pt15to7000_l1NtupleChunkyDonut":  "time python L1T_HCALL2Calib_stage1.py  --HcalPUS PFA1p   --l1ntuple '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/QCD_Pt15to7000_TuneCP5_13p6TeV-pythia8/L1TNtuple_forL1JetL2Calib_12_6_0_pre1_Run3Winter22DR-L1TPU0to99FEVT_castor_122X_mcRun3_2021_realistic_v9-v2/220922_062812/0000/L1Ntuple_*.root'   --PUrangeTag nVtxAll  --l1MatchGen --l1NtupleChunkyDonut ",

    # BDT training
    #"ChunkyDonut": "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut --l1MatchGen ",
    #"PhiRing": "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen ",

    # BDT training 13_1_0_pre2_HBZS0p5 
    #"ChunkyDonut_GenEt_0.01":          "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut --l1MatchGen --MLTarget GenEt            --fracOfDataToUse 0.01",
    #"ChunkyDonut_logGenEt_0.01":       "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut --l1MatchGen --MLTarget logGenEt         --fracOfDataToUse 0.01",
    #"ChunkyDonut_GenEtByL1Et_0.01":    "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut --l1MatchGen --MLTarget GenEtByL1Et      --fracOfDataToUse 0.01",
    #"ChunkyDonut_logGenEtByL1Et_0.01": "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut --l1MatchGen --MLTarget logGenEtByL1Et   --fracOfDataToUse 0.01",

    #"PhiRing_GenEt_0.01":              "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget GenEt            --fracOfDataToUse 0.01",
    #"PhiRing_logGenEt_0.01":           "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget logGenEt         --fracOfDataToUse 0.01",
    #"PhiRing_GenEtByL1Et_0.01":        "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget GenEtByL1Et      --fracOfDataToUse 0.01",
    #"PhiRing_logGenEtByL1Et_0.01":     "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget logGenEtByL1Et   --fracOfDataToUse 0.01",
    

    #"ChunkyDonut_GenEt_0.10":          "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut --l1MatchGen --MLTarget GenEt            --fracOfDataToUse 0.10",
    #"ChunkyDonut_logGenEt_0.10":       "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut --l1MatchGen --MLTarget logGenEt         --fracOfDataToUse 0.10",
    #"ChunkyDonut_GenEtByL1Et_0.10":    "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut --l1MatchGen --MLTarget GenEtByL1Et      --fracOfDataToUse 0.10",
    #"ChunkyDonut_logGenEtByL1Et_0.10": "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut --l1MatchGen --MLTarget logGenEtByL1Et   --fracOfDataToUse 0.10",

    #"PhiRing_GenEt_0.10":              "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --PhiRing --l1MatchGen --MLTarget GenEt            --fracOfDataToUse 0.10",
    #"PhiRing_logGenEt_0.10":           "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --PhiRing --l1MatchGen --MLTarget logGenEt         --fracOfDataToUse 0.10",
    #"PhiRing_GenEtByL1Et_0.10":        "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --PhiRing --l1MatchGen --MLTarget GenEtByL1Et      --fracOfDataToUse 0.10",
    #"PhiRing_logGenEtByL1Et_0.10":     "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --PhiRing --l1MatchGen --MLTarget logGenEtByL1Et   --fracOfDataToUse 0.10",
    

    
    #"ChunkyDonut_GenEt_0.50":          "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut --l1MatchGen --MLTarget GenEt            --fracOfDataToUse 0.50",
    #"ChunkyDonut_logGenEt_0.50":       "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut --l1MatchGen --MLTarget logGenEt         --fracOfDataToUse 0.50",
    #"ChunkyDonut_GenEtByL1Et_0.50":    "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut --l1MatchGen --MLTarget GenEtByL1Et      --fracOfDataToUse 0.50",
    #"ChunkyDonut_logGenEtByL1Et_0.50": "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut --l1MatchGen --MLTarget logGenEtByL1Et   --fracOfDataToUse 0.50",

    #"PhiRing_GenEt_0.50":              "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget GenEt            --fracOfDataToUse 0.50",
    #"PhiRing_logGenEt_0.50":           "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget logGenEt         --fracOfDataToUse 0.50",
    #"PhiRing_GenEtByL1Et_0.50":        "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget GenEtByL1Et      --fracOfDataToUse 0.50",
    #"PhiRing_logGenEtByL1Et_0.50":     "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget logGenEtByL1Et   --fracOfDataToUse 0.50",

    

    #"ChunkyDonut_GenEt_1.00":          "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut --l1MatchGen --MLTarget GenEt            --fracOfDataToUse 1.00",
    #"ChunkyDonut_logGenEt_1.00":       "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut --l1MatchGen --MLTarget logGenEt         --fracOfDataToUse 1.00",
    #"ChunkyDonut_GenEtByL1Et_1.00":    "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut --l1MatchGen --MLTarget GenEtByL1Et      --fracOfDataToUse 1.00",
    #"ChunkyDonut_logGenEtByL1Et_1.00": "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut --l1MatchGen --MLTarget logGenEtByL1Et   --fracOfDataToUse 1.00",

    #"PhiRing_GenEt_1.00":              "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget GenEt            --fracOfDataToUse 1.00",
    #"PhiRing_logGenEt_1.00":           "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget logGenEt         --fracOfDataToUse 1.00",
    #"PhiRing_GenEtByL1Et_1.00":        "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget GenEtByL1Et      --fracOfDataToUse 1.00",
    #"PhiRing_logGenEtByL1Et_1.00":     "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget logGenEtByL1Et   --fracOfDataToUse 1.00",
    
    # 13_1 : Layer1JSFs from Olivier_v2
    "ChunkyDonut_GenEt_1.00":              "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut     --l1MatchOffline --MLTarget GenEt            --fracOfDataToUse 1.00",
    "ChunkyDonut_logGenEt_1.00":           "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut     --l1MatchOffline --MLTarget logGenEt         --fracOfDataToUse 1.00",
    "ChunkyDonut_GenEtByL1Et_1.00":        "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut     --l1MatchOffline --MLTarget GenEtByL1Et      --fracOfDataToUse 1.00",
    "ChunkyDonut_logGenEtByL1Et_1.00":     "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --ChunkyDonut     --l1MatchOffline --MLTarget logGenEtByL1Et   --fracOfDataToUse 1.00",

    "PhiRing_GenEt_1.00":              "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchOffline --MLTarget GenEt            --fracOfDataToUse 1.00",
    "PhiRing_logGenEt_1.00":           "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchOffline --MLTarget logGenEt         --fracOfDataToUse 1.00",
    "PhiRing_GenEtByL1Et_1.00":        "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchOffline --MLTarget GenEtByL1Et      --fracOfDataToUse 1.00",
    "PhiRing_logGenEtByL1Et_1.00":     "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchOffline --MLTarget logGenEtByL1Et   --fracOfDataToUse 1.00",

} 

runLocally = False 
nSplitsInput = 1 # if runLocally else 1



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
            f.write("echo pwd: \n")
            f.write("pwd \n")
            f.write("cd %s \n" % pwd)
            f.write("export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch \n")
            f.write("export SCRAM_ARCH=slc6_amd64_gcc700  \n")
            f.write("source /cvmfs/cms.cern.ch/cmsset_default.sh \n\n")
            #f.write("cd ")
            f.write("export X509_USER_PROXY=/home/ssawant/x509up_u56558 \n")
            f.write("eval \n")
            f.write("cd %s \n" % (pwd))
            f.write("source /afs/cern.ch/user/s/ssawant/.bashrc \n")
            f.write("which conda \n")
            f.write("time conda env list \n")
            f.write("conda activate ana_htoaa \n")
            f.write("time conda env list \n")

            f.write("time conda list \n")
            f.write("which python3 \n")
            f.write("python3 -V \n")
            #f.write(" \n")
            f.write("conda activate ana_htoaa \n")
            #f.write("time python3 %s/%s  %s \n" % (pwd,sAnalysis, sConfig_to_use))
            #f.write('echo "GPU information:" \n')
            #f.write("nvidia-smi \n")
            
            f.write("%s \n" % (command0))
            


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
            #f.write("request_gpus = 2 \n")
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
        print("Now:  %s " % cmd1)
        if not runLocally:
            os.system(cmd1)
        else:
            os.system("time source ./%s &" % (condor_exec_file))

os.system("wait ")
print("\n\n\nJobs done\n\n")
