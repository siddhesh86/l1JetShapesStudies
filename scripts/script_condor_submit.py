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
    #"PFA1p_Run3_QCD_Pt15to7000_l1NtupleChunkyDonut":  "time python L1T_HCALL2Calib_stage1.py  --HcalPUS PFA1p   --l1ntuple '<path to produced L1TNtuples. Wildcard character * is supported>'   --PUrangeTag nVtxAll  --l1MatchGen --l1NtupleChunkyDonut ",

    "PFA1p_Run3_QCD_Pt15to7000_l1NtupleChunkyDonut":  "time python L1T_HCALL2Calib_stage1.py  --HcalPUS PFA1p   --l1ntuple '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/L1TNtuple_forHCalL2Calib_PFA1p_PhiRing_Run3Summer21DR-FlatPU0to80FEVT_castor_120X_mcRun3_2021_realistic_v6-v1/220413_162531/00*/L1Ntuple_HCAL_TP_OOTPUS_PFA1p_ITPUS_PhiRing_*.root'   --PUrangeTag nVtxAll  --l1MatchGen --l1NtupleChunkyDonut ",
} 

runLocally = False
nSplitsInput = 1 if runLocally else 60



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
            f.write("cd /afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/CMSSW_11_2_0/src/   \n")
            f.write("pwd   \n")
            f.write("cmsenv   \n")
            f.write('eval `scramv1 runtime -sh` \n') # which is alias for cmsenv

            '''
            f.write("echo CMSSW_BASE: $CMSSW_BASE  \n")
            f.write("echo PATH: $PATH  \n")
            f.write("echo PYTHONHOME: $PYTHONHOME  \n")
            f.write('export PYTHONPATH="$PYTHONPATH:/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/CMSSW_11_2_0/src:/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/CMSSW_11_2_0/lib/slc7_amd64_gcc900"   \n')
            f.write("echo PYTHONPATH: $PYTHONPATH  \n")
            '''
            
            
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
