# L1T jet layer 2 calibration
Script to make L1T JEC LUTs

## Step 1: Produce L1T ntuples
Latest recipe to make L1T ntuples can be found at [L1 Trigger Emulator Stage 2 Upgrade Instructions](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TStage2Instructions#Environment_Setup_with_Integrati). Choose the latest CMSSW release and l1t-integration tag.

Once CMSSW area is been set up, use cmsDriver command to produce config file mc.py following [recommendation](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TStage2Instructions#Workflows) command for 'Re-emulating Run-3 (MC) with Run3 era for deriving calibrations'. 
This will produce 'mc.py'. Update mc.py to add settings for HCAL PFA1p if it is recommendaded to do so.

An example mc.py used for L1T ntuple production can be found at l1ntupleMaker/l1ntuple_maker_run3_mc_PFA1p.py

Update crab job script l1ntupleMaker/crabConfig.py for 'inputDataset', 'secondaryInputDataset' and 'outLFNDirBase', and submit crab jobs to produce L1T ntuples with the following commands:
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
voms-proxy-init --voms cms
source /cvmfs/cms.cern.ch/crab3/crab.sh

crab submit -c l1ntupleMaker/crabConfig.py
```

Status of the crab jobs can be checked with
```
crab status -d <your crab projects area>
```
Details about crab commands can be found at [Running CMSSW code on the Grid using CRAB3](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial)


## Step 2: Produce BDT input .csv file
The produced L1Ntuples are read and BDT input .csv files are produced by scripts/L1T_HCALL2Calib_stage1.py script. And scripts/script_condor_submit.py is used to submit the jobs on condor on lxplus using the following command.
Set path to the produced L1TNtuples after '--l1ntuple' in scripts/script_condor_submit.py.
```
python scripts/script_condor_submit.py
```

Hadd the output root files and concatenate produced .csv files, if needed.
The .csv file is now ready to use for BDT training.

If L1Ntuples are analyzed by submitting multiple condor jobs on lxplus, each condor job produces one output .csv file, like L1T_HCALL2Calib_stage1_l1NtupleChunkyDonut_PFA1p_nVtxAll_part*_of_60.csv.
The following command concatenate multiple .csv files into a single .csv file:
```
time awk '
    FNR==1 && NR!=1 { while (/^runNumber/) getline; } 
    1 {print}    
' L1T_HCALL2Calib_stage1_l1NtupleChunkyDonut_PFA1p_nVtxAll_part*_of_60.csv   > L1T_Jet_MLInputs_Run3_QCD_Pt15to7000_PFA1p_wHCALL1Run2Scheme_nVtxAll_20220626.csv 
```

## Step 3: BDT trainning
Run l1JetLayer2Calibration_usingBDT/calculate_L1JetSFs_usingBD.ipynb by running l1JetLayer2Calibration_usingBDT/script_condor_gpu_submit.py

Run l1JetLayer2Calibration_usingBDT/compare_JEC_SFs.ipynb to get final JEC_SF.csv file.

Run check_L1JetSFs.ipynb to make quality control plots.

## Step 4: Make LUTs
cd makeLUTs

Run python3 updateSFPtEtaBins.py

Run g++ ex_bitwise_10.cpp -o ex_bitwise_10 && ./ex_bitwise_10




Scripts for l1ntuple making and
for jet shape studies


# git commit:
git remote add origin git@github.com:siddhesh86/l1JetShapesStudies.git
git push -u origin master