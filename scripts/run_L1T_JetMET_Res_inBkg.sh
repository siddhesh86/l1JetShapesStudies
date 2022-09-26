#!/bin/bash

'''
# Run3 MC
declare -a commands=(
    #"time python L1T_JetMET_Res.py def /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/TTToSemiLeptonic_TuneCP5_14TeV-powheg-pythia8/L1TNtuple_HCal_OOT_PUS_mc_effi_PFA2_wCaloOption_Run3Winter20DRPremixMiniAOD-110X_mcRun3_2021_realistic_v6-v2/210510_221640/0000/ "
    #"time python L1T_JetMET_Res.py PFA1p /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/TTToSemiLeptonic_TuneCP5_14TeV-powheg-pythia8/L1TNtuple_HCal_OOT_PUS_mc_effi_PFA1p_wCaloOption_Run3Winter20DRPremixMiniAOD-110X_mcRun3_2021_realistic_v6-v2/210510_222636/0000/ "

    # Run3 ttbar After wLUTGenFalse_PFA1pRun3ContainPhaseNSm2
    "time python L1T_JetMET_Res.py def   /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/TTToSemiLeptonic_TuneCP5_14TeV-powheg-pythia8/L1TNtuple_HCal_OOT_PUS_mc_effi_PFA2_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_Run3Winter20DRPremixMiniAOD-110X_mcRun3_2021_realistic_v6-v2/210729_133650/0000/ ",
    "time python L1T_JetMET_Res.py PFA1p /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/TTToSemiLeptonic_TuneCP5_14TeV-powheg-pythia8/L1TNtuple_HCal_OOT_PUS_mc_effi_PFA1p_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_Run3Winter20DRPremixMiniAOD-110X_mcRun3_2021_realistic_v6-v2/210729_133844/0000/ ",
)
nSplitsInput=20
'''


# 2018 data
declare -a commands=(
    # ZeroBias data
    #"time python L1T_JetMET_Res.py def   /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/ZeroBias/L1TNtuple_HCal_OOT_PUS_data_PFA2_wCaloOption_Run2018D-PromptReco-v2/210707_122906/0000/  "
    #"time python L1T_JetMET_Res.py PFA1p /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/ZeroBias/L1TNtuple_HCal_OOT_PUS_data_PFA1p_wCaloOption_Run2018D-PromptReco-v2/210707_122636/0000/ "

    # SingleMu data
    #"time python L1T_JetMET_Res.py def   /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA2_wCaloOption_Run2018D-ZMu-12Nov2019_UL2018-v4/210714_114235/0000/  "
    #"time python L1T_JetMET_Res.py PFA1p /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA1p_wCaloOption_Run2018D-ZMu-12Nov2019_UL2018-v4/210714_114456/0000/ "

    # SingleMu data: 2018 D era
    #"time python L1T_JetMET_Res.py def   '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA2_wCaloOption_Run2018D-ZMu-12Nov2019_UL2018-v4/210725_101056/*'  "
    #"time python L1T_JetMET_Res.py PFA1p /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA1p_wCaloOption_Run2018D-ZMu-12Nov2019_UL2018-v4/210725_100240/*/ "
    #"time python L1T_JetMET_Res.py def   '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA2_wCaloOption_Run2018D-ZMu-12Nov2019_UL2018-v4/210725_101056/0000'  "

    # 2018 SingleMu data: 2018 D era. After wLUTGenFalse_PFA1pRun3ContainPhaseNSm2
    #"time python L1T_JetMET_Res.py def   '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA2_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_Run2018D-ZMu-12Nov2019_UL2018-v4/210729_132736/00*/' "
    #"time python L1T_JetMET_Res.py PFA1p '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA1p_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_Run2018D-ZMu-12Nov2019_UL2018-v4/210729_133100/00*/' "
    #"time python L1T_JetMET_Res.py def   /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA2_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_Run2018D-ZMu-12Nov2019_UL2018-v4/210729_132736/0000/ "

    # 2018 SingleMu data: 2018 D era. After wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxlt25  nVts < 25
    #"time python L1T_JetMET_Res.py def   '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA2_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxlt25_Run2018D-ZMu-12Nov2019_UL2018-v4/210809_215415/00*/' "
    

    ### 2021/09/14
    # 2018 SingleMu data: 2018 D era. After wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxgt50 nVtx > 50
    #"time python L1T_JetMET_Res.py def   '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA2_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxgt50_Run2018*/210914*/00*/'   nVtxgt50 "
    #"time python L1T_JetMET_Res.py PFA1p '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA1p_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxgt50_Run2018*/210914*/00*/'  nVtxgt50 "

    # 2018 SingleMu data: 2018 D era. After wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxlt25 nVtx < 25
    #"time python L1T_JetMET_Res.py def   '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA2_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxlt25_Run2018*/210914*/00*/'   nVtxlt25 "
    #"time python L1T_JetMET_Res.py PFA1p '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ssawant/SingleMuon/L1TNtuple_HCal_OOT_PUS_data_PFA1p_wLUTGenFalse_PFA1pRun3ContainPhaseNSm2_nVtxlt25_Run2018*/210914*/00*/'  nVtxlt25 "
    
    
)
nSplitsInput=1




declare -i kCommand=0
for command0 in "${commands[@]}"
do
    #echo "command0: ${command0}"
    
    for ((kSplit=0; kSplit<${nSplitsInput}; kSplit++));
    do
	command="${command0} ${nSplitsInput} ${kSplit}"
	echo "    command: ${command}"
	${command} 2>&1 | tee cout_run_L1T_JetMET_Res_inBkg_${kCommand}_${nSplitsInput}_${kSplit}.txt &
    done

    kCommand=kCommand+1
done

wait
echo "run_L1T_JetMET_Res_inBkg.sh done"
