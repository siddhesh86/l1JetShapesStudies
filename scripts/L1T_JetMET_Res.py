#!/bin/env python

####! /usr/bin/env python

## *************************************************************************** ##
##  Look at scale and resolution of L1T jets and MET with different input TPs  ##
## *************************************************************************** ##

import os
import sys
import math 
from array import array
import ROOT as R
R.gROOT.SetBatch(False)  ## Don't print histograms to screen while processing

from collections import OrderedDict
import csv
import glob
import json

PRT_EVT  = 1000   ## Print every Nth event
MAX_EVT  = -1     ## Number of events to process per chain
VERBOSE  = False  ## Verbose print-out
PrintLevel = 0
JetClustByHand = True # True ## Run jet clustering by hand
#JetShapes = ['9x9', '8x9', '7x9', '6x9', '5x9', '4x9', '3x9', '3x9_plus_0.5_times_9x9']
JetShapes = ['Default', '9x9', '8x9', '7x9', '6x9', '5x9', '4x9', '3x9', '3x9_plus_0.5_times_9x9']
JetShapesType2 = ['L1TauDefault']
#PUSAlgosAll      = ['Raw', 'RawPUS', 'RawPUS_phiRingMin4', 'RawPUS_phiRingSide4', 'RawPUS_phiRingAdjacent', 'RawPUS_phiDefault']
PUSAlgosAll      = ['L1JDefault', 'Raw', 'RawPUS', 'RawPUS_phiRingMin4', 'RawPUS_phiRingSide4', 'RawPUS_phiRingAdjacent', 'RawPUS_phiDefault']
#PUSAlgosSelected = ['RawPUS', 'RawPUS_phiRingMin4']
PUSAlgosSelected = ['Raw', 'RawPUS', 'RawPUS_phiDefault']
PUSAlgosAllType2 = ['Et', 'RawEt']

runMode = 'CalibJetByHand' # '', 'CalCalibSF', 'CalibJetByHand', 'makeInputForML'
# 'test'           # run on L1Ntuple_*_1.root ntuple for tests
# ''               # 1st round to make jet resolution plots
# 'CalCalibSF'     # set true to fill PFjetPt vs L1jetPt histograms to calculate calibration SFs
# 'CalibJetByHand' # apply CaliSF to JetsByHand
# 'makeInputForML' # write jet information into .csv file to be input to MachineLearning
# 'makePUHisto'    # run quick to make PU histograms
sOutExt = ""

L1TEffiTurnOn_TrigThrshs = [12.0, 35.0, 60.0, 90.0, 120.0, 180.0] # L1T eT threshold for L1T efficiency turn-on curve


''' Command to concatenate files with same heaher. Header has string starting with 'PFJetEtCorr'

time awk '
    FNR==1 && NR!=1 { while (/^PFJetEtCorr/) getline; }
    1 {print}
' L1T_JetMET_Res_def_2018_SingleMu_nVtxgt50_part*_of_60.csv  > L1T_Jet_MLInputs_2018_SingleMu_PFA2_nVtxgt50_20211008.csv 


time awk '
    FNR==1 && NR!=1 { while (/^PFJetEtCorr/) getline; }
    1 {print}
' L1T_JetMET_Res_def_2018_SingleMu_nVtxlt25_part*_of_60.csv  > L1T_Jet_MLInputs_2018_SingleMu_PFA2_nVtxlt25_20211008.csv 

'''

usePUReweighting_0   = True;
#ipFilePUData = "/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/myStudies/run_1/2018_Data/l1analysis_def.root";
#ipFilePUData = "/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/myStudies/run_1/mimicMLJetRec/2018_Data/L1T_JetMET_Res_def_hadded.root"; # larger statistics
#ipFilePUData = "/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/myStudies/run_1/mimicMLJetRec/2018_SingleMu/L1T_JetMET_Res_def_hadded.root"; # larger statistics
ipFilePUData = "/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/myStudies/run_2/2018_SingleMu/IMPORTANT/L1T_JetMET_Res_def_hadded.root"
sHistoPUData = "nVtx";
ipFilePUMC   = "/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/myStudies/run_2/Run3_MC/IMPORTANT/L1T_JetMET_Res_def_hadded_v0.root"
sHistoPUMC   = "nVtx";

nPV_median = 55 # make avgLayer2SF plot for 50 < nPV <= nPV_median and nPV > nPV_median 

PURangeUsedForCalibSF = "PU50to100" # "PU50to100", "PU1to25"
PtRangeForCalibSF = "Pt25To35" #  medPt, lowPt
#sipFileCalibSF = "L1T_JetMET_Layer2CalibSF_$OOTPUS.root";
#sHistoCalibSF = "IEta$ETA/h_jet_byHand_L1JetPt_vs_PFJetPt$JETSHAPE_$PUSALGORITHM_$ETA_PtAllBins_0_$L1MODE_fixBinWidthPFJetPt_ProfileAlongL1JetPt";
#sipFileCalibSF = "/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/myStudies/run_2/2018_SingleMu/IMPORTANT/L1TJetEtSF_2018_SingleMu_20210925.root";
#sipFileCalibSF = "/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/myStudies/run_2/2018_SingleMu/IMPORTANT/L1TJetEtSF_2018_SingleMu_20211002.root";
#sipFileCalibSF = "/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/myStudies/run_3/2018_SingleMu/IMPORTANT/L1TJetEtSF_2018_SingleMu_20211124.root";
sipFileCalibSF = "/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/myStudies/run_3/run_3_1_lowPt25to35/2018_SingleMu/IMPORTANT/L1TJetEtSF_2018_SingleMu_L1JetPt25To35_20211202.root"
#sHistoCalibSF  = "h_jet_byHand_res_vs_iEta_vs_nVtx_$JETSHAPE$PUSALGORITHM_HBEF_medPt_0_emu_$OOTPUS_PU50to100_SF_PFJByL1TJ"; # medPt: SF derived for PFJetpT in [60, 90] GeV
sHistoCalibSF  = "h_jet_byHand_res_vs_iEta_vs_nVtx$JETSHAPE_$PUSALGORITHM_HBEF_%s_0_emu_$OOTPUS_%s_SF_PFJByL1TJ" % (PtRangeForCalibSF, PURangeUsedForCalibSF); # medPt: SF derived for PFJetpT in [60, 90] GeV

PT_MIN = 10 # 30   ## Minimum offline jet pT to consider
DR_MIN = 1.2  ## Minimum dR between offline jets considered for measurement
DR_MAX = 0.3  ## Maximum dR between L1T and offline jets for matching


PU_CUT = 'nVtxMin_45'  ## Pileup selection
ERA    = '2018AB'      ## Era, e.g. 2018AB or 2018D

PT_CAT = {}
#PT_CAT['lowPt'] = [30,  60,   60]  ## Low pT, turn-on threshold, high pT
#PT_CAT['medPt'] = [60,  75,   90]  ## Low pT, turn-on threshold, high pT   #[60,  90,   90]
#PT_CAT['hiPt']  = [90, 120, 9999]  ## Low pT, turn-on threshold, high pT
#PT_CAT['modPt'] = [60,  75,   90]  ## Low pT, turn-on threshold, high pT
## for lowpt jet between 25 to 35 GeV
#PT_CAT['Below25Pt'] = [0,  15,   25]  ## Low pT, turn-on threshold, high pT
#PT_CAT['lowPt'] = [25,  30,   35]  ## Low pT, turn-on threshold, high pT
#PT_CAT['medPt'] = [35,  75,   90]  ## Low pT, turn-on threshold, high pT   #[60,  90,   90]
#PT_CAT['hiPt']  = [90, 120, 9999]  ## Low pT, turn-on threshold, high pT
## for lowpt jet between 25 to 35 GeV - version 2
PT_CAT['Ptlt25']   = [ 0,  15,   25]  ## Low pT, turn-on threshold, high pT
PT_CAT['Pt25To35'] = [25,  30,   35]  ## Low pT, turn-on threshold, high pT
PT_CAT['Pt35To60'] = [35,  55,   60]  ## Low pT, turn-on threshold, high pT
PT_CAT['Pt60To90'] = [60,  75,   90]  ## Low pT, turn-on threshold, high pT   #[60,  90,   90]
PT_CAT['Ptgt90']   = [90, 120, 9999]  ## Low pT, turn-on threshold, high pT
PtTrshForTau = {'Pt35To60': 40}

if runMode in ['CalCalibSF', 'makeInputForML']:
    PT_CAT['Below30Pt'] = [0,  30,   30]  ## Low pT, turn-on threshold, high pT

ETA_CAT = {}
ETA_CAT['HBEF'] = [0.000, 5.210]  ## Whole detector, 1 - 41
ETA_CAT['HB']   = [0.000, 1.392]  ## Trigger towers  1 - 16
ETA_CAT['HE1']  = [1.392, 1.740]  ## Trigger towers 17 - 20
ETA_CAT['HE2a'] = [1.740, 2.322]  ## Trigger towers 21 - 25
ETA_CAT['HE2b'] = [2.322, 3.000]  ## Trigger towers 26 - 28
ETA_CAT['HF']   = [3.000, 5.210]  ## Trigger towers 30 - 41

IETA_CAT = {}
IETA_CAT['HBEF'] = [ 1, 41]  ## Whole detector, 1 - 41
IETA_CAT['HB']   = [ 1, 16]  ## Trigger towers  1 - 16
IETA_CAT['HE1']  = [17, 20]  ## Trigger towers 17 - 20
IETA_CAT['HE2a'] = [21, 25]  ## Trigger towers 21 - 25
IETA_CAT['HE2b'] = [26, 28]  ## Trigger towers 26 - 28
IETA_CAT['HF']   = [30, 41]  ## Trigger towers 30 - 41

HCAL_ETA_CAT = ['HB', 'HE1', 'HE2', 'HF']

useAbsEtaBins = True
ETA_Bins = []
for iEta in range(-41,42):
    if iEta in [-29, 0, 29]:        continue;
    if useAbsEtaBins and iEta < 0:  continue;
    ETA_Bins.append(str(iEta))
#ETA_Bins.append('all')
ETA_Bins.append('HBEF')

map_iEta_Eta = {
    1:  [0.000, 0.087],
    2:  [0.087, 0.174],
    3:  [0.174, 0.261],
    4:  [0.261, 0.348],
    5:  [0.348, 0.435],
    6:  [0.435, 0.522],
    7:  [0.522, 0.609],
    8:  [0.609, 0.695],
    9:  [0.695, 0.783],
    10: [0.783, 0.870],

    11: [0.870, 0.957],
    12: [0.957, 1.044],
    13: [1.044, 1.131],
    14: [1.131, 1.218],
    15: [1.218, 1.305],
    16: [1.305, 1.392],
    17: [1.392, 1.479],
    18: [1.479, 1.566], # endcap starts
    19: [1.566, 1.653],    
    20: [1.653, 1.740],
    
    21: [1.740, 1.830],
    22: [1.830, 1.930],
    23: [1.930, 2.043],
    24: [2.043, 2.172],
    25: [2.172, 2.322],
    26: [2.322, 2.500],
    27: [2.500, 2.650],
    28: [2.650, 3.000], # endcap ends   
}

map_Jet_TowerIEta_IEta = {
    #0:  X,
    1:  1,
    2:  3,
    3:  5,
    4:  7,
    5:  9,
    6: 11,
    7: 13,
    8: 15,
    9: 17,
    10: 19,
    11: 21,
    12: 23,
    13: 25,
    14: 27,
    15: 29,
    16: 31,
    17: 33,
    18: 35,
    19: 37,
    20: 39,
    21: 41,
    22: 43,
    23: 46,
    24: 48,
    25: 52,
    26: 55,
    27: 59,
    28: 63,
    #29: XX,
    30: 69,
    31: 74,
    32: 78,
    33: 82,
    34: 86,
    35: 90,
    36: 94,
    37: 98,
    38: 102,
    39: 106,
    40: 110,
    41: 116
}

def convert_jetIPhi_to_jetTowerIPhi(jetIPhi):
    jetTowerIPhi = (jetIPhi + 1) / 2
    return jetTowerIPhi


jetPtBins_forCalibration = [3, 6, 9, 12, 15, 20, 25, 30, 35, 40, 45, 55, 70, 90, 120, 160, 200, 999.0]



## Delta-phi calculation for integer phi values (1 - 72) with wrap-around
def dIPhi(iPhiA, iPhiB):
    if iPhiA > iPhiB:
        if abs((iPhiA - 72) - iPhiB) < abs(iPhiA - iPhiB):
            return (iPhiA - 72) - iPhiB
        else:
            return iPhiA - iPhiB
    else:
        if abs(iPhiA - (iPhiB - 72)) < abs(iPhiA - iPhiB):
            return iPhiA - (iPhiB - 72)
        else:
            return iPhiA - iPhiB

## Delta-eta calculation taking into account there is no iEta = 0 or +/-29
def dIEta(iEtaA, iEtaB):
    dIEtaAbs = abs(iEtaA - iEtaB)
    if (iEtaA < -29) != (iEtaB < -29):
        dIEtaAbs -= 1
    if (iEtaA <   0) != (iEtaB <   0):
        dIEtaAbs -= 1
    if (iEtaA <  29) != (iEtaB <  29):
        dIEtaAbs -= 1
    if (iEtaA < iEtaB):
        dIEtaAbs *= -1
    return dIEtaAbs


def calculateJetIEta(eta):
    jetIEta_offlineJet = -50.0 # None # abs(vOff.Eta())
    for iEta, etaBinRange in map_iEta_Eta.items():
        if abs(eta) >= etaBinRange[0] and abs(eta) < etaBinRange[1]:
            jetIEta_offlineJet = float( iEta * math.copysign(1, eta) )
    
    return jetIEta_offlineJet

def convert_jetIEta_to_jetTowerIEta(jetIEta):
    jetTowerIEta = None
    jetIEtaAbs = abs(jetIEta)
    for jetTowerIEta_1, jetIEta_1 in map_Jet_TowerIEta_IEta.items():
        if jetIEtaAbs == jetIEta_1:
            jetTowerIEta = int( jetTowerIEta_1 * math.copysign(1, jetIEta) )
            break
    
    if not jetTowerIEta:
        print "convert_jetIEta_to_jetTowerIEta(jetIEta):: jetTowerIEta not set... why ?? \t\t **** ERROR ****"
        exit(0)
    
    return jetTowerIEta

def getJetPtCategory(pt):
    iPFJetPtCat = 'None'
    for iCat in PT_CAT.keys():
        if pt >= PT_CAT[iCat][0] and pt < PT_CAT[iCat][2]:
            iPFJetPtCat = iCat
            
    #if iPFJetPtCat == 'None':
    #    print "getJetPtCategory():: pt {}, iPFJetPtCat {}, PT_CAT {}".format(pt, iPFJetPtCat, PT_CAT)
    return iPFJetPtCat

                        
def run():
    #
    # argument 1: OOT PU subtraction scheme name. For e.g. def, PFA2, PFA1p
    # argument 2: i/p L1 ntuple directory
    # argument 3: N.  Split total number of events into Nth part. For e.g. 10
    # argument 4: m. Run on mth quantile events out of spillted N events. For e.g. 0, 1, 2, 3 ..., 9

    print '\nInside L1T_JetMET_Res\n'

    in_dir = '/afs/cern.ch/work/a/abrinke1/public/L1T/CaloTP/2018/HCAL_OOT_PUS/L1TNtuples/output/'

    in_file_names = []
    # in_file_names.append(in_dir+'L1Ntuple_HCAL_TP_OOT_PUS_'+PU_CUT+'_nSamp_2_nPresamp_0_HB_1p0_1p0_HE1_1p0_1p0_HE2_1p0_1p0_'+ERA+'_100k.root')
    # in_file_names.append(in_dir+'L1Ntuple_HCAL_TP_OOT_PUS_'+PU_CUT+'_nSamp_1_nPresamp_0_HB_1p0_HE1_1p0_HE2_1p0_'+ERA+'_100k.root')
    # in_file_names.append(in_dir+'L1Ntuple_HCAL_TP_OOT_PUS_'+PU_CUT+'_nSamp_2_nPresamp_1_HB_n0p51_1p0_HE1_n0p49_1p0_HE2_n0p45_1p0_'+ERA+'_100k.root')
    # in_file_names.append(in_dir+'L1Ntuple_HCAL_TP_OOT_PUS_'+PU_CUT+'_nSamp_3_nPresamp_1_HB_n0p71_1p0_1p0_HE1_n0p71_1p0_1p0_HE2_n1p17_1p0_1p0_'+ERA+'_100k.root')
    # #in_file_names.append(in_dir+'L1Ntuple_HCAL_TP_OOT_PUS_'+PU_CUT+'_nSamp_2_nPresamp_1_HB_n1p0_1p0_HE1_n1p0_1p0_HE2_n1p0_1p0_'+ERA+'_100k.root')
    # #in_file_names.append(in_dir+'L1Ntuple_HCAL_TP_OOT_PUS_'+PU_CUT+'_nSamp_3_nPresamp_2_HB_n0p5_n0p5_1p0_HE1_n0p5_n0p5_1p0_HE2_n0p5_n0p5_1p0_'+ERA+'_100k.root')
    # #in_file_names.append(in_dir+'L1Ntuple_HCAL_TP_OOT_PUS_'+PU_CUT+'_nSamp_3_nPresamp_1_HB_n1p0_1p0_1p0_HE1_n1p0_1p0_1p0_HE2_n1p0_1p0_1p0_'+ERA+'_100k.root')
    # #in_file_names.append(in_dir+'L1Ntuple_HCAL_TP_OOT_PUS_'+PU_CUT+'_nSamp_3_nPresamp_1_HB_n2p0_1p0_1p0_HE1_n2p0_1p0_1p0_HE2_n2p0_1p0_1p0_'+ERA+'_100k.root')
    # #in_file_names.append(in_dir+'L1Ntuple_HCAL_TP_OOT_PUS_'+PU_CUT+'_nSamp_4_nPresamp_2_HB_n0p5_n0p5_1p0_1p0_HE1_n0p5_n0p5_1p0_1p0_HE2_n0p5_n0p5_1p0_1p0_'+ERA+'_100k.root')
    # #in_file_names.append(in_dir+'L1Ntuple_HCAL_TP_OOT_PUS_'+PU_CUT+'_nSamp_4_nPresamp_2_HB_n1p0_n1p0_1p0_1p0_HE1_n1p0_n1p0_1p0_1p0_HE2_n1p0_n1p0_1p0_1p0_'+ERA+'_100k.root')

    # in_file_names.append(in_dir+'L1Ntuple_HCAL_TP_OOT_PUS_Run2018D_ExpressPhysics_FEVT_PFA2.root')
    # in_file_names.append(in_dir+'L1Ntuple_HCAL_TP_OOT_PUS_Run2018D_ExpressPhysics_FEVT_PFA1p.root')

    # in_file_names.append(in_dir+'L1Ntuple_HCAL_TP_OOT_PUS_nVtxMin_45_FEVT_nSamp_2_nPresamp_1_HB_n0p51_1p0_HE1_n0p49_1p0_HE2_n0p45_1p0_2018AB_100k.root')
    # in_file_names.append(in_dir+'L1Ntuple_HCAL_TP_OOT_PUS_nVtxMin_45_FEVT2_nSamp_2_nPresamp_1_HB_n0p51_1p0_HE1_n0p49_1p0_HE2_n0p45_1p0_2018AB_100k.root')

    #in_file_names.append(in_dir+'L1Ntuple_HCAL_TP_OOT_PUS_nVtxMin_45_nSamp_2_nPresamp_0_HB_1p0_1p0_HE1_1p0_1p0_HE2_1p0_1p0_2018AB_-1k.root')
    
    in_file_names.append(in_dir+'L1Ntuple_HCAL_TP_OOT_PUS_nSamp_2_nPresamp_0_HB_1p0_1p0_HE1_1p0_1p0_HE2_1p0_1p0_2018D_5k.root')
    

    
    if not os.path.exists('plots'): os.makedirs('plots')

    out_file_str = 'L1T_JetMET_Res_'+PU_CUT+'_'+ERA+'_Comp%d' % len(in_file_names)
    out_file_str += ('_%dk' % (MAX_EVT / 1000))

    # out_file_str = 'L1T_JetMET_Res_nVtxMin_45_FEVT_2018AB_Comp2_100k'
    out_file_str = 'L1T_JetMET_Res_nVtxMin_45_nSamp_2_nPresamp_0_HB_1p0_1p0_HE1_1p0_1p0_HE2_1p0_1p0_2018AB_-1k'
    out_file_str = 'plots/'+out_file_str+'.root'

    usePUReweighting = False
    isMC = False
    
    N_parts = None
    M_quantilesIpFilesSet = None
    sipFileCalibSF_toUse = ""
    OOT_PU_scheme = None
    print "python sys.argv: {}".format(sys.argv)
    if len(sys.argv) >= 3:        
        OOT_PU_scheme = sys.argv[1]
        in_file_Dirs  = sys.argv[2]
        PUrangeTag    = sys.argv[3]
        N_parts         = None
        M_quantilesIpFilesSet = None
        if len(sys.argv) >= 5:
            N_parts         = int(sys.argv[4])
            M_quantilesIpFilesSet = int(sys.argv[5])
        
        print "in_file_Dirs: {}".format(in_file_Dirs)
        in_file_Dirs  = in_file_Dirs.replace("'","")
        in_file_Dirs  = in_file_Dirs.replace('"','')
        print "in_file_Dirs: {}".format(in_file_Dirs)
        if in_file_Dirs.endswith('/'): in_file_Dirs = in_file_Dirs[:-1]
        print "in_file_Dirs: {}".format(in_file_Dirs)
        in_file_names = [ in_file_Dirs + "/L1Ntuple_*.root" ]
        print "in_filenames_0: {}".format(in_file_names)
        if runMode in ['test'] or (MAX_EVT > 0 and MAX_EVT < 1000): in_file_names = [ in_file_Dirs + "/L1Ntuple_*_1.root" ]
        print "in_filenames: {}".format(in_file_names)
        
        out_file_str  = "L1T_JetMET_Res_%s_%s.root" % (OOT_PU_scheme, PUrangeTag)
        if runMode in ['CalCalibSF']:
            out_file_str = out_file_str.replace(".root", "_CalCalibSF.root")
        if runMode in ['CalibJetByHand']:
            out_file_str = out_file_str.replace(".root", "_wCalibJetByHand.root")
        if sOutExt:
            out_file_str = out_file_str.replace(".root", "%s.root" % (sOutExt))
        if N_parts:
            #out_file_str  = "L1T_JetMET_Res_%s_part%d_of_%d.root" % (OOT_PU_scheme, M_quantilesIpFilesSet,N_parts)
            out_file_str = out_file_str.replace(".root", "_part%d_of_%d.root" % (M_quantilesIpFilesSet,N_parts))
        print "out_file_str: {}".format(out_file_str)
        
        if usePUReweighting_0:
            isMC = False
            for in_file in in_file_names:
                if "mc" in in_file.lower():
                    isMC = True
            '''
            if not isMC:
                print "input ntuples are not for MC.... and still running with usePUReweighting == true ??? \t\t **** ERROR ****"
                exit(0)
            '''
            if isMC: usePUReweighting = True
        
        print "OOT_PU_scheme: {}".format(OOT_PU_scheme)
        '''
        if OOT_PU_scheme.lower() in ['def', 'pfa2']:            
            sipFileCalibSF_toUse = sipFileCalibSF.replace('$OOTPUS', 'PFA2')
        elif OOT_PU_scheme.lower() in ['pfa1p']:
            sipFileCalibSF_toUse = sipFileCalibSF.replace('$OOTPUS', 'PFA1p')
        '''
        sipFileCalibSF_toUse = sipFileCalibSF
    
    
    fOut_MLInputs = None
    fOut_MLInputs_writer = None
    isFirstEntry_WriteInputForML = True
    if runMode in ['makeInputForML']:        
        out_file_str_0 = out_file_str.replace('.root','.csv')
        print "runMode=makeInputForML:: output file: {}".format(out_file_str_0)
        fOut_MLInputs = open(out_file_str_0, mode='w')

        '''
        fieldnames = ['emp_name', 'dept', 'birth_month']
        writer = csv.DictWriter(fOut_MLInputs, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerow({'emp_name': 'John Smith', 'dept': 'Accounting', 'birth_month': 'November'})
        writer.writerow({'emp_name': 'Erica Meyers', 'dept': 'IT', 'birth_month': 'March'})

        fOut_MLInputs.close()
        '''

    if runMode in ['makePUHisto']:
        print "Running with runMode = makePUHisto \n"
        
                
    chains = {}
    chains['Evt'] = []  ## Event info
    chains['Vtx'] = []  ## RECO vertex info
    chains['Jet'] = []  ## RECO jet info
    chains['Unp'] = []  ## Unpacked Layer-2 info
    chains['Emu'] = []  ## Emulated Layer-2 info
    chains['uTP'] = []  ## Unpacked trigger primitives
    chains['eTP'] = []  ## Emulated trigger primitives    
    for i in range(len(in_file_names)):
        print 'Adding file %s' % in_file_names[i]
        in_file_names_ith = glob.glob(in_file_names[i])
        print "in_file_names_ith ({}): {}".format(len(in_file_names_ith), in_file_names_ith)
        
        chains['Evt'].append( R.TChain('l1EventTree/L1EventTree') )
        chains['Vtx'].append( R.TChain('l1RecoTree/RecoTree') )
        chains['Jet'].append( R.TChain('l1JetRecoTree/JetRecoTree') )
        chains['Unp'].append( R.TChain('l1UpgradeTree/L1UpgradeTree') )
        chains['Emu'].append( R.TChain('l1UpgradeEmuTree/L1UpgradeTree') )
        chains['uTP'].append( R.TChain('l1CaloTowerTree/L1CaloTowerTree') )
        chains['eTP'].append( R.TChain('l1CaloTowerEmuTree/L1CaloTowerTree') )
        '''
        chains['Evt'][i].Add( in_file_names[i] )
        chains['Vtx'][i].Add( in_file_names[i] )
        chains['Jet'][i].Add( in_file_names[i] )
        chains['Unp'][i].Add( in_file_names[i] )
        chains['Emu'][i].Add( in_file_names[i] )
        chains['uTP'][i].Add( in_file_names[i] )
        chains['eTP'][i].Add( in_file_names[i] )
        '''
        Ntot = len(in_file_names_ith)
        print '\n\nFor chain %d, total number of input ntuples: %d' % (i, Ntot); sys.stdout.flush();
        Nfirst = 0;  Nlast = Ntot # initial and last input file number
        if N_parts: # split total number of input files into parts
            Nfirst = int(1.0 * Ntot / N_parts * (M_quantilesIpFilesSet)) 
            Nlast  = int(1.0 * Ntot / N_parts * (M_quantilesIpFilesSet + 1) - 1)
        print "M_quantilesIpFilesSet: {}, N_parts: {}, Nfirst: {}, Nlast: {}".format(M_quantilesIpFilesSet, N_parts, Nfirst, Nlast)
        
        for ith_file in range(Nfirst, Nlast+1):
            in_file_name = in_file_names_ith[ith_file]
            print "   %s" %  in_file_name
            chains['Evt'][i].Add( in_file_name )
            chains['Vtx'][i].Add( in_file_name )
            chains['Jet'][i].Add( in_file_name )
            chains['Unp'][i].Add( in_file_name )
            chains['Emu'][i].Add( in_file_name )
            chains['uTP'][i].Add( in_file_name )
            chains['eTP'][i].Add( in_file_name )
        print " %d input files read. " % ( (Nlast+1) - Nfirst )
            


    ###################
    ### Book histograms
    ###################
    if PrintLevel >= 1:
        print "Booking hitograms"
    hDummy = R.TH1D("hDummy","",1,0,1)
    hDummy.SetDefaultSumw2()

    pt_bins  = [20,    0, 200]
    res_bins = [80, -1.5, 2.5]
    dR_bins  = [35,  0.0, 0.35]
    puEt_bins = [60,   0, 120]
    PUByRawPt_bins = [80, 0, 2]
    iEta_bins = [83, -41.5, 41.5]
    if useAbsEtaBins:
        iEta_bins = [41, 0.5, 41.5]

    # use nTotalEvents to normalize to min-bias events
    hnTotalEvents = R.TH1D("h_nTotalEvents", ";Channels;Events", 11,-0.5,10.5)

    hStat = R.TH1D("h_Stat", ";Conditions;Events", 51,-0.5,50.5)
    
    ## Variable distributions
    dists = ['jet_eff', 'jet_num', 'jet_den', 'jet_res', 'jet_dR']

    hist = {}
    ## Loop over all distributions
    for dist in dists:
        hist[dist] = {}
        ## Loop over L1T jet algorithms
        for algo in ['PUS','noPUS','Raw','RawPUS']:
            hist[dist][algo] = {}
            ## Loop over eta regions
            for iEta in ETA_CAT.keys():
                hist[dist][algo][iEta] = {}
                ## Loop over pT ranges
                for iPt in PT_CAT.keys()+['PtAllBins']:
                    hist[dist][algo][iEta][iPt] = {}
                    ## Loop over unpacked and emulated
                    for src in ['unp','emu']:
                        hist[dist][algo][iEta][iPt][src] = []
                        ## Separate histogram for each input file (different TP reconstructions)
                        for iTP in range(len(in_file_names)):

                            ## Pick binning for each histogram
                            h_bins = []
                            if   dist == 'jet_eff': h_bins = []
                            elif dist == 'jet_num': h_bins = pt_bins
                            elif dist == 'jet_den': h_bins = pt_bins
                            elif dist == 'jet_res': h_bins = res_bins
                            elif dist == 'jet_dR':  h_bins = dR_bins
                            else: print '\nInvalid distribution %s - no binning found!!! \nQuitting.', sys.exit()

                            if dist == 'jet_eff':
                                hist[dist][algo][iEta][iPt][src].append( iTP+1 )
                            else:
                                hist[dist][algo][iEta][iPt][src].append( R.TH1D( 'h_%s_%s_%s_%s_%d_%s' % (dist, algo, iEta, iPt, iTP, src),
                                                                                 'L1T %s %s %s in %s for %s' % (src, algo, dist, iEta, iPt),
                                                                                 h_bins[0], h_bins[1], h_bins[2] ) )

                        ## End loop: for iTP in range(len(in_file_names))
                    ## End loop: for src in ['unp','emu']
                ## End loop: for iPt in PT_CAT.keys()
            ## End loop: for iEta in ETA_CAT.keys()
        ## End loop: for algo in ['PUS','noPUS','Raw','RawPUS']
    ## End loop: for dist in dists


    dists1 = ['l1jetEt_vs_RawEtMinusPU']
    
    hist1 = {}
    for dist in dists1:
        hist1[dist] = {}
        for iEta in ETA_Bins:            
            hist1[dist][iEta] = {}
            for src in ['unp', 'emu']:
                hist1[dist][iEta][src] = []
                
                ## Separate histogram for each input file (different TP reconstructions)
                for iTP in range(len(in_file_names)):
                    hist1[dist][iEta][src] .append( R.TH2D( 'h_%s_%s_%d_%s' % (dist, iEta, iTP, src),
                                                                                 'L1T %s %s in ieta %s' % (src, dist, iEta),
                                                                                 100,0,200, 100,0,3) )
    
    
    dists2 = [
        'jet_byHand_eff', 'jet_byHand_num', 'jet_byHand_den', #'jet_byHand_res_vs_iEta', 'jet_byHand_PU_vs_iEta', 'jet_byHand_PUByRawPt_vs_iEta',
        # for calibration with PtAllBins 
        #'jet_byHand_L1JetPt_vs_PFJetPt',
        #'jet_byHand_L1JetPt_vs_DefaultL1JetPt',
        'jet_byHand_res_vs_iEta_vs_nVtx', 
    ]
    hist2 = {}
    #print "\nSetting hist2::\ndists2: {}".format(dists2)
    #for jetShape in ['Default'] + JetShapes:
    for jetShape in JetShapes:
        # JetShape = "" plots are with the first version of code for 9x9 jets
        jetShape1 = jetShape
        if jetShape == 'Default':  jetShape1 = ""
        else:                      jetShape1 = "_%s" % (jetShape)
        
        ## Loop over all distributions
        for dist_1 in dists2:
            dist = '%s%s' % (dist_1, jetShape1)
            #print "    dist: {}".format(dist)
            
            hist2[dist] = {}
            ## Loop over L1T jet algorithms
            for algo in PUSAlgosAll: # ['Raw', 'RawPUS', 'RawPUS_phiRingMin4', 'RawPUS_phiRingSide4', 'RawPUS_phiRingAdjacent']: # ['PUS','noPUS','Raw','RawPUS']:
                if (algo == 'L1JDefault') and (jetShape != 'Default'): continue
                
                if dist_1 in ['jet_byHand_L1JetPt_vs_DefaultL1JetPt'] and algo != 'RawPUS': continue
                hist2[dist][algo] = {}                
                
                #print "   algo : {}".format(algo)
                ## Loop over eta regions
                for iEta in ETA_Bins: 
                    hist2[dist][algo][iEta] = {}
                    #print "    iEta: {}".format(iEta)
                    ## Loop over pT ranges
                    for iPt in PT_CAT.keys() + ['PtAllBins']:
                        # L1JetPt_vs_PFJetPt plot to be plot with PtAllBins
                        if dist_1 in ['jet_byHand_L1JetPt_vs_PFJetPt', 'jet_byHand_L1JetPt_vs_DefaultL1JetPt'] and iPt != 'PtAllBins': continue
                        
                        hist2[dist][algo][iEta][iPt] = {}
                        #print "    iPt: {}".format(iPt)
                        ## Loop over unpacked and emulated
                        for src in ['unp','emu']: #['unp','emu']:
                            hist2[dist][algo][iEta][iPt][src] = []
                            #print "    src: {}".format(src)
                            ## Separate histogram for each input file (different TP reconstructions)
                            for iTP in range(len(in_file_names)):
                                ## Pick binning for each histogram
                                h_bins = []
                                if   'eff' in dist: h_bins = []
                                elif 'num' in dist: h_bins = pt_bins
                                elif 'den' in dist: h_bins = pt_bins
                                elif 'res' in dist: h_bins = res_bins
                                elif 'dR'  in dist: h_bins = dR_bins
                                elif '_PU_'  in dist: h_bins = puEt_bins
                                elif '_PUByRawPt_'  in dist: h_bins = PUByRawPt_bins
                                elif 'jet_byHand_L1JetPt_vs_PFJetPt' in dist: h_bins = pt_bins
                                elif 'jet_byHand_L1JetPt_vs_DefaultL1JetPt' in dist: h_bins = res_bins
                                else: print '\nInvalid distribution %s - no binning found!!! \nQuitting.', sys.exit()
                                
                                #print "      dist {}, algo {}, iEta {}, iPt {}, src {}, iTP {}".format(dist, algo, iEta, iPt, src, iTP)
                                
                                if   'eff' in dist:
                                    hist2[dist][algo][iEta][iPt][src].append( iTP+1 )
                                #elif dist in ['jet_byHand_num', 'jet_byHand_den']:
                                elif 'jet_byHand_num' in dist or 'jet_byHand_den' in dist:
                                    # TH1D
                                    hist2[dist][algo][iEta][iPt][src].append( R.TH1D( 'h_%s_%s_%s_%s_%d_%s' % (dist, algo, iEta, iPt, iTP, src),
                                                                                     'L1T %s %s %s in %s for %s' % (src, algo, dist, iEta, iPt),
                                                                                     h_bins[0], h_bins[1], h_bins[2] ) )
                                elif dist_1 in ['jet_byHand_L1JetPt_vs_PFJetPt']:
                                    # TH1D,  for calibration with PtAllBins 
                                    hist2[dist][algo][iEta][iPt][src].append( R.TH2D( 'h_%s_%s_%s_%s_%d_%s_fixBinWidthPFJetPt' % (dist, algo, iEta, iPt, iTP, src),
                                                                                      'L1T %s %s %s in %s for %s' % (src, algo, dist, iEta, iPt),
                                                                                      len(jetPtBins_forCalibration)-1, array('d',jetPtBins_forCalibration),
                                                                                      h_bins[0], h_bins[1], h_bins[2] ) )
                                elif dist_1 in ['jet_byHand_L1JetPt_vs_DefaultL1JetPt'] and jetShape == 'Default':
                                    hist2[dist][algo][iEta][iPt][src].append( R.TH2D( 'h_%s_%s_%s_%s_%d_%s' % (dist, algo, iEta, iPt, iTP, src),
                                                                                      'L1T %s %s %s in %s for %s' % (src, algo, dist, iEta, iPt),
                                                                                      len(jetPtBins_forCalibration)-1, array('d',jetPtBins_forCalibration),
                                                                                      h_bins[0], h_bins[1], h_bins[2]) )
                                elif dist_1 in ['jet_byHand_res_vs_iEta', 'jet_byHand_PU_vs_iEta', 'jet_byHand_PUByRawPt_vs_iEta']:
                                    # TH2D
                                    if iEta == 'HBEF':
                                        hist2[dist][algo][iEta][iPt][src].append( R.TH2D( 'h_%s_%s_%s_%s_%d_%s' % (dist, algo, iEta, iPt, iTP, src),
                                                                                          'L1T %s %s %s in %s for %s' % (src, algo, dist, iEta, iPt),
                                                                                          iEta_bins[0],iEta_bins[1],iEta_bins[2], h_bins[0], h_bins[1], h_bins[2] ) )
                                    else:
                                        continue
                                elif dist_1 in ['jet_byHand_res_vs_iEta_vs_nVtx']:
                                    if iEta == 'HBEF':
                                        hist2[dist][algo][iEta][iPt][src].append( R.TH3D( 'h_%s_%s_%s_%s_%d_%s' % (dist, algo, iEta, iPt, iTP, src),
                                                                                          'L1T %s %s %s in %s for %s' % (src, algo, dist, iEta, iPt),
                                                                                          iEta_bins[0],iEta_bins[1],iEta_bins[2], 101,-0.5,100.5, h_bins[0], h_bins[1], h_bins[2] ) )

                                #print "      dist {}, algo {}, iEta {}, iPt {}, src {}, iTP {}".format(dist, algo, iEta, iPt, src, iTP)    
    #print "\nhist2.keys(): {}".format(hist2.keys())


    hist_PFJetPt_iEtawise = {}
    for iEta in ETA_Bins:
        hist_PFJetPt_iEtawise[iEta] = R.TH1D('h_PFJetPt_%s' % (iEta), 'h_PFJetPt_%s' % (iEta), pt_bins[0], pt_bins[1], pt_bins[2])


    #hist_L1Jet_unp_TowerIEta_vs_IEta = R.TH2D('h_L1Jet_unp_TowerIEta_vs_IEta', 'h_L1Jet_unp_TowerIEta_vs_IEta', 101,-50.5,50.5, 241,-120.5,120.5)
    #hist_L1Jet_emu_TowerIEta_vs_IEta = R.TH2D('h_L1Jet_emu_TowerIEta_vs_IEta', 'h_L1Jet_emu_TowerIEta_vs_IEta', 101,-50.5,50.5, 241,-120.5,120.5)

    #hist_L1Jet_unp_TowerIPhi_vs_IPhi = R.TH2D('h_L1Jet_unp_TowerIPhi_vs_IPhi', 'h_L1Jet_unp_TowerIPhi_vs_IPhi', 161,-0.5,160.5, 161,-0.5,160.5)
    #hist_L1Jet_emu_TowerIPhi_vs_IPhi = R.TH2D('h_L1Jet_emu_TowerIPhi_vs_IPhi', 'h_L1Jet_emu_TowerIPhi_vs_IPhi', 161,-0.5,160.5, 161,-0.5,160.5)
    
    '''
    hist_nPV_vs_L1JetDefaultRAW_SF = {}
    hist_nPV_vs_L1JetDefaultPUS_SF = {}
    ## Loop over unpacked and emulated
    for src in ['unp','emu']:
        hist_nPV_vs_L1JetDefaultRAW_SF[src] = {}
        hist_nPV_vs_L1JetDefaultPUS_SF[src] = {}
        for iEta in ETA_Bins:
            hist_nPV_vs_L1JetDefaultRAW_SF[src][iEta] = R.TH2D('h_nPV_vs_L1JetDefaultRAW_SF_%s_%s' % (src,iEta), 'h_nPV_vs_L1JetDefaultRAW_SF_%s_%s' % (src,iEta),
                                                               101,-0.5,100.5, 150,0,3.8)
            hist_nPV_vs_L1JetDefaultPUS_SF[src][iEta] = R.TH2D('h_nPV_vs_L1JetDefaultPUS_SF_%s_%s' % (src,iEta), 'h_nPV_vs_L1JetDefaultPUS_SF_%s_%s' % (src,iEta),
                                                               101,-0.5,100.5, 150,0,3.8)
    '''
    

    jetRate_bins = [400, 0., 400.0];   jetRate_binWidth = (jetRate_bins[2] - jetRate_bins[1]) / jetRate_bins[0]; 
    eff_bins = [500, 0, 500]
    PU_bins  = [101, -0.5, 100.5]
    sAxexName_JetRate = ";Threshold E_{T} (GeV);nPV;Rate (Hz)"
    sAxexName_Eff     = ";E_{T} (GeV);nPV;Events / bin"

    dists3 = [
        'jet_byHand_eff_num_vs_PU', 'jet_byHand_eff_den_vs_PU'
    ]
    dists4 = [
        'jet_byHand_rates_singleJet', 'jet_byHand_rates_doubleJet', 'jet_byHand_rates_trippleJet', 'jet_byHand_rates_quadJet'
    ]
    
    hist3 = {}
    #print "\nSetting hist3::\ndists3: {}".format(dists3)
    #for jetShape in ['Default'] + JetShapes:
    for jetShape in JetShapes + JetShapesType2:
        # JetShape = "" plots are with the first version of code for 9x9 jets
        jetShape1 = jetShape
        if jetShape == 'Default':  jetShape1 = ""
        else:                      jetShape1 = "_%s" % (jetShape)
        
        ## Loop over all distributions
        for dist_1 in dists3:
            dist = '%s%s' % (dist_1, jetShape1)
            #print "    dist: {}".format(dist)
            
            hist3[dist] = {}
            ## Loop over L1T jet algorithms
            for algo in PUSAlgosAll + PUSAlgosAllType2: # ['Raw', 'RawPUS', 'RawPUS_phiRingMin4', 'RawPUS_phiRingSide4', 'RawPUS_phiRingAdjacent']: # ['PUS','noPUS','Raw','RawPUS']:
                # read proper jetShape and PUSAlgo conbination
                if (jetShape in JetShapes      and algo not in PUSAlgosAll) or \
                   (jetShape in JetShapesType2 and algo not in PUSAlgosAllType2 ):
                    continue
                
                if (algo == 'L1JDefault') and (jetShape != 'Default'): continue
                
                hist3[dist][algo] = {}                
                
                #print "   algo : {}".format(algo)
                ## Loop over eta categories
                for iEta in ETA_CAT.keys():
                    hist3[dist][algo][iEta] = {}
                    #print "    iEta: {}".format(iEta)
                    
                    ## Loop over pT ranges
                    for iPt in PT_CAT.keys():                        
                        hist3[dist][algo][iEta][iPt] = {}
                        #print "    iPt: {}".format(iPt)
                        
                        ## Loop over unpacked and emulated
                        for src in ['unp','emu']: #['unp','emu']:
                            hist3[dist][algo][iEta][iPt][src] = []
                            #print "    src: {}".format(src)
                            
                            ## Separate histogram for each input file (different TP reconstructions)
                            for iTP in range(len(in_file_names)):
                                ## Pick binning for each histogram
                                h_bins = []
                                sAxexName = None
                                if   'eff_num' in dist or 'eff_den' in dist:
                                    h_bins = eff_bins
                                    sAxexName = sAxexName_Eff
                                else: print '\nInvalid distribution %s - no binning found!!! \nQuitting.', sys.exit()
                                
                                #print "      dist {}, algo {}, iEta {}, iPt {}, src {}, iTP {}".format(dist, algo, iEta, iPt, src, iTP)
                                
                                if 'jet_byHand_eff_num_vs_PU' in dist or 'jet_byHand_eff_den_vs_PU' in dist:
                                    # TH2D
                                    hist3[dist][algo][iEta][iPt][src].append( R.TH2D( 'h_%s_%s_%s_TrgTrsh%d_%d_%s' % (dist, algo, iEta, PT_CAT[iPt][1], iTP, src),
                                                                                      sAxexName,
                                                                                      int(h_bins[0]), h_bins[1], h_bins[2],
                                                                                      int(PU_bins[0]),PU_bins[1],PU_bins[2] ) )
    

    hist4 = {}
    #print "\nSetting hist4::\ndists4: {}".format(dists4)
    #for jetShape in ['Default'] + JetShapes:
    for jetShape in JetShapes + JetShapesType2:
        # JetShape = "" plots are with the first version of code for 9x9 jets
        jetShape1 = jetShape
        if jetShape == 'Default':  jetShape1 = ""
        else:                      jetShape1 = "_%s" % (jetShape)
        
        ## Loop over all distributions
        for dist_1 in dists4:
            dist = '%s%s' % (dist_1, jetShape1)
            #print "    dist: {}".format(dist)
            
            hist4[dist] = {}
            ## Loop over L1T jet algorithms
            for algo in PUSAlgosAll + PUSAlgosAllType2: # ['Raw', 'RawPUS', 'RawPUS_phiRingMin4', 'RawPUS_phiRingSide4', 'RawPUS_phiRingAdjacent']: # ['PUS','noPUS','Raw','RawPUS']:
                # read proper jetShape and PUSAlgo conbination
                if (jetShape in JetShapes      and algo not in PUSAlgosAll) or \
                   (jetShape in JetShapesType2 and algo not in PUSAlgosAllType2 ):
                    continue
                
                if (algo == 'L1JDefault') and (jetShape != 'Default'): continue
                
                hist4[dist][algo] = {}                
                
                #print "   algo : {}".format(algo)
                ## Loop over eta categories
                for iEta in ETA_CAT.keys():
                    hist4[dist][algo][iEta] = {}
                    #print "    iEta: {}".format(iEta)
                    
                    ## Loop over unpacked and emulated
                    for src in ['unp','emu']: #['unp','emu']:
                        hist4[dist][algo][iEta][src] = []
                        #print "    src: {}".format(src)

                        ## Separate histogram for each input file (different TP reconstructions)
                        for iTP in range(len(in_file_names)):
                            ## Pick binning for each histogram
                            h_bins = []
                            sAxexName = None
                            if 'jet_byHand_rates_' in dist:
                                h_bins = jetRate_bins
                                sAxexName = sAxexName_JetRate
                            else: print '\nInvalid distribution %s - no binning found!!! \nQuitting.', sys.exit()

                            #print "      dist {}, algo {}, iEta {}, src {}, iTP {}".format(dist, algo, iEta, src, iTP)

                            if 'jet_byHand_rates_' in dist:
                                # TH2D
                                hist4[dist][algo][iEta][src].append( R.TH2D( 'h_%s_%s_%s_%d_%s' % (dist, algo, iEta, iTP, src),
                                                                             sAxexName,
                                                                             int(h_bins[0]), h_bins[1], h_bins[2],
                                                                             int(PU_bins[0]),PU_bins[1],PU_bins[2] ) )


    # JetShapesType2, PUSAlgosAllType2
    dists5 = [
        #'l1TauMatchingPFJet_res_vs_iEta_vs_nVtx',
        #'jet_byHand_res_vs_iEta_vs_nVtx_l1TauDefault',
        'jet_byHand_res_vs_iEta_vs_nVtx', 
    ]
    hist5 = {}
    for jetShape in JetShapesType2:
        jetShape1 = "_%s" % (jetShape)        
    
        for dist_1 in dists5:
            dist = '%s%s' % (dist_1, jetShape1)
            print "hist5: dist {}".format(dist)
            hist5[dist] = {}

            for algo in PUSAlgosAllType2: # ['Et', 'RawEt']
                hist5[dist][algo] = {}

                ## Loop over eta categories
                for iEta in ETA_CAT.keys():
                    hist5[dist][algo][iEta] = {}

                    ## Loop over pT ranges
                    for iPt in PT_CAT.keys()+['PtAllBins']:                        
                        hist5[dist][algo][iEta][iPt] = {}

                        ## Loop over unpacked and emulated
                        for src in ['unp','emu']: #['unp','emu']:
                            hist5[dist][algo][iEta][iPt][src] = []

                            ## Separate histogram for each input file (different TP reconstructions)
                            for iTP in range(len(in_file_names)):
                                ## Pick binning for each histogram
                                h_bins = []
                                sAxexName = None
                                if 'res' in dist:
                                    h_bins = res_bins
                                else: print '\nInvalid distribution %s - no binning found!!! \nQuitting.', sys.exit()

                                if dist_1 in ['jet_byHand_res_vs_iEta_vs_nVtx']:
                                    if iEta == 'HBEF':
                                        hist5[dist][algo][iEta][iPt][src].append( R.TH3D( 'h_%s_%s_%s_%s_%d_%s' % (dist, algo, iEta, iPt, iTP, src),
                                                                                              'L1T %s %s %s in %s for %s' % (src, algo, dist, iEta, iPt),
                                                                                              iEta_bins[0],iEta_bins[1],iEta_bins[2], 101,-0.5,100.5, h_bins[0], h_bins[1], h_bins[2] ) )
    
    #print "hist5: {}".format(hist5)
    
    if PrintLevel >= 1: 
        print "Hitograms booked"
    
    # nVtx
    hnVtx       = R.TH1D("nVtx",       "nVtx",       201,-0.5,200.5);
    hnVtx_ReWtd = R.TH1D("nVtx_ReWtd", "nVtx_ReWtd", 201,-0.5,200.5);
    hPUWt       = None
    hnVtxData   = None
    hnVtxMC     = None
    hnVtxMCWtd  = None
    
    if usePUReweighting:
        print "ipFilePUData: %s \nipFilePUMC: %s" % (ipFilePUData, ipFilePUMC)
        fInPUData = R.TFile(ipFilePUData)
        fInPUMC   = R. TFile(ipFilePUMC)
        if  ( not fInPUData.IsOpen() ):
            print "PUReweighting: input file %s couldn't open" % (ipFilePUData)
            exit(0)
        if ( not fInPUMC.IsOpen() ):
            print "PUReweighting: input file %s couldn't open" % (ipFilePUMC)
            exit(0)

        hnVtxData = fInPUData.Get(sHistoPUData);
        hnVtxMC   = fInPUMC.Get(sHistoPUMC);
        if ( not hnVtxData):
            print "PUReweighting: Couldn't fetch histogram %s  from %s  \t\t ******** ERROR *********\n" % (sHistoPUData, ipFilePUData)
            exit(0)

        if ( not hnVtxMC):
            print "PUReweighting: Couldn't fetch histogram %s  from %s  \t\t ******** ERROR *********\n" % (sHistoPUMC, ipFilePUMC)
            exit(0)

        if (hnVtxData.GetNbinsX() != hnVtxMC.GetNbinsX()): 
            printf("PUReweighting: hnVtxData->GetNbinsX() != hnVtxMC->GetNbinsX() \t\t ******** ERROR *********\n");
            exit(0)

        hnVtxData.SetName("%s_Data" % (hnVtxData.GetName()));
        hnVtxMC.SetName("%s_MC" % (hnVtxMC.GetName()));

        # scale hnVtxData, hnVtxMC to have area = 1
        print "Area before scaling: hnVtxData %g, \t hnVtxMC %g " % (hnVtxData.Integral(1,hnVtxData.GetNbinsX()), hnVtxMC.Integral(1,hnVtxMC.GetNbinsX()))
        hnVtxData.Scale(1. / hnVtxData.Integral(1,hnVtxData.GetNbinsX()))
        hnVtxMC.Scale(1. / hnVtxMC.Integral(1,hnVtxMC.GetNbinsX()))    
        print "Area after scaling: hnVtxData %g, \t hnVtxMC %g " % (hnVtxData.Integral(1,hnVtxData.GetNbinsX()), hnVtxMC.Integral(1,hnVtxMC.GetNbinsX()))

        # calculate PU weights
        hPUWt = R.TH1D("hPUWeight","PU weight", hnVtxData.GetNbinsX(), hnVtxData.GetXaxis().GetXmin(), hnVtxData.GetXaxis().GetXmax())
        print "PUWeights::\n"
        for iBin in range(1, hPUWt.GetNbinsX()+1):
            nVtx = hPUWt.GetBinCenter(iBin);
            nCountsData = hnVtxData.GetBinContent(iBin);
            nCountsMC   = hnVtxMC.GetBinContent(iBin);
            PUwt = 0; 

            if nCountsData > 0 and nCountsMC > 0:
                PUwt = nCountsData / nCountsMC

            hPUWt.SetBinContent(iBin, PUwt);
            print "iBin %3d, nVtx %3g, PUwt %g, \t nCountsData %g,  nCountsMC %g  " % (iBin,nVtx,PUwt, nCountsData,nCountsMC)

        hnVtxMCWtd = hnVtxMC.Clone("%s_wtd" % (hnVtxMC.GetName()))
        for iBin in range(1, hnVtxMCWtd.GetNbinsX()+1):
            wt = hPUWt.GetBinContent(iBin);
            hnVtxMCWtd.SetBinContent(iBin, hnVtxMCWtd.GetBinContent(iBin) * wt);
    
    if PrintLevel >= 1:
        print "Read PUwgts. Now read Layer2CalibSF"

    calibSFs = OrderedDict()
    if runMode in ['CalibJetByHand']:
        fInFileCalibSF = R.TFile(sipFileCalibSF_toUse)
        if not fInFileCalibSF.IsOpen():
            print "CalibJetByHand: input file %s couldn't open" % (sipFileCalibSF_toUse)
            exit(0)
        print "CalibJetByHand: input calibration SF reading from %s " % (sipFileCalibSF_toUse)
        
        if PrintLevel >= 1:
            print "CalibJetByHand:: Calibration SFs::"
        #for jetShape in ['Default'] + JetShapes:
        for jetShape in JetShapes + JetShapesType2:
            # JetShape = "" plots are with the first version of code for 9x9 jets
            jetShape1 = jetShape
            if jetShape == 'Default':  jetShape1 = ""
            else:                      jetShape1 = "_%s" % (jetShape)
            if PrintLevel >= 1:
                print "jetShape: %s" % (jetShape)
                
            
            #calibSFs[jetShape] = {}
            calibSFs_jetShape_tmp1 = OrderedDict()
            print "0: jetShape: {}, calibSFs_jetShape_tmp1: {}".format(jetShape, calibSFs_jetShape_tmp1)
            for PUSAlgo in PUSAlgosSelected + PUSAlgosAllType2:
                
                # read proper jetShape and PUSAlgo conbination
                if (jetShape in JetShapes      and PUSAlgo not in PUSAlgosSelected ) or \
                   (jetShape in JetShapesType2 and PUSAlgo not in PUSAlgosAllType2 ):
                    continue
                
                #calibSFs[jetShape][PUSAlgo] = {}                
                if PrintLevel >= 1:
                    print "%4s%s" % (" ", PUSAlgo)

                sHistoName = sHistoCalibSF
                sHistoName = sHistoName.replace('$JETSHAPE',     jetShape1)
                sHistoName = sHistoName.replace('$PUSALGORITHM', PUSAlgo)
                if OOT_PU_scheme.lower() in ['def', 'pfa2']:
                    sHistoName = sHistoName.replace('$OOTPUS', 'PFA2')
                elif OOT_PU_scheme.lower() in ['pfa1p']:
                    sHistoName = sHistoName.replace('$OOTPUS', 'PFA1p')
               # print "CalibJetByHand: input calibration SF histogram %s" % (sHistoName)
                
                hHisto_SF_vs_iEta = fInFileCalibSF.Get(sHistoName)
                if not hHisto_SF_vs_iEta: # L1T calib SF are measured only for a few jetShape and PUSAlgorithms
                    continue
                
                print "CalibJetByHand: input calibration SF histogram %s" % (sHistoName)
                calibSFs_jetShape_tmp1[PUSAlgo] = OrderedDict()
                
                for iEta in ETA_Bins:
                    if iEta == 'HBEF': continue
                    if PrintLevel >= 1:
                        print "%8s%3s" % (" ", iEta)

                    binNumber_iEta = hHisto_SF_vs_iEta.FindBin(float(iEta))
                    SF_tmp         = hHisto_SF_vs_iEta.GetBinContent(binNumber_iEta)
                    calibSFs_jetShape_tmp1[PUSAlgo][iEta] = []
                    calibSFs_jetShape_tmp1[PUSAlgo][iEta].append([0., 999., SF_tmp]) # ['L1T pT bin low-edge',  'L1T pT bin high-edge',  'L1T layer2 calib SF']
                    
                    '''
                    sHistoName = sHistoCalibSF
                    sHistoName = sHistoName.replace('$JETSHAPE',     jetShape1)
                    sHistoName = sHistoName.replace('$PUSALGORITHM', PUSAlgo)
                    sHistoName = sHistoName.replace('$ETA',          iEta)
                    sHistoName = sHistoName.replace('$L1MODE',       'emu')
                    hProfile = fInFileCalibSF.Get(sHistoName)
                    if not hProfile:
                        print "CalibJetByHand: hProfile %s couldn't read" % (sHistoName)
                        exit(0)
                        
                    
                    calibSFs[jetShape][PUSAlgo][iEta] = []
                    for iBin in range(1, hProfile.GetNbinsX()+1):
                        binLowEdge = hProfile.GetXaxis().GetBinLowEdge(iBin)
                        binUpEdge  = hProfile.GetXaxis().GetBinUpEdge(iBin)
                        binCenter  = hProfile.GetXaxis().GetBinCenter(iBin)
                        AvgPFJetPt = hProfile.GetBinContent(iBin)
                        calibSF    = AvgPFJetPt / binCenter
                        if calibSF < 1e-3: calibSF = 1 # 
                        calibSFs[jetShape][PUSAlgo][iEta].append([binLowEdge, binUpEdge, calibSF])
                        if PrintLevel >= 1:
                            print "%12s[%g, %g]: SF = %g / %g = %g" % (" ", binLowEdge,binUpEdge, AvgPFJetPt,binCenter,calibSF)
                    '''
            
            print "1: jetShape: {}, calibSFs_jetShape_tmp1: {}".format(jetShape, calibSFs_jetShape_tmp1)
            if calibSFs_jetShape_tmp1: # SF for the corresponding jetShape exist
                calibSFs[jetShape] = calibSFs_jetShape_tmp1
                
        print "\n\nCalibJetByHand: calibSFs: {}".format(json.dumps(calibSFs, sort_keys=True, indent=4))
    
        
        
    if PrintLevel >= 1:
        print "Start event loop"

    nTotalEvents_byChains = []
    iEvt = 0
    print '\nEntering loop over %d chains' % len(chains['Unp'])
    for iCh in range(len(chains['Unp'])):

        ## Faster tecnhique, inspired by https://github.com/thomreis/l1tMuonTools/blob/master/L1Analysis.py
        Evt_br = R.L1Analysis.L1AnalysisEventDataFormat()
        Vtx_br = R.L1Analysis.L1AnalysisRecoVertexDataFormat()
        Jet_br = R.L1Analysis.L1AnalysisRecoJetDataFormat()
        Unp_br = R.L1Analysis.L1AnalysisL1UpgradeDataFormat()
        Emu_br = R.L1Analysis.L1AnalysisL1UpgradeDataFormat()
        uTP_br = R.L1Analysis.L1AnalysisCaloTPDataFormat()
        eTP_br = R.L1Analysis.L1AnalysisCaloTPDataFormat()
        uTT_br = R.L1Analysis.L1AnalysisL1CaloTowerDataFormat()
        eTT_br = R.L1Analysis.L1AnalysisL1CaloTowerDataFormat()
        uTC_br = R.L1Analysis.L1AnalysisL1CaloClusterDataFormat()
        eTC_br = R.L1Analysis.L1AnalysisL1CaloClusterDataFormat()
        
        chains['Evt'][iCh].SetBranchAddress('Event',     R.AddressOf(Evt_br))
        chains['Vtx'][iCh].SetBranchAddress('Vertex',    R.AddressOf(Vtx_br))
        chains['Jet'][iCh].SetBranchAddress('Jet',       R.AddressOf(Jet_br))
        chains['Unp'][iCh].SetBranchAddress('L1Upgrade', R.AddressOf(Unp_br))
        chains['Emu'][iCh].SetBranchAddress('L1Upgrade', R.AddressOf(Emu_br))
        chains['uTP'][iCh].SetBranchAddress('CaloTP',      R.AddressOf(uTP_br))
        chains['eTP'][iCh].SetBranchAddress('CaloTP',      R.AddressOf(eTP_br))
        chains['uTP'][iCh].SetBranchAddress('L1CaloTower', R.AddressOf(uTT_br))
        chains['eTP'][iCh].SetBranchAddress('L1CaloTower', R.AddressOf(eTT_br))
        #chains['uTP'][iCh].SetBranchAddress('L1CaloCluster', R.AddressOf(uTC_br))
        chains['eTP'][iCh].SetBranchAddress('L1CaloCluster', R.AddressOf(eTC_br))

        Ntot = chains['Unp'][iCh].GetEntries()
        print "Chain %d: nEntries %d" % (iCh, Ntot)
        nTotalEvents_byChains.append( 0 )
        for jEvt in range(Ntot):

            iEvt += 1
            nTotalEvents_byChains[iCh] += 1
            
            if jEvt > MAX_EVT and MAX_EVT > 0: break
            if iEvt % PRT_EVT is 0: print '\nEvent # %d (%dth in chain)' % (iEvt, jEvt+1); sys.stdout.flush();

            hStat.Fill(0)
            
            chains['Evt'][iCh].GetEntry(jEvt)
            chains['Vtx'][iCh].GetEntry(jEvt)
            chains['Jet'][iCh].GetEntry(jEvt)
            chains['Unp'][iCh].GetEntry(jEvt)
            chains['Emu'][iCh].GetEntry(jEvt)
            chains['uTP'][iCh].GetEntry(jEvt)
            chains['eTP'][iCh].GetEntry(jEvt)

            # ## Use these lines if you don't explicitly define the DataFormat and then do SetBranchAddress above
            # Evt_br = chains['Evt'][iCh].Event
            # Vtx_br = chains['Vtx'][iCh].Vertex
            # Jet_br = chains['Jet'][iCh].Vertex
            # Unp_br = chains['Unp'][iCh].L1Upgrade
            # Emu_br = chains['Emu'][iCh].L1Upgrade

            if VERBOSE and iEvt % PRT_EVT is 0: print '  * Run %d, LS %d, event %d, nVtx %d' % (int(Evt_br.run), int(Evt_br.lumi), int(Evt_br.event), int(Vtx_br.nVtx))

            
            puWeight = 1.0
            if usePUReweighting:
                bin_puWeight = hPUWt.FindBin(Vtx_br.nVtx);
                puWeight = hPUWt.GetBinContent(bin_puWeight);
                
                if VERBOSE:
                    print "Vtx_br.nVtx: {} puWeight: {}".format(Vtx_br.nVtx, puWeight)
            hnVtx.Fill(Vtx_br.nVtx);    
            hnVtx_ReWtd.Fill(Vtx_br.nVtx, puWeight);

            if runMode in ['makePUHisto']: continue # run quick to make PU histograms
            
            nOffJets = int(Jet_br.nJets)
            nUnpJets = int(Unp_br.nJets)
            nEmuJets = int(Emu_br.nJets)
            nUnpHTPs = int(uTP_br.nHCALTP)
            nEmuHTPs = int(eTP_br.nHCALTP)
            nUnpETPs = int(uTP_br.nECALTP)
            nEmuETPs = int(eTP_br.nECALTP)
            nUnpTTs  = int(uTT_br.nTower)
            nEmuTTs  = int(eTT_br.nTower)
            nUnpTCs  = int(uTC_br.nCluster)
            nEmuTCs  = int(eTC_br.nCluster)


            if VERBOSE or PrintLevel >= 1:
                print 'Number of jets: RECO = %d, \t L1T Unpacked = %d, Emulated = %d. \t Unpack TC: %d, TT: %d, ECALTP: %d, HCALTP: %d. \t Emulated TC: %d, TT: %d, ECALTP: %d, HCALTP: %d.' % \
                (nOffJets,
                 nUnpJets, nEmuJets,
                 nUnpTCs,nUnpTTs,nUnpETPs,nUnpHTPs,
                 nEmuTCs,nEmuTTs,nEmuETPs,nEmuHTPs)

            ### save l1jets reconstructed with different JetShape+PUS for single/double/tripple/qud-jet trigger rates ---------------------
            l1JetCollection = OrderedDict()
            for src in ['unp','emu']:
                l1JetCollection[src] = OrderedDict()
                
                #for jetShape in ['Default'] + JetShapes:
                for jetShape in JetShapes + JetShapesType2:
                    l1JetCollection[src][jetShape] = OrderedDict()
                    
                    for algo1 in PUSAlgosAll + PUSAlgosAllType2: # ['Raw', 'RawPUS', 'RawPUS_phiRingMin4', 'RawPUS_phiRingSide4', 'RawPUS_phiRingAdjacent']:
                        # read proper jetShape and PUSAlgo conbination
                        if (jetShape in JetShapes      and algo1 not in PUSAlgosAll) or \
                           (jetShape in JetShapesType2 and algo1 not in PUSAlgosAllType2 ):
                            continue
                
                        if (algo1 == 'L1JDefault') and (jetShape != 'Default'): continue
                        l1JetCollection[src][jetShape][algo1] = []
            ### ---------------------------------------------------------------------------------------------------------------------------

            
            hStat.Fill(1)
            
            ## Create list of offfline RECO jets which are too close to other RECO jets
            bad_off_jets = []
            ## Loop over all offline RECO jets
            for iOff in range(nOffJets):
                iOff_vec = R.TLorentzVector()
                iOff_vec.SetPtEtaPhiM(Jet_br.et[iOff], Jet_br.eta[iOff], Jet_br.phi[iOff], 0)
                ## Loop over all offline RECO jets with higher pT
                for jOff in range(iOff):
                    jOff_vec = R.TLorentzVector()
                    jOff_vec.SetPtEtaPhiM(Jet_br.et[jOff], Jet_br.eta[jOff], Jet_br.phi[jOff], 0)

                    if iOff_vec.DeltaR(jOff_vec) < DR_MIN:
                        # print '\n  * Removing offline jet pT = %.1f, eta = %.2f, phi = %.2f' % (iOff_vec.Pt(), iOff_vec.Eta(), iOff_vec.Phi())
                        # print '  * Has dR = %.2f to jet pT = %.1f, eta = %.2f, phi = %.2f' % (iOff_vec.DeltaR(jOff_vec), jOff_vec.Pt(), jOff_vec.Eta(), jOff_vec.Phi())
                        bad_off_jets.append(iOff)
                        break

            ## Loop over all offline RECO jets
            for iOff in range(nOffJets):

                hStat.Fill(2)
                
                ## Remove offline jets which overlap other jets
                if iOff in bad_off_jets: continue

                hStat.Fill(3)
                
                selectPFJet = True
                ## PF jet filters as recommended by Aaron: https://github.com/cms-l1t-offline/cms-l1t-analysis/blob/master/cmsl1t/filters/jets.py
                abs_eta = abs(Jet_br.eta[iOff])
                isInnerJet = abs_eta <= 2.4
                isCentralJet = abs_eta <= 2.7
                isForwardCentralJet = (abs_eta > 2.7 and abs_eta <= 3.0)
                isForwardJet = abs_eta > 3.0
                reject_if = [
                    Jet_br.muMult[iOff] != 0,
                    isCentralJet and Jet_br.nhef[iOff] >= 0.9,
                    isCentralJet and Jet_br.nemef[iOff] >= 0.9,
                    isCentralJet and (Jet_br.cMult[iOff] + Jet_br.nMult[iOff]) <= 1,
                    isCentralJet and Jet_br.mef[iOff] >= 0.8,
                    isInnerJet and Jet_br.chef[iOff] <= 0,
                    isInnerJet and Jet_br.cMult[iOff] <= 0,
                    isInnerJet and Jet_br.cemef[iOff] >= 0.9,
                    isForwardCentralJet and Jet_br.nhef[iOff] >= 0.98,
                    isForwardCentralJet and Jet_br.nemef[iOff] <= 0.01,
                    isForwardCentralJet and Jet_br.nMult[iOff] <= 2,
                    isForwardJet and Jet_br.nemef[iOff] >= 0.9,
                    isForwardJet and Jet_br.nMult[iOff] <= 10
                ]
                if any(reject_if):
                    selectPFJet = False
                    
                if PrintLevel >= 10:
                    print "    PFjet {}: selectPFJet: {}    e {}, et {}, etCorr {}, corrFactor {}, eta {}, phi {}, nCaloJets {}, caloE {}, caloEt {}, caloCorrFactor {}, caloEta {}, caloPhi {},    eEMF {}, eHadHB {}, eHadHE {}, eHadHO {}, eHadHF {},  eEmEB {}, eEmEE {}, eEmHF {}, eMaxEcalTow {}, eMaxHcalTow {}, towerArea {},   n60 {}, chef {}, nhef {}, pef {}, eef {}, mef {}, hfhef {}, hfemef {},      chMult {}, nhMult {}, phMult {}, elMult {}, muMult {}, hfhMult {}, hfemMult {}, cemef {}, cmef {}, nemef {}, cMult {}, nMult {}, ".format(\
                        iOff, selectPFJet, \
                        Jet_br.e[iOff], Jet_br.et[iOff], Jet_br.etCorr[iOff], Jet_br.corrFactor[iOff], Jet_br.eta[iOff], Jet_br.phi[iOff], \
                        Jet_br.nCaloJets, Jet_br.caloE[iOff], Jet_br.caloEt[iOff], Jet_br.caloCorrFactor[iOff], Jet_br.caloEta[iOff], Jet_br.caloPhi[iOff], \
                        Jet_br.eEMF[iOff], Jet_br.eHadHB[iOff], Jet_br.eHadHE[iOff], Jet_br.eHadHO[iOff], Jet_br.eHadHF[iOff], \
                        Jet_br.eEmEB[iOff], Jet_br.eEmEE[iOff], Jet_br.eEmHF[iOff], Jet_br.eMaxEcalTow[iOff], Jet_br.eMaxHcalTow[iOff], Jet_br.towerArea[iOff], \
                        Jet_br.n60[iOff], Jet_br.chef[iOff], Jet_br.nhef[iOff], Jet_br.pef[iOff], Jet_br.eef[iOff], Jet_br.mef[iOff], Jet_br.hfhef[iOff], Jet_br.hfemef[iOff], \
                        Jet_br.chMult[iOff], Jet_br.nhMult[iOff], Jet_br.phMult[iOff], Jet_br.elMult[iOff], Jet_br.muMult[iOff], Jet_br.hfhMult[iOff], Jet_br.hfemMult[iOff], Jet_br.cemef[iOff], Jet_br.cmef[iOff], Jet_br.nemef[iOff], Jet_br.cMult[iOff], Jet_br.nMult[iOff]
                    )
                
                if not selectPFJet: continue

                hStat.Fill(4)

                # 2018 data: HE- dead zone (HEM15/16)
                if not isMC:
                    if Evt_br.run > 319077 and \
                       Jet_br.eta[iOff] > -3.4  and Jet_br.eta[iOff] < -1.17 and \
                       Jet_br.phi[iOff] > -1.97 and Jet_br.phi[iOff] < -0.47:
                        continue
                
                hStat.Fill(5)
                
                vOff = R.TLorentzVector()
                vOff.SetPtEtaPhiM(Jet_br.etCorr[iOff], Jet_br.eta[iOff], Jet_br.phi[iOff], 0) # Aaron: use jet.etCorr instead of jet.et
                jetIEta_offlineJet     = calculateJetIEta(vOff.Eta())
                jetIEtaAbs_offlineJet  = abs(jetIEta_offlineJet)
                sjetIEta_offlineJet    = "%d" % (int(jetIEta_offlineJet))
                sjetIEtaAbs_offlineJet = "%d" % (int(jetIEtaAbs_offlineJet))
                if sjetIEta_offlineJet in hist_PFJetPt_iEtawise:
                    hist_PFJetPt_iEtawise[sjetIEta_offlineJet].Fill(vOff.Pt(), puWeight )

                data_dict = OrderedDict()
                data_dict['PFJetEtCorr'] = Jet_br.etCorr[iOff]
                #data_dict['PFJetEta']    = Jet_br.eta[iOff]
                #data_dict['PFJetPhi']    = Jet_br.phi[iOff]
                    
                ## Apply selection cut(s) to offline jet
                #if runMode_CalCalibSF == False and vOff.Pt() < PT_MIN: continue
                if runMode not in ['CalCalibSF', 'makeInputForML'] and vOff.Pt() < PT_MIN: continue
                
                hStat.Fill(6)
                
                ## Pick the |eta| and pT categories
                PFJetEtaCat = 'None'
                for iCat in ETA_CAT.keys():
                    if iCat == 'HBEF': continue
                    if abs(vOff.Eta()) > ETA_CAT[iCat][0] and abs(vOff.Eta()) < ETA_CAT[iCat][1]:
                        PFJetEtaCat = iCat
                if PFJetEtaCat == 'None' or PFJetEtaCat == 'HBEF':
                    print '\n\nSUPER-BIZZARE JET THAT FALLS INTO NO ETA CATEGORIES!!!  eta = %.3f\n\n' % vOff.Eta()
                    continue
                sEtaCat_PFJet = PFJetEtaCat
                
                iPFJetPtCat = 'None'
                for iCat in PT_CAT.keys():
                    if vOff.Pt() > PT_CAT[iCat][0] and vOff.Pt() < PT_CAT[iCat][2]:
                        iPFJetPtCat = iCat
                if iPFJetPtCat == 'None':
                    if vOff.Pt() > PT_CAT['lowPt'][0] and vOff.Pt() < PT_CAT['hiPt'][2]:
                        print '\n\nSUPER-BIZZARE JET THAT FALLS INTO NO PT CATEGORIES!!!  pT = %.3f\n\n' % vOff.Pt()
                    continue
                
                if VERBOSE or PrintLevel >= 1:
                    print "    offlineJet: eta: {}, phi: {}, et: {}, etCorr: {}".format(Jet_br.eta[iOff], Jet_br.phi[iOff], Jet_br.et[iOff], Jet_br.etCorr[iOff])

                hStat.Fill(7)
                    
                ## Find highest-pT Level-1 jet with good dR matching to unpacked jet
                max_pt = {}
                vMax   = {}
                matchedEmuIdx = {}
                for src in ['unp','emu']:
                    max_pt[src] = {}
                    vMax  [src] = {}
                    matchedEmuIdx[src] = {}
                    for algo in ['PUS','noPUS','Raw','RawPUS']:
                        max_pt[src][algo] = -99
                        vMax  [src][algo] = R.TLorentzVector()
                        matchedEmuIdx[src][algo] = -1

                ## Loop over all L1T unpacked jets
                if PrintLevel >= 10:
                    print "  * UnpJets ({}):: ".format(nUnpJets)
                for iUnp in range(nUnpJets):

                    if Unp_br.jetBx[iUnp] != 0: continue  ## Use only jets in BX 0

                    vUnp = {}
                    for algo in ['PUS','noPUS','Raw','RawPUS']:
                        vUnp[algo] = R.TLorentzVector()  ## Create a 4-vector of the L1T jet

                    vUnp['PUS']   .SetPtEtaPhiM(Unp_br.jetEt[iUnp],                           Unp_br.jetEta[iUnp], Unp_br.jetPhi[iUnp], 0)
                    vUnp['noPUS'] .SetPtEtaPhiM(Unp_br.jetEt[iUnp]    + Unp_br.jetPUEt[iUnp], Unp_br.jetEta[iUnp], Unp_br.jetPhi[iUnp], 0)
                    vUnp['Raw']   .SetPtEtaPhiM(Unp_br.jetRawEt[iUnp],                        Unp_br.jetEta[iUnp], Unp_br.jetPhi[iUnp], 0)
                    vUnp['RawPUS'].SetPtEtaPhiM(Unp_br.jetRawEt[iUnp] - Unp_br.jetPUEt[iUnp], Unp_br.jetEta[iUnp], Unp_br.jetPhi[iUnp], 0)

                    for algo in ['PUS','noPUS','Raw','RawPUS']:
                        if vUnp[algo].DeltaR(vOff) < DR_MAX and vUnp[algo].Pt() > max_pt['unp'][algo]:
                            max_pt['unp'][algo] = vUnp[algo].Pt()
                            vMax  ['unp'][algo] = vUnp[algo]
                            matchedEmuIdx['unp'][algo] = iUnp
                    
                    #hist_L1Jet_unp_TowerIEta_vs_IEta.Fill(Unp_br.jetTowerIEta[iUnp], Unp_br.jetIEta[iUnp])
                    #hist_L1Jet_unp_TowerIPhi_vs_IPhi.Fill(Unp_br.jetTowerIPhi[iUnp], Unp_br.jetIPhi[iUnp]) 
                    if PrintLevel >= 10:
                        print "    {}: jetEt {}, jetEta {}, jetPhi {}, jetIEt {}, jetIEta {}, jetIPhi {},  jetBx {}, jetTowerIPhi {}, jetTowerIEta {}, jetRawEt {}, jetSeedEt {},  jetPUEt {}, jetPUDonutEt0 {}, jetPUDonutEt1 {}, jetPUDonutEt2 {}, jetPUDonutEt3 {}".format(iUnp, \
                            Unp_br.jetEt[iUnp], Unp_br.jetEta[iUnp], Unp_br.jetPhi[iUnp], Unp_br.jetIEt[iUnp], Unp_br.jetIEta[iUnp], Unp_br.jetIPhi[iUnp], \
                            Unp_br.jetBx[iUnp], Unp_br.jetTowerIPhi[iUnp], Unp_br.jetTowerIEta[iUnp], Unp_br.jetRawEt[iUnp], Unp_br.jetSeedEt[iUnp], \
                            Unp_br.jetPUEt[iUnp], Unp_br.jetPUDonutEt0[iUnp], Unp_br.jetPUDonutEt1[iUnp], Unp_br.jetPUDonutEt2[iUnp], Unp_br.jetPUDonutEt3[iUnp]
                        )                    
                    
                ## End loop: for iUnp in range(nUnpJets)

                ## Loop over all L1T emulated jets
                if PrintLevel >= 10:
                    print "  * EmuJets ({}):: ".format(nEmuJets)                
                for iEmu in range(nEmuJets):

                    if Emu_br.jetBx[iEmu] != 0: continue  ## Use only jets in BX 0

                    vEmu = {}
                    for algo in ['PUS','noPUS','Raw','RawPUS']:
                        vEmu[algo] = R.TLorentzVector()  ## Create a 4-vector of the L1T jet

                    vEmu['PUS']   .SetPtEtaPhiM(Emu_br.jetEt[iEmu],                           Emu_br.jetEta[iEmu], Emu_br.jetPhi[iEmu], 0)
                    vEmu['noPUS'] .SetPtEtaPhiM(Emu_br.jetEt[iEmu]    + Emu_br.jetPUEt[iEmu], Emu_br.jetEta[iEmu], Emu_br.jetPhi[iEmu], 0)
                    vEmu['Raw']   .SetPtEtaPhiM(Emu_br.jetRawEt[iEmu],                        Emu_br.jetEta[iEmu], Emu_br.jetPhi[iEmu], 0)
                    vEmu['RawPUS'].SetPtEtaPhiM(Emu_br.jetRawEt[iEmu] - Emu_br.jetPUEt[iEmu], Emu_br.jetEta[iEmu], Emu_br.jetPhi[iEmu], 0)

                    for algo in ['PUS','noPUS','Raw','RawPUS']:
                        if vEmu[algo].DeltaR(vOff) < DR_MAX and vEmu[algo].Pt() > max_pt['emu'][algo]:
                            max_pt['emu'][algo] = vEmu[algo].Pt()
                            vMax  ['emu'][algo] = vEmu[algo]
                            matchedEmuIdx['emu'][algo] = iEmu

                    #hist_L1Jet_emu_TowerIEta_vs_IEta.Fill(Emu_br.jetTowerIEta[iEmu], Emu_br.jetIEta[iEmu])
                    #hist_L1Jet_emu_TowerIPhi_vs_IPhi.Fill(Emu_br.jetTowerIPhi[iEmu], Emu_br.jetIPhi[iEmu])
                    if PrintLevel >= 10:
                        print "    {}: jetEt {}, jetEta {}, jetPhi {}, jetIEt {}, jetIEta {}, jetIPhi {},  jetBx {}, jetTowerIPhi {}, jetTowerIEta {}, jetRawEt {}, jetSeedEt {},  jetPUEt {}, jetPUDonutEt0 {}, jetPUDonutEt1 {}, jetPUDonutEt2 {}, jetPUDonutEt3 {}".format(iUnp, \
                            Emu_br.jetEt[iEmu], Emu_br.jetEta[iEmu], Emu_br.jetPhi[iEmu], Emu_br.jetIEt[iEmu], Emu_br.jetIEta[iEmu], Emu_br.jetIPhi[iEmu], \
                            Emu_br.jetBx[iEmu], Emu_br.jetTowerIPhi[iEmu], Emu_br.jetTowerIEta[iEmu], Emu_br.jetRawEt[iEmu], Emu_br.jetSeedEt[iEmu], \
                            Emu_br.jetPUEt[iEmu], Emu_br.jetPUDonutEt0[iEmu], Emu_br.jetPUDonutEt1[iEmu], Emu_br.jetPUDonutEt2[iEmu], Emu_br.jetPUDonutEt3[iEmu]
                        )                                                
                ## End loop: for iEmu in range(nEmuJets)

                #print "    dR(emuJet, offlineJet):: PUS: {}, npPUS: {}, Raw: {}, RawPUS: {}".format(vMax['emu']['PUS'].DeltaR(vOff), vMax['emu']['noPUS'].DeltaR(vOff), vMax['emu']['Raw'].DeltaR(vOff), vMax['emu']['RawPUS'].DeltaR(vOff))

                
                ## Re-set the |eta| categories based on emulated and unpacked L1T jet eta, if there is a match
                etaCat = {}
                for src in ['unp','emu']:
                    etaCat[src] = {}
                    for algo in ['PUS','noPUS','Raw','RawPUS']:
                        etaCat[src][algo] = 'None'

                        if max_pt[src][algo] > 0:
                            for iCat in ETA_CAT.keys():
                                if abs(vMax[src][algo].Eta()) > ETA_CAT[iCat][0] and abs(vMax[src][algo].Eta()) < ETA_CAT[iCat][1]:
                                    etaCat[src][algo] = iCat
                        else:       etaCat[src][algo] = PFJetEtaCat

                        if etaCat[src][algo] == 'None':
                            print '\n\nSUPER-BIZZARE JET THAT FALLS INTO NO ETA CATEGORIES!!!  eta[%s][%s] = %.3f\n\n' % (vMax[src][algo].Eta(), src, algo)
                            continue
                        # if etaCat[src][algo] != PFJetEtaCat:
                        #     print '  * L1T jet (eta = %.3f) not in same category as RECO jet (eta = %.3f)' % (vMax[src][algo].Eta(), vOff.Eta())


                for src in ['unp','emu']:
                    for algo in ['PUS','noPUS','Raw','RawPUS']:
                        for jPt in PT_CAT.keys():
                            hist['jet_den'][algo][etaCat[src][algo]][jPt][src][iCh].Fill( vOff.Pt() )
                            hist['jet_den'][algo]['HBEF'           ][jPt][src][iCh].Fill( vOff.Pt() )
                            if vMax[src][algo].Pt() > PT_CAT[jPt][1]:
                                hist['jet_num'][algo][etaCat[src][algo]][jPt][src][iCh].Fill( vOff.Pt() )
                                hist['jet_num'][algo]['HBEF'           ][jPt][src][iCh].Fill( vOff.Pt() )

                        if max_pt[src][algo] > 0:
                            iL1JetPtCat = getJetPtCategory( vMax[src][algo].Pt() )
                            hist['jet_res'][algo][etaCat[src][algo]][iL1JetPtCat][src][iCh].Fill( (vMax[src][algo].Pt() - vOff.Pt()) / vOff.Pt(), puWeight )
                            #hist['jet_res'][algo][etaCat[src][algo]][iL1JetPtCat][src][iCh].Fill( vMax[src][algo].Pt() / vOff.Pt(), puWeight )
                            hist['jet_dR'] [algo][etaCat[src][algo]][iL1JetPtCat][src][iCh].Fill( (vMax[src][algo].DeltaR(vOff)), puWeight )
                            hist['jet_res'][algo]['HBEF'           ][iL1JetPtCat][src][iCh].Fill( (vMax[src][algo].Pt() - vOff.Pt()) / vOff.Pt(), puWeight )
                            #hist['jet_res'][algo]['HBEF'           ][iL1JetPtCat][src][iCh].Fill( vMax[src][algo].Pt() / vOff.Pt(), puWeight )
                            hist['jet_dR'] [algo]['HBEF'           ][iL1JetPtCat][src][iCh].Fill( (vMax[src][algo].DeltaR(vOff)), puWeight )

                    ## End loop: for algo in ['PUS','noPUS','Raw','RawPUS']
                ## End loop: for src in ['unp','emu']
                
                hStat.Fill(8)
                
                ### Compare PFjet with L1Taus --------------------------------------------------------------------------------------------------
                for src in ['emu']:
                    l1_br = None
                    if src is 'unp':
                        l1_br = Unp_br
                    elif src is 'emu':
                        l1_br = Emu_br
                    
                    hStat.Fill(9)
                    
                    ## Find highest-pT Level-1 jet with good dR matching to unpacked jet
                    matchedTau_pt = {}
                    matchedTauIdx   = {}
                    vMatchedTau = {}
                    for algo in PUSAlgosAllType2: # ['Et', 'RawEt']
                        matchedTau_pt[algo] = -99.0
                        matchedTauIdx[algo] = -1
                        vMatchedTau[algo]   = R.TLorentzVector()
                    
                    for iL1Tau in range(l1_br.nTaus):
                        
                        if l1_br.tauBx[iL1Tau] != 0: continue  ## Use only taus in BX 0
                        
                        vTau = {}
                        for algo in PUSAlgosAllType2: # ['Et', 'RawEt']
                            vTau[algo] = R.TLorentzVector()  ## Create a 4-vector of the L1T jet
                        
                        vTau['Et']   .SetPtEtaPhiM(l1_br.tauEt[iL1Tau],     l1_br.tauEta[iL1Tau], l1_br.tauPhi[iL1Tau], 0)
                        vTau['RawEt'].SetPtEtaPhiM(l1_br.tauRawEt[iL1Tau],  l1_br.tauEta[iL1Tau], l1_br.tauPhi[iL1Tau], 0)
                        
                        for algo in PUSAlgosAllType2: # ['Et', 'RawEt']
                            if vTau[algo].DeltaR(vOff) < DR_MAX and vTau[algo].Pt() > matchedTau_pt[algo]:
                                matchedTau_pt[algo] = vTau[algo].Pt()
                                matchedTauIdx[algo] = iL1Tau
                                vMatchedTau[algo] = vTau[algo]

                    jetShape = JetShapesType2[0]
                    jetShape1 = "_%s" % (JetShapesType2[0])
                    for algo in PUSAlgosAllType2: # ['Et', 'RawEt']
                        hStat.Fill(10)
                        
                        l1tau_idx = matchedTauIdx[algo]
                        if l1tau_idx < 0: continue

                        hStat.Fill(11)
                        
                        tauIEta = None
                        tauIPhi = None
                        if   src == 'emu':
                            tauIEta = l1_br.tauTowerIEta[l1tau_idx]
                            tauIPhi = l1_br.tauTowerIPhi[l1tau_idx]
                        elif src == 'unp':
                            tauIEta = convert_tauIEta_to_tauTowerIEta( l1_br.tauIEta[l1tau_idx] )
                            tauIPhi = convert_tauIPhi_to_tauTowerIPhi( l1_br.tauIPhi[l1tau_idx] )
                        tauIEtaAbs  = abs(tauIEta)
                        stauIEta    = str(tauIEta)
                        stauIEtaAbs = str(tauIEtaAbs)
                        stauIEta_toUse = stauIEtaAbs if useAbsEtaBins else stauIEta;
                        tauIEta_toUse  =  tauIEtaAbs if useAbsEtaBins else  tauIEta;
                        sL1TauEtaCat = None
                        for ieta_cat in IETA_CAT.keys():
                            if ieta_cat == 'HBEF': continue
                            if tauIEtaAbs >= IETA_CAT[ieta_cat][0] and tauIEtaAbs <= IETA_CAT[ieta_cat][1]:
                                sL1TauEtaCat = ieta_cat

                        l1jet_pt = vMatchedTau[algo].Pt()
                        l1jet_pt_woLayer2Calib = l1jet_pt
                        
                        algo1 = algo
                        sjetIEta_toUse = stauIEta_toUse
                        # l1jet_pt calibration
                        if runMode in ['CalibJetByHand'] and \
                           jetShape in calibSFs.keys() and algo1 in calibSFs[jetShape].keys() and sjetIEta_toUse in calibSFs[jetShape][algo1].keys():
                            for layer2CalibSF_list in calibSFs[jetShape][algo1][sjetIEta_toUse]:
                                # layer2CalibSF_list: [bin_pt_low, bin_pt_up, layer2CalibSF]
                                if l1jet_pt_woLayer2Calib >= layer2CalibSF_list[0] and l1jet_pt_woLayer2Calib < layer2CalibSF_list[1]:
                                    l1jet_pt = l1jet_pt_woLayer2Calib * layer2CalibSF_list[2]
                                    if PrintLevel >= 1:
                                        print "%8s l1tau pT = %g * %g = %g" % (' ',l1jet_pt_woLayer2Calib, layer2CalibSF_list[2], l1jet_pt)
                        
                        
                        res = (l1jet_pt - vOff.Pt()) / vOff.Pt()
                        
                        '''
                        print "tauIEta_toUse {}, Vtx_br.nVtx {}, res {}, puWeight {}  \t l1TauMatchingPFJet_res_vs_iEta_vs_nVtx algo {}, iEta {}', iPFJetPtCat {}, src {}, iCh {}".format(tauIEta_toUse, Vtx_br.nVtx, res, puWeight, algo,'HBEF', iPFJetPtCat,src,iCh)
                        print "jetShape1 {}, algo {}: {}".format(jetShape1, algo, 'jet_byHand_res_vs_iEta_vs_nVtx%s' % (jetShape1))
                        print "hist5: {}".format(hist5)
                        if 'jet_byHand_res_vs_iEta_vs_nVtx%s' % (jetShape1) in hist5:
                            print "hist5: {} exist".format('jet_byHand_res_vs_iEta_vs_nVtx%s' % (jetShape1))
                        if algo in hist5['jet_byHand_res_vs_iEta_vs_nVtx%s' % (jetShape1)]:
                            print "algo {} exist".format(algo)
                        if 'HBEF' in hist5['jet_byHand_res_vs_iEta_vs_nVtx%s' % (jetShape1)][algo]:
                            print "HBEF exist"
                        if iPFJetPtCat in hist5['jet_byHand_res_vs_iEta_vs_nVtx%s' % (jetShape1)][algo]['HBEF'        ]:
                            print "iPFJetPtCat {} exist".format(iPFJetPtCat)
                        if src in hist5['jet_byHand_res_vs_iEta_vs_nVtx%s' % (jetShape1)][algo]['HBEF'        ][iPFJetPtCat        ]:
                            print "src: {}, iCh: {} len(iCh) {}".format(src, iCh, len(hist5['jet_byHand_res_vs_iEta_vs_nVtx%s' % (jetShape1)][algo]['HBEF'        ][iPFJetPtCat        ][src]))
                        '''
                        
                        iL1JetPtCat = getJetPtCategory( l1jet_pt )
                        if iL1JetPtCat != 'None':
                            hist5['jet_byHand_res_vs_iEta_vs_nVtx%s' % (jetShape1)][algo]['HBEF'        ][iL1JetPtCat        ][src][iCh].Fill(tauIEta_toUse, Vtx_br.nVtx, res, puWeight)
                        hist5['jet_byHand_res_vs_iEta_vs_nVtx%s' % (jetShape1)][algo]['HBEF'        ]['PtAllBins'][src][iCh].Fill(tauIEta_toUse, Vtx_br.nVtx, res, puWeight)
                        
                        l1JetCollection[src][jetShape][algo1].append( [l1jet_pt, sL1TauEtaCat] )
                        
                        # trigger efficiency plots                        
                        for jPt in PT_CAT.keys():
                            hist3['jet_byHand_eff_den_vs_PU%s' % (jetShape1)][algo1][sEtaCat_PFJet][jPt][src][iCh].Fill( vOff.Pt(), Vtx_br.nVtx, puWeight )
                            hist3['jet_byHand_eff_den_vs_PU%s' % (jetShape1)][algo1]['HBEF'       ][jPt][src][iCh].Fill( vOff.Pt(), Vtx_br.nVtx, puWeight )
                            
                            PtTrshToUse = PT_CAT[jPt][1]
                            if jPt in PtTrshForTau:
                                PtTrshToUse = PtTrshForTau[jPt]
                            #print "jPt {}, PtTrshToUse {} \t instead of {} ".format(jPt, PtTrshToUse, PT_CAT[jPt][1])
                            
                            #if l1jet_pt > PT_CAT[jPt][1]:
                            if l1jet_pt > PtTrshToUse:
                                hist3['jet_byHand_eff_num_vs_PU%s' % (jetShape1)][algo1][sEtaCat_PFJet][jPt][src][iCh].Fill( vOff.Pt(), Vtx_br.nVtx, puWeight )
                                hist3['jet_byHand_eff_num_vs_PU%s' % (jetShape1)][algo1]['HBEF'       ][jPt][src][iCh].Fill( vOff.Pt(), Vtx_br.nVtx, puWeight )
                        
                    ### End compare PFjet with L1Taus -------------------------------------------------------------------------------------------                        


                hStat.Fill(20)        
                
                if not JetClustByHand:
                    continue
                
                # run jet clustring by hand ---------------------------------------------
                if VERBOSE or PrintLevel >= 1:
                    sTmp = "matchedEmuIdx: "
                    for src in ['unp','emu']:
                        for algo in ['PUS','noPUS','Raw','RawPUS']:
                            sTmp += "  %s_%s %d" % (src,algo,matchedEmuIdx[src][algo])
                    print "        {}".format(sTmp)

                    
                # check clusters/TTs of the emulated/unpacked jet that matched to the offline jet
                # emulated jet index: l1jet_idx = matchedEmuIdx[src]['PUS']
                #for src in ['unp', 'emu']: # ['unp','emu'] # 'unp' doesn't work as jetTOwerIEta=0 variables are stored in unpacked branch
                for src in ['emu']:
                    l1jet_br = None
                    l1TC_br  = None
                    l1TT_br  = None
                    if src is 'unp':
                        l1jet_br = Unp_br
                        l1TT_br  = uTT_br
                    elif src is 'emu':
                        l1jet_br = Emu_br
                        l1TC_br  = eTC_br
                        l1TT_br  = eTT_br

                    # use l1 jet, leading in pT with algo='PUS' that matches to offline jet, 
                    # as a reference (for jetToweriEta, jetTowerIPhi) to form cluster around (jetToweriEta, jetTowerIPhi).
                    l1jet_idx = matchedEmuIdx[src]['PUS']                    
                         
                    if l1jet_idx < 0: # no dR matching between emulated/unpacked jet and offline jet is found
                        res_dummy          = -1.49  # (l1jet_pt - vOff.Pt()) / vOff.Pt()
                        #jetIEta_offlineJet = calculateJetIEta(vOff.Eta())  # -50. # None # abs(vOff.Eta())
                        '''
                        for iEta, etaBinRange in map_iEta_Eta.items():
                            if abs(vOff.Eta()) >= etaBinRange[0] and abs(vOff.Eta()) < etaBinRange[1]:
                                jetIEta_offlineJet = float( iEta * math.copysign(1, vOff.Eta()) )
                        #print "    * Matched emulated jet not find. Offline jet eta: {}, iEta: {}".format(vOff.Eta(), jetIEta_offlineJet)
                        '''
                        
                        # fill dummy value in jet resolution histograms
                        if jetIEta_offlineJet >= -41:  # skip jetIEta_offlineJet in HF, hence not set
                            #for jetShape in ['Default'] + JetShapes:
                            for jetShape in JetShapes:
                                # JetShape = "" plots are with the first version of code for 9x9 jets
                                jetShape1 = jetShape
                                if jetShape == 'Default':  jetShape1 = ""
                                else:                      jetShape1 = "_%s" % (jetShape)
                            
                                for algo1 in PUSAlgosAll: # ['Raw', 'RawPUS', 'RawPUS_phiRingMin4', 'RawPUS_phiRingSide4', 'RawPUS_phiRingAdjacent']: # ['PUS','noPUS','Raw','RawPUS']:
                                    if (algo1 == 'L1JDefault') and (jetShape != 'Default'): continue
                                    
                                    if jetIEta_offlineJet < -41: break # jetIEta_offlineJet is in HF, hence not set

                                    jetIEta_offlineJet_tmp = jetIEtaAbs_offlineJet if useAbsEtaBins else jetIEta_offlineJet
                                    #hist2['jet_byHand_res_vs_iEta%s' % (jetShape1)][algo1]['HBEF'      ][iPFJetPtCat        ][src][iCh].Fill(jetIEta_offlineJet_tmp, res_dummy, puWeight)
                                    #hist2['jet_byHand_res_vs_iEta%s' % (jetShape1)][algo1]['HBEF'      ]['PtAllBins'][src][iCh].Fill(jetIEta_offlineJet_tmp, res_dummy, puWeight)
                        
                        continue 
                    
                    #print "src: {}, l1jet_idx: {}, l1jet_br.nJets: {}".format(src,l1jet_idx,l1jet_br.nJets)
                    jetEtPUS_L1JetDefault = l1jet_br.jetEt[l1jet_idx]
                    jetIEta = None
                    jetIPhi = None
                    if src == 'emu':
                        jetIEta = l1jet_br.jetTowerIEta[l1jet_idx]
                        jetIPhi = l1jet_br.jetTowerIPhi[l1jet_idx]
                    elif src == 'unp':
                        jetIEta = convert_jetIEta_to_jetTowerIEta( l1jet_br.jetIEta[l1jet_idx] )
                        jetIPhi = convert_jetIPhi_to_jetTowerIPhi( l1jet_br.jetIPhi[l1jet_idx] )
                    jetEt       = l1jet_br.jetEt[l1jet_idx]
                    jetIEtaAbs  = abs(jetIEta)
                    sjetIEta    = str(jetIEta)
                    sjetIEtaAbs = str(jetIEtaAbs)
                    sjetIEta_toUse = sjetIEtaAbs if useAbsEtaBins else sjetIEta;
                    jetIEta_toUse  =  jetIEtaAbs if useAbsEtaBins else  jetIEta;
                    sL1JetEtaCat = None
                    for ieta_cat in IETA_CAT.keys():
                        if ieta_cat == 'HBEF': continue
                        if jetIEtaAbs >= IETA_CAT[ieta_cat][0] and jetIEtaAbs <= IETA_CAT[ieta_cat][1]:
                            sL1JetEtaCat = ieta_cat
                    
                    
                    if VERBOSE:
                        print "        l1jet {}: jetTowerIEta: {}, jetTowerIPhi: {}, Et+PU: {}, PUet: {} RawEt: {}".format(src, l1jet_br.jetTowerIEta[l1jet_idx], l1jet_br.jetTowerIPhi[l1jet_idx], l1jet_br.jetEt[l1jet_idx] + l1jet_br.jetPUEt[l1jet_idx], l1jet_br.jetPUEt[l1jet_idx], l1jet_br.jetRawEt[l1jet_idx] )
                    if PrintLevel >= 1:
                        '''
                        print "    {}: jetEt {}, jetEta {}, jetPhi {}, jetIEt {}, jetIEta {}, jetIPhi {},   jetTowerIEta {}, jetTowerIPhi {},     jetIEta {},  jetIPhi {}".format(src, \
                            l1jet_br.jetEt[l1jet_idx], l1jet_br.jetEta[l1jet_idx], l1jet_br.jetPhi[l1jet_idx], l1jet_br.jetIEt[l1jet_idx], l1jet_br.jetIEta[l1jet_idx], l1jet_br.jetIPhi[l1jet_idx], \
                            l1jet_br.jetTowerIEta[l1jet_idx], l1jet_br.jetTowerIPhi[l1jet_idx], \
                            jetIEta,  jetIPhi                                                                                                                           
                        )
                        '''
                        print "    {} l1J({}): jetEt {}, jetEta {}, jetPhi {}, jetIEt {}, jetIEta {}, jetIPhi {},  jetBx {}, jetTowerIPhi {}, jetTowerIEta {}, jetRawEt {}, jetSeedEt {},  jetPUEt {}, jetPUDonutEt0 {}, jetPUDonutEt1 {}, jetPUDonutEt2 {}, jetPUDonutEt3 {}, @@ jetIEta {},  jetIPhi {}, jetEtPUS_L1JetDefault {}".format(src, \
                            l1jet_idx,\
                            l1jet_br.jetEt[l1jet_idx], l1jet_br.jetEta[l1jet_idx], l1jet_br.jetPhi[l1jet_idx], l1jet_br.jetIEt[l1jet_idx], l1jet_br.jetIEta[l1jet_idx], l1jet_br.jetIPhi[l1jet_idx], \
                            l1jet_br.jetBx[l1jet_idx], l1jet_br.jetTowerIPhi[l1jet_idx], l1jet_br.jetTowerIEta[l1jet_idx], l1jet_br.jetRawEt[l1jet_idx], l1jet_br.jetSeedEt[l1jet_idx], \
                            l1jet_br.jetPUEt[l1jet_idx], l1jet_br.jetPUDonutEt0[l1jet_idx], l1jet_br.jetPUDonutEt1[l1jet_idx], l1jet_br.jetPUDonutEt2[l1jet_idx], l1jet_br.jetPUDonutEt3[l1jet_idx],
                            jetIEta,  jetIPhi, jetEtPUS_L1JetDefault                                                                                                                           
                        )
                        
                    '''    
                    # check matched Calo cluster
                    nTCs = 0
                    if l1TC_br:
                        nTCs = l1TC_br.nCluster
                    for iTC in range(nTCs):
                        TCieta = l1TC_br.ieta[iTC]
                        TCiphi = l1TC_br.iphi[iTC]
                        
                        if (TCieta != jetIEta) or (TCiphi != jetIPhi): continue
                        #if abs( dIEta(TCieta, jetIEta) ) > 4 or abs( dIPhi(TCiphi, jetIPhi) ) > 4: continue
                        
                        if VERBOSE:
                            print "            l1TC: ieta: {}, iphi: {}, iet: {}".format(l1TC_br.ieta[iTC], l1TC_br.iphi[iTC], l1TC_br.iet[iTC] )
                    '''

                    '''
                    hist_nPV_vs_L1JetDefaultRAW_SF[src][sjetIEtaAbs].Fill(Vtx_br.nVtx,  l1jet_br.jetRawEt[l1jet_idx] * 0.5 / Jet_br.etCorr[iOff])
                    hist_nPV_vs_L1JetDefaultRAW_SF[src]['HBEF'     ].Fill(Vtx_br.nVtx,  l1jet_br.jetRawEt[l1jet_idx] * 0.5 / Jet_br.etCorr[iOff])
                    hist_nPV_vs_L1JetDefaultPUS_SF[src][sjetIEtaAbs].Fill(Vtx_br.nVtx,  l1jet_br.jetEt[l1jet_idx]    / Jet_br.etCorr[iOff])
                    hist_nPV_vs_L1JetDefaultPUS_SF[src]['HBEF'     ].Fill(Vtx_br.nVtx,  l1jet_br.jetEt[l1jet_idx]    / Jet_br.etCorr[iOff])
                    '''
                    
                    data_dict['L1JetType']          = src
                    data_dict['L1JetDefault_Et']    = jetEtPUS_L1JetDefault 
                    data_dict['L1JetTowerIEtaAbs']  = jetIEtaAbs 
                    #data_dict['L1JetTowerIEta'] = jetIEta
                    #data_dict['L1JetTowerIPhi'] = jetIPhi 
                    
                    ## Compare pT to sum from ECAL and HCAL TPs
                    Raw_HTP_ET  = 0
                    Raw_HTP_iET = 0
                    Raw_ETP_ET  = 0
                    Raw_ETP_iET = 0
                    Raw_TT_iET  = 0
                    Raw_TT_nET  = 0
                    PUS_TT_iET  = [0,0,0,0]
                    PUS_TT_ring = [0,0,0,0,0] # phi ring excluding central and adjacant clusters
                    PUS_TT_ring2 = [0,0,0,0,0,0,0] # phi ring excluding central cluster
                    
                    # Different shape jets
                    Raw_TT_iET_ByJetShape   = {}
                    PUS_TT_ring_ByJetShape  = {} # phi ring excluding central and adjacant clusters
                    PUS_TT_ring2_ByJetShape = {} # phi ring excluding central cluster
                    for jetShape in JetShapes:
                        if jetShape == 'Default': continue
                        Raw_TT_iET_ByJetShape[jetShape]   = 0
                        PUS_TT_ring_ByJetShape[jetShape]  = [0,0,0,0,0] # phi ring excluding central and adjacant clusters
                        PUS_TT_ring2_ByJetShape[jetShape] = [0,0,0,0,0,0,0] # phi ring excluding central cluster
                    
                    if PrintLevel >= 1 : print "    %s %s nTT %g" % (" " * 4, src, l1TT_br.nTower)    
                    for iTT in range(l1TT_br.nTower):
                        TTieta = l1TT_br.ieta[iTT]
                        TTiphi = l1TT_br.iphi[iTT]
                        TTiet  = l1TT_br.iet [iTT] * 0.5 # multiply each "TT.iet" value by 0.5, to convert to units of GeV
                        
                        dIEta_TT_Seed    = dIEta(TTieta, jetIEta)
                        dIPhi_TT_Seed    = dIPhi(TTiphi, jetIPhi)
                        AbsdIEta_TT_Seed = abs( dIEta_TT_Seed )
                        AbsdIPhi_TT_Seed = abs( dIPhi_TT_Seed )
                        
                        if PrintLevel >= 5:
                            print "    %s %s TTieta %g, TTiphi %g, TTiet %g, dIEta_TT_Seed %g, dIPhi_TT_Seed %g" % \
                                (" " * 6, src,
                                 TTieta,TTiphi,TTiet, dIEta_TT_Seed,dIPhi_TT_Seed )
                        
                        # variable shape jets with phi ring PU subtraction                       
                        for jetShape in JetShapes:
                            if jetShape == 'Default': continue
                            # jetShape:
                            # '3x9':                    jet Et = sum Et in 3x9 region aroud seed
                            # '3x9_plus_0.5_times_9x9': jet Et = (sum Et in 3x9 aroud seed) + 0.5*(sum Et in 9x9 aroud seed, but outside 3x9 central part)
                            outerClusterEtWgt = None # this is also used as a flag if cluster has central and outer parts with lower weights to outer parts
                            
                            # cluster's central part. For e.g. for jetShape = '3x9_plus_0.5_times_9x9', it read 3x9 
                            centralClustSizeInEta = int(jetShape.split('_')[0].split('x')[0])
                            centralClustSizeInPhi = int(jetShape.split('_')[0].split('x')[1]) 
                            # cluster's outer part, if that is the case.  For e.g. for jetShape = '3x9_plus_0.5_times_9x9', it read 9x9 
                            outerClustSizeInEta   = int(jetShape.split('_')[-1].split('x')[0])
                            outerClustSizeInPhi   = int(jetShape.split('_')[-1].split('x')[1])
                            if '_plus_' in jetShape:
                                outerClusterEtWgt = float(jetShape.split('_plus_')[1].split('_times_')[0])
                            
                            # consider full cluster size. i.e. including outer cluster if that is the case
                            # for case e.g. jetShape = '3x9'   or   central part for '3x9_plus_0.5_times_9x9'
                            clustSizeInEta        = centralClustSizeInEta
                            clustSizeInPhi        = centralClustSizeInPhi
                            if outerClusterEtWgt: # cluster's outer part with lower Et weights to consider. for case e.g. jetShape = '3x9_plus_0.5_times_9x9'
                                clustSizeInEta    = outerClustSizeInEta
                                clustSizeInPhi    = outerClustSizeInPhi
                            
                            clustSizeAroundSeedInPhi                                   = clustSizeInPhi / 2
                            clustSizeAroundSeedInEtaLow = clustSizeAroundSeedInEtaHigh = clustSizeInEta / 2
                            if (clustSizeInEta % 2) == 0: # even TT size cluster
                                if jetIEta > 0:
                                    clustSizeAroundSeedInEtaLow  = (clustSizeInEta / 2)
                                    clustSizeAroundSeedInEtaHigh = (clustSizeInEta / 2) - 1
                                else:
                                    clustSizeAroundSeedInEtaLow  = (clustSizeInEta / 2) - 1
                                    clustSizeAroundSeedInEtaHigh = (clustSizeInEta / 2)
                            
                            centralClustSizeAroundSeedInEtaLow = centralClustSizeAroundSeedInEtaHigh = centralClustSizeInEta / 2
                            if outerClusterEtWgt: # cluster's outer part with lower Et weights to consider:
                                if (centralClustSizeInEta % 2) == 0: # even TT size cluster
                                    if jetIEta > 0:
                                        centralClustSizeAroundSeedInEtaLow  = (centralClustSizeInEta / 2)
                                        centralClustSizeAroundSeedInEtaHigh = (centralClustSizeInEta / 2) - 1
                                    else:
                                        centralClustSizeAroundSeedInEtaLow  = (centralClustSizeInEta / 2) - 1
                                        centralClustSizeAroundSeedInEtaHigh = (centralClustSizeInEta / 2)


                            
                            if ( dIEta_TT_Seed >= -1*clustSizeAroundSeedInEtaLow and dIEta_TT_Seed <= clustSizeAroundSeedInEtaHigh ): # within cluster phi ring
                                TTiet_toUse = TTiet
                                if outerClusterEtWgt: # cluster's outer part with lower Et weights to consider:
                                    if not ( dIEta_TT_Seed >= -1*centralClustSizeAroundSeedInEtaLow and dIEta_TT_Seed <= centralClustSizeAroundSeedInEtaHigh ): # outside centralCluster phi ring
                                        TTiet_toUse = TTiet * outerClusterEtWgt
                                
                                ## "Phi ring" PU subtraction sides
                                ring_dPhi     = dIPhi_TT_Seed 
                                ring_dPhi_pos = ring_dPhi
                                if ring_dPhi < 0:
                                    ring_dPhi_pos = 72 - abs(ring_dPhi)  ## Shift into all-positive scale

                                if abs(ring_dPhi) > (4 + 9):  ## Not adjacent in phi
                                    PUS_TT_ring_ByJetShape[jetShape][(ring_dPhi_pos - 14) / 9] += TTiet_toUse  ## Fill 5 9x9 regions

                                if abs(ring_dPhi) > 4: ## ring starting from adjacent phi cluster
                                    PUS_TT_ring2_ByJetShape[jetShape][(ring_dPhi_pos - 5) / 9] += TTiet_toUse  ## Fill 7 9x9 regions


                                if ( AbsdIPhi_TT_Seed <= clustSizeAroundSeedInPhi ): # within cluster
                                    Raw_TT_iET_ByJetShape[jetShape] += TTiet_toUse
                                    if VERBOSE or PrintLevel >= 5:
                                        print "    %s %s: ieta %d, iphi %d, iet %g, iet_toUse %g,  dIEta_TT_Seed %g, dIPhi_TT_Seed %g" % \
                                        (" " * 8, jetShape, TTieta,TTiphi,TTiet,TTiet_toUse,
                                         dIEta_TT_Seed,dIPhi_TT_Seed )
                                        


                        
                        ## "Phi ring" PU subtraction sides
                        if abs( dIEta(TTieta, jetIEta) ) <= 4:  ## In the same 9-ring region as the jet
                            ring_dPhi     = dIPhi(TTiphi, jetIPhi)
                            ring_dPhi_pos = ring_dPhi
                            if ring_dPhi < 0:
                                ring_dPhi_pos = 72 - abs(ring_dPhi)  ## Shift into all-positive scale

                            if abs(ring_dPhi) > (4 + 9):  ## Not adjacent in phi
                                PUS_TT_ring[(ring_dPhi_pos - 14) / 9] += TTiet  ## Fill 5 9x9 regions
                                
                            if abs(ring_dPhi) > 4: ## ring starting from adjacent phi cluster
                                PUS_TT_ring2[(ring_dPhi_pos - 5) / 9] += TTiet  ## Fill 7 9x9 regions

                        ## "Chunky doughnut" PU subtraction sides
                        if abs( dIEta(TTieta, jetIEta) ) > 7 or  abs( dIPhi(TTiphi, jetIPhi) ) > 7: continue
                        if abs( dIEta(TTieta, jetIEta) ) > 4 and abs( dIPhi(TTiphi, jetIPhi) ) > 4: continue

                        if dIEta(TTieta, jetIEta) >  4 and dIEta(TTieta, jetIEta) <  8 and abs( dIPhi(TTiphi, jetIPhi) ) < 5:
                            PUS_TT_iET[0] += TTiet
                            # print '      ** Added to PUS[0]'
                        if dIEta(TTieta, jetIEta) < -4 and dIEta(TTieta, jetIEta) > -8 and abs( dIPhi(TTiphi, jetIPhi) ) < 5:
                            PUS_TT_iET[1] += TTiet
                            # print '      ** Added to PUS[1]'
                        if dIPhi(TTiphi, jetIPhi) >  4 and dIPhi(TTiphi, jetIPhi) <  8 and abs( dIEta(TTieta, jetIEta) ) < 5:
                            PUS_TT_iET[2] += TTiet
                            # print '      ** Added to PUS[2]'
                        if dIPhi(TTiphi, jetIPhi) < -4 and dIPhi(TTiphi, jetIPhi) > -8 and abs( dIEta(TTieta, jetIEta) ) < 5:
                            PUS_TT_iET[3] += TTiet
                            # print '      ** Added to PUS[3]'

                        ## Central jet sum
                        if abs( dIEta(TTieta, jetIEta) ) > 4 or abs( dIPhi(TTiphi, jetIPhi) ) > 4: continue
                        Raw_TT_iET += TTiet
                        Raw_TT_nET += 1
                        # print '    - Trigger tower iEta = %d, iPhi = %d, pT = %.1f' % (TTieta, TTiphi, TTiet)
                        # print '      dIEta(%d, %d) = %d' % (TTieta, jetIEta, dIEta(TTieta, jetIEta))
                        # print '      ** Added to central sum'
                        if VERBOSE:
                            print "    %sieta = %d, iphi = %d, iet = %g, iem = %g, ihad = %g, iratio = %g, iqual = %g, et = %g, eta = %g, phi = %g" % (" " * 8,eTT_br.ieta[iTT], eTT_br.iphi[iTT], eTT_br.iet[iTT], eTT_br.iem[iTT], eTT_br.ihad[iTT], eTT_br.iratio[iTT], eTT_br.iqual[iTT], eTT_br.et[iTT], eTT_br.eta[iTT], eTT_br.phi[iTT])

                    PUet                 =  sum(PUS_TT_iET)  - max(PUS_TT_iET)
                    PUS_TT_ring_Min4     = (sum(PUS_TT_ring) - max(PUS_TT_ring)) / 4.0
                    PUS_TT_ring_Side4    = (sum(PUS_TT_ring) - PUS_TT_ring[2])   / 4.0
                    PUS_TT_ring_Adjacent = (PUS_TT_ring2[0] + PUS_TT_ring2[-1]) / 2.0
                    PUS_TT_ring_Full     = sum(PUS_TT_ring2) + Raw_TT_iET
                    
                    if VERBOSE or PrintLevel >= 1:
                        print "%8s       Raw_TT_iET: %g, PUet: %g, PUS_TT_ring_Min4: %g, PUS_TT_ring_Side4: %g, PUS_TT_ring_Adjacent: %g" % (" ", Raw_TT_iET, PUet, PUS_TT_ring_Min4, PUS_TT_ring_Side4, PUS_TT_ring_Adjacent)
                        for jetShape in JetShapes:
                            if jetShape == 'Default': continue
                            PUS_TT_ring_Min4_tmp     = (sum(PUS_TT_ring_ByJetShape[jetShape]) - max(PUS_TT_ring_ByJetShape[jetShape])) / 4.0
                            PUS_TT_ring_Side4_tmp    = (sum(PUS_TT_ring_ByJetShape[jetShape]) - PUS_TT_ring_ByJetShape[jetShape][2])   / 4.0
                            PUS_TT_ring_Adjacent_tmp = (PUS_TT_ring2_ByJetShape[jetShape][0] + PUS_TT_ring2_ByJetShape[jetShape][-1]) / 2.0                            
                            #print "%8s %s: Raw_TT_iET: %g, PUet: %g, PUS_TT_ring_Min4: %g, PUS_TT_ring_Side4: %g, PUS_TT_ring_Adjacent: %g" % (" ", jetShape, Raw_TT_iET_ByJetShape[jetShape], -1, PUS_TT_ring_Min4_tmp, PUS_TT_ring_Side4_tmp, PUS_TT_ring_Adjacent_tmp)


                    #print "jetIEta {}, src {}, iCh {}".format(jetIEta, src, iCh)
                    hist1['l1jetEt_vs_RawEtMinusPU'][sjetIEta_toUse][src][iCh].Fill(jetEt, (Raw_TT_iET - PUet) / jetEt, puWeight)
                    hist1['l1jetEt_vs_RawEtMinusPU']['HBEF'        ][src][iCh].Fill(jetEt, (Raw_TT_iET - PUet) / jetEt, puWeight)
                    
                    
                    
                    
                    # jet resolution plots for different jet shapes   
                    #for jetShape in ['Default'] + JetShapes:
                    for jetShape in JetShapes:
                        # JetShape = "" plots are with the first version of code for 9x9 jets
                        jetShape1 = jetShape
                        if jetShape == 'Default':  jetShape1 = ""
                        else:                      jetShape1 = "_%s" % (jetShape)
                        
                        # assign cluster and PU energy for jetShape under consideration
                        Raw_TT_iET_tmp           = 0
                        PUet_tmp                 = 0
                        PUS_TT_ring_Min4_tmp     = 0
                        PUS_TT_ring_Side4_tmp    = 0
                        PUS_TT_ring_Adjacent_tmp = 0
                        PUS_TT_ring_Full_tmp     = 0
                        if jetShape == 'Default':
                            Raw_TT_iET_tmp           = Raw_TT_iET
                            PUet_tmp                 = PUet
                            PUS_TT_ring_Min4_tmp     = PUS_TT_ring_Min4
                            PUS_TT_ring_Side4_tmp    = PUS_TT_ring_Side4
                            PUS_TT_ring_Adjacent_tmp = PUS_TT_ring_Adjacent
                            PUS_TT_ring_Full_tmp     = PUS_TT_ring_Full
                        else:
                            Raw_TT_iET_tmp           = Raw_TT_iET_ByJetShape[jetShape]
                            PUet_tmp                 = 0
                            PUS_TT_ring_Min4_tmp     = (sum(PUS_TT_ring_ByJetShape[jetShape]) - max(PUS_TT_ring_ByJetShape[jetShape])) / 4.0 
                            PUS_TT_ring_Side4_tmp    = (sum(PUS_TT_ring_ByJetShape[jetShape]) - PUS_TT_ring_ByJetShape[jetShape][2])   / 4.0 
                            PUS_TT_ring_Adjacent_tmp = (PUS_TT_ring2_ByJetShape[jetShape][0] + PUS_TT_ring2_ByJetShape[jetShape][-1]) / 2.0
                            PUS_TT_ring_Full_tmp     = sum(PUS_TT_ring2_ByJetShape[jetShape]) + Raw_TT_iET_ByJetShape[jetShape]      
                        
                        if jetShape == 'Default':
                            data_dict['L1JetDefault_RawEtPUS']                  = Raw_TT_iET_tmp - PUet_tmp
                            data_dict['L1JetDefault_PU']                     = PUet_tmp
                        elif jetShape in ['9x9', '7x9','5x9','3x9']: # for machine learning training 
                            data_dict['L1Jet%s_RawEt' % (jetShape)]          = Raw_TT_iET_ByJetShape[jetShape]
                            data_dict['L1Jet%s_EtSum7PUTowers' % (jetShape)] = sum(PUS_TT_ring2_ByJetShape[jetShape])
                            
                        if jetShape == 'Default' and ((Raw_TT_iET_tmp - PUet_tmp) < 0.):
                        #if jetShape == 'Default':
                            #if PrintLevel >= 1:
                            #print "".format(jEvt,)
                            if jetShape == 'Default' and ((Raw_TT_iET_tmp - PUet_tmp) < 0.):
                                print "\t\t\t (L1JetDefault_RawEtPUS < 0.) for default jet \t\t\t *** ATTENTION *** "
                            print "         Event: {} ({}:{}:{}), nVtx {}, {}: Raw_TT_iET_tmp {}, PUet_tmp {}, (Raw_TT_iET_tmp - PUet_tmp) {},  jetEtPUS_L1JetDefault {}, jetIEta {}, jetIPhi {},  PFjePt {}".format(jEvt, \
                                int(Evt_br.run), int(Evt_br.lumi), int(Evt_br.event), int(Vtx_br.nVtx), \
                                src, Raw_TT_iET_tmp,PUet_tmp,(Raw_TT_iET_tmp - PUet_tmp), jetEtPUS_L1JetDefault, jetIEta,jetIPhi, vOff.Pt() )
                            print "           l1jet_idx {}: jetEt {}, jetEta {}, jetPhi {}, jetIEt {}, jetIEta {}, jetIPhi {},  jetBx {}, jetTowerIPhi {}, jetTowerIEta {}, jetRawEt/2 {}, jetSeedEt {},  jetPUEt/2 {}, jetPUDonutEt0 {}, jetPUDonutEt1 {}, jetPUDonutEt2 {}, jetPUDonutEt3 {}\n".format(l1jet_idx, \
                                l1jet_br.jetEt[l1jet_idx], l1jet_br.jetEta[l1jet_idx], l1jet_br.jetPhi[l1jet_idx], l1jet_br.jetIEt[l1jet_idx], l1jet_br.jetIEta[l1jet_idx], l1jet_br.jetIPhi[l1jet_idx], \
                                l1jet_br.jetBx[l1jet_idx], l1jet_br.jetTowerIPhi[l1jet_idx], l1jet_br.jetTowerIEta[l1jet_idx], l1jet_br.jetRawEt[l1jet_idx]/2., l1jet_br.jetSeedEt[l1jet_idx], \
                                l1jet_br.jetPUEt[l1jet_idx]/2., l1jet_br.jetPUDonutEt0[l1jet_idx], l1jet_br.jetPUDonutEt1[l1jet_idx], l1jet_br.jetPUDonutEt2[l1jet_idx], l1jet_br.jetPUDonutEt3[l1jet_idx]
                            ) 
                        
                        
                        for algo1 in PUSAlgosAll: # ['Raw', 'RawPUS', 'RawPUS_phiRingMin4', 'RawPUS_phiRingSide4', 'RawPUS_phiRingAdjacent']: # ['PUS','noPUS','Raw','RawPUS']:
                            if (algo1 == 'L1JDefault') and (jetShape != 'Default'): continue
                            
                            if Raw_TT_iET_tmp <= 0: continue # skip cluster with eT=0

                            l1jet_PU_pt = 0
                            if algo1 == 'Raw':                     l1jet_PU_pt = 0.
                            if algo1 == 'RawPUS':                  l1jet_PU_pt = PUet_tmp
                            if algo1 == 'RawPUS_phiRingMin4':      l1jet_PU_pt = PUS_TT_ring_Min4_tmp
                            if algo1 == 'RawPUS_phiRingSide4':     l1jet_PU_pt = PUS_TT_ring_Side4_tmp
                            if algo1 == 'RawPUS_phiRingAdjacent':  l1jet_PU_pt = PUS_TT_ring_Adjacent_tmp
                            if algo1 == 'RawPUS_phiDefault':       l1jet_PU_pt = Raw_TT_iET_tmp - (8.0/7*Raw_TT_iET_tmp - 1./7*PUS_TT_ring_Full_tmp) #  l1jet_pt = (8.0/7*JetRaw - 1./7*SumFullPhiRing)                            
                            l1jet_pt = Raw_TT_iET_tmp - l1jet_PU_pt
                            l1jet_pt_woLayer2Calib = l1jet_pt
                            
                            if PrintLevel >= 5 and algo1 in PUSAlgosSelected and runMode in ['CalibJetByHand']:
                                print "jetShape {}, algo1 {}, jetIEta {}".format(jetShape,algo1,sjetIEta)
                                print "  calibSFs.keys(): {}".format(calibSFs.keys())
                                if jetShape in calibSFs.keys():
                                    print "  calibSFs[jetShape].keys(): {}".format(calibSFs[jetShape].keys())
                                    if algo1 in calibSFs[jetShape].keys():
                                        print "  calibSFs[jetShape][algo1].keys(): {}".format(calibSFs[jetShape][algo1].keys())
                                        if sjetIEta_toUse in calibSFs[jetShape][algo1].keys():
                                            print "  calibSFs[jetShape][algo1][jetIEta]: {}".format(calibSFs[jetShape][algo1][sjetIEta_toUse])
                            
                            # l1jet_pt calibration
                            if runMode in ['CalibJetByHand'] and \
                               jetShape in calibSFs.keys() and algo1 in calibSFs[jetShape].keys() and sjetIEta_toUse in calibSFs[jetShape][algo1].keys():
                                for layer2CalibSF_list in calibSFs[jetShape][algo1][sjetIEta_toUse]:
                                    # layer2CalibSF_list: [bin_pt_low, bin_pt_up, layer2CalibSF]
                                    if l1jet_pt_woLayer2Calib >= layer2CalibSF_list[0] and l1jet_pt_woLayer2Calib < layer2CalibSF_list[1]:
                                        l1jet_pt = l1jet_pt_woLayer2Calib * layer2CalibSF_list[2]
                                        if PrintLevel >= 1:
                                            print "%8s l1j pT = %g * %g = %g" % (' ',l1jet_pt_woLayer2Calib, layer2CalibSF_list[2], l1jet_pt)

                            if PrintLevel >= 10 and algo1 in PUSAlgosSelected and \
                               ( (jetShape == 'Default' and algo1 == 'RawPUS') or  (jetShape != 'Default' and algo1 != 'RawPUS') ):
                                print "%4s%8s, %24s, %3s, pT: Raw: %7.2f, PU: %7.2f, l1j: %7.2f,  %7.2f,    diff: %e " % \
                                    (' ',jetShape,algo1,sjetIEta_toUse,  \
                                     Raw_TT_iET_tmp, l1jet_PU_pt, l1jet_pt_woLayer2Calib, l1jet_pt, (l1jet_pt - l1jet_pt_woLayer2Calib))


                            if PrintLevel >= 5 :
                                print "%4s%8s, %24s, %3s, pT: Raw: %7.2f, PU: %7.2f, l1j: %7.2f,  %7.2f,    diff: %e " % \
                                    (' ',jetShape,algo1,sjetIEta_toUse,  \
                                     Raw_TT_iET_tmp, l1jet_PU_pt, l1jet_pt_woLayer2Calib, l1jet_pt, (l1jet_pt - l1jet_pt_woLayer2Calib))
                            
                            
                            
                            # L1JetDefault -------------------------------------------------------
                            if (algo1 == 'L1JDefault') and (jetShape == 'Default'):
                                l1jet_pt = jetEtPUS_L1JetDefault
                            # --------------------------------------------------------------------
                            
                            iL1JetPtCat = getJetPtCategory( l1jet_pt )
                                
                                
                            '''
                            print "jetShape1: {}, algo1 {}, str(jetIEta) {}, jPt {}, src {}, iCh {}".format('jet_byHand_den%s' % (jetShape1), algo1, str(jetIEta), jPt, src, iCh)
                            if 'jet_byHand_den%s' % (jetShape1) not in hist2:
                                print "  {} doesn't exist".format('jet_byHand_den%s' % (jetShape1))
                            if algo1 not in hist2['jet_byHand_den%s' % (jetShape1)]:
                                print "  {} doesn't exist".format(algo1)
                            if str(jetIEta) not in hist2['jet_byHand_den%s' % (jetShape1)][algo1]:
                                print "  {} doesn't exist".format(str(jetIEta))
                            if 'HBEF' not in hist2['jet_byHand_den%s' % (jetShape1)][algo1]:
                                print "  {} doesn't exist".format()
                            if jPt not in hist2['jet_byHand_den%s' % (jetShape1)][algo1][str(jetIEta)]:
                                print "  {} doesn't exist".format(jPt)
                            if src not in hist2['jet_byHand_den%s' % (jetShape1)][algo1][str(jetIEta)][jPt]:
                                print "  {} doesn't exist".format(src)
                            if iCh not in hist2['jet_byHand_den%s' % (jetShape1)][algo1][str(jetIEta)][jPt][src]:
                                print "  {} doesn't exist".format(iCh)
                            print "hist2['jet_byHand_den%s' % (jetShape1)][algo1][str(jetIEta)][jPt][src]: {}".format(hist2['jet_byHand_den%s' % (jetShape1)][algo1][str(jetIEta)][jPt][src])
                            '''
                            
                            l1JetCollection[src][jetShape][algo1].append( [l1jet_pt, sL1JetEtaCat] )
                            
                            for jPt in PT_CAT.keys():
                                hist2['jet_byHand_den%s' % (jetShape1)][algo1][sjetIEta_toUse][jPt][src][iCh].Fill( vOff.Pt() )
                                hist2['jet_byHand_den%s' % (jetShape1)][algo1]['HBEF'        ][jPt][src][iCh].Fill( vOff.Pt() )
                                
                                hist3['jet_byHand_eff_den_vs_PU%s' % (jetShape1)][algo1][sEtaCat_PFJet][jPt][src][iCh].Fill( vOff.Pt(), Vtx_br.nVtx, puWeight )
                                hist3['jet_byHand_eff_den_vs_PU%s' % (jetShape1)][algo1]['HBEF'       ][jPt][src][iCh].Fill( vOff.Pt(), Vtx_br.nVtx, puWeight )
                                
                                if l1jet_pt > PT_CAT[jPt][1]:
                                    hist2['jet_byHand_num%s' % (jetShape1)][algo1][sjetIEta_toUse][jPt][src][iCh].Fill( vOff.Pt() )
                                    hist2['jet_byHand_num%s' % (jetShape1)][algo1]['HBEF'        ][jPt][src][iCh].Fill( vOff.Pt() )
                                    
                                    hist3['jet_byHand_eff_num_vs_PU%s' % (jetShape1)][algo1][sEtaCat_PFJet][jPt][src][iCh].Fill( vOff.Pt(), Vtx_br.nVtx, puWeight )
                                    hist3['jet_byHand_eff_num_vs_PU%s' % (jetShape1)][algo1]['HBEF'       ][jPt][src][iCh].Fill( vOff.Pt(), Vtx_br.nVtx, puWeight )

                            
                                    
                            res = (l1jet_pt - vOff.Pt()) / vOff.Pt()
                            #hist2['jet_byHand_res_vs_iEta%s' % (jetShape1)][algo1]['HBEF'        ][iPFJetPtCat        ][src][iCh].Fill(jetIEta_toUse, res, puWeight)
                            #hist2['jet_byHand_res_vs_iEta%s' % (jetShape1)][algo1]['HBEF'        ]['PtAllBins'][src][iCh].Fill(jetIEta_toUse, res, puWeight)

                            #hist2['jet_byHand_PU_vs_iEta%s' % (jetShape1)][algo1]['HBEF'        ][iPFJetPtCat        ][src][iCh].Fill(jetIEta_toUse, l1jet_PU_pt, puWeight)
                            #hist2['jet_byHand_PU_vs_iEta%s' % (jetShape1)][algo1]['HBEF'        ]['PtAllBins'][src][iCh].Fill(jetIEta_toUse, l1jet_PU_pt, puWeight)

                            #hist2['jet_byHand_PUByRawPt_vs_iEta%s' % (jetShape1)][algo1]['HBEF'        ][iPFJetPtCat        ][src][iCh].Fill(jetIEta_toUse, l1jet_PU_pt/Raw_TT_iET_tmp, puWeight)
                            #hist2['jet_byHand_PUByRawPt_vs_iEta%s' % (jetShape1)][algo1]['HBEF'        ]['PtAllBins'][src][iCh].Fill(jetIEta_toUse, l1jet_PU_pt/Raw_TT_iET_tmp, puWeight)
                            
                            ## Calibration plots with fixBinWidthPFJetPt
                            #hist2['jet_byHand_L1JetPt_vs_PFJetPt%s' % (jetShape1)][algo1][sjetIEta_toUse]['PtAllBins'][src][iCh].Fill(l1jet_pt, vOff.Pt(), puWeight)
                            #hist2['jet_byHand_L1JetPt_vs_PFJetPt%s' % (jetShape1)][algo1]['HBEF'        ]['PtAllBins'][src][iCh].Fill(l1jet_pt, vOff.Pt(), puWeight)

                            if runMode in ['CalibJetByHand'] and jetShape == 'Default' and algo1 == 'RawPUS':
                                # validate Layer2Calibration by hand
                                #print "jetShape1: {}, algo1 {}, str(jetIEta) {}, jPt {}, src {}, iCh {}".format('jet_byHand_den%s' % (jetShape1), algo1, str(jetIEta), jPt, src, iCh)
                                #hist2['jet_byHand_L1JetPt_vs_DefaultL1JetPt%s' % (jetShape1)][algo1][sjetIEta_toUse]['PtAllBins'][src][iCh].Fill(l1jet_pt_woLayer2Calib, (l1jet_pt - jetEt) / jetEt, puWeight)
                                #hist2['jet_byHand_L1JetPt_vs_DefaultL1JetPt%s' % (jetShape1)][algo1]['HBEF'        ]['PtAllBins'][src][iCh].Fill(l1jet_pt_woLayer2Calib, (l1jet_pt - jetEt) / jetEt, puWeight)
                                pass
                            
                            # w.r.t. nVts
                            if iL1JetPtCat != 'None':
                                hist2['jet_byHand_res_vs_iEta_vs_nVtx%s' % (jetShape1)][algo1]['HBEF'        ][iL1JetPtCat        ][src][iCh].Fill(jetIEta_toUse, Vtx_br.nVtx, res, puWeight)
                            hist2['jet_byHand_res_vs_iEta_vs_nVtx%s' % (jetShape1)][algo1]['HBEF'        ]['PtAllBins'][src][iCh].Fill(jetIEta_toUse, Vtx_br.nVtx, res, puWeight)
                    

                    #if runMode in ['makeInputForML'] and data_dict['L1JetDefault_RawEtPUS'] < 0.:
                    #    print "-ve RawEtPUS: data_dict: {}".format(data_dict)
                    if runMode in ['makeInputForML'] and isFirstEntry_WriteInputForML:
                        fOut_MLInputs_writer = csv.DictWriter(fOut_MLInputs, fieldnames=data_dict.keys())
                        fOut_MLInputs_writer.writeheader()
                        isFirstEntry_WriteInputForML = False
                        if PrintLevel >= 14:
                            print "WriteInputForML: data_dict.keys(): {}".format(data_dict.keys())
                        
                    if runMode in ['makeInputForML']: fOut_MLInputs_writer.writerow( data_dict )
                    if PrintLevel >= 14:
                        print "WriteInputForML: data_dict: {}".format(data_dict)
                    
            ## End loop: for iOff in range(nOffJets):
            
            
            ### Fill trigger rates plots per event level --------------------------------------------------------
            for src in ['unp','emu']:
                
                #for jetShape in ['Default'] + JetShapes:
                for jetShape in JetShapes + JetShapesType2:
                    # JetShape = "" plots are with the first version of code for 9x9 jets
                    jetShape1 = jetShape
                    if jetShape == 'Default':  jetShape1 = ""
                    else:                      jetShape1 = "_%s" % (jetShape)
                    
                    for algo1 in PUSAlgosAll + PUSAlgosAllType2:
                        # read proper jetShape and PUSAlgo conbination
                        if (jetShape in JetShapes      and algo1 not in PUSAlgosAll) or \
                           (jetShape in JetShapesType2 and algo1 not in PUSAlgosAllType2 ):
                            continue
                        
                        if (algo1 == 'L1JDefault') and (jetShape != 'Default'): continue
                        
                        l1JetsInEvent_sortedByL1JPt = sorted(l1JetCollection[src][jetShape][algo1], key=lambda x: x[0], reverse=True)
                        if len(l1JetsInEvent_sortedByL1JPt) == 0: continue
                        
                        if PrintLevel >= 5:
                            print "    l1JetCollection[{}][{}][{}]: {} \t\t sorted".format(src, jetShape, algo1,  l1JetCollection[src][jetShape][algo1], l1JetsInEvent_sortedByL1JPt)

                        # leading jet
                        l1JetPt_toConsider      = l1JetsInEvent_sortedByL1JPt[0][0]
                        l1JetIEtaCat_toConsider = l1JetsInEvent_sortedByL1JPt[0][1]
                        for jetBin in range(jetRate_bins[0]):
                            pT_thrsh = jetRate_bins[1] + (jetBin * jetRate_binWidth)
                            if l1JetPt_toConsider <= pT_thrsh: continue
                            hist4['jet_byHand_rates_singleJet%s' % (jetShape1)][algo1][l1JetIEtaCat_toConsider][src][iCh].Fill( pT_thrsh, Vtx_br.nVtx, puWeight )
                            hist4['jet_byHand_rates_singleJet%s' % (jetShape1)][algo1]['HBEF'                 ][src][iCh].Fill( pT_thrsh, Vtx_br.nVtx, puWeight )

                        if len(l1JetsInEvent_sortedByL1JPt) <= 1: continue
                        # 2nd leading jet
                        l1JetPt_toConsider      = l1JetsInEvent_sortedByL1JPt[1][0]
                        l1JetIEtaCat_toConsider = l1JetsInEvent_sortedByL1JPt[1][1]
                        for jetBin in range(jetRate_bins[0]):
                            pT_thrsh = jetRate_bins[1] + (jetBin * jetRate_binWidth)
                            if l1JetPt_toConsider <= pT_thrsh: continue
                            hist4['jet_byHand_rates_doubleJet%s' % (jetShape1)][algo1][l1JetIEtaCat_toConsider][src][iCh].Fill( pT_thrsh, Vtx_br.nVtx, puWeight )
                            hist4['jet_byHand_rates_doubleJet%s' % (jetShape1)][algo1]['HBEF'                 ][src][iCh].Fill( pT_thrsh, Vtx_br.nVtx, puWeight )

                        if len(l1JetsInEvent_sortedByL1JPt) <= 2: continue
                        # 3rd leading jet
                        l1JetPt_toConsider      = l1JetsInEvent_sortedByL1JPt[2][0]
                        l1JetIEtaCat_toConsider = l1JetsInEvent_sortedByL1JPt[2][1]
                        for jetBin in range(jetRate_bins[0]):
                            pT_thrsh = jetRate_bins[1] + (jetBin * jetRate_binWidth)
                            if l1JetPt_toConsider <= pT_thrsh: continue
                            hist4['jet_byHand_rates_trippleJet%s' % (jetShape1)][algo1][l1JetIEtaCat_toConsider][src][iCh].Fill( pT_thrsh, Vtx_br.nVtx, puWeight )
                            hist4['jet_byHand_rates_trippleJet%s' % (jetShape1)][algo1]['HBEF'                 ][src][iCh].Fill( pT_thrsh, Vtx_br.nVtx, puWeight )

                        if len(l1JetsInEvent_sortedByL1JPt) <= 3: continue
                        # 4th leading jet
                        l1JetPt_toConsider      = l1JetsInEvent_sortedByL1JPt[3][0]
                        l1JetIEtaCat_toConsider = l1JetsInEvent_sortedByL1JPt[3][1]
                        for jetBin in range(jetRate_bins[0]):
                            pT_thrsh = jetRate_bins[1] + (jetBin * jetRate_binWidth)
                            if l1JetPt_toConsider <= pT_thrsh: continue
                            hist4['jet_byHand_rates_quadJet%s' % (jetShape1)][algo1][l1JetIEtaCat_toConsider][src][iCh].Fill( pT_thrsh, Vtx_br.nVtx, puWeight )
                            hist4['jet_byHand_rates_quadJet%s' % (jetShape1)][algo1]['HBEF'                 ][src][iCh].Fill( pT_thrsh, Vtx_br.nVtx, puWeight )
                ###  trigger rates plots per event level --------------------------------------------------------                
                

        ## End loop: for jEvt in range(chains['Unp'][iCh].GetEntries()):
        
        print "\n\n nTotalEvents_byChains[iCh {}]: {} ".format(iCh, nTotalEvents_byChains[iCh])
        hnTotalEvents.SetBinContent(iCh+1, nTotalEvents_byChains[iCh])
    ## End loop: for iCh in range(len(chains['Unp'])):

    print '\nFinished loop over chains'
    if runMode in ['makeInputForML']: fOut_MLInputs.close()
    
    #out_file_str = "tmp.root"
    out_file = R.TFile(out_file_str,'recreate')
    out_file.cd()

    colors = [R.kGray, R.kViolet, R.kBlue, R.kCyan+3, R.kTeal, R.kSpring, R.kRed, R.kOrange, R.kMagenta+3, R.kMagenta]
    while len(colors) < len(in_file_names):
        colors += colors

    if usePUReweighting:
        dir1 = out_file.mkdir("PUWeights")
        dir1.cd()
        hPUWt.Write();
        hnVtxData.Write();
        hnVtxMC.Write();
        hnVtxMCWtd.Write();
    
    out_file.cd()
    hnVtx.Write();
    hnVtx_ReWtd.Write();
    hnTotalEvents.Write();
    
    
    '''
    hist_L1Jet_unp_TowerIEta_vs_IEta.Write()
    hist_L1Jet_emu_TowerIEta_vs_IEta.Write()
    hist_L1Jet_unp_TowerIPhi_vs_IPhi.Write()
    hist_L1Jet_emu_TowerIPhi_vs_IPhi.Write()
    
    for iEta in ETA_Bins:
        hist_PFJetPt_iEtawise[iEta].Write()

    for src in ['unp','emu']:    
        for iEta in ETA_Bins:
            hist_nPV_vs_L1JetDefaultRAW_SF[src][iEta].Write()
            hist_nPV_vs_L1JetDefaultPUS_SF[src][iEta].Write()
    '''


    '''
    ## Loop over all histograms
    if runMode not in ['CalCalibSF', 'CalibJetByHand', 'makeInputForML']:
        out_file.cd()
        for algo in ['PUS','noPUS','Raw','RawPUS']:
            for iEta in ETA_CAT.keys():
                for iPt in PT_CAT.keys():
                    for src in ['unp','emu']:
                        for iTP in range(len(in_file_names)):

                            hist['jet_eff'][algo][iEta][iPt][src][iTP] = R.TEfficiency( hist['jet_num'][algo][iEta][iPt][src][iTP],
                                                                                        hist['jet_den'][algo][iEta][iPt][src][iTP] )
                            hist['jet_eff'][algo][iEta][iPt][src][iTP].SetName( hist['jet_num'][algo][iEta][iPt][src][iTP].GetName().replace('num','eff') )

                            for dist in dists:

                                hist[dist][algo][iEta][iPt][src][iTP].SetLineWidth(2)
                                if src == 'unp': hist[dist][algo][iEta][iPt][src][iTP].SetLineColor(R.kBlack)
                                if src == 'emu': hist[dist][algo][iEta][iPt][src][iTP].SetLineColor(colors[iTP])
                                hist[dist][algo][iEta][iPt][src][iTP].Write()

                            ## End loop: for dist in dists+['jet_eff']
                        ## End loop: for iTP in range(len(in_file_names))
                    ## End loop: for src in ['unp','emu']
                ## End loop: for iPt in PT_CAT.keys()
            ## End loop: for iEta in ETA_CAT.keys()
        ## End loop: for algo in ['PUS','noPUS','Raw','RawPUS']
        
    if runMode not in ['CalCalibSF', 'CalibJetByHand', 'makeInputForML']:
        out_file.cd()
        for dist in dists1:
            for iEta in ETA_Bins:
                for src in ['unp', 'emu']:
                    for iTP in range(len(in_file_names)):
                        hist1[dist][iEta][src][iTP].SetLineWidth(2)
                        if src == 'unp': hist1[dist][iEta][src][iTP].SetLineColor(R.kBlack)
                        if src == 'emu': hist1[dist][iEta][src][iTP].SetLineColor(colors[iTP])
                        hist1[dist][iEta][src][iTP].Write()
    '''
    
    ## Loop over all histograms
    # hist2
    if runMode not in ['CalibJetByHand', 'makeInputForML']:
        out_file.cd()        
        for algo in PUSAlgosAll: # ['Raw', 'RawPUS', 'RawPUS_phiRingMin4', 'RawPUS_phiRingSide4', 'RawPUS_phiRingAdjacent']: # ['PUS','noPUS','Raw','RawPUS']:        
            for iEta in ETA_Bins:
                for iPt in PT_CAT.keys() + ['PtAllBins']:
                    for src in ['unp', 'emu']: #['unp','emu']:
                        for iTP in range(len(in_file_names)):
                            #for jetShape in ['Default'] + JetShapes:
                            for jetShape in JetShapes:
                                if (algo == 'L1JDefault') and (jetShape != 'Default'): continue

                                # JetShape = "" plots are with the first version of code for 9x9 jets
                                jetShape1 = jetShape
                                if jetShape == 'Default':  jetShape1 = ""
                                else:                      jetShape1 = "_%s" % (jetShape)
                                 # JetShape = "" plots are with the first version of code for 9x9 jets                        

                                hist2['jet_byHand_eff%s' % (jetShape1)][algo][iEta][iPt][src][iTP] = R.TEfficiency( hist2['jet_byHand_num%s' % (jetShape1)][algo][iEta][iPt][src][iTP],
                                                                                                    hist2['jet_byHand_den%s' % (jetShape1)][algo][iEta][iPt][src][iTP] )
                                hist2['jet_byHand_eff%s' % (jetShape1)][algo][iEta][iPt][src][iTP].SetName( hist2['jet_byHand_num%s' % (jetShape1)][algo][iEta][iPt][src][iTP].GetName().replace('num','eff') )

                                for dist_1 in dists2:
                                    dist = '%s%s' % (dist_1, jetShape1)
                                    if '_vs_iEta' in dist and iEta != 'HBEF': continue
                                    if dist_1 in ['jet_byHand_L1JetPt_vs_DefaultL1JetPt'] and algo != 'RawPUS': continue
                                    if dist_1 in ['jet_byHand_L1JetPt_vs_PFJetPt', 'jet_byHand_L1JetPt_vs_DefaultL1JetPt'] and iPt != 'PtAllBins': continue
                                    if dist_1 in ['jet_byHand_L1JetPt_vs_DefaultL1JetPt'] and jetShape != 'Default': continue

                                    hist2[dist][algo][iEta][iPt][src][iTP].SetLineWidth(2)
                                    if src == 'unp': hist2[dist][algo][iEta][iPt][src][iTP].SetLineColor(R.kBlack)
                                    if src == 'emu': hist2[dist][algo][iEta][iPt][src][iTP].SetLineColor(colors[iTP])
                                    hist2[dist][algo][iEta][iPt][src][iTP].Write()
    
    
    
    if runMode in ['CalibJetByHand']:
        # hist3
        out_file.cd()        
        for algo in PUSAlgosAll + PUSAlgosAllType2: # ['Raw', 'RawPUS', 'RawPUS_phiRingMin4', 'RawPUS_phiRingSide4', 'RawPUS_phiRingAdjacent']: # ['PUS','noPUS','Raw','RawPUS']:
            for iEta in ETA_CAT.keys():
                for iPt in PT_CAT.keys():
                    for src in ['unp', 'emu']: #['unp','emu']:
                        for iTP in range(len(in_file_names)):
                            #for jetShape in ['Default'] + JetShapes:
                            for jetShape in JetShapes + JetShapesType2:
                                # read proper jetShape and PUSAlgo conbination
                                if (jetShape in JetShapes      and algo not in PUSAlgosAll) or \
                                   (jetShape in JetShapesType2 and algo not in PUSAlgosAllType2 ):
                                    continue
                
                                if (algo == 'L1JDefault') and (jetShape != 'Default'): continue

                                # JetShape = "" plots are with the first version of code for 9x9 jets
                                jetShape1 = jetShape
                                if jetShape == 'Default':  jetShape1 = ""
                                else:                      jetShape1 = "_%s" % (jetShape)
                                 # JetShape = "" plots are with the first version of code for 9x9 jets   

                                for dist_1 in dists3:
                                    dist = '%s%s' % (dist_1, jetShape1)

                                    hist3[dist][algo][iEta][iPt][src][iTP].SetLineWidth(2)
                                    if src == 'unp': hist3[dist][algo][iEta][iPt][src][iTP].SetLineColor(R.kBlack)
                                    if src == 'emu': hist3[dist][algo][iEta][iPt][src][iTP].SetLineColor(colors[iTP])
                                    hist3[dist][algo][iEta][iPt][src][iTP].Write()

        # hist4
        out_file.cd()        
        for algo in PUSAlgosAll + PUSAlgosAllType2: # ['Raw', 'RawPUS', 'RawPUS_phiRingMin4', 'RawPUS_phiRingSide4', 'RawPUS_phiRingAdjacent']: # ['PUS','noPUS','Raw','RawPUS']:
            for iEta in ETA_CAT.keys():
                for src in ['unp', 'emu']: #['unp','emu']:
                    for iTP in range(len(in_file_names)):
                        #for jetShape in ['Default'] + JetShapes:
                        for jetShape in JetShapes + JetShapesType2:
                            # read proper jetShape and PUSAlgo conbination
                            if (jetShape in JetShapes      and algo not in PUSAlgosAll) or \
                               (jetShape in JetShapesType2 and algo not in PUSAlgosAllType2 ):
                                continue
                            
                            if (algo == 'L1JDefault') and (jetShape != 'Default'): continue
                            
                            # JetShape = "" plots are with the first version of code for 9x9 jets
                            jetShape1 = jetShape
                            if jetShape == 'Default':  jetShape1 = ""
                            else:                      jetShape1 = "_%s" % (jetShape)
                             # JetShape = "" plots are with the first version of code for 9x9 jets   

                            for dist_1 in dists4:
                                dist = '%s%s' % (dist_1, jetShape1)

                                hist4[dist][algo][iEta][src][iTP].SetLineWidth(2)
                                if src == 'unp': hist4[dist][algo][iEta][src][iTP].SetLineColor(R.kBlack)
                                if src == 'emu': hist4[dist][algo][iEta][src][iTP].SetLineColor(colors[iTP])
                                hist4[dist][algo][iEta][src][iTP].Write()


    # JetShapesType2, PUSAlgosAllType2
    # hist5
    if runMode not in ['makeInputForML']:
        out_file.cd()
        hStat.Fill(50)
        
        for iEta in ETA_CAT.keys():
            for iPt in PT_CAT.keys()+['PtAllBins']:
                for src in ['emu']: # ['unp', 'emu']:
                    for iTP in range(len(in_file_names)):
                        for jetShape in JetShapesType2:
                            jetShape1 = "_%s" % (jetShape)      
                            
                            for algo in PUSAlgosAllType2: # ['Et', 'RawEt']
                            
                                for dist_1 in dists5:
                                    dist = '%s%s' % (dist_1, jetShape1)
                                    

                                    if dist_1 in ['jet_byHand_res_vs_iEta_vs_nVtx'] and \
                                       iEta != 'HBEF':
                                        continue

                                    #print "dist5 {} algo {}, iEta {}', iPt {}, src {}, iTP {}".format(dist, algo,iEta, iPt,src,iTP)
                                    hist5[dist][algo][iEta][iPt][src][iTP].Write()

    hStat.Write()
    
    out_file.Close()
    del chains

    #print '\nWrote out file:  plots/'+out_file_str+'.root'
    print '\nWrote out file:  '+out_file_str


if __name__ == '__main__':
    run()
