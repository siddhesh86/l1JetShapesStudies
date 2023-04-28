#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
import argparse
from collections import OrderedDict
from collections import OrderedDict as OD
import csv

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib
import hist
import mplhep as hep

from xgboost import XGBRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from hyperopt import fmin, tpe, hp, STATUS_OK, Trials, space_eval # hyperparameter optimization
import pickle


from CommonTools import (
    convert_Et_to_logEt, 
    convert_logEt_to_Et, 
    convert_GenEt_to_GenEtByL1Et, 
    convert_GenEtByL1Et_to_GenEt, 
    convert_GenEt_to_logGenEtByL1Et, 
    convert_logGenEtByL1Et_to_GenEt,
    convert_CaloToolMPEta_to_IEta,
    GaussianFunction,
    calculate_errorOfRatio
)


parser = argparse.ArgumentParser()
parser.add_argument('--MLTarget',           type=str,   dest='MLTarget',  default='logGenEt', choices=['GenEt', 'logGenEt', 'GenEtByL1Et', 'logGenEtByL1Et'])
parser.add_argument('--PUForSFComputation', type=int,   dest='PUForSFComputation', help="PU at which SFs to compute", default='48')
parser.add_argument('--fracOfDataToUse',    type=float, dest='fracOfDataToUse', help="fraction of data to use", default='1.0')
parseGroup1 = parser.add_mutually_exclusive_group(required=True)
parseGroup1.add_argument('--ChunkyDonut',    action='store_true')
parseGroup1.add_argument('--PhiRing',        action='store_true')
parseGroup2 = parser.add_mutually_exclusive_group(required=True)
parseGroup2.add_argument('--l1MatchOffline', action='store_true')
parseGroup2.add_argument('--l1MatchGen',     action='store_true')

runLocally = False

args = None
if not runLocally:
    matplotlib.use('Agg') # use for condor jobs to disable display of plots
    args = parser.parse_args()
else:
    #args = parser.parse_args("--ChunkyDonut --l1MatchGen --MLTarget logGenEt --fracOfDataToUse 0.01".split()) # to run in jupyter-notebook     
    args = parser.parse_args("--PhiRing --l1MatchGen --MLTarget logGenEtByL1Et --fracOfDataToUse 1.00 --PUForSFComputation 33 ".split()) # to run in jupyter-notebook     
    from IPython.display import display, HTML
    display(HTML("<style>.container { width:100% !important; }</style>"))

l1Jet_ChunkyDonut = args.ChunkyDonut
l1Jet_PhiRing     = args.PhiRing
l1MatchOffline    = args.l1MatchOffline
l1MatchGen        = args.l1MatchGen
fracOfDataToUse   = args.fracOfDataToUse
MLTarget          = args.MLTarget
PUForSFComputation = args.PUForSFComputation



printLevel = PrintLevel = 5
iEtaBins = [i for i in range(1, 42) if i!=29]
sL1JetEt_PUS_ChunkyDonut = 'L1JetEt_PUS_ChunkyDonut'
sL1JetEt_PUS_PhiRing     = 'L1JetEt_PUS_PhiRing'
sOfflineJetEt            = 'PFJetEtCorr'
sGenJetEt                = 'GenJetEt'
sL1JetTowerIEtaAbs       = 'L1JetTowerIEtaAbs'
L1JetPtThrsh             = 10.0 # GeV
L1JetPtMax               = 255.0 # GeV
RefJetPtThrsh            = 10.0 # GeV
snVtx                    = 'nVertexReco'

NCompPtBins = 16 # 16 # No. of compressed pT bins
calibSF_L1JetPtRange = [15., 255., 1.] # [<lowest pT>,  <hightest pT>,  <pT bin width>] # pT range for SFs to read from Syed's SF.csv file
LUT_PtRange = [0., 255., 1.] # pT range for SFs for LUT
SF_forZeroPt = 1.0



sL1JetEt  = sL1JetEt_PUS_ChunkyDonut if l1Jet_ChunkyDonut else sL1JetEt_PUS_PhiRing
sRefJetEt = sOfflineJetEt if l1MatchOffline else sGenJetEt 


#version         = "v%s_%s_MLTarget_%s_dataFrac%.2f_20220925" % (sL1JetEt, sRefJetEt, MLTarget, fracOfDataToUse) 
version         = "v%s_%s_MLTarget_%s_dataFrac%.2f_20220925_wRefJetPtHighThrsh400GeV" % (sL1JetEt, sRefJetEt, MLTarget, fracOfDataToUse) 
sIpFileName     = "../data/L1T_Jet_MLInputs_Run3_QCD_Pt15to7000_PFA1p_CMSSW12_6_0_pre1_nVtxAll_20220925.csv"
sOpFileName_SFs = "../data/L1T_Jet_SFs_2023_QCD_122X_mcRun3_2021_realistic_v9_12_6_0_pre1_20220925_%s.csv" % (version)
sOutDir         = "./plots_check_BDT_performance_%s" % (version)


#sOpFileName_SFs = sOpFileName_SFs.replace('.csv', '_%s.csv' % (sL1JetEt))
#sOutDir = '%s_%s' % (sOutDir, sL1JetEt)


BDTFileNames_dict = {
    sL1JetEt_PUS_ChunkyDonut: {
        'GenEt': {
            'HBEF': {
                #'PtAll': '../data/BDTModel_vL1JetEt_PUS_ChunkyDonut_GenJetEt_MLTarget_GenEt_dataFrac1.00_20230403_wRefJetPtHighThrsh999999GeV_wOptimizedHyperparams_HBEF_PtAll.pkl',                
                'PtAll': '../data/BDTModel_vL1JetEt_PUS_ChunkyDonut_GenJetEt_MLTarget_GenEt_dataFrac1.00_20220925_wRefJetPtHighThrsh400GeV_wOptimizedHyperparams_HBEF_PtAll.pkl',                
            },            
        },
        
        'logGenEt': {
            'HBEF': {
                #'PtAll': '../data/BDTModel_vL1JetEt_PUS_ChunkyDonut_GenJetEt_MLTarget_logGenEt_dataFrac1.00_20230403_wRefJetPtHighThrsh999999GeV_wOptimizedHyperparams_HBEF_PtAll.pkl',                
                'PtAll': '../data/BDTModel_vL1JetEt_PUS_ChunkyDonut_GenJetEt_MLTarget_logGenEt_dataFrac1.00_20220925_wRefJetPtHighThrsh400GeV_wOptimizedHyperparams_HBEF_PtAll.pkl',                
            },            
        },
        
        'GenEtByL1Et': {
            'HBEF': {
                #'PtAll': '../data/BDTModel_vL1JetEt_PUS_ChunkyDonut_GenJetEt_MLTarget_GenEtByL1Et_dataFrac1.00_20230403_wRefJetPtHighThrsh999999GeV_wOptimizedHyperparams_HBEF_PtAll.pkl',                
                'PtAll': '../data/BDTModel_vL1JetEt_PUS_ChunkyDonut_GenJetEt_MLTarget_GenEtByL1Et_dataFrac1.00_20220925_wRefJetPtHighThrsh400GeV_wOptimizedHyperparams_HBEF_PtAll.pkl',                
            },            
        },
        
        'logGenEtByL1Et': {
            'HBEF': {
                #'PtAll': '../data/BDTModel_vL1JetEt_PUS_ChunkyDonut_GenJetEt_MLTarget_logGenEtByL1Et_dataFrac1.00_20230403_wRefJetPtHighThrsh999999GeV_wOptimizedHyperparams_HBEF_PtAll.pkl',                
                'PtAll': '../data/BDTModel_vL1JetEt_PUS_ChunkyDonut_GenJetEt_MLTarget_logGenEtByL1Et_dataFrac1.00_20220925_wRefJetPtHighThrsh400GeV_wOptimizedHyperparams_HBEF_PtAll.pkl',                
            },            
        },        
    },
    
    sL1JetEt_PUS_PhiRing: {
        'GenEt': {
            'HBEF': {
                #'PtAll': '../data/BDTModel_vL1JetEt_PUS_PhiRing_GenJetEt_MLTarget_GenEt_dataFrac1.00_20230403_wRefJetPtHighThrsh999999GeV_wOptimizedHyperparams_HBEF_PtAll.pkl',                
                'PtAll': '.../data/BDTModel_vL1JetEt_PUS_PhiRing_GenJetEt_MLTarget_GenEt_dataFrac1.00_20220925_wRefJetPtHighThrsh400GeV_wOptimizedHyperparams_HBEF_PtAll.pkl',                
            },            
        },
        
        'logGenEt': {
            'HBEF': {
                #'PtAll': '../data/BDTModel_vL1JetEt_PUS_PhiRing_GenJetEt_MLTarget_logGenEt_dataFrac1.00_20230403_wRefJetPtHighThrsh999999GeV_wOptimizedHyperparams_HBEF_PtAll.pkl',                
                'PtAll': '../data/BDTModel_vL1JetEt_PUS_PhiRing_GenJetEt_MLTarget_logGenEt_dataFrac1.00_20220925_wRefJetPtHighThrsh400GeV_wOptimizedHyperparams_HBEF_PtAll.pkl',                
            },            
        },
        
        'GenEtByL1Et': {
            'HBEF': {
                #'PtAll': '../data/BDTModel_vL1JetEt_PUS_PhiRing_GenJetEt_MLTarget_GenEtByL1Et_dataFrac1.00_20230403_wRefJetPtHighThrsh999999GeV_wOptimizedHyperparams_HBEF_PtAll.pkl',                
                'PtAll': '../data/BDTModel_vL1JetEt_PUS_PhiRing_GenJetEt_MLTarget_GenEtByL1Et_dataFrac1.00_20220925_wRefJetPtHighThrsh400GeV_wOptimizedHyperparams_HBEF_PtAll.pkl',                
            },            
        },
        
        'logGenEtByL1Et': {
            'HBEF': {
                #'PtAll': '../data/BDTModel_vL1JetEt_PUS_PhiRing_GenJetEt_MLTarget_logGenEtByL1Et_dataFrac1.00_20230403_wRefJetPtHighThrsh999999GeV_wOptimizedHyperparams_HBEF_PtAll.pkl',                
                'PtAll': '../data/BDTModel_vL1JetEt_PUS_PhiRing_GenJetEt_MLTarget_logGenEtByL1Et_dataFrac1.00_20220925_wRefJetPtHighThrsh400GeV_wOptimizedHyperparams_HBEF_PtAll.pkl',                
            },            
        },        
    },    
    
}
PUPointsForMLEvaluation = [33, 48, 63]
PUPointForMLEvaluation_selected = PUForSFComputation #48

plotPerformancePlots = False

PT_CAT = OD()
PT_CAT['Ptlt25']   = [ 0,  15,   25]  ## Low pT, turn-on threshold, high pT
PT_CAT['Pt25To35'] = [25,  30,   35]  ## Low pT, turn-on threshold, high pT
PT_CAT['Pt35To60'] = [35,  55,   60]  ## Low pT, turn-on threshold, high pT
PT_CAT['Pt60To90'] = [60,  75,   90]  ## Low pT, turn-on threshold, high pT   #[60,  90,   90]
PT_CAT['Ptgt90']   = [90, 120, 9999]  ## Low pT, turn-on threshold, high pT

ETA_CAT = OD()
ETA_CAT['HBEF'] = [0.000, 5.210]  ## Whole detector, 1 - 41
ETA_CAT['HB']   = [0.000, 1.392]  ## Trigger towers  1 - 16
ETA_CAT['HE1']  = [1.392, 1.740]  ## Trigger towers 17 - 20
ETA_CAT['HE2a'] = [1.740, 2.322]  ## Trigger towers 21 - 25
ETA_CAT['HE2b'] = [2.322, 3.000]  ## Trigger towers 26 - 28
ETA_CAT['HF']   = [3.000, 5.210]  ## Trigger towers 30 - 41

IETA_CAT = OD()
IETA_CAT['HBEF'] = [ 1, 41]  ## Whole detector, 1 - 41
IETA_CAT['HB']   = [ 1, 16]  ## Trigger towers  1 - 16
IETA_CAT['HE1']  = [17, 20]  ## Trigger towers 17 - 20
IETA_CAT['HE2a'] = [21, 25]  ## Trigger towers 21 - 25
IETA_CAT['HE2b'] = [26, 28]  ## Trigger towers 26 - 28
IETA_CAT['HF']   = [30, 41]  ## Trigger towers 30 - 41

CaloToolMPEtaBinsMerge_forEtaCompressedLUT_2018 = OrderedDict([
    ( 0, [*range(1, 5+1)]),  # 0
    ( 1, [*range(6, 9+1)]),  # 1
    ( 2, [*range(10,13+1)]), # 2
    ( 3, [*range(14,15+1)]), # 3
    ( 4, [*range(16,17+1)]), # 4
    ( 5, [*range(18,19+1)]), # 5
    ( 6, [*range(20,21+1)]), # 6
    ( 7, [22]), # 7
    ( 8, [23]), # 8
    ( 9, [24]), # 9
    (10, [25]), # 10
    (11, [26]), # 11
    (12, [*range(27,28+1)]), # 12
    (13, [*range(29,31+1)]), # 13
    (14, [*range(32,35+1)]), # 14
    (15, [*range(36,40+1)]), # 15
])
CaloToolMPEtaBinsMerge_forEtaCompressedLUT = CaloToolMPEtaBinsMerge_forEtaCompressedLUT_2018

useAbsEtaBins = True
ETA_Bins = []
for iEta in range(-41,42):
    if iEta in [-29, 0, 29]:        continue;
    if useAbsEtaBins and iEta < 0:  continue;
    ETA_Bins.append(str(iEta))


sOutDirBeforeJEC    = "%s/beforeJEC" % (sOutDir)
sOutDirAfterJEC     = "%s/afterJEC" % (sOutDir)
sOutDirBeforeJEC_1D = '%s/1D' % (sOutDirBeforeJEC)
sOutDirAfterJEC_1D  = '%s/1D' % (sOutDirAfterJEC)
if not os.path.exists(sOutDir):             os.makedirs( sOutDir , exist_ok=True )
if not os.path.exists(sOutDirBeforeJEC):    os.makedirs( sOutDirBeforeJEC , exist_ok=True )
if not os.path.exists(sOutDirAfterJEC):     os.makedirs( sOutDirAfterJEC , exist_ok=True )    
if not os.path.exists(sOutDirBeforeJEC_1D): os.makedirs( sOutDirBeforeJEC_1D , exist_ok=True )    
if not os.path.exists(sOutDirAfterJEC_1D):  os.makedirs( sOutDirAfterJEC_1D , exist_ok=True )
if not os.path.exists("../data"):           os.makedirs("../data", exist_ok=True )
    
print("Input file: %s" % (sIpFileName))
print(f"{fracOfDataToUse = }")
print(f"{PUForSFComputation = }")
print(f"{version = }")
print("iEtaBins ({}): {}".format(len(iEtaBins), iEtaBins))
print("sRefJetEt: {}, \t sL1Jet: {}, \t L1JetPtThrsh: {}".format(sRefJetEt, sL1JetEt, L1JetPtThrsh))
print("l1Jet_ChunkyDonut {}, l1Jet_PhiRing {}, l1MatchOffline {}, l1MatchGen {}".format(
    l1Jet_ChunkyDonut, l1Jet_PhiRing, l1MatchOffline, l1MatchGen)); sys.stdout.flush();    


# In[2]:


data_all = pd.read_csv(sIpFileName)


# In[3]:


data_all[sL1JetEt_PUS_ChunkyDonut] = data_all['L1Jet9x9_RawEt'] - data_all['L1Jet9x9_PUEt_ChunkyDonut']

data_all[sL1JetEt_PUS_PhiRing]     = data_all['L1Jet9x9_RawEt'] - (data_all['L1Jet9x9_EtSum7PUTowers'] / 7.0 )


# In[4]:


## data cleaning--------

# Drop entries with L1JetEt < L1JetPtThrsh
data_all_L1EtBelowThrsh = data_all[(
    (data_all[sL1JetEt]  < L1JetPtThrsh) | 
    (data_all[sRefJetEt] < RefJetPtThrsh)
)]
if printLevel >= 8:
    print("data_all[ data_all['{}'] < {} ]: \n{}".format(sL1JetEt, L1JetPtThrsh, data_all_L1EtBelowThrsh))
data_all.drop(index=data_all_L1EtBelowThrsh.index, inplace=True)

print("\nDoes any of the columns have NaN entries: \ndata_all.isna().sum(): \n{}".format(data_all.isna().sum()))
if printLevel >= 5:
    print("\nAfter cleaning, data_all.describe(): \n{}".format(data_all.describe()))


# In[5]:


print("After data cleaning: data_all.columns: {}, \ndata_all.shape: {}".format(data_all.columns, data_all.shape))


# In[6]:


data_all = data_all.sample(frac=fracOfDataToUse, random_state=1)


# In[7]:


print("Sample to use: data_all.columns: {}, \ndata_all.shape: {}".format(data_all.columns, data_all.shape))


# In[8]:


# set trainning and target variables

sL1JetEt_forML  = None
sRefJetEt_forML = None

if   MLTarget == 'GenEt':
    sL1JetEt_forML  = sL1JetEt
    sRefJetEt_forML = sRefJetEt    
    
elif MLTarget == 'logGenEt':    
    sL1JetEt_forML  = 'log%s' % (sL1JetEt)
    sRefJetEt_forML = 'log%s' % (sRefJetEt)
    
    data_all[sL1JetEt_forML] = convert_Et_to_logEt( data_all[sL1JetEt] )
    data_all[sRefJetEt_forML] = convert_Et_to_logEt( data_all[sRefJetEt] )
    
elif MLTarget == 'GenEtByL1Et':    
    sL1JetEt_forML  = sL1JetEt
    sRefJetEt_forML = '%sBy%s' % (sRefJetEt, sL1JetEt)
        
    data_all[sRefJetEt_forML] = convert_GenEt_to_GenEtByL1Et( data_all[sRefJetEt], data_all[sL1JetEt] )   
    
elif MLTarget == 'logGenEtByL1Et':    
    sL1JetEt_forML  = 'log%s' % (sL1JetEt)
    sRefJetEt_forML = 'log%sBy%s' % (sRefJetEt, sL1JetEt)
    
    data_all[sL1JetEt_forML] = convert_Et_to_logEt( data_all[sL1JetEt] )
    data_all[sRefJetEt_forML] = convert_GenEt_to_logGenEtByL1Et( data_all[sRefJetEt],  data_all[sL1JetEt] )
    
    
train_vars = [sL1JetTowerIEtaAbs, sL1JetEt_forML, snVtx]
target_var = sRefJetEt_forML
print("\nsL1JetEt_forML: {}, sRefJetEt_forML: {}".format(sL1JetEt_forML, sRefJetEt_forML))
print("train_vars: {}, \ntarget_var: {}\n".format(train_vars, target_var))

if printLevel >= 1:
    print("data_all.describe(): \n{}".format(data_all.describe()))


# In[9]:


#%%time


sL1JetEt_forML_predict = "%s_predict" % (sL1JetEt_forML)
sL1JetEt_predict       = "%s_predict" % (sL1JetEt)
sSF                    = "SF"

IEta_Cat_forML = OD()
#IEta_Cat_forML['HB'] = [ 1, 16]
#IEta_Cat_forML['HE12a'] = [17, 26]
#IEta_Cat_forML['HE2b'] = [27, 28]
#IEta_Cat_forML['HF30to32'] = [30, 32]
#IEta_Cat_forML['HF33to36'] = [33, 36]
#IEta_Cat_forML['HF37to41'] = [37, 41]
#IEta_Cat['HBEF'] = [ 1, 41]
IEta_Cat_forML['HBEF'] = [ 1, 41]

Pt_Cat_forML = OD()
#Pt_Cat_forML['Ptlt60'] = [ 0, 60]
#Pt_Cat_forML['Ptgt60'] = [60, L1JetPtMax]
#Pt_Cat_forML['Ptlt25']   = [ 0, 25]
#Pt_Cat_forML['Pt25to35'] = [25, 35]
#Pt_Cat_forML['Pt35to60'] = [35, 60]
#Pt_Cat_forML['Pt60to90'] = [60, 90]
#Pt_Cat_forML['Ptgt90']   = [90, L1JetPtMax]
Pt_Cat_forML['PtAll'] = [0, L1JetPtMax]

if printLevel >= 5:
    print("train_vars: {}, \ntarget_var: {}, \nsL1JetEt_forML_predict: {}, \nsL1JetEt_predict: {}, \nsSF: {}".format(
        train_vars, target_var, sL1JetEt_forML_predict, sL1JetEt_predict, sSF))
    
varsOfInterest = train_vars.copy()
#varsOfInterest.extend([target_var, sL1JetEt, sRefJetEt])
for var_ in [target_var, sL1JetEt, sRefJetEt]:
    if var_ not in varsOfInterest:
        varsOfInterest.append(var_)
if printLevel >= 0:
    print("Going for BDT training: varsOfInterest: {}\n".format(varsOfInterest))
    print("After train_vars: {}, \ntarget_var: {}, \nsL1JetEt_forML_predict: {}, \nsL1JetEt_predict: {}, \nsSF: {}".format(
        train_vars, target_var, sL1JetEt_forML_predict, sL1JetEt_predict, sSF))

'''    
# ML training ----------------------------------------------------------------------
def train_MLModel_wHyperopt(X, y):
    X_train, X_valid, y_train, y_valid = train_test_split(X, y, test_size=0.3, random_state=0)
    
    hyperparameter_space = { 
        #'n_estimators': hp.choice('n_estimators', np.arange(500, 2001, 100, dtype=int)),
        'n_estimators': hp.choice('n_estimators', np.arange(700, 701, 100, dtype=int)),
        #'learning_rate':hp.quniform('learning_rate', 0.01, 0.2, 0.01),
        'learning_rate':hp.quniform('learning_rate', 0.05, 0.055, 0.01),
        'early_stopping_rounds': 10
    }    
    max_evals = 1 # 30
    
    def ML_score(params):
        model = XGBRegressor(**params)
        model.fit(
            X_train, y_train,
            eval_set=[(X_valid, y_valid)],
            verbose=False       
        )
        score = mean_squared_error(y_valid, model.predict(X_valid), squared=False)
        if printLevel >= 3:
            print("score: valid {}, train {}. params: {}".format(
                score,
                mean_squared_error(y_train, model.predict(X_train), squared=False),
                params))
        return {'loss': score, 'status': STATUS_OK, 'ML_model': model}
            
    
    def getBestMLModel(trials):
        # https://stackoverflow.com/questions/54273199/how-to-save-the-best-hyperopt-optimized-keras-models-and-its-weights
        valid_trial_list = [trial for trial in trials  if STATUS_OK == trial['result']['status']]
        losses = [float(trial['result']['loss']) for trial in valid_trial_list]
        index_having_minimum_loss = np.argmin(losses)
        best_trial_obj = valid_trial_list[index_having_minimum_loss]
        return best_trial_obj['result']['ML_model']
        
        
        
    # ref: 
    # https://sites.google.com/view/raybellwaves/blog/using-xgboost-and-hyperopt-in-a-kaggle-comp
    # https://www.kaggle.com/code/prashant111/a-guide-on-xgboost-hyperparameters-tuning/notebook
    trials = Trials()
    best_params = fmin(
        fn=ML_score,
        space=hyperparameter_space, 
        algo=tpe.suggest,
        max_evals=max_evals,
        trials=trials)
    print("best_params: {}".format(best_params))
    print("space_eval(hyperparameter_space, best_params): {}".format(space_eval(hyperparameter_space, best_params)))
    
    return getBestMLModel(trials)
# ----------------------------------------------------------------------------        
    
    
    
    
    
BDTModel_dict = OD([])
data_SFs = None
for iEta_category, iEtaBinRange in IEta_Cat_forML.items():
    iEtaBins_i = range(iEtaBinRange[0], iEtaBinRange[-1]+1)

    BDTModel_dict[iEta_category] = OD([])
    for Pt_category, PtRange in Pt_Cat_forML.items():
        PtRangeMin = PtRange[0]
        PtRangeMax = PtRange[1]
        
        data_all_iEtaBins = data_all[
            (data_all[sL1JetTowerIEtaAbs] >= iEtaBinRange[0]) & 
            (data_all[sL1JetTowerIEtaAbs] <= iEtaBinRange[-1]) &
            (data_all[sL1JetEt] >= PtRangeMin) &
            (data_all[sL1JetEt] <  PtRangeMax)
        ][varsOfInterest]
        if printLevel >= 0:
            print("\nbefore training iEta_category {}, iEtaBinRange {}, Pt_category {}, PtRange {}, data_all_iEtaBins.describe(): \n{}".format(
                iEta_category, iEtaBinRange, Pt_category, PtRange, data_all_iEtaBins.describe()))


        X = data_all_iEtaBins[train_vars]
        y = data_all_iEtaBins[target_var]

        xgb_rg = train_MLModel_wHyperopt(X, y)


        BDTModel_dict[iEta_category][Pt_category] = xgb_rg
        
        data_all_iEtaBins[sL1JetEt_forML_predict] = xgb_rg.predict(data_all_iEtaBins[train_vars])
        if   MLTarget == 'GenEt':
            data_all_iEtaBins[sL1JetEt_predict] = data_all_iEtaBins[sL1JetEt_forML_predict]
            
        elif MLTarget == 'logGenEt':
            data_all_iEtaBins[sL1JetEt_predict] = convert_logEt_to_Et( data_all_iEtaBins[sL1JetEt_forML_predict] )

        elif MLTarget == 'GenEtByL1Et':
            data_all_iEtaBins[sL1JetEt_predict] = convert_GenEtByL1Et_to_GenEt( data_all_iEtaBins[sL1JetEt_forML_predict], data_all_iEtaBins[sL1JetEt] )
            
        elif MLTarget == 'logGenEtByL1Et':
            data_all_iEtaBins[sL1JetEt_predict] = convert_logGenEtByL1Et_to_GenEt( data_all_iEtaBins[sL1JetEt_forML_predict], data_all_iEtaBins[sL1JetEt] )
                       
        if printLevel >= 0:
            print("\nafter training iEta_category {}, iEtaBinRange {}, Pt_category {}, PtRange {}, data_all_iEtaBins.describe(): \n{}".format(
                iEta_category, iEtaBinRange, Pt_category, PtRange, data_all_iEtaBins.describe()))
            print(f"\n{data_all_iEtaBins = }")
            
            
            
            
        # save BDT model version
        #sTrain_vars_ = '_'.join()
        #sBDTModel_fileName = '../data/BDTModel_%s_vs_%s__%s_%s.pkl' % ('_'.join(train_vars), target_var, iEta_category, Pt_category)
        sBDTModel_fileName = '../data/BDTModel_%s_%s_%s.pkl' % (version, iEta_category, Pt_category)
        pickle.dump(BDTModel_dict[iEta_category][Pt_category], open(sBDTModel_fileName, "wb"))   
        print(f"\n\nWrote BDT model to {sBDTModel_fileName = }"); sys.stdout.flush()
'''


# In[10]:


def prepareDataframeForSFs(iEtaBinRange, PtRangeMin=L1JetPtThrsh, PtRangeMax=L1JetPtMax, nVtx=48):
    dict_iEta_Et = OD([ (sL1JetTowerIEtaAbs, []), (sL1JetEt, []) ])
    if snVtx in train_vars:
        dict_iEta_Et[snVtx] = []
        
    for iEta in iEtaBinRange:
        list_pt      = np.arange(PtRangeMin, PtRangeMax+1.0)
        list_ietabin = [iEta] * len(list_pt)
        dict_iEta_Et[sL1JetTowerIEtaAbs].extend(list_ietabin) 
        dict_iEta_Et[sL1JetEt].extend(list_pt) 
        if snVtx in train_vars:
            list_nVtx = [nVtx] * len(list_pt)
            dict_iEta_Et[snVtx].extend(list_nVtx) 
            
          
    data_SFs = pd.DataFrame(dict_iEta_Et)
    
    if   MLTarget == 'GenEt':
        data_SFs[sL1JetEt_forML] = data_SFs[sL1JetEt]

    elif MLTarget == 'logGenEt':    
        data_SFs[sL1JetEt_forML] = convert_Et_to_logEt( data_SFs[sL1JetEt] )

    elif MLTarget == 'GenEtByL1Et':    
        data_SFs[sL1JetEt_forML] = data_SFs[sL1JetEt]

    elif MLTarget == 'logGenEtByL1Et':    
        data_SFs[sL1JetEt_forML] = convert_Et_to_logEt( data_SFs[sL1JetEt] )
            
    return data_SFs


# In[11]:


# read BDT models
BDTModel_dict = OD([])
for iEta_category, iEtaBinRange in IEta_Cat_forML.items():
    iEtaBins_i = range(iEtaBinRange[0], iEtaBinRange[-1]+1)
    BDTModel_dict[iEta_category] = OD([])
    
    for Pt_category, PtRange in Pt_Cat_forML.items():
        PtRangeMin = PtRange[0]
        PtRangeMax = PtRange[1]
        
        #xgb_rg = train_MLModel_wHyperopt(X, y)
        #sBDTModel_fileName = '../data/BDTModel_%s_vs_%s__%s_%s.pkl' % ('_'.join(train_vars), target_var, iEta_category, Pt_category)
        sBDTModel_fileName = BDTFileNames_dict[sL1JetEt][MLTarget][iEta_category][Pt_category]
        print(f"{sL1JetEt = }, {MLTarget =}, {iEta_category = }, {Pt_category = }, {sBDTModel_fileName = }")
        xgb_rg = pickle.load(open(sBDTModel_fileName, "rb"))

        BDTModel_dict[iEta_category][Pt_category] = xgb_rg

# compute data_SFs
data_SFs_atPUPoints = OD()

for PUPoint_ForMLEvaluation in PUPointsForMLEvaluation:
    data_SFs = None
    
    for iEta_category, iEtaBinRange in IEta_Cat_forML.items():
        iEtaBins_i = range(iEtaBinRange[0], iEtaBinRange[-1]+1)

        for Pt_category, PtRange in Pt_Cat_forML.items():
            PtRangeMin = PtRange[0]
            PtRangeMax = PtRange[1]

            #xgb_rg = train_MLModel_wHyperopt(X, y)
            #sBDTModel_fileName = '../data/BDTModel_%s_vs_%s__%s_%s.pkl' % ('_'.join(train_vars), target_var, iEta_category, Pt_category)
            sBDTModel_fileName = BDTFileNames_dict[sL1JetEt][MLTarget][iEta_category][Pt_category]
            print(f"{sL1JetEt = }, {MLTarget =}, {iEta_category = }, {Pt_category = }, {sBDTModel_fileName = }")
            xgb_rg = pickle.load(open(sBDTModel_fileName, "rb"))

            BDTModel_dict[iEta_category][Pt_category] = xgb_rg

            xgb_rg = BDTModel_dict[iEta_category][Pt_category]

            if snVtx in train_vars:
                data_SFs_i = prepareDataframeForSFs(iEtaBins_i, PtRangeMin=PtRangeMin, PtRangeMax=PtRangeMax, nVtx=PUPoint_ForMLEvaluation)
            else:
                data_SFs_i = prepareDataframeForSFs(iEtaBins_i, PtRangeMin=PtRangeMin, PtRangeMax=PtRangeMax)

            data_SFs_i[sL1JetEt_forML_predict] = xgb_rg.predict(data_SFs_i[train_vars])        

            if   MLTarget == 'GenEt':
                data_SFs_i[sL1JetEt_predict] = data_SFs_i[sL1JetEt_forML_predict]

            elif MLTarget == 'logGenEt':
                data_SFs_i[sL1JetEt_predict] = convert_logEt_to_Et( data_SFs_i[sL1JetEt_forML_predict] )

            elif MLTarget == 'GenEtByL1Et':
                data_SFs_i[sL1JetEt_predict] = convert_GenEtByL1Et_to_GenEt( data_SFs_i[sL1JetEt_forML_predict], data_SFs_i[sL1JetEt] )

            elif MLTarget == 'logGenEtByL1Et':
                data_SFs_i[sL1JetEt_predict] = convert_logGenEtByL1Et_to_GenEt( data_SFs_i[sL1JetEt_forML_predict], data_SFs_i[sL1JetEt] )

            data_SFs_i[sSF]                    = data_SFs_i[sL1JetEt_predict] / data_SFs_i[sL1JetEt]
            if printLevel >= 11:
                print("iEtaBins_i: {}".format(iEtaBins_i))
                print("data_SFs_i: {}".format(data_SFs_i.describe()))

            if data_SFs is None:
                data_SFs = data_SFs_i
            else:
                data_SFs = pd.concat([data_SFs, data_SFs_i])    
                
    data_SFs_atPUPoints[PUPoint_ForMLEvaluation] = data_SFs

'''    
if snVtx in train_vars:
    sOpFileName_SFs = sOpFileName_SFs.replace('.csv', '_forPU%d.csv' %(nVtx_forSF))
    
#print("Hello1")    
#print("\n\ndata_SFs: \n{}".format(data_SFs.to_string()))
data_SFs.to_csv(sOpFileName_SFs, index=False)
print("Wrote {}".format(sOpFileName_SFs))                
'''


# In[12]:


sL1JetEt_calib = '%s_calib' % (sL1JetEt)
dataVars_forL1JetResponsePlots = [sL1JetTowerIEtaAbs, sL1JetEt, sGenJetEt]
if snVtx in train_vars: 
    dataVars_forL1JetResponsePlots.append(snVtx)


# In[13]:


def calibrateJets_wBDTAtPUPoint(data_all_, data_SFs_):
    data_copy1_     = data_all_[dataVars_forL1JetResponsePlots].copy()
    data_SFs_copy1_ = data_SFs_[[sL1JetTowerIEtaAbs, sL1JetEt, sSF]].copy()
    data_SFs_copy1_ = data_SFs_copy1_.set_index([sL1JetTowerIEtaAbs, sL1JetEt])
    SFs_dict_       = data_SFs_copy1_.to_dict()[sSF]

    def calibrateJet(Et_0, iEta):
        Et = round(Et_0)
        if Et < calibSF_L1JetPtRange[0]: Et = round(calibSF_L1JetPtRange[0])
        if Et > calibSF_L1JetPtRange[1]: Et = round(calibSF_L1JetPtRange[1])
        #print(f"iEta {iEta}, Et {Et}")
        sf = SFs_dict_[(iEta, Et)] if Et >= 1 else SF_forZeroPt
        return Et_0 * sf

    data_copy1_[sL1JetEt_calib] = data_copy1_.apply(lambda row: calibrateJet(row[sL1JetEt], row[sL1JetTowerIEtaAbs]), axis=1)
    #data_copy1[sL1JetEt_calib] = np.vectorize(calibrateJet)(data_copy1[sL1JetEt], data_copy1[sL1JetTowerIEtaAbs])
    if printLevel >= 5:
        print("calibrateJets_wBDTAtPUPoint():: data_copy1_: {}".format(data_copy1_))
        
    return data_copy1_[sL1JetEt_calib]


# In[14]:


data_copy1     = data_all[dataVars_forL1JetResponsePlots].copy()
for PUPoint_ForMLEvaluation in PUPointsForMLEvaluation:
    data_copy1['%s_wBDTAtPU%d' % (sL1JetEt_calib, PUPoint_ForMLEvaluation)] = calibrateJets_wBDTAtPUPoint(
        data_all, data_SFs_atPUPoints[PUPoint_ForMLEvaluation]
    )
if printLevel >= 5:
    print("data_copy1: {}".format(data_copy1))    


# In[15]:


# L1JetResponse vs Eta after JEC

sOutDir_toUse   = '%s/L1JetResponse_vs_Eta_perPU_perPtEtaCat' % (sOutDirAfterJEC)
sOutDir1D_toUse = '%s/L1JetResponse_vs_Eta_perPU_perPtEtaCat/1D' % (sOutDirAfterJEC)
if not os.path.exists(sOutDir_toUse):             os.makedirs( sOutDir_toUse , exist_ok=True )
if not os.path.exists(sOutDir1D_toUse):           os.makedirs( sOutDir1D_toUse , exist_ok=True )    

Pt_Cat_forResolutionPlots = OD()
Pt_Cat_forResolutionPlots['Pt25to35']   = [ 25,  35]
Pt_Cat_forResolutionPlots['Pt40to55']   = [ 40,  55]
Pt_Cat_forResolutionPlots['Pt60to90']   = [ 60,  90]
Pt_Cat_forResolutionPlots['Pt100to200'] = [100, 200]
#Pt_Cat_forResolutionPlots['PtAll'] = [0, L1JetPtMax]

PU_Cat_forResolutionPlots = OD()
PU_Cat_forResolutionPlots['PUlt38'] = [ 0, 38]
PU_Cat_forResolutionPlots['PUgt58'] = [58, 99]
#PU_Cat_forResolutionPlots['PUAll'] = [0, 99]

sL1JetResponse = 'L1JetEt/%s' % (sRefJetEt)

marker_color_list = ['r', 'b', 'darkviolet', 'c', 'orange', 'green']
marker_style_list = ["o", "X", '>', '^', 'v', "s", "+", 'x', '*']
marker_size_list  = [5, 5, 3, 3, 3, 3, 2, 2, 2]


for PU_category, PURange in PU_Cat_forResolutionPlots.items():
    PURangeMin = PURange[0]
    PURangeMax = PURange[1]

    for Pt_category, PtRange in Pt_Cat_forResolutionPlots.items():
        PtRangeMin = PtRange[0]
        PtRangeMax = PtRange[1]
        
        JES_atPUPointsForMLEvaluations = OD()
        JER_atPUPointsForMLEvaluations = OD()        
        for PUPoint_ForMLEvaluation in [-1] + PUPointsForMLEvaluation :
            if PUPoint_ForMLEvaluation == -1:
                sL1JetEt_forResolutionPlots = sL1JetEt
            else:
                sL1JetEt_forResolutionPlots = '%s_wBDTAtPU%d' % (sL1JetEt_calib, PUPoint_ForMLEvaluation)
            
            JES = OD()
            JER = OD()
            for iEtaBin in iEtaBins:
                data_toUse_ = data_copy1[
                    (data_copy1[sL1JetTowerIEtaAbs] == iEtaBin   ) & 
                    (data_copy1[sRefJetEt]               >= PtRangeMin) &
                    (data_copy1[sRefJetEt]               <  PtRangeMax) & 
                    (data_copy1[snVtx]                   >= PURangeMin) &
                    (data_copy1[snVtx]                   <  PURangeMax)
                ]#.copy()
                data_toUse_[sL1JetResponse] = data_toUse_[sL1JetEt_forResolutionPlots] / data_toUse_[sRefJetEt]
                print(f"PU: {PURange}, Pt: {PtRange}, {PUPoint_ForMLEvaluation = }, iEta: {iEtaBin}")

                h = hist.Hist.new.Regular(100,0,2.5, name=sL1JetResponse).Weight()            
                h.fill(data_toUse_[sL1JetResponse])

                x_    = h.axes[0].centers
                y_    = h.values()
                yErr_ = np.sqrt(h.variances())

                # index of bins with y>0
                idx_NonZeroY = np.nonzero(y_)

                # Give initial parameters for Gaussian function fit
                #pInitial = [y_.max(), 1, 0.3]
                pInitial = [y_.max(), np.mean(data_toUse_.loc[data_toUse_[sL1JetResponse] < 2.5][sL1JetResponse]), np.std(data_toUse_.loc[data_toUse_[sL1JetResponse] < 2.5][sL1JetResponse])]

                #popt, pcov = curve_fit(GaussianFunction, xdata=x_, ydata=y_, p0=pInitial)
                try:
                    popt, pcov = curve_fit(GaussianFunction, xdata=x_[idx_NonZeroY], ydata=y_[idx_NonZeroY], p0=pInitial, sigma=yErr_[idx_NonZeroY])
                except:
                    print(f"{PU_category = }, {Pt_category = }, {PUPoint_ForMLEvaluation = }, {iEtaBin = }, {pInitial = }:  fit did not converge *** ")
                    continue
                poptErr    = np.sqrt(np.diag(pcov))
                print(f"popt: {popt}, \npcov: \n{pcov},  \n{poptErr = }")

                Mean_  = popt[1];    errMean_  = poptErr[1]
                Sigma_ = popt[2];    errSigma_ = poptErr[2]

                x1Dense_ = np.arange(x_.min(), x_.max(), 0.01)            

                fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5,4), layout='constrained')
                axs.errorbar(x_, y_, yerr=yErr_, label='Entries', fmt='o', markersize=1)
                axs.plot(x1Dense_, GaussianFunction(x1Dense_, *popt), label=r'Gaussian fit, $\mu:\,  %.2f \pm %.2f$, $\sigma:\, %.2f \pm %.2f$'%(Mean_,errMean_, Sigma_,errSigma_))
                axs.set_xlabel(sL1JetResponse)
                axs.set_ylabel('Entries')
                axs.set_title('%s %s iEta %s' % (PU_category, Pt_category, iEtaBin))            
                axs.legend(bbox_to_anchor=(0.0, 1.0), loc='upper left', borderaxespad=0.9)
                axs.set_ylim(0, y_.max() * 1.4)
                axs.grid()

                errJER_      = calculate_errorOfRatio(N=Sigma_, D=Mean_, eN=errSigma_, eD=errMean_)            
                JES[iEtaBin] = {'value': Mean_,          'error': errMean_}
                JER[iEtaBin] = {'value': Sigma_ / Mean_, 'error': errJER_}

                sOutDir1D_toUse_ = '%s/MLAtPU%d' % (sOutDir1D_toUse, PUPoint_ForMLEvaluation)
                if not os.path.exists(sOutDir1D_toUse_):           os.makedirs( sOutDir1D_toUse_ , exist_ok=True )    
                fig.savefig('%s/L1JetResponse_1D_%s_%s_ieta_%d.png' % (sOutDir1D_toUse_, PU_category, Pt_category, iEtaBin))
                plt.close(fig)
                
            JES_atPUPointsForMLEvaluations[PUPoint_ForMLEvaluation] = JES
            JER_atPUPointsForMLEvaluations[PUPoint_ForMLEvaluation] = JER
                
            

        # plot JES vs iEta
        fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5,4), layout='constrained')
        yMin_ =  9999
        yMax_ = -9999
        for iPU, PUPoint_ForMLEvaluation in enumerate([-1] + PUPointsForMLEvaluation):
            JES = JES_atPUPointsForMLEvaluations[PUPoint_ForMLEvaluation]
            
            JES_iEtawise    = [ JES[iEtaBin]['value'] for iEtaBin in JES.keys()]
            errJES_iEtawise = [ JES[iEtaBin]['error'] for iEtaBin in JES.keys()]

            yMin1_ = np.min(np.array(JES_iEtawise) - np.array(errJES_iEtawise))
            yMax1_ = np.max(np.array(JES_iEtawise) + np.array(errJES_iEtawise))
            if yMin1_ < yMin_: yMin_ = yMin1_
            if yMax1_ > yMax_: yMax_ = yMax1_

            label_ = r'JEC SFs @PU%d' % (PUPoint_ForMLEvaluation) if PUPoint_ForMLEvaluation != -1 else 'Before JEC'
            #axs.errorbar(list(JES.keys()), JES_iEtawise, yerr=errJES_iEtawise, fmt='o', markersize=2, label=label_)
            axs.errorbar(list(JES.keys()), JES_iEtawise, yerr=errJES_iEtawise, fmt='', linestyle='', markersize=2, marker=marker_style_list[iPU], color=marker_color_list[iPU],  label=label_)
            
        axs.set_xlabel('iEta')
        axs.set_ylabel(r'$\mu$ (%s)' % (sL1JetResponse))
        axs.set_title('%s %s' % (PU_category, Pt_category))
        if yMin_ < 0.4 or yMax_ > 1.5:
            axs.set_ylim(0.4, 1.5)
        else:
            axs.set_ylim(yMin_, yMax_)
        #axs.axhline(y=1, linestyle='--')
        axs.legend(bbox_to_anchor=(0.1, 1.05), loc='upper left', borderaxespad=0.9, ncol=2)
        axs.margins(y=0.3)
        axs.grid()
        fig.savefig('%s/L1JetResponse_vs_iEta_%s_%s_Mean.png' % (sOutDir_toUse, PU_category, Pt_category))
        plt.close(fig)


        # plot JER vs iEta
        fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5,4), layout='constrained')
        yMin_ =  9999
        yMax_ = -9999
        for iPU, PUPoint_ForMLEvaluation in enumerate([-1] + PUPointsForMLEvaluation):
            JER = JER_atPUPointsForMLEvaluations[PUPoint_ForMLEvaluation]            
            
            JER_iEtawise    = [ JER[iEtaBin]['value'] for iEtaBin in JER.keys()]
            errJER_iEtawise = [ JER[iEtaBin]['error'] for iEtaBin in JER.keys()]

            yMin1_ = np.min(np.array(JER_iEtawise) - np.array(errJER_iEtawise))
            yMax1_ = np.max(np.array(JER_iEtawise) + np.array(errJER_iEtawise))
            if yMin1_ < yMin_: yMin_ = yMin1_
            if yMax1_ > yMax_: yMax_ = yMax1_

            label_ = r'JEC SFs @PU%d' % (PUPoint_ForMLEvaluation) if PUPoint_ForMLEvaluation != -1 else 'Before JEC'
            #axs.errorbar(list(JER.keys()), JER_iEtawise, yerr=errJER_iEtawise, fmt='o', markersize=2, label=label_)
            axs.errorbar(list(JER.keys()), JER_iEtawise, yerr=errJER_iEtawise, fmt='', linestyle='', markersize=2, marker=marker_style_list[iPU], color=marker_color_list[iPU],  label=label_)
            
        axs.set_xlabel('iEta')
        axs.set_ylabel(r'$\sigma/\mu$ (%s)' % (sL1JetResponse))
        axs.set_title('%s %s' % (PU_category, Pt_category))
        if yMin_ < 0. or yMax_ > 1.:
            axs.set_ylim(0., 1.)
        else:
            axs.set_ylim(yMin_, yMax_)                
        axs.legend(bbox_to_anchor=(0.1, 1.05), loc='upper left', borderaxespad=0.9, ncol=2)
        axs.margins(y=0.3)
        axs.grid()
        fig.savefig('%s/L1JetResponse_vs_iEta_%s_%s_Resolution.png' % (sOutDir_toUse, PU_category, Pt_category)) 
        plt.close(fig)


# In[18]:


# L1JetResponse vs Pt after JEC

exit(0)

sOutDir_toUse   = '%s/L1JetResponse_vs_Pt_perPU_perEtaCat_BDTAtPU%d' % (sOutDirAfterJEC, PUPointForMLEvaluation_selected)
if not os.path.exists(sOutDir_toUse):             os.makedirs( sOutDir_toUse , exist_ok=True )   

sL1JetEt_forResolutionPlots = '%s_wBDTAtPU%d' % (sL1JetEt_calib, PUPointForMLEvaluation_selected)

IETA_CAT_forResolutionPlots = OD()
#IETA_CAT_forResolutionPlots['HBEF'] = [ 1, 41]  ## Whole detector, 1 - 41
IETA_CAT_forResolutionPlots['iEta1to16']   = [ 1, 16]  ## Trigger towers  1 - 16
IETA_CAT_forResolutionPlots['iEta17to23']  = [17, 23]  ## Trigger towers 17 - 20
IETA_CAT_forResolutionPlots['iEta24to28']  = [24, 28]  ## Trigger towers 21 - 25
IETA_CAT_forResolutionPlots['iEta30to34']  = [30, 34]  ## Trigger towers 30 - 41
IETA_CAT_forResolutionPlots['iEta34to41']  = [35, 41]  ## Trigger towers 30 - 41

PU_Cat_forResolutionPlots = OD()
PU_Cat_forResolutionPlots['PUlt38']   = [ 0, 38]
PU_Cat_forResolutionPlots['PU39to57'] = [39, 57]
PU_Cat_forResolutionPlots['PUgt58']   = [58, 99]
#PU_Cat_forResolutionPlots['PUAll'] = [0, 99]

RefJetPtBins_forResolutionPlots = [
    *np.arange( 24.5,  60.5, 3),
    *np.arange( 60.5, 120.5, 5),
    *np.arange(120.5, 150.5, 10),
    *np.arange(150.5, 230.5, 20),
    *np.arange(230.5, 255.5, 25),
]
print(f"RefJetPtBins_forResolutionPlots: {RefJetPtBins_forResolutionPlots}")





for PU_category, PURange in PU_Cat_forResolutionPlots.items():
    PURangeMin = PURange[0]
    PURangeMax = PURange[1]

    JES_iEtaCatwise = OD()
    JER_iEtaCatwise = OD()
    for iEta_cat, iEtaBinRange in IETA_CAT_forResolutionPlots.items():
        iEtaBinRangeMin = iEtaBinRange[0]
        iEtaBinRangeMax = iEtaBinRange[1]                    

        JES = OD()
        JER = OD()
        for iPtBin in range(len(RefJetPtBins_forResolutionPlots) - 1):
            PtRangeMin = RefJetPtBins_forResolutionPlots[iPtBin]
            PtRangeMax = RefJetPtBins_forResolutionPlots[iPtBin + 1]
            PtMean = (PtRangeMin + PtRangeMax)/2

            data_toUse_ = data_copy1[
                (data_copy1[sL1JetTowerIEtaAbs]      >= iEtaBinRangeMin   ) & 
                (data_copy1[sL1JetTowerIEtaAbs]      <= iEtaBinRangeMax   ) & 
                (data_copy1[sRefJetEt]               >= PtRangeMin) &
                (data_copy1[sRefJetEt]               <  PtRangeMax) & 
                (data_copy1[snVtx]                   >= PURangeMin) &
                (data_copy1[snVtx]                   <  PURangeMax)
            ]#.copy()
            data_toUse_[sL1JetResponse] = data_toUse_[sL1JetEt_forResolutionPlots] / data_toUse_[sRefJetEt]

            h = hist.Hist.new.Regular(100,0,2.5, name=sL1JetResponse).Weight()            
            h.fill(data_toUse_[sL1JetResponse])

            x_    = h.axes[0].centers
            y_    = h.values()
            yErr_ = np.sqrt(h.variances())

            # index of bins with y>0
            idx_NonZeroY = np.nonzero(y_)

            # Give initial parameters for Gaussian function fit
            #pInitial = [y_.max(), 1, 0.3]
            pInitial = [y_.max(), np.mean(data_toUse_.loc[data_toUse_[sL1JetResponse] < 2.5][sL1JetResponse]), np.std(data_toUse_.loc[data_toUse_[sL1JetResponse] < 2.5][sL1JetResponse])]

            #popt, pcov = curve_fit(GaussianFunction, xdata=x_, ydata=y_, p0=pInitial)
            try:
                popt, pcov = curve_fit(GaussianFunction, xdata=x_[idx_NonZeroY], ydata=y_[idx_NonZeroY], p0=pInitial, sigma=yErr_[idx_NonZeroY])
            except:
                print(f"{PU_category = }, {iEta_cat = }, {PtRangeMin = }, {PtRangeMax = }, {pInitial = }:  fit did not converge *** ")
                continue
            poptErr    = np.sqrt(np.diag(pcov))
            print(f"popt: {popt}, \npcov: \n{pcov},  \n{poptErr = }")

            Mean_  = popt[1];    errMean_  = poptErr[1]
            Sigma_ = popt[2];    errSigma_ = poptErr[2]

            x1Dense_ = np.arange(x_.min(), x_.max(), 0.01)            

            fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5,4), layout='constrained')
            axs.errorbar(x_, y_, yerr=yErr_, label='Entries', fmt='o', markersize=1)
            axs.plot(x1Dense_, GaussianFunction(x1Dense_, *popt), label=r'Gaussian fit, $\mu:\,  %.2f \pm %.2f$, $\sigma:\, %.2f \pm %.2f$'%(Mean_,errMean_, Sigma_,errSigma_))
            axs.set_xlabel(sL1JetResponse)
            axs.set_ylabel('Entries')
            axs.set_title('%s %s pT [%.1f, %.1f]' % (PU_category, iEta_cat, PtRangeMin, PtRangeMax))
            axs.legend(bbox_to_anchor=(0.0, 1.0), loc='upper left', borderaxespad=0.9)
            axs.set_ylim(0, y_.max() * 1.4)
            axs.grid()

            errJER_      = calculate_errorOfRatio(N=Sigma_, D=Mean_, eN=errSigma_, eD=errMean_)            
            JES[PtMean] = {'value': Mean_,          'error': errMean_}
            JER[PtMean] = {'value': Sigma_ / Mean_, 'error': errJER_}

            sOutDir1D_toUse_ = '%s/%s/1D' % (sOutDir_toUse,PU_category)
            if not os.path.exists(sOutDir1D_toUse_):           os.makedirs( sOutDir1D_toUse_, exist_ok=True )    
            fig.savefig('%s/L1JetResponse_1D_%s_%s_Pt_%.1f.png' % (sOutDir1D_toUse_, PU_category, iEta_cat, PtMean))
            plt.close(fig)

        JES_iEtaCatwise[iEta_cat] = JES
        JER_iEtaCatwise[iEta_cat] = JER
        
        
            
    # plot JES vs Pt
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5,4), layout='constrained')
    yMin_ =  9999
    yMax_ = -9999
    i_ = 0
    for iEta_cat, iEtaBinRange in IETA_CAT_forResolutionPlots.items():
        JES = JES_iEtaCatwise[iEta_cat]
        
        JES_Ptwise    = [ JES[Pt]['value'] for Pt in JES.keys()]
        errJES_Ptwise = [ JES[Pt]['error'] for Pt in JES.keys()]

        yMin1_ = np.min(np.array(JES_Ptwise) - np.array(errJES_Ptwise))
        yMax1_ = np.max(np.array(JES_Ptwise) + np.array(errJES_Ptwise))
        if yMin1_ < yMin_: yMin_ = yMin1_
        if yMax1_ > yMax_: yMax_ = yMax1_
        
        axs.errorbar(list(JES.keys()), JES_Ptwise, yerr=errJES_Ptwise, fmt='', linestyle='', markersize=2, marker=marker_style_list[i_], color=marker_color_list[i_],  label=iEta_cat)
        i_ += 1
        
    axs.set_xlabel('%s [GeV]' % (sRefJetEt))
    axs.set_ylabel(r'$\mu$ (%s)' % (sL1JetResponse))
    axs.set_title('%s ' % (PU_category))
    if yMin_ < 0.4 or yMax_ > 1.4:
        axs.set_ylim(0.4, 1.4)
    else:
        axs.set_ylim(yMin_, yMax_)
    axs.legend(bbox_to_anchor=(0.1, 1.05), loc='upper left', borderaxespad=0.9, ncol=1)
    axs.margins(y=0.3)
    axs.grid()
    fig.savefig('%s/L1JetResponse_vs_Pt_%s_Mean.png' % (sOutDir_toUse, PU_category))   
    plt.close(fig)
        
        
        

    # plot JER vs iEta
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5,4), layout='constrained')
    yMin_ =  9999
    yMax_ = -9999
    i_ = 0
    for iEta_cat, iEtaBinRange in IETA_CAT_forResolutionPlots.items():
        JER = JER_iEtaCatwise[iEta_cat]
    
        JER_Ptwise    = [ JER[Pt]['value'] for Pt in JER.keys()]
        errJER_Ptwise = [ JER[Pt]['error'] for Pt in JER.keys()]

        yMin1_ = np.min(np.array(JER_Ptwise) - np.array(errJER_Ptwise))
        yMax1_ = np.max(np.array(JER_Ptwise) + np.array(errJER_Ptwise))        
        if yMin1_ < yMin_: yMin_ = yMin1_
        if yMax1_ > yMax_: yMax_ = yMax1_

        axs.errorbar(list(JER.keys()), JER_Ptwise, yerr=errJER_Ptwise, fmt='', linestyle='', markersize=2, marker=marker_style_list[i_], color=marker_color_list[i_],  label=iEta_cat)
        i_ += 1
        
    axs.set_xlabel('%s [GeV]' % (sRefJetEt))
    axs.set_ylabel(r'$\sigma/\mu$ (%s)' % (sL1JetResponse))
    axs.set_title('%s' % (PU_category))
    if yMin_ < 0. or yMax_ > 1.:
        axs.set_ylim(0., 1.)
    else:
        axs.set_ylim(yMin_, yMax_)        
    axs.legend(bbox_to_anchor=(0.1, 1.05), loc='upper left', borderaxespad=0.9, ncol=1)
    axs.margins(y=0.3)
    axs.grid()
    fig.savefig('%s/L1JetResponse_vs_Pt_%s_Resolution.png' % (sOutDir_toUse, PU_category))
    plt.close(fig)

