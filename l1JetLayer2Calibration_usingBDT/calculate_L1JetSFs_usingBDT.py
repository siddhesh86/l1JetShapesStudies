#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import argparse
from collections import OrderedDict as OD

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors

from xgboost import XGBRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error


parser = argparse.ArgumentParser()
parseGroup1 = parser.add_mutually_exclusive_group(required=True)
parseGroup1.add_argument('--ChunkyDonut',    action='store_true', default=True)
parseGroup1.add_argument('--PhiRing',        action='store_true', default=False)
parseGroup2 = parser.add_mutually_exclusive_group(required=True)
parseGroup2.add_argument('--l1MatchOffline', action='store_true', default=False)
parseGroup2.add_argument('--l1MatchGen',     action='store_true', default=True)


#args = parser.parse_args()
args = parser.parse_args("--ChunkyDonut --l1MatchGen".split()) # to run in jupyter-notebook 
l1Jet_ChunkyDonut = args.ChunkyDonut
l1Jet_PhiRing     = args.PhiRing
l1MatchOffline    = args.l1MatchOffline
l1MatchGen        = args.l1MatchGen


sIpFileName     = "../data/L1T_Jet_MLInputs_Run3_QCD_Pt15to7000_PFA1p_CMSSW12_6_0_pre1_nVtxAll_20220925.csv"
sOpFileName_SFs = "../data/L1T_Jet_SFs_Run3_QCD_Pt15to7000_PFA1p_CMSSW12_6_0_pre1_nVtxAll_20220925.csv"
sOutDir         = "./plots"

iEtaBins = [i for i in range(1, 42) if i!=29]
sL1JetEt_PUS_ChunkyDonut = 'L1JetEt_PUS_ChunkyDonut'
sL1JetEt_PUS_PhiRing     = 'L1JetEt_PUS_PhiRing'
sOfflineJetEt            = 'PFJetEtCorr'
sGenJetEt                = 'GenJetEt'
sL1JetTowerIEtaAbs       = 'L1JetTowerIEtaAbs'
L1JetPtThrsh             = 10.0 # GeV
L1JetPtMax               = 255.0 # GeV

sL1JetEt  = sL1JetEt_PUS_ChunkyDonut if l1Jet_ChunkyDonut else sL1JetEt_PUS_PhiRing
sRefJetEt = sOfflineJetEt if l1MatchOffline else sGenJetEt 


data_all = pd.read_csv(sIpFileName)
print("Input file: %s" % (sIpFileName))
print("iEtaBins ({}): {}".format(len(iEtaBins), iEtaBins))
print("sRefJetEt: {}, \t sL1Jet: {}, \t L1JetPtThrsh: {}".format(sRefJetEt, sL1JetEt, L1JetPtThrsh))


# In[2]:


print("data_all.columns: {}, \ndata_all.shape: {}".format(data_all.columns, data_all.shape))



# In[3]:


data_all[sL1JetEt_PUS_ChunkyDonut] = data_all['L1Jet9x9_RawEt'] - data_all['L1Jet9x9_PUEt_ChunkyDonut']

data_all[sL1JetEt_PUS_PhiRing]     = data_all['L1Jet9x9_RawEt'] - (data_all['L1Jet9x9_EtSum7PUTowers'] / 7.0 )

sL1JetEt_forML              = 'log_%s' % (sL1JetEt)
sRefJetEt_forML             = 'log_%s' % (sRefJetEt)
train_vars = [sL1JetTowerIEtaAbs, sL1JetEt_forML]
target_var = sRefJetEt_forML
print("\nsL1JetEt_forML: {}, sRefJetEt_forML: {}".format(sL1JetEt_forML, sRefJetEt_forML))
print("train_vars: {}, \ntarget_var: {}\n".format(train_vars, target_var))

def transform_JetEt_forML(series):
    series_new = np.log(series)
    return series_new;

def transform_back_JetEt_fromML(series):
    series_new = np.exp(series)
    return series_new;

def prepareDataframeForSFs(iEtaBinRange):
    dict_iEta_Et = OD([ (sL1JetTowerIEtaAbs, []), (sL1JetEt, []) ])
    for iEta in iEtaBinRange:
        list_pt      = np.arange(L1JetPtThrsh+1.0, L1JetPtMax+1.0)
        list_ietabin = [iEta] * len(list_pt)
        dict_iEta_Et[sL1JetTowerIEtaAbs].extend(list_ietabin) 
        dict_iEta_Et[sL1JetEt].extend(list_pt) 
          
    data_SFs = pd.DataFrame(dict_iEta_Et)
    data_SFs[sL1JetEt_forML]    = transform_JetEt_forML(data_SFs[sL1JetEt])
    return data_SFs

data_all[sL1JetEt_forML]    = transform_JetEt_forML(data_all[sL1JetEt])
data_all[sRefJetEt_forML]   = transform_JetEt_forML(data_all[sRefJetEt])

print("data_all.describe(): \n{}".format(data_all.describe()))
#print("\n\ndata_SFs.describe(): \n{}".format(data_SFs.describe()))
#print("\n\ndata_SFs.describe(): \n{}".format(data_SFs.to_string()))


# In[4]:


## data cleaning--------

# Drop entries with L1JetEt < L1JetPtThrsh
data_all_L1EtBelowThrsh = data_all[ data_all[sL1JetEt] < L1JetPtThrsh ]
print("data_all[ data_all['{}'] < {} ]: \n{}".format(sL1JetEt, L1JetPtThrsh, data_all_L1EtBelowThrsh))
data_all.drop(index=data_all_L1EtBelowThrsh.index, inplace=True)

print("\nDoes any of the columns have NaN entries: \ndata_all.isna().sum(): \n{}".format(data_all.isna().sum()))
print("\nAfter cleaning, data_all.describe(): \n{}".format(data_all.describe()))


# In[5]:


## make resolution plots before JEC

if not os.path.exists(sOutDir): 
    os.makedirs(sOutDir)

#print("".format())
for iEtaBinRange in np.array_split(iEtaBins, 8):
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5,4), layout='constrained')
        
    
    for iEtaBin in iEtaBinRange:
        data_all_iEtaBin = data_all[data_all[sL1JetTowerIEtaAbs] == iEtaBin]
        axs.hist(
            (data_all_iEtaBin[sL1JetEt]/data_all_iEtaBin[sRefJetEt]), 
            bins=100, range=(0, 2.6),
            label='iEta %d' % (iEtaBin),
            histtype='step',#, linewidth=2
            density=True
        )
    axs.set_xlabel('L1JetEt / %s' % (sRefJetEt))
    axs.set_ylabel('Normalized entries')
    axs.set_title('%s' % (sL1JetEt))
    axs.legend()
        
    fig.savefig('%s/beforeJEC_%s_ieta_%d_to_%d.png' % (sOutDir, sL1JetEt, iEtaBinRange[0], iEtaBinRange[-1]))
        


# In[6]:


nEntriesPerIEtaBin = [len(data_all[data_all[sL1JetTowerIEtaAbs] == iEtaBin].index) for iEtaBin in iEtaBins]   
nEntriesPerIEtaBin_1 = { iEtaBin: len(data_all[data_all[sL1JetTowerIEtaAbs] == iEtaBin].index) for iEtaBin in iEtaBins}   
fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5,4), layout='constrained')
axs.plot(iEtaBins, nEntriesPerIEtaBin)
axs.set_xlabel('iEta bins')
axs.set_ylabel('Entries')
axs.set_title('%s' % (sL1JetEt))
#axs.legend()
fig.savefig('%s/%s_nEntriesPerIEtaBin.png' % (sOutDir, sL1JetEt))
print("nEntriesPerIEtaBin: {}".format(nEntriesPerIEtaBin))
print("nEntriesPerIEtaBin_1: {}".format(nEntriesPerIEtaBin_1))


# In[ ]:


get_ipython().run_cell_magic('time', '', '\nIEta_Cat = OD()\n#IEta_Cat[\'HB\'] = [ 1, 16]\n#IEta_Cat[\'HE\'] = [17, 28]\n#IEta_Cat[\'HF\'] = [30, 41]\nIEta_Cat[\'HBEF\'] = [ 1, 41]\n\nvarsOfInterest = train_vars\nvarsOfInterest.extend([target_var, sL1JetEt, sRefJetEt])\nprint("varsOfInterest: {}\\n".format(varsOfInterest))\n\nfor iEta_category, iEtaBinRange in IEta_Cat.items():\n    data_all_iEtaBins = data_all[\n        (data_all[sL1JetTowerIEtaAbs] >= iEtaBinRange[0]) & \n        (data_all[sL1JetTowerIEtaAbs] <= iEtaBinRange[-1])\n    ][varsOfInterest]\n    print("\\niEta_category {}, iEtaBinRange {}, data_all_iEtaBins.describe(): \\n{}".format(iEta_category, iEtaBinRange, data_all_iEtaBins.describe()))\n    \n    X = data_all_iEtaBins[train_vars]\n    y = data_all_iEtaBins[target_var]\n\n    X_train, X_valid, y_train, y_valid = train_test_split(X, y, test_size=0.3, random_state=0)\n    \n    # declare parameters\n    params_i = {\n        \'n_estimators\': 1000, \n        \'learning_rate\': 0.01, \n        \'early_stopping_rounds\': 5\n    }\n    xgb_rg = XGBRegressor(**params_i)\n    xgb_rg.fit(X_train, y_train,\n               eval_set=[(X_valid, y_valid)],\n               verbose=False\n             )\n    print("\\niEta_category {}, iEtaBinRange {}, params: {}, mean_squared_error: {}; {}".format(\n        iEta_category, iEtaBinRange, \n        params_i,\n        mean_squared_error(y_valid, xgb_rg.predict(X_valid)),\n        mean_squared_error(y_valid, xgb_rg.predict(X_valid), squared=False)\n    ))\n    \n    \n    \n    \nprint("Hello1")    \n')


# In[ ]:




