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


args = parser.parse_args()
#args = parser.parse_args("--ChunkyDonut --l1MatchGen".split()) # to run in jupyter-notebook 
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

sOpFileName_SFs = sOpFileName_SFs.replace('.csv', '_%s.csv' % (sL1JetEt))
sOutDir = '%s_%s' % (sOutDir, sL1JetEt)


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


# In[7]:


#%%time


sL1JetEt_forML_predict = "%s_predict" % (sL1JetEt_forML)
sL1JetEt_predict       = "%s_predict" % (sL1JetEt)
sSF                    = "SF"

IEta_Cat = OD()
IEta_Cat['HB'] = [ 1, 16]
IEta_Cat['HE'] = [17, 28]
IEta_Cat['HF'] = [30, 41]
#IEta_Cat['HBEF'] = [ 1, 41]

print("train_vars: {}, \ntarget_var: {}, \nsL1JetEt_forML_predict: {}, \nsL1JetEt_predict: {}, \nsSF: {}".format(
        train_vars, target_var, sL1JetEt_forML_predict, sL1JetEt_predict, sSF
))
varsOfInterest = train_vars.copy()
varsOfInterest.extend([target_var, sL1JetEt, sRefJetEt])
print("varsOfInterest: {}\n".format(varsOfInterest))

print("After train_vars: {}, \ntarget_var: {}, \nsL1JetEt_forML_predict: {}, \nsL1JetEt_predict: {}, \nsSF: {}".format(
        train_vars, target_var, sL1JetEt_forML_predict, sL1JetEt_predict, sSF
))

dafa_SFs = None
for iEta_category, iEtaBinRange in IEta_Cat.items():
    iEtaBins_i = range(iEtaBinRange[0], iEtaBinRange[-1]+1)
    data_all_iEtaBins = data_all[
        (data_all[sL1JetTowerIEtaAbs] >= iEtaBinRange[0]) & 
        (data_all[sL1JetTowerIEtaAbs] <= iEtaBinRange[-1])
    ][varsOfInterest]
    print("\niEta_category {}, iEtaBinRange {}, data_all_iEtaBins.describe(): \n{}".format(
        iEta_category, iEtaBinRange, data_all_iEtaBins.describe()
    ))
    print("train_vars: {}, target_var: {}, sL1JetEt_forML_predict: {}, sL1JetEt_predict: {}, sSF: {}".format(
        train_vars, target_var, sL1JetEt_forML_predict, sL1JetEt_predict, sSF
    ))
    
        
    X = data_all_iEtaBins[train_vars]
    y = data_all_iEtaBins[target_var]

    X_train, X_valid, y_train, y_valid = train_test_split(X, y, test_size=0.3, random_state=0)
    
    # declare parameters
    params_i = {
        'n_estimators': 1000, 
        'learning_rate': 0.01, 
        'early_stopping_rounds': 5
    }
    xgb_rg = XGBRegressor(**params_i)
    xgb_rg.fit(X_train, y_train,
               eval_set=[(X_valid, y_valid)],
               verbose=False
    )
    print("\niEta_category {}, iEtaBinRange {}, params: {}, mean_squared_error: valid {}; train {}".format(
        iEta_category, iEtaBinRange, 
        params_i,
        #mean_squared_error(y_valid, xgb_rg.predict(X_valid)),
        mean_squared_error(y_valid, xgb_rg.predict(X_valid), squared=False),
        mean_squared_error(y_train, xgb_rg.predict(X_train), squared=False)
        
    ))
    
    print("{} here1".format(iEta_category))
    dafa_SFs_i = prepareDataframeForSFs(iEtaBins_i)
    print("{} here2".format(iEta_category))
    dafa_SFs_i[sL1JetEt_forML_predict] = xgb_rg.predict(dafa_SFs_i[train_vars])
    print("{} here3".format(iEta_category))
    dafa_SFs_i[sL1JetEt_predict]       = transform_back_JetEt_fromML( dafa_SFs_i[sL1JetEt_forML_predict] )
    print("{} here4".format(iEta_category))
    dafa_SFs_i[sSF]                    = dafa_SFs_i[sL1JetEt_predict] / dafa_SFs_i[sL1JetEt]
    print("iEtaBins_i: {}".format(iEtaBins_i))
    print("dafa_SFs_i: {}".format(dafa_SFs_i.describe()))
    
    if dafa_SFs is None:
        dafa_SFs = dafa_SFs_i
    else:
        dafa_SFs = pd.concat([dafa_SFs, dafa_SFs_i])    
    
    
    
print("Hello1")    
print("\n\ndafa_SFs: \n{}".format(dafa_SFs.to_string()))
dafa_SFs.to_csv(sOpFileName_SFs, index=False)
print("Wrote {}".format(sOpFileName_SFs))


# In[40]:


sL1JetEt_calib = '%s_calib' % (sL1JetEt)
data_copy1     = data_all[[sL1JetTowerIEtaAbs, sL1JetEt, sGenJetEt]].copy()
dafa_SFs_copy1 = dafa_SFs[[sL1JetTowerIEtaAbs, sL1JetEt, sSF]].copy()
dafa_SFs_copy1 = dafa_SFs_copy1.set_index([sL1JetTowerIEtaAbs, sL1JetEt])
SFs_dict       = dafa_SFs_copy1.to_dict()[sSF]

def calibrateJet(Et_0, iEta):
    Et = round(Et_0)
    if Et < 11: Et = 11
    if Et > 255: Et = 255
    return Et_0 * SFs_dict[(iEta, Et)]

data_copy1[sL1JetEt_calib] = data_copy1.apply(lambda row: calibrateJet(row[sL1JetEt], row[sL1JetTowerIEtaAbs]), axis=1)
#data_copy1[sL1JetEt_calib] = np.vectorize(calibrateJet)(data_copy1[sL1JetEt], data_copy1[sL1JetTowerIEtaAbs])
print("data_copy1: {}".format(data_copy1))


# In[42]:


# SF vs Et plots ----
for iEtaBinRange in np.array_split(iEtaBins, 8):
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5,4), layout='constrained')        
    
    for iEtaBin in iEtaBinRange:
        dafa_SFs_iEtaBin = dafa_SFs[dafa_SFs[sL1JetTowerIEtaAbs] == iEtaBin]
        axs.plot(
            dafa_SFs_iEtaBin[sL1JetEt],
            dafa_SFs_iEtaBin[sSF],
            label='iEta %d' % (iEtaBin)
        )
    axs.set_xlabel('L1JetEt [GeV]')
    axs.set_ylabel('Scale factor')
    axs.set_title('%s' % (sL1JetEt))
    axs.legend()
        
    fig.savefig('%s/SF_vs_Et_%s_ieta_%d_to_%d.png' % (sOutDir, sL1JetEt, iEtaBinRange[0], iEtaBinRange[-1]))
 


# In[43]:


# Resolution plots 
for iEtaBinRange in np.array_split(iEtaBins, 8):
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5,4), layout='constrained')        
    
    for iEtaBin in iEtaBinRange:
        data_copy1_iEtaBin = data_copy1[data_copy1[sL1JetTowerIEtaAbs] == iEtaBin]
        axs.hist(
            (data_copy1_iEtaBin[sL1JetEt_calib]/data_copy1_iEtaBin[sRefJetEt]), 
            bins=100, range=(0, 2.6),
            label='iEta %d' % (iEtaBin),
            histtype='step',#, linewidth=2
            density=True
        )
    axs.set_xlabel('L1JetEt / %s' % (sRefJetEt))
    axs.set_ylabel('Normalized entries')
    axs.set_title('%s' % (sL1JetEt))
    axs.legend()
        
    fig.savefig('%s/AfterJEC_%s_ieta_%d_to_%d.png' % (sOutDir, sL1JetEt, iEtaBinRange[0], iEtaBinRange[-1]))

