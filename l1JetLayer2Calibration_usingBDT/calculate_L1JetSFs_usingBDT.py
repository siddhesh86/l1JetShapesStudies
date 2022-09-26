#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parseGroup1 = parser.add_mutually_exclusive_group(required=False)
parseGroup1.add_argument('--ChunkyDonut', action='store_true', default=True)
parseGroup1.add_argument('--PhiRing',     action='store_true', default=False)

#args = parser.parse_args()
args = parser.parse_args("--ChunkyDonut".split()) # to run in jupyter-notebook 
l1Jet_ChunkyDonut = args.ChunkyDonut
l1Jet_PhiRing     = args.PhiRing

sIpFileName = "../data/L1T_Jet_MLInputs_Run3_QCD_Pt15to7000_PFA1p_CMSSW12_6_0_pre1_nVtxAll_20220925.csv"


data_all = pd.read_csv(sIpFileName)
print("Input file: %s" % (sIpFileName))


# In[2]:


print("data_all.columns: {}, \ndata_all.shape: {}".format(data_all.columns, data_all.shape))


# In[8]:


data_all['L1JetEt_PUS_ChunkyDonut'] = data_all['L1Jet9x9_RawEt'] - data_all['L1Jet9x9_PUEt_ChunkyDonut']

data_all['L1JetEt_PUS_PhiRing']     = data_all['L1Jet9x9_RawEt'] - (data_all['L1Jet9x9_EtSum7PUTowers'] / 7.0 )

print("data_all.describe(): {}".format(data_all.describe()))


# In[ ]:




