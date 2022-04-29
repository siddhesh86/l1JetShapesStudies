#!/bin/env python


import os
import sys
import math 
from array import array
import ROOT as R
R.gROOT.SetBatch(False)  ## Don't print histograms to screen while processing
from collections import OrderedDict
import math


sDir = "plots_troubleshootPhiRingPU"


def printHistoBinContent(h, sSaveAs):
    if not os.path.exists(sDir): os.mkdir(sDir)
    with open("%s/%s.csv" % (sDir, sSaveAs), 'w') as fOut:
        for binY in range(h.GetNbinsY(), 0, -1):
            for binX in range(1, h.GetNbinsX()+1):
                n = h.GetBinContent(binX, binY)
                sTmp = str(n) if n>0 else " "
                fOut.write( "%s\t" % (sTmp) )
            fOut.write( "\n" )


def plotHisto2D(h, sSaveAs):

    c1 = R.TCanvas("c1","c1",1700,1000)
    c1.cd()

    h.GetXaxis().SetTitle("iEta")
    h.GetYaxis().SetTitle("iPhi")
    
    #hEffi.SetMaximum(0.5)
    #hEffi.SetMinimum(0.0)

    #R.gStyle.SetPaintTextFormat("4.2f");
    #h.SetMarkerColor(2)
    h.SetMarkerSize(0.5)
    h.Draw('colz TEXT')
    #h.Draw('colz')



    line = R.TLine()
    line.SetLineStyle(3)
    line.SetLineColor(R.kGray+1)
    #line.DrawLineNDC(-20, 5, -20, 70)
    #line.DrawLine(-20.5, 5.5, -20.5, 70.5)
    #line.DrawLine(-30.5, 5.5, -30.5, 70.5)
    iEtaRange = [-41.5, 41.5]
    iPhiRange = [0.5, 72]
    i = iEtaRange[0]
    while i <= iEtaRange[1]:
        line.DrawLine(i, iPhiRange[0], i, iPhiRange[1])
        i += 1
    i = iPhiRange[0]
    while i <= iPhiRange[1]:
        line.DrawLine(iEtaRange[0], i, iEtaRange[1], i)
        i += 1

    c1.Update()
    
    if not os.path.exists(sDir): os.mkdir(sDir)
    c1.SaveAs("%s/%s.png" % (sDir, sSaveAs))
    
    

if __name__ == '__main__':


    
    sInFile = "L1T_HCALL2Calib_stage1_PFA1p_nVtxAll_part0_of_1.root"
    sHistoDir = "caloTTs/"

    evts = [
        '1:9983:9982010',
        '1:9983:9982116',
        '1:9983:9982086',
        '1:9983:9982118', 
    ]

    sHisto_dict = OrderedDict()
    for evt0 in evts:
        evt = evt0.replace(":", "_")
        sHisto_dict[evt] = [
            'hCaloTowers_iEta_vs_iPhi_%s' % evt,
            'hCaloTTs_iEta_vs_iPhi_%s' % evt,            
        ]
    

    R.gStyle.SetOptStat(0)
    
    fIn = R.TFile(sInFile)
    if not fIn.IsOpen():
        print("{} could not open".format(sInFile))


    for evt0, sHistos in sHisto_dict.items():
        evt = evt0.replace(":", "_")
        histos = []
        for sHisto0 in sHistos:
            sHisto = "%s%s" % (sHistoDir, sHisto0)
            h = fIn.Get(sHisto)
            if not h:
                print("%s<<histogram not found" % (sHisto0))
                continue
            #histos.append( h )

            sSaveAs = sHisto0.replace(evt, "")
            sSaveAs = "%s_%s" % (evt, sSaveAs)

                       
            plotHisto2D(h, sSaveAs)

            printHistoBinContent(h, sSaveAs)
