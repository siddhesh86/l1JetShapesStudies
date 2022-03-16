#! /usr/bin/env python

## ******************************************************************************************** ##
##  Look at scale and resolution of L1T jets and MET with different input TPs - draw extension  ##
## ******************************************************************************************** ##

import os
import sys
from collections import OrderedDict
from array import array
from ctypes import *
import time

import ROOT as R
R.gROOT.SetBatch(False)  ## Don't print histograms to screen while processing


kNullLargeNeg = -9999.0

def getHists1DYRange(histos_lits):
    yMax = -99999.0
    yMin =  99999.0
    for iHisto in range(len(histos_lits)):
        h = histos_lits[iHisto]
        for iBin in range(1, h.GetNbinsX()+1):
            # ignore bins with >N% error
            if h.GetBinError(iBin) >= abs(h.GetBinContent(iBin)) * 0.30:
                continue
            
            yplus  = h.GetBinContent(iBin) + h.GetBinError(iBin)
            yminus = h.GetBinContent(iBin) - h.GetBinError(iBin)
            if yplus  > yMax: yMax = yplus
            if yminus < yMin: yMin = yminus

    if yMax > 0: yMax *= 1.3 # 
    else:        yMax *= 0.7
    if yMin > 0: yMin *= 0.8
    else:        yMin *= 1.3
        
    return yMin, yMax;    

#def plotHistos1DAndRatioPlot(h_list, legend, canvas, sSaveAs, hRatio_list = []):
def plotHistos1DAndRatioPlot(h_list, hRatio_list , legend, canvas, sSaveAs, drawOption="PE", pad1SetLogY=False, yRange=[], setLogX=False):
    print "type(h_list[0]): {} ".format(type(h_list[0]))
    # yRange=[yMinPad1, yMaxPad1,   yMinPad2, yMaxPad2]
    print "plotHistos1DAndRatioPlot():: h_list ({}), hRatio_list ({}) ".format(len(h_list), len(hRatio_list))
    for ih in range(0,len(h_list)):
        h = h_list[ih]
        #print "hist %s: %g " % (h.GetName(), h.GetBinContent(5))
        
    
    if len(hRatio_list) == 0: 
        for ih in range(1,len(h_list)):
            hRatio = h_list[ih].Clone("%s_Ratio" % (h_list[ih].GetName()))
            hRatio.Divide(h_list[0])

            '''
            for iBin in range(1, hRatio.GetNbinsX()+1):
                binContent = hRatio.GetBinContent(iBin)
                if abs(binContent - 0) < 1e-6: hRatio.SetBinContent(iBin, kNullLargeNeg)
                else:                          hRatio.SetBinContent(iBin, binContent - 1.0)
            '''
            hRatio_list.append(hRatio)
        else:
            if (len(h_list) - 1) != len(hRatio_list):
                print "plotHistos1DAndRatioPlot():: (len(h_list) - 1) != len(hRatio_list) \t *** ERROR ***"
                exit(0)

    drawOption0 = drawOption1 = drawOption
    if 'TGraph' in str(type(h_list[0])): drawOption1 = drawOption1.replace('A','') # don't use option 'A' for grpoh ploting on the same canvas
    
    canvas.cd()
    pad1 = R.TPad("pad1","pad2",0,0.3, 1,1)
    if pad1SetLogY: pad1.SetLogy()
    if setLogX: pad1.SetLogx()
    pad1.Draw()
    pad1.SetGridx()
    pad1.SetGridy()

    pad2 = R.TPad("pad2","pad2",0,0, 1,0.3)
    if setLogX: pad2.SetLogx()
    pad2.Draw()
    pad2.SetGridx()
    pad2.SetGridy()
    
    pad1.SetBottomMargin(0.01);
    
    pad2.SetTopMargin(0.03);
    pad2.SetBottomMargin(0.3);
    
    pad1.SetRightMargin(0.05);
    pad2.SetRightMargin(0.05);
    
    pad1.SetLeftMargin(0.13);
    pad2.SetLeftMargin(0.13);

    ymin = ymax = yRatiomin = yRatiomax = None
    if len(yRange) < 2: ymin, ymax           = getHists1DYRange(h_list)
    if len(yRange) < 4: yRatiomin, yRatiomax = getHists1DYRange(hRatio_list)
    print "plotHistos1DAndRatioPlot:: ymin {}, ymax {}, yRatiomin {}, yRatiomax {}".format(ymin, ymax, yRatiomin, yRatiomax)
    if len(yRange) < 4 and yRatiomin < 0: yRatiomin = 0
    for ih in range(0,len(h_list)):
        pad1.cd()
        h = h_list[ih]
        if (ih == 0):
            if len(yRange) >= 2: h.GetYaxis().SetRangeUser(yRange[0], yRange[1]) 
            else:                h.GetYaxis().SetRangeUser(ymin, ymax)     
            h.Draw("%s" % (drawOption0))
            #print "hist %s - 1: %g " % (h.GetName(), h.GetBinContent(5))
        else:
            h.Draw("same %s" % (drawOption1))
            #print "hist %s - 2: %g " % (h.GetName(), h.GetBinContent(5))

    for ihRatio in range(0, len(hRatio_list)):
        # ratio plots
        pad2.cd()
        hRatio = hRatio_list[ihRatio]
        if (ihRatio == 0):
            #hRatio.GetYaxis().SetRangeUser(yRatiomin, yRatiomax)
            #hRatio.GetYaxis().SetRangeUser(0, 2)
            if len(yRange) >= 4: hRatio.GetYaxis().SetRangeUser(yRange[2], yRange[3])
            #else:                hRatio.GetYaxis().SetRangeUser(0.5, 1.5)
            else:                hRatio.GetYaxis().SetRangeUser(0, 2)
            
            #hRatio.GetYaxis().SetTitle("Ratio w.r.t. 1st");
            
            xLabSize   = hRatio.GetXaxis().GetLabelSize();
            xTitleSize = hRatio.GetXaxis().GetTitleSize();
            xOffSet    = hRatio.GetXaxis().GetTitleOffset();
            yLabSize   = hRatio.GetYaxis().GetLabelSize();
            yTitleSize = hRatio.GetYaxis().GetTitleSize();
            #yOffSet = hRatio.GetYaxis().GetTitleOffset();
            
            hRatio.GetXaxis().SetTitleSize(3.2*xTitleSize);
            hRatio.GetXaxis().SetLabelSize(2.3*xLabSize);
            #hRatio.GetXaxis().SetLabelSize(1.7*xLabSize);
            hRatio.GetXaxis().SetTitleOffset(0.8*xOffSet);
            
            hRatio.GetYaxis().SetTitleSize(1.5*yTitleSize); # 3.0 
            hRatio.GetYaxis().SetLabelSize(2.3*yLabSize);
            
            #hRatio.GetYaxis().SetTitleOffset(0.32);
            hRatio.GetYaxis().SetTitleOffset(0.72);
            
            hRatio.GetYaxis().SetNdivisions(505);
            
            hRatio.Draw("%s" % (drawOption0))
        else:
            hRatio.Draw("same %s" % (drawOption1))
    
    pad1.cd()
    legend.Draw()
    canvas.Update()
    canvas.SaveAs(sSaveAs)
    #input("Enter something")
    return canvas


def cloneHistogram(histo, sSubstringForCloneName):
    h = histo.Clone("%s_%s" % (histo.GetName(), sSubstringForCloneName))
    return h


def DivideGraphAsymmErr(grN, grD):
    if grN.GetN() != grD.GetN():
        print "DivideGraphAsymmErr():: grN.GetN() ({}) != grD.GetN() ({}) \t *** ERROR *** \n".format(grN.GetN(), grD.GetN())
        #input("\nEnter anything")
        return None
    
    nPts = min(grN.GetN(), grD.GetN())
    x = []
    y = []
    for iPt in range(nPts):
        #xN = yN = xD = yD = None
        #xN = yN = xD = yD = ctypes.c_double(1.)
        #xN = yN = xD = yD = 1.0
        xN = grN.GetPointX(iPt);   yN = grN.GetPointY(iPt);
        xD = grD.GetPointX(iPt);   yD = grD.GetPointY(iPt);
        '''
        if (grN.GetPoint(iPt, xN, yN) == -1) or (grD.GetPoint(iPt, xD, yD) == -1):
            print "DivideGraphAsymmErr():: error in reading graph point {} \t *** ERROR *** \n".format(iPt)
            return
        '''
        
        x.append(xN)
        if abs(yD - 0) > 1e-5: y.append( yN / yD )
        else:                  y.append( 0. )
    
    gr = R.TGraphAsymmErrors(nPts, array('d', x), array('d', y))
    
    '''
    print "print gr::"
    for iPt in range(nPts):
        xN = grN.GetPointX(iPt);   yN = grN.GetPointY(iPt);
        xD = grD.GetPointX(iPt);   yD = grD.GetPointY(iPt);
        x1 = gr.GetPointX(iPt);    y1 = gr.GetPointY(iPt);
        if yD > 0:
            print "grN: ({}, {}), grD: ({}, {}), yN/yD: {}, gr: ({}, {})".format(xN,yN, xD,yD, yN/yD, x1,y1)
    '''
    
    return gr


def run(resPtBin):
    #in_fileNames = ["L1T_JetMET_Res_nVtxMin_45_nSamp_2_nPresamp_0_HB_1p0_1p0_HE1_1p0_1p0_HE2_1p0_1p0_2018AB_-1k.root"]
    #in_fileNames = ["L1T_JetMET_Res_def.root", "L1T_JetMET_Res_PFA1p.root"]
    #in_fileNames = ["L1T_JetMET_Res_def.root"]
    #in_fileNames = ["L1T_JetMET_Res_def_hadded.root", "L1T_JetMET_Res_PFA1p_hadded.root"]; sInFileVersion = "_wPFJetFilters";
    #in_fileNames = ["L1T_JetMET_Res_def_hadded.root", "L1T_JetMET_Res_PFA1p_hadded.root"]; sInFileVersion = "_wLUTGenModeTrue_PFA1pRun3ContainPhaseNSm2";
    #in_fileNames = ["L1T_JetMET_Res_def_hadded.root"]
    #in_fileNames = ["L1T_JetMET_Res_def_CalibJetByHand_hadded.root", ]; sInFileVersion = "_afterLayer2Calib";
    #in_fileNames = ["L1T_JetMET_Res_def_woPUS_hadded.root", "L1T_JetMET_Res_PFA1p_woPUS_hadded.root"]; sInFileVersion = "_woPUS";
    #in_fileNames = ["L1T_JetMET_Res_def_hadded.root"]; sInFileVersion = "_wLUTGenModeTrue_PFA1pRun3ContainPhaseNSm2";
    #in_fileNames = ["L1T_JetMET_Res_def_LowAndHighPU_hadded.root"]; sInFileVersion = "_wLUTGenModeTrue_PFA1pRun3ContainPhaseNSm2_v5_1";

    # 'run_TrgRatesAndEffi'
    #in_fileNames = ["L1T_JetMET_Res_def_2018_ZeroBias_nVts_Lt25_and_Gt50_wCalibJetByHand.root"]; sInFileVersion = "_wLUTGenModeTrue_PFA1pRun3ContainPhaseNSm2_Rates_2018_ZeroBias";
    #in_fileNames = ["L1T_JetMET_Res_def_2018_ZeroBias_nVts_Lt25_and_Gt50_wCalibJetByHand_calibSFAtPU50to100_v2.root"]; sInFileVersion = "_wLUTGenModeTrue_PFA1pRun3ContainPhaseNSm2_Rates_2018_ZeroBias_calibSFAtPU50to100_set2"; # need "set1" or "set2" string in sInFileVersion
   # in_fileNames = ["L1T_JetMET_Res_def_2018_ZeroBias_nVts_Lt25_and_Gt50_wCalibJetByHand_calibSFAtPU1to25_v2.root"]; sInFileVersion = "_wLUTGenModeTrue_PFA1pRun3ContainPhaseNSm2_Rates_2018_ZeroBias_calibSFAtPU1to25_set1"; # need "set1" or "set2" string in sInFileVersion
    
    # 'run_TrgEffi'
    #in_fileNames = ["L1T_JetMET_Res_def_2018_SingleMu_nVts_Lt25_and_Gt50_wCalibJetByHand.root"]; sInFileVersion = "_wLUTGenModeTrue_PFA1pRun3ContainPhaseNSm2_Effi_2018_SingleMu";
    #in_fileNames = ["L1T_JetMET_Res_def_2018_SingleMu_nVts_Lt25_and_Gt50_wCalibJetByHand_calibSFAtPU50to100_v2.root"]; sInFileVersion = "_wLUTGenModeTrue_PFA1pRun3ContainPhaseNSm2_Effi_2018_SingleMu_calibSFAtPU50to100_set2"; # need "set1" or "set2" string in sInFileVersion
    #in_fileNames = ["L1T_JetMET_Res_def_2018_SingleMu_nVts_Lt25_and_Gt50_wCalibJetByHand_calibSFAtPU1to25_v2.root"]; sInFileVersion = "_wLUTGenModeTrue_PFA1pRun3ContainPhaseNSm2_Effi_2018_SingleMu_calibSFAtPU1to25_set2"; # need "set1" or "set2" string in sInFileVersion


    # 'L1TauMatchedToPFJetVersions'
    # 'run_TrgEffi'
    #in_fileNames = ["L1T_JetMET_Res_def_2018_SingleMu_nVts_Lt25_and_Gt50_wCalibJetByHand_calibSFAtPU50to100_v5.root"]; sInFileVersion = "_wLUTGenModeTrue_PFA1pRun3ContainPhaseNSm2_L1TauMatchedToPFJet_v5_set3"; # need "set1" or "set2" or "set3" string in sInFileVersion
    # 'run_TrgRate'
    #in_fileNames = ["L1T_JetMET_Res_def_2018_ZeroBias_nVts_Lt25_and_Gt50_wCalibJetByHand_calibSFAtPU50to100_v5.root"]; sInFileVersion = "_wLUTGenModeTrue_PFA1pRun3ContainPhaseNSm2_L1TauMatchedToPFJet_v5_set3"; # need "set1" or "set2" or "set3" string in sInFileVersion


    #
    #
    #in_fileNames = ["L1T_JetMET_Res_def_2018_SingleMu_nVts_Lt25_and_Gt50_v6.root"]; sInFileVersion = "_wLUTGenModeTrue_PFA1pRun3ContainPhaseNSm2_L1TauMatchedToPFJet_v6_set3"; # need "set1" or "set2" or "set3" string in sInFileVersion
    # 'run_TrgEffi'
    #in_fileNames = ["L1T_JetMET_Res_def_2018_SingleMu_nVts_Lt25_and_Gt50_wCalibJetByHand_calibSFAtPU50to100_Pt25to35_v7.root"]; sInFileVersion = "_wLUTGenModeTrue_PFA1pRun3ContainPhaseNSm2_L1TauMatchedToPFJet_v7_set3"; # need "set1" or "set2" or "set3" string in sInFileVersion
    # 'run_TrgRate'
    #in_fileNames = ["L1T_JetMET_Res_def_2018_ZeroBias_nVts_Lt25_and_Gt50_wCalibJetByHand_calibSFAtPU50to100_PFPt25to35_v7.root"]; sInFileVersion = "_wLUTGenModeTrue_PFA1pRun3ContainPhaseNSm2_L1TauMatchedToPFJet_v7_set3"; # need "set1" or "set2" or "set3" string in sInFileVersion
    
    #
    #
    #in_fileNames = ["L1T_JetMET_Res_def_2018_SingleMu_nVts_Lt25_and_Gt50_v8.root"]; sInFileVersion = "_wLUTGenModeTrue_PFA1pRun3ContainPhaseNSm2_L1TauMatchedToPFJet_v8_set3"; # need "set1" or "set2" or "set3" string in sInFileVersion
    # 'run_TrgEffi'
    #in_fileNames = ["L1T_JetMET_Res_def_2018_SingleMu_nVts_Lt25_and_Gt50_wCalibJetByHand_calibSFAtPU50to100_Pt25to35_v9.root"]; sInFileVersion = "_wLUTGenModeTrue_PFA1pRun3ContainPhaseNSm2_L1TauMatchedToPFJet_v9_set3"; # need "set1" or "set2" or "set3" string in sInFileVersion
    # 'run_TrgRate'
    #in_fileNames = ["L1T_JetMET_Res_def_2018_ZeroBias_nVts_Lt25_and_Gt50_wCalibJetByHand_calibSFAtPU50to100_PFPt25to35_v9.root"]; sInFileVersion = "_wLUTGenModeTrue_PFA1pRun3ContainPhaseNSm2_L1TauMatchedToPFJet_v9_set3"; # need "set1" or "set2" or "set3" string in sInFileVersion
    
    #
    #
    in_fileNames = ["L1T_JetMET_Res_def_2018_SingleMu_nVts_Lt25_and_Gt50_v10.root"]; sInFileVersion = "_wLUTGenModeTrue_PFA1pRun3ContainPhaseNSm2_L1TauMatchedToPFJet_v10_set3"; # need "set1" or "set2" or "set3" string in sInFileVersion
    # 'run_TrgEffi'
    #in_fileNames = ["L1T_JetMET_Res_def_2018_SingleMu_nVts_Lt25_and_Gt50_wCalibJetByHand_calibSFAtPU50to100_Pt25to35_v11.root"]; sInFileVersion = "_wLUTGenModeTrue_PFA1pRun3ContainPhaseNSm2_L1TauMatchedToPFJet_v11_set3"; # need "set1" or "set2" or "set3" string in sInFileVersion
    # 'run_TrgRate'
    #in_fileNames = ["L1T_JetMET_Res_def_2018_ZeroBias_nVts_Lt25_and_Gt50_wCalibJetByHand_calibSFAtPU50to100_PFPt25to35_v11.root"]; sInFileVersion = "_wLUTGenModeTrue_PFA1pRun3ContainPhaseNSm2_L1TauMatchedToPFJet_v11_set3"; # need "set1" or "set2" or "set3" string in sInFileVersion

    
    sInFileTags  = ["PFA2", "PFA1p"]
    
    #JetShapes = ['9x9', '8x9', '7x9', '6x9', '5x9', '4x9', '3x9']
    #JetShapes = ['9x9', '8x9', '7x9', '6x9', '5x9', '4x9', '3x9', '3x9_plus_0.5_times_9x9']
    JetShapes = ['Default'] + ['9x9', '8x9', '7x9', '6x9', '5x9', '4x9', '3x9', '3x9_plus_0.5_times_9x9']    

    # now pasing as main() function argument
    #resPtBin = ["medPt"] # "medPt"  "PtAllBins"

    resTypes   = [
        "h_jet_byHand_res_vs_iEta$JETSHAPE_Raw_HBEF_%s_0" % (resPtBin),
        "h_jet_byHand_res_vs_iEta$JETSHAPE_RawPUS_HBEF_%s_0" % (resPtBin),
        "h_jet_byHand_res_vs_iEta$JETSHAPE_RawPUS_phiRingMin4_HBEF_%s_0" % (resPtBin),
        "h_jet_byHand_res_vs_iEta$JETSHAPE_RawPUS_phiRingSide4_HBEF_%s_0" % (resPtBin),
        "h_jet_byHand_res_vs_iEta$JETSHAPE_RawPUS_phiRingAdjacent_HBEF_%s_0" % (resPtBin),
        "h_jet_byHand_res_vs_iEta$JETSHAPE_RawPUS_phiDefault_HBEF_%s_0" % (resPtBin),
    ]
    sresTypes_hRatio_Denominator_forDefault = "_RawPUS_HBEF_" # string to identify denominator histogram to take ratio for other histograms
    sresTypes_hRatio_Denominator_forPhiRing = "_RawPUS_phiRingMin4_" # string to identify denominator histogram to take ratio for other histograms
    
    sResLegend = [
        "Raw",
        "RawPUS",
        "RawPUS_phiRingMin4",
        "RawPUS_phiRingSide4",
        "RawPUS_phiRingAdjacent",
        "RawPUS_phiRingDefault"
    ];

    effiTypes   = [
        "h_jet_byHand_eff$JETSHAPE_Raw_IETA_PTBIN_0",
        "h_jet_byHand_eff$JETSHAPE_RawPUS_IETA_PTBIN_0",
        "h_jet_byHand_eff$JETSHAPE_RawPUS_phiRingMin4_IETA_PTBIN_0",
        "h_jet_byHand_eff$JETSHAPE_RawPUS_phiRingSide4_IETA_PTBIN_0",
        "h_jet_byHand_eff$JETSHAPE_RawPUS_phiRingAdjacent_IETA_PTBIN_0",
    ]    
    

    l1Modes = ["emu"] # ["emu", "hw"]


    resType_selected = "h_jet_byHand_res_vs_iEta$JETSHAPE_RawPUS_phiRingMin4_HBEF_%s_0" % (resPtBin)
    sResLegend_selected = "RawPUS_phiRingMin4"
    #resType_selected = "h_jet_byHand_res_vs_iEta$JETSHAPE_RawPUS_phiDefault_HBEF_%s_0" % (resPtBin)
    #sResLegend_selected = "RawPUS_phiRingDefault"
    resType_default = "h_jet_byHand_res_vs_iEta_RawPUS_HBEF_%s_0" % (resPtBin)
    sResLegend_default = "Default"
    
    '''
    JetShapesAndPUSs = OrderedDict([
        ('9x9', OrderedDict([
            ("Raw",                   "h_jet_byHand_res_vs_iEta_vs_nVtx$JETSHAPE_Raw_HBEF_%s_0" % (resPtBin)),
            ("RawPUS",                "h_jet_byHand_res_vs_iEta_vs_nVtx$JETSHAPE_RawPUS_HBEF_%s_0" % (resPtBin)),
            ("RawPUS_phiRingDefault", "h_jet_byHand_res_vs_iEta_vs_nVtx$JETSHAPE_RawPUS_phiDefault_HBEF_%s_0" % (resPtBin)),
            ("RawPUS_phiRingMin4",    "h_jet_byHand_res_vs_iEta_vs_nVtx$JETSHAPE_RawPUS_phiRingSide4_HBEF_%s_0" % (resPtBin))
        ])),
        
        ('3x9_plus_0.5_times_9x9', OrderedDict([
            ("RawPUS_phiRingDefault", "h_jet_byHand_res_vs_iEta_vs_nVtx$JETSHAPE_RawPUS_phiDefault_HBEF_%s_0" % (resPtBin)),
            ("RawPUS_phiRingMin4",    "h_jet_byHand_res_vs_iEta_vs_nVtx$JETSHAPE_RawPUS_phiRingSide4_HBEF_%s_0" % (resPtBin))
        ]))
        
    ])
    JetShapesAndPUSs_Denom_forRatioPlot = ['9x9', "RawPUS", "lowPU"] # dict key names used to refer histogram to be used as a denominator histogram in ratio plot
    '''

    JetShapesAndPUSs = OrderedDict([
       ('Default', OrderedDict([
            ("Raw",                   "h_jet_byHand_res_vs_iEta_vs_nVtx$JETSHAPE_Raw_HBEF_%s_0" % (resPtBin)),
            ("RawPUS",                "h_jet_byHand_res_vs_iEta_vs_nVtx$JETSHAPE_RawPUS_HBEF_%s_0" % (resPtBin)),
            ("L1 jets",               "h_jet_byHand_res_vs_iEta_vs_nVtx$JETSHAPE_L1JDefault_HBEF_%s_0" % (resPtBin)),
        ])),        
        
        ('9x9', OrderedDict([
            ("RawPUS_phiRingDefault", "h_jet_byHand_res_vs_iEta_vs_nVtx$JETSHAPE_RawPUS_phiDefault_HBEF_%s_0" % (resPtBin)),
            #("RawPUS_phiRingMin4",    "h_jet_byHand_res_vs_iEta_vs_nVtx$JETSHAPE_RawPUS_phiRingSide4_HBEF_%s_0" % (resPtBin))
        ])),

        ('3x9', OrderedDict([
            ("RawPUS_phiRingDefault", "h_jet_byHand_res_vs_iEta_vs_nVtx$JETSHAPE_RawPUS_phiDefault_HBEF_%s_0" % (resPtBin)),
#            ("RawPUS_phiRingMin4",    "h_jet_byHand_res_vs_iEta_vs_nVtx$JETSHAPE_RawPUS_phiRingSide4_HBEF_%s_0" % (resPtBin))
        ])),
        
        ('3x9_plus_0.5_times_9x9', OrderedDict([
            ("RawPUS_phiRingDefault", "h_jet_byHand_res_vs_iEta_vs_nVtx$JETSHAPE_RawPUS_phiDefault_HBEF_%s_0" % (resPtBin)),
#            ("RawPUS_phiRingMin4",    "h_jet_byHand_res_vs_iEta_vs_nVtx$JETSHAPE_RawPUS_phiRingSide4_HBEF_%s_0" % (resPtBin))
        ])),

       ('L1TauDefault', OrderedDict([
            ("Et",                   "h_jet_byHand_res_vs_iEta_vs_nVtx$JETSHAPE_Et_HBEF_%s_0" % (resPtBin)),
            ("RawEt",                "h_jet_byHand_res_vs_iEta_vs_nVtx$JETSHAPE_RawEt_HBEF_%s_0" % (resPtBin)),
        ])),        
         
        
    ])    


    JetShapesAndPUSs_histoName_TrgRates_1 = OrderedDict([
       ('Default', OrderedDict([
           ("RawPUS",                "h_jet_byHand_rates_singleJet$JETSHAPE_RawPUS_$ETACAT_0"),
           ("Raw",                   "h_jet_byHand_rates_singleJet$JETSHAPE_Raw_$ETACAT_0"),
           ("jet (PUS)",                "h_jet_byHand_rates_singleJet$JETSHAPE_L1JDefault_$ETACAT_0"),
        ])),        
        
        ('9x9', OrderedDict([
            ("RawPUS_phiRingDefault", "h_jet_byHand_rates_singleJet$JETSHAPE_RawPUS_phiDefault_$ETACAT_0"),
            ##("RawPUS_phiRingMin4",    "h_jet_byHand_rates_singleJet$JETSHAPE_Raw_HBEF_0")
        ]))
    ])
    JetShapesAndPUSs_histoName_TrgRates_2 = OrderedDict([
        ('9x9', OrderedDict([
            ("RawPUS_phiRingDefault", "h_jet_byHand_rates_singleJet$JETSHAPE_RawPUS_phiDefault_$ETACAT_0"),
            ##("RawPUS_phiRingMin4",    "h_jet_byHand_rates_singleJet$JETSHAPE_Raw_HBEF_0")
        ])),

        ('3x9', OrderedDict([
            ("RawPUS_phiRingDefault", "h_jet_byHand_rates_singleJet$JETSHAPE_RawPUS_phiDefault_$ETACAT_0"),
            ##("RawPUS_phiRingMin4",    "h_jet_byHand_res_vs_iEta_vs_nVtx$JETSHAPE_RawPUS_phiRingSide4_HBEF_%s_0" % (resPtBin))
        ])),
        
        ('3x9_plus_0.5_times_9x9', OrderedDict([
            ("RawPUS_phiRingDefault", "h_jet_byHand_rates_singleJet$JETSHAPE_RawPUS_phiDefault_$ETACAT_0"),
            ##("RawPUS_phiRingMin4",    "h_jet_byHand_rates_singleJet$JETSHAPE_Raw_HBEF_0")
        ]))
    ])
    JetShapesAndPUSs_histoName_TrgRates_3 = OrderedDict([
       ('Default', OrderedDict([
           ("RawPUS",                 "h_jet_byHand_rates_singleJet$JETSHAPE_RawPUS_$ETACAT_0"),
           ("Jet (PUS)",              "h_jet_byHand_rates_singleJet$JETSHAPE_L1JDefault_$ETACAT_0"),
        ])),        
        
        ('9x9', OrderedDict([
            ("RawPUS_phiRingDefault", "h_jet_byHand_rates_singleJet$JETSHAPE_RawPUS_phiDefault_$ETACAT_0"),
            ##("RawPUS_phiRingMin4",    "h_jet_byHand_rates_singleJet$JETSHAPE_Raw_HBEF_0")
        ])),
        
       ('L1TauDefault', OrderedDict([
           ("Et",                   "h_jet_byHand_rates_singleJet$JETSHAPE_Et_$ETACAT_0"),
           ("RawEt",                "h_jet_byHand_rates_singleJet$JETSHAPE_RawEt_$ETACAT_0"),
        ])),        
        
    ])
    JetShapesAndPUSs_histoName_TrgRates = None
    if   "set1" in sInFileVersion:
        JetShapesAndPUSs_histoName_TrgRates = JetShapesAndPUSs_histoName_TrgRates_1
    elif "set2" in sInFileVersion:
        JetShapesAndPUSs_histoName_TrgRates = JetShapesAndPUSs_histoName_TrgRates_2
    elif "set3" in sInFileVersion:
        JetShapesAndPUSs_histoName_TrgRates = JetShapesAndPUSs_histoName_TrgRates_3
    
    TrgRates_Types = ['singleJet', 'doubleJet', 'trippleJet', 'quadJet']
    '''
    TrgEffi_HistNameReplacements = [
        'jet_byHand_rates_singleJet', # replace this string from JetShapesAndPUSs_histoName_TrgRates by the following strings
        'jet_byHand_eff_num_vs_PU'
    ]
    '''
    sHistoName_nTotalEvents_forTrgRates = "h_nTotalEvents"
    #// Run3 LHC parameters for normalizing rates
    instLumi = 2e34; #// Hz/cm^2, https://indico.cern.ch/event/880508/contributions/3720014/attachments/1980197/3297287/CMS-Week_20200203_LHCStatus_Schaumann_v2.pdf
    mbXSec   = 6.92e-26; #// cm^2, minimum bias cross section from Run2: https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData#Recommended_cross_section
    axisRange_EtTrsh_forTrgRates = [10, 50] # [50, 100]

    JetShapesAndPUSs_histoName_TrgEffi_1 = OrderedDict([
       ('Default', OrderedDict([
           ("RawPUS",                "h_jet_byHand_eff_num_vs_PU$JETSHAPE_RawPUS_$ETACAT_TrgTrsh$TRGTRSH_0"),
           ("Raw",                   "h_jet_byHand_eff_num_vs_PU$JETSHAPE_Raw_$ETACAT_TrgTrsh$TRGTRSH_0"),
           ("jet (PUS)",             "h_jet_byHand_eff_num_vs_PU$JETSHAPE_L1JDefault_$ETACAT_TrgTrsh$TRGTRSH_0"),
        ])),        
        
        ('9x9', OrderedDict([
            ("RawPUS_phiRingDefault", "h_jet_byHand_eff_num_vs_PU$JETSHAPE_RawPUS_phiDefault_$ETACAT_TrgTrsh$TRGTRSH_0"),
            #("RawPUS_phiRingMin4",    "h_jet_byHand_rates_singleJet$JETSHAPE_Raw_HBEF_0")
        ])),
    ])
    JetShapesAndPUSs_histoName_TrgEffi_2 = OrderedDict([
        ('9x9', OrderedDict([
            ("RawPUS_phiRingDefault", "h_jet_byHand_eff_num_vs_PU$JETSHAPE_RawPUS_phiDefault_$ETACAT_TrgTrsh$TRGTRSH_0"),
            #("RawPUS_phiRingMin4",    "h_jet_byHand_rates_singleJet$JETSHAPE_Raw_HBEF_0")
        ])),
        
        ('3x9', OrderedDict([
            ("RawPUS_phiRingDefault", "h_jet_byHand_eff_num_vs_PU$JETSHAPE_RawPUS_phiDefault_$ETACAT_TrgTrsh$TRGTRSH_0"),
            #("RawPUS_phiRingMin4",    "h_jet_byHand_rates_singleJet$JETSHAPE_Raw_HBEF_0")
        ])),
        
        ('3x9_plus_0.5_times_9x9', OrderedDict([
            ("RawPUS_phiRingDefault", "h_jet_byHand_eff_num_vs_PU$JETSHAPE_RawPUS_phiDefault_$ETACAT_TrgTrsh$TRGTRSH_0"),
#            ("RawPUS_phiRingMin4",    "h_jet_byHand_rates_singleJet$JETSHAPE_Raw_HBEF_0")
        ]))
        
    ])
    JetShapesAndPUSs_histoName_TrgEffi_3 = OrderedDict([
       ('Default', OrderedDict([
           ("RawPUS",            "h_jet_byHand_eff_num_vs_PU$JETSHAPE_RawPUS_$ETACAT_TrgTrsh$TRGTRSH_0"),
           ("Jet (PUS)",             "h_jet_byHand_eff_num_vs_PU$JETSHAPE_L1JDefault_$ETACAT_TrgTrsh$TRGTRSH_0"),
        ])),
        
        ('9x9', OrderedDict([
            ("RawPUS_phiRingDefault", "h_jet_byHand_eff_num_vs_PU$JETSHAPE_RawPUS_phiDefault_$ETACAT_TrgTrsh$TRGTRSH_0"),
            #("RawPUS_phiRingMin4",    "h_jet_byHand_rates_singleJet$JETSHAPE_Raw_HBEF_0")
        ])),
        
        ('L1TauDefault', OrderedDict([
           ("Et",                "h_jet_byHand_eff_num_vs_PU$JETSHAPE_Et_$ETACAT_TrgTrsh$TRGTRSH_0"),
           ("RawEt",             "h_jet_byHand_eff_num_vs_PU$JETSHAPE_RawEt_$ETACAT_TrgTrsh$TRGTRSH_0"),
        ])),
        
    ])
    JetShapesAndPUSs_histoName_TrgEffi = None
    if   "set1" in sInFileVersion:
        JetShapesAndPUSs_histoName_TrgEffi = JetShapesAndPUSs_histoName_TrgEffi_1
    elif "set2" in sInFileVersion:
        JetShapesAndPUSs_histoName_TrgEffi = JetShapesAndPUSs_histoName_TrgEffi_2
    elif "set3" in sInFileVersion:
        JetShapesAndPUSs_histoName_TrgEffi = JetShapesAndPUSs_histoName_TrgEffi_3
    TrgEffi_TrgTrshs = [55] #  [75] #[60, 90]  [30]
    axisRange_PFJetPt_forTrgEffi = [10, 100] # [40, 130]
    #axisRebin_PFJetPt_forTrgEffi = int( (axisRange_PFJetPt_forTrgEffi[1] - axisRange_PFJetPt_forTrgEffi[0]) / 18 )
    axisRebin_PFJetPt_forTrgEffi = int( (axisRange_PFJetPt_forTrgEffi[1] - axisRange_PFJetPt_forTrgEffi[0]) / 18 )
    
    nPVRanges = OrderedDict([
        #("lowPU",  [ 1,  25]),
        #("nVtx #in [1, 25]",  [ 1,  25]),
        ("PU1to25",  [ 1,  25]),
        #"medPU":  [50,  54],
        #("highPU", [50, 100]), # [55, 100]
        #("nVtx #in [50, 100]", [50, 100]), # [55, 100]
        ("PU50to100", [50, 100]), # [55, 100]
    ])

    #JetShapesAndPUSs_Denom_forRatioPlot = ['Default', "RawPUS", "nVtx #in [1, 25]"] # dict key names used to refer histogram to be used as a denominator histogram in ratio plot
    JetShapesAndPUSs_Denom_forRatioPlot = ["PU1to25"] # dict key names used to refer histogram to be used as a denominator histogram in ratio plot

    JetShapesAndPUSs_Denom_forRatioPlot_forTrgRates = None
    if JetShapesAndPUSs_histoName_TrgRates == JetShapesAndPUSs_histoName_TrgRates_1:
        JetShapesAndPUSs_Denom_forRatioPlot_forTrgRates = ['Default', "RawPUS"]
    elif JetShapesAndPUSs_histoName_TrgRates == JetShapesAndPUSs_histoName_TrgRates_2:
        JetShapesAndPUSs_Denom_forRatioPlot_forTrgRates = ['9x9', "RawPUS_phiRingDefault"]
    elif JetShapesAndPUSs_histoName_TrgRates == JetShapesAndPUSs_histoName_TrgRates_3:
        JetShapesAndPUSs_Denom_forRatioPlot_forTrgRates = ['Default', "RawPUS"]
        
    
    sOpL1JetEnergySF = "L1TJetEtSF.root"




    L1TauMatchedToPFJetVersions = OrderedDict([
        #('DefaultL1JetEt',  "h_jet_byHand_res_vs_iEta_vs_nVtx_RawPUS_HBEF_%s_0" % (resPtBin)),
        #('DefaultL1TauEt', "h_l1TauMatchingPFJet_res_vs_iEta_vs_nVtx_Et_HBEF_%s_0" % (resPtBin)),
        #('DefaultL1TauRawEt', "h_l1TauMatchingPFJet_res_vs_iEta_vs_nVtx_RawEt_HBEF_%s_0" % (resPtBin)),
        ##
        ('DefaultL1Jet Et',        "h_jet_byHand_res_vs_iEta_vs_nVtx_L1JDefault_HBEF_%s_0" % (resPtBin)),
        ('DefaultL1Jet RAWPUS',    "h_jet_byHand_res_vs_iEta_vs_nVtx_RawPUS_HBEF_%s_0" % (resPtBin)),
        ('DefaultL1Tau Et',        "h_jet_byHand_res_vs_iEta_vs_nVtx_L1TauDefault_Et_HBEF_%s_0" % (resPtBin)),
        ('DefaultL1Tau RawEt',     "h_jet_byHand_res_vs_iEta_vs_nVtx_L1TauDefault_RawEt_HBEF_%s_0" % (resPtBin)),
    ])    




    

    
    run_list =  [
        #'run_res1DPlots',
        #'run_res_vs_ieta',
        #'run_effi_periEta',
        #'run_res_vs_jetShape',
        #'run_res1DPlotsForDiffJetShapes',
        #'run_L1JetPt_vs_DefaultL1JetPt'
        #'run_AvgSF_vs_nPV',
        'run_res_vs_jetShapeAndPUS'
        
        #'run_TrgRates'
        #'run_TrgEffi'

        #'run_res_vs_tauJets'
    ]

    
    useAbsEtaBins = True
    ETA_Bins = []
    for iEta in range(-41,42):
        if iEta in [-29, 0, 29]:        continue;
        if useAbsEtaBins and iEta < 0:  continue;
        ETA_Bins.append(str(iEta))
    #ETA_Bins.append('all')
    ETA_Bins.append('HBEF')
    
    IETA_CAT = {}
    IETA_CAT['HBEF'] = [ 1, 41]  ## Whole detector, 1 - 41
    IETA_CAT['HB']   = [ 1, 16]  ## Trigger towers  1 - 16
    IETA_CAT['HE1']  = [17, 20]  ## Trigger towers 17 - 20
    IETA_CAT['HE2a'] = [21, 25]  ## Trigger towers 21 - 25
    IETA_CAT['HE2b'] = [26, 28]  ## Trigger towers 26 - 28
    IETA_CAT['HF']   = [30, 41]  ## Trigger towers 30 - 41

    canvasDims_iEta = [1200,600]
    if useAbsEtaBins: canvasDims_iEta = [700,500]
    
    R.gStyle.SetPadTopMargin(0.055);
    R.gStyle.SetOptTitle(0)
    
    R.gStyle.SetPadTopMargin(0.055);
    R.gStyle.SetPadRightMargin(0.05);
    R.gStyle.SetPadBottomMargin(0.12);
    R.gStyle.SetPadLeftMargin(0.12);
    
    #// use large Times-Roman fonts
    R.gStyle.SetTextFont(132);
    R.gStyle.SetTextSize(0.05);
    
    R.gStyle.SetLabelFont(132,"x");
    R.gStyle.SetLabelFont(132,"y");
    R.gStyle.SetLabelFont(132,"z");
    
    R.gStyle.SetLabelSize(0.045,"x");
    R.gStyle.SetLabelSize(0.045,"y");
    R.gStyle.SetLabelSize(0.045,"z");
    
    R.gStyle.SetTitleSize(0.05,"x");
    R.gStyle.SetTitleSize(0.05,"y");
    R.gStyle.SetTitleSize(0.05,"z");
    
    R.gStyle.SetTitleOffset(0.85,"x");
    R.gStyle.SetTitleOffset(1.0,"y");
    
    #R.gStyle.SetNdivisions(520, "x");
    R.gStyle.SetNdivisions(505, "y");
    
    #// legend attributes
    R.gStyle.SetLegendFillColor(-1);
    R.gStyle.SetLegendBorderSize(1);
    
    R.gStyle.SetOptTitle(0);
    R.gStyle.SetOptStat(0);

    R.gStyle.SetLegendFillColor(0);

    
    if not os.path.exists("./plots"):
        os.mkdir("plots")
    
    #colors   = [R.kBlack, R.kRed, R.kBlue, R.kMagenta, R.kCyan+1, R.kOrange+2, R.kGray+22, R.kGreen+2, R.kTeal, R.kViolet, R.kSpring, R.kBlack, ]
    #markers  = [24, 23, 20]
    colors   = [R.kRed, R.kBlue, R.kMagenta, R.kCyan+1, R.kOrange+2, R.kGreen+2, R.kGray+22, R.kTeal, R.kViolet, R.kSpring, R.kBlack, ]
    markers  = [20, 24, 23, 20]
    colors2  = [R.kBlack, R.kRed, R.kBlue, R.kMagenta, R.kCyan+1, R.kOrange+2, R.kGray+2, R.kGreen+2, R.kBlue-10, R.kViolet, R.kTeal, R.kSpring, R.kBlack, ]
    markers2 = [20, 22, 24, 23, 21]
    
    
    fGaus = R.TF1("g1", "gaus");
    fGaus.SetRange(-1.4, 2.5)
    

    inFiles = []
    for iInFile in range(len(in_fileNames)):
        inFileName  = in_fileNames[iInFile]
        sInFileTag  = sInFileTags[iInFile]            
        
        print "inFileName: {}".format(inFileName); sys.stdout.flush()
        fIn = R.TFile(inFileName)
        if not fIn.IsOpen():
            print "Input file %s couldn't open." % (inFileName)
            return
        inFiles.append(fIn)
    


    if 'run_res1DPlots' in run_list:
        # plot resolution 1D for each iEta
        c5 = R.TCanvas("c5","c5",600,500)
        for jetShape in JetShapes:
            # JetShape = "" plots are with the first version of code for 9x9 jets
            jetShape1 = jetShape
            if jetShape == 'Default':  jetShape1 = ""
            else:                      jetShape1 = "_%s" % (jetShape)
            print "  jetShape1: {}".format(jetShape1)
        
            for plot1DName in ['res', 'PU', 'PUByRawPt']: # ['jet_byHand_res_vs_iEta', 'jet_byHand_PU_vs_iEta', 'jet_byHand_PUByRawPt_vs_iEta']
                for l1Mode in l1Modes:

                    hres2D = {}
                    hres1D = {} 
                    for iInFile in range(len(in_fileNames)):
                        fIn         = inFiles[iInFile]
                        sInFileTag  = sInFileTags[iInFile]

                        hres2D[iInFile] = {}
                        hres1D[iInFile] = {}
                        for iResType in range(len(resTypes)):
                            resType = resTypes[iResType]
                            resType = resType.replace('res', plot1DName)
                            resType = resType.replace('$JETSHAPE', jetShape1)
                            res_histoName = "%s_%s" % (resType, l1Mode)
                            print "res histo: %s" % (res_histoName); sys.stdout.flush()
                            h2D = fIn.Get(res_histoName)
                            hres2D[iInFile][res_histoName] = h2D
                            #hres2D.append(fIn.Get(res_histoName))
                            sNameTmp = "%s_InFile%d" % (res_histoName, iInFile)
                            hres2D[iInFile][res_histoName].SetNameTitle(sNameTmp, sNameTmp)
                            #hres2D[-1].GetXaxis().SetTitle("iEta")
                            #print "res_histoName new: {}".format(hres2D[-1].GetName()); sys.stdout.flush()

                            hres1D[iInFile][res_histoName] = {}
                            for idxEtaBin in range(0, h2D.GetNbinsX()+1):
                                Eta = None
                                hres = None
                                if idxEtaBin == 0: # inclusive iEta
                                    Eta = "HBEF"
                                    hres = h2D.ProjectionY("%s_ProjY_%s" % (h2D.GetName(), Eta), 1, h2D.GetNbinsX())
                                else:
                                    Eta = str( int( h2D.GetXaxis().GetBinCenter(idxEtaBin) ) )
                                    hres = h2D.ProjectionY("%s_ProjY_%s" % (h2D.GetName(), Eta), idxEtaBin, idxEtaBin)

                                if plot1DName == 'PUByRawPt':
                                    hres.Rebin(3)

                                hres1D[iInFile][res_histoName][Eta] = hres


                        resType0 = resTypes[0]
                        resType0 = resType0.replace('res', plot1DName)
                        resType0 = resType0.replace('$JETSHAPE', jetShape1)
                        list_etaBins = hres1D[iInFile]["%s_%s" % (resType0, l1Mode)].keys()
                        #print "list_etaBins: {}".format(list_etaBins)
                        for eta in list_etaBins:
                            c5.Clear()

                            sJetShape = "JetShape 9x9 default"
                            if jetShape1: sJetShape = "JetShape %s" % (jetShape1)
                            legend = R.TLegend(0.3,0.8,1,1)
                            legend.SetNColumns(2)
                            if not useAbsEtaBins:
                                legend.SetHeader("i#eta = %s, %s, %s, %s" % (eta, l1Mode,sInFileTags[iInFile], sJetShape ), "C")
                            else:
                                legend.SetHeader("|i#eta| = %s, %s, %s, %s" % (eta, l1Mode,sInFileTags[iInFile], sJetShape ), "C")

                            # make list of histograms that will be plotted on the same canvas
                            # This is doen to set histo-y-range
                            h_group = []
                            for iResType in range(len(resTypes)):
                                resType = resTypes[iResType]
                                resType = resType.replace('res', plot1DName)
                                resType = resType.replace('$JETSHAPE', jetShape1)
                                res_histoName = "%s_%s" % (resType, l1Mode)

                                # Don't consider PU_RAW histograms to decide yaxis range as PU_RAW histograms are filled with '0' entries  
                                if plot1DName != 'res' and '_Raw_' in resType:
                                    print "Skip h_group.append: plot1DName {}, resType {}".format(plot1DName, resType)
                                    continue

                                h = hres1D[iInFile][res_histoName][eta]
                                h_group.append(h)
                                print "h_append: iInFile {}, res_histoName {}, eta {}".format(iInFile, res_histoName, eta)

                            for iResType in range(len(resTypes)):
                                resType = resTypes[iResType]
                                resType = resType.replace('res', plot1DName)
                                resType = resType.replace('$JETSHAPE', jetShape1)
                                res_histoName = "%s_%s" % (resType, l1Mode)
                                h = hres1D[iInFile][res_histoName][eta]

                                xaxisName = ""
                                plotName = ""
                                if plot1DName == 'res':
                                    xaxisName = "(p_{T}^{L1 jet} - p_{T}^{PF jet}) / p_{T}^{PF jet}"
                                    plotName = "Res1D"
                                elif plot1DName == 'PU':
                                    xaxisName = "L1 jet PU E_{T} [GeV]"
                                    plotName = "PU1D"
                                elif plot1DName == 'PUByRawPt':
                                    xaxisName = "(L1 jet PU E_{T})/(L1 jet raw E_{T})"
                                    plotName = "PUByRawPt1D"

                                h.GetXaxis().SetTitle(xaxisName)
                                h.GetYaxis().SetTitle("Entries")
                                #h.GetYaxis().SetRangeUser(0, h.GetMaximum() * 1.4)
                                h.GetXaxis().SetTitleOffset(1.3)
                                h.SetLineColor(colors[iResType])
                                h.SetMarkerColor(colors[iResType])
                                h.SetMarkerStyle(markers[iInFile])
                                h.SetMarkerSize(0.7)
                                c5.cd()
                                if iResType == 0:
                                    ymin, ymax = getHists1DYRange(h_group)
                                    print "ymin {}, ymax {}, h_group ({}) : {}".format(ymin, ymax, len(h_group), h_group)
                                    h.GetYaxis().SetRangeUser(0, ymax * 1.3)
                                    h.Draw()
                                else:
                                    h.Draw("same")

                                legend.AddEntry(h, sResLegend[iResType], "lep")

                            c5.cd()
                            legend.Draw()
                            c5.Update()
                            sDir1 = "plots%s/JetShape%s/%s_wPFJet%s/%s_%s" % (sInFileVersion, jetShape1, plot1DName,resPtBin, plotName, sInFileTags[iInFile])
                            if not os.path.exists(sDir1):
                                os.makedirs(sDir1)                                
                            c5.SaveAs("%s/jetByHand_%s_%s_%s_%s_%s.png" % (sDir1, l1Mode, plotName, eta, resPtBin, sInFileTags[iInFile]))


            c5.Update()            
            #sTmp = input("enter something")



    if 'run_res_vs_ieta' in run_list:
        c1 = R.TCanvas("c1","c1",500,400)
        c2 = R.TCanvas("c2","c2",canvasDims_iEta[0],canvasDims_iEta[1])
        c3 = R.TCanvas("c3","c3",canvasDims_iEta[0],canvasDims_iEta[1])
        c4 = R.TCanvas("c4","c4",canvasDims_iEta[0],canvasDims_iEta[1])

        #for jetShape in [""] + JetShapes:
        for jetShape in JetShapes:
            # JetShape = "" plots are with the first version of code for 9x9 jets
            jetShape1 = jetShape
            if jetShape == 'Default':  jetShape1 = ""
            else:                      jetShape1 = "_%s" % (jetShape)
            print "  jetShape1: {}".format(jetShape1)
        
            #for plot1DName in ['res', 'PU', 'PUByRawPt']: # ['jet_byHand_res_vs_iEta', 'jet_byHand_PU_vs_iEta', 'jet_byHand_PUByRawPt_vs_iEta']
            for plot1DName in ['res']: # ['jet_byHand_res_vs_iEta', 'jet_byHand_PU_vs_iEta', 'jet_byHand_PUByRawPt_vs_iEta']
                plotName = ""
                if   plot1DName == 'res':       plotName = "Res"
                elif plot1DName == 'PU':        plotName = "PU"
                elif plot1DName == 'PUByRawPt': plotName = "PUByRawPt"

                for l1Mode in l1Modes:
                    legend = R.TLegend(0.3,0.8,1,1)
                    legend.SetNColumns(2)

                    hDummy = R.TH1D("hDummy", "", 1,0,1)

                    hTmps = []
                    #hres2D = []
                    #hres2D_1 = []
                    #hres2D_2 = []
                    h1D_1_OD = OrderedDict()
                    h1D_2_OD = OrderedDict()
                    h1D_3_OD = OrderedDict()
                    hRatio1D_1_OD = OrderedDict()
                    hRatio1D_2_OD = OrderedDict()
                    hRatio1D_3_OD = OrderedDict()
                    for iInFile in range(len(in_fileNames)):
                        sInFileTag  = sInFileTags[iInFile]
                        hTmp = R.TH1D("hTmp%d"%iInFile, "", 1,0,1)
                        hTmp.SetMarkerStyle(markers[iInFile])
                        hTmps.append(hTmp)
                        legend.AddEntry(hTmps[iInFile], sInFileTag, "P")

                    for iInFile in range(len(in_fileNames)):
                        fIn         = inFiles[iInFile]
                        sInFileTag  = sInFileTags[iInFile]                 

                        iHisto = 0
                        h1D_1_OD[sInFileTag] = OrderedDict()
                        h1D_2_OD[sInFileTag] = OrderedDict()
                        h1D_3_OD[sInFileTag] = OrderedDict() 
                        shRatio_Denom = None
                        hRatio1D_1_OD[sInFileTag] = OrderedDict()
                        hRatio1D_2_OD[sInFileTag] = OrderedDict()
                        hRatio1D_3_OD[sInFileTag] = OrderedDict()                         
                        for iResType in range(len(resTypes)):
                            resType = resTypes[iResType]
                            resType = resType.replace('res', plot1DName)
                            resType = resType.replace('$JETSHAPE', jetShape1)
                            res_histoName = "%s_%s" % (resType, l1Mode)

                            # Don't consider PU_RAW histograms to decide yaxis range as PU_RAW histograms are filled with '0' entries 
                            if   plot1DName != 'res' and '_Raw_' in resType: continue 

                            print "res histo: %s" % (res_histoName); sys.stdout.flush()
                            hres2D = fIn.Get(res_histoName)
                            sNameTmp = "%s_InFile%d" % (res_histoName, iInFile)
                            hres2D.SetNameTitle(sNameTmp, sNameTmp)
                            if not useAbsEtaBins:
                                hres2D.GetXaxis().SetTitle("i#eta")
                            else:
                                hres2D.GetXaxis().SetTitle("|i#eta|")
                            print "res_histoName new: {}".format(hres2D.GetName()); sys.stdout.flush()

                            yaxisName = ""
                            if   plot1DName == 'res':
                                yaxisName_mean  = "#mu(L1T/PF - 1)"
                                yaxisName_sigma = "#sigma(L1T/PF - 1)"
                                yaxisName_sigmaByMean = "#frac{#sigma(L1T/PF - 1)}{#mu(L1T/PF - 1) + 1}"
                                fGaus.SetRange(-1.4, 2.5)
                                axisRange_mean  = [-0.8, 1]
                                axisRange_sigma = [0.1, 0.9]
                            elif plot1DName == 'PU':
                                yaxisName_mean  = "#mu(L1 jet PU E_{T}) [GeV]"
                                yaxisName_sigma = "#sigma(L1 jet PU E_{T}) [GeV]"
                                fGaus.SetRange(0, 60)                        
                                axisRange_mean  = [0, 60]
                                axisRange_sigma = [0, 20]
                            elif plot1DName == 'PUByRawPt':
                                yaxisName_mean  = "#mu(#frac{L1 jet PU E_{T}}{L1 jet raw E_{T}}) [GeV]"
                                yaxisName_sigma = "#sigma(#frac{L1 jet PU E_{T}}{L1 jet raw E_{T}}) [GeV]"
                                fGaus.SetRange(0.0, 1.3)
                                axisRange_mean  = [0, 0.9]
                                axisRange_sigma = [0., 0.4]
                                hres2D.RebinY(2)

                            c1.cd()
                            hres2D.Draw("colz")
                            
                            c2.cd()
                            hres2D.FitSlicesY(fGaus, 1,hres2D.GetNbinsX(), 1)
                            
                            
                            c2.cd()
                            #hres2D_1 = R.gDirectory.Get("%s_1" % (res_histoName))
                            hres2D_1 = R.gDirectory.Get("%s_1" % (sNameTmp))
                            hres2D_1.SetNameTitle("%s_1" % (sNameTmp), "%s_1" % (sNameTmp))
                            hres2D_1.SetLineColor(colors[iResType])
                            hres2D_1.SetMarkerColor(colors[iResType])
                            hres2D_1.SetMarkerStyle(markers[iInFile])
                            hres2D_1.SetMarkerSize(0.7)
                            hres2D_1.GetYaxis().SetTitle(yaxisName_mean)
                            #hres2D_1[-1].GetYaxis().SetRangeUser(-0.6,0.8)
                            '''
                            if iHisto == 0 and iInFile == 0:
                                ymin, ymax = getHists1DYRange([hres2D_1])
                                #hres2D_1[-1].GetYaxis().SetRangeUser(ymin, ymax * 1.2)
                                hres2D_1.GetYaxis().SetRangeUser(axisRange_mean[0], axisRange_mean[1])
                                hres2D_1.Draw("PE")
                                c2.Update()
                                c2.SaveAs("tmp.png")
                            else:
                                hres2D_1.Draw("PE same")
                            '''
                            
                            c3.cd()
                            #hres2D_2 = R.gDirectory.Get("%s_2" % (res_histoName))
                            hres2D_2 = R.gDirectory.Get("%s_2" % (sNameTmp)) 
                            hres2D_2.SetNameTitle("%s_2" % (sNameTmp), "%s_2" % (sNameTmp))
                            hres2D_2.SetLineColor(colors[iResType])
                            hres2D_2.SetMarkerColor(colors[iResType])
                            hres2D_2.SetMarkerStyle(markers[iInFile])
                            hres2D_2.SetMarkerSize(0.7)
                            hres2D_2.GetYaxis().SetTitle(yaxisName_sigma)
                            #hres2D_2[-1].GetYaxis().SetRangeUser(0.1,0.7)
                            '''
                            if iHisto == 0 and iInFile == 0:
                                ymin, ymax = getHists1DYRange([hres2D_2])
                                #hres2D_2[-1].GetYaxis().SetRangeUser(ymin, ymax * 1.2)
                                hres2D_2.GetYaxis().SetRangeUser(axisRange_sigma[0], axisRange_sigma[1])                        
                                hres2D_2.Draw("PE")
                            else:
                                hres2D_2.Draw("PE same")
                            '''

                            if iInFile == 0:
                                legend.AddEntry(hres2D_2, sResLegend[iResType], "le")
                                #legend.AddEntry(hDummy, sResLegend[iResType], "le")
                                print "    legend {}: {}, histo: {}".format(iResType, sResLegend[iResType], hres2D_2.GetName())
                            
                            hres2D_MeanPlus1 = hres2D_1.Clone("%s_MeanPlus1" % (sNameTmp))
                            hres2D_MeanPlus1.SetNameTitle("%s_MeanPlus1" % (sNameTmp), "%s_MeanPlus1" % (sNameTmp))
                            for iBin in range(1, hres2D_MeanPlus1.GetNbinsX()+1):
                                binContent = hres2D_MeanPlus1.GetBinContent(iBin)
                                hres2D_MeanPlus1.SetBinContent(iBin, binContent + 1)
                            hEffectiveResolution = hres2D_2.Clone("%s_EffectiveResolution" % (hres2D.GetName()))
                            hEffectiveResolution.Divide(hres2D_2, hres2D_MeanPlus1)
                            hEffectiveResolution.GetYaxis().SetTitle(yaxisName_sigmaByMean)
                            
                            
                            h1D_1_OD[sInFileTag][res_histoName] = hres2D_1
                            h1D_2_OD[sInFileTag][res_histoName] = hres2D_2
                            h1D_3_OD[sInFileTag][res_histoName] = hEffectiveResolution
                            if jetShape in ['Default', '']:
                                if sresTypes_hRatio_Denominator_forDefault in res_histoName:
                                    shRatio_Denom = res_histoName
                            else:
                                if sresTypes_hRatio_Denominator_forPhiRing in res_histoName:
                                    shRatio_Denom = res_histoName
                            
                            iHisto += 1

                            h_dict = h1D_1_OD
                            print "h1D_3_OD::"
                            for sInFileTag in h_dict.keys():
                                print "    h1D_3_OD[%s]" % (sInFileTag)
                                for res_histoName in h_dict[sInFileTag].keys():
                                    h = h_dict[sInFileTag][res_histoName]
                                    print "    h1D_3_OD[%s][%s]: %s, %g" % \
                                        (sInFileTag, res_histoName, h.GetName(),
                                         h.GetBinContent(5) )
                                    
                        print "shRatio_Denom: {}".format(shRatio_Denom)
                        for iResType in range(len(resTypes)):
                            resType = resTypes[iResType]
                            resType = resType.replace('res', plot1DName)
                            resType = resType.replace('$JETSHAPE', jetShape1)
                            res_histoName = "%s_%s" % (resType, l1Mode)
                            
                            # Don't consider PU_RAW histograms to decide yaxis range as PU_RAW histograms are filled with '0' entries 
                            if   plot1DName != 'res' and '_Raw_' in resType: continue 
                            
                            if shRatio_Denom and shRatio_Denom != res_histoName:
                                hRatio1D_1_OD[sInFileTag][res_histoName] = h1D_1_OD[sInFileTag][res_histoName].Clone("%s_ratio" % (h1D_1_OD[sInFileTag][res_histoName].GetName()))
                                hRatio1D_2_OD[sInFileTag][res_histoName] = h1D_2_OD[sInFileTag][res_histoName].Clone("%s_ratio" % (h1D_2_OD[sInFileTag][res_histoName].GetName()))
                                hRatio1D_3_OD[sInFileTag][res_histoName] = h1D_3_OD[sInFileTag][res_histoName].Clone("%s_ratio" % (h1D_3_OD[sInFileTag][res_histoName].GetName()))
                                #
                                hRatio1D_1_OD[sInFileTag][res_histoName].Divide( h1D_1_OD[sInFileTag][res_histoName], h1D_1_OD[sInFileTag][shRatio_Denom] )
                                hRatio1D_2_OD[sInFileTag][res_histoName].Divide( h1D_2_OD[sInFileTag][res_histoName], h1D_2_OD[sInFileTag][shRatio_Denom] )
                                hRatio1D_3_OD[sInFileTag][res_histoName].Divide( h1D_3_OD[sInFileTag][res_histoName], h1D_3_OD[sInFileTag][shRatio_Denom] )
                                #
                                #sRatio = "#frac{h(i)}{h(RAWPUS)}"
                                sRatio = None
                                if jetShape in ['Default', '']:
                                    sRatio = "#frac{h(i)}{h(RAWPUS)}"
                                else:
                                    sRatio = "#frac{h(i)}{h(%s)}" % (sresTypes_hRatio_Denominator_forPhiRing[1:-1])
                                hRatio1D_1_OD[sInFileTag][res_histoName].GetYaxis().SetTitle(sRatio)
                                hRatio1D_2_OD[sInFileTag][res_histoName].GetYaxis().SetTitle(sRatio)
                                hRatio1D_3_OD[sInFileTag][res_histoName].GetYaxis().SetTitle(sRatio)
                                
                    '''            
                    c1.Update()
                    c2.cd()
                    legend.Draw()
                    c2.Update()
                    sDir1 = "plots%s/JetShape%s/%s_wPFJet%s" % (sInFileVersion, jetShape1, plot1DName,resPtBin)
                    if not os.path.exists(sDir1):
                        os.makedirs(sDir1)                                
                    
                    c2.SaveAs("%s/jetByHand_%s_%s_%s_mean_1.png" % (sDir1, l1Mode, resPtBin, plotName))
                    c3.SaveAs("%s/jetByHand_%s_%s_%s_sigma_1.png" % (sDir1, l1Mode, resPtBin, plotName))
                    '''
                    h_dict = h1D_1_OD
                    print "h_dict:: h1D_1_OD - 2"
                    for sInFileTag in h_dict.keys():
                        print "    h_dict[%s]" % (sInFileTag)
                        for res_histoName in h_dict[sInFileTag].keys():
                            h = h_dict[sInFileTag][res_histoName]
                            print "    h_dict[%s][%s]: %s, %g" % \
                                (sInFileTag, res_histoName, h.GetName(),
                                 h.GetBinContent(5) )
                    
                    h1D_1_list = []
                    h1D_2_list = []
                    h1D_3_list = []
                    hRatio1D_1_list = []
                    hRatio1D_2_list = []
                    hRatio1D_3_list = []
                    for iInFile in range(len(in_fileNames)):
                        fIn         = inFiles[iInFile]
                        sInFileTag  = sInFileTags[iInFile]                 
                        
                        for iResType in range(len(resTypes)):
                            resType = resTypes[iResType]
                            resType = resType.replace('res', plot1DName)
                            resType = resType.replace('$JETSHAPE', jetShape1)
                            res_histoName = "%s_%s" % (resType, l1Mode)
                            
                            # Don't consider PU_RAW histograms to decide yaxis range as PU_RAW histograms are filled with '0' entries 
                            if   plot1DName != 'res' and '_Raw_' in resType: continue 
                            
                            h1D_1_list.append( h1D_1_OD[sInFileTag][res_histoName] )
                            h1D_2_list.append( h1D_2_OD[sInFileTag][res_histoName] )
                            h1D_3_list.append( h1D_3_OD[sInFileTag][res_histoName] )                            
                            if shRatio_Denom and shRatio_Denom != res_histoName:
                                hRatio1D_1_list.append( hRatio1D_1_OD[sInFileTag][res_histoName] )
                                hRatio1D_2_list.append( hRatio1D_2_OD[sInFileTag][res_histoName] )
                                hRatio1D_3_list.append( hRatio1D_3_OD[sInFileTag][res_histoName] )


                    h_dict = h1D_1_OD
                    print "h_dict:: h1D_1_OD - 3"
                    for sInFileTag in h_dict.keys():
                        print "    h_dict[%s]" % (sInFileTag)
                        for res_histoName in h_dict[sInFileTag].keys():
                            h = h_dict[sInFileTag][res_histoName]
                            print "    h_dict[%s][%s]: %s, %g" % \
                                (sInFileTag, res_histoName, h.GetName(),
                                 h.GetBinContent(5) )
                            
                            
                    sDir1 = "plots%s/JetShape%s/%s_wPFJet%s" % (sInFileVersion, jetShape1, plot1DName,resPtBin)
                    if not os.path.exists(sDir1):
                        os.makedirs(sDir1)
                    c2 = plotHistos1DAndRatioPlot(h1D_1_list, hRatio1D_1_list, legend, c2, sSaveAs="%s/jetByHand_%s_%s_%s_mean.png" % (sDir1, l1Mode, resPtBin, plotName) )
                    c3 = plotHistos1DAndRatioPlot(h1D_2_list, hRatio1D_2_list, legend, c3, sSaveAs="%s/jetByHand_%s_%s_%s_sigma.png" % (sDir1, l1Mode, resPtBin, plotName) )
                    c4 = plotHistos1DAndRatioPlot(h1D_3_list, hRatio1D_3_list, legend, c4, sSaveAs="%s/jetByHand_%s_%s_%s_effectiveResolution.png" % (sDir1, l1Mode, resPtBin, plotName) )    
                    
                    '''
                    c1.Update()
                    c2.cd()
                    legend.Draw()
                    c2.Update()
                    

                    c3.cd()
                    legend.Draw()
                    c3.Update()
                    
                    sDir1 = "plots%s/JetShape%s/%s_wPFJet%s" % (sInFileVersion, jetShape1, plot1DName,resPtBin)
                    if not os.path.exists(sDir1):
                        os.makedirs(sDir1)                                
                    
                    c2.SaveAs("%s/jetByHand_%s_%s_%s_mean.png" % (sDir1, l1Mode, resPtBin, plotName))
                    c3.SaveAs("%s/jetByHand_%s_%s_%s_sigma.png" % (sDir1, l1Mode, resPtBin, plotName))
                    print "here14 "; sys.stdout.flush()
                    '''


    
    
    
    if 'run_effi_periEta' in run_list:
        for jetShape in JetShapes:
            # JetShape = "" plots are with the first version of code for 9x9 jets
            jetShape1 = jetShape
            if jetShape == 'Default':  jetShape1 = ""
            else:                      jetShape1 = "_%s" % (jetShape)
            print "  jetShape1: {}".format(jetShape1)
        
            c4 = R.TCanvas("c4","c4",650,500)

            # efficiency
            for l1Mode in l1Modes:

                for iEta in range(-41, 42):
                #for iEta in range(1, 2):
                    if iEta in [-29, 29, 0]: continue
                    siEta = str(iEta)

                    for iPt in ['lowPt', 'medPt', 'hiPt']:
                    #for iPt in ['medPt']:
                        #c4.Clear()

                        legend = R.TLegend(0.3,0.8,1,1)
                        legend.SetNColumns(2)

                        #hDummy = R.TH1D("hDummy_" % (l1Mode,), "hDummy", 1,0,1)


                        hTmps = []
                        heffi1D = []
                        for iInFile in range(len(in_fileNames)):
                            sInFileTag  = sInFileTags[iInFile]
                            hTmp = R.TH1D("hTmp1%d"%iInFile, "", 1,0,1)
                            hTmp.SetMarkerStyle(markers[iInFile])
                            hTmps.append(hTmp)
                            legend.AddEntry(hTmps[iInFile], sInFileTag, "P")

                        for iInFile in range(len(in_fileNames)):
                            fIn         = inFiles[iInFile]
                            sInFileTag  = sInFileTags[iInFile]                 


                            for iEffiType in range(len(effiTypes)):
                                effiType = effiTypes[iEffiType]
                                effiType = effiType.replace('$JETSHAPE', jetShape1)
                                effi_histoName = "%s_%s" % (effiType, l1Mode)
                                effi_histoName = effi_histoName.replace("IETA", str(iEta))
                                effi_histoName = effi_histoName.replace("PTBIN", iPt)
                                print "effi histo: %s" % (effi_histoName); sys.stdout.flush()
                                heffi1D.append( fIn.Get(effi_histoName) )
                                sNameTmp = "%s_InFile%d" % (effi_histoName, iInFile)
                                heffi1D[-1].SetTitle("%s;E_{T}(offline jet) [GeV]; Efficiency" % (sNameTmp))

                                c4.cd()
                                heffi1D[-1].SetLineColor(colors[iEffiType])
                                heffi1D[-1].SetLineWidth(1)
                                heffi1D[-1].SetMarkerColor(colors[iEffiType])
                                heffi1D[-1].SetMarkerStyle(markers[iInFile])
                                heffi1D[-1].SetMarkerSize(0.7)

                                if iEffiType == 0 and iInFile == 0:
                                    hTmp = heffi1D[-1].GetCopyTotalHisto()
                                    hTmp.GetYaxis().SetRangeUser(0, 1.3)
                                    hTmp.Draw("AXIS")
                                    heffi1D[-1].Draw("same")
                                else:                                heffi1D[-1].Draw(" same")

                                if iInFile == 0:
                                    legend.AddEntry(heffi1D[-1], sResLegend[iEffiType], "le")

                        c4.cd()
                        legend.Draw()
                        c4.Update()

                        sDir1 = "plots%s/JetShape%s/effi_wPFJet%s/%s" % (sInFileVersion, jetShape1, resPtBin, sInFileTags[iInFile])
                        if not os.path.exists(sDir1):
                            os.makedirs(sDir1)                                
                        
                        c4.SaveAs("%s/jetByHand_%s_%s_%s_Effi.png" % (sDir1, l1Mode,str(iEta),iPt))

            #sTmp = input("enter something")
            print "here15 "; sys.stdout.flush()

    



    
    
    if 'run_res1DPlotsForDiffJetShapes' in run_list:
        # plot resolution 1D for each iEta
        c5 = R.TCanvas("c5","c5",600,500)
        for plot1DName in ['res', 'PU', 'PUByRawPt']: # ['jet_byHand_res_vs_iEta', 'jet_byHand_PU_vs_iEta', 'jet_byHand_PUByRawPt_vs_iEta']
            plotName = ""
            if   plot1DName == 'res':       plotName = "Res"
            elif plot1DName == 'PU':        plotName = "PU"
            elif plot1DName == 'PUByRawPt': plotName = "PUByRawPt"

            for l1Mode in l1Modes:
                hDummy = R.TH1D("hDummy", "", 1,0,1)
        
                hres1D_dict    = {}
                for iInFile in range(len(in_fileNames)):
                    fIn         = inFiles[iInFile]
                    sInFileTag  = sInFileTags[iInFile]                 
                   
                    hres1D_dict[sInFileTag]    = {}
                    iResType = 0
                    for jetShape in JetShapes:
                        # JetShape = "" plots are with the first version of code for 9x9 jets
                        jetShape1 = jetShape
                        if jetShape == 'Default':  jetShape1 = ""
                        else:                      jetShape1 = "_%s" % (jetShape)

                        resType = resType_selected
                        resType = resType.replace('res', plot1DName)
                        resType = resType.replace('$JETSHAPE', jetShape1)
                        res_histoName = "%s_%s" % (resType, l1Mode)

                        # Don't consider PU_RAW histograms to decide yaxis range as PU_RAW histograms are filled with '0' entries 
                        if   plot1DName != 'res' and '_Raw_' in resType: continue 

                        print "res histo: %s" % (res_histoName); sys.stdout.flush()
                        h2D = (fIn.Get(res_histoName))
                        sNameTmp = "%s_InFile%d" % (res_histoName, iInFile)
                        h2D.SetNameTitle(sNameTmp, sNameTmp)
                        h2D.GetXaxis().SetTitle("iEta")
                        print "res_histoName new: {}".format(h2D.GetName()); sys.stdout.flush()

                        yaxisName = ""
                        if   plot1DName == 'res':
                            yaxisName_mean  = "#mu(#frac{online - offline}{offline})"
                            yaxisName_sigma = "#sigma(#frac{online - offline}{offline})"
                            fGaus.SetRange(-1.4, 2.5)
                            axisRange_mean  = [-0.8, 1]
                            axisRange_sigma = [0.1, 0.9]
                        elif plot1DName == 'PU':
                            yaxisName_mean  = "#mu(L1 jet PU E_{T}) [GeV]"
                            yaxisName_sigma = "#sigma(L1 jet PU E_{T}) [GeV]"
                            fGaus.SetRange(0, 60)                        
                            axisRange_mean  = [0, 60]
                            axisRange_sigma = [0, 20]
                        elif plot1DName == 'PUByRawPt':
                            yaxisName_mean  = "#mu(#frac{L1 jet PU E_{T}}{L1 jet raw E_{T}}) [GeV]"
                            yaxisName_sigma = "#sigma(#frac{L1 jet PU E_{T}}{L1 jet raw E_{T}}) [GeV]"
                            fGaus.SetRange(0.0, 1.3)
                            axisRange_mean  = [0, 0.9]
                            axisRange_sigma = [0., 0.4]
                            #h2D.RebinY(2)

                        hres1D_dict[sInFileTag][jetShape] = {}
                        for idxEtaBin in range(1, h2D.GetNbinsX()+1):
                            Eta = int( h2D.GetXaxis().GetBinCenter(idxEtaBin) )

                            hres = h2D.ProjectionY("%s_ProjY_%d" % (h2D.GetName(), Eta), idxEtaBin, idxEtaBin)

                            if plot1DName == 'PUByRawPt':
                                hres.Rebin(3)

                            hres1D_dict[sInFileTag][jetShape][Eta] = hres

                        iResType += 1


                    list_etaBins = hres1D_dict[sInFileTag][jetShape].keys()                    
                    for eta in list_etaBins:
                        c5.Clear()
                        legend = R.TLegend(0.3,0.85,1,1)
                        legend.SetNColumns(3)
                        legend.SetHeader("i#eta = %d, %s, %s" % (eta, l1Mode,sInFileTag), "C")
                        
                        hres1D_list = []
                        for jetShape in JetShapes:
                            hres1D_list.append(hres1D_dict[sInFileTag][jetShape][eta])

                        for iJetShape in range(len(JetShapes)):
                            jetShape = JetShapes[iJetShape]
                            h = hres1D_dict[sInFileTag][jetShape][eta]
                            
                            xaxisName = ""
                            plotName = ""
                            if plot1DName == 'res':
                                xaxisName = "(p_{T}^{L1 jet} - p_{T}^{PF jet}) / p_{T}^{PF jet}"
                                plotName = "Res1D"
                            elif plot1DName == 'PU':
                                xaxisName = "L1 jet PU E_{T} [GeV]"
                                plotName = "PU1D"
                            elif plot1DName == 'PUByRawPt':
                                xaxisName = "(L1 jet PU E_{T})/(L1 jet raw E_{T})"
                                plotName = "PUByRawPt1D"

                            h.GetXaxis().SetTitle(xaxisName)
                            h.GetYaxis().SetTitle("Entries")
                            #h.GetYaxis().SetRangeUser(0, h.GetMaximum() * 1.4)
                            h.GetXaxis().SetTitleOffset(1.3)
                            h.SetLineColor(colors2[iJetShape])
                            h.SetMarkerColor(colors2[iJetShape])
                            h.SetMarkerStyle(markers2[iJetShape % 4])
                            h.SetMarkerSize(0.7)
                            c5.cd()
                            if jetShape == JetShapes[0]:
                                ymin, ymax = getHists1DYRange(hres1D_list)
                                h.GetYaxis().SetRangeUser(0, ymax)
                                h.Draw()
                            else:
                                h.Draw("same")

                            legend.AddEntry(h, jetShape, "lep")

                        c5.cd()
                        legend.Draw()
                        c5.Update()
                        sDir1 = "plots%s/CompareDiffJetShapes/%s/PFJet%s/%s/%s" % (sInFileVersion, plot1DName, resPtBin, sInFileTags[iInFile],plotName)
                        if not os.path.exists(sDir1):
                            os.makedirs(sDir1)                                
                        c5.SaveAs("%s/%s_%s_%d.png" % (sDir1, l1Mode, plotName, eta))






            
    if 'run_res_vs_jetShape' in run_list:
        c1 = R.TCanvas("c1","c1",500,400)
        #c2 = R.TCanvas("c2","c2",1200,600)
        #c3 = R.TCanvas("c3","c3",1200,600)
        c2_1 = R.TCanvas("c2_1","c2_1",canvasDims_iEta[0],canvasDims_iEta[1])
        c3_1 = R.TCanvas("c3_1","c3_1",canvasDims_iEta[0],canvasDims_iEta[1])
        c4_1 = R.TCanvas("c4_1","c4_1",canvasDims_iEta[0],canvasDims_iEta[1])

        for plot1DName in ['res', 'PU', 'PUByRawPt']: # ['jet_byHand_res_vs_iEta', 'jet_byHand_PU_vs_iEta', 'jet_byHand_PUByRawPt_vs_iEta']
            plotName = ""
            if   plot1DName == 'res':       plotName = "Res"
            elif plot1DName == 'PU':        plotName = "PU"
            elif plot1DName == 'PUByRawPt': plotName = "PUByRawPt"

            for l1Mode in l1Modes:
                legend = R.TLegend(0.3,0.8,1,1)
                legend.SetNColumns(2)

                hDummy = R.TH1D("hDummy", "", 1,0,1)
        
                hres2D_dict    = OrderedDict()
                hresMean_dict  = OrderedDict()
                hresSigma_dict = OrderedDict()
                hEffectiveResolution_dict = OrderedDict()
                hRatiores2D_dict    = OrderedDict()
                hRatioresMean_dict  = OrderedDict()
                hRatioresSigma_dict = OrderedDict()
                hRatioEffectiveResolution_dict = OrderedDict()
                for iInFile in range(len(in_fileNames)):
                    fIn         = inFiles[iInFile]
                    sInFileTag  = sInFileTags[iInFile]                 

                    legend = R.TLegend(0.3,0.85,1,1)
                    legend.SetNColumns(3)
                    
                    hres2D_dict[sInFileTag] = OrderedDict()
                    hresMean_dict[sInFileTag] = OrderedDict()
                    hresSigma_dict[sInFileTag] = OrderedDict()
                    hEffectiveResolution_dict[sInFileTag] = OrderedDict()
                    hRatiores2D_dict[sInFileTag] = OrderedDict()
                    hRatioresMean_dict[sInFileTag] = OrderedDict()
                    hRatioresSigma_dict[sInFileTag] = OrderedDict()
                    hRatioEffectiveResolution_dict[sInFileTag] = OrderedDict()
                    shRatio_Denom_JetShape = None
                    iResType = 0
                    for jetShape in JetShapes:
                        # JetShape = "" plots are with the first version of code for 9x9 jets
                        jetShape1 = jetShape
                        if jetShape == 'Default':  jetShape1 = ""
                        else:                      jetShape1 = "_%s" % (jetShape)
                        #print "  jetShape1: {}".format(jetShape1)

                        resType = None
                        if jetShape == 'Default': resType = resType_default
                        else:                     resType = resType_selected
                        resType = resType.replace('res', plot1DName)
                        resType = resType.replace('$JETSHAPE', jetShape1)
                        res_histoName = "%s_%s" % (resType, l1Mode)

                        # Don't consider PU_RAW histograms to decide yaxis range as PU_RAW histograms are filled with '0' entries 
                        if   plot1DName != 'res' and '_Raw_' in resType: continue 

                        print "res histo: %s" % (res_histoName); sys.stdout.flush()
                        h2D = (fIn.Get(res_histoName))
                        sNameTmp = "%s_InFile%d" % (res_histoName, iInFile)
                        h2D.SetNameTitle(sNameTmp, sNameTmp)
                        if not useAbsEtaBins:
                            h2D.GetXaxis().SetTitle("iEta")
                        else:
                            h2D.GetXaxis().SetTitle("|i#eta|")
                        print "res_histoName new: {}".format(h2D.GetName()); sys.stdout.flush()

                        yaxisName = ""
                        if   plot1DName == 'res':
                            yaxisName_mean  = "#mu(L1T/PF - 1)"
                            yaxisName_sigma = "#sigma(L1T/PF - 1)"
                            yaxisName_sigmaByMean = "#frac{#sigma(L1T/PF - 1)}{#mu(L1T/PF - 1) + 1}"
                            fGaus.SetRange(-1.4, 2.5)
                            axisRange_mean  = [-0.8, 1]
                            axisRange_sigma = [0.1, 0.9]
                        elif plot1DName == 'PU':
                            yaxisName_mean  = "#mu(L1 jet PU E_{T}) [GeV]"
                            yaxisName_sigma = "#sigma(L1 jet PU E_{T}) [GeV]"
                            fGaus.SetRange(0, 60)                        
                            axisRange_mean  = [0, 60]
                            axisRange_sigma = [0, 20]
                        elif plot1DName == 'PUByRawPt':
                            yaxisName_mean  = "#mu(#frac{L1 jet PU E_{T}}{L1 jet raw E_{T}}) [GeV]"
                            yaxisName_sigma = "#sigma(#frac{L1 jet PU E_{T}}{L1 jet raw E_{T}}) [GeV]"
                            fGaus.SetRange(0.0, 1.3)
                            axisRange_mean  = [0, 0.9]
                            axisRange_sigma = [0., 0.4]
                            #h2D.RebinY(2)

                        c2_1.cd()
                        h2D.FitSlicesY(fGaus, 1,h2D.GetNbinsX(), 1)
                        

                        c2_1.cd()
                        #hres2D_1 = R.gDirectory.Get("%s_1" % (res_histoName))
                        h2D_1 = (R.gDirectory.Get("%s_1" % (sNameTmp)))
                        h2D_1.SetNameTitle("%s_1" % (sNameTmp), "%s_1" % (sNameTmp))
                        h2D_1.SetLineColor(colors2[iResType])
                        h2D_1.SetMarkerColor(colors2[iResType])
                        #h2D_1.SetMarkerStyle(markers[iInFile])
                        h2D_1.SetMarkerStyle(markers2[iResType % 4])
                        h2D_1.SetMarkerSize(0.7)
                        h2D_1.GetYaxis().SetTitle(yaxisName_mean)

                        c3_1.cd()
                        #hres2D_2 = R.gDirectory.Get("%s_2" % (res_histoName))
                        h2D_2 = ( R.gDirectory.Get("%s_2" % (sNameTmp)) )
                        h2D_2.SetNameTitle("%s_2" % (sNameTmp), "%s_2" % (sNameTmp))
                        h2D_2.SetLineColor(colors2[iResType])
                        h2D_2.SetMarkerColor(colors2[iResType])
                        #h2D_2.SetMarkerStyle(markers[iInFile])
                        h2D_2.SetMarkerStyle(markers2[iResType % 4])
                        h2D_2.SetMarkerSize(0.7)
                        h2D_2.GetYaxis().SetTitle(yaxisName_sigma)

                        c4_1.cd()
                        h2D_MeanPlus1 = h2D_1.Clone("%s_MeanPlus1" % (h2D.GetName()))
                        for iBin in range(1, h2D_MeanPlus1.GetNbinsX()+1):
                            binContent = h2D_MeanPlus1.GetBinContent(iBin)
                            h2D_MeanPlus1.SetBinContent(iBin, binContent + 1)
                        hEffectiveResolution = h2D_2.Clone("%s_EffectiveResolution" % (h2D.GetName()))
                        hEffectiveResolution.Divide(h2D_2, h2D_MeanPlus1)
                        hEffectiveResolution.GetYaxis().SetTitle(yaxisName_sigmaByMean)
                        
                        hres2D_dict[sInFileTag][jetShape] = h2D
                        hresMean_dict[sInFileTag][jetShape] = h2D_1
                        hresSigma_dict[sInFileTag][jetShape] = h2D_2
                        hEffectiveResolution_dict[sInFileTag][jetShape] = hEffectiveResolution
                        if jetShape == 'Default':
                            shRatio_Denom_JetShape = jetShape
                        iResType += 1
                        
                        legend.AddEntry(hresSigma_dict[sInFileTag][jetShape], jetShape, "lep")


                    hresMean_list = []
                    hresSigma_list = []
                    hEffectiveResolution_list = []
                    hRatioresMean_list = []
                    hRatioresSigma_list = []
                    hRatioEffectiveResolution_list = []
                    for jetShape in JetShapes:
                        hresMean_list.append(hresMean_dict[sInFileTag][jetShape])
                        hresSigma_list.append(hresSigma_dict[sInFileTag][jetShape])
                        hEffectiveResolution_list.append(hEffectiveResolution_dict[sInFileTag][jetShape])
                        
                        if jetShape == shRatio_Denom_JetShape: continue
                        
                        hRatioresMean_dict[sInFileTag][jetShape] = hresMean_dict[sInFileTag][jetShape].Clone("%s_Ratio" % hresMean_dict[sInFileTag][jetShape].GetName())
                        hRatioresSigma_dict[sInFileTag][jetShape] = hresSigma_dict[sInFileTag][jetShape].Clone("%s_Ratio" % hresSigma_dict[sInFileTag][jetShape].GetName())
                        hRatioEffectiveResolution_dict[sInFileTag][jetShape] = hEffectiveResolution_dict[sInFileTag][jetShape].Clone("%s_Ratio" % hEffectiveResolution_dict[sInFileTag][jetShape].GetName())
                        
                        hRatioresMean_dict[sInFileTag][jetShape].Divide(hresMean_dict[sInFileTag][jetShape], hresMean_dict[sInFileTag][shRatio_Denom_JetShape])
                        hRatioresSigma_dict[sInFileTag][jetShape].Divide(hresSigma_dict[sInFileTag][jetShape], hresSigma_dict[sInFileTag][shRatio_Denom_JetShape])
                        hRatioEffectiveResolution_dict[sInFileTag][jetShape].Divide(hEffectiveResolution_dict[sInFileTag][jetShape], hEffectiveResolution_dict[sInFileTag][shRatio_Denom_JetShape])
                        
                        sYName = "#frac{h(i)}{h(%s)}" % (shRatio_Denom_JetShape)
                        hRatioresMean_dict[sInFileTag][jetShape].GetYaxis().SetTitle(sYName)
                        hRatioresSigma_dict[sInFileTag][jetShape].GetYaxis().SetTitle(sYName)
                        hRatioEffectiveResolution_dict[sInFileTag][jetShape].GetYaxis().SetTitle(sYName)
                        
                        hRatioresMean_list.append(hRatioresMean_dict[sInFileTag][jetShape])
                        hRatioresSigma_list.append(hRatioresSigma_dict[sInFileTag][jetShape])
                        hRatioEffectiveResolution_list.append(hRatioEffectiveResolution_dict[sInFileTag][jetShape])
                        
                    
                    sDir1 = "plots%s/CompareDiffJetShapes/%s/PFJet%s/%s/%s" % (sInFileVersion, plot1DName, resPtBin, sInFileTags[iInFile],sResLegend_selected)
                    if not os.path.exists(sDir1):
                        os.makedirs(sDir1)                                
                    c2_1 = plotHistos1DAndRatioPlot(hresMean_list, hRatioresMean_list, legend, c2_1, sSaveAs="%s/DiffJetShapes_%s_%s_mean.png" % (sDir1, l1Mode, plotName))
                    c3_1 = plotHistos1DAndRatioPlot(hresSigma_list, hRatioresSigma_list, legend, c3_1, sSaveAs="%s/DiffJetShapes_%s_%s_sigma.png" % (sDir1, l1Mode, plotName))
                    c4_1 = plotHistos1DAndRatioPlot(hEffectiveResolution_list, hRatioEffectiveResolution_list, legend, c4_1, sSaveAs="%s/DiffJetShapes_%s_%s_effectiveResolution.png" % (sDir1, l1Mode, plotName))
                    '''
                    c1.Update()
                    c2_1.Update()
                    c3_1.Update()
                    '''
                    #c2_1.SaveAs("%s/DiffJetShapes_%s_%s_mean.png" % (sDir1, l1Mode, plotName))
                    #c3_1.SaveAs("%s/DiffJetShapes_%s_%s_sigma.png" % (sDir1, l1Mode, plotName))
                    


    if 'run_L1JetPt_vs_DefaultL1JetPt' in run_list:
        sHistoName0 = "h_jet_byHand_L1JetPt_vs_DefaultL1JetPt_RawPUS_$IETA_PtAllBins_0_emu"

        R.gStyle.SetPadRightMargin(0.12);
        c1 = R.TCanvas("c1","c1",600,500)
        #R.gStyle.SetOptStat(111111)
        
        ETA_Bins = ['HBEF']
        for iInFile in range(len(in_fileNames)):
            fIn         = inFiles[iInFile]
            sInFileTag  = sInFileTags[iInFile]
            
            for iEta in ETA_Bins:
                sHistoName = sHistoName0
                sHistoName = sHistoName.replace('$IETA', iEta)
            
                h2D = fIn.Get(sHistoName)
                h2D.GetXaxis().SetTitle("pT(L1 jet before layer2 calibration) [GeV]")
                h2D.GetYaxis().SetTitle("#frac{pT(L1 jet by hand) - pT(default L1 jet)}{pT(default L1 jet)}")
                
                c1.cd()
                R.gPad.SetLogx()
                h2D.Draw("colz")

                #line = R.TLine(0.1,0.5,0.95,0.5)
                line = R.TLine(h2D.GetXaxis().GetXmin(), 0, h2D.GetXaxis().GetXmax(), 0)
                line.SetLineStyle(2)
                #line.DrawLineNDC(0.1,0.5,0.95,0.5)
                line.Draw()
                
                c1.Update()
                sDir1 = "plots%s/ValidateLayer2Calibration/%s" % (sInFileVersion, sInFileTags[iInFile])
                if not os.path.exists(sDir1):
                    os.makedirs(sDir1)
                c1.SaveAs("%s/L1JetPt_New_vs_Default_ieta%s.png" % (sDir1, iEta))

                #input('Enter anything')



    if 'run_AvgSF_vs_nPV' in run_list:
        l1Mode = 'emu'
        sHistoName0_list = [
            "h_nPV_vs_L1JetDefaultRAW_SF_%s_$IETA" % l1Mode,
            'h_nPV_vs_L1JetDefaultPUS_SF_%s_$IETA' % l1Mode
        ]
        sLegends = ["RAW #times 0.5", "RAWPUS"]
        nPV_Median = 54
        
        nPVRanges = {
            "lowPU":  [ 1,  25],
            #"medPU":  [50,  54],
            "highPU": [50, 100],  # [55, 100]
        }
        
        
        c1 = R.TCanvas("c1","c1",500,400)
        
        R.gStyle.SetPadRightMargin(0.12);
        c1 = R.TCanvas("c1","c1",600,500)
        #R.gStyle.SetOptStat(111111)
                
        iInFile = 0
        fIn         = inFiles[iInFile]
        sInFileTag  = sInFileTags[iInFile]
        
        for iEta in ETA_Bins:
            c1.cd()
            legend = R.TLegend(0.1,0.8,1,1)
            legend.SetHeader("|i#eta| = %s, %s, 9x9 default L1T jets " % (iEta, sInFileTags[iInFile] ), "C")
            legend.SetNColumns(2)

            hTmps = []
            for iHisto in range(0, len(sHistoName0_list)):
                hTmp = R.TH1D("hTmp%d"%iHisto, "", 1,0,1)
                hTmp.SetMarkerStyle(markers[iHisto])
                hTmp.SetLineColor(R.kBlack)
                hTmps.append(hTmp)
                legend.AddEntry(hTmps[iHisto], sLegends[iHisto], "lep")
            
            hTmp = R.TH1D("hTmp%d"%iHisto, "", 1,0,1)    
            #legend.AddEntry(hTmp, "", "")    
                
            h_list = []
            for iHisto in range(0, len(sHistoName0_list)):
                sHistoName = sHistoName0_list[iHisto]
                sHistoName = sHistoName.replace('$IETA', iEta)
                #print "sHistoName: %s" % (sHistoName)
                h2D = fIn.Get(sHistoName)
                h2D.GetYaxis().SetTitle("E_{T}(L1T)/E_{T}(PF)")
                
                #kBin_nPV_Median = h2D.GetXaxis().FindBin(nPV_Median)
                #h1 = h2D.ProjectionY("%s_ProjY_1"%(h2D.GetName()), 1,kBin_nPV_Median)
                #h2 = h2D.ProjectionY("%s_ProjY_2"%(h2D.GetName()), kBin_nPV_Median+1,h2D.GetNbinsX())
                kBin_nPVRanges = {}
                kBin_nPVRanges["lowPU"] = [
                    h2D.GetXaxis().FindBin(nPVRanges["lowPU"][0]),
                    h2D.GetXaxis().FindBin(nPVRanges["lowPU"][1])
                ]
                if "medPU" in nPVRanges.keys():
                    kBin_nPVRanges["medPU"] = [
                        h2D.GetXaxis().FindBin(nPVRanges["medPU"][0]),
                        h2D.GetXaxis().FindBin(nPVRanges["medPU"][1])
                    ]
                kBin_nPVRanges["highPU"] = [
                    h2D.GetXaxis().FindBin(nPVRanges["highPU"][0]),
                    h2D.GetXaxis().FindBin(nPVRanges["highPU"][1])
                ]
                h1 = h2 = h3 = None
                h1 = h2D.ProjectionY("%s_ProjY_1"%(h2D.GetName()), kBin_nPVRanges["lowPU"][0],  kBin_nPVRanges["lowPU"][1])
                if "medPU" in nPVRanges.keys():
                    h2 = h2D.ProjectionY("%s_ProjY_2"%(h2D.GetName()), kBin_nPVRanges["medPU"][0],  kBin_nPVRanges["medPU"][1])
                h3 = h2D.ProjectionY("%s_ProjY_3"%(h2D.GetName()), kBin_nPVRanges["highPU"][0], kBin_nPVRanges["highPU"][1])
                #print "iEta {}, sHistoName {}, kBin_nPV_Median {}, ".format(iEta, sHistoName, kBin_nPV_Median)
                print "iEta {}, sHistoName {},   kBin_nPVRanges {}".format(iEta, sHistoName, kBin_nPVRanges)
                h1.Rebin(2)
                if "medPU" in nPVRanges.keys(): h2.Rebin(2)
                h3.Rebin(2)

                h1.GetXaxis().SetRangeUser(0,2.5)
                
                h1.Scale(0.04)
                
                h1.GetYaxis().SetTitle("Entries")
                h1.GetXaxis().SetTitleOffset(1.3)
                h1.SetLineColor(colors[0]) # (colors[iHisto*2])
                h1.SetMarkerColor(colors[0])
                h1.SetMarkerStyle(markers[iHisto])
                h1.SetMarkerSize(0.7)

                if "medPU" in nPVRanges.keys():
                    h2.GetYaxis().SetTitle("Entries")
                    h2.GetXaxis().SetTitleOffset(1.3)
                    h2.SetLineColor(colors[1])
                    h2.SetMarkerColor(colors[1])
                    h2.SetMarkerStyle(markers[iHisto])
                    h2.SetMarkerSize(0.7)

                h3.GetYaxis().SetTitle("Entries")
                h3.GetXaxis().SetTitleOffset(1.3)
                h3.SetLineColor(colors[2])
                h3.SetMarkerColor(colors[2])
                h3.SetMarkerStyle(markers[iHisto])
                h3.SetMarkerSize(0.7)

                fGaus.SetRange(0, 2.5)
                fitResults = {}
                fGaus_1 = R.TF1("g1", "gaus");
                if "medPU" in nPVRanges.keys(): fGaus_2 = R.TF1("g1", "gaus");
                fGaus_3 = R.TF1("g1", "gaus");
                fGaus_1.SetRange(0, 2.5)
                if "medPU" in nPVRanges.keys(): fGaus_2.SetRange(0, 2.5)
                fGaus_3.SetRange(0, 2.5)

                h1.Fit(fGaus_1, "R0")
                if "medPU" in nPVRanges.keys(): h2.Fit(fGaus_2, "R0")
                h3.Fit(fGaus_3, "R0")
                
                
                '''
                legend.AddEntry(h1, "%s, %d #leq nPV #leq %d" % (sLegends[iHisto], nPVRanges["lowPU"][0],  nPVRanges["lowPU"][1]), "lep")
                legend.AddEntry(h2, "%s, %d #leq nPV #leq %d" % (sLegends[iHisto], nPVRanges["medPU"][0],  nPVRanges["medPU"][1]), "lep")
                legend.AddEntry(h3, "%s, %d #leq nPV #leq %d" % (sLegends[iHisto], nPVRanges["highPU"][0], nPVRanges["highPU"][1]), "lep")
                '''
                legend.AddEntry(h1, "nPV#in[%d, %d], #mu=%5.4f, #sigma=%.2e" % (nPVRanges["lowPU"][0],  nPVRanges["lowPU"][1], fGaus_1.GetParameter(1),fGaus_1.GetParameter(2)), "lep")
                if "medPU" in nPVRanges.keys(): legend.AddEntry(h2, "nPV#in[%d, %d], #mu=%5.4f, #sigma=%.2e" % (nPVRanges["medPU"][0],  nPVRanges["medPU"][1], fGaus_2.GetParameter(1),fGaus_2.GetParameter(2)), "lep")
                legend.AddEntry(h3, "nPV#in[%d, %d], #mu=%5.4f, #sigma=%.2e" % (nPVRanges["highPU"][0], nPVRanges["highPU"][1], fGaus_3.GetParameter(1),fGaus_3.GetParameter(2)), "lep")
                
                h_list.append(h1)
                if "medPU" in nPVRanges.keys(): h_list.append(h2)
                h_list.append(h3)
                
            ymin, ymax           = getHists1DYRange(h_list)
            for ih in range(0, len(h_list)):
                h = h_list[ih]
                c1.cd()
                if ih == 0:
                    h.GetYaxis().SetRangeUser(ymin, ymax)
                    h.Draw("PE")
                else:
                    h.Draw("PE same")

            legend.Draw()
            c1.Update()
            
            sDir1 = "plots%s/AvgL1TLayer2SF/wPFJet%s/%s" % (sInFileVersion, resPtBin, sInFileTags[iInFile])
            if not os.path.exists(sDir1):
                os.makedirs(sDir1)                                
            c1.SaveAs("%s/AvgL1TLayer2SF_%s_%s.png" % (sDir1, l1Mode, iEta))
            
            

            
    if 'run_res_vs_jetShapeAndPUS' in run_list:
        c1 = R.TCanvas("c1","c1",500,400)
        #c2 = R.TCanvas("c2","c2",1200,600)
        #c3 = R.TCanvas("c3","c3",1200,600)
        c2_1 = R.TCanvas("c2_1","c2_1",canvasDims_iEta[0],canvasDims_iEta[1])
        c3_1 = R.TCanvas("c3_1","c3_1",canvasDims_iEta[0],canvasDims_iEta[1])
        c4_1 = R.TCanvas("c4_1","c4_1",canvasDims_iEta[0],canvasDims_iEta[1])
        
        
        fOpL1JetEnergySF = R.TFile(sOpL1JetEnergySF, "RECREATE")
        
        
        for JetShape, JetPUSs in JetShapesAndPUSs.items():
            print "JetShape {}, JetPUS {}".format(JetShape, JetPUSs)
            for JetPUSName, JetPUSHistName in JetPUSs.items():
                print "    JetShape {}, JetPUSName {}, JetPUSHistName {}".format(JetShape, JetPUSName, JetPUSHistName)
            
        #exit(0)
        
        for plot1DName in ['res']: #['res', 'PU', 'PUByRawPt']: # ['jet_byHand_res_vs_iEta', 'jet_byHand_PU_vs_iEta', 'jet_byHand_PUByRawPt_vs_iEta']
            plotName = ""
            if   plot1DName == 'res':       plotName = "Res"
            elif plot1DName == 'PU':        plotName = "PU"
            elif plot1DName == 'PUByRawPt': plotName = "PUByRawPt"

            for l1Mode in l1Modes:
                hDummy = R.TH1D("hDummy", "", 1,0,1)
        
                hres2D_dict    = OrderedDict()
                hresMean_dict  = OrderedDict()
                hresSigma_dict = OrderedDict()
                hEffectiveResolution_dict = OrderedDict()
                hRatiores2D_dict    = OrderedDict()
                hRatioresMean_dict  = OrderedDict()
                hRatioresSigma_dict = OrderedDict()
                hRatioEffectiveResolution_dict = OrderedDict()
                for iInFile in range(len(in_fileNames)):
                    fIn         = inFiles[iInFile]
                    sInFileTag  = sInFileTags[iInFile]                 

                    legend = R.TLegend(0.2,0.8,1,1)
                    legend.SetNColumns(2)
                    
                    hres2D_dict[sInFileTag] = OrderedDict()
                    hresMean_dict[sInFileTag] = OrderedDict()
                    hresSigma_dict[sInFileTag] = OrderedDict()
                    hEffectiveResolution_dict[sInFileTag] = OrderedDict()
                    hRatiores2D_dict[sInFileTag] = OrderedDict()
                    hRatioresMean_dict[sInFileTag] = OrderedDict()
                    hRatioresSigma_dict[sInFileTag] = OrderedDict()
                    hRatioEffectiveResolution_dict[sInFileTag] = OrderedDict()
                    shRatio_Denom_JetShape = None
                    iPURange = 0
                    hTmps = []
                    for PURangeName, PURange in nPVRanges.items():
                        hTmp = R.TH1D("hTmp%d"%(iPURange), "", 1,0,1)
                        hTmp.SetMarkerStyle(markers[iPURange])
                        hTmp.SetLineColor(R.kBlack)
                        hTmp.SetMarkerSize(0.7)
                        hTmps.append(hTmp)
                        legend.AddEntry(hTmps[iPURange], PURangeName, "lep")
                        iPURange += 1
                    
                    iResType = 0                    
                    for jetShape, JetPUSs in JetShapesAndPUSs.items():
                        # JetShape = "" plots are with the first version of code for 9x9 jets
                        jetShape1 = jetShape
                        if jetShape == 'Default':  jetShape1 = ""
                        else:                      jetShape1 = "_%s" % (jetShape)
                        #print "  jetShape1: {}".format(jetShape1)
                        
                        for JetPUSName, JetPUSHistName in JetPUSs.items():
                            resType = JetPUSHistName
                            resType = resType.replace('res', plot1DName)
                            resType = resType.replace('$JETSHAPE', jetShape1)
                            res_histoName = "%s_%s" % (resType, l1Mode)
                            print "res histo: %s" % (res_histoName); sys.stdout.flush()
                            h3D = (fIn.Get(res_histoName))
                            #sNameTmp0 = "%s_InFile%d" % (res_histoName, iInFile)
                            sNameTmp0 = "%s_%s" % (res_histoName, sInFileTags[iInFile])
                            h3D.SetNameTitle(sNameTmp0, sNameTmp0)
                            
                            iPURange = 0
                            for PURangeName, PURange in nPVRanges.items():
                                # h3D histogram:: X-axis: iEta, Y-axis: nVts, Z-axis: resolution
                                kBin_PURangeMin = h3D.GetYaxis().FindBin(PURange[0])
                                kBin_PURangeMax = h3D.GetYaxis().FindBin(PURange[1])
                                kBin_iEtaMin    = 1
                                kBin_iEtaMax    = h3D.GetNbinsX()
                                print "H3D {}, PURangeName: {}, PURange {}".format(h3D.GetName(), PURangeName, PURange)
                                #h2D = h3D.ProjectionZ("%s_%s" % (h3D.GetName(), PURangeName), kBin_iEtaMin,kBin_iEtaMax, kBin_PURangeMin,kBin_PURangeMax)
                                h3D.GetYaxis().SetRange(kBin_PURangeMin, kBin_PURangeMax)
                                h2D = h3D.Project3D("zx") # 1st lable gets plot along y-axis
                                sNameTmp = "%s_%s" % (sNameTmp0, PURangeName)
                                h2D.SetNameTitle(sNameTmp, sNameTmp)
                                
                                jetShape_PUS_PUrange = "%s_%s_%s" % (jetShape,JetPUSName,PURangeName)
                                
                                if not useAbsEtaBins:
                                    h2D.GetXaxis().SetTitle("iEta")
                                else:
                                    h2D.GetXaxis().SetTitle("|i#eta|")
                                print "res_histoName new: {}".format(h2D.GetName()); sys.stdout.flush()

                                yaxisName = ""
                                if   plot1DName == 'res':
                                    yaxisName_mean  = "#mu(L1T/PF - 1)"
                                    yaxisName_sigma = "#sigma(L1T/PF - 1)"
                                    yaxisName_sigmaByMean = "#frac{#sigma(L1T/PF - 1)}{#mu(L1T/PF - 1) + 1}"
                                    fGaus.SetRange(-1.4, 2.5)
                                    axisRange_mean  = [-0.8, 1]
                                    axisRange_sigma = [0.1, 0.9]
                                elif plot1DName == 'PU':
                                    yaxisName_mean  = "#mu(L1 jet PU E_{T}) [GeV]"
                                    yaxisName_sigma = "#sigma(L1 jet PU E_{T}) [GeV]"
                                    fGaus.SetRange(0, 60)                        
                                    axisRange_mean  = [0, 60]
                                    axisRange_sigma = [0, 20]
                                elif plot1DName == 'PUByRawPt':
                                    yaxisName_mean  = "#mu(#frac{L1 jet PU E_{T}}{L1 jet raw E_{T}}) [GeV]"
                                    yaxisName_sigma = "#sigma(#frac{L1 jet PU E_{T}}{L1 jet raw E_{T}}) [GeV]"
                                    fGaus.SetRange(0.0, 1.3)
                                    axisRange_mean  = [0, 0.9]
                                    axisRange_sigma = [0., 0.4]
                                    #h2D.RebinY(2)

                                c2_1.cd()
                                h2D.FitSlicesY(fGaus, 1,h2D.GetNbinsX(), 1)


                                c2_1.cd()
                                #hres2D_1 = R.gDirectory.Get("%s_1" % (res_histoName))
                                h2D_1 = (R.gDirectory.Get("%s_1" % (sNameTmp)))
                                #h2D_1.SetNameTitle("%s_1" % (sNameTmp), "%s_1" % (sNameTmp))
                                h2D_1.SetNameTitle("%s_Mean" % (sNameTmp), "%s_Mean" % (sNameTmp))
                                h2D_1.SetLineColor(colors[iResType])
                                h2D_1.SetMarkerColor(colors[iResType])
                                #h2D_1.SetMarkerStyle(markers[iInFile])
                                h2D_1.SetMarkerStyle(markers[iPURange])
                                h2D_1.SetMarkerSize(0.7)
                                h2D_1.GetYaxis().SetTitle(yaxisName_mean)

                                c3_1.cd()
                                #hres2D_2 = R.gDirectory.Get("%s_2" % (res_histoName))
                                h2D_2 = ( R.gDirectory.Get("%s_2" % (sNameTmp)) )
                                #h2D_2.SetNameTitle("%s_2" % (sNameTmp), "%s_2" % (sNameTmp))
                                h2D_2.SetNameTitle("%s_Sigma" % (sNameTmp), "%s_Sigma" % (sNameTmp))
                                h2D_2.SetLineColor(colors[iResType])
                                h2D_2.SetMarkerColor(colors[iResType])
                                #h2D_2.SetMarkerStyle(markers[iInFile])
                                h2D_2.SetMarkerStyle(markers[iPURange])
                                h2D_2.SetMarkerSize(0.7)
                                h2D_2.GetYaxis().SetTitle(yaxisName_sigma)

                                c4_1.cd()
                                h2D_MeanPlus1    = h2D_1.Clone("%s_MeanPlus1" % (sNameTmp))
                                h2D_SF_PFJByL1TJ = h2D_1.Clone("%s_SF_PFJByL1TJ" % (sNameTmp))
                                h2D_SF_PFJByL1TJ.GetYaxis().SetTitle("#frac{1}{#mu(L1T/PF - 1) + 1}")
                                for iBin in range(1, h2D_MeanPlus1.GetNbinsX()+1):
                                    binContent    = h2D_MeanPlus1.GetBinContent(iBin)
                                    binError      = h2D_MeanPlus1.GetBinError(iBin)
                                    MeanPlus1     = binContent + 1
                                    binContent_SF = 1 / (MeanPlus1)
                                    binError_SF   = binError / (MeanPlus1 * MeanPlus1)
                                    h2D_MeanPlus1.SetBinContent(iBin,    MeanPlus1)
                                    h2D_SF_PFJByL1TJ.SetBinContent(iBin, binContent_SF)
                                    h2D_SF_PFJByL1TJ.SetBinError(iBin,   binError_SF)
                                print "h2D_SF_PFJByL1TJ: {}".format(h2D_SF_PFJByL1TJ.GetName())    
                                hEffectiveResolution = h2D_2.Clone("%s_EffectiveResolution" % (sNameTmp))
                                hEffectiveResolution.Divide(h2D_2, h2D_MeanPlus1)
                                hEffectiveResolution.GetYaxis().SetTitle(yaxisName_sigmaByMean)

                                if   plot1DName == 'res':
                                    fOpL1JetEnergySF.cd();
                                    h2D_SF_PFJByL1TJ.Write()

                                hres2D_dict[sInFileTag][jetShape_PUS_PUrange] = h2D
                                hresMean_dict[sInFileTag][jetShape_PUS_PUrange] = h2D_1
                                hresSigma_dict[sInFileTag][jetShape_PUS_PUrange] = h2D_2
                                hEffectiveResolution_dict[sInFileTag][jetShape_PUS_PUrange] = hEffectiveResolution

                                # make ratio plot of "highPU / lowPU"
                                if PURangeName not in JetShapesAndPUSs_Denom_forRatioPlot: # JetShapesAndPUSs_Denom_forRatioPlot contains label for denomintor histogram
                                    jetShape_PUS_PUrange_forRatioPlot = "%s_%s_%s" % (jetShape,JetPUSName,JetShapesAndPUSs_Denom_forRatioPlot[0])
                                    sYName  = "#frac{h(%s)}{h(%s)}" % (PURangeName, JetShapesAndPUSs_Denom_forRatioPlot[0])
                                    sYName1 = "#frac{h(%s) + 1}{h(%s) + 1}" % (PURangeName, JetShapesAndPUSs_Denom_forRatioPlot[0])
                                    
                                    hRatio2D_1       = cloneHistogram(hresMean_dict[sInFileTag][jetShape_PUS_PUrange],              "ratio")
                                    hRatio2D_1_Denom = cloneHistogram(hresMean_dict[sInFileTag][jetShape_PUS_PUrange_forRatioPlot], "ratio")
                                    # make "mean + 1" so as to get L1TJetEt scale factor
                                    for iBin in range(1, hRatio2D_1.GetNbinsX()):
                                        hRatio2D_1      .SetBinContent(iBin, (hRatio2D_1      .GetBinContent(iBin) + 1) )
                                        hRatio2D_1_Denom.SetBinContent(iBin, (hRatio2D_1_Denom.GetBinContent(iBin) + 1) )
                                    hRatio2D_1.Divide(hRatio2D_1_Denom)
                                    
                                    hRatio2D_2       = cloneHistogram(hresSigma_dict[sInFileTag][jetShape_PUS_PUrange],              "ratio")
                                    hRatio2D_2_Denom = cloneHistogram(hresSigma_dict[sInFileTag][jetShape_PUS_PUrange_forRatioPlot], "ratio")
                                    hRatio2D_2.Divide(hRatio2D_2_Denom)

                                    hRatio2D_3       = cloneHistogram(hEffectiveResolution_dict[sInFileTag][jetShape_PUS_PUrange],              "ratio")
                                    hRatio2D_3_Denom = cloneHistogram(hEffectiveResolution_dict[sInFileTag][jetShape_PUS_PUrange_forRatioPlot], "ratio")
                                    hRatio2D_3.Divide(hRatio2D_3_Denom)

                                    hRatioresMean_dict[sInFileTag][jetShape_PUS_PUrange]             = hRatio2D_1
                                    hRatioresSigma_dict[sInFileTag][jetShape_PUS_PUrange]            = hRatio2D_2
                                    hRatioEffectiveResolution_dict[sInFileTag][jetShape_PUS_PUrange] = hRatio2D_3
                                    
                                    hRatioresMean_dict[sInFileTag][jetShape_PUS_PUrange].GetYaxis().SetTitle(sYName1)
                                    hRatioresSigma_dict[sInFileTag][jetShape_PUS_PUrange].GetYaxis().SetTitle(sYName)
                                    hRatioEffectiveResolution_dict[sInFileTag][jetShape_PUS_PUrange].GetYaxis().SetTitle(sYName)

                                if iPURange == 0:
                                    jetShape_PUS_PUrange_nice = "%s %s" % (jetShape,JetPUSName)
                                    jetShape_PUS_PUrange_nice = jetShape_PUS_PUrange_nice.replace('_plus_', ' + ')
                                    jetShape_PUS_PUrange_nice = jetShape_PUS_PUrange_nice.replace('_times_', ' * ')
                                    legend.AddEntry(hresSigma_dict[sInFileTag][jetShape_PUS_PUrange], jetShape_PUS_PUrange_nice, "lep")
                                iPURange += 1
                            
                            # Different colors for different JetSpaeAndPUS, but same color for PUrange
                            iResType += 1
                            
                    '''            
                    shRatio_Denom_JetShape = "%s_%s_%s" % ( 
                        JetShapesAndPUSs_Denom_forRatioPlot[0],
                        JetShapesAndPUSs_Denom_forRatioPlot[1],
                        JetShapesAndPUSs_Denom_forRatioPlot[2],
                    )
                    '''




                    hresMean_list = []
                    hresSigma_list = []
                    hEffectiveResolution_list = []
                    hRatioresMean_list = []
                    hRatioresSigma_list = []
                    hRatioEffectiveResolution_list = []
                    #for jetShape in JetShapes:
                    for jetShape in hEffectiveResolution_dict[sInFileTag].keys():
                        hresMean_list.append(hresMean_dict[sInFileTag][jetShape])
                        hresSigma_list.append(hresSigma_dict[sInFileTag][jetShape])
                        hEffectiveResolution_list.append(hEffectiveResolution_dict[sInFileTag][jetShape])
                        
                    for jetShape in hRatioEffectiveResolution_dict[sInFileTag].keys():
                        hRatioresMean_list.append(hRatioresMean_dict[sInFileTag][jetShape])
                        hRatioresSigma_list.append(hRatioresSigma_dict[sInFileTag][jetShape])
                        hRatioEffectiveResolution_list.append(hRatioEffectiveResolution_dict[sInFileTag][jetShape])
                        
                        '''
                        if jetShape == shRatio_Denom_JetShape: continue
                        
                        hRatioresMean_dict[sInFileTag][jetShape] = hresMean_dict[sInFileTag][jetShape].Clone("%s_Ratio" % hresMean_dict[sInFileTag][jetShape].GetName())
                        hRatioresSigma_dict[sInFileTag][jetShape] = hresSigma_dict[sInFileTag][jetShape].Clone("%s_Ratio" % hresSigma_dict[sInFileTag][jetShape].GetName())
                        hRatioEffectiveResolution_dict[sInFileTag][jetShape] = hEffectiveResolution_dict[sInFileTag][jetShape].Clone("%s_Ratio" % hEffectiveResolution_dict[sInFileTag][jetShape].GetName())
                        
                        hRatioresMean_dict[sInFileTag][jetShape].Divide(hresMean_dict[sInFileTag][jetShape], hresMean_dict[sInFileTag][shRatio_Denom_JetShape])
                        hRatioresSigma_dict[sInFileTag][jetShape].Divide(hresSigma_dict[sInFileTag][jetShape], hresSigma_dict[sInFileTag][shRatio_Denom_JetShape])
                        hRatioEffectiveResolution_dict[sInFileTag][jetShape].Divide(hEffectiveResolution_dict[sInFileTag][jetShape], hEffectiveResolution_dict[sInFileTag][shRatio_Denom_JetShape])
                        
                        sYName = "#frac{h(i)}{h(%s)}" % (shRatio_Denom_JetShape)
                        hRatioresMean_dict[sInFileTag][jetShape].GetYaxis().SetTitle(sYName)
                        hRatioresSigma_dict[sInFileTag][jetShape].GetYaxis().SetTitle(sYName)
                        hRatioEffectiveResolution_dict[sInFileTag][jetShape].GetYaxis().SetTitle(sYName)
                        
                        hRatioresMean_list.append(hRatioresMean_dict[sInFileTag][jetShape])
                        hRatioresSigma_list.append(hRatioresSigma_dict[sInFileTag][jetShape])
                        hRatioEffectiveResolution_list.append(hRatioEffectiveResolution_dict[sInFileTag][jetShape])
                        '''
                    
                    
                    print "here1"
                    print "iInFile: {}".format(iInFile)
                    print "sInFileTags[iInFile]: {}".format(sInFileTags[iInFile])
                    print "sResLegend_selected: {}".format(sResLegend_selected)
                    sDir1 = "plots%s/CompareDiffJetShapesAndPUS/%s/%s/%s" % (sInFileVersion, plot1DName, resPtBin, sInFileTags[iInFile])
                    if not os.path.exists(sDir1):
                        os.makedirs(sDir1)                                
                    c2_1 = plotHistos1DAndRatioPlot(hresMean_list, hRatioresMean_list, legend, c2_1, sSaveAs="%s/DiffJetShapesAndPUS_%s_%s_mean.png" % (sDir1, l1Mode, plotName))
                    c3_1 = plotHistos1DAndRatioPlot(hresSigma_list, hRatioresSigma_list, legend, c3_1, sSaveAs="%s/DiffJetShapesAndPUS_%s_%s_sigma.png" % (sDir1, l1Mode, plotName))
                    c4_1 = plotHistos1DAndRatioPlot(hEffectiveResolution_list, hRatioEffectiveResolution_list, legend, c4_1, sSaveAs="%s/DiffJetShapesAndPUS_%s_%s_effectiveResolution.png" % (sDir1, l1Mode, plotName))
                    '''
                    c1.Update()
                    c2_1.Update()
                    c3_1.Update()
                    '''
                    #c2_1.SaveAs("%s/DiffJetShapes_%s_%s_mean.png" % (sDir1, l1Mode, plotName))
                    #c3_1.SaveAs("%s/DiffJetShapes_%s_%s_sigma.png" % (sDir1, l1Mode, plotName))
        
        
        fOpL1JetEnergySF.Close()            



            
    if 'run_TrgRates' in run_list:
        c1 = R.TCanvas("c1","c1",500,400)
        c2 = R.TCanvas("c2","c2",600,800)
        #c3 = R.TCanvas("c3","c3",1200,600)
        #c2_1 = R.TCanvas("c2_1","c2_1",canvasDims_iEta[0],canvasDims_iEta[1])
        #c3_1 = R.TCanvas("c3_1","c3_1",canvasDims_iEta[0],canvasDims_iEta[1])
        #c4_1 = R.TCanvas("c4_1","c4_1",canvasDims_iEta[0],canvasDims_iEta[1])
        
        
        
        
        ### Trigger rates
        for iInFile in range(len(in_fileNames)):
            fIn         = inFiles[iInFile]
            sInFileTag  = sInFileTags[iInFile]
            
            if 'ZeroBias' not in in_fileNames[iInFile]:
                print "i/p file {} for trigger rate plots is not ZeroBias dataset ????? *** ERROR ??? ***".format(in_fileNames[iInFile])
                return
            
            # read histograms for nTotalEvents
            hnTotalEvents = fIn.Get(sHistoName_nTotalEvents_forTrgRates)
            if not hnTotalEvents:
                print "Histogram {} doesn't exist in {} \t\t *** ERROR ***".format(sHistoName_nTotalEvents_forTrgRates, in_fileNames[iInFile])
                continue
            
            nTotalEvents = hnTotalEvents.GetBinContent(1)            
            print "{}: nTotalEvents: {}".format(in_fileNames[iInFile], nTotalEvents)
            norm_TrgRates = instLumi * mbXSec / nTotalEvents
            print "instLumi {}, mbXSec {}, nTotalEvents {}, norm_TrgRates {}".format(instLumi, mbXSec, nTotalEvents, norm_TrgRates)
            
            
            
            for l1Mode in l1Modes:
                
                for eta_cat in IETA_CAT.keys():
                    iPURange = 0
                    for PURangeName, PURange in nPVRanges.items():
                        
                        for TrgRate_Type in TrgRates_Types:
                        
                            #legend = R.TLegend(0.2,0.8,1,1)
                            legend = R.TLegend(0.5,0.85,1,1)
                            legend.SetNColumns(2)

                            iResType = 0
                            h_OD = OrderedDict()
                            hRatio_OD = OrderedDict() 
                            for jetShape, JetPUSs in JetShapesAndPUSs_histoName_TrgRates.items():
                                # JetShape = "" plots are with the first version of code for 9x9 jets
                                jetShape1 = jetShape
                                if jetShape == 'Default':  jetShape1 = ""
                                else:                      jetShape1 = "_%s" % (jetShape)
                                #print "  jetShape1: {}".format(jetShape1)


                                for JetPUSName, JetPUSHistName in JetPUSs.items():
                                    resType = JetPUSHistName
                                    resType = resType.replace('$JETSHAPE', jetShape1)
                                    resType = resType.replace('$ETACAT', eta_cat)
                                    resType = resType.replace(TrgRates_Types[0], TrgRate_Type)
                                    res_histoName = "%s_%s" % (resType, l1Mode)
                                    print "res histo: %s" % (res_histoName); sys.stdout.flush()
                                    h2D = (fIn.Get(res_histoName))
                                    #sNameTmp0 = "%s_InFile%d" % (res_histoName, iInFile)
                                    sNameTmp0 = "%s_%s" % (res_histoName, sInFileTags[iInFile])
                                    h2D.SetNameTitle(sNameTmp0, sNameTmp0)

                                    h2D.Scale( norm_TrgRates )

                                    c1.cd()
                                    h2D.Draw("colz")

                                    # take projection along the selected PU range
                                    kBin_PURangeMin = h2D.GetYaxis().FindBin(PURange[0])
                                    kBin_PURangeMax = h2D.GetYaxis().FindBin(PURange[1])

                                    h1D = h2D.ProjectionX("%s_ProjX%s"%(h2D.GetName(), PURangeName),
                                                          kBin_PURangeMin, kBin_PURangeMax )

                                    h1D.GetXaxis().SetTitle("Threshold E_{T} (GeV)")
                                    h1D.GetYaxis().SetTitle("Rates (Hz)")
                                    h1D.GetXaxis().SetRangeUser(axisRange_EtTrsh_forTrgRates[0], axisRange_EtTrsh_forTrgRates[1])
                                    
                                    h1D.SetLineWidth(2)
                                    h1D.SetLineColor(colors[iResType])
                                    h1D.SetMarkerColor(colors[iResType])
                                    #h1D.SetMarkerStyle(markers[iInFile])
                                    #h1D.SetMarkerStyle(markers[iPURange])
                                    #h1D.SetMarkerSize(0.7)
                                    #h1D.GetYaxis().SetRangeUser(1e3, 1e10)

                                    sJetShapeAndPUS_current = "%s_%s" % (jetShape, JetPUSName)
                                    sJetShapeAndPUS_Denom = "%s_%s" % (JetShapesAndPUSs_Denom_forRatioPlot_forTrgRates[0], JetShapesAndPUSs_Denom_forRatioPlot_forTrgRates[1])
                                    h_OD[sJetShapeAndPUS_current] = h1D

                                    # legend
                                    jetShape_PUS_PUrange_nice = "%s %s" % (jetShape,JetPUSName)
                                    jetShape_PUS_PUrange_nice = jetShape_PUS_PUrange_nice.replace('_plus_', ' + ')
                                    jetShape_PUS_PUrange_nice = jetShape_PUS_PUrange_nice.replace('_times_', ' * ')
                                    legend.AddEntry(h_OD[sJetShapeAndPUS_current], jetShape_PUS_PUrange_nice, 'l')

                                    iResType += 1

                                    if sJetShapeAndPUS_current == sJetShapeAndPUS_Denom: continue                                

                                    # ratio plot
                                    hRatio       = cloneHistogram(h_OD[sJetShapeAndPUS_current],        "ratio")
                                    hRatio_Denom = cloneHistogram(h_OD[sJetShapeAndPUS_Denom],          "ratio")
                                    hRatio.Divide(hRatio_Denom)
                                    #hRatio.GetYaxis().SetTitle("#frac{New}{Default RAWPUS}")
                                    hRatio.GetYaxis().SetTitle("#frac{New}{%s %s}" % (JetShapesAndPUSs_Denom_forRatioPlot_forTrgRates[0], JetShapesAndPUSs_Denom_forRatioPlot_forTrgRates[1]))

                                    hRatio_OD[sJetShapeAndPUS_current] = hRatio



                            h_list = []
                            hRatio_list = []
                            for jetShape in h_OD.keys():
                                h_list.append(h_OD[jetShape])

                            for jetShape in hRatio_OD.keys():
                                hRatio_list.append(hRatio_OD[jetShape])

                            sDir1 = "plots%s/CompareDiffJetShapesAndPUS/rates/%s" % (sInFileVersion, sInFileTags[iInFile])
                            if not os.path.exists(sDir1):
                                os.makedirs(sDir1)
                            yRange1=[1e1,1e8, 0,2]
                            if TrgRate_Type == 'singleJet':
                                yRange1=[1e3,1e8, 0,2]
                            if TrgRate_Type == 'doubleJet':
                                yRange1=[1e3,1e7, 0,2]
                            if TrgRate_Type == 'trippleJet':
                                yRange1=[1e2,1e6, 0,2]
                            if TrgRate_Type == 'quadJet':
                                yRange1=[1e1,1e5, 0,2]
                            #yRange1 = [yRange1[0], yRange1[1]]
                            yRange1=[]
                            c2 = plotHistos1DAndRatioPlot(h_list, hRatio_list, legend, c2, sSaveAs="%s/rates_%s_%s_%s_%s.png" % (sDir1, TrgRate_Type, eta_cat, l1Mode, PURangeName), drawOption="HIST", pad1SetLogY=True, yRange=yRange1)

                    
                    iPURange += 1







    if 'run_TrgEffi' in run_list:
        #c1 = R.TCanvas("c1","c1",500,400)
        c1 = R.TCanvas("c1","c1",1200,1400); c1.Divide(2,3);
        c2 = R.TCanvas("c2","c2",600,800)
        c3 = R.TCanvas("c3","c3",1200,600); c3.Divide(2,1)
        #c2_1 = R.TCanvas("c2_1","c2_1",canvasDims_iEta[0],canvasDims_iEta[1])
        #c3_1 = R.TCanvas("c3_1","c3_1",canvasDims_iEta[0],canvasDims_iEta[1])
        #c4_1 = R.TCanvas("c4_1","c4_1",canvasDims_iEta[0],canvasDims_iEta[1])
        
        


        ### Trigger efficiency
        for iInFile in range(len(in_fileNames)):
            fIn         = inFiles[iInFile]
            sInFileTag  = sInFileTags[iInFile]
            
            if 'SingleMu' not in in_fileNames[iInFile]:
                print "i/p file {} for trigger efficiency plots is not SingleMu dataset ????? *** ERROR ??? ***".format(in_fileNames[iInFile])
                return
                       
            
            
            for l1Mode in l1Modes:
                
                for eta_cat in IETA_CAT.keys():
                    iPURange = 0
                    for PURangeName, PURange in nPVRanges.items():

                        for TrgTrsh in TrgEffi_TrgTrshs:

                            legend = R.TLegend(0.2,0.85,1,1)
                            legend.SetNColumns(2)

                            iResType = 0
                            h_OD = OrderedDict()
                            hRatio_OD = OrderedDict()
                            c1.Clear(); c1.Divide(2,3);
                            for jetShape, JetPUSs in JetShapesAndPUSs_histoName_TrgEffi.items():
                                # JetShape = "" plots are with the first version of code for 9x9 jets
                                jetShape1 = jetShape
                                if jetShape == 'Default':  jetShape1 = ""
                                else:                      jetShape1 = "_%s" % (jetShape)
                                #print "  jetShape1: {}".format(jetShape1)


                                for JetPUSName, JetPUSHistName in JetPUSs.items():
                                    resType = JetPUSHistName
                                    resType = resType.replace('$JETSHAPE', jetShape1)
                                    resType = resType.replace('$ETACAT', eta_cat)
                                    resType = resType.replace('$TRGTRSH', str(TrgTrsh))
                                    res_histoName = "%s_%s" % (resType, l1Mode)
                                    res_histoName_den = res_histoName.replace('num', 'den')

                                    print "res histo_num: %s \t histo_num: %s" % (res_histoName, res_histoName_den); sys.stdout.flush()
                                    h2D_num = (fIn.Get(res_histoName))
                                    sNameTmp0 = "%s_%s" % (res_histoName, sInFileTags[iInFile])
                                    h2D_num.SetNameTitle(sNameTmp0, sNameTmp0)
                                    h2D_den = (fIn.Get(res_histoName_den))
                                    sNameTmp0 = "%s_%s" % (res_histoName_den, sInFileTags[iInFile])
                                    h2D_den.SetNameTitle(sNameTmp0, sNameTmp0)


                                    # take projection along the selected PU range
                                    kBin_PURangeMin = h2D_num.GetYaxis().FindBin(PURange[0])
                                    kBin_PURangeMax = h2D_num.GetYaxis().FindBin(PURange[1])

                                    h1D_num = h2D_num.ProjectionX("%s_ProjX%s"%(h2D_num.GetName(), PURangeName),
                                                                  kBin_PURangeMin, kBin_PURangeMax )
                                    h1D_den = h2D_den.ProjectionX("%s_ProjX%s"%(h2D_den.GetName(), PURangeName),
                                                                  kBin_PURangeMin, kBin_PURangeMax )

                                    rebin = axisRebin_PFJetPt_forTrgEffi # 4
                                    h1D_num.Rebin(rebin)
                                    h1D_den.Rebin(rebin)
                                    h1D_num.GetXaxis().SetRangeUser(axisRange_PFJetPt_forTrgEffi[0], axisRange_PFJetPt_forTrgEffi[1])
                                    h1D_den.GetXaxis().SetRangeUser(axisRange_PFJetPt_forTrgEffi[0], axisRange_PFJetPt_forTrgEffi[1])
                                    
                                    '''
                                    h1D = R.TEfficiency(h1D_num, h1D_den)
                                    #h1D.GetXaxis().SetTitle("Offline Jet E_{T} (GeV)")
                                    #h1D.GetYaxis().SetTitle("Efficieny")
                                    h1D.SetName("%s_Efficiency_%s"%(h2D_num.GetName(), PURangeName))
                                    '''
                                    
                                    h1D = R.TGraphAsymmErrors()
                                    h1D.BayesDivide(h1D_num, h1D_den)

                                    h1D.GetXaxis().SetTitle("Offline jet E_{T} (GeV)")
                                    h1D.GetYaxis().SetTitle("Efficiency")
                                    #h1D.GetXaxis().SetRangeUser(20,500)
                                    h1D.GetXaxis().SetRangeUser(axisRange_PFJetPt_forTrgEffi[0], axisRange_PFJetPt_forTrgEffi[1])
                                    
                                    #h1D.SetLineWidth(2)
                                    h1D.SetLineColor(colors[iResType])
                                    h1D.SetMarkerColor(colors[iResType])
                                    #h1D.SetMarkerStyle(markers[iInFile])
                                    h1D.SetMarkerStyle(markers[iPURange])
                                    h1D.SetMarkerSize(0.7)
                                    #h1D.GetYaxis().SetRangeUser(1e3, 1e10)

                                    
                                    
                                    
                                    c1.cd(1)
                                    h2D_num.Draw("colz")
                                    c1.cd(2)
                                    h2D_den.Draw("colz")
                                    
                                    c1.cd(3)
                                    h1D_num.Draw()
                                    c1.cd(4)
                                    h1D_den.Draw()
                                    c1.cd(5)
                                    h1D.Draw()
                                    
                                    c1.Update()
                                    

                                    sJetShapeAndPUS_current = "%s_%s" % (jetShape, JetPUSName)
                                    sJetShapeAndPUS_Denom = "%s_%s" % (JetShapesAndPUSs_Denom_forRatioPlot_forTrgRates[0], JetShapesAndPUSs_Denom_forRatioPlot_forTrgRates[1])
                                    print "sJetShapeAndPUS_current {}, sJetShapeAndPUS_Denom {}".format(sJetShapeAndPUS_current, sJetShapeAndPUS_Denom)
                                    h_OD[sJetShapeAndPUS_current] = h1D

                                    # legend
                                    jetShape_PUS_PUrange_nice = "%s %s" % (jetShape,JetPUSName)
                                    jetShape_PUS_PUrange_nice = jetShape_PUS_PUrange_nice.replace('_plus_', ' + ')
                                    jetShape_PUS_PUrange_nice = jetShape_PUS_PUrange_nice.replace('_times_', ' * ')
                                    legend.AddEntry(h_OD[sJetShapeAndPUS_current], jetShape_PUS_PUrange_nice, 'lep')

                                    iResType += 1

                                    if sJetShapeAndPUS_current == sJetShapeAndPUS_Denom: continue                                

                                    # ratio plot
                                    hRatio       = cloneHistogram(h_OD[sJetShapeAndPUS_current],        "ratio")
                                    hRatio_Denom = cloneHistogram(h_OD[sJetShapeAndPUS_Denom],          "ratio")
                                    c3.Clear(); c3.Divide(2,1)
                                    c3.cd(1)
                                    hRatio.Draw()
                                    c3.cd(2)
                                    hRatio_Denom.Draw()
                                    c3.Update()
                                    #hRatio.Divide(hRatio_Denom)
                                    hRatio = DivideGraphAsymmErr(hRatio, hRatio_Denom)
                                    #hRatio.GetYaxis().SetTitle("#frac{New}{Default RAWPUS}")

                                    if hRatio == None: continue
                                    
                                    hRatio.GetXaxis().SetTitle("Offline jet E_{T} (GeV)")
                                    #hRatio.GetYaxis().SetTitle("#frac{New}{Default RAWPUS}")
                                    hRatio.GetYaxis().SetTitle("#frac{New}{%s %s}" % (JetShapesAndPUSs_Denom_forRatioPlot_forTrgRates[0], JetShapesAndPUSs_Denom_forRatioPlot_forTrgRates[1]))
                                    #hRatio.GetXaxis().SetRangeUser(20,500)
                                    hRatio.GetXaxis().SetRangeUser(axisRange_PFJetPt_forTrgEffi[0], axisRange_PFJetPt_forTrgEffi[1])
                                    
                                    #hRatio.SetLineWidth(2)
                                    hRatio.SetLineColor(colors[iResType - 1])
                                    hRatio.SetMarkerColor(colors[iResType - 1])
                                    #hRatio.SetMarkerStyle(markers[iInFile])
                                    hRatio.SetMarkerStyle(markers[iPURange])
                                    hRatio.SetMarkerSize(0.7)
                                    #hRatio.GetYaxis().SetRangeUser(1e3, 1e10)

                                    hRatio_OD[sJetShapeAndPUS_current] = hRatio



                            h_list = []
                            hRatio_list = []
                            for jetShape in h_OD.keys():
                                h_list.append(h_OD[jetShape])

                            for jetShape in hRatio_OD.keys():
                                hRatio_list.append(hRatio_OD[jetShape])

                            sDir1 = "plots%s/CompareDiffJetShapesAndPUS/efficiency/%s" % (sInFileVersion, sInFileTags[iInFile])
                            if not os.path.exists(sDir1):
                                os.makedirs(sDir1)                                
                            #c2 = plotHistos1DAndRatioPlot(h_list, hRatio_list, legend, c2, sSaveAs="%s/effi_%s_%s_%s.png" % (sDir1, eta_cat, l1Mode, PURangeName), drawOption="HIST", pad1SetLogY=True, yRange=[0,1.1, 0,2], setLogX=True)
                            c2 = plotHistos1DAndRatioPlot(h_list, hRatio_list, legend, c2, sSaveAs="%s/effi_%s_%s_TrgTrsh%s_%s.png" % (sDir1, eta_cat, l1Mode, str(TrgTrsh), PURangeName), drawOption="PA", pad1SetLogY=False, yRange=[0,1.2, 0,2], setLogX=False)

                    
                    iPURange += 1









    if 'run_res_vs_tauJets' in run_list:
        c1 = R.TCanvas("c1","c1",500,400)
        #c2 = R.TCanvas("c2","c2",1200,600)
        #c3 = R.TCanvas("c3","c3",1200,600)
        c2_1 = R.TCanvas("c2_1","c2_1",canvasDims_iEta[0],canvasDims_iEta[1])
        c3_1 = R.TCanvas("c3_1","c3_1",canvasDims_iEta[0],canvasDims_iEta[1])
        c4_1 = R.TCanvas("c4_1","c4_1",canvasDims_iEta[0],canvasDims_iEta[1])
        
        
        fOpL1JetEnergySF = R.TFile(sOpL1JetEnergySF, "RECREATE")
        
        
        
        for plot1DName in ['res']: #['res', 'PU', 'PUByRawPt']: # ['jet_byHand_res_vs_iEta', 'jet_byHand_PU_vs_iEta', 'jet_byHand_PUByRawPt_vs_iEta']
            plotName = ""
            if   plot1DName == 'res':       plotName = "Res"
            elif plot1DName == 'PU':        plotName = "PU"
            elif plot1DName == 'PUByRawPt': plotName = "PUByRawPt"

            for l1Mode in l1Modes:
                hDummy = R.TH1D("hDummy", "", 1,0,1)
        
                hres2D_dict    = OrderedDict()
                hresMean_dict  = OrderedDict()
                hresSigma_dict = OrderedDict()
                hEffectiveResolution_dict = OrderedDict()
                hRatiores2D_dict    = OrderedDict()
                hRatioresMean_dict  = OrderedDict()
                hRatioresSigma_dict = OrderedDict()
                hRatioEffectiveResolution_dict = OrderedDict()
                for iInFile in range(len(in_fileNames)):
                    fIn         = inFiles[iInFile]
                    sInFileTag  = sInFileTags[iInFile]                 

                    legend = R.TLegend(0.35,0.8,1,1)
                    legend.SetNColumns(2)
                    
                    hres2D_dict[sInFileTag] = OrderedDict()
                    hresMean_dict[sInFileTag] = OrderedDict()
                    hresSigma_dict[sInFileTag] = OrderedDict()
                    hEffectiveResolution_dict[sInFileTag] = OrderedDict()
                    hRatiores2D_dict[sInFileTag] = OrderedDict()
                    hRatioresMean_dict[sInFileTag] = OrderedDict()
                    hRatioresSigma_dict[sInFileTag] = OrderedDict()
                    hRatioEffectiveResolution_dict[sInFileTag] = OrderedDict()
                    shRatio_Denom_JetShape = None
                    iPURange = 0
                    hTmps = []
                    for PURangeName, PURange in nPVRanges.items():
                        hTmp = R.TH1D("hTmp%d"%(iPURange), "", 1,0,1)
                        hTmp.SetMarkerStyle(markers[iPURange])
                        hTmp.SetLineColor(R.kBlack)
                        hTmp.SetMarkerSize(0.7)
                        hTmps.append(hTmp)
                        legend.AddEntry(hTmps[iPURange], PURangeName, "lep")
                        iPURange += 1
                    
                    iResType = 0
                    for jetShape, JetPUSHistName in L1TauMatchedToPFJetVersions.items():
                        resType = JetPUSHistName
                        resType = resType.replace('res', plot1DName)
                        res_histoName = "%s_%s" % (resType, l1Mode)
                        print "res histo: %s" % (res_histoName); sys.stdout.flush()
                        h3D = (fIn.Get(res_histoName))
                        #sNameTmp0 = "%s_InFile%d" % (res_histoName, iInFile)
                        sNameTmp0 = "%s_%s" % (res_histoName, sInFileTags[iInFile])
                        print "h3D {}, sNameTmp0 {}".format(h3D, sNameTmp0)
                        h3D.SetNameTitle(sNameTmp0, sNameTmp0)

                        iPURange = 0
                        for PURangeName, PURange in nPVRanges.items():
                            # h3D histogram:: X-axis: iEta, Y-axis: nVts, Z-axis: resolution
                            kBin_PURangeMin = h3D.GetYaxis().FindBin(PURange[0])
                            kBin_PURangeMax = h3D.GetYaxis().FindBin(PURange[1])
                            kBin_iEtaMin    = 1
                            kBin_iEtaMax    = h3D.GetNbinsX()
                            print "H3D {}, PURangeName: {}, PURange {}".format(h3D.GetName(), PURangeName, PURange)
                            #h2D = h3D.ProjectionZ("%s_%s" % (h3D.GetName(), PURangeName), kBin_iEtaMin,kBin_iEtaMax, kBin_PURangeMin,kBin_PURangeMax)
                            print "H3D: GetEntries {}, Xaxis: {}, {}, {}, Yaxis: {}, {}, {}, Zaxis: {}, {}, {},".format(
                                h3D.GetEntries(),
                                h3D.GetXaxis().GetNbins(),h3D.GetXaxis().GetXmin(),h3D.GetXaxis().GetXmin(),
                                h3D.GetYaxis().GetNbins(),h3D.GetYaxis().GetXmin(),h3D.GetYaxis().GetXmin(),
                                h3D.GetZaxis().GetNbins(),h3D.GetZaxis().GetXmin(),h3D.GetZaxis().GetXmin(),
                                
                            )
                            
                            h3D.GetYaxis().SetRange(kBin_PURangeMin, kBin_PURangeMax)
                            h2D = h3D.Project3D("zx") # 1st lable gets plot along y-axis
                            sNameTmp = "%s_%s" % (sNameTmp0, PURangeName)
                            h2D.SetNameTitle(sNameTmp, sNameTmp)

                            c1.cd()
                            h2D.Draw('colz')

                            jetShape_PUS_PUrange = "%s_%s" % (jetShape,PURangeName)

                            if not useAbsEtaBins:
                                h2D.GetXaxis().SetTitle("iEta")
                            else:
                                h2D.GetXaxis().SetTitle("|i#eta|")
                            print "res_histoName new: {}".format(h2D.GetName()); sys.stdout.flush()

                            yaxisName = ""
                            if   plot1DName == 'res':
                                yaxisName_mean  = "#mu(L1T/PF - 1)"
                                yaxisName_sigma = "#sigma(L1T/PF - 1)"
                                yaxisName_sigmaByMean = "#frac{#sigma(L1T/PF - 1)}{#mu(L1T/PF - 1) + 1}"
                                fGaus.SetRange(-1.4, 2.5)
                                axisRange_mean  = [-0.8, 1]
                                axisRange_sigma = [0.1, 0.9]
                            elif plot1DName == 'PU':
                                yaxisName_mean  = "#mu(L1 jet PU E_{T}) [GeV]"
                                yaxisName_sigma = "#sigma(L1 jet PU E_{T}) [GeV]"
                                fGaus.SetRange(0, 60)                        
                                axisRange_mean  = [0, 60]
                                axisRange_sigma = [0, 20]
                            elif plot1DName == 'PUByRawPt':
                                yaxisName_mean  = "#mu(#frac{L1 jet PU E_{T}}{L1 jet raw E_{T}}) [GeV]"
                                yaxisName_sigma = "#sigma(#frac{L1 jet PU E_{T}}{L1 jet raw E_{T}}) [GeV]"
                                fGaus.SetRange(0.0, 1.3)
                                axisRange_mean  = [0, 0.9]
                                axisRange_sigma = [0., 0.4]
                                #h2D.RebinY(2)

                            c2_1.cd()
                            h2D.FitSlicesY(fGaus, 1,h2D.GetNbinsX(), 1)


                            c2_1.cd()
                            #hres2D_1 = R.gDirectory.Get("%s_1" % (res_histoName))
                            h2D_1 = (R.gDirectory.Get("%s_1" % (sNameTmp)))
                            #h2D_1.SetNameTitle("%s_1" % (sNameTmp), "%s_1" % (sNameTmp))
                            h2D_1.SetNameTitle("%s_Mean" % (sNameTmp), "%s_Mean" % (sNameTmp))
                            h2D_1.SetLineColor(colors[iResType])
                            h2D_1.SetMarkerColor(colors[iResType])
                            #h2D_1.SetMarkerStyle(markers[iInFile])
                            h2D_1.SetMarkerStyle(markers[iPURange])
                            h2D_1.SetMarkerSize(0.7)
                            h2D_1.GetYaxis().SetTitle(yaxisName_mean)

                            c3_1.cd()
                            #hres2D_2 = R.gDirectory.Get("%s_2" % (res_histoName))
                            h2D_2 = ( R.gDirectory.Get("%s_2" % (sNameTmp)) )
                            #h2D_2.SetNameTitle("%s_2" % (sNameTmp), "%s_2" % (sNameTmp))
                            h2D_2.SetNameTitle("%s_Sigma" % (sNameTmp), "%s_Sigma" % (sNameTmp))
                            h2D_2.SetLineColor(colors[iResType])
                            h2D_2.SetMarkerColor(colors[iResType])
                            #h2D_2.SetMarkerStyle(markers[iInFile])
                            h2D_2.SetMarkerStyle(markers[iPURange])
                            h2D_2.SetMarkerSize(0.7)
                            h2D_2.GetYaxis().SetTitle(yaxisName_sigma)

                            c4_1.cd()
                            h2D_MeanPlus1    = h2D_1.Clone("%s_MeanPlus1" % (sNameTmp))
                            h2D_SF_PFJByL1TJ = h2D_1.Clone("%s_SF_PFJByL1TJ" % (sNameTmp))
                            h2D_SF_PFJByL1TJ.GetYaxis().SetTitle("#frac{1}{#mu(L1T/PF - 1) + 1}")
                            for iBin in range(1, h2D_MeanPlus1.GetNbinsX()+1):
                                binContent    = h2D_MeanPlus1.GetBinContent(iBin)
                                binError      = h2D_MeanPlus1.GetBinError(iBin)
                                MeanPlus1     = binContent + 1
                                binContent_SF = 1 / (MeanPlus1)
                                binError_SF   = binError / (MeanPlus1 * MeanPlus1)
                                h2D_MeanPlus1.SetBinContent(iBin,    MeanPlus1)
                                h2D_SF_PFJByL1TJ.SetBinContent(iBin, binContent_SF)
                                h2D_SF_PFJByL1TJ.SetBinError(iBin,   binError_SF)
                            print "h2D_SF_PFJByL1TJ: {}".format(h2D_SF_PFJByL1TJ.GetName())    
                            hEffectiveResolution = h2D_2.Clone("%s_EffectiveResolution" % (sNameTmp))
                            hEffectiveResolution.Divide(h2D_2, h2D_MeanPlus1)
                            hEffectiveResolution.GetYaxis().SetTitle(yaxisName_sigmaByMean)

                            if   plot1DName == 'res':
                                fOpL1JetEnergySF.cd();
                                h2D_SF_PFJByL1TJ.Write()

                            hres2D_dict[sInFileTag][jetShape_PUS_PUrange] = h2D
                            hresMean_dict[sInFileTag][jetShape_PUS_PUrange] = h2D_1
                            hresSigma_dict[sInFileTag][jetShape_PUS_PUrange] = h2D_2
                            hEffectiveResolution_dict[sInFileTag][jetShape_PUS_PUrange] = hEffectiveResolution

                            # make ratio plot of "highPU / lowPU"
                            if PURangeName not in JetShapesAndPUSs_Denom_forRatioPlot: # JetShapesAndPUSs_Denom_forRatioPlot contains label for denomintor histogram
                                jetShape_PUS_PUrange_forRatioPlot = "%s_%s" % (jetShape,JetShapesAndPUSs_Denom_forRatioPlot[0])
                                sYName  = "#frac{h(%s)}{h(%s)}" % (PURangeName, JetShapesAndPUSs_Denom_forRatioPlot[0])
                                sYName1 = "#frac{h(%s) + 1}{h(%s) + 1}" % (PURangeName, JetShapesAndPUSs_Denom_forRatioPlot[0])

                                hRatio2D_1       = cloneHistogram(hresMean_dict[sInFileTag][jetShape_PUS_PUrange],              "ratio")
                                hRatio2D_1_Denom = cloneHistogram(hresMean_dict[sInFileTag][jetShape_PUS_PUrange_forRatioPlot], "ratio")
                                # make "mean + 1" so as to get L1TJetEt scale factor
                                for iBin in range(1, hRatio2D_1.GetNbinsX()):
                                    hRatio2D_1      .SetBinContent(iBin, (hRatio2D_1      .GetBinContent(iBin) + 1) )
                                    hRatio2D_1_Denom.SetBinContent(iBin, (hRatio2D_1_Denom.GetBinContent(iBin) + 1) )
                                hRatio2D_1.Divide(hRatio2D_1_Denom)

                                hRatio2D_2       = cloneHistogram(hresSigma_dict[sInFileTag][jetShape_PUS_PUrange],              "ratio")
                                hRatio2D_2_Denom = cloneHistogram(hresSigma_dict[sInFileTag][jetShape_PUS_PUrange_forRatioPlot], "ratio")
                                hRatio2D_2.Divide(hRatio2D_2_Denom)

                                hRatio2D_3       = cloneHistogram(hEffectiveResolution_dict[sInFileTag][jetShape_PUS_PUrange],              "ratio")
                                hRatio2D_3_Denom = cloneHistogram(hEffectiveResolution_dict[sInFileTag][jetShape_PUS_PUrange_forRatioPlot], "ratio")
                                hRatio2D_3.Divide(hRatio2D_3_Denom)

                                hRatioresMean_dict[sInFileTag][jetShape_PUS_PUrange]             = hRatio2D_1
                                hRatioresSigma_dict[sInFileTag][jetShape_PUS_PUrange]            = hRatio2D_2
                                hRatioEffectiveResolution_dict[sInFileTag][jetShape_PUS_PUrange] = hRatio2D_3

                                hRatioresMean_dict[sInFileTag][jetShape_PUS_PUrange].GetYaxis().SetTitle(sYName1)
                                hRatioresSigma_dict[sInFileTag][jetShape_PUS_PUrange].GetYaxis().SetTitle(sYName)
                                hRatioEffectiveResolution_dict[sInFileTag][jetShape_PUS_PUrange].GetYaxis().SetTitle(sYName)

                            if iPURange == 0:
                                jetShape_PUS_PUrange_nice = "%s" % (jetShape)
                                jetShape_PUS_PUrange_nice = jetShape_PUS_PUrange_nice.replace('_plus_', ' + ')
                                jetShape_PUS_PUrange_nice = jetShape_PUS_PUrange_nice.replace('_times_', ' * ')
                                legend.AddEntry(hresSigma_dict[sInFileTag][jetShape_PUS_PUrange], jetShape_PUS_PUrange_nice, "lep")
                            iPURange += 1

                            c1.Update()
                            c2_1.Update()
                            c3_1.Update()
                            c4_1.Update()
                            time.sleep( 1 )

                        # Different colors for different JetSpaeAndPUS, but same color for PUrange
                        iResType += 1

                    '''            
                    shRatio_Denom_JetShape = "%s_%s_%s" % ( 
                        JetShapesAndPUSs_Denom_forRatioPlot[0],
                        JetShapesAndPUSs_Denom_forRatioPlot[1],
                        JetShapesAndPUSs_Denom_forRatioPlot[2],
                    )
                    '''




                    hresMean_list = []
                    hresSigma_list = []
                    hEffectiveResolution_list = []
                    hRatioresMean_list = []
                    hRatioresSigma_list = []
                    hRatioEffectiveResolution_list = []
                    #for jetShape in JetShapes:
                    for jetShape in hEffectiveResolution_dict[sInFileTag].keys():
                        hresMean_list.append(hresMean_dict[sInFileTag][jetShape])
                        hresSigma_list.append(hresSigma_dict[sInFileTag][jetShape])
                        hEffectiveResolution_list.append(hEffectiveResolution_dict[sInFileTag][jetShape])
                        
                    for jetShape in hRatioEffectiveResolution_dict[sInFileTag].keys():
                        hRatioresMean_list.append(hRatioresMean_dict[sInFileTag][jetShape])
                        hRatioresSigma_list.append(hRatioresSigma_dict[sInFileTag][jetShape])
                        hRatioEffectiveResolution_list.append(hRatioEffectiveResolution_dict[sInFileTag][jetShape])
                        
                        '''
                        if jetShape == shRatio_Denom_JetShape: continue
                        
                        hRatioresMean_dict[sInFileTag][jetShape] = hresMean_dict[sInFileTag][jetShape].Clone("%s_Ratio" % hresMean_dict[sInFileTag][jetShape].GetName())
                        hRatioresSigma_dict[sInFileTag][jetShape] = hresSigma_dict[sInFileTag][jetShape].Clone("%s_Ratio" % hresSigma_dict[sInFileTag][jetShape].GetName())
                        hRatioEffectiveResolution_dict[sInFileTag][jetShape] = hEffectiveResolution_dict[sInFileTag][jetShape].Clone("%s_Ratio" % hEffectiveResolution_dict[sInFileTag][jetShape].GetName())
                        
                        hRatioresMean_dict[sInFileTag][jetShape].Divide(hresMean_dict[sInFileTag][jetShape], hresMean_dict[sInFileTag][shRatio_Denom_JetShape])
                        hRatioresSigma_dict[sInFileTag][jetShape].Divide(hresSigma_dict[sInFileTag][jetShape], hresSigma_dict[sInFileTag][shRatio_Denom_JetShape])
                        hRatioEffectiveResolution_dict[sInFileTag][jetShape].Divide(hEffectiveResolution_dict[sInFileTag][jetShape], hEffectiveResolution_dict[sInFileTag][shRatio_Denom_JetShape])
                        
                        sYName = "#frac{h(i)}{h(%s)}" % (shRatio_Denom_JetShape)
                        hRatioresMean_dict[sInFileTag][jetShape].GetYaxis().SetTitle(sYName)
                        hRatioresSigma_dict[sInFileTag][jetShape].GetYaxis().SetTitle(sYName)
                        hRatioEffectiveResolution_dict[sInFileTag][jetShape].GetYaxis().SetTitle(sYName)
                        
                        hRatioresMean_list.append(hRatioresMean_dict[sInFileTag][jetShape])
                        hRatioresSigma_list.append(hRatioresSigma_dict[sInFileTag][jetShape])
                        hRatioEffectiveResolution_list.append(hRatioEffectiveResolution_dict[sInFileTag][jetShape])
                        '''
                    
                    
                    print "here1"
                    print "iInFile: {}".format(iInFile)
                    print "sInFileTags[iInFile]: {}".format(sInFileTags[iInFile])
                    print "sResLegend_selected: {}".format(sResLegend_selected)
                    sDir1 = "plots%s/CompareDiffJetShapesAndPUS/%s/%s/%s" % (sInFileVersion, plot1DName, resPtBin, sInFileTags[iInFile])
                    if not os.path.exists(sDir1):
                        os.makedirs(sDir1)                                
                    #c2_1 = plotHistos1DAndRatioPlot(hresMean_list, hRatioresMean_list, legend, c2_1, sSaveAs="%s/L1TauMatchedToPFJet_%s_%s_mean.png" % (sDir1, l1Mode, plotName), yRange=[-0.6, 0.4])
                    c2_1 = plotHistos1DAndRatioPlot(hresMean_list, hRatioresMean_list, legend, c2_1, sSaveAs="%s/L1TauMatchedToPFJet_%s_%s_mean.png" % (sDir1, l1Mode, plotName))
                    c3_1 = plotHistos1DAndRatioPlot(hresSigma_list, hRatioresSigma_list, legend, c3_1, sSaveAs="%s/L1TauMatchedToPFJet_%s_%s_sigma.png" % (sDir1, l1Mode, plotName))
                    c4_1 = plotHistos1DAndRatioPlot(hEffectiveResolution_list, hRatioEffectiveResolution_list, legend, c4_1, sSaveAs="%s/L1TauMatchedToPFJet_%s_%s_effectiveResolution.png" % (sDir1, l1Mode, plotName))
                    
                    '''
                    c1.Update()
                    c2_1.Update()
                    c3_1.Update()
                    '''
                    #c2_1.SaveAs("%s/DiffJetShapes_%s_%s_mean.png" % (sDir1, l1Mode, plotName))
                    #c3_1.SaveAs("%s/DiffJetShapes_%s_%s_sigma.png" % (sDir1, l1Mode, plotName))
        
        
        fOpL1JetEnergySF.Close()            




















        
if __name__ == '__main__':

    #resPtBins = ["lowPt"] # "medPt"  "PtAllBins"  "modPt"
    resPtBins = ["Pt60To90"] # "medPt"  "PtAllBins"  "modPt", "Pt25To35", "Pt35To60"
    
    for resPtBin in resPtBins:
        run(resPtBin)
