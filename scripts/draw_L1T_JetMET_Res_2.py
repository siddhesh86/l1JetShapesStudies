import os
import sys
from collections import OrderedDict as OD

import ROOT as R
R.gROOT.SetBatch(False)  ## Don't print histograms to screen while processing


sIpFile = 'L1T_HCALL2Calib_stage1_l1NtupleChunkyDonut_PFA1p_nVtxAll_QCDMC_hadded_JEC2022v5.root'; sRefJet = 'GEN'
sOutDir = 'plot_JEC2022v5_QCD'

PT_CAT = OD()
PT_CAT['Ptlt25']   = [ 0,  15,   25]  ## Low pT, turn-on threshold, high pT
PT_CAT['Pt25To35'] = [25,  30,   35]  ## Low pT, turn-on threshold, high pT
PT_CAT['Pt35To60'] = [35,  55,   60]  ## Low pT, turn-on threshold, high pT
PT_CAT['Pt60To90'] = [60,  75,   90]  ## Low pT, turn-on threshold, high pT   #[60,  90,   90]
PT_CAT['Ptgt90']   = [90, 120, 9999]  ## Low pT, turn-on threshold, high pT
sPtAll = 'PtAll'

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

sHistograms = OD([
    ('L1JetEt / %sJetEt - 1'%(sRefJet), 'h_jet_res_PUS_$ETACAT_$PTCAT_0_emu'),
])







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





if __name__ == '__main__':

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



    colors   = [R.kRed, R.kBlue, R.kMagenta, R.kCyan+1, R.kOrange+2, R.kGreen+2, R.kGray+22, R.kTeal, R.kViolet, R.kSpring, R.kBlack, ]
    markers  = [20, 24, 23, 20]
    colors2  = [R.kBlack, R.kRed, R.kBlue, R.kMagenta, R.kCyan+1, R.kOrange+2, R.kGray+2, R.kGreen+2, R.kBlue-10, R.kViolet, R.kTeal, R.kSpring, R.kBlack, ]
    markers2 = [20, 22, 24, 23, 21]



    fIp = R.TFile(sIpFile)
    print("Input file %s " % (sIpFile))
    if not fIp.IsOpen():
        print("Input file %s couldn't open." % (sIpFile))
        exit(0)


    if not os.path.exists(sOutDir):
        os.mkdir(sOutDir)
    

    for sHistoName_nice, sHistoName in sHistograms.items():
        histos = OD()
        for eta_cat in ETA_CAT.keys():
            if eta_cat == 'HBEF': continue

            histos[eta_cat] = OD()
            for pt_cat in PT_CAT.keys():
                sHistoName_toUse = sHistoName
                sHistoName_toUse = sHistoName_toUse.replace('$ETACAT', eta_cat)
                sHistoName_toUse = sHistoName_toUse.replace('$PTCAT',  pt_cat)
                print(f"histogram {sHistoName_toUse}")
                histo = fIp.Get(sHistoName_toUse)
                #print(f"sHisto {sHistoName_toUse}, {type(histo)}: {histo}")
                histos[eta_cat][pt_cat] = histo

            if sPtAll not in PT_CAT.keys():
                histos[eta_cat][sPtAll] = None
                for pt_cat in PT_CAT:
                    if histos[eta_cat][sPtAll] is None: histos[eta_cat][sPtAll] = histos[eta_cat][pt_cat].Clone('%s_clone'%(histos[eta_cat][pt_cat].GetName()))
                    else:                               histos[eta_cat][sPtAll].Add( histos[eta_cat][pt_cat] )


        # normalize all histograms
        for i_eta_cat, eta_cat in enumerate(list(histos.keys())):
            for i_pt_cat, pt_cat in enumerate(list(histos[eta_cat].keys())):
                h = histos[eta_cat][pt_cat]

                nEntries = h.Integral()
                scaleFactor = 1./nEntries
                #print(f"{eta_cat} {pt_cat}: nEntries {nEntries}, scaleFactor {scaleFactor}, h.GetMaximum {h.GetMaximum()}")
                h.Scale( scaleFactor )

                
        # plot different pt_cat for each eta_cat
        for i_eta_cat, eta_cat in enumerate(list(histos.keys())):
            c1 = R.TCanvas("c1","c1",600,500)
            legend = R.TLegend(0.3,0.75,1,1)
            legend.SetHeader('%s (%.2f < #eta < %.2f)' % (eta_cat, ETA_CAT[eta_cat][0],ETA_CAT[eta_cat][1]), 'C')

            for i_pt_cat, pt_cat in enumerate(list(histos[eta_cat].keys())):
                h = histos[eta_cat][pt_cat]

                '''
                nEntries = h.Integral()
                scaleFactor = 1./nEntries
                print(f"{eta_cat} {pt_cat}: nEntries {nEntries}, scaleFactor {scaleFactor}, h.GetMaximum {h.GetMaximum()}")
                h.Scale( scaleFactor )
                '''

                fGaus = R.TF1('gaus','gaus', -1, 1)
                h.Fit(fGaus, 'R0')
                peak_mu        = fGaus.GetParameter(1)
                peak_mu_err    = fGaus.GetParError(1)
                peak_sigma     = fGaus.GetParameter(2)
                peak_sigma_err = fGaus.GetParError(2)

                h.GetXaxis().SetTitle(sHistoName_nice)
                h.GetYaxis().SetTitle("Entries")
                h.GetXaxis().SetTitleOffset(1.3)
                h.SetLineColor(colors[i_pt_cat])
                h.SetMarkerColor(colors[i_pt_cat])
                h.SetMarkerStyle(markers[0])
                h.SetMarkerSize(0.7)

                if i_pt_cat == 0:
                    h_group = list(list(histos[eta_cat].values()))
                    ymin, ymax = getHists1DYRange(h_group)
                    #ymin *= scaleFactor
                    #ymax *= scaleFactor
                    #print("ymin {}, ymax {}, h_group ({}) : {}".format(ymin, ymax, len(h_group), h_group))
                    h.GetYaxis().SetRangeUser(0, ymax * 1.1)
                    #h.GetYaxis().SetRangeUser(0, 0.25)
                    h.Draw()
                else:
                    h.Draw('same')

                legend.AddEntry(h, '%s %-10s #mu=%4.2f #pm %.2f, #sigma=%.2f #pm %.2f'%(sRefJet, pt_cat, peak_mu,peak_mu_err, peak_sigma,peak_sigma_err), "lep")
                    
            legend.Draw()            
            c1.Update()
            c1.SaveAs('%s/L1Jet_res_%s.png' % (sOutDir, eta_cat))



        # plot different eta_cat
        c1 = R.TCanvas("c1","c1",600,500)
        legend = R.TLegend(0.3,0.75,1,1)
        legend.SetHeader('%s' % (sPtAll), 'C')
        for i_eta_cat, eta_cat in enumerate(list(histos.keys())):
            h = histos[eta_cat][sPtAll]

            fGaus = R.TF1('gaus','gaus', -1, 1)
            h.Fit(fGaus, 'R0')
            peak_mu        = fGaus.GetParameter(1)
            peak_mu_err    = fGaus.GetParError(1)
            peak_sigma     = fGaus.GetParameter(2)
            peak_sigma_err = fGaus.GetParError(2)

            h.GetXaxis().SetTitle(sHistoName_nice)
            h.GetYaxis().SetTitle("Entries")
            h.GetXaxis().SetTitleOffset(1.3)
            h.SetLineColor(colors[i_eta_cat])
            h.SetMarkerColor(colors[i_eta_cat])
            h.SetMarkerStyle(markers[0])
            h.SetMarkerSize(0.7)

            if i_eta_cat == 0:
                h_group = [histos[eta_cat_1][sPtAll] for eta_cat_1 in histos.keys()]
                list(list(histos[eta_cat].values()))
                ymin, ymax = getHists1DYRange(h_group)
                #ymin *= scaleFactor
                #ymax *= scaleFactor
                #print("ymin {}, ymax {}, h_group ({}) : {}".format(ymin, ymax, len(h_group), h_group))
                h.GetYaxis().SetRangeUser(0, ymax * 1.1)
                #h.GetYaxis().SetRangeUser(0, 0.25)
                h.Draw()
            else:
                h.Draw('same')

            legend.AddEntry(h, '%.2f < #eta < %.2f, #mu=%4.2f #pm %.2f, #sigma=%.2f #pm %.2f'%(ETA_CAT[eta_cat][0],ETA_CAT[eta_cat][1], peak_mu,peak_mu_err, peak_sigma,peak_sigma_err), "lep")

        legend.Draw()            
        c1.Update()
        c1.SaveAs('%s/L1Jet_res_HBEF.png' % (sOutDir))


