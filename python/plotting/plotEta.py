#!/usr/bin/python
from ROOT import ROOT,gROOT,gStyle,TCanvas,TFile,TH1D,TH2D,TLegend,THStack,TPaveText,TBox
import sys
import os
import PlotStyle
import numpy as np

PlotStyle.setPlotStyle()
if(len(sys.argv) < 2):
    print 'Error! Filename as first argument needed.'
    sys.exit(1)

if( not os.path.exists('plots')):
    os.mkdir('plots')

filename = sys.argv[1]
print 'Opening file:',filename

file = TFile.Open(filename)

h1eta = file.Get("hoMuonAnalyzer/etaPhi/horeco_Eta")
h1MatchedEtaHo = file.Get("hoMuonAnalyzer/etaPhi/L1MuonWithHoMatch_Eta")

canv = TCanvas("canvasEta",'canvasEta',1200,1200)
canv.SetLogy()

h1eta.GetXaxis().SetTitle("#eta")
h1eta.GetYaxis().SetTitle("N")
h1eta.Rebin(8)
h1eta.GetYaxis().SetRangeUser(0.1,3e5)

h1eta.Draw()


h1MatchedEtaHo.SetLineColor(ROOT.kRed)
h1MatchedEtaHo.Rebin(4)
h1MatchedEtaHo.Draw('Same')

canv.Update()

stats = h1eta.GetListOfFunctions().FindObject("stats")
stats.SetX1NDC(.1)
stats.SetX2NDC(.3)
stats.SetY1NDC(.9)
stats.SetY2NDC(1.)


canv.Update()

canv.SaveAs("plots/eta.png")
canv.SaveAs("plots/eta.pdf")

f = TFile.Open("plots/eta.root","RECREATE")
canv.Write()
f.Close()
