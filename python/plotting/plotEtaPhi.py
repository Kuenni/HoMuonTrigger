#!/usr/bin/python
from ROOT import ROOT,gROOT,gStyle,TCanvas,TFile,TH1D,TH2D,TLegend,THStack,TPaveText,TBox
import sys
import os
import PlotStyle
import numpy as np

PlotStyle.setPlotStyle()

DEBUG = 1
prefix = '[plotEtaPhi] '


def drawHoBoxes(canvas):
    canvas.cd()
    boxes = []
    for i in range(1,8):
        for j in range(1,8):
            box = TBox(-0.261 + (i-1)*0.087,-0.261 + (j-1)*0.087,-0.261 + i*0.087,-0.261 + j*0.087)
            box.SetFillStyle(0)
            box.SetLineColor(ROOT.kBlack)
            box.SetLineWidth(1)
            box.Draw()
            boxes.append(box)
    return boxes

def plotEtaPhi(folder):
	
	if(DEBUG):
		print prefix + 'was called'
	
	if(folder == None):
		print prefix + 'Error! Filename as first argument needed.'
		return	
	if( not os.path.exists('plots')):
	    os.mkdir('plots')
   	if( not os.path.exists('plots/' + folder)):
		os.mkdir('plots/' + folder)
	
	filename = folder + '/L1MuonHistogram.root'
	if( not os.path.exists(filename)):
		print 'Error! File ' + filename + ' does not exist!'
		return
	print prefix + 'Opening file:',filename
	
	file = TFile.Open(filename)
	
	
	h2dDeltaEtaDeltaPhi = file.Get("hoMuonAnalyzer/etaPhi/horeco_EtaPhi")
	
	canv = TCanvas("canvasEtaPhi",'canvasEtaPhi',1200,1200)
	canv.SetLogz()
	#h2dDeltaEtaDeltaPhi.Rebin2D(8,8)
	#h2dDeltaEtaDeltaPhi.GetXaxis().SetRangeUser(-.45,.45)
	h2dDeltaEtaDeltaPhi.GetXaxis().SetTitle("#Delta#eta")
	#h2dDeltaEtaDeltaPhi.GetYaxis().SetRangeUser(-.45,.45)
	h2dDeltaEtaDeltaPhi.GetYaxis().SetTitle("#Delta#phi")
	h2dDeltaEtaDeltaPhi.GetZaxis().SetTitle("N")
	h2dDeltaEtaDeltaPhi.Draw("colz")
	
	boxList = drawHoBoxes(canv)
	
	canv.Update()
	pal = h2dDeltaEtaDeltaPhi.GetListOfFunctions().FindObject("palette")
	pal.SetX2NDC(0.92)
	
	stats = h2dDeltaEtaDeltaPhi.GetListOfFunctions().FindObject("stats")
	stats.SetX1NDC(.1)
	stats.SetX2NDC(.2)
	stats.SetY1NDC(.1)
	stats.SetY2NDC(.25)
	
	legend = TLegend(0.7,0.8,0.9,0.9)
	legend.AddEntry(boxList[0],"HO tile dimensions","le")
	legend.Draw()
	
	canv.Update();
	
	canv.SaveAs("plots/" + folder + "/etaPhi.png")
	canv.SaveAs("plots/" + folder + "/etaPhi.pdf")
	
	f = TFile.Open("plots/" + folder + "/etaPhi.root","RECREATE")
	canv.Write()
	f.Close()
