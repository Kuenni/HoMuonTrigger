#!/usr/bin/python
import os,sys
from math import sqrt
from ROOT import TCanvas,ROOT,TFile,TLegend,TF1,TLine,gROOT,TPaveText,TH1D,Double,TH2D,THStack,gStyle
from plotting.PlotStyle import setPlotStyle,calcSigma,getLabelCmsPrivateSimulation,colorRwthDarkBlue,colorRwthMagenta,setupAxes,convertToHcalCoords,chimney1,chimney2
from plotting.RootFileHandler import RootFileHandler

DEBUG = 1

gROOT.Reset()

setPlotStyle()

if len(sys.argv) < 2:
	print 'First argument has to be the file name scheme!'
fileHandler = RootFileHandler(sys.argv[1])
fileHandler.printStatus()

if( not os.path.exists('plots')):
	os.mkdir('plots')
if( not os.path.exists('plots/timing')):
	os.mkdir('plots/timing')

# Plot the delta timing distribution for Ho
# and the L1MuonObject
filename = 'L1MuonHistogram.root'
if(DEBUG):
	print 'Opening file:',filename
file = TFile.Open(filename)
if(file == None):
	print 'Error opening file:',filename

def plotDeltaTime():
	hDeltaTAllHo = fileHandler.getHistogram('hoMuonAnalyzer/L1MuonPresentHoMatch_DeltaTime')
	hDeltaTCleanHo = fileHandler.getHistogram('hoMuonAnalyzer/L1MuonAboveThr_DeltaTime')
	
	c = TCanvas("c","Delta Time",1200,1200)
	c.SetLogy()
	
	hDeltaTAllHo.SetLineColor(ROOT.kBlue)
	hDeltaTAllHo.SetLineWidth(3)
	hDeltaTAllHo.SetFillColor(ROOT.kBlue)
	hDeltaTAllHo.SetFillStyle(3017)
	hDeltaTAllHo.SetTitle("#Delta time")
	hDeltaTAllHo.SetStats(0)
	
	hDeltaTCleanHo.SetLineColor(8)
	hDeltaTCleanHo.SetFillColor(8)
	hDeltaTCleanHo.SetLineWidth(3)
	hDeltaTCleanHo.SetFillStyle(3002)
	
	#hDeltaTAllHo.Scale(1/hDeltaTAllHo.Integral())
	#hDeltaTCleanHo.Scale(1/hDeltaTCleanHo.Integral())
	
	print hDeltaTCleanHo.Integral(),hDeltaTAllHo.Integral()
	
	fitFirstMin = TF1("fitFirstMin","[0]+x*[1]+[2]*x**2")
	fitSecondMin = TF1("fitsecondMin","[0]+x*[1]+[2]*x**2",10,20)
	
	hDeltaTCleanHo.Fit(fitFirstMin,"+q","",-20,-10)
	hDeltaTCleanHo.Fit(fitSecondMin,"R+q","")
	
	hDeltaTAllHo.Draw()
	legend = TLegend(0.6,0.75,0.9,0.9)
	legend.AddEntry(hDeltaTAllHo,"L1Muon matched to any HO","le")
	legend.Draw()
	
	label = getLabelCmsPrivateSimulation()
	label.Draw()
	c.Update()
	
	c.SaveAs("plots/timing/deltaTimeAllHo.png")
	c.SaveAs("plots/timing/deltaTimeAllHo.pdf")
	
	hDeltaTCleanHo.Draw('same')
	
	fitFirstMin.SetRange(-50,50)
	fitSecondMin.SetRange(-50,50)
	
	#fitFirstMin.Draw('lSame')
	#fitSecondMin.Draw('lSame')
	
	lineFirstMin = TLine(fitFirstMin.GetMinimumX(-20,-10),hDeltaTAllHo.GetMinimum(),fitFirstMin.GetMinimumX(-20,-10),hDeltaTAllHo.GetMaximum())
	lineFirstMin.SetLineWidth(3)
	lineFirstMin.SetLineColor(ROOT.kRed)
	lineFirstMin.Draw()
	
	lineSecondMin = TLine(fitSecondMin.GetMinimumX(10,20),hDeltaTAllHo.GetMinimum(),fitSecondMin.GetMinimumX(10,20),hDeltaTAllHo.GetMaximum())
	lineSecondMin.SetLineWidth(3)
	lineSecondMin.SetLineColor(ROOT.kRed)
	lineSecondMin.Draw()
	
	
	legend.AddEntry(hDeltaTCleanHo,"L1Muon matched to HO > 0.2 GeV","le")
	legend.AddEntry(lineFirstMin,"Integral boundaries","e")
	legend.Draw()
	
	integralCenter = hDeltaTCleanHo.Integral(hDeltaTCleanHo.FindBin(fitFirstMin.GetMinimumX(-20,-10)),hDeltaTCleanHo.FindBin(fitSecondMin.GetMinimumX(10,20)))
	integralCenterAll = hDeltaTAllHo.Integral(hDeltaTAllHo.FindBin(fitFirstMin.GetMinimumX(-20,-10)),hDeltaTAllHo.FindBin(fitSecondMin.GetMinimumX(10,20)))
	print 80*'#'
	print 'Integral of center area in clean histogram :',integralCenter
	print '==> %.2f%% +/- %.2f%%' % (integralCenter/hDeltaTCleanHo.Integral()*100,calcSigma(integralCenter, hDeltaTCleanHo.Integral())*100)
	print 'Integral of center area in all matched HO events:',integralCenterAll
	print '==> %.2f%% +/- %.2f%%' % (integralCenterAll/hDeltaTAllHo.Integral()*100,calcSigma(integralCenterAll, hDeltaTAllHo.Integral())*100)
	print 80*'#'
	
	paveText = TPaveText(0.6,0.7,0.9,0.75,'NDC')
	paveText.AddText('%s' % ('Central peak contains (filtered hist.)'))
	paveText.AddText('%.2f%% +/- %.2f%%' % (integralCenter/hDeltaTCleanHo.Integral()*100,calcSigma(integralCenter, hDeltaTCleanHo.Integral())*100))
	paveText.SetBorderSize(1)
	paveText.Draw()
	
	label = getLabelCmsPrivateSimulation()
	label.Draw()
	c.Update()
	
	
	c.SaveAs("plots/timing/deltaTime.png")
	c.SaveAs("plots/timing/deltaTime.pdf")

def plotL1BxId():
	c2 = TCanvas("c2","BX ID",1200,1200)
	c2.SetLogy()
	
	hBxId = fileHandler.getHistogram('hoMuonAnalyzer/L1MuonPresent_BxId')
	hBxIdAboveThr = fileHandler.getHistogram('hoMuonAnalyzer/L1MuonAboveThr_BxId')
	hBxId.SetLineColor(colorRwthDarkBlue)
	hBxId.SetLineWidth(3)
	hBxId.SetStats(0)
	hBxId.SetTitle("BX ID distribution")
	hBxId.GetXaxis().SetRangeUser(-2,5)
	#hBxId.SetFillColor(colorRwthDarkBlue)
	#hBxId.SetFillStyle(3017)
	hBxIdAboveThr.SetFillStyle(3002)
	hBxIdAboveThr.SetFillColor(colorRwthMagenta)
	hBxIdAboveThr.SetLineColor(colorRwthMagenta)
	hBxIdAboveThr.SetLineWidth(2)
	
	hBxId.Scale(1/hBxId.Integral())
	hBxIdAboveThr.Scale(1/hBxIdAboveThr.Integral())
	
	hBxId.GetYaxis().SetTitle("rel. fraction")
	
	hBxId.Draw()
	#hBxIdAboveThr.Draw("same")
	
	legend2 = TLegend(0.45,0.8,0.9,0.9)
	legend2.AddEntry(hBxId,"L1Muon","f")
	#legend2.AddEntry(hBxIdAboveThr,"L1Muon matched to HO > 0.2 GeV","f")
	legend2.Draw()
	
	label = getLabelCmsPrivateSimulation()
	label.Draw()
	
	c2.SaveAs("plots/timing/bxId.png")
	c2.SaveAs("plots/timing/bxId.pdf")
	
	c4 = TCanvas("c4","BX L1Muon in ns",1200,1200)
	c4.SetLogy()
	hL1InNs = TH1D("hL1InNs","L1 Muon time in ns;time / ns;#",201,-100.5,100.5)
	for i in range(0,hBxId.GetNbinsX()):
		x = hBxId.GetBinCenter(i)*25
		y = hBxId.GetBinContent(i)
		hL1InNs.SetBinContent(hL1InNs.FindBin(x),y)
	hL1InNs.SetStats(0)
	hL1InNs.Draw()
	c4.SaveAs("plots/timing/timeL1Only.pdf")

def plotHoTimeLog():
	c3 = TCanvas("c3Log","HO Time Log",1200,1200)
	skipNoisePlot = True
	if not skipNoisePlot:
		c3.Divide(1,2)
		c3.cd(1).SetLogy()
		label = getLabelCmsPrivateSimulation()
		label.Draw()
		hHoTime = fileHandler.getHistogram('hoMuonAnalyzer/hoRecHits_Time')
		
		hHoTime.SetStats(0)
		hHoTime.SetTitle("Time distribution for all HO Rec Hits")
		hHoTime.SetLineColor(colorRwthDarkBlue)
		hHoTime.SetLineWidth(3)
		hHoTime.Draw()
		label = getLabelCmsPrivateSimulation()
		label.Draw()
		
	hHoTimeAboveThr = fileHandler.getHistogram('hoMuonAnalyzer/hoRecHitsAboveThr_Time')
	c3.cd(2).SetLogy()
	hHoTimeAboveThr.SetStats(0)
	hHoTimeAboveThr.SetTitle("Time distribution for HO Rec Hits > 0.2 GeV")
	hHoTimeAboveThr.SetLineColor(colorRwthDarkBlue)
	hHoTimeAboveThr.SetLineWidth(3)
	setupAxes(hHoTimeAboveThr)
	hHoTimeAboveThr.Draw()
	
	label = getLabelCmsPrivateSimulation()
	label.Draw()

	fit = TF1("fit","gaus",-10,10)
	fit.SetParameter(1,0)
	fit.SetParameter(2,1)
	hHoTimeAboveThr.Fit(fit,'','R',-12.5,12.5)
	
	pText = TPaveText(0.7,0.8,0.9,0.9,'NDC')
	pText.AddText('Mean: %.2f ns' % (fit.GetParameter(1)))
	pText.AddText('#sigma: %.2f ns' % (fit.GetParameter(2)))
	pText.SetBorderSize(1)
	pText.SetFillColor(0)
	pText.Draw()
	
	c3.Update()
	c3.SaveAs("plots/timing/hoTimeLog.png")
	c3.SaveAs("plots/timing/hoTimeLog.pdf")
	
	return c3,pText,hHoTimeAboveThr

def plotHoTime():
	c3 = TCanvas("c3","HO Time",1200,1200)
	skipNoisePlot = True
	if not skipNoisePlot:
		c3.Divide(1,2)
		c3.cd(1).SetLogy()
		label = getLabelCmsPrivateSimulation()
		label.Draw()
		hHoTime = fileHandler.getHistogram('hoMuonAnalyzer/hoRecHits_Time')
		
		hHoTime.SetStats(0)
		hHoTime.SetTitle("Time distribution for all HO Rec Hits")
		hHoTime.SetLineColor(colorRwthDarkBlue)
		hHoTime.SetLineWidth(3)
		hHoTime.Draw()
		label = getLabelCmsPrivateSimulation()
		label.Draw()
		
	hHoTimeAboveThr = fileHandler.getHistogram('hoMuonAnalyzer/hoRecHitsAboveThr_Time')
	c3.cd(2).SetLogy()
	hHoTimeAboveThr.SetStats(0)
	hHoTimeAboveThr.SetTitle("Time distribution for HO Rec Hits > 0.2 GeV")
	hHoTimeAboveThr.SetLineColor(colorRwthDarkBlue)
	hHoTimeAboveThr.SetLineWidth(3)
	setupAxes(hHoTimeAboveThr)
	hHoTimeAboveThr.Draw()
	
	label = getLabelCmsPrivateSimulation()
	label.Draw()
	
	print 80*'#'
	print 'Integral of HO > 0.2 GeV time histogram:'
	print hHoTimeAboveThr.Integral()
	print
	
	xLow = -5
	xHigh = 5
	histogramBetween = hHoTimeAboveThr.Integral(hHoTimeAboveThr.FindBin(xLow),hHoTimeAboveThr.FindBin(xHigh))
	histogramTotal = float(hHoTimeAboveThr.Integral())
	print 'Histogram integral between %.f ns and %.f ns' % (xLow,xHigh)
	print '%d/%d => %.2f +/- %f' % (histogramBetween,histogramTotal,histogramBetween/histogramTotal,calcSigma(histogramBetween, histogramTotal))  
	print 80*'#'
	
	fit = TF1("fit","gaus",-10,10)
	hHoTimeAboveThr.Fit(fit)
	
	pText = TPaveText(0.7,0.8,0.9,0.9,'NDC')
	pText.AddText('Mean: %.2f ns' % (fit.GetParameter(1)))
	pText.AddText('#sigma: %.2f ns' % (fit.GetParameter(2)))
	pText.SetBorderSize(1)
	pText.SetFillColor(0)
	pText.Draw()
	
	c3.Update()
	c3.SaveAs("plots/timing/hoTime.png")
	c3.SaveAs("plots/timing/hoTime.pdf")
	
	return c3,label,hHoTimeAboveThr,pText

def plotFractionsOfBxId():
	##BX right plotting pt
	canvasBxRightPt = TCanvas("cBxRightPt","cBxRightPt",1200,1200)
	canvasBxRightPt.cd().SetLeftMargin(0.15)
	hBxRightPt = fileHandler.getHistogram('hoMuonAnalyzer/BxRightGen_Pt').Clone()
	setPlotStyle()
	hBxRightPt.Rebin(50)
	hBxRightPt.GetXaxis().SetRangeUser(0,200)
	hBxRightPt.GetYaxis().SetTitle("normalized Entries / 5 GeV")
	hBxRightPt.GetXaxis().SetTitle("p_{T} Gen")
	hBxRightPt.GetYaxis().SetTitleOffset(2)
	hBxRightPt.SetTitle("Events with right BX ID vs. p_{T}")
	hBxRightPt.SetStats(0)
	hBxRightPt.SetLineWidth(2)
	hBxRightPt.Scale(1/hBxRightPt.Integral())
	hBxRightPt.Draw()
	label = getLabelCmsPrivateSimulation()
	label.Draw()
	canvasBxRightPt.SaveAs('plots/timing/bxRightPt.pdf')
	
	##BX wrong plotting pt
	canvasBxWrongPt = TCanvas("cBxWrongPt","cBxWrongPt",1200,1200)
	canvasBxWrongPt.cd().SetLeftMargin(0.15)
	hBxWrongPt = fileHandler.getHistogram('hoMuonAnalyzer/BxWrongGen_Pt').Clone()
	setPlotStyle()
	hBxWrongPt.Rebin(50)
	hBxWrongPt.GetXaxis().SetRangeUser(0,200)
	hBxWrongPt.GetYaxis().SetTitle("normalized Entries / 5 GeV")
	hBxWrongPt.GetXaxis().SetTitle("p_{T} Gen")
	hBxWrongPt.GetYaxis().SetTitleOffset(2)
	hBxWrongPt.SetTitle("Events with wrong BX ID vs. p_{T}")
	hBxWrongPt.SetStats(0)
	hBxWrongPt.SetLineWidth(2)
	hBxWrongPt.Scale(1/hBxWrongPt.Integral())
	hBxWrongPt.DrawCopy()
	label = getLabelCmsPrivateSimulation()
	label.Draw()
	canvasBxWrongPt.SaveAs('plots/timing/bxWrongPt.pdf')
	
	#Plot the histogram stack
	canvasStack = TCanvas("cStacked","cStacked",1200,1200)
	canvasStack.cd().SetLeftMargin(0.15)
	hWrong = fileHandler.getHistogram('hoMuonAnalyzer/BxWrongGen_Pt')
	hRight = fileHandler.getHistogram('hoMuonAnalyzer/BxRightGen_Pt')
	hRightFraction = TH1D('hRightFraction','',100,0,500)
	hWrongFraction = TH1D('hWrongFraction','',100,0,500)
	hWrong.Rebin(50)
	hRight.Rebin(50)
	#Fill the histograms with the bin wise fractions
	for i in range(0,hRight.GetNbinsX()):
		nRight = hRight.GetBinContent(i+1)
		nWrong = hWrong.GetBinContent(i+1)
		if(nRight + nWrong == 0):
			continue
		hRightFraction.SetBinContent(i+1,nRight/(nRight+nWrong))
		hWrongFraction.SetBinContent(i+1,nWrong/(nRight+nWrong))
	
	#Create the stack
	stack = THStack("hstack","Fractions of events for BX ID correct and wrong")
	nRight = hRight.Integral()
	nWrong = hWrong.Integral()
	nTotal = nRight + nWrong
	hRightFraction.SetLineColor(colorRwthDarkBlue)
	hRightFraction.SetFillColor(colorRwthDarkBlue)
	hRightFraction.SetFillStyle(3002)
	hWrongFraction.SetLineColor(colorRwthMagenta)
	hWrongFraction.SetFillColor(colorRwthMagenta)
	hWrongFraction.SetFillStyle(3002)
	stack.Add(hRightFraction)
	stack.Add(hWrongFraction)
	stack.Draw()
	stack.GetXaxis().SetRangeUser(0,201)
	stack.GetYaxis().SetTitle('rel. fraction / 5 GeV')
	stack.GetYaxis().SetTitleOffset(2)
	stack.GetXaxis().SetTitle('p_{T} Gen')
	stack.SetMinimum(0.9)
	stack.SetMaximum(1)
	legend = TLegend(0.6,0.75,0.9,0.9)
	legend.AddEntry(hRightFraction,"BX ID right","f")
	legend.AddEntry(hWrongFraction,"BX ID wrong","f")
	legend.Draw()
	label = getLabelCmsPrivateSimulation()
	label.Draw()
	canvasStack.Update()
	canvasStack.SaveAs('plots/timing/bxStacked.pdf')
	canvasStack.SaveAs('plots/timing/bxStacked.root')

#Plot eta and phi of wrong bx ids
def plotEtaPhiOfWrongBxId():
	canvasEtaPhiBxWrong = TCanvas("canvasEtaPhiBxWrong","canvasEtaPhiBxWrong",1200,1200)
	
	hWrong = fileHandler.getHistogram('hoMuonAnalyzer/BxWrongGen_Pt')
	hRight = fileHandler.getHistogram('hoMuonAnalyzer/BxRightGen_Pt')
	nRight = hRight.Integral()
	nWrong = hWrong.Integral()
	
	nEventsWithL1 = fileHandler.getHistogram('hoMuonAnalyzer/count/L1MuonPresent_Count').GetBinContent(2)
	
	etaPhiBxWrongNC = fileHandler.getGraph("hoMuonAnalyzer/graphs/BxWrongGen")
	etaPhiBxWrong = convertToHcalCoords(etaPhiBxWrongNC)
	etaPhiBxWrong.GetXaxis().SetTitle("i#eta / a.u.")
	etaPhiBxWrong.GetYaxis().SetTitle("i#phi / a.u.")
	etaPhiBxWrong.SetMarkerStyle(6)
	etaPhiBxWrong.SetMarkerColor(colorRwthDarkBlue)
	etaPhiBxWrong.SetTitle("#eta #phi plot of events with BX ID wrong")
	etaPhiBxWrong.Draw("AP")
	
	pText = TPaveText(0.7,0.85,0.9,0.9,'NDC')
	pText.AddText('Events with L1 objects: %d' % (nEventsWithL1))
	pText.AddText('Events in Plot: %d' % (etaPhiBxWrong.GetN()))
	pText.SetBorderSize(1)
	pText.Draw()
	
	chimney1Converted = convertToHcalCoords(chimney1)
	chimney2Converted = convertToHcalCoords(chimney2)
	chimney1Converted.SetLineColor(colorRwthMagenta)
	chimney2Converted.SetLineColor(colorRwthMagenta)
	chimney1Converted.Draw("same,l")
	chimney2Converted.Draw("same,l")
	
	labelCmsPrivateSimulation = getLabelCmsPrivateSimulation()
	labelCmsPrivateSimulation.Draw()
	
	legend = TLegend(0.1,0.87,0.3,0.9)
	legend.AddEntry(chimney2Converted,"chimney","l")
	legend.Draw()
	
	canvasEtaPhiBxWrong.Update()
	canvasEtaPhiBxWrong.SaveAs("plots/timing/bxWrongEtaPhi.pdf")
	canvasEtaPhiBxWrong.SaveAs("plots/timing/bxWrongEtaPhi.png")

def plotEtaOfWrongBxId():
	#Make eta histogram of the graph before
	etaPhiBxWrongNC = fileHandler.getGraph("hoMuonAnalyzer/graphs/BxWrongGen")
	etaPhiBxWrong = convertToHcalCoords(etaPhiBxWrongNC)
	canvasEtaBxWrong = TCanvas("canvasEtaBxWrong","canvasEtaBxWrong",1200,1200)
	histEtaBxWrong = TH1D("histEtaBxWrong","histEtaBxWrong;#eta Gen;# Events",20,-0.8,0.8)
	x = Double(0)
	y = Double(0)
	for i in range(0,etaPhiBxWrong.GetN()):
		etaPhiBxWrongNC.GetPoint(i,x,y)
		histEtaBxWrong.Fill(x)
	histEtaBxWrong.Draw()
	canvasEtaBxWrong.Update()
	canvasEtaBxWrong.SaveAs("plots/timing/bxWrongEta.pdf")
	canvasEtaBxWrong.SaveAs("plots/timing/bxWrongEta.png")
	
	#
	# Create a binwise normalized histogram of eta
	#
	canvasEtaBxTotal = TCanvas("canvasEtaBxTotal","canvasEtaBxtotal",1200,1200)
	etaPhiTotalNC = fileHandler.getGraph("hoMuonAnalyzer/graphs/L1ToGen")
	histEtaBxTotal= TH1D("histEtaBxTotal","histEtaBxTotal;#eta Gen;entries / 0.08 #eta",20,-0.8,0.8)
	x = Double(0)
	y = Double(0)
	for i in range(0,etaPhiTotalNC.GetN()):
		etaPhiTotalNC.GetPoint(i,x,y)
		histEtaBxTotal.Fill(x)
	histEtaBxTotal.Draw()
	canvasEtaBxTotal.SaveAs("plots/timing/bxEtaTotal.png")
	
	canvasEtaBxWrongNorm = TCanvas("canvasEtaBxWrongNorm","canvasEtaBxWrongNorm",1200,1200)
	histEtaBxWrongNorm = TH1D("histEtaBxWrongNorm","Fraction of L1 with BX ID Wrong;#eta Gen;fraction / 0.08 #eta (%)",20,-0.8,0.8)
	histEtaBxWrongNorm.SetStats(0)
	histEtaBxWrongNorm.SetLineColor(colorRwthDarkBlue)
	#fill the histogram bins
	for i in range(1,21):
		w = histEtaBxWrong.GetBinContent(i)
		t = histEtaBxTotal.GetBinContent(i)
		error = sqrt(w/float(t**2) + w**2/float(t**3))*100
		print '%5d / %d = %.2f +/- %.2f' % (w,t,w/t*100,error)
		histEtaBxWrongNorm.SetBinContent(i,w/t*100)
		histEtaBxWrongNorm.SetBinError(i,error)
	histEtaBxWrongNorm.Draw("ehist")
	canvasEtaBxWrongNorm.Update()
	canvasEtaBxWrongNorm.SaveAs("plots/timing/bxWrongEtaNorm.pdf")
	canvasEtaBxWrongNorm.SaveAs("plots/timing/bxWrongEtaNorm.png")

def plotDetectorContributionsToTiming():
	#Prepare canvas
	canvas = TCanvas("canvasDetectorContributions","detectorContributions",1200,600)
	canvas.Divide(2,1)
	canvas.cd(1).SetGridx(0)
	canvas.cd(1).SetGridy(0)
	canvas.cd(1).SetLogy()

	#prepare histogram
	hist = fileHandler.getHistogram("hoMuonAnalyzer/detectorIndexBxWrong_Multiplicity")
	hist.GetXaxis().SetRangeUser(0.5,5.5)
	hist.SetLineColor(colorRwthDarkBlue)
	hist.SetFillColor(colorRwthDarkBlue)
	#Prepare the bin labels
	x = ['RPC','DT','DT/RPC','CSC','CSC/RPC']
	for i in range(1,6):
		hist.GetXaxis().SetBinLabel(i+1,x[i-1])
	hist.GetXaxis().SetLabelFont(62)
	hist.GetYaxis().SetLabelFont(62)
	hist.GetYaxis().SetTitleFont(62)
	hist.Scale(1/hist.Integral())
	hist.GetYaxis().SetTitle('rel. fraction')
	hist.SetStats(0)
	hist.SetTitle('Subdetectors in wrong L1 BX ID')
	hist.Draw()
	
	#add label
	label = getLabelCmsPrivateSimulation()
	label.Draw()
	
	#Make a histogram to grey out the CSC parts
	histGrey = TH1D('greyHist','grey hist',2,3.5,5.5)
	histGrey.SetBinContent(1,1)
	histGrey.SetBinContent(2,1)
	histGrey.SetFillStyle(3004)
	histGrey.SetLineWidth(2)
	histGrey.SetLineColor(1)
	histGrey.SetFillColor(1)
	histGrey.Draw('same')
	
	pText = TPaveText(0.7,0.33,0.79,0.62,"NDC")
	tTextPointer = pText.AddText('Not in #eta range')
	tTextPointer.SetTextAngle(90)
	tTextPointer.SetTextSize(0.039)
	pText.SetBorderSize(1)
	pText.SetFillColor(0)
	pText.Draw()
	
	canvas.cd(2).SetGridx(0)
	canvas.cd(2).SetGridy(0)
	canvas.cd(2).SetLogy()

	#prepare second histogram
	hist2 = fileHandler.getHistogram("hoMuonAnalyzer/detectorIndexBxRight_Multiplicity")
	hist2.GetXaxis().SetRangeUser(0.5,5.5)
	hist2.SetLineColor(colorRwthDarkBlue)
	hist2.SetFillColor(colorRwthDarkBlue)
	#Prepare the bin labels
	x = ['RPC','DT','DT/RPC','CSC','CSC/RPC']
	for i in range(1,6):
		hist2.GetXaxis().SetBinLabel(i+1,x[i-1])
	hist2.GetXaxis().SetLabelFont(62)
	hist2.GetYaxis().SetLabelFont(62)
	hist2.GetYaxis().SetTitleFont(62)
	hist2.Scale(1/hist2.Integral())
	hist2.GetYaxis().SetTitle('rel. fraction')
	hist2.SetStats(0)
	hist2.SetTitle('Subdetectors in right L1 BX ID')
	hist2.Draw()
	
	#add label
	label2 = getLabelCmsPrivateSimulation()
	label2.Draw()
	
	histGrey.DrawCopy('same')
	pText2 = pText.Clone("pText2")
	pText2.Draw()
	
	canvas.Update()
	canvas.SaveAs('plots/timing/bxWrongDetectorContributions.pdf')
	canvas.SaveAs('plots/timing/bxWrongDetectorContributions.png')
	
	
	
	return canvas,hist,label,histGrey,pText,pText2,label2

def plotPtAndEtaOfWrongBxId():
	#Prepare canvas
	canvas = TCanvas("canvasEtaPtBxWrong","EtaPtBxWrong",1200,1200)
	canvas.cd().Draw()
	#prepare histogram
	hist = fileHandler.getHistogram("hoMuonAnalyzer/etaPhi/3D/BxWrongGen_EtaPhiPt")

	stack = THStack(hist,"zx","2dStack","",-1,-1,-1,-1,"zx","")

	#Create new histogram and add the histograms from the stack
	histNew = TH2D("histPtEtaBxWrong","p_{T} vs. #eta distribution for wrong BX ID;#eta;p_{T} / 5 GeV;#",40,-1.6,1.6,40,0,200)
	histNew.GetYaxis().SetTitleOffset(1.2)
	for i in stack.GetHists():
		histNew.Add(i)

	gStyle.SetPalette(1)
	histNew.SetStats(0)
	histNew.Draw('COLZ')
	canvas.Update()

	palette = histNew.FindObject("palette")
	palette.SetX1NDC(0.9)
	palette.SetX2NDC(0.92)
	#add label
	label = getLabelCmsPrivateSimulation()
	label.Draw()
	
	canvas.Update()
	canvas.SaveAs('plots/timing/bxWrongEtaPt.pdf')
	canvas.SaveAs('plots/timing/bxWrongEtaPt.png')
	return canvas,hist,stack,histNew,label

def plotPtAndPhiOfWrongBxId():
	#Prepare canvas
	canvas = TCanvas("canvasPtPhiBxWrong","PtPhiBxWrong",1200,1200)
	canvas.cd().Draw()
	#prepare histogram
	hist = fileHandler.getHistogram("hoMuonAnalyzer/etaPhi/3D/BxWrongGen_EtaPhiPt")

	stack = THStack(hist,"zy","2dStack","",-1,-1,-1,-1,"zy","")

	#Create new histogram and add the histograms from the stack
	histNew = TH2D("histPtPhiBxWrong","p_{T} vs. #phi distribution for wrong BX ID;#phi;p_{T} / 5 GeV;#",80,-3.2,3.2,40,0,200)
	histNew.GetYaxis().SetTitleOffset(1.2)
	for i in stack.GetHists():
		histNew.Add(i)

	gStyle.SetPalette(1)
	histNew.SetStats(0)
	setupAxes(histNew)
	histNew.Draw('COLZ')
	canvas.Update()

	palette = histNew.FindObject("palette")
	palette.SetX1NDC(0.9)
	palette.SetX2NDC(0.92)
	#add label
	label = getLabelCmsPrivateSimulation()
	label.Draw()
	
	canvas.Update()
	canvas.SaveAs('plots/timing/bxWrongPtPhi.pdf')
	canvas.SaveAs('plots/timing/bxWrongPtPhi.png')
	return canvas,hist,stack,histNew,label


def plotImprovementInDt():
	#Prepare canvas
	setPlotStyle()
	canvas = TCanvas("canvasDtImprovement","DT improvement",1200,1200)
	canvas.SetLogy()
	histDt = fileHandler.getHistogram("hoMuonAnalyzer/BxDtOnly_BxId")
	
	#Define variables for integrals
	histHoTime = fileHandler.getHistogram('hoMuonAnalyzer/hoRecHitsAboveThr_Time')
	integralHoCorrect = histHoTime.Integral(histHoTime.FindBin(-12.5),histHoTime.FindBin(12.5))
	integralHoTotal = histHoTime.Integral()
	integralHoOutside = integralHoTotal - integralHoCorrect
	hoFractionWrong = integralHoOutside/float(integralHoTotal)
	hoFractionRight = integralHoCorrect/float(integralHoTotal)
	
	#Print some information
	heading = 'Integrals of the Ho timing:'
	print 80*'#'
	print heading
	print len(heading)*'-'
	print 'Timing correct:\t%d\t=>\t%5.2f%% +/- %5.2f%%'%(integralHoCorrect,hoFractionRight*100,100*calcSigma(integralHoCorrect, integralHoTotal))
	print 'Timing outside:\t%d\t=>\t%5.2f%% +/- %f%%'%(integralHoOutside,hoFractionWrong*100,100*calcSigma(integralHoOutside, integralHoTotal))
	print 'Timing total:%d'%(integralHoTotal)
	print
	
	#Define Variables for bx id counts
	dtBx0 = histDt.GetBinContent(6)
	dtBxM1 = histDt.GetBinContent(5)
	dtBxP1 = histDt.GetBinContent(7)
	dtBxTotal = dtBx0 + dtBxM1 + dtBxP1
	dtFractionWrongM1 = dtBxM1/float(dtBxTotal)
	dtFractionWrongP1 = dtBxP1/float(dtBxTotal)

	#Print some information
	heading = 'Bin contents for DT timing:'
	print heading
	print len(heading)*'-'
	print 'BX ID 0:\t%d\t=>\t%5.2f%% +/- %5.2f%%'%(dtBx0,dtBx0/float(dtBxTotal)*100,calcSigma(dtBx0, dtBxTotal))
	print 'BX ID -1:\t%d\t=>\t%5.2f%% +/- %5.2f%%'%(dtBxM1,dtFractionWrongM1*100,calcSigma(dtBxM1, dtBxTotal))
	print 'BX ID +1:\t%d\t=>\t%5.2f%% +/- %5.2f%%'%(dtBxP1,dtFractionWrongP1*100,calcSigma(dtBxP1, dtBxTotal))
	print 'BX ID total:\t%d' % (dtBxTotal)
	print
	
	#Calculate corrected numbers
	correctedBxIdM1 = dtBxM1 + hoFractionWrong*dtBx0/2. - hoFractionRight*dtBxM1
	correctedBxId0 = dtBx0 - hoFractionWrong*dtBx0 + hoFractionRight*dtBxM1 + hoFractionRight*dtBxP1
	correctedBxIdP1 = dtBxP1 + hoFractionWrong*dtBx0/2. - hoFractionRight*dtBxP1
	correctedTotal = correctedBxIdM1 + correctedBxId0 + correctedBxIdP1
	correctedRightFraction = correctedBxId0/float(correctedTotal)
	
	heading = 'DT After correction:'
	print heading
	print len(heading)*'-'
	print 'BX -1:\t',int(correctedBxIdM1)
	print 'BX 0:\t',int(correctedBxId0)
	print 'BX +1:\t',int(correctedBxIdP1)
	print 
	#Fill corrected histogram
	histNew = TH1D("histNew","BX ID in DT only triggers;BX ID;rel. fraction",6,-2.5,3.5)
	histNew.SetBinContent(histNew.FindBin(-1),correctedBxIdM1)
	histNew.SetBinContent(histNew.FindBin(0),correctedBxId0)
	histNew.SetBinContent(histNew.FindBin(1),correctedBxIdP1)
	histNew.SetLineColor(colorRwthMagenta)
	histNew.SetStats(0)
	histNew.Scale(1/histNew.Integral())
	histNew.SetLineStyle(9)
	setupAxes(histNew)
	
	histDt.GetXaxis().SetRangeUser(-3,3)
	histDt.SetLineWidth(3)
	histDt.Scale(1/histDt.Integral())
	histDt.SetLineColor(colorRwthDarkBlue)
	
	
	histNew.Draw()
	histDt.Draw('same')
	histNew.Draw('same')
	

	#Add label
	label = getLabelCmsPrivateSimulation()
	label.Draw()
	
	#Add legend
	legend = TLegend(0.7,0.65,0.9,0.8)
	legend.AddEntry(histDt,"DT Only","l")
	legend.AddEntry(histNew,"DT shifted with HO","l")
	legend.SetBorderSize(1)
	legend.Draw()
	
	#Add text object
	pText = TPaveText(0.52,0.8,0.9,0.9,'NDC')
	pText.AddText('Fraction in BX ID 0: %6.3f%% #pm %6.3f%%' % (dtBx0/float(dtBxTotal)*100,calcSigma(dtBx0, dtBxTotal)))
	pText.AddText('Fraction in BX ID 0 (HO corr.): %6.3f%% #pm %6.3f%%' % (correctedRightFraction*100,calcSigma(correctedBxId0, correctedTotal)))
	pText.SetBorderSize(1)
	pText.SetFillColor(0)
	pText.Draw()
	
	pText2 = TPaveText(0.7,0.6,0.9,0.65,'NDC')
	pText2.AddText('Entries: %d' % (histDt.GetEntries()))
	pText2.SetBorderSize(1)
	pText2.SetFillColor(0)
	pText2.Draw()
	
	#Print again some information
	heading = 'Fraction of correct BXID:'
	print heading
	print len(heading)*'-'
	print 'Uncorrected:\t%5.2f%% #pm %f%%' % (dtBx0/float(dtBxTotal)*100,calcSigma(dtBx0, dtBxTotal))
	print 'Corrected\t%5.2f%% #pm %f%%' % (correctedRightFraction*100,calcSigma(correctedBxId0, correctedTotal))
	print 80*'#'
	
	canvas.Update()
	canvas.SaveAs('plots/timing/correctedDt.pdf')
	canvas.SaveAs('plots/timing/correctedDt.png')
	canvas.SaveAs('plots/timing/correctedDt.root')
	return canvas, histDt,histNew,label,legend,pText2,pText
	
plotDeltaTime()
plotEtaOfWrongBxId()
plotEtaPhiOfWrongBxId()
plotFractionsOfBxId()
plotL1BxId()
res = plotHoTime()
res2 = plotDetectorContributionsToTiming()
res3 = plotPtAndEtaOfWrongBxId()
res4 = plotImprovementInDt()
res5 = plotPtAndPhiOfWrongBxId()
res6 = plotHoTimeLog()
raw_input('-->')