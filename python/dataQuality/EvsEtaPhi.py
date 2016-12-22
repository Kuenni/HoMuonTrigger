#!/usr/bin/python
from ROOT import TCanvas,ROOT,TFile,TF1,TLine,gROOT,TPaveText,TH1D,Double,TH2D,THStack,gStyle,TMarker
from plotting.PlotStyle import getLabelCmsPrivateSimulation,setupPalette,\
	getTH2D
from plotting.PlotStyle import setupAxes,drawHoBoxes
from plotting.PlotStyle import setStatBoxOptions,setStatBoxPosition,pyplotCmsPrivateLabel
from plotting.Utils import setupEAvplot, L1_PHI_BIN, L1_ETA_BIN, calcPercent, calcSigma,\
	fillGraphIn2DHist

from plotting.Plot import Plot

import numpy as np
import matplotlib.pyplot as plt
import math
from OfflineAnalysis.HoMatcher import HO_BIN

class EvsEtaPhi(Plot):
	
	def __init__(self,filename,data = False,debug = False):
		Plot.__init__(self,filename,data,debug)
		self.createPlotSubdir('averageEnergy')
		self.fileHandler.printNEvents()
	'''
	Plots the average energy seen in in the tiles around the direction
	of the L1 muons
	'''
	def plotAverageEnergyAroundL1(self):
		canvas = TCanvas('canvasAverageEnergy','Average energy',900,800)
		canvas.cd().SetLogz()
		canvas.cd().SetRightMargin(.15)

		hSum = self.fileHandler.getHistogram('averageEnergy/averageEnergyAroundPoint' + self.key + '_SummedEnergy')
		hCounter = self.fileHandler.getHistogram('averageEnergy/averageEnergyAroundPoint' + self.key + '_Counter')
	
		for i in range(0,hSum.GetNbinsX()):
			for j in range(0,hSum.GetNbinsY()):
				if hCounter.GetBinContent(hCounter.GetBin(i,j)) != 0:
					hSum.SetBinContent(hSum.GetBin(i,j),hSum.GetBinContent(hSum.GetBin(i,j))/hCounter.GetBinContent(hCounter.GetBin(i,j)))
					pass
		hSum.GetXaxis().SetRangeUser(-0.6,0.6)
		hSum.GetYaxis().SetRangeUser(-0.6,0.6)
		hSum.GetZaxis().SetTitle('Reconstructed Energy / GeV')
		hSum.SetTitle(';#Delta#eta;#Delta#phi / rad;Reconstructed Energy / GeV')#'Average Energy in HO tiles around L1 direction
		hSum.Draw('colz')
						
		canvas.Update()
		
		#Setup plot style
		setupAxes(hSum)	
		setStatBoxOptions(hSum,1100)
		setStatBoxPosition(hSum,x1=.65,x2=.85)
		setupPalette(hSum,shiftBy=.05)
	
		canvas.Update()
		self.storeCanvas(canvas,'averageEnergy',marginRight=.15)
		return canvas,hSum,hCounter,drawHoBoxes(canvas)
	
	
	def calculateCentralFractionInTight(self):
		hSum = self.fileHandler.getHistogram('deltaEtaDeltaPhiEnergy/averageEMaxAroundPointpatTightToL1Muons_2dSummedWeights')
		hCounter = self.fileHandler.getHistogram('deltaEtaDeltaPhiEnergy/averageEMaxAroundPointpatTightToL1Muons_2dCounter')

		points = []
		sum = 0
		for x in range(0,5):
			for y in range(0,5):
				sum += hCounter.GetBinContent(hSum.FindBin(x*0.0435 - 0.087,y*0.0435 - 0.087))
				p = TMarker(x*0.0435 - 0.087,y*0.0435 - 0.087,20)
		#		p.Draw('same')
				points.append(p)
		print sum
		print hCounter.Integral()
		print '-----'
		print sum/hCounter.Integral()
	
	def plotAverageEMaxAroundL1(self):
		canvas = TCanvas('canvasAverageEMax','Average EMax',1200,1200)
		canvas.cd().SetLogz()
		
		hSum = self.fileHandler.getHistogram('deltaEtaDeltaPhiEnergy/averageEMaxAroundPoint' + self.key + '_2dSummedWeights')
		hCounter = self.fileHandler.getHistogram('deltaEtaDeltaPhiEnergy/averageEMaxAroundPoint' + self.key + '_2dCounter')
		
		hSum = setupEAvplot(hSum, hCounter,same=True,limitForAll=0.3)
		hSum.SetTitle('Mean E_{Max} in HO tiles around L1 direction')
		hSum.SetMaximum(2)
		hSum.Draw('colz')
		setupEAvplot(hCounter,same=True,limitForAll=0.3).Draw('same,text')
		
		canvas.Update()		
				
		setupPalette(hSum)
		
		canvas.Update()
		self.storeCanvas(canvas, 'averageEmax')
		hCounter.SaveAs('histogramEMaxCounter.root')
		
		return canvas,hSum,hCounter
	
	def plot1DEnergyAroundL1(self):	
		'''
			eta[P,M][2,1,0]phi[P,M][2,1,0]_averageEnergyAroundPoint
			Central tile is central
		'''
		histList = []
		fitList = []
		labelList = []
		canvas = TCanvas('canvas1DEnergy','1D energy',1200,1200)
		for p in reversed(range(-2,3)):
			for e in range(-2,3):
				if e == 0 and p == 0:
					histList.append(self.fileHandler.getHistogram('etaPhi/energy1D/central_averageEnergyAroundPoint' + self.key))
				else:
					histName = ('etaPhi/energy1D/eta%s%dPhi%s%d_averageEnergyAroundPoint' + self.key) % ('P' if e >= 0 else 'M',abs(e),'P' if p >= 0 else 'M',abs(p))
					histList.append(self.fileHandler.getHistogram(histName))
		canvas.Divide(5,5)
		for i,hist in enumerate(histList):
			canvas.cd(i+1).SetLogy()
			hist.GetXaxis().SetRangeUser(-0.5,4)
			hist.SetLineWidth(3)
			setupAxes(hist)
			hist.Draw()
			fit = TF1('fit%d' % (i),'landau',0.5,2)
			hist.Fit(fit,'RQ')
			label = TPaveText(0.6,0.7,0.9,0.9,"NDC")
			label.AddText('MPV: %5.2f' % (fit.GetParameter(1)))
			label.Draw()
			labelList.append(label)
			fitList.append(fit)
		canvas.Update()
		self.storeCanvas(canvas, '1DPlots')
		return histList,canvas,fitList,labelList
	
	def plot1DEMaxAroundL1(self):	
		'''
			eta[P,M][2,1,0]phi[P,M][2,1,0]_averageEnergyAroundPoint
			Central tile is central
		'''
		histList = []
		fitList = []
		labelList = []
		canvas = TCanvas('canvas1DEMax','1D EMax',1200,1200)
		for p in reversed(range(-2,3)):
			for e in range(-2,3):
				if e == 0 and p == 0:
					histList.append(self.fileHandler.getHistogram('etaPhi/energy1D/central_averageEMaxAroundPoint' + self.key))
				else:
					histName = ('etaPhi/energy1D/eta%s%dPhi%s%d_averageEMaxAroundPoint' + self.key) % ('P' if e >= 0 else 'M',abs(e),'P' if p >= 0 else 'M',abs(p))
					histList.append(self.fileHandler.getHistogram(histName))
		canvas.Divide(5,5)
		for i,hist in enumerate(histList):
			canvas.cd(i+1).SetLogy()
			hist.GetXaxis().SetRangeUser(-0.5,4)
			hist.SetLineWidth(3)
			setupAxes(hist)
			hist.Draw()
			fit = TF1('fit%d' % (i),'landau',0.5,2)
			hist.Fit(fit,'RQ')
			label = TPaveText(0.6,0.7,0.9,0.9,"NDC")
			label.AddText('MPV: %5.2f' % (fit.GetParameter(1)))
			label.Draw()
			labelList.append(label)
			fitList.append(fit)
		canvas.Update()
		self.storeCanvas(canvas, '1DEMaxPlots')
		return histList,canvas,fitList,labelList
	
	def plotMPVs(self,fitList):
		x = np.arange(-2,3)
		y = np.arange(-2,3)
		z = []
		for fit in fitList:
			z.append(fit.GetParameter(1))
		
		z = np.reshape(z, [len(x), len(y)])
		z = np.flipud(z)
		
		fig, ax = plt.subplots()
		im = ax.pcolor(z)
		ax.set_xticks(np.arange(z.shape[0])+0.5, minor=False)
		ax.set_yticks(np.arange(z.shape[1])+0.5, minor=False)
		ax.set_xticklabels(x, minor=False)
		ax.set_yticklabels(y, minor=False)
		
		colorbar = fig.colorbar(im)
		colorbar.set_label('MPV reconstructed energy / GeV')
		ax.axis('tight')
		pyplotCmsPrivateLabel(ax,y=1)
		plt.xlabel(r'$\Delta i\eta$')
		plt.ylabel(r'$\Delta i\phi$')
		plt.savefig('plots/averageEnergy/mpv.png')
		plt.show()
	
	def compareHistogramMethods(self):
		canvas = TCanvas('cComparison','Comparison btween histograms')
		
	#	canvas.Divide(2,1)
		
		histNormal = self.fileHandler.getHistogram('averageEnergy/averageEnergyAroundPoint' + self.key + '_SummedEnergy')
		histNormalCounter = self.fileHandler.getHistogram('averageEnergy/averageEnergyAroundPoint' + self.key + '_Counter')
		
		histNormal = setupEAvplot(histNormal, histNormalCounter,same=True,limitForAll=0.6)
		
	#	histNew = self.fileHandler.getHistogram('deltaEtaDeltaPhiEnergy/averageEnergyAroundPoint_2dSummedWeightsIEtaIPhi')
	#	histNewCounter = self.fileHandler.getHistogram('deltaEtaDeltaPhiEnergy/averageEnergyAroundPoint_2dCounterIEtaIPhi')
		
		canvas.cd(1).SetLogz()
		
		histNormal.SetTitle('Average Energy in HO tiles around L1 direction, i#eta by binning')
		histNormal.SetStats(1)
		histNormal.Draw('colz')
		
		label = getLabelCmsPrivateSimulation()
		label.Draw()
		
	#	canvas.cd(2).SetLogz()
		
	#	histNew = average2DHistogramBinwise(histNew, histNewCounter)
	#	histNew.GetXaxis().SetRangeUser(-8,8)
	#	histNew.GetYaxis().SetRangeUser(-8,8)
	#	histNew.GetXaxis().SetTitle('#Delta#eta')
	#	histNew.GetYaxis().SetTitle('#Delta#phi')
	#	histNew.GetZaxis().SetTitle('Reconstructed Energy / GeV')
	#	histNew.SetTitle('Mean Energy in HO tiles around L1 direction, i#eta by binning')
	#	histNew.Draw('colz')
			
	#	label2 = getLabelCmsPrivateSimulation()
	#	label2.Draw()
		
		canvas.Update()
		
		#Setup plot style
		setStatBoxOptions(histNormal,1100)
		setStatBoxPosition(histNormal)
		setupPalette(histNormal)
		
	#	setupAxes(histNew)	
	#	setStatBoxOptions(histNew,1100)
	#	setStatBoxPosition(histNew)
	#	setupPalette(histNew)
	
		canvas.Update()
		
		#TODO: Print the bin contents subtracted
		
		return canvas, histNormal,label#,histNew,label2
	
	def plotEavForTightMuons(self):
		canvas = TCanvas('canvasEavTightMuons','EAv Tight muons',800,800)
		canvas.cd().SetLogz()
			
		hSum = self.fileHandler.getHistogram('averageEnergy/averageEnergyAroundPointpatTightToL1Muons_SummedEnergy')
		hCounter = self.fileHandler.getHistogram('averageEnergy/averageEnergyAroundPointpatTightToL1Muons_Counter')
	
		hSum = setupEAvplot(hSum, hCounter,same=True,limitForAll=0.3)
		hSum.SetTitle('Average E_{Rec} in HO tiles around tight L1 direction')
		hSum.SetMaximum(1.2)
		hSum.SetMinimum(5e-3)
		hSum.Draw('colz')
		label = self.drawLabel()
		canvas.Update()		
		setupPalette(hSum)
		canvas.Update()
		#boxes = drawHoBoxes(canvas)
		self.storeCanvas(canvas,'eAverageTightMuons')
		return canvas,hSum,label#,boxes
	
	### ##########################################################
	### Put data of fraction of events inside tile grid on the CLI
	### ##########################################################
	def outputFractionsInTileGrid(self,hist):
		#Calculate fraction in grids
		integralCentral = hist.Integral(hist.GetXaxis().FindBin(-.0435),hist.GetXaxis().FindBin(.0435),
									hist.GetYaxis().FindBin(-.0435),hist.GetYaxis().FindBin(.0435))
		integral3x3 = hist.Integral(hist.GetXaxis().FindBin(-.1305),hist.GetXaxis().FindBin(.1305),
									hist.GetYaxis().FindBin(-.1305),hist.GetYaxis().FindBin(.1305))
		integral5x5 = hist.Integral(hist.GetXaxis().FindBin(-.2175),hist.GetXaxis().FindBin(.2175),
									hist.GetYaxis().FindBin(-.2175),hist.GetYaxis().FindBin(.2175))
		integralTotal = hist.Integral()
		
		self.debug(80*'#')
		self.debug('%20s:%5.2f%% +/- %5.2f%%' % ('Central Fraction',calcPercent(integralCentral,integralTotal),
												calcSigma(integralCentral,integralTotal)*100))
		self.debug('%20s:%5.2f%% +/- %5.2f%%' % ('3x3 Fraction',calcPercent(integral3x3,integralTotal),
												calcSigma(integral3x3,integralTotal)*100))
		self.debug('%20s:%5.2f%% +/- %5.2f%%' % ('5x5 Fraction',calcPercent(integral5x5,integralTotal),
												calcSigma(integral5x5,integralTotal)*100))
		self.debug(80*'#')
	
	#########################################################
	# Generalized plotting script: Emax around L1 direction #
	#########################################################
	def makeEmaxPlot(self,source, title = ""):
		if title == "":
			title = '# of E_{Max} in HO tiles around L1 direction ' + source
			
		canvas = TCanvas('canvasEmaxcounts' + source,'E max counts' + source,900,800)
		canvas.cd().SetLogz()
		canvas.cd().SetRightMargin(.15)
		hCounter = self.fileHandler.getHistogram('deltaEtaDeltaPhiEnergy/averageEMaxAroundPoint' + source + '_2dCounter')
		hCounter.GetXaxis().SetRangeUser(-.6,.6)
		hCounter.GetYaxis().SetRangeUser(-.6,.6)
		hCounter.SetTitle(';#Delta#eta;#Delta#phi / rad;Entries')#title +
		hCounter.SetStats(0)
		hCounter.Draw('colz')
		canvas.Update()		
		setupAxes(hCounter)
		setupPalette(hCounter,x2ndc=.87)
		canvas.Update()
		#boxes = drawHoBoxes(canvas)
		self.storeCanvas(canvas,'eMaxCounts' + source,marginRight=.15)

		#Output fractions in grid		
		self.debug('Emax fraction for ' + source)
		self.outputFractionsInTileGrid(hCounter)
		return canvas,hCounter
	
	####################
	# Predefined plots #
	####################
	def plotEMaxCountsForTightMuons(self):
		return self.makeEmaxPlot('patTightToL1Muons')
	
	def plotEMaxCounts(self):
		return self.makeEmaxPlot('L1MuonPresent')
		
	def plotEavPerWheelForTightMuons(self):
		canvas = TCanvas('canvasEavPerWheelTightMuons','EAv Per Wheel Tight muons',1800,800)
		canvas.Divide(3,1)
				
		hM1Energy = self.fileHandler.getHistogram('averageEnergy/wh1m/averageEnergyAroundPointpatTightToL1Muons_wh-1SummedEnergy')
		hM1Counter = self.fileHandler.getHistogram('averageEnergy/wh1m/averageEnergyAroundPointpatTightToL1Muons_wh-1Counter')
		hM1Energy = setupEAvplot(hM1Energy, hM1Counter)
		hM1Energy.SetStats(0)
	
		h0Energy = self.fileHandler.getHistogram('averageEnergy/wh0/averageEnergyAroundPointpatTightToL1Muons_wh0SummedEnergy')
		h0Counter = self.fileHandler.getHistogram('averageEnergy/wh0/averageEnergyAroundPointpatTightToL1Muons_wh0Counter')
		h0Energy = setupEAvplot(h0Energy, h0Counter)
		h0Energy.SetStats(0)
	
		hP1Energy = self.fileHandler.getHistogram('averageEnergy/wh1p/averageEnergyAroundPointpatTightToL1Muons_wh1SummedEnergy')
		hP1Counter = self.fileHandler.getHistogram('averageEnergy/wh1p/averageEnergyAroundPointpatTightToL1Muons_wh1Counter')
		hP1Energy = setupEAvplot(hP1Energy, hP1Counter)
		hP1Energy.SetStats(0)
	
		canvas.cd(1).SetLogz()
		setupAxes(hM1Energy)
		hM1Energy.SetMaximum(1.2)
		hM1Energy.SetMinimum(5e-3)
		hM1Energy.Draw('colz')
		canvas.Update()
		setupPalette(hM1Energy)
		label1 = self.drawLabel()
		
		canvas.cd(2).SetLogz()
		setupAxes(h0Energy)
		h0Energy.SetMaximum(1.2)
		h0Energy.SetMinimum(5e-3)
		h0Energy.Draw('colz')
		canvas.Update()
		setupPalette(h0Energy)
		label2 = self.drawLabel()
		
		canvas.cd(3).SetLogz()
		setupAxes(hP1Energy)
		hP1Energy.SetMaximum(1.2)
		hP1Energy.SetMinimum(5e-3)
		hP1Energy.Draw('colz')
		canvas.Update()
		setupPalette(hP1Energy)
		label3 = self.drawLabel()
	
		canvas.Update()
		
		self.storeCanvas(canvas,'eAveragePerWheelTightMuons')
		
		return hM1Energy,canvas,h0Energy,hP1Energy,h0Counter,label1,label2,label3

	def plotEAveragePerWheel(self):
		canvas = TCanvas('cEAvPerWheel',"E Average per Wheel",1800,800)
		canvas.Divide(3,1)
	
		hM1Energy = self.fileHandler.getHistogram('averageEnergy/wh1m/averageEnergyAroundPoint' + self.key + '_wh-1SummedEnergy')
		hM1Counter = self.fileHandler.getHistogram('averageEnergy/wh1m/averageEnergyAroundPoint' + self.key + '_wh-1Counter')
		hM1Energy = setupEAvplot(hM1Energy, hM1Counter)
		hM1Energy.SetStats(0)
	
		h0Energy = self.fileHandler.getHistogram('averageEnergy/wh0/averageEnergyAroundPoint' + self.key + '_wh0SummedEnergy')
		h0Counter = self.fileHandler.getHistogram('averageEnergy/wh0/averageEnergyAroundPoint' + self.key + '_wh0Counter')
		h0Energy = setupEAvplot(h0Energy, h0Counter)
		h0Energy.SetStats(0)
	
		hP1Energy = self.fileHandler.getHistogram('averageEnergy/wh1p/averageEnergyAroundPoint' + self.key + '_wh1SummedEnergy')
		hP1Counter = self.fileHandler.getHistogram('averageEnergy/wh1p/averageEnergyAroundPoint' + self.key + '_wh1Counter')
		hP1Energy = setupEAvplot(hP1Energy, hP1Counter)
		hP1Energy.SetStats(0)
	
		canvas.cd(1).SetLogz()
		setupAxes(hM1Energy)
		hM1Energy.SetMaximum(1.2)
		hM1Energy.SetMinimum(5e-3)
		hM1Energy.Draw('colz')
		canvas.Update()
		setupPalette(hM1Energy)
		label1 = self.drawLabel()

		canvas.cd(2).SetLogz()
		setupAxes(h0Energy)
		h0Energy.SetMaximum(1.2)
		h0Energy.SetMinimum(5e-3)
		h0Energy.Draw('colz')
		#h0Counter.Draw('same,text')
		canvas.Update()
		setupPalette(h0Energy)
		label2 = self.drawLabel()
		
		canvas.cd(3).SetLogz()
		setupAxes(hP1Energy)
		hP1Energy.SetMaximum(1.2)
		hP1Energy.SetMinimum(5e-3)
		hP1Energy.Draw('colz')
		canvas.Update()
		setupPalette(hP1Energy)
		label3 = self.drawLabel()
	
		canvas.Update()
		self.storeCanvas(canvas,'eAveragePerWheel')
		
		return hM1Energy,canvas,h0Energy,hP1Energy,h0Counter,label1,label2,label3
	
	def plotEtaPhiForTightL1(self):
		canvas = TCanvas("cEtaPhi","Eta Phi",1200,900)
		graphAll = self.fileHandler.getGraph('graphs/patTightToL1Muons')
		graphWithHo = self.fileHandler.getGraph('graphs/patTightToL1Muons3x3')
				
		halfPhiBinwidth = L1_PHI_BIN/2.
		l1BinOffset = L1_PHI_BIN*3/4.
		
		histAll = TH2D('hEtaPhiAll',"",30,-15*L1_ETA_BIN	,15*L1_ETA_BIN,
					144, -math.pi,math.pi)
		histWithHo = TH2D('hEtaPhiWithHO',"",30,-15*L1_ETA_BIN,15*L1_ETA_BIN,
					144, -math.pi,math.pi)
		
		x = Double(0)
		y = Double(0)
		
		for i in range(0,graphAll.GetN()):
			graphAll.GetPoint(i,x,y)
			histAll.Fill(x,y)
			
		for i in range(0,graphWithHo.GetN()):
			graphWithHo.GetPoint(i,x,y)
			histWithHo.Fill(x,y)
		
		canvas.cd().SetRightMargin(.15)
		histAll.SetStats(0)
		histAll.GetXaxis().SetRangeUser(-1,1)
		histAll.SetTitle(histAll.GetTitle() + ';#eta_{L1};#phi_{L1} / rad;Entries')
		histAll.Draw('colz')
		canvas.Update()
		setupAxes(histAll)
		setupPalette(histAll,shiftBy=.05)
		#label1 = self.drawLabel(x1ndc=.55,x2ndc=.85)
		histAll.GetZaxis().SetTitleOffset(1.)
		histAll.GetZaxis().SetRangeUser(0,1250)
		canvas.Update()
		
		canvas2 = TCanvas("cEtaPhiAndHo","Eta Phi And HO",1200,900)
		
		canvas2.cd().SetRightMargin(.15)
		histWithHo.SetStats(0)
		histWithHo.GetXaxis().SetRangeUser(-1,1)
		histWithHo.SetTitle(histWithHo.GetTitle() + ';#eta_{L1};#phi_{L1} / rad;Entries')
		histWithHo.Draw('colz')
		#label2 = self.drawLabel(x1ndc=.55,x2ndc=.85)
		histWithHo.GetZaxis().SetRangeUser(0,1250)
		
		canvas2.Update()
		setupAxes(histWithHo)
		setupPalette(histWithHo,shiftBy=.05)
		histWithHo.GetZaxis().SetTitleOffset(1.)
		
		canvas2.Update()
		
		self.storeCanvas(canvas, 'etaPhiForTightL1',marginRight=.15)
		self.storeCanvas(canvas2, 'etaPhiForTightL1AndHo',marginRight=.15)
		return canvas,histAll,histWithHo,canvas2
	
	def printFractionsAboveEthr(self):
		c = TCanvas('cFractions','fractionsAbveThr')

		gStyle.SetPaintTextFormat('4.2f') 
		lateralSpreadForPlot = 5;
		nBinsPerHoPhi = 1;
		phiBinOffsetFactor = 1/2.;
		nTotalBinsPhi = 2*nBinsPerHoPhi*lateralSpreadForPlot + 1;
				
		histTotal = TH2D('histTightTotal','',nTotalBinsPhi,-lateralSpreadForPlot*HO_BIN - HO_BIN*phiBinOffsetFactor,
			lateralSpreadForPlot*HO_BIN + HO_BIN*phiBinOffsetFactor,
			nTotalBinsPhi,-lateralSpreadForPlot*HO_BIN - HO_BIN*phiBinOffsetFactor,
			lateralSpreadForPlot*HO_BIN + HO_BIN*phiBinOffsetFactor)
		
		fillGraphIn2DHist(self.fileHandler.getGraph('graphs/eAvAboveThrCounterL1MuonPresent'), histTotal)
		
		hist3x3 = TH2D('histTight5x5',';#Delta#eta;#Delta#phi / rad;Fraction / %',nTotalBinsPhi,-lateralSpreadForPlot*HO_BIN - HO_BIN*phiBinOffsetFactor,
			lateralSpreadForPlot*HO_BIN + HO_BIN*phiBinOffsetFactor,
			nTotalBinsPhi,-lateralSpreadForPlot*HO_BIN - HO_BIN*phiBinOffsetFactor,
			lateralSpreadForPlot*HO_BIN + HO_BIN*phiBinOffsetFactor)
		
		fillGraphIn2DHist(self.fileHandler.getGraph('graphs/eAvAboveThr5x5L1MuonPresent'), hist3x3)
		hist3x3.Divide(histTotal)
		hist3x3.SetStats(0)
		hist3x3.GetXaxis().SetRangeUser(-.3,.3)
		hist3x3.GetYaxis().SetRangeUser(-.3,.3)
		hist3x3.Draw('colz')
		c.Update()
		setupAxes(hist3x3)
		setupPalette(hist3x3, shiftBy=.05)
		hCopy = hist3x3.Clone('hCloneTight')
		hCopy.Scale(100)
		hCopy.SetMarkerSize(2)
		hCopy.Draw('text,same')
		self.storeCanvas(c,'fractionAboveThr',marginRight=.15)
		return hist3x3,hCopy,c
	
	def printFractionsAboveEthrTight(self):
		c = TCanvas('cFractionsTight','fractionsTightAbveThr')

		gStyle.SetPaintTextFormat('4.2f') 
		lateralSpreadForPlot = 5;
		nBinsPerHoPhi = 1;
		phiBinOffsetFactor = 1/2.;
		nTotalBinsPhi = 2*nBinsPerHoPhi*lateralSpreadForPlot + 1;
				
		histTotal = TH2D('histTightTotal','',nTotalBinsPhi,-lateralSpreadForPlot*HO_BIN - HO_BIN*phiBinOffsetFactor,
			lateralSpreadForPlot*HO_BIN + HO_BIN*phiBinOffsetFactor,
			nTotalBinsPhi,-lateralSpreadForPlot*HO_BIN - HO_BIN*phiBinOffsetFactor,
			lateralSpreadForPlot*HO_BIN + HO_BIN*phiBinOffsetFactor)
		
		fillGraphIn2DHist(self.fileHandler.getGraph('graphs/eAvAboveThrCounterpatTightToL1Muons'), histTotal)
		
		hist3x3 = TH2D('histTight5x5',';#Delta#eta;#Delta#phi / rad;Fraction / %',nTotalBinsPhi,-lateralSpreadForPlot*HO_BIN - HO_BIN*phiBinOffsetFactor,
			lateralSpreadForPlot*HO_BIN + HO_BIN*phiBinOffsetFactor,
			nTotalBinsPhi,-lateralSpreadForPlot*HO_BIN - HO_BIN*phiBinOffsetFactor,
			lateralSpreadForPlot*HO_BIN + HO_BIN*phiBinOffsetFactor)
		
		fillGraphIn2DHist(self.fileHandler.getGraph('graphs/eAvAboveThr5x5patTightToL1Muons'), hist3x3)
		hist3x3.Divide(histTotal)
		hist3x3.SetStats(0)
		hist3x3.GetXaxis().SetRangeUser(-.3,.3)
		hist3x3.GetYaxis().SetRangeUser(-.3,.3)
		hist3x3.Draw('colz')
		c.Update()
		setupAxes(hist3x3)
		setupPalette(hist3x3, shiftBy=.05)
		hCopy = hist3x3.Clone('hCloneTight')
		hCopy.Scale(100)		
		hCopy.SetMarkerSize(2)
		hCopy.Draw('text,same')
		self.storeCanvas(c,'fractionAboveThrTight',marginRight=.15)

		return hist3x3,hCopy,c
