#!/usr/bin/python
from ROOT import TCanvas,TFile,ROOT,TPaveText,TGraph,TLegend
from plotting.Plot import Plot
#from efficiency.QualityCodes import plotQualityCodes,plotQualityCodesVsPt,plotQualityCodesStacked
import matplotlib.pyplot as plt
from cmath import sqrt
from plotting.PlotStyle import colorRwthMagenta, colorRwthGruen,\
	colorRwthDarkBlue, setupAxes, pyplotCmsPrivateLabel
from plotting.Colors import rwthDarkBlue, rwthGruen, rwthRot
from libPyROOT import Double

	
class GridMatching(Plot):
	def __init__(self,filename,data =False,debug=False):
		Plot.__init__(self,filename,data)
		self.createPlotSubdir('efficiency')
		plt.ioff()
		
	def get3x3MatchingPlot(self,hist3x3,hist3x3Tight,histTitle):
		eff3x3 = self.fileHandler.getHistogram(hist3x3)
		eff3x3Tight = self.fileHandler.getHistogram(hist3x3Tight)
		c = TCanvas(histTitle,histTitle,900,900)
		c.SetBottomMargin(.12)

		eff3x3.SetMarkerColor(colorRwthDarkBlue)
		eff3x3.SetLineColor(colorRwthDarkBlue)
		eff3x3.SetMarkerStyle(20)
		eff3x3.SetTitle(histTitle)

		eff3x3Tight.SetMarkerColor(colorRwthMagenta)
		eff3x3Tight.SetLineColor(colorRwthMagenta)
		eff3x3Tight.SetMarkerStyle(21)

		eff3x3.Draw('')
		eff3x3Tight.Draw('same')

		c.Update()
		
		setupAxes(eff3x3)
		setupAxes(eff3x3Tight)
		
		eff3x3.GetPaintedGraph().GetXaxis().SetTitleOffset(.8)
		eff3x3.GetPaintedGraph().GetYaxis().SetTitleOffset(.68)
		
		legend = TLegend(0.25,0.2,0.83,0.33)
		
		return c,eff3x3,eff3x3Tight,legend
		
	def plotL13x3AndL1Tight3x3(self):
		c,effL1Muon3x3,effL1TightMuon3x3,legend = self.get3x3MatchingPlot('efficiency/gridMatching_loose3x3_Efficiency',
			 'efficiency/gridMatching_tight3x3_Efficiency',
			  ';p_{T,RECO} / GeV;Efficiency')
				
		legend.AddEntry(effL1Muon3x3,'Matches in 3x3 grid','ep')
		legend.AddEntry(effL1TightMuon3x3,'Matches from tight in 3x3 grid','ep')
		legend.Draw()

		c.Update()
		self.storeCanvas(c,'efficiencyNormalAndTightVsPt_fullRange',marginRight=.02)
		
		effL1Muon3x3.GetPaintedGraph().GetXaxis().SetRangeUser(0,20)
		effL1Muon3x3.GetPaintedGraph().GetYaxis().SetRangeUser(0,1)
		
		
		for i in range(3,20):
			yL1 = effL1Muon3x3.GetPaintedGraph().Eval(i)
			yL1Tight = effL1TightMuon3x3.GetPaintedGraph().Eval(i)
			self.debug('L1: (%5.2f)\tTight: (%5.2f)\t Diff: %5.2f'% (yL1,yL1Tight,(yL1Tight-yL1)*100))
		
		c.Update()
		self.storeCanvas(c,'efficiencyNormalAndTightVsPt_zoom',marginRight=.02)
		
		return c, legend, effL1Muon3x3, effL1TightMuon3x3
	
	def plotL13x3AndL1Tight3x3L1Coordinates(self):
		
		c,effL1Muon3x3,effL1TightMuon3x3,legend = self.get3x3MatchingPlot('efficiency/gridMatching_L1pT_loose3x3_Efficiency',
			 'efficiency/gridMatching_L1pT_tight3x3_Efficiency',
			 ';p_{T,L1} / GeV;Efficiency')
				
		c.Update()
		self.storeCanvas(c,'efficiencyNormalAndTightVsPtL1_fullRange',marginRight=.02)	
			
		effL1Muon3x3.GetPaintedGraph().GetXaxis().SetRangeUser(0,20)
		effL1Muon3x3.GetPaintedGraph().GetYaxis().SetRangeUser(0,1)
		
		legend.AddEntry(effL1Muon3x3,'Matches in 3x3 grid','ep')
		legend.AddEntry(effL1TightMuon3x3,'Matches from tight in 3x3 grid','ep')
		legend.Draw()

		c.Update()
		self.storeCanvas(c,'efficiencyNormalAndTightVsPtL1_zoom',marginRight=.02)
		return c, legend, effL1Muon3x3, effL1TightMuon3x3

	def plotL13x3AndL1Tight3x3FromPatL1Pt(self):
		
		c,effL1Muon3x3,effL1TightMuon3x3,legend = self.get3x3MatchingPlot('efficiency/patToL1Muons_L1pT3x3_Efficiency',
			 'efficiency/patTightToL1Muons_L1pT3x3_Efficiency',
			 ';p_{T,L1} / GeV;Efficiency')
				
		c.Update()
		self.storeCanvas(c,'efficiencyFromPatNormalAndTightVsPtL1_fullRange',marginRight=.02)		
		
		effL1Muon3x3.GetPaintedGraph().GetXaxis().SetRangeUser(0,20)
		effL1Muon3x3.GetPaintedGraph().GetYaxis().SetRangeUser(0,1)
		
		legend.AddEntry(effL1Muon3x3,'Matches in 3x3 grid','ep')
		legend.AddEntry(effL1TightMuon3x3,'Matches from tight in 3x3 grid','ep')
		legend.Draw()

		c.Update()
		self.storeCanvas(c,'efficiencyFromPatNormalAndTightVsPtL1_zoom',marginRight=.02)
		return c, legend, effL1Muon3x3, effL1TightMuon3x3
	
	def plotL13x3AndL1Tight3x3FromPat(self):
		c,effL1Muon3x3,effL1TightMuon3x3,legend = self.get3x3MatchingPlot('efficiency/patToL1Muons3x3_Efficiency',
			 'efficiency/patTightToL1Muons3x3_Efficiency',
			 ';p_{T,RECO} / GeV;Efficiency')
				
		c.Update()
		self.storeCanvas(c,'efficiencyFromPatNormalAndTightVsPt_fullRange',marginRight=.02)		
		
		effL1Muon3x3.GetPaintedGraph().GetXaxis().SetRangeUser(0,20)
		effL1Muon3x3.GetPaintedGraph().GetYaxis().SetRangeUser(0,1)
		
		legend.AddEntry(effL1Muon3x3,'Matches in 3x3 grid','ep')
		legend.AddEntry(effL1TightMuon3x3,'Matches from tight in 3x3 grid','ep')
		legend.Draw()

		c.Update()
		self.storeCanvas(c,'efficiencyFromPatNormalAndTightVsPt_zoom',marginRight=.02)
		return c, legend, effL1Muon3x3, effL1TightMuon3x3
	
	def plotL1GridMatchingEfficiency(self):
		effL1MuonCentral = self.fileHandler.getHistogram('efficiency/L1MuonCentral_Efficiency')
		effL1Muon3x3 = self.fileHandler.getHistogram('efficiency/L1Muon3x3_Efficiency')
		effL1Muon5x5 = self.fileHandler.getHistogram('efficiency/L1Muon5x5_Efficiency')
		
		c = TCanvas("c1","c1",1200,1200)
		
		effL1Muon5x5.SetMarkerColor(colorRwthMagenta)
		effL1Muon5x5.SetLineColor(colorRwthMagenta)
		effL1Muon5x5.SetMarkerStyle(27)
		effL1Muon5x5.SetTitle('Efficiency for Matching HO to L1;p_{T} / GeV;Efficiency')
		
		effL1Muon5x5.Draw('')
		
		c.Update()
		effL1Muon5x5.GetPaintedGraph().GetYaxis().SetRangeUser(0.4,1)
		
		effL1Muon3x3.SetMarkerColor(colorRwthGruen)
		effL1Muon3x3.SetLineColor(colorRwthGruen)
		effL1Muon3x3.SetMarkerStyle(26)
		effL1Muon3x3.Draw("same")
		
		effL1MuonCentral.SetMarkerColor(colorRwthDarkBlue)
		effL1MuonCentral.SetLineColor(colorRwthDarkBlue)
		effL1MuonCentral.SetMarkerStyle(4)
		
		effL1MuonCentral.Draw("same")
		
		c.Update()
		
		setupAxes(effL1MuonCentral)
		setupAxes(effL1Muon3x3)
		setupAxes(effL1Muon5x5)
		
		legend = TLegend(0.55,0.1,0.9,0.3)
		legend.AddEntry(effL1MuonCentral,'Matches in central tile','ep')
		legend.AddEntry(effL1Muon3x3,'Matches in 3x3 grid','ep')
		legend.AddEntry(effL1Muon5x5,'Matches in 5x5 grid','ep')
		legend.Draw()
		
		effCentral 	= []
		eff3x3		= []
		eff5x5		= []
		xValues 	= []
		xErrLow 	= []
		xErrHigh 	= []
		yErrCentralHigh	= []
		yErrCentralLow	= []
		yErr3x3High	= []
		yErr3x3Low	= []
		yErr5x5High	= []
		yErr5x5Low	= []
		
		for i in range(effL1MuonCentral.GetPassedHistogram().GetNbinsX()):
			if effL1MuonCentral.GetTotalHistogram().GetBinContent(i) != 0:
				effCentral.append(effL1MuonCentral.GetPassedHistogram().GetBinContent(i)/effL1MuonCentral.GetTotalHistogram().GetBinContent(i)*100)
				xValues.append(effL1MuonCentral.GetTotalHistogram().GetBinCenter(i))
			if effL1Muon3x3.GetTotalHistogram().GetBinContent(i) != 0:
				eff3x3.append(effL1Muon3x3.GetPassedHistogram().GetBinContent(i)/effL1Muon3x3.GetTotalHistogram().GetBinContent(i)*100)
			if effL1Muon5x5.GetTotalHistogram().GetBinContent(i) != 0:
				eff5x5.append(effL1Muon5x5.GetPassedHistogram().GetBinContent(i)/effL1Muon5x5.GetTotalHistogram().GetBinContent(i)*100)
		
		for i in range(effL1MuonCentral.GetPaintedGraph().GetN()):
			yErrCentralHigh.append(effL1MuonCentral.GetPaintedGraph().GetErrorYhigh(i)*100)
			yErrCentralLow.append(effL1MuonCentral.GetPaintedGraph().GetErrorYlow(i)*100)
			yErr3x3High.append(effL1Muon3x3.GetPaintedGraph().GetErrorYhigh(i)*100)
			yErr3x3Low.append(effL1Muon3x3.GetPaintedGraph().GetErrorYlow(i)*100)
			yErr5x5High.append(effL1Muon5x5.GetPaintedGraph().GetErrorYhigh(i)*100)
			yErr5x5Low.append(effL1Muon5x5.GetPaintedGraph().GetErrorYlow(i)*100)
		
		xErrLow.append(0.25)
		xErrHigh.append(.25)
		for i in range(1,len(xValues) -1 ):
			xErrLow.append((xValues[i] - xValues[i-1])/2.)
			xErrHigh.append((xValues[i+1] - xValues[i])/2.)
		xErrLow.append(xErrLow[-1])
		xErrHigh.append(xErrHigh[-1])
		
		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		plt.errorbar(xValues,effCentral,yerr = [yErrCentralLow,yErrCentralHigh] ,xerr=[xErrLow,xErrHigh],fmt='s',color='#00549F',label='Match in central tile')
		plt.errorbar(xValues,eff3x3,yerr = [yErr3x3Low,yErr3x3High] ,xerr=[xErrLow,xErrHigh],fmt='o',color='#E30066',label='Match in 3x3')
		plt.errorbar(xValues,eff5x5,yerr = [yErr5x5Low,yErr5x5High] ,xerr=[xErrLow,xErrHigh],fmt='^',color='#57AB27',label='Match in 5x5')
		plt.gca().xaxis.grid(True)
		plt.gca().yaxis.grid(True)
		plt.axis([0.,155,40,101])
		plt.xlabel(r'p$_{\mathrm{T}}$ / GeV')
		plt.ylabel(r'$\epsilon$ / %')
		plt.title(r'Efficiency $\epsilon$ of matching L1 to HO')
		plt.legend(loc='lower right')
		pyplotCmsPrivateLabel(ax)
		plt.savefig('plots/efficiency.png')
	
		plt.show(block=False)
	
	def plotL1TruthGridMatchingPlot(self):
		effL1MuonCentral = self.fileHandler.getHistogram('efficiency/L1MuonTruthCentral_Efficiency')
		effL1Muon3x3 = self.fileHandler.getHistogram('efficiency/L1MuonTruth3x3_Efficiency')
		effL1Muon5x5 = self.fileHandler.getHistogram('efficiency/L1MuonTruth5x5_Efficiency')
		
		c = TCanvas("c12","c1",1200,1200)
		
		effL1Muon5x5.SetMarkerColor(colorRwthMagenta)
		effL1Muon5x5.SetLineColor(colorRwthMagenta)
		effL1Muon5x5.SetMarkerStyle(27)
		effL1Muon5x5.SetTitle('Efficiency for Matching HO to L1 truth;p_{T} / GeV;Efficiency')
		
		effL1Muon5x5.Draw('')
		
		c.Update()
		effL1Muon5x5.GetPaintedGraph().GetYaxis().SetRangeUser(0.4,1)
		
		effL1Muon3x3.SetMarkerColor(colorRwthGruen)
		effL1Muon3x3.SetLineColor(colorRwthGruen)
		effL1Muon3x3.SetMarkerStyle(26)
		effL1Muon3x3.Draw("same")
		
		effL1MuonCentral.SetMarkerColor(colorRwthDarkBlue)
		effL1MuonCentral.SetLineColor(colorRwthDarkBlue)
		effL1MuonCentral.SetMarkerStyle(4)
		
		effL1MuonCentral.Draw("same")
		
		c.Update()
		
		setupAxes(effL1MuonCentral)
		setupAxes(effL1Muon3x3)
		setupAxes(effL1Muon5x5)
		
		legend = TLegend(0.55,0.1,0.9,0.3)
		legend.AddEntry(effL1MuonCentral,'Matches in central tile','ep')
		legend.AddEntry(effL1Muon3x3,'Matches in 3x3 grid','ep')
		legend.AddEntry(effL1Muon5x5,'Matches in 5x5 grid','ep')
		legend.Draw()
		
		effCentral 	= []
		eff3x3		= []
		eff5x5		= []
		xValues 	= []
		xErrLow 	= []
		xErrHigh 	= []
		yErrCentralHigh	= []
		yErrCentralLow	= []
		yErr3x3High	= []
		yErr3x3Low	= []
		yErr5x5High	= []
		yErr5x5Low	= []
		
		for i in range(effL1MuonCentral.GetPassedHistogram().GetNbinsX()):
			if effL1MuonCentral.GetTotalHistogram().GetBinContent(i) != 0:
				effCentral.append(effL1MuonCentral.GetPassedHistogram().GetBinContent(i)/effL1MuonCentral.GetTotalHistogram().GetBinContent(i)*100)
				xValues.append(effL1MuonCentral.GetTotalHistogram().GetBinCenter(i))
			if effL1Muon3x3.GetTotalHistogram().GetBinContent(i) != 0:
				eff3x3.append(effL1Muon3x3.GetPassedHistogram().GetBinContent(i)/effL1Muon3x3.GetTotalHistogram().GetBinContent(i)*100)
			if effL1Muon5x5.GetTotalHistogram().GetBinContent(i) != 0:
				eff5x5.append(effL1Muon5x5.GetPassedHistogram().GetBinContent(i)/effL1Muon5x5.GetTotalHistogram().GetBinContent(i)*100)
		
		for i in range(effL1MuonCentral.GetPaintedGraph().GetN()):
			yErrCentralHigh.append(effL1MuonCentral.GetPaintedGraph().GetErrorYhigh(i)*100)
			yErrCentralLow.append(effL1MuonCentral.GetPaintedGraph().GetErrorYlow(i)*100)
			yErr3x3High.append(effL1Muon3x3.GetPaintedGraph().GetErrorYhigh(i)*100)
			yErr3x3Low.append(effL1Muon3x3.GetPaintedGraph().GetErrorYlow(i)*100)
			yErr5x5High.append(effL1Muon5x5.GetPaintedGraph().GetErrorYhigh(i)*100)
			yErr5x5Low.append(effL1Muon5x5.GetPaintedGraph().GetErrorYlow(i)*100)
		
		xErrLow.append(0.25)
		xErrHigh.append(.25)
		for i in range(1,len(xValues) -1 ):
			xErrLow.append((xValues[i] - xValues[i-1])/2.)
			xErrHigh.append((xValues[i+1] - xValues[i])/2.)
		xErrLow.append(xErrLow[-1])
		xErrHigh.append(xErrHigh[-1])
		
		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		plt.errorbar(xValues,effCentral,yerr = [yErrCentralLow,yErrCentralHigh] ,xerr=[xErrLow,xErrHigh],fmt='s',color='#00549F',label='Match in central tile')
		plt.errorbar(xValues,eff3x3,yerr = [yErr3x3Low,yErr3x3High] ,xerr=[xErrLow,xErrHigh],fmt='o',color='#E30066',label='Match in 3x3')
		plt.errorbar(xValues,eff5x5,yerr = [yErr5x5Low,yErr5x5High] ,xerr=[xErrLow,xErrHigh],fmt='^',color='#57AB27',label='Match in 5x5')
		plt.gca().xaxis.grid(True)
		plt.gca().yaxis.grid(True)
		plt.axis([0.,155,40,101])
		plt.xlabel(r'p$_{\mathrm{T}}$ / GeV')
		plt.ylabel(r'$\epsilon$ / %')
		plt.title(r'Efficiency $\epsilon$ of matching "true" L1 to HO')
		plt.legend(loc='lower right')
		pyplotCmsPrivateLabel(ax)
		plt.savefig('plots/efficiencyTruth.png')
		plt.show(block=False)
	
	
	def plot5x5GridTogether(self):
		c = TCanvas()
		effL1Muon5x5Truth = self.fileHandler.getHistogram('efficiency/L1MuonTruth5x5_Efficiency')
		effL1Muon5x5 = self.fileHandler.getHistogram('efficiency/L1Muon5x5_Efficiency')
	
		effL1Muon5x5.Draw()
		effL1Muon5x5Truth.Draw("same")
		
		c.Update()
		
		xValues 	= []
		xErrLow 	= []
		xErrHigh 	= []
		eff5x5		= []
		yErr5x5High	= []
		yErr5x5Low	= []
		effTruth5x5		= []
		yErr5x5TruthHigh	= []
		yErr5x5TruthLow	= []
		
		diff = []
		diffErrHigh = []
		diffErrLow = []
		
		for i in range(effL1Muon5x5.GetPassedHistogram().GetNbinsX()):
			if effL1Muon5x5.GetTotalHistogram().GetBinContent(i) != 0:
				eff5x5.append(effL1Muon5x5.GetPassedHistogram().GetBinContent(i)/effL1Muon5x5.GetTotalHistogram().GetBinContent(i)*100)
				xValues.append(effL1Muon5x5.GetTotalHistogram().GetBinCenter(i))
			if effL1Muon5x5Truth.GetTotalHistogram().GetBinContent(i) != 0:
				effTruth5x5.append(effL1Muon5x5Truth.GetPassedHistogram().GetBinContent(i)/effL1Muon5x5Truth.GetTotalHistogram().GetBinContent(i)*100)
	
		for i in range(effL1Muon5x5.GetPaintedGraph().GetN()):
			yErr5x5High.append(effL1Muon5x5.GetPaintedGraph().GetErrorYhigh(i)*100)
			yErr5x5Low.append(effL1Muon5x5.GetPaintedGraph().GetErrorYlow(i)*100)
			yErr5x5TruthHigh.append(effL1Muon5x5Truth.GetPaintedGraph().GetErrorYhigh(i)*100)
			yErr5x5TruthLow.append(effL1Muon5x5Truth.GetPaintedGraph().GetErrorYlow(i)*100)
		
		
		xErrLow.append(0.25)
		xErrHigh.append(.25)
		
		for i in range(0,len(eff5x5)):
			diff.append(effTruth5x5[i] - eff5x5[i])
			diffErrHigh.append(sqrt(pow(yErr5x5High[i], 2) + pow(yErr5x5TruthHigh[i], 2)))
			diffErrLow.append(sqrt(pow(yErr5x5Low[i], 2) + pow(yErr5x5TruthLow[i], 2)))
		
		for i in range(1,len(xValues) -1 ):
			xErrLow.append((xValues[i] - xValues[i-1])/2.)
			xErrHigh.append((xValues[i+1] - xValues[i])/2.)
		xErrLow.append(xErrLow[-1])
		xErrHigh.append(xErrHigh[-1])
		
		fig = plt.figure()
		ax1 = fig.add_subplot(211)
		ax2 = fig.add_subplot(212,sharex=ax1)
	
		ax1.errorbar(xValues,eff5x5,yerr = [yErr5x5Low,yErr5x5High] ,xerr=[xErrLow,xErrHigh],fmt='^',color='#57AB27',label='L1')
		ax1.errorbar(xValues,effTruth5x5,yerr = [yErr5x5TruthLow,yErr5x5TruthHigh] ,xerr=[xErrLow,xErrHigh],fmt='^',color='#E30066',label='Truth L1')
		ax1.legend(loc='lower right')
		ax1.grid(True)
		ax1.set_ylabel(r'$\epsilon$ / %')
		ax1.set_title(r'Efficiency $\epsilon$ for 5x5 grid in direct comparison')
	
		pyplotCmsPrivateLabel(ax1)
	
		ax2.errorbar(xValues,diff,yerr=[diffErrLow,diffErrHigh], xerr=[xErrLow,xErrHigh],fmt='o', color='#00549F', label=r'$\epsilon$(Truth L1) - $\epsilon$(L1)')
		ax2.legend()
		ax2.grid(True)
		ax2.set_ylabel(r'$\Delta$ $\epsilon$ / %')
		ax2.set_xlabel(r'p$_\mathrm{T}$ / GeV')
	
		#make x labels of first plot invisible
		plt.setp(ax1.get_xticklabels(),visible=False)
		
		#avoid overlapping labels on y axis
		plt.setp(ax1.get_yticklabels()[0],visible=False)
		plt.setp(ax2.get_yticklabels()[-1],visible=False)
	
		fig.subplots_adjust(hspace=0)
		plt.savefig('plots/efficiency5x5.png')
		plt.show(block=False)
	
	def plotDeltaNL1ComparedGridMatching(self):
		effL1Muon3x3Truth = self.fileHandler.getHistogram('efficiency/L1MuonTruth3x3_Efficiency')
		effL1Muon3x3 = self.fileHandler.getHistogram('efficiency/L1Muon3x3_Efficiency')
	
		xValues = []
		xErrLow = []
		xErrHigh= []
		yValues = []
		
		for i in range(effL1Muon3x3.GetPassedHistogram().GetNbinsX()):
			if(effL1Muon3x3.GetTotalHistogram().GetBinContent(i) != 0):
				yValues.append(effL1Muon3x3Truth.GetTotalHistogram().GetBinContent(i) - effL1Muon3x3.GetTotalHistogram().GetBinContent(i))
				xValues.append(effL1Muon3x3.GetTotalHistogram().GetBinCenter(i))
		
		xErrLow.append(0.25)
		xErrHigh.append(.25)
		
		for i in range(1,len(xValues) -1 ):
			xErrLow.append((xValues[i] - xValues[i-1])/2.)
			xErrHigh.append((xValues[i+1] - xValues[i])/2.)
		xErrLow.append(xErrLow[-1])
		xErrHigh.append(xErrHigh[-1])
		
		fig = plt.figure()
		ax1 = fig.add_subplot(111)
		ax1.errorbar(xValues,yValues,xerr=[xErrLow,xErrHigh],fmt='^',color='#57AB27',label='# Truth - # All')
		ax1.legend(loc='lower right')
		ax1.grid(True)
		ax1.set_ylabel('# L1')
		ax1.set_xlabel(r'p$_\mathrm{T}$ / GeV')
		ax1.set_title(r'# L1 in direct comparison for Grid Matching (3x3)')
		pyplotCmsPrivateLabel(ax1)
		plt.savefig('plots/nL13x3.png')
		plt.show(block=False)
	
	def plotDeltaNL1ComparedMatchingByEMax(self):
		effL1Muon = self.fileHandler.getHistogram("efficiency/GenAndL1MuonPt15_Efficiency")
		effL1MuonTruth = self.fileHandler.getHistogram("efficiency/GenAndL1MuonAndHoAboveThrPt15_Efficiency")
		xValues = []
		xErrLow = []
		xErrHigh= []
		yValues = []
		
		for i in range(effL1Muon.GetPassedHistogram().GetNbinsX()):
			if(effL1Muon.GetTotalHistogram().GetBinContent(i) != 0):
				yValues.append(effL1MuonTruth.GetTotalHistogram().GetBinContent(i) - effL1Muon.GetTotalHistogram().GetBinContent(i))
				xValues.append(effL1Muon.GetTotalHistogram().GetBinCenter(i))
		
		xErrLow.append(0.25)
		xErrHigh.append(.25)
		
		for i in range(1,len(xValues) -1 ):
			xErrLow.append((xValues[i] - xValues[i-1])/2.)
			xErrHigh.append((xValues[i+1] - xValues[i])/2.)
		xErrLow.append(xErrLow[-1])
		xErrHigh.append(xErrHigh[-1])
		
		fig = plt.figure()
		ax1 = fig.add_subplot(111)
		ax1.errorbar(xValues,yValues,xerr=[xErrLow,xErrHigh],fmt='^',color='#00549F',label='(L1 + HO) - L1')
		ax1.legend(loc='lower right')
		ax1.grid(True)
		ax1.set_ylabel('# L1')
		ax1.set_xlabel(r'p$_\mathrm{T}$ / GeV')
		ax1.set_title(r'# L1 in direct comparison for $\Delta$R matching by E$_{Max}$')
		pyplotCmsPrivateLabel(ax1)
		plt.savefig('plots/nL1ByEmax.png')
		plt.show(block=False)
		
	def plotNtotalGridMatching3x3(self):
		effL1Muon3x3 = self.fileHandler.getHistogram("efficiency/L1Muon3x3_Efficiency")
		effL1Muon3x3Truth = self.fileHandler.getHistogram("efficiency/L1MuonTruth3x3_Efficiency")
		effL1Muon5x5Truth = self.fileHandler.getHistogram("efficiency/L1MuonTruth5x5_Efficiency")
	
		genPt = self.fileHandler.getHistogram("gen_Pt")
		
		yTruth = []
		yTruthErr = []
		
		xValues = []
		xErrLow = []
		xErrHigh= []
		yValues = []
		
		xGen = []
		xGenErr = []
		yGen = []
		yGenSummed = []
		yGenSummedErr = []
		
		for i in range(1,genPt.GetNbinsX()):
			if(genPt.GetBinCenter(i) > 200):
				break
			xGen.append(genPt.GetBinCenter(i))
			xGenErr.append(genPt.GetBinWidth(i)/2.)
			yGen.append( genPt.GetBinContent(i) )
			
		for i in range(1,effL1Muon3x3.GetPassedHistogram().GetNbinsX()):
			if(effL1Muon3x3.GetTotalHistogram().GetBinContent(i) != 0):
				yValues.append(effL1Muon3x3.GetTotalHistogram().GetBinContent(i))
				xValues.append(effL1Muon3x3.GetTotalHistogram().GetBinCenter(i))
				yTruth.append(effL1Muon3x3Truth.GetTotalHistogram().GetBinContent(i))
				yTruthErr.append(sqrt(yTruth[-1]))
		
		print effL1Muon3x3.GetTotalHistogram().Integral(),effL1Muon3x3Truth.GetTotalHistogram().Integral()\
			,genPt.Integral(),effL1Muon5x5Truth.GetTotalHistogram().Integral()
	
		
		xErrLow.append(0.25)
		xErrHigh.append(.25)
		
		for i in range(1,len(xValues) -1 ):
			xErrLow.append((xValues[i] - xValues[i-1])/2.)
			xErrHigh.append((xValues[i+1] - xValues[i])/2.)
		xErrLow.append(xErrLow[-1])
		xErrHigh.append(60)#xErrHigh[-1])
	
		counter = 0
		genSum = 0
		for i in range(len(xGen)):
			if(xGen[i]>(xValues[counter] + xErrHigh[counter])):
				counter+=1
				yGenSummed.append(genSum)
				yGenSummedErr.append(sqrt(yGenSummed[-1]))
	
				genSum = 0
			genSum += yGen[i]
		yGenSummed.append(genSum)
		yGenSummedErr.append(sqrt(yGenSummed[-1]))
	
		for i in range(0,len(yValues)):
			yValues[i] /= (xErrHigh[i] + xErrLow[i])
		for i in range(0,len(yTruth)):
			yTruth[i] /= (xErrHigh[i] + xErrLow[i])
			yTruthErr[i] /= (xErrHigh[i] + xErrLow[i])
		for i in range(0,len(yGenSummed)):
			yGenSummed[i] /= (xErrHigh[i] + xErrLow[i])
			yGenSummedErr[i] /= (xErrHigh[i] + xErrLow[i])
		
		fig = plt.figure()
		ax1 = fig.add_subplot(111)
		ax1.errorbar(xValues,yValues,xerr=[xErrLow,xErrHigh],fmt='^',color=rwthDarkBlue,label='L1')
		ax1.set_yscale('log')
		ax1.errorbar(xValues,yTruth,xerr=[xErrLow,xErrHigh],yerr=yTruthErr,fmt='^',color=rwthGruen,label='L1 Truth')
	#	ax1.errorbar(xGen,yGen,xerr=xGenErr,fmt='^',color=rwthMagenta,label='Gen')
		ax1.errorbar(xValues,yGenSummed,xerr=[xErrLow,xErrHigh],yerr = yGenSummedErr,fmt='^',color=rwthRot,label='Gen Summed')
		ax1.legend(loc='lower right')
		ax1.grid(True)
		ax1.set_ylabel('# / bin width (1/GeV)')
		ax1.set_xlabel(r'p$_\mathrm{T}$ / GeV')
		ax1.set_title('# L1')
		pyplotCmsPrivateLabel(ax1)
		plt.savefig('plots/nL1ByEmax3x3Absolute.png')
		plt.show(block=False)
		
	def plot3x3GridTogether(self):
		c = TCanvas()
		effL1Muon3x3Truth = self.fileHandler.getHistogram('efficiency/L1MuonTruth3x3_Efficiency')
		effL1Muon3x3 = self.fileHandler.getHistogram('efficiency/L1Muon3x3_Efficiency')
	
		effL1Muon3x3.Draw()
		effL1Muon3x3Truth.Draw("same")
		
		c.Update()
		
		xValues 	= []
		xErrLow 	= []
		xErrHigh 	= []
		eff5x5		= []
		yErr3x3High	= []
		yErr3x3Low	= []
		effTruth3x3		= []
		yErr3x3TruthHigh	= []
		yErr3x3TruthLow	= []
		
		diff = []
		diffErrHigh = []
		diffErrLow = []
		
		for i in range(effL1Muon3x3.GetPassedHistogram().GetNbinsX()):
			if effL1Muon3x3.GetTotalHistogram().GetBinContent(i) != 0:
				eff5x5.append(effL1Muon3x3.GetPassedHistogram().GetBinContent(i)/effL1Muon3x3.GetTotalHistogram().GetBinContent(i)*100)
				xValues.append(effL1Muon3x3.GetTotalHistogram().GetBinCenter(i))
			if effL1Muon3x3Truth.GetTotalHistogram().GetBinContent(i) != 0:
				effTruth3x3.append(effL1Muon3x3Truth.GetPassedHistogram().GetBinContent(i)/effL1Muon3x3Truth.GetTotalHistogram().GetBinContent(i)*100)
	
		for i in range(effL1Muon3x3.GetPaintedGraph().GetN()):
			yErr3x3High.append(effL1Muon3x3.GetPaintedGraph().GetErrorYhigh(i)*100)
			yErr3x3Low.append(effL1Muon3x3.GetPaintedGraph().GetErrorYlow(i)*100)
			yErr3x3TruthHigh.append(effL1Muon3x3Truth.GetPaintedGraph().GetErrorYhigh(i)*100)
			yErr3x3TruthLow.append(effL1Muon3x3Truth.GetPaintedGraph().GetErrorYlow(i)*100)
		
		
		xErrLow.append(0.25)
		xErrHigh.append(.25)
		
		for i in range(0,len(eff5x5)):
			diff.append(effTruth3x3[i] - eff5x5[i])
			diffErrHigh.append(sqrt(pow(yErr3x3High[i], 2) + pow(yErr3x3TruthHigh[i], 2)))
			diffErrLow.append(sqrt(pow(yErr3x3Low[i], 2) + pow(yErr3x3TruthLow[i], 2)))
		
		for i in range(1,len(xValues) -1 ):
			xErrLow.append((xValues[i] - xValues[i-1])/2.)
			xErrHigh.append((xValues[i+1] - xValues[i])/2.)
		xErrLow.append(xErrLow[-1])
		xErrHigh.append(xErrHigh[-1])
		
		fig = plt.figure()
		ax1 = fig.add_subplot(211)
		ax2 = fig.add_subplot(212,sharex=ax1)
	
		ax1.errorbar(xValues,eff5x5,yerr = [yErr3x3Low,yErr3x3High] ,xerr=[xErrLow,xErrHigh],fmt='^',color='#57AB27',label='L1')
		ax1.errorbar(xValues,effTruth3x3,yerr = [yErr3x3TruthLow,yErr3x3TruthHigh] ,xerr=[xErrLow,xErrHigh],fmt='^',color='#E30066',label='Truth L1')
		ax1.legend(loc='lower right')
		ax1.grid(True)
		ax1.set_ylabel(r'$\epsilon$ / %')
		ax1.set_title(r'Efficiency $\epsilon$ for 3x3 grid in direct comparison')
	
		pyplotCmsPrivateLabel(ax1)
	
		ax2.errorbar(xValues,diff,yerr=[diffErrLow,diffErrHigh], xerr=[xErrLow,xErrHigh],fmt='o', color='#00549F', label=r'$\epsilon$(Truth L1) - $\epsilon$(L1)')
		ax2.legend()
		ax2.grid(True)
		ax2.set_ylabel(r'$\Delta$ $\epsilon$ / %')
		ax2.set_xlabel(r'p$_\mathrm{T}$ / GeV')
	
		#make x labels of first plot invisible
		plt.setp(ax1.get_xticklabels(),visible=False)
		
		#avoid overlapping labels on y axis
		plt.setp(ax1.get_yticklabels()[0],visible=False)
		plt.setp(ax2.get_yticklabels()[-1],visible=False)
	
		fig.subplots_adjust(hspace=0)
		plt.savefig('plots/efficiency3x3.png')
		plt.show(block=False)
	
