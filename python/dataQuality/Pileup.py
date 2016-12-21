import os
from plotting.Plot import Plot
from plotting.OutputModule import CliColors
from plotting.PlotStyle import setupPalette, colorRwthDarkBlue, colorRwthMagenta,\
	setupAxes, getTH1D, setStatBoxPosition
from ROOT import TCanvas,TPaveText,TFile,gApplication

class Pileup(Plot):
	def __init__(self,filename = "MyDataPileupHistogram.root", data=True, debug = False):
		Plot.__init__(self,filename,data,debug)
		self.fileHandler.setModuleName('')
		self.createPlotSubdir('pileup')
		
	def plotPileup(self):
		c = TCanvas('cPu')
		rootfile = TFile(self.filename)
		if self.DEBUG:
			self.warning(rootfile)
		puOriginal = rootfile.Get('pileup').Clone()
		puNew = getTH1D('histPu',';Interactions per BX;a.u.',26,-.5,25.5)
		
		for x in range(1,puNew.GetNbinsX()+1):
			puNew.SetBinContent(x,puOriginal.GetBinContent(x))
		
		puNew.SetLineColor(colorRwthDarkBlue)	
		puNew.SetFillColor(colorRwthDarkBlue)	
		puNew.Scale(1/puNew.Integral())
		puNew.Draw()
		
		c.Update()
		setStatBoxPosition(puNew, x1=.78,x2=.98,y2=.95)
		self.storeCanvas(c,'pileup',drawMark=False, marginRight=.02)
		
		gApplication.Run()
		return c, puNew