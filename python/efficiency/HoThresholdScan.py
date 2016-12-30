from plotting.Plot import Plot

from ROOT import TCanvas
from plotting.Utils import getTGraphErrors, getLegend,getMedian
import math
from plotting.PlotStyle import colorRwthDarkBlue, setupAxes, colorRwthMagenta

class HoThresholdScan(Plot):
	
	def __init__(self,filename,data = False, debug = False):
		Plot.__init__(self,filename,data,debug)
		self.createPlotSubdir('hoThresholdScan')
		
	def plotHoThresholdScan(self):
		canvas = TCanvas('cThresholdScan')
		canvas.SetLogy()
		xVals = []
		yVals = []
		yErr = []
		
		meanValues = []
		rmsValues = []
		
		for i in range(0,199):
			hist = self.fileHandler.getHistogram('multiplicity/recHitThrScan%d_Multiplicity' % i)
			hist.GetXaxis().SetRangeUser(0,50)
			hist.Draw()
			#self.storeCanvas(canvas,'plots/hoThresholdScan/hist%fGeV.gif' % ((i+1)*0.025))
			yVals.append(hist.GetBinCenter(hist.GetMaximumBin()))
			yErr.append(math.sqrt(yVals[-1]))
			xVals.append((i+1)*0.025)
			meanValues.append(hist.GetMean())
			rmsValues.append(hist.GetRMS())
			#getMedian(hist)
		graph = getTGraphErrors(xVals, yVals)#, ey=yErr)
		graph.SetMarkerStyle(20)
		graph.SetMarkerColor(colorRwthDarkBlue)
		graph.SetLineColor(colorRwthDarkBlue)
		graph.SetTitle(';E_{Thr} / GeV;# per Event')#Number of HORecHits per Event;E_{Thr} / GeV;# per Event')
		setupAxes(graph)
		graph.GetXaxis().SetRangeUser(0,2)
		graph.SetMinimum(6e-2)
		graph.Draw('ap')
		
		graphMean = getTGraphErrors(xVals, meanValues)#, ey=rmsValues)
		graphMean.SetMarkerStyle(21)
		graphMean.SetMarkerColor(colorRwthMagenta)
		graphMean.SetLineColor(colorRwthMagenta)
		graphMean.Draw('samep')
		
		legend = getLegend(x1=.48,x2=.98,y1=.8,y2=.95)
		legend.AddEntry(graph,'Most frequent # per Evt.','ep')
		legend.AddEntry(graphMean,'Mean # per Evt.','ep')
		legend.Draw()
				
		canvas.Update()
		self.storeCanvas(canvas,'hoThresholdScan',marginRight=.02)
		return canvas, graph,legend,graphMean