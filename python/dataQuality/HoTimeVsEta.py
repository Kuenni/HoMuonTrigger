'''
Created on Jul 25, 2016

@author: kuensken
'''
from plotting.Plot import Plot
from plotting.Utils import getMedian, calcPercent, fillGraphIn2DHist, calcSigma,\
	getLegend, fill2DGraphIn2DHist
from ROOT import Double, TH1D, TCanvas, TEfficiency, TF1, TH2D, TBox, TGaxis,gPad,gStyle,TPaveText
from plotting.PlotStyle import setupAxes, colorRwthMagenta, colorRwthRot,\
	colorRwthDarkBlue, colorRwthDarkGray, setupPalette
from math import fabs, sqrt
import numpy as np

class HoTimeVsEta(Plot):
	'''
	classdocs
	'''
	def __init__(self,filename,data,debug):
		Plot.__init__(self,filename,data,debug)
		self.createPlotSubdir('timing')
		
	def printFractionsPerIEta(self,graph):
		self.output("Fractions of HO time in [-12.5,12.5] ns")
		counterDict = [{'total':0,'inside':0, 'hist':TH1D('hist' + str(i-10),'',201,-100.5,100.5)} for i in range(0,21)]
		x = Double(0)
		y = Double(0)
		nTotal = graph.GetN()
		
		###
		# Fill histograms for later calculation of mean
		###
		for i in range(0,nTotal):
			graph.GetPoint(i,x,y)
			indexHelper = int(x+10)
			counterDict[indexHelper]['total'] += 1
			counterDict[indexHelper]['hist'].Fill(y)
			
		
		###
		# Once all histograms are filled, calculate median
		###
		for index,item in enumerate(counterDict):
			item['median'] = getMedian(item['hist'])
		
		###
		# Fill the number of objects in interval depending on the
		# Histogram median
		###
		nDone = 0
		if self.DEBUG:
			self.debug('Filling interval counters')
		for i in range(0,nTotal):
			graph.GetPoint(i,x,y)
			indexHelper = int(x+10)
			median = counterDict[indexHelper]['median']
			if( fabs(median - y) < 12.5 ):
				counterDict[indexHelper]['inside'] += 1
			nDone += 1
			self.printProgress(nDone, nTotal)
			
		#Graph for results
		graph = TEfficiency(graph.GetName(),"",21,-10.5,10.5)
		
		###
		# Plot the results of calculations on CLI
		###
		for index,item in enumerate(counterDict):
			if(index - 10 == 0):
				continue
			total = item['total']
			inside = item['inside']
			self.output("iEta: %3d\tTotal: %5d\tInside:%5d\tFraction:%6.2f +/- %6.2f" % 
					(index - 10, total, inside, calcPercent(inside,total),calcSigma(inside, total)*100))
			graph.SetTotalEvents(graph.FindFixBin(index -10),total)
			graph.SetPassedEvents(graph.FindFixBin(index -10),inside)
		return graph,counterDict

	
	def plotFractionsVsEta(self,graph,title,saveName):
		c2 = TCanvas(saveName + '_EtaFractions',"saveName",0,600,800,600)
		c2.cd().SetLeftMargin(.15)
		graph.SetTitle(';i#eta;rel. Fraction in [-12.5,12.5]ns')
		graph.Draw('ap')
		c2.Update()
		setupAxes(graph)
		self.storeCanvas(c2,saveName,marginLeft=.15,marginRight=.02)
		return c2,graph
		
	### ========================
	### Plots of HO time vs iEta
	### ========================
	def makeTimeVsEtaPlot(self,source,title = ""):
		canvas = TCanvas(source)
		canvas.SetRightMargin(.15)
		if title == "":
			title = source
		hist = TH2D(source,";i#eta;Time / ns;Entries",33,-16.5,16.5,
				201,-100.5,100.5)
		graph = self.fileHandler.getGraph('graphs/timingSupport_' + source)
		fillGraphIn2DHist(graph, hist)
		hist.SetStats(0)
		hist.GetXaxis().SetRangeUser(-11,11)
		hist.GetYaxis().SetRangeUser(-60,60)
		hist.Draw('colz')
		canvas.Update()
		setupAxes(hist)
		setupPalette(hist,shiftBy=.05)
		canvas.Update()
		fractionGraph,counterDict = self.printFractionsPerIEta(graph)
		self.debug('Integral of plot %20s: %d' % (source,hist.Integral()))
		medianTofZero = (counterDict[9]['median'] + counterDict[11]['median'])/2.
		tofFunction = TF1('f','4*sqrt(1+ 1/(tan(2*atan(exp(-x*0.087/2.)))**2))/300000000.*1e9 - 13.3 + [0]',-10,10)
		tofFunction.SetParameter(0,medianTofZero)
		tofFunction.Draw('same')

		###
		# Draw Boxes for the Intervals
		###
		intervalBoxes = []
		for index,item in enumerate(counterDict):
			iEta = index - 10
			# Skip iEta 0
			if iEta == 0:
				continue
			box = TBox(iEta - 0.5, item['median'] - 12.5, iEta + 0.5, item['median'] + 12.5)
			box.SetFillStyle(0)
			box.SetLineColor(colorRwthMagenta)
			box.SetLineWidth(2)
			box.Draw()
			intervalBoxes.append(box)

		return canvas,hist,fractionGraph,intervalBoxes,tofFunction
	
	### ========================
	### Plots of HO time vs iPhi
	### ========================
	def makeTimeVsPhiPlot(self,source,iEta,title = ""):
		canvas = TCanvas(source + str(iEta),str(iEta))
		if title == "":
			title = source
		hist = TH2D(source + str(iEta),title + str(iEta) + ";i#phi;Time / ns;#",72,.5,72.5,
				201,-100.5,100.5)
		graph = self.fileHandler.getGraph('graphs/iEta/timingSupport_' + source + 'Ieta' + str(iEta))
		if not graph: return
		fillGraphIn2DHist(graph, hist)
		hist.SetStats(0)
		hist.Draw('colz')
		hist.GetYaxis().SetRangeUser(-60,60)
		canvas.Update()
		setupAxes(hist)
		canvas.Update()

		self.storeCanvas(canvas, 'timeVsPhi/' +  source + '_iEta' + str(iEta),
						labelPosition={'x1ndc' : .6, 'y1ndc' : 0.92, 'x2ndc' : .9, 'y2ndc' : 0.95} )

		return canvas,hist,graph
	
	#
	# Plot all iEta bins, time vs iphi
	#
	def plotTimeVsPhi(self):
		results = []
		for i in range(-12,13):
			results.append(self.makeTimeVsPhiPlot("UnmatchedDtHoIphiTime",i,'Unmatched DT i#eta '))
		return results
	
	#
	# Plot all iEta bins, time vs iphi tight
	#
	def plotTimeVsPhiTight(self):
		results = []
		for i in range(-12,13):
			results.append(self.makeTimeVsPhiPlot("tight_UnmatchedDtHoIphiTime",i,'#verb{tight} Unmatched DT i#eta '))
		return results
	
	#
	# Plot all iEta bins, time vs iphi DT/RPC
	#
	def plotTimeVsPhiDtRpc(self):
		results = []
		for i in range(-12,13):
			results.append(self.makeTimeVsPhiPlot("MatchedDtRpcHoIphiTime",i,'DT/RPC i#eta '))
		return results	
	
	### ====================
	### Plot L1 BX ID vs eta
	### ====================
	def makeL1TimeVsEtaPlot(self,source):
		canvas = TCanvas(source,source)
		canvas.SetLogz()
		hist = TH2D(source,source + ";#eta_{L1};BXID;#",20,-1,1,
				7,-3.5,3.5)
		graph = self.fileHandler.getGraph('graphs/timingSupport_' + source)
		fillGraphIn2DHist(graph, hist)
		hist.SetStats(0)
		hist.Draw('colz')
		canvas.Update()
		setupAxes(hist)
		label = self.drawLabel()
		canvas.Update()
		return canvas,label,hist

	### ============================================
	### Plot fraction vs. Eta for HO and L1 together
	### ============================================
	def makeCombinedEtaPlot(self, tight = False):
		hist = self.makeL1TimeVsEtaPlot(('tight_' if tight else '') + 'dtOnly_bxidVsEta')[2]
		countsInL1 = []
		for x in np.arange(-.95,1.05,0.1):
			totalCounter = 0
			zeroCount = 0
			for y in range(-2,3):
				totalCounter += hist.GetBinContent(hist.FindBin(x,y))
				if y == 0:
					zeroCount = hist.GetBinContent(hist.FindBin(x,y))
			countsInL1.append({'total':totalCounter,'zero':zeroCount,'eta':x})
	
		#Graph for results
		graph1 = TEfficiency(hist.GetName(),"",8,-9.195,-.5)
		graph2 = TEfficiency(hist.GetName(),"",8,.5,9.195)
	
		for item in countsInL1:
			if item['total'] == 0:
				continue
			print item['total'],item['zero'],item['eta']
			if item['eta'] < 0:
				graph1.SetTotalEvents(graph1.FindFixBin(-0.5 + item['eta']/0.087),int(item['total']))
				graph1.SetPassedEvents(graph1.FindFixBin(-0.5 + item['eta']/0.087),int(item['zero']))
			else:
				graph2.SetTotalEvents(graph2.FindFixBin(0.5 + item['eta']/0.087),int(item['total']))
				graph2.SetPassedEvents(graph2.FindFixBin(0.5 + item['eta']/0.087),int(item['zero']))
		
		histHo = None
		if tight:
			histHo = self.plotTightHoTimeVsEta()[2][1]
		else:
			histHo = self.plotHoTimeVsEta()[2][1]
			
		histHo.SetTitle('')
		
		canvas = TCanvas('combinedPlot' + ('Tight ' if tight else '') + hist.GetName(),'combinedPlot')
		histHo.Draw('ap')
		canvas.Update()
		canvas.cd().SetTicks(0,0)
		
		histHo.SetMarkerStyle(2)
		histHo.SetLineColor(colorRwthDarkBlue)
		histHo.SetMarkerColor(colorRwthDarkBlue)
		histHo.GetPaintedGraph().GetXaxis().SetRangeUser(-12,12)
		histHo.GetPaintedGraph().GetXaxis().SetLabelColor(colorRwthDarkBlue)
		histHo.GetPaintedGraph().GetXaxis().SetTitleColor(colorRwthDarkBlue)
		histHo.GetPaintedGraph().GetXaxis().SetAxisColor(colorRwthDarkBlue)
		yMax = gPad.GetFrame().GetY2()
		yMin = gPad.GetFrame().GetY1()
		
		#Print average Fraction excluding iEta +/-2
		x = Double(0)
		y = Double(0)
		mean = 0
		var = 0
		for i in range(0,histHo.GetPaintedGraph().GetN()):
			histHo.GetPaintedGraph().GetPoint(i,x,y)
			if abs(x) == 2:
				continue
			mean += y
			var += histHo.GetPaintedGraph().GetErrorY(i)*histHo.GetPaintedGraph().GetErrorY(i)
			
		mean /= histHo.GetPaintedGraph().GetN() - 2
		sigma = sqrt(var / (histHo.GetPaintedGraph().GetN() - 2))

		self.debug("Average fraction excluding iEta +/- 2 %s: %5.2f%% +/- %5.2f%%" % ('[Tight]' if tight else '',mean*100,sigma*100))

		nTotal = 0
		nPassed = 0
		for item in countsInL1:
			if fabs(item['eta']) == 2 or fabs(item['eta'] == 0):
				continue
			nTotal += item['total']
			nPassed += item['zero']

		#Print again with ClopperPearson uncertainty, the counts are for L1!!!
		mean = nPassed/nTotal*100
		sigmaPlus = TEfficiency.ClopperPearson(int(nTotal),int(nPassed),.68,1)*100 - mean
		sigmaMinus = mean - TEfficiency.ClopperPearson(int(nTotal),int(nPassed),.68,0)*100
		#self.debug("Average fraction excluding iEta +/- 2 %s with Clop.Pear.: %5.2f%% +%5.2f%% -%5.2f%%" 
		#		% ('[Tight]' if tight else '',mean,sigmaPlus,sigmaMinus))
	
		
		#Left axis part
		f1 = TF1("f1","x",-0.87,0)
		A1 = TGaxis(-10,yMax,-0.5,yMax,"f1",010,"-")
		A1.SetLineColor(colorRwthRot)
		A1.SetLabelColor(colorRwthRot)
		A1.Draw()
		
		#Right axis part
		f2 = TF1("f2","x",0,0.87)
		A2 = TGaxis(0.5,yMax,10,yMax,"f2",010,"-")
		A2.SetLineColor(colorRwthRot)
		A2.SetLabelColor(colorRwthRot)
		A2.Draw()
		
		#Box for shading out 0
		box = TBox(-.5,yMin,0.5, yMax)
		box.SetLineColor(colorRwthDarkGray)
		box.SetFillColor(colorRwthDarkGray)
		box.SetFillStyle(3013)
		box.Draw('same')
		
		#Left L1 eta
		graph1.SetMarkerColor(colorRwthRot)
		graph1.SetLineColor(colorRwthRot)
		graph1.SetMarkerStyle(20)
		graph1.Draw('same,p')
		
		#Right L1Eta
		graph2.SetMarkerColor(colorRwthRot)
		graph2.SetLineColor(colorRwthRot)
		graph2.SetMarkerStyle(20)
		graph2.Draw('same,p')
		
		#Label for extra axis
		axisLabel = TPaveText( 0.83,0.94,0.9,0.98,"NDC")
		axisLabel.AddText('#eta_{L1}')
		axisLabel.SetBorderSize(0)
		axisLabel.SetFillStyle(0)
		axisLabel.SetTextColor(colorRwthRot)
		axisLabel.SetTextSize(0.06)
		axisLabel.Draw()
		
		#Legend
		legend = getLegend(x1=0.1,y1=0.1,x2=0.4,y2=.35)
		legend.AddEntry(histHo,'HO #in #pm12.5 ns','pe')
		legend.AddEntry(graph1,('Tight ' if tight else '') + 'L1 BXID = 0','pe')
		legend.Draw()
		
		label = self.drawLabel({'x1ndc' : .1, 'y1ndc' : 0.35, 'x2ndc' : .4, 'y2ndc' : 0.4})
		
		canvas.Update()
		self.storeCanvas(canvas, "combinedFractionL1AndHo" + ('Tight' if tight else ''),drawMark=False,drawLabel=False,marginTop=0.08)
			
		return histHo, graph1,canvas,A1,f1,A2,f2,box,graph2,axisLabel,legend
	
	### ============================
	### Predefined plotter functions
	### ============================
	
	def plotHoTimeVsEta(self):
		canvas, hist, graph, boxes, tofFunction = self.makeTimeVsEtaPlot('UnmatchedDtHoTimeGraph', "Unmatched DT")
		self.storeCanvas(canvas,'UnmatchedDtHo_TimeVsEta',marginRight=.15)
		fractionsPerEtaData = self.plotFractionsVsEta(graph,"Fraction of HORecHits at BXID 0",'UnmatchedDtHoTimeGraph_fractionVsEta')
		return canvas,hist,fractionsPerEtaData, boxes, tofFunction
	
	def plotHoTimeVsEtaBxWrong(self):
		canvas, hist, graph, boxes, tofFunction = self.makeTimeVsEtaPlot('UnmatchedDtHoBxNot0TimeGraph',"Unmatched DT BX Wrong")
		self.storeCanvas(canvas,'UnmatchedDtHoBxNot0TimeGraph_TimeVsEta',marginRight=.15)
		fractionsPerEtaData = self.plotFractionsVsEta(graph,"Fraction of HORecHits at BXID 0,L1 BX wrong",'UnmatchedDtHoBxNot0TimeGraph_fractionVsEta')
		return canvas,hist,fractionsPerEtaData, boxes, tofFunction
	
	def plotTightHoTimeVsEta(self):
		canvas, hist, graph, boxes, tofFunction = self.makeTimeVsEtaPlot('tight_UnmatchedDtHoTimeGraph',"Tight unmatched DT")
		self.storeCanvas(canvas,'tight_UnmatchedDtHoTimeGraph_TimeVsEta',marginRight=.15)
		fractionsPerEtaData = self.plotFractionsVsEta(graph,"Fraction of HORecHits at BXID 0,Tight L1",'tight_UnmatchedDtHoTimeGraph_fractionVsEta')
		return canvas,hist,fractionsPerEtaData, boxes, tofFunction
		
	def plotTightHoTimeVsEtaBxWrong(self):
		canvas, hist, graph, boxes, tofFunction = self.makeTimeVsEtaPlot('tight_UnmatchedDtHoBxNot0TimeGraph',"Tight unmatched DT BX Wrong")
		self.storeCanvas(canvas,'tight_UnmatchedDtHoBxNot0TimeGraph_TimeVsEta',marginRight=.15)
		fractionsPerEtaData = self.plotFractionsVsEta(graph,"Fraction of HORecHits at BXID 0,Tight L1 BX wrong",'tight_UnmatchedDtHoBxNot0TimeGraph_fractionVsEta')
		return canvas,hist,fractionsPerEtaData, boxes, tofFunction
	
	def plotHoTimeVsEtaDtRpcTight(self):
		canvas, hist, graph, boxes, tofFunction = self.makeTimeVsEtaPlot('tight_MatchedDtRpcHoTimeGraph')
		self.storeCanvas(canvas,'tight_MatchedDtRpcHoTimeGraph_TimeVsEta',marginRight=.15)
		fractionsPerEtaData = self.plotFractionsVsEta(graph,"Fraction of HORecHits at BXID 0,Tight L1 DT/RPC tight",'tight_MatchedDtRpcHoTimeGraph_fractionVsEta')
		return canvas,hist,fractionsPerEtaData, boxes, tofFunction
	
	def plotL1TimeVsEta(self):
		canvas, hist, label = self.makeL1TimeVsEtaPlot('bxidVsEta')
		self.storeCanvas(canvas,'bxidVsEta')
		return canvas,hist,label
	
	def plotCombined(self):
		return self.makeCombinedEtaPlot()
	
	def plotCombinedTight(self):
		return self.makeCombinedEtaPlot(True)