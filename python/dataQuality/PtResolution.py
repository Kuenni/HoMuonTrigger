from plotting.Utils import getTGraphErrors,getLegend

from ROOT import TCanvas,TF1,TGraphErrors,Double,TMath,TGraph
from plotting.Plot import Plot
import math
from plotting.PlotStyle import setupAxes, colorRwthDarkBlue, colorRwthMagenta,\
	colorRwthGruen, colorRwthTuerkis, colorRwthRot, colorRwthOrange, colorRwthLila
import numpy as np

class PtResolution(Plot):
	
	cached_data = {
		'x':{'values' : [],'errors' : []},
		'L1':{'rms':[],
			'rmsError':[],
			'mean':[],
			'meanError':[],
			'median':[],
			'q16':[],
			'q84':[],
			'q50':[]},
		'L1 Tight':{
			'rms':[],
			'rmsError':[],
			'mean':[],
			'meanError':[],
			'median':[],
			'q16':[],
			'q84':[],
			'q50':[]},
		'L1 And HO':{
			'rms':[],
			'rmsError':[],
			'mean':[],
			'meanError':[],
			'median':[],
			'q16':[],
			'q84':[],
			'q50':[]},
		'L1 Tight And HO':{
			'rms':[],
			'rmsError':[],
			'mean':[],
			'meanError':[],
			'median':[],
			'q16':[],
			'q84':[],
			'q50':[]},
		'L1 !HO':{
			'rms':[],
			'rmsError':[],
			'mean':[],
			'meanError':[],
			'median':[],
			'q16':[],
			'q84':[],
			'q50':[]},
		'L1 Tight !HO':{
			'rms':[],
			'rmsError':[],
			'mean':[],
			'meanError':[],
			'median':[],
			'q16':[],
			'q84':[],
			'q50':[]},
		}
	
	def __init__(self,filename,truth,data = False,debug=False):
		Plot.__init__(self,filename,data=data,debug=debug)
		self.createPlotSubdir('ptResolution')
		self.truthTag = 'patTo' if truth else ''
		self.loadData()
		
	def getMedian(self,h):
	#compute the median for 1-d histogram h1
		nbins = h.GetXaxis().GetNbins()
		xList = []
		yList = []
		for i in range(0,nbins):
			xList.append(h.GetBinCenter(i+1))
			yList.append(h.GetBinContent(i+1))

		return TMath.Median(len(xList),np.array(xList,'d'),np.array(yList,'d'))

	def getQuantiles(self, h ):
		q = np.zeros(5)
		if not h: return q

		probs = np.array([0.025, 0.25, 0.5, .75, 0.975 ],'d')
		h.GetQuantiles( 5, q, probs )
	
		r = np.zeros(5);
		for i in range(0,5): 
			r[i] = q[i]
		return r

	def fillDataIntoCache(self,hist,dataSet):
		self.cached_data[dataSet]['rms'].append(0 if not hist else hist.GetRMS())
		self.cached_data[dataSet]['rmsError'].append(0 if not hist else hist.GetRMSError())
		self.cached_data[dataSet]['mean'].append(0 if not hist else hist.GetMean())
		self.cached_data[dataSet]['meanError'].append(0 if not hist else hist.GetMeanError())
		self.cached_data[dataSet]['median'].append(0 if not hist else self.getMedian(hist))
		quantiles = self.getQuantiles(0 if not hist else hist)
		self.cached_data[dataSet]['q16'].append(quantiles[1])
		self.cached_data[dataSet]['q50'].append(quantiles[2])
		self.cached_data[dataSet]['q84'].append(quantiles[3])
		
	def loadData(self):
		for i in range(0,121):
			#calculate pt range from bin number
			histPt = self.fileHandler.getHistogram('l1PtResolution/' + self.truthTag + 'L1MuonBin%d' % i)
			histPtTight = self.fileHandler.getHistogram('l1PtResolution/' + self.truthTag + 'L1MuonTightBin%d' % i)
			histPtMatch = self.fileHandler.getHistogram('l1PtResolution/' + self.truthTag + 'L1MuonHoMatchBin%d' % i)
			histPtTightMatch = self.fileHandler.getHistogram('l1PtResolution/' + self.truthTag + 'L1MuonTightHoMatchBin%d' % i)
			histNoMatch = self.fileHandler.getHistogram('l1PtResolution/' + self.truthTag + 'L1MuonNotHoMatchBin%d' % i)
			histTightNoMatch = self.fileHandler.getHistogram('l1PtResolution/' + self.truthTag + 'L1MuonTightNotHoMatchBin%d' % i)
			if i < 40:
				self.cached_data['x']['values'].append(i + 0.5)
				self.cached_data['x']['errors'].append(0.5)
			else:
				self.cached_data['x']['values'].append(i*2 - 40 + 1)
				self.cached_data['x']['errors'].append(1)
			
			self.fillDataIntoCache(histPt,'L1')
			self.fillDataIntoCache(histPtMatch,'L1 And HO')
			self.fillDataIntoCache(histPtTight,'L1 Tight')
			self.fillDataIntoCache(histPtTightMatch,'L1 Tight And HO')
			self.fillDataIntoCache(histNoMatch,'L1 !HO')
			self.fillDataIntoCache(histTightNoMatch,'L1 Tight !HO')

			self.printProgress(i+1, 121,1)
	
	###
	#	Get the intrgrals of the source histograms
	###
	def getHistogramIntegralsAsList(self,sourceName):
		values 		= []
		valuesErr 	= []
		
		for i in range(0,121):
			hist = self.fileHandler.getHistogram(sourceName + 'Bin%d' % i)
			if hist != None:
				integral = hist.Integral()
				values.append(integral)
				valuesErr.append(math.sqrt(integral))
			else:
				values.append(0)
				valuesErr.append(0)
				self.printProgress(i+1, 121,1)
		return values,valuesErr

	
	def plotL1(self):
		return self.makeQuantilePlot('L1')
	
	def plotL1Tight(self):
		return self.makeQuantilePlot('L1 Tight')
	
	def plotL1AndHo(self):
		return self.makeQuantilePlot('L1 And HO')
	
	def plotL1TightAndHo(self):
		return self.makeQuantilePlot('L1 Tight And HO')
	
	def plotL1NotHo(self):
		return self.makeQuantilePlot('L1 !HO')
	
	def plotL1TightNotHo(self):
		return self.makeQuantilePlot('L1 Tight !HO')
	
	def plotPtResolutionHistograms(self):

		c = TCanvas()
		graphL1 = getTGraphErrors(self.cached_data['x']['values'],
								self.cached_data['L1']['rms'],
								ey = self.cached_data['L1']['rmsError'],
								ex=self.cached_data['x']['errors'])
		graphL1.SetMarkerStyle(20)
		graphL1.SetMarkerColor(colorRwthDarkBlue)
		graphL1.SetLineColor(colorRwthDarkBlue)
		graphL1.SetTitle("RMS of L1 Objects;p_{T} / GeV;RMS / GeV")
		graphL1.GetYaxis().SetRangeUser(0,75)
		graphL1.Draw('ap')
		
		graphL1Tight = getTGraphErrors(self.cached_data['x']['values'],
								self.cached_data['L1 Tight']['rms'],
								ey = self.cached_data['L1 Tight']['rmsError'],
								ex=self.cached_data['x']['errors'])#rmsL1Tight
		graphL1Tight.SetMarkerStyle(21)
		graphL1Tight.SetMarkerColor(colorRwthGruen)
		graphL1Tight.SetLineColor(colorRwthGruen)
		graphL1Tight.Draw('samep')
		
		graphL1AndHo = getTGraphErrors(self.cached_data['x']['values'],
								self.cached_data['L1 And HO']['rms'],
								ey = self.cached_data['L1 And HO']['rmsError'],
								ex=self.cached_data['x']['errors'])#rmsL1AndHo
		graphL1AndHo.SetMarkerStyle(26)
		graphL1AndHo.SetMarkerColor(colorRwthMagenta)
		graphL1AndHo.SetLineColor(colorRwthMagenta)
		graphL1AndHo.Draw('samep')
		
		graphL1TightAndHo = getTGraphErrors(self.cached_data['x']['values'],
								self.cached_data['L1 Tight And HO']['rms'],
								ey = self.cached_data['L1 Tight And HO']['rmsError'],
								ex=self.cached_data['x']['errors'])#rmsL1TightAndHo
		graphL1TightAndHo.SetMarkerStyle(27)
		graphL1TightAndHo.SetMarkerColor(colorRwthRot)
		graphL1TightAndHo.SetLineColor(colorRwthRot)
		graphL1TightAndHo.Draw('samep')
		
		graphL1NotHo = getTGraphErrors(self.cached_data['x']['values'],
								self.cached_data['L1 !HO']['rms'],
								ey = self.cached_data['L1 !HO']['rmsError'],
								ex=self.cached_data['x']['errors'])#rmsL1NotHo
		graphL1NotHo.SetMarkerStyle(29)
		graphL1NotHo.SetMarkerColor(colorRwthOrange)
		graphL1NotHo.SetLineColor(colorRwthOrange)
		graphL1NotHo.Draw('samep')
		
		graphL1TightNotHo = getTGraphErrors(self.cached_data['x']['values'],
								self.cached_data['L1 Tight !HO']['rms'],
								ey = self.cached_data['L1 Tight !HO']['rmsError'],
								ex=self.cached_data['x']['errors'])#rmsL1TightNotHo
		graphL1TightNotHo.SetMarkerStyle(34)
		graphL1TightNotHo.SetMarkerColor(colorRwthLila)
		graphL1TightNotHo.SetLineColor(colorRwthLila)
		graphL1TightNotHo.Draw('samep')
		
		setupAxes(graphL1)
		
		legend = getLegend(y2 = .9)
		legend.AddEntry(graphL1,'RMS L1','ep')
		legend.AddEntry(graphL1AndHo,'RMS L1 and HO','ep')
		legend.AddEntry(graphL1Tight,'RMS L1 Tight','ep')
		legend.AddEntry(graphL1TightAndHo,'RMS L1 Tight and HO','ep')
		legend.AddEntry(graphL1NotHo,'RMS L1 & !HO','ep')
		legend.AddEntry(graphL1TightNotHo,'RMS L1 Tight & !HO','ep')
		legend.Draw()
		
		label = self.drawLabel()
		
		c.Update()
		
		self.storeCanvas(c, 'rmsVsPt' + self.truthTag)
		#c2 = TCanvas('cfitResults','fitResults',800,0,800,600)
		#graphL1Fit.Draw('AP')
				
		return graphL1,graphL1AndHo,legend,label,graphL1TightAndHo,graphL1Tight, graphL1NotHo, graphL1TightNotHo#,c2,c

	def makeQuantilePlot(self,dataSet):
		c = TCanvas(dataSet + 'quantiles')

		graphMean = getTGraphErrors(self.cached_data['x']['values'],
								self.cached_data[dataSet]['mean'],
								ey = self.cached_data[dataSet]['meanError'],
								ex=self.cached_data['x']['errors'])
		setupAxes(graphMean)
		graphMean.SetMarkerStyle(20)
		graphMean.SetMarkerColor(colorRwthDarkBlue)
		graphMean.SetLineColor(colorRwthDarkBlue)
		
		graphMedian = getTGraphErrors(self.cached_data['x']['values'],
								self.cached_data[dataSet]['median'])
		graphMedian.SetMarkerStyle(21)
		graphMedian.SetMarkerColor(colorRwthMagenta)
		graphMedian.SetLineColor(colorRwthMagenta)
		graphMedian.SetFillColor(colorRwthMagenta)
		
		graphQ50 = getTGraphErrors(self.cached_data['x']['values'],
								self.cached_data[dataSet]['q50'])
		graphQ50.SetMarkerStyle(22)
		graphQ50.SetMarkerColor(colorRwthLila)
		graphQ50.SetLineColor(colorRwthLila)
		
		graphShape68 = TGraph()
		graphShape68.SetName(dataSet)
		graphShape68.SetFillStyle(1001)
		graphShape68.SetFillColorAlpha(colorRwthGruen,.75)
		graphShape68.SetLineColorAlpha(colorRwthGruen,.75)

		for i in zip(self.cached_data['x']['values'],self.cached_data[dataSet]['q16']):
			graphShape68.SetPoint(graphShape68.GetN(),i[0],i[1])
		
		for i in zip(reversed(self.cached_data['x']['values']),reversed(self.cached_data[dataSet]['q84'])):
			graphShape68.SetPoint(graphShape68.GetN(),i[0],i[1])
		
		graphShape68.GetYaxis().SetRangeUser(0,200)
		graphShape68.SetTitle(";p_{T,Reco} / GeV;p_{T,L1} / GeV")#"Mean, Median, and Quantiles of " + dataSet + " p_{T};p_{T,Reco} / GeV;p_{T,L1} / GeV")
		setupAxes(graphShape68)

		graphShape68.Draw('a f')
		graphMean.Draw('same p')
		graphMedian.Draw('same,l')
		graphQ50.Draw('same,l')

		
		legend = getLegend(y1=.75,y2 = .95,x1=.1,x2=.5)
		legend.AddEntry(graphMean,dataSet + ' p_{T} Mean','ep')
		legend.AddEntry(graphMedian,dataSet + ' p_{T} Median','l')
		legend.AddEntry(graphShape68,dataSet + ' p_{T} Q_{25} - Q_{75}','f')
		legend.AddEntry(graphQ50,dataSet + ' p_{T} Q_{50}','l')
		legend.Draw()
				
		c.Update()
		self.storeCanvas(c,'quantiles_' + self.truthTag + '_' + dataSet.replace(' ','_').replace('!','Not'), marginRight = .02)
		return c,graphMean,legend,graphMedian,graphShape68,graphQ50

	def plotTightPtResolution(self):
		c = TCanvas('cTightResolution','cTightResolution')
		xData = [self.cached_data['x']['values'],self.cached_data['x']['errors']]
		tight = [self.cached_data['L1 Tight']['rms'],
				self.cached_data['L1 Tight']['rmsError']]
		tightAndHo = [self.cached_data['L1 Tight And HO']['rms'],
				self.cached_data['L1 Tight And HO']['rmsError']]
		tightAndNotHo = [self.cached_data['L1 Tight !HO']['rms'],
				self.cached_data['L1 Tight !HO']['rmsError']]
		
		graphL1Tight = getTGraphErrors(xData[0], tight[0], ex=xData[1], ey=tight[1])
		graphL1TightHo = getTGraphErrors(xData[0], tightAndHo[0], ex=xData[1], ey=tightAndHo[1])
		graphL1TightNotHo = getTGraphErrors(xData[0], tightAndNotHo[0], ex=xData[1], ey=tightAndNotHo[1])
		
		graphL1Tight.SetMarkerStyle(25)
		graphL1Tight.SetMarkerColor(colorRwthGruen)
		graphL1Tight.SetLineColor(colorRwthGruen)
		
		graphL1TightHo.SetMarkerStyle(20)
		graphL1TightHo.SetTitle('RMS of tight L1 Objects;p_{T,RECO} / GeV;L1 p_{T} RMS / GeV')
		graphL1TightHo.SetMarkerColor(colorRwthDarkBlue)
		graphL1TightHo.SetLineColor(colorRwthDarkBlue)
		
		graphL1TightNotHo.SetMarkerStyle(34)
		graphL1TightNotHo.SetMarkerColor(colorRwthMagenta)
		graphL1TightNotHo.SetLineColor(colorRwthMagenta)
		
		graphL1TightHo.Draw('ap')
		graphL1TightNotHo.Draw('samep')
		graphL1Tight.Draw('samep')
		
		label = self.drawLabel()
		
		setupAxes(graphL1TightHo)
		legend = getLegend(y2 = .9)
		legend.AddEntry(graphL1Tight,'L1 Tight','ep')
		legend.AddEntry(graphL1TightHo,'L1 Tight & HO','ep')
		legend.AddEntry(graphL1TightNotHo,'L1 Tight & !HO','ep')
		legend.Draw()
		c.Update()
		
		self.storeCanvas(c, 'rmsVsPt_tight' + self.truthTag)
		
		c2 = self.makeIntegralPlot('patToL1MuonTight')
		
		return  c,graphL1TightHo,graphL1TightNotHo,legend,graphL1Tight,c2,label
		
		##
		# Plot the integral for each histogram as a control plot
		##
	def makeIntegralPlot(self,datasetName):
		cControlPlot = TCanvas('cControlPlots','control plot integral')
		xData = [self.cached_data['x']['values'],self.cached_data['x']['errors']]
		dataAllL1,errorAllL1 = self.getHistogramIntegralsAsList('l1PtResolution/' + datasetName)
		dataMatch,errorMatch = self.getHistogramIntegralsAsList('l1PtResolution/' + datasetName + "HoMatch")
		dataNoMatch,errorNoMatch = self.getHistogramIntegralsAsList('l1PtResolution/' + datasetName + "NotHoMatch")
		
		graphIntegralsAllL1 = getTGraphErrors(xData[0], dataAllL1, ex=xData[1], ey=errorAllL1)
		graphIntegralsMatch = getTGraphErrors(xData[0], dataMatch, ex=xData[1], ey=errorMatch)
		graphIntegralsNoMatch = getTGraphErrors(xData[0], dataNoMatch, ex=xData[1], ey=errorNoMatch)
		
		graphIntegralsAllL1.SetMarkerStyle(21)
		graphIntegralsAllL1.SetTitle('Integrals of ' + datasetName + ';p_{T} / GeV;# entries')
		graphIntegralsAllL1.SetMarkerColor(colorRwthDarkBlue)
		graphIntegralsAllL1.SetLineColor(colorRwthDarkBlue)
		
		graphIntegralsMatch.SetMarkerStyle(20)
		graphIntegralsMatch.SetMarkerColor(colorRwthMagenta)
		graphIntegralsMatch.SetLineColor(colorRwthMagenta)
		
		graphIntegralsNoMatch.SetMarkerStyle(34)
		graphIntegralsNoMatch.SetMarkerColor(colorRwthTuerkis)
		graphIntegralsNoMatch.SetLineColor(colorRwthTuerkis)
		
		graphIntegralsAllL1.Draw('ap')
		graphIntegralsMatch.Draw('samep')
		graphIntegralsNoMatch.Draw('samep')
		
		label = self.drawLabel()
		
		setupAxes(graphIntegralsAllL1)
		legend2 = getLegend(y2 = .9)
		legend2.AddEntry(graphIntegralsAllL1,datasetName,'ep')
		legend2.AddEntry(graphIntegralsMatch,datasetName + " And HO",'ep')
		legend2.AddEntry(graphIntegralsNoMatch,datasetName + " !HO",'ep')
		legend2.Draw()
		cControlPlot.Update()
		
		self.storeCanvas(cControlPlot, 'rmsVsPt_integrals' + datasetName)
		return label, legend2, graphIntegralsMatch, graphIntegralsNoMatch, cControlPlot,graphIntegralsAllL1
	
	def plotLoosePtResolution(self):
		c = TCanvas('cLooseResolution','cLooseResolution')

		xData = [self.cached_data['x']['values'],self.cached_data['x']['errors']]
		tight = [self.cached_data['L1']['rms'],
				self.cached_data['L1']['rmsError']]
		tightAndHo = [self.cached_data['L1 And HO']['rms'],
				self.cached_data['L1 And HO']['rmsError']]
		tightAndNotHo = [self.cached_data['L1 !HO']['rms'],
				self.cached_data['L1 !HO']['rmsError']]
		
		graphL1Tight = getTGraphErrors(xData[0], tight[0], ex=xData[1], ey=tight[1])
		graphL1TightHo = getTGraphErrors(xData[0], tightAndHo[0], ex=xData[1], ey=tightAndHo[1])
		graphL1TightNotHo = getTGraphErrors(xData[0], tightAndNotHo[0], ex=xData[1], ey=tightAndNotHo[1])
		
		graphL1TightHo.SetMarkerStyle(20)
		graphL1TightHo.SetTitle('RMS of loose L1 Objects;p_{T,RECO} / GeV;L1 p_{T} RMS / GeV')
		graphL1TightHo.SetMarkerColor(colorRwthDarkBlue)
		graphL1TightHo.SetLineColor(colorRwthDarkBlue)
		
		graphL1TightNotHo.SetMarkerStyle(34)
		graphL1TightNotHo.SetMarkerColor(colorRwthMagenta)
		graphL1TightNotHo.SetLineColor(colorRwthMagenta)
		
		graphL1Tight.SetMarkerStyle(25)
		graphL1Tight.SetMarkerColor(colorRwthGruen)
		graphL1Tight.SetLineColor(colorRwthGruen)
		
		graphL1TightHo.Draw('ap')
		graphL1TightNotHo.Draw('samep')
		graphL1Tight.Draw('samep')
		
		label = self.drawLabel()
		
		setupAxes(graphL1TightHo)
		legend = getLegend(y2 = .9)
		legend.AddEntry(graphL1Tight,'L1','ep')
		legend.AddEntry(graphL1TightHo,'L1 & HO','ep')
		legend.AddEntry(graphL1TightNotHo,'L1 & !HO','ep')
		legend.Draw()
		c.Update()

		self.storeCanvas(c, 'rmsVsPt_loose' + self.truthTag)
	
		c2 = self.makeIntegralPlot('patToL1Muon')
		return c,graphL1TightHo,graphL1TightNotHo,legend,label,graphL1Tight,c2
