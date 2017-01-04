from plotting import OutputModule
from plotting.PlotStyle import setupAxes

from ROOT import Double,TGraphErrors,TLegend,vector,TMath,gPad,TPad
from array import array

import math
from numpy import nan
import numpy as np
from math import sqrt

commandLine = OutputModule.CommandLineHandler('[Utils.py] ')

L1_PHI_BIN = math.pi/72.
L1_ETA_BIN = 0.1

def average2DHistogramBinwise(histWeights,histCounter):
	for i in range(0,histWeights.GetNbinsX()):
		for j in range(0,histWeights.GetNbinsY()):
			if histCounter.GetBinContent(histCounter.GetBin(i,j)) != 0:
				histWeights.SetBinContent(histWeights.GetBin(i,j),histWeights.GetBinContent(histWeights.GetBin(i,j))
										/histCounter.GetBinContent(histCounter.GetBin(i,j)))
	return histWeights

#Set axis range and labels for the 2D histograms showing E Average around L1 direction
def setupEAvplot(histE,histC = None,xmin = -0.4, xmax = 0.4, ymin = -0.4, ymax = 0.4,same = False, limitForAll = None):
	if histC != None:
		histE = average2DHistogramBinwise(histE,histC)
	if same:
		if limitForAll == None:
			commandLine.warning('WARNING: Requested same histogram borders for all ranges but '
							'did not give limitForAll parameter. Using default values instead!')
		else:
			xmin = ymin = -limitForAll
			xmax = ymax = limitForAll
	histE.GetXaxis().SetRangeUser(xmin,xmax)
	histE.GetYaxis().SetRangeUser(ymin,ymax)
	histE.SetStats(0)
	histE.GetXaxis().SetTitle('#Delta#eta')
	histE.GetYaxis().SetTitle('#Delta#phi')
	histE.GetZaxis().SetTitle('Reconstructed Energy / GeV')
	setupAxes(histE)
	return histE

def normalizeTH2DAlongXAxis(hist):
	#First count all entries along axis
	for yBin in range(hist.GetYaxis().GetNbins()+1):
		nEntries = 0
		#First count all entries in x-direction
		for xBin in range(hist.GetXaxis().GetNbins()):
			nEntries += hist.GetBinContent(xBin,yBin)
		#Then normalize
		for xBin in range(hist.GetXaxis().GetNbins()):
			if nEntries:
				hist.SetBinContent(xBin,yBin,hist.GetBinContent(xBin,yBin)/nEntries)

def fillGraphIn2DHist(graph,hist):
	x = Double(0)
	y = Double(0)
	commandLine.output('Filling graph in 2D histogram:')
	if not graph:
		commandLine.error('Graph is null!')
	nTotal = graph.GetN()
	for i in range(0,nTotal):
		graph.GetPoint(i,x,y)
		hist.Fill(x,y)
		if(not i%10000):
			commandLine.printProgress(i,nTotal)
		if(i == nTotal - 1):
			commandLine.printProgress(nTotal, nTotal)
	return hist

def fillGraphIn1DHist(graph,hist,axis = 'y'):
	x = Double(0)
	y = Double(0)
	if( not axis == 'y' or not axis == 'x'):
		commandLine.warning('Unsupported axis: %s' % axis)
		axis = y
	commandLine.output('Filling %s axis of graph in 1D histogram:' % axis)
	nTotal = graph.GetN()
	for i in range(0,nTotal):
		graph.GetPoint(i,x,y)
		if axis == 'x':
			hist.Fill(x)
		else:
			hist.Fill(y)
		if(not i%10000):
			commandLine.printProgress(i,nTotal)
		if(i == nTotal - 1):
			commandLine.printProgress(nTotal, nTotal)
	return hist

#
# Fill a slice of a 2D ROOT graph in a 2D ROOT histogram
#
def fill2DGraphIn2DHist(graph,hist,sliceOnValue):
	x = Double(0)
	y = Double(0)
	z = Double(0)
	
	# sliceOnValue can either be a value or a range 
	minSlice = maxSlice = sliceOnValue
	
	if sliceOnValue is list:
		minSlice = sliceOnValue[0]
		maxSlice = sliceOnValue[1]
	
	commandLine.output('Filling 2D graph in 2D histogram:')
	nTotal = graph.GetN()
	for i in range(0,nTotal):
		graph.GetPoint(i,x,y,z)
		if(x >= minSlice or x <= maxSlice):
			hist.Fill(y,z)
		if(not i%10000):
			commandLine.printProgress(i,nTotal)
		if(i == nTotal - 1):
			commandLine.printProgress(nTotal, nTotal)
	return hist

#Returns a 2D hisotgram containing the binwise difference of both objects
#Python histogram means the entries matrix coming from numpy histogram 2d 
def comparePythonAndRoot2DHist(pythonHist, rootHist):
	contourLevels = [-1,-.2,-.1,.1,.2,1]
	
	if not (len(pythonHist) == rootHist.GetNbinsY() and len(pythonHist[0]) == rootHist.GetNbinsX()):
		commandLine.output('Error! Cannot compare python and root histogram with different number of bins!')
		return
	
	comparisonHist = rootHist.Clone('comparison' + rootHist.GetName())
	comparisonHist.Reset()
	comparisonHist.SetStats(0)
	comparisonHist.SetContour(len(contourLevels),array('d',contourLevels))
	for j,y in enumerate(reversed(pythonHist)):
		for i,x in enumerate(y):
			comparisonHist.SetBinContent(i+1,len(pythonHist)-j, x - rootHist.GetBinContent( i+1 , len(pythonHist) - j))
			#sys.stdout.write( str(x - histNormalBins.GetBinContent( i+1 , len(pythonHist) - j)) + '\t' )
			pass
	return comparisonHist

#Returns a 2D hisotgram containing the binwise difference of both objects
def compareTwoRoot2DHists(rootHist1, rootHist2):
	contourLevels = [-3.4e48,-.2,-.1,.1,.2,1]
	if not (rootHist1.GetNbinsX() == rootHist2.GetNbinsX() and rootHist1.GetNbinsY() == rootHist2.GetNbinsY()):
		commandLine.output('Error! Cannot compare two root histograms with different number of bins!')
		return
	comparisonHist = rootHist1.Clone('comparison' + rootHist1.GetName())
	comparisonHist.Reset()
	comparisonHist.SetStats(0)
	comparisonHist.SetContour(len(contourLevels),array('d',contourLevels))
	for x in range(0,rootHist1.GetNbinsX()):
		for y in range(0,rootHist1.GetNbinsY()):
			comparisonHist.SetBinContent(x,y,rootHist1.GetBinContent(x,y) - rootHist2.GetBinContent(x,y))
			pass
	return comparisonHist

def extractTEfficiencyToList(tEffObject):
	xVals = []
	yVals = []
	yErrLow = []
	yErrUp = []
	for i in range(tEffObject.GetPassedHistogram().GetNbinsX()):
		if tEffObject.GetTotalHistogram().GetBinContent(i) != 0:
			yVals.append(tEffObject.GetPassedHistogram().GetBinContent(i)/tEffObject.GetTotalHistogram().GetBinContent(i)*100)
			xVals.append(tEffObject.GetTotalHistogram().GetBinCenter(i))
			yErrLow.append(tEffObject.GetEfficiencyErrorLow(i))
			yErrUp.append(tEffObject.GetEfficiencyErrorUp(i))
	return xVals, yVals,yErrLow,yErrUp
	
def getTGraphErrors(x,y,ex = None,ey = None):
	if (ex == None):
		ex = [0]*len(x)
	if (ey == None):
		ey = [0]*len(x)
	return TGraphErrors(len(x),array('f',x),array('f',y),array('f',ex),array('f',ey))
	
def getLegend(x1=.6,y1=.8,x2=.9,y2=.85):
	l =  TLegend(x1,y1,x2,y2)
#	l.SetTextFont(62)
	return l

def makeResidualsPad(pad):
	pad.Divide(1,2)
	pad.cd(1).SetBottomMargin(0)
	pad.cd(1).SetPad(0,0.3,1,1)
	pad.cd(2).SetTopMargin(0)
	pad.cd(2).SetBottomMargin(0.15)
	pad.cd(2).SetPad(0,0,1,0.3)
	pad.cd(1)
	return pad
	
def calcPercent(numerator, denominator):
	if(denominator == 0):
		commandLine.error('Tried to divide by 0')
		return nan
	return numerator/float(denominator)*100

def calcSigma(num,denom):
	if(denom == 0):
		commandLine.error('Tried to divide by 0')
		return nan
	return sqrt(num/float(denom*denom) + num*num/float(pow(denom, 3)))

def getMedian(h):
	#compute the median for 1-d histogram h1
		nbins = h.GetXaxis().GetNbins()
		xList = []
		yList = []
		for i in range(0,nbins):
			xList.append(h.GetBinCenter(i+1))
			yList.append(h.GetBinContent(i+1))
		# For update to ROOT 6, explicit type is needed. Seems to be in connection with
		# only templated method available in ROOT 6. Write:
		#
		# return TMath.Median('double')(len(xList),np.array(xList,'d'),np.array(yList,'d'))
		#
		return TMath.Median(Double)(len(xList),np.array(xList,'d'),np.array(yList,'d'))

def getXinNDC(x):
	gPad.Update()
	return (x - gPad.GetX1())/(gPad.GetX2()-gPad.GetX1())

def phiWrapCheck(phi2,phi1):
	delta_phi = phi2 - phi1;
	if(delta_phi < -math.pi):
		return (2*math.pi + delta_phi)
	if(delta_phi > math.pi):
		return (delta_phi - 2*math.pi)
	return delta_phi

