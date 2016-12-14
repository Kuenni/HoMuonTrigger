from ROOT import gROOT,gStyle, TColor, TPaveText, TGraph, ROOT, Double, TBox, TH2D, TH1D, TEfficiency, TGraphErrors, THStack
from math import pi
import matplotlib.pyplot as plt
import sys
from array import array

colorRwthMagenta 	= TColor.GetColor("#E30066")
colorRwthLightBlue 	= TColor.GetColor("#8EBAE5")
colorRwthDarkBlue 	= TColor.GetColor("#00549F")
colorRwthTuerkis	= TColor.GetColor("#0098A1")
colorRwthGruen		= TColor.GetColor("#57AB27")
colorRwthRot		= TColor.GetColor("#CC071E")
colorRwthOrange 	= TColor.GetColor("#F6A800")
colorRwthViolett	= TColor.GetColor("#612158")
colorRwthLila		= TColor.GetColor("#7A6FAC")
colorRwthLightGray	= TColor.GetColor("#CFD1D2")
colorRwthDarkGray	= TColor.GetColor("#646567")

def setColorPalette():
	nGradients = 100
	palette = []
	r = [0, 0/255., 227./255., 255/255., 244/255.]
	g = [0,0/255.,0/255.,237/255.,244/255.]
	b = [0,159/255.,102/255.,0/255.,244/255.]
	stop = [0,.25,.55,.9,1]
	FI = TColor.CreateGradientColorTable(5,array('d',stop),array('d',r),array('d',g),array('d',b),nGradients)
	for i in range(0,nGradients):
		palette.append(FI+i)
	gStyle.SetPalette(nGradients,array('i',palette))
	
def getWaterMark(markPosition):
	watermark = TPaveText(markPosition['x1ndc'], markPosition['y1ndc'], markPosition['x2ndc'], markPosition['y2ndc'],"NDC")
	markText = watermark.AddText("A. Kuensken")
	markText.SetTextAlign(12)
	markText.SetTextSize(0.0316456)
	markText.SetTextColor(15)
	watermark.SetBorderSize(0)
	watermark.SetFillStyle(0)
	return watermark

def drawWaterMark(markPosition = {'x1ndc' : 0.093, 'y1ndc' : 0.898, 'x2ndc' : 0.319, 'y2ndc' : 0.940}):
	watermark = getWaterMark(markPosition = markPosition)
	watermark.Draw()
	return watermark

def getLabelCmsPrivateSimulation( x1ndc = 0.6, y1ndc = 0.90, x2ndc = 0.9, y2ndc = 0.93):
	labelCmsPrivateSimulation = TPaveText(x1ndc,y1ndc,x2ndc,y2ndc,"NDC")
	labelCmsPrivateSimulation.AddText("#font[62]{CMS private work}, #font[72]{simulation}")
	labelCmsPrivateSimulation.SetBorderSize(1)
	labelCmsPrivateSimulation.SetFillColor(0) # 0 == White
	return labelCmsPrivateSimulation

def drawLabelCmsPrivateSimulation( x1ndc = 0.6, y1ndc = 0.90, x2ndc = 0.9, y2ndc = 0.93):
	label = getLabelCmsPrivateSimulation(x1ndc,y1ndc,x2ndc,y2ndc)
	label.Draw()
	return label

def getLabelCmsPrivateData( x1ndc = 0.6, y1ndc = 0.90, x2ndc = 0.9, y2ndc = 0.93):
	labelCmsPrivateSimulation = TPaveText(x1ndc,y1ndc,x2ndc,y2ndc,"NDC")
	labelCmsPrivateSimulation.AddText("#font[62]{CMS private work}, #font[72]{2015 data}")
	labelCmsPrivateSimulation.SetBorderSize(1)
	labelCmsPrivateSimulation.SetFillColor(0) # 0 == White
	return labelCmsPrivateSimulation

def drawLabelCmsPrivateData( x1ndc = 0.6, y1ndc = 0.90, x2ndc = 0.9, y2ndc = 0.93):
	label = getLabelCmsPrivateData(x1ndc,y1ndc,x2ndc,y2ndc)
	label.Draw()
	return label

chimney1 = TGraph()
chimney1.SetPoint(0,.3425,1.48)
chimney1.SetPoint(1,.3425,1.65)
chimney1.SetPoint(2,.435,1.65)
chimney1.SetPoint(3,.435,1.48)
chimney1.SetPoint(4,.3425,1.48)
chimney1.SetLineColor(colorRwthMagenta)
chimney1.SetLineWidth(2)

chimney2 = TGraph()
chimney2.SetPoint(0,-.3425,0.87)
chimney2.SetPoint(1,-.3425,1.22)
chimney2.SetPoint(2,-.435,1.22)
chimney2.SetPoint(3,-.435,0.87)
chimney2.SetPoint(4,-.3425,0.87)
chimney2.SetLineColor(colorRwthMagenta)
chimney2.SetLineWidth(2)

#Converts a graph in eta and phi to hcal coordinates with ieta and iphi
def convertToHcalCoords(inputGraph):
	gConvertedToiEtaiPhi = TGraph()
	gConvertedToiEtaiPhi.SetMarkerStyle(6)
	for i in range(0, inputGraph.GetN()):
		x = Double(0)
		y = Double(0)
		inputGraph.GetPoint(i,x,y)
		if( y < 0 ):
			y += 2*pi
		if( y > 2*pi ):
			y -= 2*pi
		xIEta = x/0.087 + ( 0.5 if x > 0 else -0.5 )
		yIPhi = y/0.087 + 1
		gConvertedToiEtaiPhi.SetPoint(i,xIEta,yIPhi)
	gConvertedToiEtaiPhi.GetXaxis().SetTitle("i#eta")
	gConvertedToiEtaiPhi.GetYaxis().SetTitle("i#phi")
	return gConvertedToiEtaiPhi

def drawHcalBoxesHcalCoords(canvas):
	canvas.cd()
	boxes = []
	for i in range(-15,15):
		for j in range(0,72):
			box = TBox(i,j,(i+1),(j+1))
			box.SetFillStyle(0)
			box.SetLineColor(colorRwthMagenta)
			box.SetLineWidth(2)
			box.Draw()
			boxes.append(box)
	return boxes

def drawHoBoxes(canvas):
	canvas.cd()
	boxes = []
	for i in range(-1,2):
		for j in range(-1,2):
			box = TBox(i*0.087 - 0.0435,j*0.087 - 0.0435,(i*0.087+0.0435),(j*0.087+0.0435))
			box.SetFillStyle(0)
			box.SetLineColor(colorRwthMagenta)
			box.SetLineWidth(3)
			box.Draw()
			boxes.append(box)
	return boxes

def setPlotStyle():
	gStyle.SetPadGridX(1)
	gStyle.SetPadGridY(1)
	gStyle.SetOptStat(1110)
	gStyle.SetLineWidth(2)
	gStyle.SetHistLineWidth(3)
	gStyle.SetPadTickX(1)
	gStyle.SetPadTickY(1)
	gStyle.SetTextFont(62)
	gStyle.SetLabelFont(62)
	gStyle.SetTitleFont(62)
	gStyle.SetLegendBorderSize(1)
	gStyle.SetLegendFont(62)
#	setColorPalette()

#Set all Axes to bold font
def setupAxes(plot):
	
	if plot.__class__ == TEfficiency:
		plot = plot.GetPaintedGraph()
	plot.GetXaxis().SetTitleFont(62)
	plot.GetYaxis().SetTitleFont(62)
	plot.GetXaxis().SetLabelFont(62)
	plot.GetYaxis().SetLabelFont(62)
	# Check for the function. Otherwise it crashes with TGraphs
	if hasattr(plot, 'GetZaxis'):
		plot.GetZaxis().SetTitleFont(62)
		plot.GetZaxis().SetLabelFont(62)
	#Do not continue if you're working with TGraphs
	if (plot.__class__.__name__.find('TGraph') != -1):
		return
	
	if(plot.__class__ == THStack):
		plot = plot.GetHistogram()
	plot.SetTitleFont(62)
	plot.SetTitleSize(0.06)
	setBigAxisTitles(plot)
	setupPalette(plot)

def setBigAxisTitles(plot):
	if plot.__class__ == TEfficiency:
		plot = plot.GetPaintedGraph()
	plot.GetXaxis().SetTitleSize(0.06)
	plot.GetXaxis().SetTitleOffset(0.7)
	plot.GetYaxis().SetTitleSize(0.06)
	plot.GetYaxis().SetTitleOffset(0.8)
	plot.GetZaxis().SetTitleSize(0.06)
	plot.GetZaxis().SetTitleOffset(0.6)
	
#Set the stat box display Options
def setStatBoxOptions(plot,option):
	stats = plot.GetListOfFunctions().FindObject("stats")
	stats.SetOptStat(option)

#Set the stat box position
def setStatBoxPosition(plot,x1 = 0.7, x2 = 0.9, y1 = 0.8, y2 = 0.95):
	stats = plot.GetListOfFunctions().FindObject("stats")
	stats.SetX1NDC(x1)
	stats.SetX2NDC(x2)
	stats.SetY1NDC(y1)
	stats.SetY2NDC(y2)

'''
The x2ndc is for now kept for downwards compatibility
'''
def setupPalette(plot,x2ndc=.92,shiftBy=0.0):
	#Set as many color palette divisions as possible
	plot.SetContour(99)
	#make the palette as small as possible
	pal = plot.GetListOfFunctions().FindObject("palette")
	if not pal:
		return
	pal.SetX2NDC(x2ndc - shiftBy)
	pal.SetX1NDC(0.905 - shiftBy)
	pal.SetY2NDC(.95)
	plot.GetZaxis().SetTitleOffset(1.)

#Function that returns a new TH2D with the axes already set up
def getTH2D(name,title,nBinsX,xLow,xHigh,nBinsY,yLow,yHigh):
	hist = TH2D(name,title,nBinsX,xLow,xHigh,nBinsY,yLow,yHigh)
	setupAxes(hist)
	return hist

#Function that returns a new TH2D with the axes already set up
def getTH1D(name,title,nBinsX,xLow,xHigh):
	hist = TH1D(name,title,nBinsX,xLow,xHigh)
	setupAxes(hist)
	return hist

#Output function for the progress in a python script
def printProgress(done,total):
	s = getProgressString(done, total)
	sys.stdout.write(s)
	sys.stdout.flush()
	pass

#Output function for the progress in a python script
def getProgressString(done,total):
	nHashes = int(done/float(total)*80)
	progressbar = '\r[%s%s] %5.2f%% done.' % (nHashes*'#',(80-nHashes)*' ',done*100/float(total))
	return progressbar

def pyplotCmsPrivateLabel(ax,x=0.995,y=0.945):
	plt.text(x, y, r'$\mathbf{CMS}$ private, $\mathit{simulation}$',
        horizontalalignment='right',
        verticalalignment='bottom',
        transform=ax.transAxes,
        bbox=dict(facecolor='white'))
