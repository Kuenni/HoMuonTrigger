import os
from plotting.Plot import Plot
from plotting.OutputModule import CliColors
from plotting.PlotStyle import setupPalette, colorRwthDarkBlue, colorRwthMagenta,\
	setupAxes
from ROOT import TCanvas,TBox,THStack,TLegend,TPaveText

class QualityCode(Plot):
	def __init__(self,filename,data=False, debug = False):
		Plot.__init__(self,filename,data,debug)
		self.createPlotSubdir('qualityCodes')

		self.qualityCodeDict = {
				1:'Halo Muon',
				2:'Very low qual. Type 1',
				3:'Very low qual. Type 2',
				4:'Very low qual. Type 3',
				5:'unmatched RPC',
				6:'unmatched DT or CSC',
				7:'matched DT-RPC or CSC-RPC'
				}

		self.gridSizeDict = {
			0:'Central',
			1:'3x3',
			2:'5x5'}

		
	def plot3x3MatchQualityCodesVsPt(self):
		c = TCanvas('cMatchQC3x3','Match QC 3x3',0,0,900,700)
		c.SetLogz()
		hist = self.fileHandler.getHistogram('qualityCode/L1Muon3x3Match_QcVsPt')
		hist.SetStats(0)
		hist.Scale(1,'width')
		hist.Draw('colz')
		c.Update()
		setupPalette(hist)
		c.Update()
				
		label = self.drawLabel()
		c.Update()
		
		box = TBox(0,6.5,180,7.5)
		box.SetLineColor(3)
		box.SetLineWidth(2)
		box.Draw()
		
		c.Update()
		
		return hist,c,label,box
	
	def plot3x3FailQualityCodesVsPt(self):
		c = TCanvas('cFailQC3x3','Fail QC 3x3',910,0,900,700)
		c.SetLogz()
		hist = self.fileHandler.getHistogram('qualityCode/L1Muon3x3Fail_QcVsPt')
		hist.SetStats(0)
		hist.Scale(1,'width')
		hist.Draw('colz')
		c.Update()
		setupPalette(hist)
		c.Update()

		label = self.drawLabel()
		c.Update()

		return hist,c,label
		
	def plotQualityCodesStacked(self,gridSize):
		gridString = self.gridSizeDict.get(gridSize)
		c = TCanvas('cQualityCodes' + gridString + 'Stacked','Stacked QC ' + gridString,600,0,800,600)
		c.cd().SetBottomMargin(0.15)
		c.cd().SetRightMargin(0.20)
		qualityCodes = self.fileHandler.getHistogram('multiplicity/L1MuonQualityCodes' + gridString + '_Multiplicity')
		qualityCodesFail = self.fileHandler.getHistogram('multiplicity/L1MuonQualityCodes' + gridString + 'Fail_Multiplicity')
		
		countQualityCodes = self.fileHandler.getHistogram('multiplicity/L1MuonAllQualityCodes_Multiplicity')
		
		self.commandLine.output('Sanity check for quality code counts')
		for i in range(1,8):
			nTotalHistogram = countQualityCodes.GetBinContent(countQualityCodes.FindBin(i))
			nFail = qualityCodesFail.GetBinContent(qualityCodesFail.FindBin(i))
			nPass = qualityCodes.GetBinContent(qualityCodes.FindBin(i))
			nSummed = nFail + nPass
			print
			self.commandLine.output('NTotal: %d\t\tNSummed: %d' % (nTotalHistogram,nSummed))
			self.commandLine.output('Sanity check: %s'% (CliColors.OKBLUE + 'OK' + CliColors.ENDC if nTotalHistogram == nSummed else CliColors.FAIL + 'FAIL' + CliColors.ENDC) )
			print
			if nTotalHistogram:
				qualityCodes.SetBinContent(qualityCodes.FindBin(i),nPass/float(nTotalHistogram))
				qualityCodesFail.SetBinContent(qualityCodesFail.FindBin(i),nFail/float(nTotalHistogram))
		
		stack = THStack("hstack",";;rel. fraction")
		
		
		qualityCodes.SetLineColor(colorRwthDarkBlue)
		qualityCodes.SetFillColor(colorRwthDarkBlue)
		qualityCodes.SetFillStyle(3002)
	
		qualityCodesFail.SetFillColor(colorRwthMagenta)
		qualityCodesFail.SetLineColor(colorRwthMagenta)
		qualityCodesFail.SetFillStyle(3002)
	
		stack.Add(qualityCodes)
		stack.Add(qualityCodesFail)
		
		stack.Draw()
		stack.GetXaxis().SetRangeUser(0,8)

			#Label the bins with the meaning of the quality code
		for i in range(1,8):
			stack.GetXaxis().SetBinLabel(stack.GetXaxis().FindBin(i),self.qualityCodeDict.get(i))
			
		legend = TLegend(0.1,0.8,0.3,0.95)
		legend.AddEntry(qualityCodes,"Passed","f")
		legend.AddEntry(qualityCodesFail,"Failed","f")
		legend.Draw()
				
		setupAxes(stack)
		stack.SetMinimum(0.75)

		c.Update()
		
		self.storeCanvas(c, 'qualityCodesStacked' + gridString)
		
		return stack,c,qualityCodes,qualityCodesFail,legend
	
	def plotAllQualitiyCodes(self):	
		canvas = TCanvas('allQCCodes',"All QC")
		canvas.SetLogy()
		canvas.cd().SetBottomMargin(0.15)
	
		histAllCodes = self.fileHandler.getHistogram('multiplicity/L1MuonAllQualityCodes_Multiplicity')
	
		setupAxes(histAllCodes)
	
		histAllCodes.SetLineWidth(3)
		histAllCodes.SetLineColor(colorRwthDarkBlue)
		histAllCodes.Scale(1/histAllCodes.Integral())
		histAllCodes.GetXaxis().SetRangeUser(0,8)	
		histAllCodes.SetStats(0)
		histAllCodes.SetTitle(';;rel. fraction')
		#Label the bins with the meaning of the quality code
		for i in range(1,8):
			histAllCodes.GetXaxis().SetBinLabel(histAllCodes.GetXaxis().FindBin(i),self.qualityCodeDict.get(i))
		
		histAllCodes.Draw()
		
		histAllCodesTruth = None 
		if not self.data:
			histAllCodesTruth = self.fileHandler.getHistogram('multiplicity/L1MuonTruthAllQualityCodes_Multiplicity')
			histAllCodesTruth.Scale(1/histAllCodesTruth.Integral())
			histAllCodesTruth.SetLineWidth(3)
			histAllCodesTruth.SetLineColor(colorRwthMagenta)
			histAllCodesTruth.Draw('Same')
			
		label = self.drawLabel()
		
		legend = TLegend(0.1,0.8,0.3,0.95)
		legend.AddEntry(histAllCodes,"All L1","l")
		if not self.data: legend.AddEntry(histAllCodesTruth,"L1 Truth","l")
		legend.Draw()
		
		canvas.Update()
		self.storeCanvas(canvas,'allQualityCodes')
				
		return canvas,legend, label, histAllCodes,histAllCodesTruth