from plotting.Plot import Plot

from ROOT import TCanvas,TFile,TPad
from plotting.Utils import getLegend
from plotting.PlotStyle import colorRwthDarkBlue, colorRwthTuerkis, colorRwthRot,\
	colorRwthGruen, setupAxes, colorRwthLightBlue, colorRwthMagenta, colorRwthLila

MIP_PEAK_POSITION = 1.1

class Energy(Plot):
	#Initialize
	def __init__(self,filename,data,debug):
		Plot.__init__(self,filename,data,debug)
		self.createPlotSubdir('energy')
	
	def plotHoEnergyPerWheel(self):	
		res = self.plotEnergyPerWheel('horeco')
		res[-2].SetName('cEPerWheelHo')
		res[0].SetTitle('Reconstructed Energy per Wheel')
		res[-2].Update()
		self.storeCanvas(res[-2],'energyPerWheel')
		return res
	
	def plotMatchedHoEnergyPerWheel(self):
		res = self.plotEnergyPerWheel('L1MuonWithHoMatchAboveThr')
		res[-2].SetName('cEPerWheelMatchedHo')
		res[0].SetTitle('Reconstructed Energy per Wheel for matched HO')
		res[-2].Update()
		self.storeCanvas(res[-2],'energyPerWheelMatched')
		return res
	
	def plotMatchedAndNotMatchedPerWheel(self):
		resHo = self.plotEnergyPerWheel('horeco')
		resHoMatched = self.plotEnergyPerWheel('L1MuonWithHoMatchAboveThr')
		
		cTogether = TCanvas('cTogether','Matched and not Matched',1800,500)
		cTogether.Divide(4,1)
		
		plotTitles = [
					'Wheel -1',
					'Wheel 0 (-)',
					'Wheel 0 (+)',
					'Wheel +1',
					]
		
		for i in range(0,4):
			cTogether.cd(i+1).SetLogy()
			resHo[i].SetTitle(plotTitles[i] + ';E_{Rec} / GeV;fraction of MIP peak')
			#resHo[i].GetXaxis().SetRangeUser(-1,6)
			maxBin = resHoMatched[i].GetMaximumBin()#resHo[i].FindBin()
			resHo[i].Scale(1/resHo[i].GetBinContent(maxBin))
			resHo[i].SetStats(0)
			resHo[i].Draw()
			setupAxes(resHo[i])
			#resHo[i].GetYaxis().SetRangeUser(4e-3,2)
			resHoMatched[i].Scale(1/resHoMatched[i].GetBinContent(maxBin))
			resHoMatched[i].Draw('same')
		cTogether.Update()
		cTogether.SaveAs('plots/energyPerWheelTogether.root')

		return cTogether,resHo, resHoMatched
	
	def plotEnergyPerWheel(self,sourceName):	
		hoM1 = self.fileHandler.getHistogram('energy/perWheel/' + sourceName + '_Energy_M1')
		hoM0 = self.fileHandler.getHistogram('energy/perWheel/' + sourceName + '_Energy_M0')
		hoP0 = self.fileHandler.getHistogram('energy/perWheel/' + sourceName + '_Energy_P0')
		hoP1 = self.fileHandler.getHistogram('energy/perWheel/' + sourceName + '_Energy_P1')
		
		c = TCanvas('cEPerWheel','E Per Wheel')
		c.SetLogy()
		hoM1.SetLineColor(colorRwthDarkBlue)
#		hoM1.GetXaxis().SetRangeUser(-.8,6)
		hoM1.SetStats(0)
		hoM1.SetTitle('Reconstructed Energy per Wheel;E_{Rec} / GeV;# Entries')
		
		hoM0.SetLineColor(colorRwthTuerkis)
		hoP0.SetLineColor(colorRwthGruen)
		hoP1.SetLineColor(colorRwthMagenta)
		
		setupAxes(hoM1)
		
		hoM1.Draw('')
		hoM0.Draw('same')
		hoP0.Draw('same')
		hoP1.Draw('same')
		
		legend = getLegend(x1=.7,y2=0.9,y1=.6)
		legend.AddEntry(hoM1,'Wheel M1','l')
		legend.AddEntry(hoM0,'Wheel M0','l')
		legend.AddEntry(hoP0,'Wheel P0','l')
		legend.AddEntry(hoP1,'Wheel P1','l')
		legend.Draw()
		
		label = self.drawLabel()
		
		c.Update()
				
		return hoM1, hoM0, hoP0, hoP1, legend, c, label
		
	def plotEnergyNormalizedToMip(self):
		ho = self.fileHandler.getHistogram("energy/horeco_Energy")
		L1MuonAndHoMatch = self.fileHandler.getHistogram('energy/L1MuonWithHoMatch_Energy')
		L1MuonAndHoMatchAboveThr = self.fileHandler.getHistogram('energy/L1MuonWithHoMatchAboveThr_Energy')
		L1MuonAndHoMatchAboveThrFilt = self.fileHandler.getHistogram('energy/L1MuonWithHoMatchAboveThrFilt_Energy')
	
		canv = TCanvas("cEnergyNormToMip",'Energy Norm To MIP',1200,1200)
		canv.SetLogy()
	
		ho.SetStats(0)
		ho.SetTitle('Normalized energy distribution of HO hits (to MIP Maximum)')
		ho.GetXaxis().SetTitle('Reconstructed HO energy / GeV')
		ho.GetYaxis().SetTitle('rel. fraction')
		#ho.GetXaxis().SetRangeUser(-2,6)
	
		ho.SetLineColor(colorRwthDarkBlue)
		ho.SetLineWidth(3)
		ho.Scale(1/ho.GetBinContent(ho.FindBin(MIP_PEAK_POSITION)))
		ho.Draw()
		
		label = self.drawLabel()
		
		legend = getLegend(0.5,0.65,0.9,0.9)
		legend.AddEntry(ho,'All HO hits','l')
		legend.Draw()
	
		if(L1MuonAndHoMatch):
			L1MuonAndHoMatch.SetLineColor(colorRwthTuerkis)
			L1MuonAndHoMatch.SetLineWidth(3)
			L1MuonAndHoMatch.Scale(1/L1MuonAndHoMatch.GetBinContent(L1MuonAndHoMatch.FindBin(MIP_PEAK_POSITION)))
			L1MuonAndHoMatch.Draw('same')
			legend.AddEntry(L1MuonAndHoMatch,'L1Muon + HO match','l')
			
		if(L1MuonAndHoMatchAboveThr):
			L1MuonAndHoMatchAboveThr.SetLineColor(colorRwthRot)
			L1MuonAndHoMatchAboveThr.SetLineWidth(3)
			L1MuonAndHoMatchAboveThr.Scale(1/L1MuonAndHoMatchAboveThr.GetBinContent(L1MuonAndHoMatchAboveThr.FindBin(MIP_PEAK_POSITION)))
			L1MuonAndHoMatchAboveThr.Draw('same')
			legend.AddEntry(L1MuonAndHoMatchAboveThr,'L1Muon + HO match > 0.2 GeV','l')
	
	
		if(L1MuonAndHoMatchAboveThrFilt):
			L1MuonAndHoMatchAboveThrFilt.SetLineColor(colorRwthGruen)
			L1MuonAndHoMatchAboveThrFilt.SetLineWidth(3)
			L1MuonAndHoMatchAboveThrFilt.Scale(1/L1MuonAndHoMatchAboveThrFilt.GetBinContent(L1MuonAndHoMatchAboveThrFilt.FindBin(MIP_PEAK_POSITION)))
			L1MuonAndHoMatchAboveThrFilt.Draw('same')
			legend.AddEntry(L1MuonAndHoMatchAboveThrFilt,'L1Muon + HO match > 0.2 GeV (In Ho Geom.)','l')
	
		
	
		self.storeCanvas(canv,'energyNormToMip')
		canv.SaveAs("plots/energy/energyNormToMip.pdf")
	
		return canv,ho,L1MuonAndHoMatch, L1MuonAndHoMatchAboveThr,L1MuonAndHoMatchAboveThrFilt,label, legend
			
	def plotEnergyNormalized(self):

		ho = self.fileHandler.getHistogram("energy/horeco_Energy")
		hoNoise = self.fileHandler.getHistogram("energy/L1MuonPlusPi_Energy")
		L1MuonAndHoMatch = self.fileHandler.getHistogram('energy/L1MuonWithHoMatch_Energy')
		L1MuonAndHoMatchAboveThr = self.fileHandler.getHistogram('energy/L1MuonWithHoMatchAboveThr_Energy')
		L1MuonAndHoMatchAboveThrFilt = self.fileHandler.getHistogram('energy/L1MuonWithHoMatchAboveThrFilt_Energy')
	
		L1MuonAndHoMatch = None
		L1MuonAndHoMatchAboveThr = None
		L1MuonAndHoMatchAboveThrFilt = None
		
		canv = TCanvas("energieNormCanvas",'Energy Norm canvas',800,800)
		canv.SetLogy()
	
		ho.SetStats(0)
		ho.SetTitle('')#Normalized energy distribution of HO hits')
		ho.GetXaxis().SetTitle('Reconstructed HO energy / GeV')
		ho.GetYaxis().SetTitle('rel. fraction / 0.05 GeV')
		ho.GetXaxis().SetRangeUser(-0.2,2)
	
		ho.SetLineColor(colorRwthDarkBlue)
		ho.SetLineWidth(3)
		ho.Scale(1/ho.Integral())
		
		hoNoise.SetLineColor(colorRwthMagenta)
		hoNoise.SetLineWidth(3)
		hoNoise.Scale(1/hoNoise.Integral())
		
		setupAxes(ho)		
		ho.Draw()
		hoNoise.Draw('same')
				
		legend = getLegend(0.4,0.8,0.98,0.95)
		legend.AddEntry(ho,'All HO hits','l')
		legend.AddEntry(hoNoise,'HO hits, no #mu expected','l')
		legend.Draw()
	
		if(L1MuonAndHoMatch):
			L1MuonAndHoMatch.SetLineColor(colorRwthTuerkis)
			L1MuonAndHoMatch.SetLineWidth(3)
			L1MuonAndHoMatch.Scale(1/L1MuonAndHoMatch.Integral())
			L1MuonAndHoMatch.Draw('same')
			legend.AddEntry(L1MuonAndHoMatch,'L1Muon + HO match','l')
			
		if(L1MuonAndHoMatchAboveThr):
			L1MuonAndHoMatchAboveThr.SetLineColor(colorRwthRot)
			L1MuonAndHoMatchAboveThr.SetLineWidth(3)
			L1MuonAndHoMatchAboveThr.Scale(1/L1MuonAndHoMatchAboveThr.Integral())
			L1MuonAndHoMatchAboveThr.Draw('same')
			legend.AddEntry(L1MuonAndHoMatchAboveThr,'L1Muon + HO match > 0.2 GeV','l')
	
	
		if(L1MuonAndHoMatchAboveThrFilt):
			L1MuonAndHoMatchAboveThrFilt.SetLineColor(colorRwthGruen)
			L1MuonAndHoMatchAboveThrFilt.SetLineWidth(3)
			L1MuonAndHoMatchAboveThrFilt.Scale(1/L1MuonAndHoMatchAboveThrFilt.Integral())
			L1MuonAndHoMatchAboveThrFilt.Draw('same')
			legend.AddEntry(L1MuonAndHoMatchAboveThrFilt,'L1Muon + HO match > 0.2 GeV (In Ho Geom.)','l')
	
		setupAxes(ho)
		canv.Update()
		
		ho.GetYaxis().SetTitleOffset(1.)
		canv.Update()
		
		self.storeCanvas(canv,'energyNorm',marginRight=.02,marginLeft=.15)

		return [canv,ho,L1MuonAndHoMatch, L1MuonAndHoMatchAboveThr,L1MuonAndHoMatchAboveThrFilt,hoNoise,legend]
	
	def plotEnergy(self):

		ho = self.fileHandler.getHistogram("energy/horeco_Energy")
		L1MuonAndHoMatch = self.fileHandler.getHistogram('energy/L1MuonWithHoMatch_Energy')
		L1MuonAndHoMatchAboveThr = self.fileHandler.getHistogram('energy/L1MuonWithHoMatchAboveThr_Energy')
		L1MuonAndHoMatchAboveThrFilt = self.fileHandler.getHistogram('energy/L1MuonWithHoMatchAboveThrFilt_Energy')
	
		canv = TCanvas("energieCanvas",'Energy canvas',1200,1200)
		canv.SetLogy()
	
		ho.SetStats(0)
		ho.SetTitle('Energy distribution of HO hits')
		ho.GetXaxis().SetTitle('Reconstructed HO energy / GeV')
		ho.GetYaxis().SetTitle('N')
		ho.GetXaxis().SetRangeUser(-2,6)
	
		ho.SetLineColor(colorRwthDarkBlue)
		ho.SetLineWidth(3)
		ho.Draw()
		
		label = self.drawLabel()
		
		legend = getLegend(0.5,0.65,0.9,0.9)
		legend.AddEntry(ho,'All HO hits','l')
		legend.Draw()
	
		if(L1MuonAndHoMatch):
			L1MuonAndHoMatch.SetLineColor(colorRwthTuerkis)
			L1MuonAndHoMatch.SetLineWidth(3)
			L1MuonAndHoMatch.Draw('same')
	#		legend.AddEntry(L1MuonAndHoMatch,'L1Muon + HO match','l')
			
		if(L1MuonAndHoMatchAboveThr):
			L1MuonAndHoMatchAboveThr.SetLineColor(colorRwthRot)
			L1MuonAndHoMatchAboveThr.SetLineWidth(3)
			L1MuonAndHoMatchAboveThr.Draw('same')
			legend.AddEntry(L1MuonAndHoMatchAboveThr,'L1Muon + HO match > 0.2 GeV','l')
	
	
		if(L1MuonAndHoMatchAboveThrFilt):
			L1MuonAndHoMatchAboveThrFilt.SetLineColor(colorRwthGruen)
			L1MuonAndHoMatchAboveThrFilt.SetLineWidth(3)
			L1MuonAndHoMatchAboveThrFilt.Draw('same')
			legend.AddEntry(L1MuonAndHoMatchAboveThrFilt,'L1Muon + HO match > 0.2 GeV (In Ho Geom.)','l')
	
		
		self.storeCanvas(canv,'energy')
	
		#f = TFile.Open("plots/energy/energy.root","RECREATE")
		#canv.Write()
		#f.Close()
		return [canv,ho,L1MuonAndHoMatch, L1MuonAndHoMatchAboveThr,L1MuonAndHoMatchAboveThrFilt,label,legend]
	
	def plotEnergyVsEta(self,sourceHistogram = 'L1MuonWithHoMatch_EnergyVsEta'):
	
		canv = TCanvas("energieCanvas",'Energy canvas',1200,1200)
	
		energyVsEta = self.fileHandler.getHistogram("energy/" + sourceHistogram)
		energyVsEta.Rebin2D(10,1)
		energyVsEta.GetXaxis().SetRangeUser(-1.1,1.1)
		energyVsEta.GetYaxis().SetRangeUser(0,2.5)
		energyVsEta.SetStats(0)
	
		energyVsEta.GetZaxis().SetTitle('Entries / 0.05 GeV')
		
		energyVsEta.Draw('colz')
		
		canv.Update()
		pal = energyVsEta.GetListOfFunctions().FindObject("palette")
		pal.SetX2NDC(0.92)
		
		self.storeCanvas(canv,'energyVsEta')
	
# 		f = TFile.Open("plots/energy/energyVsEta.root","RECREATE")
# 		canv.Write()
# 		f.Close()
		return canv
	
	def plotEnergyVsPhi(self,sourceHistogram = 'L1MuonWithHoMatch_EnergyVsPhi'):
	
		canv = TCanvas("energieCanvas",'Energy canvas',1200,1200)
	
		energyVsEta = self.fileHandler.getHistogram("energy/" + sourceHistogram)
		energyVsEta.Rebin2D(10,1)
		energyVsEta.GetXaxis().SetRangeUser(-3.17,3.17)
		energyVsEta.GetYaxis().SetRangeUser(0,2.5)
		energyVsEta.SetStats(0)
	
		energyVsEta.GetZaxis().SetTitle('Entries / 0.05 GeV')
		
		energyVsEta.Draw('colz')
		
		canv.Update()
		pal = energyVsEta.GetListOfFunctions().FindObject("palette")
		pal.SetX2NDC(0.92)
		

		self.storeCanvas(canv,'energyVsPhi')
	
# 		f = TFile.Open("plots/energy/energyVsPhi.root","RECREATE")
# 		canv.Write()
# 		f.Close()
		return canv
	
	# Generate the plot based on the energy vs eta and phin histogram
	def plotEnergyVsEtaPhi(self,sourceHistogram = 'L1MuonWithHoMatch_EnergyVsEtaPhi'):

		canv = TCanvas("energieVsPositionCanvas",'Energy canvas',1200,1200)

		energyVsPos = self.fileHandler.getHistogram("energy/" + sourceHistogram)
		projection = energyVsPos.Project3DProfile()
		projection.GetXaxis().SetTitle(energyVsPos.GetXaxis().GetTitle())
		projection.GetYaxis().SetTitle(energyVsPos.GetYaxis().GetTitle())
		projection.GetZaxis().SetTitle(energyVsPos.GetZaxis().GetTitle())
		projection.Draw('colz')
		
		canv.Update()
		pal = projection.GetListOfFunctions().FindObject("palette")
		pal.SetX2NDC(0.92)
		
		stats = projection.GetListOfFunctions().FindObject("stats")
		stats.SetX1NDC(.1)
		stats.SetX2NDC(.3)
		stats.SetY1NDC(.7)
		stats.SetY2NDC(.9)
		
		canv.SaveAs("plots/energy/energyVsPosition.gif")
		canv.SaveAs("plots/energy/energyVsPosition.pdf")
	
		f = TFile.Open("plots/energy/energyVsPosition.root","RECREATE")
		canv.Write()
		f.Close()
		return canv
	